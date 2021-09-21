# %%

# Validation

import csv
import math
#pysam: https://pysam.readthedocs.io/en/latest/usage.html
import pysam
#print(pysam.__version__)

import numpy as np 
import matplotlib.pyplot as plt  
import sys


#build centromere position dictionary
def build_centro_dict(centromere_file):
    #centromere file
    centromere_raw = []
    with open(centromere_file) as f:
        reader = csv.reader(f, delimiter="\t")
        centromere_raw = list(reader)
    f.close() 

    #build dictionay
    dict_centromere = dict()
    pre_centro_chr = ""
    start_centro_pos = ""
    for record in centromere_raw:
        if record[0] == pre_centro_chr:
            end_centro_pos = str(record[2])
            dict_centromere[record[0]] = (start_centro_pos, end_centro_pos)
        else:
            start_centro_pos = str(record[1])
            pre_centro_chr = record[0]
    return dict_centromere

#build lists for excluded SV positions
def get_filtered_sv_pos(exclude_assem1_non_cover_file, exclude_assem2_non_cover_file):

    with open(exclude_assem1_non_cover_file) as f:
        reader = csv.reader(f, delimiter="\t")
        exclude_assem1_non_cover = list(reader)
    f.close()

    with open(exclude_assem2_non_cover_file) as f:
        reader = csv.reader(f, delimiter="\t")
        exclude_assem2_non_cover = list(reader)
    f.close()

    '''
    with open(output_dir + "exclude_assem1_short_reads_250000.bed") as f:
        reader = csv.reader(f, delimiter="\t")
        exclude_assem1_short_reads = list(reader)
    f.close()

    with open(output_dir + "exclude_assem2_short_reads_250000.bed") as f:
        reader = csv.reader(f, delimiter="\t")
        exclude_assem2_short_reads = list(reader)
    f.close()
    '''

#     with open(exclude_high_depth_file) as f:
#     #with open("/panfs/qcb-panasas/jianzhiy/illuVeri/use_mcutils/output/v4.6.1/lumpy/HG002/exclude_high_depth.bed") as f:
#         reader = csv.reader(f, delimiter="\t")
#         exclude_high_depth = list(reader)
#     f.close()
    
    return exclude_assem1_non_cover, exclude_assem2_non_cover

#check if current cases is in the excluded list
def check_exclude(list_to_check, exclude_assem1_non_cover, exclude_assem2_non_cover):
	if (list_to_check in exclude_assem1_non_cover) or (list_to_check in exclude_assem2_non_cover):
		return True
    #With overlapping contigs trimmed, 
	#elif (list_to_check in exclude_assem1_short_reads) or (list_to_check in exclude_assem2_short_reads):
	#	return True
# 	elif list_to_check in exclude_high_depth:
# 		return True
	else:
		return False
    
#check if current cases is in the excluded list for male chr X
def check_exclude_chrx(list_to_check, exclude_assem1_non_cover, exclude_assem2_non_cover):
	if (list_to_check in exclude_assem1_non_cover) and (list_to_check in exclude_assem2_non_cover):
		return True
    #With overlapping contigs trimmed, 
	#elif (list_to_check in exclude_assem1_short_reads) or (list_to_check in exclude_assem2_short_reads):
	#	return True
# 	elif list_to_check in exclude_high_depth:
# 		return True
	else:
		return False

#get depth of coverage given position of a chromosome
def get_depth(ref_name, ref_pos, bam_file):
    pos_str = ref_name + ':' + str(int(ref_pos) - 1) + '-' + str(ref_pos)
    res = pysam.depth("-r", pos_str, bam_file)
    if res=='':
        return 0
    start = 0
    end = len(res) - 1
    for i in range(len(res) - 1, -1, -1):
        if res[i] == '\t':
            start = i + 1
            break
    return int(res[start:end])

#check if true positive or not
def check_tp(rela_len, rela_score, sv_type):
    result = True
    if sv_type in ['DEL', 'DUP', 'DUP:TANDEM']:
        if rela_score >= 0 and rela_score <= 2.5:
            if rela_len >= -0.05*rela_score + 0.8 and rela_len <= 0.05*rela_score + 1.2:
                result = True
            else:
                result = False
        elif rela_score > 2.5:
            if rela_len >= 0.675 and rela_len <= 1.325:
                result = True
            else:
                result = False
        else:
            result = False
    elif sv_type == 'INS':
        if rela_len < 0.675 or rela_len > 1.325:
            result = False
    elif sv_type == 'INV':
        if rela_score <= 0:
            result = False
    return result

# build dictionary for validation
def updateDict(dict_score, align_info, exclude_assem1_non_cover, exclude_assem2_non_cover, dict_centromere, chr_list, if_hg38):
    for record in align_info:
        #len > 30
        if int(record[3]) > 0 and abs(int(record[8])) > 30 and str(record[9]) in chr_list:
#         if int(record[3]) > 0 and abs(int(record[8])) >= 10 and abs(int(record[8])) <= 30 and str(record[9]) in chr_list:
            #filter out centromere cases
            index = str(record[0])
            ref_name = str(record[9])
            #interval start and end
            ref_start = int(record[12])
            ref_end = int(record[13])

            if if_hg38:
                centro_start = int(dict_centromere[ref_name][0])
                centro_end = int(dict_centromere[ref_name][1])
            else:
                centro_start = int(dict_centromere['chr'+ref_name][0])
                centro_end = int(dict_centromere['chr'+ref_name][1])

            #if SV in the exclude_list: start or end of a contig
            sv_pos = int(record[10])
            sv_end = int(record[11])
            list_to_check = [str(ref_name), str(sv_pos), str(sv_end)]
            #if sv in high-depth regions or non-covered regions, skip
            if check_exclude(list_to_check, exclude_assem1_non_cover, exclude_assem2_non_cover):
                continue

            #if ref start or ref end in centromere, skip
            if (ref_start > centro_start and ref_start < centro_end) or (ref_end > centro_start and ref_end < centro_end):
                continue

            #start to build a dictionay: merge validation information from 2 haplotypes
            #TODO: solve the score/length = 0 problem
            #skip before score/length = 0
            if float(record[1]) == 0 or float(record[8]) == 0:
                continue

            if record[0] not in dict_score:
                dict_score[record[0]] = record[1:len(record)]
            else:
                #TODO: solve the score = 0 problem
                #choose the better relative score one:
                #if one of the before score is 0, choose the better relative length one:
                '''
                if float(dict_score[record[0]][0]) == 0 or float(record[1]) == 0:
                    old_rela_len = (float(dict_score[record[0]][5]) - float(dict_score[record[0]][6]))/abs(float(dict_score[record[0]][7]))
                    new_rela_len = (float(record[6]) - float(record[7]))/abs(float(record[8]))
                    if new_rela_len < old_rela_len:
                        dict_score.update({record[0]: record[1:len(record)]})
                    continue
                old_rela_change = (float(dict_score[record[0]][1]) - float(dict_score[record[0]][0]))/abs(float(dict_score[record[0]][0]))
                new_rela_change = (float(record[2]) - float(record[1]))/abs(float(record[1]))
                if new_rela_change > old_rela_change:
                    dict_score.update({record[0]: record[1:len(record)]})
                '''
                #if INS or DEL choose the tp one
                #if both tp/fp, choose the better relative length one
                if record[4] in ['INS', 'DEL', 'DUP', 'DUP:TANDEM']:
                    old_rela_score = (float(dict_score[record[0]][1]) - float(dict_score[record[0]][0]))/abs(float(dict_score[record[0]][0]))
                    new_rela_score = (float(record[2]) - float(record[1]))/abs(float(record[1]))
                    old_rela_len = (float(dict_score[record[0]][5]) - float(dict_score[record[0]][6]))/(float(dict_score[record[0]][7]))
                    new_rela_len = (float(record[6]) - float(record[7]))/(float(record[8]))

                    old_res = check_tp(old_rela_len, old_rela_score, record[4])
                    new_res = check_tp(new_rela_len, new_rela_score, record[4])

                    if new_res and not old_res:
                        dict_score.update({record[0]: record[1:len(record)]})
                    elif old_res and not new_res:
                        continue
                    else:
                        if abs(new_rela_len - 1) < abs(old_rela_len - 1):
                            dict_score.update({record[0]: record[1:len(record)]})
                #if INV choose the better relative score one:
                elif record[4] == 'INV':
                    if float(dict_score[record[0]][0]) == 0 or float(record[1]) == 0:
                        continue
                    else:
                        old_rela_score = (float(dict_score[record[0]][1]) - float(dict_score[record[0]][0]))/abs(float(dict_score[record[0]][0]))
                        new_rela_score = (float(record[2]) - float(record[1]))/abs(float(record[1]))

                        if old_rela_score < new_rela_score:
                            dict_score.update({record[0]: record[1:len(record)]})
    return dict_score

#validate by both haplotypes
def vali_info(output_dir, exclude_assem1_non_cover, exclude_assem2_non_cover, assem1_info_file, assem2_info_file, dict_centromere, chr_list, if_hg38):
    with open(output_dir + assem1_info_file) as f:
        reader = csv.reader(f, delimiter="\t")
        align_info_assem1 = list(reader)
    f.close()

    with open(output_dir + assem2_info_file) as f:
        reader = csv.reader(f, delimiter="\t")
        align_info_assem2 = list(reader)
    f.close()

    dict_comb = dict()
    dict_comb = updateDict(dict_comb, align_info_assem1, exclude_assem1_non_cover, exclude_assem2_non_cover, dict_centromere, chr_list, if_hg38)
    dict_comb = updateDict(dict_comb, align_info_assem2, exclude_assem1_non_cover, exclude_assem2_non_cover, dict_centromere, chr_list, if_hg38)

    return dict_comb


#TODO: output tp/fp as vcf files
#output a text file
#CHR POS END SVTYPE rela_len rela_score validation_res

#write output
def write_output(output_dir, dict_comb):
    g = open(output_dir + "ttmars_res.txt", "w")
    for record in dict_comb:
        #TODO: solve 0 score problem
        #if zero score
        if float(dict_comb[record][0]) == 0 or float(dict_comb[record][7]) == 0:
            #rela_score = 0.01
            continue
        rela_score = round((float(dict_comb[record][1]) - float(dict_comb[record][0]))/abs(float(dict_comb[record][0])), 2)
        rela_len = round((float(dict_comb[record][5]) - float(dict_comb[record][6]))/float(dict_comb[record][7]), 2)

        g.write(str(dict_comb[record][8]) + "\t")
        g.write(str(dict_comb[record][9]) + "\t")
        g.write(str(dict_comb[record][10]) + "\t")
        g.write(str(dict_comb[record][3]) + "\t")
        g.write(str(rela_len) + "\t")
        g.write(str(rela_score) + "\t")
        g.write(str(check_tp(rela_len, rela_score, str(dict_comb[record][3]))))
        g.write("\n")
    g.close()

#main function
def main():
    #input
    output_dir = sys.argv[1] + "/"
    if_hg38_input = sys.argv[2]
    centromere_file = sys.argv[3]
    #exclude_assem1_non_cover_file = sys.argv[4]
    #exclude_assem2_non_cover_file = sys.argv[5]
    #exclude_high_depth_file = sys.argv[6]
    
    #constants
    interval = 20
    if_hg38 = False
    if if_hg38_input == "True":
        if_hg38 = True
    chr_list = []
    if if_hg38:
        chr_list = ["chr1", "chr2", "chr3", "chr4", "chr5",
                    "chr6", "chr7", "chr8", "chr9", "chr10",
                    "chr11", "chr12", "chr13", "chr14", "chr15",
                    "chr16", "chr17", "chr18", "chr19", "chr20",
                    "chr21", "chr22", "chrX", "chrY"]
    else:
        chr_list = ["1", "2", "3", "4", "5",
                    "6", "7", "8", "9", "10",
                    "11", "12", "13", "14", "15",
                    "16", "17", "18", "19", "20",
                    "21", "22", "X", "Y"]
        
    #build centromere dictionary
    dict_centromere = build_centro_dict(centromere_file)
    
    #build lists for excluded SV positions
    exclude_assem1_non_cover, exclude_assem2_non_cover = get_filtered_sv_pos(output_dir + "exclude_assem1_non_cover.bed", 
                                                                             output_dir + "exclude_assem2_non_cover.bed")
    
    #validate by both haplotypes
    dict_comb = vali_info(output_dir, 
                          exclude_assem1_non_cover, 
                          exclude_assem2_non_cover, 
                          "align_info_assem1_chrall.txt",
                          "align_info_assem2_chrall.txt", 
                          dict_centromere, 
                          chr_list, 
                          if_hg38)
    
    #write output
    write_output(output_dir, dict_comb)

if __name__ == "__main__":
    main()
