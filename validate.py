# %%

# Validation

import csv
import math
#pysam: https://pysam.readthedocs.io/en/latest/usage.html
import pysam
#print(pysam.__version__)

import numpy as np 
import matplotlib.pyplot as plt  

# %%
#input

output_dir = sys.argv[1] + "/"
vcf_file = sys.argv[2]
if_hg38_input = sys.argv[3]
centromere_file = sys.argv[4]
exclude_assem1_non_cover_file = sys.argv[5]
exclude_assem2_non_cover_file = sys.argv[6]
exclude_high_depth_file = sys.argv[7]

# %%
#constants

interval = 20

if_hg38 = False
if if_hg38_input == "True":
    if_hg38 = True

# %%
#build centromere position dictionary

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


# %%
#build lists for excluded SV positions

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

with open(exclude_high_depth_file) as f:
#with open("/panfs/qcb-panasas/jianzhiy/illuVeri/use_mcutils/output/v4.6.1/lumpy/HG002/exclude_high_depth.bed") as f:
    reader = csv.reader(f, delimiter="\t")
    exclude_high_depth = list(reader)
f.close()

#check if current cases is in the excluded list
def check_exclude(list_to_check):
	if (list_to_check in exclude_assem1_non_cover) or (list_to_check in exclude_assem2_non_cover):
		return True
    #With overlapping contigs trimmed, 
	#elif (list_to_check in exclude_assem1_short_reads) or (list_to_check in exclude_assem2_short_reads):
	#	return True
	elif list_to_check in exclude_high_depth:
		return True
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

def check_tp (rela_len, rela_score):
    result = True
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
    return result

#in align_info:
# SV_index score_before score_after no_found SV_type hap:1 or 2
# (counter alignment_beforeSV[0][2] alignment_afterSV[0][2] len(query_start_dic) sv_type hap)
# contig_seg_len ref_seg_len SV_len
# (len(query_frag) len(ref_frag) sv_end-sv_pos+1)
# chr_name SV_start SV_end interval_start interval_end
# (ref_name sv_pos sv_end ref_start ref_end)
# error message
# (err_mes)

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

def updateDict(dict_score, align_info):
	for record in align_info:
		#some cases are not aligned successfully due to function's (memory) issues: give -1 in that column
		if int(record[3]) > 0 and abs(int(record[8])) > 30 and str(record[9]) in chr_list:
			#filter out centromere cases
			index = str(record[0])
			ref_name = str(record[9])
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
			if check_exclude(list_to_check):
				continue

			#if ref start or ref end in centromere, skip this case
			if (ref_start > centro_start and ref_start < centro_end) or (ref_end > centro_start and ref_end < centro_end):
				continue
            
			#filter out call at high coverage location
			#ref_start_depth = get_depth(ref_name, ref_start, bam_file)
			#ref_end_depth = get_depth(ref_name, ref_end, bam_file)
			#if ((ref_start_depth > 2 * avg_depth) or (ref_end_depth > 2 * avg_depth)):
			#	continue

			if record[0] not in dict_score:
				dict_score[record[0]] = record[1:len(record)]
			else:
				#TODO: may change
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
				#if INS or DEL choose the better relative length one:
				if record[4] == 'INS' or record[4] == 'DEL':
					if float(dict_score[record[0]][7]) == 0 or float(record[8]) == 0:
						continue
					else:
						old_rela_len = (float(dict_score[record[0]][5]) - float(dict_score[record[0]][6]))/(float(dict_score[record[0]][7]))
						new_rela_len = (float(record[6]) - float(record[7]))/(float(record[8]))
                    
              	      #test
						#if record[0] == "90":
						#	print(new_rela_len, old_rela_len)
                        
						if abs(new_rela_len - 1) < abs(old_rela_len - 1):
							dict_score.update({record[0]: record[1:len(record)]})
                	#if INV choose the better relative score one:
				elif record[4] == 'INV':
					if float(dict_score[record[0]][0]) == 0 or float(record[1]) == 0:
						continue
					else:
						old_rela_score = (float(dict_score[record[0]][1]) - float(dict_score[record[0]][0]))/abs(float(dict_score[record[0]][0]))
						new_rela_score = (float(record[2]) - float(record[1]))/abs(float(record[1]))
                    
              	      #test
						#if record[0] == "90":
						#	print(new_rela_len, old_rela_len)
                        
						if old_rela_score < new_rela_score:
							dict_score.update({record[0]: record[1:len(record)]})
	return dict_score



# %%

#validate by both haplotypes
with open(output_dir + "align_info_assem1_chrall.txt") as f:
	reader = csv.reader(f, delimiter="\t")
	align_info_assem1 = list(reader)
f.close()

with open(output_dir + "align_info_assem2_chrall.txt") as f:
	reader = csv.reader(f, delimiter="\t")
	align_info_assem2 = list(reader)
f.close()

#dicts for single haplotypes and combined
dict_comb = dict()
dict_comb = updateDict(dict_comb, align_info_assem1)
dict_comb = updateDict(dict_comb, align_info_assem2)

align_plot_info_1 = []
for record in align_info_assem1:
	if record[0] in dict_1:
		#TODO: handle this case
		#if zero score_before
		if float(dict_1[record[0]][0]) == 0 or float(record[8]) == 0:
			continue
		#round to 2 decimal points
		#print(float(record[8]))
		rela_len = round((float(record[6]) - float(record[7]))/float(record[8]), 2)
		rela_score = round((float(dict_1[record[0]][1]) - float(dict_1[record[0]][0]))/abs(float(dict_1[record[0]][0])), 2)
		align_plot_info_1.append([rela_len, rela_score, int(record[0])])

align_plot_info_2 = []
for record in align_info_assem2:
	if record[0] in dict_2:
		#TODO: handle this case
		#if zero score_before
		if float(dict_2[record[0]][0]) == 0 or float(record[8]) == 0:
			continue
		#round to 2 decimal points
		#print(float(record[8]))
		rela_len = round((float(record[6]) - float(record[7]))/float(record[8]), 2)
		rela_score = round((float(dict_2[record[0]][1]) - float(dict_2[record[0]][0]))/abs(float(dict_2[record[0]][0])), 2)
		align_plot_info_2.append([rela_len, rela_score, int(record[0])])

#combine two haplotypes
#TODO: debug this part: check the results!!!
align_plot_info_comb = []
for record in dict_comb:
		#TODO: solve 0 score problem: should be included!!!
		#if zero score
		if float(dict_comb[record][0]) == 0 or float(dict_comb[record][7]) == 0:
			#rela_score = 0.01
			continue
		else:
			#round to 2 decimal points
			rela_score = round((float(dict_comb[record][1]) - float(dict_comb[record][0]))/abs(float(dict_comb[record][0])), 2)
		rela_len = round((float(dict_comb[record][5]) - float(dict_comb[record][6]))/float(dict_comb[record][7]), 2)

		align_plot_info_comb.append([rela_len, rela_score, int(record)])

        
        
        
#output tp/fp as vcf files

