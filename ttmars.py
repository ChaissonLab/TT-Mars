import sys
import csv
import pysam
import numpy as np 
import math

import get_conf_int
import validate
import get_align_info

##########################################################
##########################################################
##input
output_dir = sys.argv[1] + "/"
if_hg38_input = sys.argv[2]
centromere_file = sys.argv[3]
#assembly bam files
assem1_non_cov_regions_file = sys.argv[4]
assem2_non_cov_regions_file = sys.argv[5]
#avg_read_depth = sys.argv[6]
#read_bam_file = sys.argv[6]
vcf_file = sys.argv[6]
#ref fasta file
ref_file = sys.argv[7]
#assembly fasta files
query_file1 = sys.argv[8]
query_file2 = sys.argv[9]
liftover_file1 = sys.argv[10]
liftover_file2 = sys.argv[11]
tandem_file = sys.argv[12]
if_passonly_input = sys.argv[13]
seq_resolved_input = sys.argv[14]

# liftover_file1_0 = sys.argv[12]
# liftover_file2_0 = sys.argv[13]

##########################################################
##########################################################
#constants

#liftover interval
interval = 20
if_hg38 = False
if if_hg38_input == "True":
    if_hg38 = True
#if pass_only
if_pass_only = False
if if_passonly_input == "True":
    if_pass_only = True
#if seq_resolved
seq_resolved = False
if seq_resolved_input == "True":
    seq_resolved = True
#chr names
chr_list = []
if if_hg38:
    chr_list = ["chr1", "chr2", "chr3", "chr4", "chr5",
                "chr6", "chr7", "chr8", "chr9", "chr10",
                "chr11", "chr12", "chr13", "chr14", "chr15",
                "chr16", "chr17", "chr18", "chr19", "chr20",
                "chr21", "chr22", "chrX"]
else:
    chr_list = ["1", "2", "3", "4", "5",
                "6", "7", "8", "9", "10",
                "11", "12", "13", "14", "15",
                "16", "17", "18", "19", "20",
                "21", "22", "X"]
#approximate length of chromosomes
chr_len = [250000000, 244000000, 199000000, 192000000, 182000000, 
            172000000, 160000000, 147000000, 142000000, 136000000, 
            136000000, 134000000, 116000000, 108000000, 103000000, 
            90400000, 83300000, 80400000, 59200000, 64500000, 
            48200000, 51400000, 157000000, 59400000]

#max/min length of allowed SV not DUP
memory_limit = 100000
memory_min = 10
#max length of allowed DUP
dup_memory_limit = 50000
dup_memory_min = 10

#valid types
valid_types = ['DEL', 'INS', 'INV', 'DUP:TANDEM', 'DUP']

#tandem repeats regions file
with open(tandem_file) as f:
    reader = csv.reader(f, delimiter="\t")
    tandem_info = list(reader)
f.close()  

#get tandem start and end list
tandem_start_list, tandem_end_list = get_align_info.get_chr_tandem_shart_end_list(tandem_info, if_hg38)

##########################################################
##########################################################

#build lists for excluded SV positions

#Output regions on ref where its not covered by at least one of the assembly
# get_conf_int.get_non_cover_regions(output_dir, bam_file1, 1, chr_list)
# get_conf_int.get_non_cover_regions(output_dir, bam_file2, 2, chr_list)

#Get regions where read depth > 2 * avg_read_depth
#get_conf_int.get_high_depth_calls_info(output_dir, read_bam_file, vcf_file, avg_read_depth)

#Output sv positions
get_conf_int.get_sv_positions(output_dir, vcf_file)

#Output filtered calls in non-covered regions
SV_positions_file = output_dir + "SV_positions.bed"
# assem1_non_cov_regions_file = output_dir + "assem1_non_cov_regions.bed"
# assem2_non_cov_regions_file = output_dir + "assem2_non_cov_regions.bed"
get_conf_int.output_non_cov_call_info(output_dir, SV_positions_file, assem1_non_cov_regions_file, assem2_non_cov_regions_file)

#get filtered sv info, using results from get_conf_int.py 
exclude_assem1_non_cover, exclude_assem2_non_cover = validate.get_filtered_sv_pos(output_dir + "exclude_assem1_non_cover.bed", 
                                                                                  output_dir + "exclude_assem2_non_cover.bed")

#build centromere dictionary
dict_centromere = validate.build_centro_dict(centromere_file)
    


#return False if not filtered
#first_filter: type, PASS, chr_name
def first_filter(sv, sv_type):
    #type filter
    if sv_type not in valid_types:
        return True
    #PASS filter
    if if_pass_only:
        if 'PASS' not in sv.filter.keys():
            return True
    chr_name = sv.chrom
    #chr filter
    if chr_name not in chr_list:
        return True
    return False

#second_filter: centromere, non-cov
def second_filter(sv):
    index = sv.idx
    ref_name = sv.ref_name
    sv_pos = sv.sv_pos
    sv_stop = sv.sv_stop

    if if_hg38:
        centro_start = int(dict_centromere[ref_name][0])
        centro_end = int(dict_centromere[ref_name][1])
    else:
        centro_start = int(dict_centromere['chr'+ref_name][0])
        centro_end = int(dict_centromere['chr'+ref_name][1])

    #centromere
    if (sv_pos > centro_start and sv_pos < centro_end) or (sv_stop > centro_start and sv_stop < centro_end):
        sv.is_sec_fil = True
        return True
        
    #non-cov
    list_to_check = [str(ref_name), str(sv_pos), str(sv_stop)]
    #if sv in high-depth regions or non-covered regions, skip
    if validate.check_exclude(list_to_check, exclude_assem1_non_cover, exclude_assem2_non_cover):
        sv.is_sec_fil = True
        return True
    
#third_filter: size
def third_filter(sv):
    #size
    if sv.sv_type not in ['DUP:TANDEM', 'DUP']:
        if abs(sv.length) < memory_min or abs(sv.length) > memory_limit:
            sv.is_third_fil = True
            return True
    else:
        if abs(sv.length) < dup_memory_min or abs(sv.length) > dup_memory_limit:
            sv.is_third_fil = True
            return True
    

#get validation info
def write_vali_info(sv_list):
    g = open(output_dir + "ttmars_res.txt", "w")
    for sv in sv_list:
        #skip if not analyzed
        if (not sv.analyzed_hap1) or (not sv.analyzed_hap2):
            continue
        
        res = sv.get_vali_res()
        
        g.write(str(sv.ref_name) + "\t")
        g.write(str(sv.sv_pos) + "\t")
        g.write(str(sv.sv_stop) + "\t")
        g.write(str(sv.sv_type) + "\t")
        g.write(str(res[1]) + "\t")
        g.write(str(res[2]) + "\t")
        g.write(str(res[0]))
        g.write("\n")
    g.close()
    
#define class
class struc_var:
    def __init__(self, idx, ref_name, sv_type, sv_pos, sv_stop, length, gt):
        self.idx = idx
        self.ref_name = ref_name
        self.sv_pos = sv_pos
        self.sv_stop = sv_stop
        self.sv_type = sv_type
        self.length = length
        self.gt = gt
        #if the call is part of an aggregate SV
        self.is_agg = False
        #if second filtered out
        self.is_sec_fil = False
        self.is_third_fil = False
        
        self.query_name_hap1 = "NA"
        self.query_name_hap2 = "NA"
        
        self.ref_start_best_hap1 = -1
        self.ref_end_best_hap1 = -1
        self.query_start_best_hap1 = -1
        self.query_end_best_hap1 = -1
        
        self.ref_start_best_hap2 = -1
        self.ref_end_best_hap2 = -1
        self.query_start_best_hap2 = -1
        self.query_end_best_hap2 = -1
        
        self.analyzed_hap1 = False
        self.analyzed_hap2 = False
        
        self.len_query_hap1 = -1
        self.len_ref_hap1 = -1
        self.len_query_hap2 = -1
        self.len_ref_hap2 = -1
        
        self.score_before_hap1 = -1
        self.score_after_hap1 = -1
        self.score_before_hap2 = -1
        self.score_after_hap2 = -1
        
        self.neg_strand_hap1 = False
        self.neg_strand_hap2 = False
        
        self.ins_seq = ""
        self.if_seq_resolved = False
        
    def check_tp(self, rela_len, rela_score):
        result = True
        if self.sv_type in ['DEL', 'DUP', 'DUP:TANDEM']:
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
        elif self.sv_type == 'INS':
            #not seq-resolved
            #if len(self.ins_seq) == 0:
            if not self.if_seq_resolved:
                if rela_len < 0.675 or rela_len > 1.325:
                    result = False
            #seq-resolved
            else:
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
                
        elif self.sv_type == 'INV':
            if rela_score <= 0:
                result = False
        return result
        
    def print_info(self):
        print(self.idx, self.ref_name, self.sv_pos, self.sv_stop, self.sv_type, self.length, self.gt, self.is_agg, self.is_sec_fil, self.is_third_fil)
        
    def cal_rela_score(self, score_before, score_after):
        if score_before > -1 and score_before < 0:
            tmp_score_before = -1
            tmp_score_after = score_after + (tmp_score_before - score_before)
            return round((tmp_score_after - tmp_score_before) / abs(tmp_score_before), 2)
        
        elif score_before >= 0 and score_before < 1:
            tmp_score_before = 1
            tmp_score_after = score_after + (tmp_score_before - score_before)
            return round((tmp_score_after - tmp_score_before) / abs(tmp_score_before), 2)
        
        else:
            return round((score_after - score_before) / abs(score_before), 2)
        
    def cal_rela_len(self, query_len, ref_len):
        return round((query_len - ref_len) / self.length, 2)
        
    def get_vali_res(self):
        if (not self.analyzed_hap1) or (not self.analyzed_hap2):
            return -1
        
        if self.analyzed_hap1 and self.analyzed_hap2:
            rela_len_1 = self.cal_rela_len(self.len_query_hap1, self.len_ref_hap1)
            rela_len_2 = self.cal_rela_len(self.len_query_hap2, self.len_ref_hap2)
            
            rela_score_1 = self.cal_rela_score(self.score_before_hap1, self.score_after_hap1)
            rela_score_2 = self.cal_rela_score(self.score_before_hap2, self.score_after_hap2)
            
            res_hap1 = self.check_tp(rela_len_1, rela_score_1)
            res_hap2 = self.check_tp(rela_len_2, rela_score_2)
            
            if res_hap1 and res_hap2:
                if abs(rela_len_1 - 1) <= abs(rela_len_2 - 1):
                    return (res_hap1, rela_len_1, rela_score_1)
                else:
                    return (res_hap2, rela_len_2, rela_score_2)
            elif res_hap1:
                return (res_hap1, rela_len_1, rela_score_1)
            elif res_hap2:
                return (res_hap2, rela_len_2, rela_score_2)
            else:
                if abs(rela_len_1 - 1) <= abs(rela_len_2 - 1):
                    return (res_hap1, rela_len_1, rela_score_1)
                else:
                    return (res_hap2, rela_len_2, rela_score_2)
            
        
class alignment:
    def __init__(self, idx, agt_rec, hap, query_length):
        self.idx = idx
        self.ref_name = 'NA'
        self.ref_start = -1
        self.ref_end = -1
        self.contig_name = agt_rec.reference_name
        self.contig_start = agt_rec.reference_start
        self.contig_end = agt_rec.reference_end
        self.query_name = agt_rec.query_name
        #This the index of the first base in seq that is not soft-clipped
        self.query_start = agt_rec.query_alignment_start
        self.query_end = agt_rec.query_alignment_end
        #the index of the last base in seq that is not soft-clipped - the index of the first base in seq that is not soft-clipped
        self.aligned_length = agt_rec.query_alignment_length
        
        #use query length from the fasta file instead!!!
        self.query_length = query_length
        self.hap = hap
    
    def cal_aligned_portion(self):
        return self.aligned_length/self.query_length
    
    def cal_ins_portion(self):
        return 1 - (self.ref_end - self.ref_start)/self.aligned_length
    
    def set_ref_info(self, ref_name, ref_start, ref_end):
        self.ref_name = ref_name
        self.ref_start = ref_start
        self.ref_end = ref_end
        
    def print_info(self):
        print(self.idx, self.ref_name, self.ref_start, self.ref_end, self.contig_name, self.contig_start, self.contig_end, self.query_name, self.query_start, self.query_end,\
             self.aligned_length, self.query_length, self.hap)

#main function
def main():

    #get validation info files
    
    #build map and get validation info on hap1
    contig_name_list_1, contig_pos_list_1, contig_name_dict_1 = get_align_info.build_map_compress(chr_len, interval, liftover_file1, if_hg38)
    #build map and get validation info on hap2
    contig_name_list_2, contig_pos_list_2, contig_name_dict_2 = get_align_info.build_map_compress(chr_len, interval, liftover_file2, if_hg38)
    
    
    #index SVs
    f = pysam.VariantFile(vcf_file,'r')
    sv_list = []
    for count, rec in enumerate(f.fetch()):
        #get sv_type
        try:
            sv_type = rec.info['SVTYPE']
        except:
            print("invalid sv type info")
            continue

        if first_filter(rec, sv_type):
            continue

        #get sv length
        if sv_type == 'INV':
            sv_len = abs(rec.stop - rec.pos + 1)
        else:
            try:
                sv_len = rec.info['SVLEN'][0]
            except:
                try:
                    sv_len = rec.info['SVLEN']
                except:
                    sv_len = abs(rec.stop - rec.pos + 1)
                    #print("invalid sv length info")
    #         try:
    #             sv_len = rec.info['SVLEN'][0]
    #         except:
    #             sv_len = rec.info['SVLEN']
        #handle del length > 0:
        if sv_type == 'DEL':
            sv_len = -abs(sv_len)
            
        if abs(sv_len) < memory_min:
            continue

        #get gt
        #only taking the first sample genotype
    #     if len(rec.samples.values()) != 1:
    #         raise Exception("Wrong number of sample genotype(s)")
    #     gts = [s['GT'] for s in rec.samples.values()]    
        sv_gt = ""
#         for s in rec.samples.values():
#             sv_gt = s['GT']
#             break
        
        sv_list.append(struc_var(count, rec.chrom, sv_type, rec.pos, rec.stop, sv_len, sv_gt))   
        
        #add ins seq for seq-resolved insertion
        #no multi-allelic considered
        if (sv_type == 'INS') and seq_resolved:
            sv_list[len(sv_list)-1].ins_seq = rec.alts[0]
            sv_list[len(sv_list)-1].if_seq_resolved = True
        
    f.close()
    #index sv: second_filter: centromere, non-cov
    #third_filter: size

    for sv in sv_list:
        second_filter(sv)
        third_filter(sv)
    
    get_align_info.get_vali_info(output_dir, vcf_file, query_file1, 1, ref_file, interval, 
              contig_name_list_1, contig_pos_list_1, contig_name_dict_1, memory_limit, if_hg38, chr_list,
              tandem_start_list, tandem_end_list, tandem_info, sv_list, seq_resolved)
    
    get_align_info.get_vali_info(output_dir, vcf_file, query_file2, 2, ref_file, interval, 
              contig_name_list_2, contig_pos_list_2, contig_name_dict_2, memory_limit, if_hg38, chr_list,
              tandem_start_list, tandem_end_list, tandem_info, sv_list, seq_resolved)
    
    #get validation info
    write_vali_info(sv_list)

if __name__ == "__main__":
    main()
    
