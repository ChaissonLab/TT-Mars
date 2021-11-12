##################################################################################
##################################################################################
from Bio import SeqIO
import csv
import math
#pysam: https://pysam.readthedocs.io/en/latest/usage.html
import pysam
#print(pysam.__version__)

from Bio.Seq import Seq
from Bio import pairwise2

import sys 
import numpy as np

import random

import mappy
import os

import random

sys.path.insert(0, '../')
import get_align_info
import get_conf_int
import validate

##################################################################################
##################################################################################
# #hg38 chr length
# chr_len_list = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
#  159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
#  114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 
#  58617616, 64444167, 46709983, 50818468, 156040895]

reg_dup_upper_len = 10000000

output_dir = sys.argv[1] + "/"
if_hg38_input = sys.argv[2]
centromere_file = sys.argv[3]
#exclude_assem1_non_cover_file = sys.argv[4]
#exclude_assem2_non_cover_file = sys.argv[5]
#exclude_high_depth_file = sys.argv[6]
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

liftover_file1_0 = sys.argv[13]
liftover_file2_0 = sys.argv[14]

if_passonly_input = sys.argv[15]

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

#max/min length of allowed SV
memory_limit = 50000
memory_min = 10

#valid types
valid_types = ['DEL', 'INS', 'INV', 'DUP:TANDEM', 'DUP']

#tandem repeats regions file
with open(tandem_file) as f:
    reader = csv.reader(f, delimiter="\t")
    tandem_info = list(reader)
f.close()  

#get tandem start and end list
tandem_start_list, tandem_end_list = get_align_info.get_chr_tandem_shart_end_list(tandem_info, if_hg38)

#query asm files and ref fasta file
query_fasta_file1 = pysam.FastaFile(query_file1)
query_fasta_file2 = pysam.FastaFile(query_file2)
ref_fasta_file = pysam.FastaFile(ref_file)


##################################################################################
##################################################################################
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
    if abs(sv.length) < memory_min or abs(sv.length) > memory_limit:
        sv.is_third_fil = True
        return True
    
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
        
        #for dup validation
        self.valid_non_ins = False
        
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
            if rela_len < 0.675 or rela_len > 1.325:
                result = False
        elif self.sv_type == 'INV':
            if rela_score <= 0:
                result = False
        return result
    
    #TP when wrong length flag presents -- looser rules for TP
    def check_tp_wlen(self, rela_len, rela_score):
        result = True
        if self.sv_type in ['DEL', 'DUP', 'DUP:TANDEM']:
            if rela_score >= 0 and rela_score <= 2.5:
                if rela_len >= -0.05*rela_score + 0.6 and rela_len <= 0.05*rela_score + 1.4:
                    result = True
                else:
                    result = False
            elif rela_score > 2.5:
                if rela_len >= 0.475 and rela_len <= 1.525:
                    result = True
                else:
                    result = False
            else:
                result = False
        elif self.sv_type == 'INS':
            #not seq-resolved
            #if len(self.ins_seq) == 0:
            if not self.if_seq_resolved:
                if rela_len < 0.475 or rela_len > 1.525:
                    result = False
            #seq-resolved
            else:
                if rela_score >= 0 and rela_score <= 2.5:
                    if rela_len >= -0.05*rela_score + 0.6 and rela_len <= 0.05*rela_score + 1.4:
                        result = True
                    else:
                        result = False
                elif rela_score > 2.5:
                    if rela_len >= 0.475 and rela_len <= 1.525:
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
        if self.sv_type == 'DUP':
            #if there's no valid non ins for a DUP call, won't validate it
            if not self.valid_non_ins:
                self.analyzed_hap1 = False
                self.analyzed_hap2 = False
                    
        if (not self.analyzed_hap1) or (not self.analyzed_hap2):
            return -1
        
        if self.analyzed_hap1 and self.analyzed_hap2:
            rela_len_1 = self.cal_rela_len(self.len_query_hap1, self.len_ref_hap1)
            rela_len_2 = self.cal_rela_len(self.len_query_hap2, self.len_ref_hap2)
            
            rela_score_1 = self.cal_rela_score(self.score_before_hap1, self.score_after_hap1)
            rela_score_2 = self.cal_rela_score(self.score_before_hap2, self.score_after_hap2)
            
            if not wrong_len:
                res_hap1 = self.check_tp(rela_len_1, rela_score_1)
                res_hap2 = self.check_tp(rela_len_2, rela_score_2)
            else:
                res_hap1 = self.check_tp_wlen(rela_len_1, rela_score_1)
                res_hap2 = self.check_tp_wlen(rela_len_2, rela_score_2)
            
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
        
#get validation info

def write_vali_info(sv_list):
    g = open(output_dir + "ttmars_regdup_res.txt", "w")
    for sv in sv_list:
        res = sv.get_vali_res()
        
        #skip if not analyzed
        if (not sv.analyzed_hap1) or (not sv.analyzed_hap2):
            continue
        
        g.write(str(sv.ref_name) + "\t")
        g.write(str(sv.sv_pos) + "\t")
        g.write(str(sv.sv_stop) + "\t")
        g.write(str(sv.sv_type) + "\t")
        g.write(str(res[1]) + "\t")
        g.write(str(res[2]) + "\t")
        g.write(str(res[0]))
        g.write("\n")
    g.close()
    
#test validate DUP

def build_map_asm_to_ref(query_fasta_file, liftover_file):
    ref_name_list = []
    ref_pos_list = []
    
    #dict stores: {contig name: [idx, length]}
    contig_idx_len = {}
    ctr = 0
    for contig in query_fasta_file.references:
        contig_len = query_fasta_file.get_reference_length(contig)
        contig_idx_len[contig] = [ctr, contig_len]
        ctr += 1
        ref_name_list.append(np.zeros(int((contig_len+1)//interval) + 1, dtype='int16') - 1)
        ref_pos_list.append(np.zeros(int((contig_len+1)//interval) + 1, dtype='uint32'))
        #print(query_fasta_file.get_reference_length(chrom))

#     #build a dictionary for contig names: to avoid store too many str
#     contig_name_dict = dict()

    with open(liftover_file) as f:
        contig_name_ctr = -1
        pre_contig_name = ""
        for line in f:
            record = line.strip().split()
            int_ref_name = get_align_info.get_int_chr_name(record[0], if_hg38)
            contig_pos = int(record[5])
            ref_pos = int(record[1])

            #store contig name in a dict to save memory
            contig_name = record[4]
                
            contig_idx = contig_idx_len[contig_name][0]
            contig_list_pos = int(contig_pos//interval)

            ref_name_list[contig_idx][contig_list_pos] = int_ref_name
            ref_pos_list[contig_idx][contig_list_pos] = ref_pos
    f.close()
    return ref_name_list, ref_pos_list, contig_idx_len

def build_map_asm_to_ref_compress(query_fasta_file, liftover_file):
    ref_name_list = []
    ref_pos_list = []
    
    #dict stores: {contig name: [idx, length]}
    contig_idx_len = {}
    ctr = 0
    for contig in query_fasta_file.references:
        contig_len = query_fasta_file.get_reference_length(contig)
        contig_idx_len[contig] = [ctr, contig_len]
        ctr += 1
        ref_name_list.append(np.zeros(int((contig_len+1)//interval) + 1, dtype='int16') - 1)
        ref_pos_list.append(np.zeros(int((contig_len+1)//interval) + 1, dtype='uint32'))
        #print(query_fasta_file.get_reference_length(chrom))

    #build a dictionary for contig names: to avoid store too many str
    with open(liftover_file) as f:
        contig_name_ctr = -1
        pre_contig_name = ""
        for line in f:
            record = line.strip().split()
            int_ref_name = get_align_info.get_int_chr_name(record[1], if_hg38)
            contig_name = record[0]
            contig_idx = contig_idx_len[contig_name][0]

            lo_info = record[2].strip().split(";")
            for info in lo_info[:len(lo_info)-1]:
                info_list = info.strip().split(":")
                ctr = int(info_list[2])
                strand = info_list[3]
                if strand == "+":
                    forward = True
                elif strand == "-":
                    forward = False
                for i in range(0, ctr):
                    contig_pos = int(info_list[0]) + i*interval
                    if forward:
                        ref_pos = int(info_list[1]) + i*interval
                    else:
                        ref_pos = int(info_list[1]) - i*interval
                    #pos in the chr list
                    contig_list_pos = int(contig_pos//interval)

                    ref_name_list[contig_idx][contig_list_pos] = int_ref_name
                    ref_pos_list[contig_idx][contig_list_pos] = ref_pos       
    f.close()
    return ref_name_list, ref_pos_list, contig_idx_len


#check overlap of two ref interval
def check_ol(list1, list2):
    #return False: no ol, True: does ol
    #list: [ref_name, start, end]
    
    max_valid_ol_len = 0.5 * min(abs(list1[2] - list1[1]), abs(list2[2] - list2[1]))
    
    if list1[0] != list2[0] or list1[1] > list2[2] - max_valid_ol_len  or list1[2] < list2[1] + max_valid_ol_len:
        return False
    return True

#     ol_start = max(list1[1], list2[1])
#     ol_end = min(list1[2], list2[2])

#check duplicated alignment
def check_duplicate(dup_list, cur_dup_info):
    #cur_dup_info: hap, query_name, ref_name, start, end
    #dup_list[i][1:6]
    for dup in dup_list:
        dup_info = dup[1:7]
        
        first_hap = dup_info[0]
        second_hap = cur_dup_info[0]

        if first_hap == 1:
            contig_idx_len_first = contig_idx_len_1
            ref_name_list_first = ref_name_list_1
            ref_pos_list_first = ref_pos_list_1
        else:
            contig_idx_len_first = contig_idx_len_2
            ref_name_list_first = ref_name_list_2
            ref_pos_list_first = ref_pos_list_2

        if second_hap == 1:
            contig_idx_len_second = contig_idx_len_1
            ref_name_list_second = ref_name_list_1
            ref_pos_list_second = ref_pos_list_1
        else:
            contig_idx_len_second = contig_idx_len_2
            ref_name_list_second = ref_name_list_2
            ref_pos_list_second = ref_pos_list_2

        first_contig_idx = contig_idx_len_first[dup_info[2]][0]
        second_contig_idx = contig_idx_len_second[cur_dup_info[2]][0]

        first_start_ref_name = ref_name_list_first[first_contig_idx][int(dup_info[3]//interval)]
        first_end_ref_name = ref_name_list_first[first_contig_idx][int(dup_info[4]//interval)]

        second_start_ref_name = ref_name_list_second[second_contig_idx][int(cur_dup_info[3]//interval)]
        second_end_ref_name = ref_name_list_second[second_contig_idx][int(cur_dup_info[4]//interval)]

#         if first_start_ref_name != first_end_ref_name or second_start_ref_name != second_end_ref_name:
#             continue
#         if -1 in [first_start_ref_name, first_end_ref_name, second_start_ref_name, second_end_ref_name]:
#             continue

        first_start = ref_pos_list_first[first_contig_idx][int(dup_info[3]//interval)]
        first_end = ref_pos_list_first[first_contig_idx][int(dup_info[4]//interval)]

        second_start = ref_pos_list_second[second_contig_idx][int(cur_dup_info[3]//interval)]
        second_end = ref_pos_list_second[second_contig_idx][int(cur_dup_info[4]//interval)]

        if first_end < first_start:
            tmp = first_start
            first_start = first_end
            first_end = tmp

        if second_end < second_start:
            tmp = second_start
            second_start = second_end
            second_end = tmp
            
        first_aligned_len = dup_info[5]
        second_aligned_len = cur_dup_info[5]
        
        if (first_end - first_start)/first_aligned_len < (1 - valid_ins_ratio) or (second_end - second_start)/second_aligned_len < (1 - valid_ins_ratio):
            continue

        if check_ol([first_start_ref_name, first_start, first_end], [second_start_ref_name, second_start, second_end]):
            return True

#         print([first_start_ref_name, first_start, first_end],  [second_start_ref_name, second_start, second_end])
#         print(potentail_dup[0][1], potentail_dup[1][1], first_end - first_start, second_end - second_start, potentail_dup[0][5], potentail_dup[1][5], potentail_dup[0][6], potentail_dup[1][6])

    return False
    
#get seq function

def getSeqRec(fasta_file, seq_name):
    seq = fasta_file.fetch(seq_name)
    return seq



##################################################################################
##################################################################################
#build lists for excluded SV positions

#Output regions on ref where its not covered by at least one of the assembly
# get_conf_int.get_non_cover_regions(output_dir, bam_file1, 1, chr_list)
# get_conf_int.get_non_cover_regions(output_dir, bam_file2, 2, chr_list)

#Get regions where read depth > 2 * avg_read_depth
#get_conf_int.get_high_depth_calls_info(output_dir, read_bam_file, vcf_file, avg_read_depth)

#Output sv positions
get_conf_int.get_sv_positions(output_dir, vcf_file)

##################################################################################
##################################################################################
#Output filtered calls in non-covered regions
SV_positions_file = output_dir + "SV_positions.bed"
# assem1_non_cov_regions_file = output_dir + "assem1_non_cov_regions.bed"
# assem2_non_cov_regions_file = output_dir + "assem2_non_cov_regions.bed"
get_conf_int.output_non_cov_call_info(output_dir, SV_positions_file, assem1_non_cov_regions_file, assem2_non_cov_regions_file)

#get filtered sv info, using results from get_conf_int.py 
exclude_assem1_non_cover, exclude_assem2_non_cover = validate.get_filtered_sv_pos(output_dir + "exclude_assem1_non_cover.bed", 
                                                                                  output_dir + "exclude_assem2_non_cover.bed")

dict_centromere = validate.build_centro_dict(centromere_file)

##################################################################################
##################################################################################
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
    for s in rec.samples.values():
        sv_gt = s['GT']
        break
    
    sv_list.append(struc_var(count, rec.chrom, sv_type, rec.pos, rec.stop, sv_len, sv_gt))
f.close()
#index sv: second_filter: centromere, non-cov
#third_filter: size

for sv in sv_list:
    second_filter(sv)
    third_filter(sv)

##################################################################################
##################################################################################
#test validate DUP

g = open(output_dir+"all_reg_dup.fasta", "w")
ref_name = ""
ref_rec = ""
for test_sv in sv_list:
    if test_sv.is_sec_fil:
        continue
    
    if test_sv.length > reg_dup_upper_len:
        continue
    
    #test
    if test_sv.sv_type != 'DUP':
        continue

    #test
    if test_sv.idx % 100 == 0:
        print(test_sv.idx)
    
    #test
    #test_sv = sv_list[10]
    
    if test_sv.ref_name != ref_name:
        ref_name = test_sv.ref_name
        ref_rec = ref_fasta_file.fetch(test_sv.ref_name)
    ref_frag = ref_rec[test_sv.sv_pos:test_sv.sv_stop+1]
    ref_frag = ref_frag.upper()
    
    g.write('>' + str(test_sv.idx) + "\n")
    g.write(str(ref_frag) + "\n")
    
g.close()

##################################################################################
##################################################################################
class mappy_alignment:
    def __init__(self, ctr, agt_rec, hap, qname, qlen):
        self.idx = ctr
        self.ref_name = 'NA'
        self.ref_start = -1
        self.ref_end = -1
        self.contig_name = agt_rec.ctg
        self.contig_start = agt_rec.r_st
        self.contig_end = agt_rec.r_en
        self.query_name = qname
        #This the index of the first base in seq that is not soft-clipped
        self.query_start = agt_rec.q_st
        self.query_end = agt_rec.q_en
        #the index of the last base in seq that is not soft-clipped - the index of the first base in seq that is not soft-clipped
        self.aligned_length = agt_rec.mlen
        
        #use query length from the fasta file instead!!!
        self.query_length = qlen
        self.hap = hap
        
        self.contig_int_start = max(int(agt_rec.r_st - 0.75*qlen), 1)
        self.contig_int_end = min(int(agt_rec.r_en + 0.75*qlen), agt_rec.ctg_len)
        
        self.ref_int_start = -1
        self.ref_int_end = -1
    
    def cal_aligned_portion(self):
        return self.aligned_length/self.query_length
    
    def cal_ins_portion(self):
        return 1 - (self.ref_end - self.ref_start)/self.aligned_length
    
    def set_ref_int_info(self, ref_int_start, ref_int_end):
        self.ref_int_start = ref_int_start
        self.ref_int_end = ref_int_end
    
    def cal_ins_rela_len(self):
        return ((self.contig_int_end - self.contig_int_start) - (self.ref_int_end - self.ref_int_start))/self.query_length
    
    def set_ref_info(self, ref_name, ref_start, ref_end):
        self.ref_name = ref_name
        self.ref_start = ref_start
        self.ref_end = ref_end
        
    def print_info(self):
        print(self.idx, self.ref_name, self.ref_start, self.ref_end, self.contig_name, self.contig_start, self.contig_end, self.query_name, self.query_start, self.query_end,\
             self.aligned_length, self.query_length, "hap:", self.hap, self.contig_int_start, self.contig_int_end, self.ref_int_start, self.ref_int_end, 
             self.cal_ins_rela_len()) 

##################################################################################
##################################################################################
#index alignment

aligner_hap1 = mappy.Aligner(fn_idx_in=query_file1)
aligner_hap2 = mappy.Aligner(fn_idx_in=query_file2)

##################################################################################
##################################################################################

ref_name_list_1, ref_pos_list_1, contig_idx_len_1 = build_map_asm_to_ref_compress(query_fasta_file1, liftover_file1_0)

#get asm to ref mapping
ref_name_list_2, ref_pos_list_2, contig_idx_len_2 = build_map_asm_to_ref_compress(query_fasta_file2, liftover_file2_0)

##################################################################################
##################################################################################
alignment_list = []

dup_alm_fasta_file_name = output_dir+"all_reg_dup.fasta"
try:
    dup_alm_fasta_file = pysam.FastaFile(dup_alm_fasta_file_name)
except:
    print("failed open interspersed duplication sequences, may not exist")
    g = open(output_dir + "ttmars_regdup_res.txt", "w")
    g.close()
    
ctr = 0

for seq_name in dup_alm_fasta_file.references:
    dup_seq = getSeqRec(dup_alm_fasta_file, seq_name)
    #if not aligner: raise Exception("ERROR: failed to load/build index")
    dup_aligner_hap1 = aligner_hap1.map(dup_seq, seq2=None, cs=False, MD=False)
    dup_aligner_hap2 = aligner_hap2.map(dup_seq, seq2=None, cs=False, MD=False)
    
    alignment_list.append([])
    for agt in dup_aligner_hap1:
        alignment_list[len(alignment_list)-1].append(mappy_alignment(ctr, agt, 1, seq_name, len(dup_seq)))
        ctr += 1
    for agt in dup_aligner_hap2:
        alignment_list[len(alignment_list)-1].append(mappy_alignment(ctr, agt, 2, seq_name, len(dup_seq)))
        ctr += 1

##################################################################################
##################################################################################
def get_ref_info(alm):
    if alm.hap == 1:
        contig_idx_len = contig_idx_len_1
        ref_name_list = ref_name_list_1
        ref_pos_list = ref_pos_list_1
    else:
        contig_idx_len = contig_idx_len_2
        ref_name_list = ref_name_list_2
        ref_pos_list = ref_pos_list_2

    contig_idx = contig_idx_len[alm.contig_name][0]

    start_ref_name = ref_name_list[contig_idx][int(alm.contig_start//interval)]
    end_ref_name = ref_name_list[contig_idx][int(alm.contig_end//interval)]

    #TODO: find the best interval!

    if start_ref_name != end_ref_name:
        return 'NA', -1, -1
    if -1 in [start_ref_name, end_ref_name]:
        return 'NA', -1, -1

    start = ref_pos_list[contig_idx][int(alm.contig_start//interval)]
    end = ref_pos_list[contig_idx][int(alm.contig_end//interval)]

    if end < start:
        tmp = start
        start = end
        end = tmp  
    
    return start_ref_name, start, end

def find_best_ref_int(alm, contig_idx_len, ref_name_list, ref_pos_list):
    contig_int_start = alm.contig_int_start
    contig_int_end = alm.contig_int_end
    
    best_rela_len = 100
    best_contig_int_start = contig_int_start
    best_contig_int_end = contig_int_end
    best_ref_int_start = -1
    best_ref_int_end = -1
    best_ref_name = 'NA'
    
    for i in range(0, (contig_int_end - contig_int_start)//interval * interval, interval):
        cur_contig_int_start = contig_int_start + i
        cur_contig_int_end = contig_int_end - i
        if cur_contig_int_end <= cur_contig_int_start:
            break
        contig_idx = contig_idx_len[alm.contig_name][0]

        cur_start_ref_int_name = ref_name_list[contig_idx][int(cur_contig_int_start//interval)]
        cur_end_ref_int_name = ref_name_list[contig_idx][int(cur_contig_int_end//interval)]

        #TODO: find the best interval!

        if cur_start_ref_int_name != cur_end_ref_int_name:
            continue
        if -1 in [cur_start_ref_int_name, cur_end_ref_int_name]:
            continue

        cur_ref_int_start = ref_pos_list[contig_idx][int(cur_contig_int_start//interval)]
        cur_ref_int_end = ref_pos_list[contig_idx][int(cur_contig_int_end//interval)]

        if cur_ref_int_end < cur_ref_int_start:
            tmp = cur_ref_int_start
            cur_ref_int_start = cur_ref_int_end
            cur_ref_int_end = tmp
        
        cur_rela_len = ((cur_contig_int_end - cur_contig_int_start) - (cur_ref_int_end - cur_ref_int_start))/alm.query_length
        if abs(cur_rela_len-1) < abs(best_rela_len-1):
            best_rela_len = cur_rela_len
            best_ref_int_start = cur_ref_int_start
            best_ref_int_end = cur_ref_int_end
            best_ref_name = cur_start_ref_int_name
            alm.contig_int_start = cur_contig_int_start
            alm.contig_int_end = cur_contig_int_end
    
    return best_ref_name, best_ref_int_start, best_ref_int_end

def get_ref_int_info(alm):
    if alm.hap == 1:
        contig_idx_len = contig_idx_len_1
        ref_name_list = ref_name_list_1
        ref_pos_list = ref_pos_list_1
    else:
        contig_idx_len = contig_idx_len_2
        ref_name_list = ref_name_list_2
        ref_pos_list = ref_pos_list_2
    
    #test
#     print(alm.hap, alm.contig_int_start, alm.contig_int_end)
    
    start_ref_int_name, start, end = find_best_ref_int(alm, contig_idx_len, ref_name_list, ref_pos_list)

#     contig_idx = contig_idx_len[alm.contig_name][0]

#     start_ref_int_name = ref_name_list[contig_idx][int(alm.contig_int_start//interval)]
#     end_ref_int_name = ref_name_list[contig_idx][int(alm.contig_int_end//interval)]

#     #TODO: find the best interval!

#     if start_ref_int_name != end_ref_int_name:
#         return 'NA', -1, -1
#     if -1 in [start_ref_int_name, end_ref_int_name]:
#         return 'NA', -1, -1

#     start = ref_pos_list[contig_idx][int(alm.contig_int_start//interval)]
#     end = ref_pos_list[contig_idx][int(alm.contig_int_end//interval)]

#     if end < start:
#         tmp = start
#         start = end
#         end = tmp  
    
    return start_ref_int_name, start, end

valid_ins_ratio  = 0.6
valid_aligned_portion = 0.9
ins_rela_len_lb = 0.7
ins_rela_len_ub = 1.3

non_ins_rela_len_ub = 0.4

##################################################################################
##################################################################################
def fake_tp_sv(sv, hap):
    if hap == 1:
        sv.analyzed_hap1 = True
        sv.len_query_hap1 = 1 + sv.length
        sv.len_ref_hap1 = 1
        sv.score_before_hap1 = 1
        sv.score_after_hap1 = 2
    elif hap == 2:
        sv.analyzed_hap2 = True
        sv.len_query_hap2 = 1 + sv.length
        sv.len_ref_hap2 = 1
        sv.score_before_hap2 = 1
        sv.score_after_hap2 = 2
        
def get_sv_list_idx(sv_list, sv_idx):
    sv_list_idx = -1
    for ctr,sv in enumerate(sv_list):
        if sv.idx == sv_idx:
            sv_list_idx = ctr
            break
    return sv_list_idx

##################################################################################
##################################################################################
for agt_list in alignment_list:
    if len(agt_list) < 2:
        continue
    
    for agt in agt_list:
        
        ref_name, start, end = get_ref_info(agt)

    #     if check_ol([first_start_ref_name, first_start, first_end], [second_start_ref_name, second_start, second_end]):
    #         continue
        agt.set_ref_info(ref_name, start, end)
        
        ref_name, int_start, int_end = get_ref_int_info(agt)
        
        agt.set_ref_int_info(int_start, int_end)
        
#         agt.print_info()
        
        if agt.cal_aligned_portion() < valid_aligned_portion:
            continue
        
        if int_start == -1 or int_end == -1:
            sv_idx = int(agt.query_name)
            
            sv_list_idx = get_sv_list_idx(sv_list, sv_idx)
            
            sv_list[sv_list_idx].analyzed_hap1 = False
            sv_list[sv_list_idx].analyzed_hap2 = False
            break
        
        if int_start != -1 and int_end != -1:
            sv_idx = int(agt.query_name)
            
            sv_list_idx = get_sv_list_idx(sv_list, sv_idx)
            
            hap = agt.hap
            sv_list[sv_list_idx].analyzed_hap1 = True
            sv_list[sv_list_idx].analyzed_hap2 = True
            
            if agt.cal_ins_rela_len() > ins_rela_len_lb and agt.cal_ins_rela_len() < ins_rela_len_ub:
#                 print(sv_list_idx)
                fake_tp_sv(sv_list[sv_list_idx], hap)
            elif agt.cal_ins_rela_len() < non_ins_rela_len_ub:
                sv_list[sv_list_idx].valid_non_ins = True   
#             if agt.cal_ins_portion() > valid_ins_ratio:
#                 fake_tp_sv(sv_list[sv_idx], hap)
#                 print(sv_idx)

##################################################################################
##################################################################################
write_vali_info(sv_list)

##################################################################################
##################################################################################
#not analyzed:
#1. no valid non-ins
#2. <2 alignment
#3. at least one alignment failed to be located on ref
##################################################################################
##################################################################################
