#code for simulating and testing DUP

##################################################################################
##################################################################################
#import

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

sys.path.insert(0, '../')
import get_align_info
import get_align_info_rewrite
import get_conf_int
import validate

##################################################################################
##################################################################################
#random.seed(5)
len_lower_bd = 100
len_upper_bd = 100000
no_simu_call = 100

#hg38 chr length
chr_len_list = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 
 58617616, 64444167, 46709983, 50818468, 156040895]

##################################################################################
##################################################################################
file_ctr = sys.argv[1]
output_dir = "./simu_dup/100dup_"+ str(file_ctr) +"/"
altered_query_file1 = output_dir + "altered_query_fasta_file1.fasta"
altered_liftover_file1_0 = output_dir + "lo_pos_assem1_0_result.bed"

if_hg38_input = "True"
centromere_file = "./TT-Mars/centromere_hg38.txt"
#assembly bam files
bam_file1 = "./assem1_nool_sort.bam"
bam_file2 = "./assem2_nool_sort.bam"
vcf_file = "./vcf/file.vcf"
#ref fasta file
ref_file = "./hg38.no_alts.fasta"
#assembly fasta files
query_file1 = "./h1.fa"
query_file2 = "./h2.fa"
liftover_file1 = "./lo_pos_assem1_result.bed"
liftover_file2 = "./lo_pos_assem2_result.bed"
tandem_file = "./TT-Mars/hg38_tandem_repeats.bed"

liftover_file1_0 = "./lo_pos_assem1_0_result.bed"
liftover_file2_0 = "./lo_pos_assem2_0_result.bed"

##################################################################################
##################################################################################
#constants

#liftover interval
interval = 20
if_hg38 = False
if if_hg38_input == "True":
    if_hg38 = True
#chr names
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
tandem_start_list, tandem_end_list = get_align_info_rewrite.get_chr_tandem_shart_end_list(tandem_info, if_hg38)

##################################################################################
##################################################################################
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

##################################################################################
##################################################################################
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

##################################################################################
##################################################################################
#get validation info

def write_vali_info(sv_list):
    g = open(output_dir + "ttmars_res.txt", "w")
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

##################################################################################
##################################################################################
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
            
            #test
#             if contig_name == "cluster18_000019F":
#                 if contig_pos > 2045482 and contig_pos < 2216341:
#                     print(contig_pos)
#             else:
#                 continue

#             if ref_pos > 118697800 and ref_pos < 118698200:
#                 print(record[0], contig_pos, contig_name)
                
                
            contig_idx = contig_idx_len[contig_name][0]
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

##################################################################################
##################################################################################
#get seq function

def getSeqRec(fasta_file, seq_name):
    seq = fasta_file.fetch(seq_name)
    return seq

##################################################################################
##################################################################################
#build map and get validation info on hap1

contig_name_list_1, contig_pos_list_1, contig_name_dict_1 = get_align_info_rewrite.build_map(chr_len, interval, liftover_file1, if_hg38)
#build map and get validation info on hap2

contig_name_list_2, contig_pos_list_2, contig_name_dict_2 = get_align_info_rewrite.build_map(chr_len, interval, liftover_file2, if_hg38)

##################################################################################
##################################################################################
no_of_contigs = len(query_fasta_file1.references)
dup_seq_info = [[] for _ in range(no_of_contigs)]

#note: chr1 only for now

from pysam import VariantFile

vcf_in = VariantFile(vcf_file)
vcfh = vcf_in.header
vcf_out = VariantFile(output_dir + "/simu_dup_1000.vcf", 'w', header=vcfh)

contig_name_dict_rev_1 = {}

#g = open("./HG00514/HG00514_filtered.vcf", "a")
counter = 0
for rec in vcf_in.fetch():
    
    #get random length
    rm_len = random.randrange(len_lower_bd, len_upper_bd, 1)
    
    #get random ref locus
    #chr1 - chrX: 1-23
    rm_ref_chr = random.randrange(1, 23, 1)
    ref_idx = rm_ref_chr - 1
    rm_chr_pos = random.randrange(1, chr_len_list[ref_idx] - rm_len, 1)
    
    ref_chr = "chr" + str(rm_ref_chr)
    
    rec.chrom = ref_chr
    rec.info['SVTYPE'] = 'DUP'
    rec.pos = rm_chr_pos
    rec.stop = rm_chr_pos + rm_len - 1
    rec.info["SVLEN"]=rm_len,
    
    ref_chr_seq = getSeqRec(ref_fasta_file, ref_chr)
    dup_seq = ref_chr_seq[rm_chr_pos:rm_chr_pos + rm_len - 1]
    
    rm_contig_idx = contig_name_list_1[ref_idx][rm_chr_pos//interval]
    rm_contig_pos = contig_pos_list_1[ref_idx][rm_chr_pos//interval]
    if rm_contig_idx == -1 or rm_contig_pos == -1:
        continue
    rm_contig_name = contig_name_dict_1[rm_contig_idx]
    contig_name_dict_rev_1[rm_contig_name] = rm_contig_idx
    #contig name: query_fasta_file1.references[rm_contig_idx]
    rm_contig_len = len(getSeqRec(query_fasta_file1, rm_contig_name))
    
    vcf_out.write(rec)
    
    dup_seq_info[int(rm_contig_idx)].append([rm_contig_idx, rm_contig_name, rm_contig_pos, rm_contig_len, dup_seq])
    
    counter += 1
    if counter >= no_simu_call:
        break
        
vcf_out.close()
vcf_in.close()

##################################################################################
##################################################################################
altered_query_fasta_file1_name = altered_query_file1
#len(chr1_seq)
h = open(altered_query_fasta_file1_name, "w")
for contig in query_fasta_file1.references:
    original_contig_seq = getSeqRec(query_fasta_file1, contig)
    if contig in contig_name_dict_rev_1:
#     try:
        contig_idx = contig_name_dict_rev_1[contig]
    else:
#     except:
        h.write('>' + contig + "\n")
        h.write(original_contig_seq + "\n")
        continue
    if len(dup_seq_info[int(contig_idx)]) > 0:
        temp_seq = original_contig_seq
        for cur_dup_seq_info in dup_seq_info[int(contig_idx)]:
            temp_seq = temp_seq[:cur_dup_seq_info[2]] + cur_dup_seq_info[4] + temp_seq[cur_dup_seq_info[2]:]
        h.write('>' + contig + "\n")
        h.write(temp_seq + "\n")
#     else:
#         h.write('>' + contig + "\n")
#         h.write(original_contig_seq + "\n")
h.close()

##################################################################################
##################################################################################


##################################################################################
##################################################################################


##################################################################################
##################################################################################


##################################################################################
##################################################################################


##################################################################################
##################################################################################


##################################################################################
##################################################################################