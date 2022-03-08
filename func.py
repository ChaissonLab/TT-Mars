##################################################################################
##################################################################################
#import

import sys
import csv
import pysam
import numpy as np 
import math

import get_conf_int
import validate
import get_align_info

##################################################################################
##################################################################################
#ttmars.py

########################################
########################################
#define class
class struc_var:
    def __init__(self, idx, ref_name, sv_type, sv_pos, sv_stop, length, gt, wrong_len):
        self.idx = idx
        self.ref_name = ref_name
        self.sv_pos = sv_pos
        self.sv_stop = sv_stop
        self.sv_type = sv_type
        self.length = length
        self.gt = gt
        self.wrong_len = wrong_len
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
        if (not self.analyzed_hap1) or (not self.analyzed_hap2):
            return -1
        
        if self.analyzed_hap1 and self.analyzed_hap2:
            rela_len_1 = self.cal_rela_len(self.len_query_hap1, self.len_ref_hap1)
            rela_len_2 = self.cal_rela_len(self.len_query_hap2, self.len_ref_hap2)
            
            rela_score_1 = self.cal_rela_score(self.score_before_hap1, self.score_after_hap1)
            rela_score_2 = self.cal_rela_score(self.score_before_hap2, self.score_after_hap2)
            
            if not self.wrong_len:
                res_hap1 = self.check_tp(rela_len_1, rela_score_1)
                res_hap2 = self.check_tp(rela_len_2, rela_score_2)
            else:
                res_hap1 = self.check_tp_wlen(rela_len_1, rela_score_1)
                res_hap2 = self.check_tp_wlen(rela_len_2, rela_score_2)
            
            gt_validate = False
            if self.gt:
                if res_hap1 and res_hap2:
                    if self.gt == (1,1):
                        gt_validate = True
                elif res_hap1 or res_hap2:
                    if self.gt == (1,0) or self.gt == (0,1):
                        gt_validate = True
                
            if res_hap1 and res_hap2:
                if abs(rela_len_1 - 1) <= abs(rela_len_2 - 1):
                    return (res_hap1, rela_len_1, rela_score_1, gt_validate)
                else:
                    return (res_hap2, rela_len_2, rela_score_2, gt_validate)
            elif res_hap1:
                return (res_hap1, rela_len_1, rela_score_1, gt_validate)
            elif res_hap2:
                return (res_hap2, rela_len_2, rela_score_2, gt_validate)
            else:
                if abs(rela_len_1 - 1) <= abs(rela_len_2 - 1):
                    return (res_hap1, rela_len_1, rela_score_1, gt_validate)
                else:
                    return (res_hap2, rela_len_2, rela_score_2, gt_validate)
                
    def get_vali_res_reg_dup(self):
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
            
            if not self.wrong_len:
                res_hap1 = self.check_tp(rela_len_1, rela_score_1)
                res_hap2 = self.check_tp(rela_len_2, rela_score_2)
            else:
                res_hap1 = self.check_tp_wlen(rela_len_1, rela_score_1)
                res_hap2 = self.check_tp_wlen(rela_len_2, rela_score_2)
            
            gt_validate = False
            if self.gt:
                if res_hap1 and res_hap2:
                    if self.gt == (1,1):
                        gt_validate = True
                elif res_hap1 or res_hap2:
                    if self.gt == (1,0) or self.gt == (0,1):
                        gt_validate = True
                
            if res_hap1 and res_hap2:
                if abs(rela_len_1 - 1) <= abs(rela_len_2 - 1):
                    return (res_hap1, rela_len_1, rela_score_1, gt_validate)
                else:
                    return (res_hap2, rela_len_2, rela_score_2, gt_validate)
            elif res_hap1:
                return (res_hap1, rela_len_1, rela_score_1, gt_validate)
            elif res_hap2:
                return (res_hap2, rela_len_2, rela_score_2, gt_validate)
            else:
                if abs(rela_len_1 - 1) <= abs(rela_len_2 - 1):
                    return (res_hap1, rela_len_1, rela_score_1, gt_validate)
                else:
                    return (res_hap2, rela_len_2, rela_score_2, gt_validate)  

    def get_vali_res_chrx(self):
        if (not self.analyzed_hap1) and (not self.analyzed_hap2):
            return -1
        
        elif self.analyzed_hap1 and self.analyzed_hap2:
            rela_len_1 = self.cal_rela_len(self.len_query_hap1, self.len_ref_hap1)
            rela_len_2 = self.cal_rela_len(self.len_query_hap2, self.len_ref_hap2)
            
            rela_score_1 = self.cal_rela_score(self.score_before_hap1, self.score_after_hap1)
            rela_score_2 = self.cal_rela_score(self.score_before_hap2, self.score_after_hap2)
            
            if not self.wrong_len:
                res_hap1 = self.check_tp(rela_len_1, rela_score_1)
                res_hap2 = self.check_tp(rela_len_2, rela_score_2)
            else:
                res_hap1 = self.check_tp_wlen(rela_len_1, rela_score_1)
                res_hap2 = self.check_tp_wlen(rela_len_2, rela_score_2)
            
            gt_validate = False
            if self.gt:
                if res_hap1 and res_hap2:
                    if self.gt == (1,1):
                        gt_validate = True
                elif res_hap1 or res_hap2:
                    if self.gt == (1,0) or self.gt == (0,1):
                        gt_validate = True
                
            if res_hap1 and res_hap2:
                if abs(rela_len_1 - 1) <= abs(rela_len_2 - 1):
                    return (res_hap1, rela_len_1, rela_score_1, gt_validate)
                else:
                    return (res_hap2, rela_len_2, rela_score_2, gt_validate)
            elif res_hap1:
                return (res_hap1, rela_len_1, rela_score_1, gt_validate)
            elif res_hap2:
                return (res_hap2, rela_len_2, rela_score_2, gt_validate)
            else:
                if abs(rela_len_1 - 1) <= abs(rela_len_2 - 1):
                    return (res_hap1, rela_len_1, rela_score_1, gt_validate)
                else:
                    return (res_hap2, rela_len_2, rela_score_2, gt_validate)
                
        elif self.analyzed_hap1:
            rela_len_1 = self.cal_rela_len(self.len_query_hap1, self.len_ref_hap1)
            
            rela_score_1 = self.cal_rela_score(self.score_before_hap1, self.score_after_hap1)
            
            res_hap1 = self.check_tp(rela_len_1, rela_score_1)
            
            gt_validate = False
            if self.gt:
                if res_hap1:
                    if self.gt == (1,0) or self.gt == (0,1):
                        gt_validate = True
            
            return (res_hap1, rela_len_1, rela_score_1, gt_validate)
        elif self.analyzed_hap2:
            rela_len_2 = self.cal_rela_len(self.len_query_hap2, self.len_ref_hap2)
            
            rela_score_2 = self.cal_rela_score(self.score_before_hap2, self.score_after_hap2)
            
            res_hap2 = self.check_tp(rela_len_2, rela_score_2)   
            
            gt_validate = False
            if self.gt:
                if res_hap2:
                    if self.gt == (1,0) or self.gt == (0,1):
                        gt_validate = True
            
            return (res_hap2, rela_len_2, rela_score_2, gt_validate)
            
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
        
########################################
########################################
#define functions

#return False if not filtered
#first_filter: type, PASS, chr_name
def first_filter(sv, sv_type, valid_types, if_pass_only, chr_list):
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
def second_filter(sv, if_hg38, dict_centromere, exclude_assem1_non_cover, exclude_assem2_non_cover):
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
    
#second_filter for chrX: centromere, non-cov
def second_filter_chrx(sv, if_hg38, dict_centromere, exclude_assem1_non_cover, exclude_assem2_non_cover):
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
    if validate.check_exclude_chrx(list_to_check, exclude_assem1_non_cover, exclude_assem2_non_cover):
        sv.is_sec_fil = True
        return True
    
#third_filter: size
def third_filter(sv, memory_min, memory_limit, dup_memory_min, dup_memory_limit):
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
def write_vali_info(sv_list, output_dir):
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
        
        if sv.gt:
            g.write("\t" + str(res[3]))
        
        g.write("\n")
    g.close()

    
##################################################################################
##################################################################################
#reg_dup.py

##############################
##############################
#define class

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

        
def build_map_asm_to_ref_compress(query_fasta_file, liftover_file, interval, if_hg38):
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

#check duplicated alignment
def check_duplicate(dup_list, 
                    cur_dup_info,
                    contig_idx_len_1,
                    ref_name_list_1,
                    ref_pos_list_1,
                    contig_idx_len_2,
                    ref_name_list_2,
                    ref_pos_list_2,
                    interval,
                    valid_ins_ratio):
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

def get_ref_info(alm,
                contig_idx_len_1,
                ref_name_list_1, 
                ref_pos_list_1, 
                contig_idx_len_2, 
                ref_name_list_2,
                ref_pos_list_2,
                interval):
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

def find_best_ref_int(alm, contig_idx_len, ref_name_list, ref_pos_list, interval):
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

def get_ref_int_info(alm,
                    contig_idx_len_1,
                    ref_name_list_1, 
                    ref_pos_list_1, 
                    contig_idx_len_2, 
                    ref_name_list_2,
                    ref_pos_list_2,
                    interval):
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
    
    start_ref_int_name, start, end = find_best_ref_int(alm, contig_idx_len, ref_name_list, ref_pos_list, interval)

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
    for ctr, sv in enumerate(sv_list):
        if sv.idx == sv_idx:
            sv_list_idx = ctr
            break
    return sv_list_idx

def write_vali_info_reg_dup(sv_list, output_dir):
    g = open(output_dir + "ttmars_regdup_res.txt", "w")
    for sv in sv_list:
        res = sv.get_vali_res_reg_dup()
        
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
        
        if sv.gt:
            g.write("\t" + str(res[3]))
            
        g.write("\n")
    g.close()
    
def write_vali_info_chrx(sv_list, output_dir):
    g = open(output_dir + "ttmars_chrx_res.txt", "w")
    for sv in sv_list:
        #skip if not analyzed
        if (not sv.analyzed_hap1) and (not sv.analyzed_hap2):
            continue
        
        res = sv.get_vali_res_chrx()
        
        g.write(str(sv.ref_name) + "\t")
        g.write(str(sv.sv_pos) + "\t")
        g.write(str(sv.sv_stop) + "\t")
        g.write(str(sv.sv_type) + "\t")
        g.write(str(res[1]) + "\t")
        g.write(str(res[2]) + "\t")
        g.write(str(res[0]))
                
        if sv.gt:
            g.write("\t" + str(res[3]))
        
        g.write("\n")
    g.close()