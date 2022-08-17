import ttmars

output_dir = ttmars.output_dir
vcf_file = ttmars.vcf_file
ref_file = ttmars.ref_file
#assembly fasta files
query_file1 = ttmars.query_file1
query_file2 = ttmars.query_file2
liftover_file1 = ttmars.liftover_file1
liftover_file2 = ttmars.liftover_file2

#liftover interval
if_hg38 = ttmars.if_hg38
#if pass_only
if_pass_only = ttmars.if_pass_only
#if seq_resolved
seq_resolved = ttmars.seq_resolved
#if include wrong length as TP
wrong_len = ttmars.wrong_len
#if take GT info
if_gt_info = ttmars.if_gt_info
#if consider phased
if_phased = ttmars.if_phased
#if validate GT
if_gt = ttmars.if_gt

# #######################################
# #######################################
import sys
import csv
import pysam
import numpy as np 
import math

import get_conf_int
import validate
import get_align_info

import func
import help_func

import heapq

import mappy
import os

#chr names
chr_list = ttmars.chr_list
#approximate length of chromosomes
chr_len = ttmars.chr_len

#max/min length of allowed SV not DUP
memory_limit = ttmars.memory_limit
memory_min = ttmars.memory_min
#max length of allowed DUP
dup_memory_limit = ttmars.dup_memory_limit
dup_memory_min = ttmars.dup_memory_min
#max length of allowed interspersed DUP
reg_dup_upper_len = ttmars.reg_dup_upper_len
#flanking regions for searching
# region_len_m = 1000
region_len_m = ttmars.region_len_m

#valid types
valid_types = ttmars.valid_types

#CONST for interspersed DUP
valid_ins_ratio  = ttmars.valid_ins_ratio
valid_aligned_portion = ttmars.valid_aligned_portion
ins_rela_len_lb = ttmars.ins_rela_len_lb
ins_rela_len_ub = ttmars.ins_rela_len_ub
non_ins_rela_len_ub = ttmars.non_ins_rela_len_ub

#alt/ref length threshold for abn SV
alt_len_lb = ttmars.alt_len_lb

#group SVs
max_btw_dist = ttmars.max_btw_dist

#max no of sv in a comb
max_no_of_sv = ttmars.max_no_of_sv

#max no of sv in a group
max_sv_group_size = ttmars.max_sv_group_size

#interval length for mapping
interval = ttmars.interval

#######################################
#######################################

#index SVs
def idx_sv(vcf_file):
    f = pysam.VariantFile(vcf_file,'r')
    sv_list = []

    for count, rec in enumerate(f.fetch()):
        #get sv_type
        try:
            sv_type = rec.info['SVTYPE']
        except:
            print("invalid sv type info")
            continue

        if func.first_filter(rec, sv_type, valid_types, if_pass_only, chr_list):
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
        #handle del length > 0:
        if sv_type == 'DEL':
            sv_len = -abs(sv_len)

        if abs(sv_len) < memory_min:
            continue

        #get gt
        #only taking the first sample genotype 
        if if_gt_info:
            sv_gt = rec.samples[0]["GT"]
            #bad genotype
            if sv_gt not in [(1, 1), (1, 0), (0, 1), (None, 1), (1, None)]:
                #test
    #             print("not valid GT", sv_gt, rec.pos, rec.stop)
                sv_gt = None
                continue
        else:
            sv_gt = None

        ref_len = len(rec.ref)
        alt_len = len(rec.alts[0])

        sv_list.append(func.struc_var(count, rec.chrom, sv_type, rec.pos, rec.stop, sv_len, sv_gt, wrong_len, ref_len, alt_len))   

        #add ins seq for seq-resolved insertion
        #no multi-allelic considered
        if (sv_type == 'INS') and seq_resolved:
            sv_list[len(sv_list)-1].ins_seq = rec.alts[0]
            sv_list[len(sv_list)-1].if_seq_resolved = True
            
        #add alt seq for abn DEL
        if (sv_type == 'DEL') and alt_len > alt_len_lb:
            sv_list[len(sv_list)-1].alt_seq = rec.alts[0]

    f.close()

    for sv in sv_list:
        #TODO: use second filter, and get_large_intervals()
    #     func.second_filter(sv, if_hg38, dict_centromere, exclude_assem1_non_cover, exclude_assem2_non_cover)
        func.third_filter(sv, memory_min, memory_limit, dup_memory_min, dup_memory_limit)
    
    return sv_list
    

#######################################
#######################################

#build sv dict by idx
def build_sv_idx_dict(sv_list):
    sv_idx_dict = dict()
    for sv in sv_list:
        sv_idx_dict[sv.idx] = sv
    return sv_idx_dict

#build sv comb dict by idx

def build_sv_comb_idx_dict(sv_groups_combs):
    comb_dict = dict()
    for comb_sv_list in sv_groups_combs:
        for comb_sv in comb_sv_list:
            comb_dict[tuple(comb_sv.idx)] = comb_sv
    return comb_dict

#check if two SV overlapping
def check_sv_ol(sv1, sv2):
    x1 = sv1.sv_pos
    x2 = sv1.sv_stop
    
    y1 = sv2.sv_pos
    y2 = sv2.sv_stop
    
    if x1 <= y2 and y1 <= x2:
        return True
    
    return False

#get sv groups: a group contains SVs that are closed to each other by at most max_btw_dist

#groups allow overlapping

class sv_group_class:
    def __init__(self):
        self.sv_idx_list = []
        
        self.start = -1
        self.end = -1
        
        self.length = 0
        
        #store all combs by idx as list
        self.comb_list = []
        
        #heap of top rela len
        #store the top k comb in terms of rela len
        self.no_top_comb = 5
        #list of tuple: (-dist_2_one, comb_key)
        #can be empty
        self.top_rela_len_comb_idx = []
        
    def add_sv(self, sv):
        self.sv_idx_list.append(sv.idx)
        
        if self.start < 0:
            self.start = sv.sv_pos
        
        self.end = max(self.end, sv.sv_stop)
        self.length += 1
    
    def add_comb(self, comb_sv):
        self.comb_list.append(comb_sv.idx)
    
    def cal_dist_2_one(self, comb_sv):
        rela_len_1 = comb_sv.cal_rela_len(comb_sv.len_query_hap1, comb_sv.len_ref_hap1)
        rela_len_2 = comb_sv.cal_rela_len(comb_sv.len_query_hap2, comb_sv.len_ref_hap2)

        dist_2_one_1 = abs(rela_len_1 - 1)
        dist_2_one_2 = abs(rela_len_2 - 1)

        dist_2_one = min(dist_2_one_1, dist_2_one_2)
        return dist_2_one
        
    def get_top_comb(self, comb_dict):
        for comb_idx_list in self.comb_list:
            comb_key = tuple(comb_idx_list)
            comb_sv = comb_dict[comb_key]
            if (not comb_sv.analyzed_hap1) or (not comb_sv.analyzed_hap2):
                continue
            
            dist_2_one = round(self.cal_dist_2_one(comb_sv), 2)
            
            heapq.heappush(self.top_rela_len_comb_idx, (-dist_2_one, comb_key))
            
            if len(self.top_rela_len_comb_idx) > self.no_top_comb:
                heapq.heappop(self.top_rela_len_comb_idx)
            
def get_sv_groups(sv_list, max_btw_dist, max_sv_group_size):
    sv_groups = []
    cur_group = sv_group_class()
    pre_chr = ""
    
    for sv in sv_list:
        #allow INS and DEL to aggregate for now
        if sv.sv_type not in ['INS', 'DEL']:
            continue
            
        if sv.is_third_fil:
            continue
        
        #if move to the next chr
        cur_chr = sv.ref_name
        if cur_chr != pre_chr:
            if cur_group.length > 1 and cur_group.length <= max_sv_group_size:
                sv_groups.append(cur_group)
            cur_group = sv_group_class()
            cur_group.add_sv(sv)
            pre_chr = cur_chr
            continue
        
        if sv.sv_pos < cur_group.end + max_btw_dist:
            cur_group.add_sv(sv)
        else:
            if cur_group.length > 1 and cur_group.length <= max_sv_group_size:
                sv_groups.append(cur_group)
            cur_group = sv_group_class()
            cur_group.add_sv(sv)
        # test
#         print(sv.sv_pos, cur_group.end, sv.sv_stop, cur_group.length, cur_group.sv_idx_list, len(sv_groups))
        
    if cur_group.length > 1 and cur_group.length <= max_sv_group_size:
        sv_groups.append(cur_group)
    
    return sv_groups


#for each sv group, get all the valid combinations (of <= n SVs)

def get_sv_comb(sv_group, if_gt_info, if_phased, sv_idx_dict, max_no_of_sv):
    sv_group_list = []
    for idx in sv_group.sv_idx_list:
        sv_group_list.append(sv_idx_dict[idx])
    
#     max_no_of_sv = 3
    
    #dfs
    cur_idx = 0
    cur_comb = []
    combs = []
    
    if if_gt_info and if_phased:
        #add combs of phase 0
        get_combs_dfs_phased(sv_group_list, cur_idx, cur_comb, combs, max_no_of_sv, 0)
        
        cur_idx = 0
        cur_comb = []
        
        #add combs of phase 1
        get_combs_dfs_phased(sv_group_list, cur_idx, cur_comb, combs, max_no_of_sv, 1)
        
    if (not if_gt_info) or (not if_phased):
        get_combs_dfs(sv_group_list, cur_idx, cur_comb, combs, max_no_of_sv)
    
    return combs

def get_combs_dfs(sv_group, cur_idx, cur_comb, combs, max_no_of_sv):
    if len(cur_comb) == len(sv_group) or len(cur_comb) == max_no_of_sv:
        if check_comb_valid(cur_comb):
            combs.append(help_func.idx_comb(list(cur_comb)))
#             combs.append(list(cur_comb))
        return
    else:
        if len(cur_comb) > 1 and check_comb_valid(cur_comb):
            combs.append(help_func.idx_comb(list(cur_comb)))
#             combs.append(list(cur_comb))
    
    for i in range(cur_idx, len(sv_group)):
        cur_comb.append(sv_group[i])
        get_combs_dfs(sv_group, i+1, cur_comb, combs, max_no_of_sv)
        cur_comb.pop()

def get_combs_dfs_phased(sv_group, cur_idx, cur_comb, combs, max_no_of_sv, phase):
    if len(cur_comb) == len(sv_group) or len(cur_comb) == max_no_of_sv:
        if check_comb_valid(cur_comb):
            combs.append(help_func.idx_comb(list(cur_comb)))
#             combs.append(list(cur_comb))
        return
    else:
        if len(cur_comb) > 1 and check_comb_valid(cur_comb):
            combs.append(help_func.idx_comb(list(cur_comb)))
#             combs.append(list(cur_comb))
    
    for i in range(cur_idx, len(sv_group)):
        #check phasing
        if sv_group[i].gt[phase] != 1:
            continue
        
        cur_comb.append(sv_group[i])
        get_combs_dfs_phased(sv_group, i+1, cur_comb, combs, max_no_of_sv, phase)
        cur_comb.pop()

def comb_overlap(comb):
    cur_sv = comb[0]
    for next_sv in comb[1:]:
        if check_sv_ol(cur_sv, next_sv):
            return True
        else:
            cur_sv = next_sv

    return False

#check if length valid, if overlapping
def check_comb_valid(comb):
    #check length requirments
    len_min = 50
    len_max = 500000
    
    #validate length
    agg_len = 0
    for sv in comb:
        agg_len += sv.length
        
    #comb sv length
    if abs(agg_len) < len_min or abs(agg_len) > len_max:
        #test
#         for sv in comb:
#             print(sv.idx)
        return False

    #length of sequence that comb sv spans
    if comb[-1].sv_stop - comb[0].sv_pos > len_max:
        return False

    #check overlapping
    if comb_overlap(comb):
        return False
        
    return True

#return combs that include all possible SVs on one hap
#require phased SVs
def get_hap_all_comb(sv_group):
    #gready
    haps = [0, 1]
    combs = []
    
    for hap in haps:
        cur_comb = []
        for sv in sv_group:
            if sv.gt[hap] == 1:
                cur_comb.append(sv)
        if len(cur_comb) > 1:
            if check_comb_valid(cur_comb):
                combs.append(cur_comb)
    
    return combs

#match each info needed to find res
def match_sv_with_comb_res(sv, comb_sv):
    sv.length = comb_sv.length
    
    sv.analyzed_hap1 = comb_sv.analyzed_hap1
    sv.analyzed_hap2 = comb_sv.analyzed_hap2
    sv.len_query_hap1 = comb_sv.len_query_hap1
    sv.len_query_hap2 = comb_sv.len_query_hap2 
    sv.len_ref_hap1 = comb_sv.len_ref_hap1
    sv.len_ref_hap2 = comb_sv.len_ref_hap2 
    sv.score_before_hap1 = comb_sv.score_before_hap1
    sv.score_after_hap1 = comb_sv.score_after_hap1
    sv.score_before_hap2 = comb_sv.score_before_hap2 
    sv.score_after_hap2 = comb_sv.score_after_hap2
    
#check if the comb res is better than an sv
def check_comb_better_than_sv(sv, comb_sv, if_gt):
    #return a tuple: (res, rela_len, rela_score, gt_validate)
    comb_res = comb_sv.get_vali_res(if_gt)
    sv_res = sv.get_vali_res(if_gt)
    
    if comb_res[0] and (not sv_res[0]):
        return True
    elif sv_res[0] and (not comb_res[0]):
        return False
    else:
        comb_rela_len = comb_res[1]
        sv_rela_len = sv_res[1]
        
        if abs(comb_rela_len - 1) <= abs(sv_rela_len - 1):
            return True
        else:
            return False

#update sv info with comb res
def update_sv_with_comb_res(sv, comb_sv, if_gt):
    if (not comb_sv.analyzed_hap1) or (not comb_sv.analyzed_hap2):
        return
        
    if (not sv.analyzed_hap1) or (not sv.analyzed_hap2):
        match_sv_with_comb_res(sv, comb_sv)
    else:
        if check_comb_better_than_sv(sv, comb_sv, if_gt):
            match_sv_with_comb_res(sv, comb_sv)

#######################################
#######################################
    
    
# #preprocessing: liftover mapping

# contig_name_list_1, contig_pos_list_1, contig_name_dict_1 = get_align_info.build_map_compress(chr_len, interval, liftover_file1, if_hg38)
# contig_name_list_2, contig_pos_list_2, contig_name_dict_2 = get_align_info.build_map_compress(chr_len, interval, liftover_file2, if_hg38)


#######################################
#######################################
#analysis data

# #build sv groups
# sv_list = idx_sv(vcf_file)
# sv_idx_dict = build_sv_idx_dict(sv_list)

# sv_groups = get_sv_groups(sv_list, max_btw_dist, max_sv_group_size)

# assert len(sv_groups) > 0

# #build sv combs
# sv_groups_combs = []
# for sv_group in sv_groups:
#     comb_sv_list = get_sv_comb(sv_group, if_gt_info, if_phased, sv_idx_dict, max_no_of_sv)
    
#     for comb_sv in comb_sv_list:
#         sv_group.add_comb(comb_sv)
    
#     sv_groups_combs.append(comb_sv_list)

# #build comb dict based on [idx]
# comb_dict = build_sv_comb_idx_dict(sv_groups_combs)

# query_fasta_file_1 = pysam.FastaFile(query_file1)
# query_fasta_file_2 = pysam.FastaFile(query_file2)
# ref_fasta_file = pysam.FastaFile(ref_file)
# cur_ref_name = ""

# #test
# counter = 0

# for comb_sv_list in sv_groups_combs:
#     for comb_sv in comb_sv_list:
#         #test
#         counter += 1
#         if counter % 500 == 1:
#             print(counter)
        
#         if cur_ref_name != comb_sv.ref_name:
#             cur_ref_name = comb_sv.ref_name
#             ref_rec = ref_fasta_file.fetch(cur_ref_name)
        
#         help_func.get_comb_vali_info_len_only(comb_sv, 1, interval, contig_name_list_1, contig_pos_list_1, 
#                                      contig_name_dict_1, if_hg38, ref_rec, query_fasta_file_1, sv_idx_dict, region_len_m)
        
#         help_func.get_comb_vali_info_len_only(comb_sv, 2, interval, contig_name_list_2, contig_pos_list_2, 
#                                      contig_name_dict_2, if_hg38, ref_rec, query_fasta_file_2, sv_idx_dict, region_len_m)
        
#         if (not comb_sv.analyzed_hap1) or (not comb_sv.analyzed_hap2):
#             continue

# #         help_func.update_sv_res_len_only(comb_sv)


# #for each sv group, find the top k comb in terms for the rela length (for either haps)
# for cur_group in sv_groups:
#     #test
# #     print(cur_group.comb_list)
#     cur_group.get_top_comb(comb_dict)



# cur_ref_name = ""

# #test
# counter = 0

# for cur_group in sv_groups:
#     #can be empty if every comb failed to be analyzed
#     if len(cur_group.top_rela_len_comb_idx) > 0:
            
#         for _, comb_key in cur_group.top_rela_len_comb_idx:
#             #test
#             counter += 1
#             if counter % 100 == 1:
#                 print(counter)
                
#             comb_sv = comb_dict[comb_key]
            
#             if cur_ref_name != comb_sv.ref_name:
#                 cur_ref_name = comb_sv.ref_name
#                 ref_rec = ref_fasta_file.fetch(cur_ref_name)

#             help_func.get_comb_vali_info_align_only(comb_sv, 1, interval, contig_name_list_1, contig_pos_list_1, 
#                                          contig_name_dict_1, if_hg38, ref_rec, query_fasta_file_1, sv_idx_dict)
#             help_func.get_comb_vali_info_align_only(comb_sv, 2, interval, contig_name_list_2, contig_pos_list_2, 
#                                          contig_name_dict_2, if_hg38, ref_rec, query_fasta_file_2, sv_idx_dict)

#             if (not comb_sv.analyzed_hap1) or (not comb_sv.analyzed_hap2):
#                 continue
                

# #update SV info

# for cur_group in sv_groups:
#     #don't consider comb not in the top list
#     #can be empty if every comb failed to be analyzed
#     if len(cur_group.top_rela_len_comb_idx) > 0:
#         for _, comb_key in cur_group.top_rela_len_comb_idx:
#             comb_sv = comb_dict[comb_key]
    
#             #skip NA comb_sv
#             if (not comb_sv.analyzed_hap1) or (not comb_sv.analyzed_hap2):
#                 continue

#             for sv_idx in comb_key:
#                 update_sv_with_comb_res(sv_idx_dict[sv_idx], comb_sv, if_gt)            

# #write SV info
# func.write_vali_info_agg(sv_list, output_dir, if_gt)