#######################################
#######################################


import help_func

import sys
import csv
import pysam
import numpy as np 
import math

import get_conf_int
import validate
import get_align_info

import func

import mappy
import os

#######################################
#######################################


#liftover interval
interval = 20
if_hg38 = True
if_pass_only = False
seq_resolved = True
wrong_len = False

if_gt_info = True
if_phased = True

if_gt = False

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
#max length of allowed interspersed DUP
reg_dup_upper_len = 10000000

#valid types
valid_types = ['DEL', 'INS', 'INV', 'DUP:TANDEM', 'DUP']

#CONST for interspersed DUP
valid_ins_ratio  = 0.6
valid_aligned_portion = 0.9
ins_rela_len_lb = 0.7
ins_rela_len_ub = 1.3
non_ins_rela_len_ub = 0.4

#file and dir
vcf_file = "/scratch2/jianzhiy/elber_vcf/test_vcf/HG00438.minigraph.chr10.vcf"
output_dir = "/scratch2/jianzhiy/ttmars/output/0413_elber_agg_test/elber_minigraph/HG00438_chr10"

#liftover files
liftover_file1 = "/scratch2/jianzhiy/data/assemblies/hprc/HG00438/lra/lo_pos_assem1_result_compressed.bed"
liftover_file2 = "/scratch2/jianzhiy/data/assemblies/hprc/HG00438/lra/lo_pos_assem2_result_compressed.bed"

#ref file
ref_file = "/panfs/qcb-panasas/jianzhiy/data/reference/hg38.no_alts.fasta"

#query file2
query_file_1 = "/scratch2/jianzhiy/data/assemblies/hprc/HG00438/h1.fa"
query_file_2 = "/scratch2/jianzhiy/data/assemblies/hprc/HG00438/h2.fa"

#######################################
#######################################

#build sv dict by idx
def build_sv_idx_dict(sv_list):
    sv_idx_dict = dict()
    for sv in sv_list:
        sv_idx_dict[sv.idx] = sv
    return sv_idx_dict

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

# class sv_group:
#     def __init__(self, sv):
#         self.comb_idx = [sv.idx]
#         self.start = sv.sv_pos
#         self.end = sv.sv_stop
        
#     def add(self, sv):
#         self.comb_idx.append(sv.idx)
#         self.end = sv.sv_stop

def get_sv_groups(sv_list, max_btw_dist):
    sv_groups = []
    cur_group = []
    cur_stop = 0
    for sv in sv_list:
        #allow INS and DEL to aggregate for now
        if sv.sv_type not in ['INS', 'DEL']:
            continue
            
        if sv.is_third_fil:
            continue
        
        if sv.sv_pos < cur_stop + max_btw_dist:
            cur_group.append(sv)
            cur_stop = max(cur_stop, sv.sv_stop)
        else:
            if len(cur_group) > 1:
                sv_groups.append(list(cur_group))
            cur_stop = sv.sv_stop
            cur_group = [sv]
        # test
#         print(sv.sv_pos, cur_stop, sv.sv_stop, len(cur_group), len(sv_groups))
        
    if len(cur_group) > 1:
        sv_groups.append(list(cur_group))
    
    return sv_groups

#for each sv group, get all the valid combinations (of <= 5 SVs)

def get_sv_comb(sv_group, if_gt_info, if_phased):
    max_no_of_sv = 3
    
    #dfs
    cur_idx = 0
    cur_comb = []
    combs = []
    
    if if_gt_info and if_phased:
        #add combs of phase 0
        get_combs_dfs_phased(sv_group, cur_idx, cur_comb, combs, max_no_of_sv, 0)
        
        cur_idx = 0
        cur_comb = []
        
        #add combs of phase 1
        get_combs_dfs_phased(sv_group, cur_idx, cur_comb, combs, max_no_of_sv, 1)
        
    if (not if_gt_info) and (not if_phased):
        get_combs_dfs(sv_group, cur_idx, cur_comb, combs, max_no_of_sv)
    
    return combs

def get_combs_dfs(sv_group, cur_idx, cur_comb, combs, max_no_of_sv):
    if len(cur_comb) == len(sv_group) or len(cur_comb) == max_no_of_sv:
        if check_comb_valid(cur_comb):
            combs.append(list(cur_comb))
        return
    else:
        if len(cur_comb) > 1 and check_comb_valid(cur_comb):
            combs.append(list(cur_comb))
    
    for i in range(cur_idx, len(sv_group)):
        cur_comb.append(sv_group[i])
        get_combs_dfs(sv_group, i+1, cur_comb, combs, max_no_of_sv)
        cur_comb.pop()

def get_combs_dfs_phased(sv_group, cur_idx, cur_comb, combs, max_no_of_sv, phase):
    if len(cur_comb) == len(sv_group) or len(cur_comb) == max_no_of_sv:
        if check_comb_valid(cur_comb):
            combs.append(list(cur_comb))
        return
    else:
        if len(cur_comb) > 1 and check_comb_valid(cur_comb):
            combs.append(list(cur_comb))
    
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
        
    if abs(agg_len) < len_min or abs(agg_len) > len_max:
        #test
#         for sv in comb:
#             print(sv.idx)
        return False

    #check overlapping
    if comb_overlap(comb):
        return False
        
    return True


#######################################
#######################################


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
#main function
def get_agg_align_info(contig_name_list_1, contig_pos_list_1, contig_name_dict_1,
                       contig_name_list_2, contig_pos_list_2, contig_name_dict_2,
                       sv_list, if_gt_info, if_phased, query_file_1, query_file_2,
                       ref_file, interval, if_hg38, if_gt):

    #group SVs
    max_btw_dist = 1000

    sv_idx_dict = build_sv_idx_dict(sv_list)

    sv_groups = get_sv_groups(sv_list, max_btw_dist)

    sv_groups_combs = []
    for sv_group in sv_groups:
        combs = get_sv_comb(sv_group, if_gt_info, if_phased)
        sv_groups_combs.append(combs)

    #######################################
    #######################################
    #preprocessing: liftover mapping

    # contig_name_list_1, contig_pos_list_1, contig_name_dict_1 = get_align_info.build_map_compress(chr_len, interval, liftover_file1, if_hg38)
    # contig_name_list_2, contig_pos_list_2, contig_name_dict_2 = get_align_info.build_map_compress(chr_len, interval, liftover_file2, if_hg38)

    query_fasta_file_1 = pysam.FastaFile(query_file_1)
    query_fasta_file_2 = pysam.FastaFile(query_file_2)

    ref_fasta_file = pysam.FastaFile(ref_file)
    cur_ref_name = ""
    ref_rec_len = 0

    comb_dict = {}

    for combs in sv_groups_combs:
        for comb in combs:
            comb_sv = help_func.idx_comb(comb)
            #test
    #         comb_sv.print_info()

            if cur_ref_name != comb_sv.ref_name:
                cur_ref_name = comb_sv.ref_name
                ref_rec = ref_fasta_file.fetch(cur_ref_name)
                ref_rec_len = len(ref_rec)

            help_func.get_comb_vali_info(comb_sv, 1, interval, contig_name_list_1, contig_pos_list_1, 
                                         contig_name_dict_1, if_hg38, ref_rec, query_fasta_file_1, sv_idx_dict)
            help_func.get_comb_vali_info(comb_sv, 2, interval, contig_name_list_2, contig_pos_list_2, 
                                         contig_name_dict_2, if_hg38, ref_rec, query_fasta_file_2, sv_idx_dict)

            if (not comb_sv.analyzed_hap1) or (not comb_sv.analyzed_hap2):
                #test
                print("NA")
                continue

            #test
            #print res
    #         gt = False
    #         print(comb_sv.get_vali_res(gt))

            comb_dict[tuple(comb_sv.idx)] = comb_sv

    #         help_func.update_sv_res_len_only(comb_sv)

    #update SV info

    for idx_tuple in comb_dict:
        comb_sv = comb_dict[idx_tuple]

        #skip NA comb_sv
        if (not comb_sv.analyzed_hap1) or (not comb_sv.analyzed_hap2):
            continue

        for sv_idx in idx_tuple:
            update_sv_with_comb_res(sv_idx_dict[sv_idx], comb_sv, if_gt)

    #write SV info
    func.write_vali_info_agg(sv_list, output_dir, if_gt)