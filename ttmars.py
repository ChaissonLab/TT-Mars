##########################################################
##########################################################
#arguments

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("output_dir",
                    help="output directory")
parser.add_argument("centromere_file",
                    help="centromere file")
parser.add_argument("assem1_non_cov_regions_file",
                    help="Regions that are not covered on hap1")
parser.add_argument("assem2_non_cov_regions_file",
                    help="Regions that are not covered on hap2")
parser.add_argument("vcf_file",
                    help="input vcf file")
parser.add_argument("ref_file",
                    help="reference file")
parser.add_argument("query_file1",
                    help="assembly fasta file hap1")
parser.add_argument("query_file2",
                    help="assembly fasta file hap2")
parser.add_argument("liftover_file1",
                    help="liftover file hap1")
parser.add_argument("liftover_file2",
                    help="liftover file hap2")

#needed in interspersed dup validation
parser.add_argument("liftover_file1_0",
                    help="liftover file hap1 asm to ref")
parser.add_argument("liftover_file2_0",
                    help="liftover file hap2 asm to ref")

parser.add_argument("tandem_file",
                    help="tandem repeats regions")

# parser.add_argument("if_hg38_input",
#                     help="if reference is hg38 or not")
parser.add_argument("-n",
                    "--not_hg38",
                    help="if reference is NOT hg38 (hg19)",
                    action="store_true")
# parser.add_argument("if_passonly_input",
#                     help="if consider PASS calls only or not")
parser.add_argument("-p",
                    "--passonly",
                    help="if consider PASS calls only",
                    action="store_true")
# parser.add_argument("seq_resolved_input",
#                     help="if consider sequence resolved calls (INS) or not")
parser.add_argument("-s",
                    "--seq_resolved",
                    help="if consider sequence resolved calls (INS)",
                    action="store_true")
# parser.add_argument("wrong_len_input",
#                     help="if count wrong length calls as True")
parser.add_argument("-w",
                    "--wrong_len",
                    help="if count wrong length calls as True",
                    action="store_true")
parser.add_argument("-g",
                    "--gt_vali",
                    help="conduct genotype validation",
                    action="store_true")



args = parser.parse_args()

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
sys.path.insert(0, '../')

##########################################################
##########################################################

#add reg_dup.py code and test

##########################################################
##########################################################


output_dir = args.output_dir + "/"
# if_hg38_input = args.if_hg38_input
centromere_file = args.centromere_file
#assembly bam files
assem1_non_cov_regions_file = args.assem1_non_cov_regions_file
assem2_non_cov_regions_file = args.assem2_non_cov_regions_file
#avg_read_depth = sys.argv[6]
#read_bam_file = sys.argv[6]
vcf_file = args.vcf_file
#ref fasta file
ref_file = args.ref_file
#assembly fasta files
query_file1 = args.query_file1
query_file2 = args.query_file2
liftover_file1 = args.liftover_file1
liftover_file2 = args.liftover_file2
tandem_file = args.tandem_file
# if_passonly_input = args.if_passonly_input
# seq_resolved_input = args.seq_resolved_input
# wrong_len_input = args.wrong_len_input
liftover_file1_0 = args.liftover_file1_0
liftover_file2_0 = args.liftover_file2_0

##########################################################
##########################################################
#constants

#liftover interval
interval = 20
if_hg38 = not args.not_hg38
# if if_hg38_input == "True":
#     if_hg38 = True
#if pass_only
if_pass_only = args.passonly
# if if_passonly_input == "True":
#     if_pass_only = True
#if seq_resolved
seq_resolved = args.seq_resolved
# if seq_resolved_input == "True":
#     seq_resolved = True
#if include wrong length as TP
wrong_len = args.wrong_len
# if wrong_len_input == "True":
#     wrong_len = True
if_gt = args.gt_vali
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


#main function
def main():
    ##################################################################################
    ##################################################################################
    #validate ordinary SVs
    
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
        if args.gt_vali:
            sv_gt = rec.samples[0]["GT"]
            #bad genotype
            if sv_gt not in [(1, 1), (1, 0), (0, 1)]:
                sv_gt = None
        else:
            sv_gt = None
        
    #     if len(rec.samples.values()) != 1:
    #         raise Exception("Wrong number of sample genotype(s)")
    #     gts = [s['GT'] for s in rec.samples.values()] 
        
        sv_list.append(func.struc_var(count, rec.chrom, sv_type, rec.pos, rec.stop, sv_len, sv_gt, wrong_len))   
        
        #add ins seq for seq-resolved insertion
        #no multi-allelic considered
        if (sv_type == 'INS') and seq_resolved:
            sv_list[len(sv_list)-1].ins_seq = rec.alts[0]
            sv_list[len(sv_list)-1].if_seq_resolved = True
        
    f.close()
    #index sv: second_filter: centromere, non-cov
    #third_filter: size

    for sv in sv_list:
        func.second_filter(sv, if_hg38, dict_centromere, exclude_assem1_non_cover, exclude_assem2_non_cover)
        func.third_filter(sv, memory_min, memory_limit, dup_memory_min, dup_memory_limit)
    
    get_align_info.get_vali_info(output_dir, vcf_file, query_file1, 1, ref_file, interval, 
              contig_name_list_1, contig_pos_list_1, contig_name_dict_1, memory_limit, if_hg38, chr_list,
              tandem_start_list, tandem_end_list, tandem_info, sv_list, seq_resolved)
    
    get_align_info.get_vali_info(output_dir, vcf_file, query_file2, 2, ref_file, interval, 
              contig_name_list_2, contig_pos_list_2, contig_name_dict_2, memory_limit, if_hg38, chr_list,
              tandem_start_list, tandem_end_list, tandem_info, sv_list, seq_resolved)
    
    #get validation info
    func.write_vali_info(sv_list, output_dir, if_gt)
    
    ##################################################################################
    ##################################################################################
    #on chrX, which is a special case
            
    #index chrX SVs
    f = pysam.VariantFile(vcf_file,'r')
    chrx_sv_list = []
    for count, rec in enumerate(f.fetch()):
        #get sv_type
        try:
            sv_type = rec.info['SVTYPE']
        except:
            print("invalid sv type info")
            continue
        
        if not (rec.chrom == 'X' or rec.chrom == 'chrX'):
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
        if args.gt_vali:
            sv_gt = rec.samples[0]["GT"]
            #bad genotype
            if sv_gt not in [(1, 1), (1, 0), (0, 1)]:
                sv_gt = None
        else:
            sv_gt = None
        
    #     if len(rec.samples.values()) != 1:
    #         raise Exception("Wrong number of sample genotype(s)")
    #     gts = [s['GT'] for s in rec.samples.values()] 
        
        chrx_sv_list.append(func.struc_var(count, rec.chrom, sv_type, rec.pos, rec.stop, sv_len, sv_gt, wrong_len))   
        
        #add ins seq for seq-resolved insertion
        #no multi-allelic considered
        if (sv_type == 'INS') and seq_resolved:
            chrx_sv_list[len(chrx_sv_list)-1].ins_seq = rec.alts[0]
            chrx_sv_list[len(chrx_sv_list)-1].if_seq_resolved = True
        
    f.close()
    #index sv: second_filter: centromere, non-cov
    #third_filter: size

    for sv in chrx_sv_list:
        func.second_filter_chrx(sv, if_hg38, dict_centromere, exclude_assem1_non_cover, exclude_assem2_non_cover)
        func.third_filter(sv, memory_min, memory_limit, dup_memory_min, dup_memory_limit)
    
    get_align_info.get_vali_info(output_dir, vcf_file, query_file1, 1, ref_file, interval, 
              contig_name_list_1, contig_pos_list_1, contig_name_dict_1, memory_limit, if_hg38, chr_list,
              tandem_start_list, tandem_end_list, tandem_info, chrx_sv_list, seq_resolved)
    
    get_align_info.get_vali_info(output_dir, vcf_file, query_file2, 2, ref_file, interval, 
              contig_name_list_2, contig_pos_list_2, contig_name_dict_2, memory_limit, if_hg38, chr_list,
              tandem_start_list, tandem_end_list, tandem_info, chrx_sv_list, seq_resolved)            
    
    if len(chrx_sv_list) > 0:
        func.write_vali_info_chrx(chrx_sv_list, output_dir, if_gt)
    
    ##################################################################################
    ##################################################################################
    #validate interspersed DUP
    
    #re-index sv
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
        if args.gt_vali:
            sv_gt = rec.samples[0]["GT"]
            #bad genotype
            if sv_gt not in [(1, 1), (1, 0), (0, 1)]:
                sv_gt = None
        else:
            sv_gt = None
        
    #     if len(rec.samples.values()) != 1:
    #         raise Exception("Wrong number of sample genotype(s)")
    #     gts = [s['GT'] for s in rec.samples.values()] 
        
        sv_list.append(func.struc_var(count, rec.chrom, sv_type, rec.pos, rec.stop, sv_len, sv_gt, wrong_len))   
        
        #add ins seq for seq-resolved insertion
        #no multi-allelic considered
        if (sv_type == 'INS') and seq_resolved:
            sv_list[len(sv_list)-1].ins_seq = rec.alts[0]
            sv_list[len(sv_list)-1].if_seq_resolved = True
        
    f.close()
    #index sv: second_filter: centromere, non-cov
    #third_filter: size

    for sv in sv_list:
        func.second_filter(sv, if_hg38, dict_centromere, exclude_assem1_non_cover, exclude_assem2_non_cover)
        func.third_filter(sv, memory_min, memory_limit, dup_memory_min, dup_memory_limit)

        
    ####################################
    ####################################
    #get duplicated seq
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
    #     if test_sv.idx % 100 == 0:
    #         print(test_sv.idx)

        if test_sv.ref_name != ref_name:
            ref_name = test_sv.ref_name
            ref_rec = ref_fasta_file.fetch(test_sv.ref_name)
        ref_frag = ref_rec[test_sv.sv_pos:test_sv.sv_stop+1]
        ref_frag = ref_frag.upper()

        g.write('>' + str(test_sv.idx) + "\n")
        g.write(str(ref_frag) + "\n")

    g.close()
    
    
    ####################################
    ####################################
    #index alignment

    aligner_hap1 = mappy.Aligner(fn_idx_in=query_file1)
    aligner_hap2 = mappy.Aligner(fn_idx_in=query_file2)

    ####################################
    ####################################

    ref_name_list_1, ref_pos_list_1, contig_idx_len_1 = func.build_map_asm_to_ref_compress(query_fasta_file1, liftover_file1_0, interval, if_hg38)

    #get asm to ref mapping
    ref_name_list_2, ref_pos_list_2, contig_idx_len_2 = func.build_map_asm_to_ref_compress(query_fasta_file2, liftover_file2_0, interval, if_hg38)

    ####################################
    ####################################
    alignment_list = []

    dup_alm_fasta_file_name = output_dir+"all_reg_dup.fasta"
    if_dup = True
    try:
        dup_alm_fasta_file = pysam.FastaFile(dup_alm_fasta_file_name)
    except:
        print("failed open interspersed duplication sequences, may not exist")
        g = open(output_dir + "ttmars_regdup_res.txt", "w")
        g.close()
        if_dup = False

    if if_dup:
        ctr = 0

        for seq_name in dup_alm_fasta_file.references:
            dup_seq = func.getSeqRec(dup_alm_fasta_file, seq_name)
            #if not aligner: raise Exception("ERROR: failed to load/build index")
            dup_aligner_hap1 = aligner_hap1.map(dup_seq, seq2=None, cs=False, MD=False)
            dup_aligner_hap2 = aligner_hap2.map(dup_seq, seq2=None, cs=False, MD=False)

            alignment_list.append([])
            for agt in dup_aligner_hap1:
                alignment_list[len(alignment_list)-1].append(func.mappy_alignment(ctr, agt, 1, seq_name, len(dup_seq)))
                ctr += 1
            for agt in dup_aligner_hap2:
                alignment_list[len(alignment_list)-1].append(func.mappy_alignment(ctr, agt, 2, seq_name, len(dup_seq)))
                ctr += 1


        ###################################
        ###################################
        for agt_list in alignment_list:
            if len(agt_list) < 2:
                continue

            for agt in agt_list:

                ref_name, start, end = func.get_ref_info(agt,
                                                    contig_idx_len_1,
                                                    ref_name_list_1, 
                                                    ref_pos_list_1, 
                                                    contig_idx_len_2, 
                                                    ref_name_list_2,
                                                    ref_pos_list_2,
                                                    interval)

            #     if check_ol([first_start_ref_name, first_start, first_end], [second_start_ref_name, second_start, second_end]):
            #         continue
                agt.set_ref_info(ref_name, start, end)

                ref_name, int_start, int_end = func.get_ref_int_info(agt,
                                                                    contig_idx_len_1,
                                                                    ref_name_list_1, 
                                                                    ref_pos_list_1, 
                                                                    contig_idx_len_2, 
                                                                    ref_name_list_2,
                                                                    ref_pos_list_2,
                                                                    interval)

                agt.set_ref_int_info(int_start, int_end)

        #         agt.print_info()

                if agt.cal_aligned_portion() < valid_aligned_portion:
                    continue

                if int_start == -1 or int_end == -1:
                    sv_idx = int(agt.query_name)

                    sv_list_idx = func.get_sv_list_idx(sv_list, sv_idx)

                    sv_list[sv_list_idx].analyzed_hap1 = False
                    sv_list[sv_list_idx].analyzed_hap2 = False
                    break

                if int_start != -1 and int_end != -1:
                    sv_idx = int(agt.query_name)

                    sv_list_idx = func.get_sv_list_idx(sv_list, sv_idx)

                    hap = agt.hap
                    sv_list[sv_list_idx].analyzed_hap1 = True
                    sv_list[sv_list_idx].analyzed_hap2 = True

                    if agt.cal_ins_rela_len() > ins_rela_len_lb and agt.cal_ins_rela_len() < ins_rela_len_ub:
        #                 print(sv_list_idx)
                        func.fake_tp_sv(sv_list[sv_list_idx], hap)
                    elif agt.cal_ins_rela_len() < non_ins_rela_len_ub:
                        sv_list[sv_list_idx].valid_non_ins = True   
        #             if agt.cal_ins_portion() > valid_ins_ratio:
        #                 fake_tp_sv(sv_list[sv_idx], hap)
        #                 print(sv_idx)

        #################################
        #################################
        func.write_vali_info_reg_dup(sv_list, output_dir, if_gt)
    
    
if __name__ == "__main__":
    main()