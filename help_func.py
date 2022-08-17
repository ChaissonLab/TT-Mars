import pysam
import csv
from Bio import Align

import func
import get_align_info
from Bio.Seq import Seq

#constants
import ttmars

memory_min = ttmars.memory_min

if_hg38 = ttmars.if_hg38
if_pass_only = ttmars.if_hg38

valid_types = ttmars.valid_types

#not consider chrY calls
chr_list = ttmars.chr_list

#modified from ttmars.py to take sv dict object as input
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


#functions

#build sv dict
def build_dict(vcf):
    sv_dict = dict()
    f = pysam.VariantFile(vcf)
    for counter, rec in enumerate(f.fetch()):
        #test
        #if counter < 91:
        #    continue

        ref_name = rec.chrom

        sv_type = rec.info['SVTYPE']
    #     sv_len = rec.rlen
    
        if first_filter(rec, sv_type, valid_types, if_pass_only, chr_list):
            continue

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

        if sv_type == 'DEL':
            sv_len = -abs(sv_len)
            
        if abs(sv_len) < memory_min:
            continue

        #TODOL double check the start for different types
        sv_pos = rec.pos
        sv_end = rec.stop
        
        n = 5
        if_abnormal = False
        if len(rec.ref) > n:
            if sv_type == 'INS':
                if_abnormal = True
        if len(rec.alts[0]) > n:
            if sv_type == 'DEL':
                if_abnormal = True

#         if (ref_name, int(sv_pos), int(sv_end), sv_type, int(sv_len)) in sv_dict:
#             print((ref_name, int(sv_pos), int(sv_end), sv_type, int(sv_len)))
        #dict to store sv: {(chr, start, end, type): [vapor score, vapor gt, vapor result, truvari result, ttmars result, rela_len, rela_score, len, assem_score, wlen result]}
        sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = ['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', int(sv_len), 'NA', 'NA', rec.filter.keys(), if_abnormal]
        #test
        #else:
        #    print(sv_type)
    f.close()  
    return sv_dict



def add_ttmars_res(sv_dict, ttmars_res):
    with open(ttmars_res) as f:
        reader = csv.reader(f, delimiter="\t")
        ttmars_result = list(reader)
    f.close()

    for rec in ttmars_result:
        key = (rec[0], int(rec[1]), int(rec[2]), rec[3])
        try:
            sv_dict[key]
        except:
            continue
        if rec[6] == 'True':
            sv_dict[key][4] = 'tp'
        else:
            sv_dict[key][4] = 'fp'
        sv_dict[key][5] = float(rec[4])
        sv_dict[key][6] = float(rec[5])
        
def tuple_2_list(input_tup):
    output_list = []
    for ele in input_tup:
        output_list.append(ele)
    return output_list


#in ttmars.py:
#call get_align_info.get_vali_info for both haplotypes
#then func.write_vali_info(sv_list, output_dir, if_gt)

#idx comb sv
#requires seq resolved here
def idx_comb(comb):
    ref_name = comb[0].ref_name
    wrong_len = False
    
    sv_len = 0
    ref_len = 0
    alt_len = 0
    
    sv_pos = float('inf')
    sv_stop = 0
    
    sv_idx = []
    
    for sv in comb:
        sv_idx.append(sv.idx)
        
        sv_len += sv.length
        
        ref_len += sv.ref_len
        alt_len += sv.alt_len

        sv_pos = min(sv_pos, sv.sv_pos)
        sv_stop = max(sv_stop, sv.sv_stop)
        
    sv_type = ""
    if sv_len > 0:
        sv_type = "INS"
    else:
        sv_type = "DEL"
        
    sv_gt = None
    
    comb_sv = func.struc_var(sv_idx, ref_name, sv_type, sv_pos, sv_stop, sv_len, sv_gt, wrong_len, ref_len, alt_len)
    comb_sv.if_seq_resolved = True
    
    return comb_sv

# def idx_comb_len_only(comb):
#     ref_name = comb[0].ref_name
#     wrong_len = False
    
#     sv_len = 0
#     ref_len = 0
#     alt_len = 0
    
#     sv_pos = float('inf')
#     sv_stop = 0
    
#     sv_idx = []
    
#     for sv in comb:
#         sv_idx.append(sv.idx)
        
#         sv_len += sv.length
        
#         ref_len += sv.ref_len
#         alt_len += sv.alt_len

#         sv_pos = min(sv_pos, sv.sv_pos)
#         sv_stop = max(sv_stop, sv.sv_stop)
        
#     sv_type = ""
#     if sv_len > 0:
#         sv_type = "INS"
#     else:
#         sv_type = "DEL"
        
#     sv_gt = None
    
#     comb_sv = func.struc_var(sv_idx, ref_name, sv_type, sv_pos, sv_stop, sv_len, sv_gt, wrong_len, ref_len, alt_len)
    
#     return comb_sv

def get_comb_vali_info_align_only(comb_sv, hap, interval, contig_name_list, contig_pos_list, contig_name_dict, if_hg38, ref_rec, query_fasta_file, sv_idx_dict):    
    if hap == 1:
        query_rec = query_fasta_file.fetch(comb_sv.query_name_hap1)
        ref_start = comb_sv.ref_start_best_hap1
        ref_end = comb_sv.ref_end_best_hap1
        query_start = comb_sv.query_start_best_hap1
        query_end = comb_sv.query_end_best_hap1
        neg_strand = comb_sv.neg_strand_hap1
    elif hap == 2:
        query_rec = query_fasta_file.fetch(comb_sv.query_name_hap2)
        ref_start = comb_sv.ref_start_best_hap2
        ref_end = comb_sv.ref_end_best_hap2
        query_start = comb_sv.query_start_best_hap2
        query_end = comb_sv.query_end_best_hap2
        neg_strand = comb_sv.neg_strand_hap2
    
    if query_start >= len(query_rec) or query_end >= len(query_rec):
        message = "bad_query_pos"
        #write_err(output_file_name, message, g)
        if hap == 1:
            comb_sv.analyzed_hap1 = False
        if hap == 2:
            comb_sv.analyzed_hap2 = True
        
        return False    
    
    #definitely need alignment score for comb_sv
    
    query_frag = query_rec[query_start:query_end]
    ref_frag = ref_rec[ref_start:ref_end]
    ref_after_sv_frag = get_comb_ref_frag_after_sv(comb_sv, ref_rec, sv_idx_dict, ref_start, ref_end)
    
    #neg strand
    if neg_strand:
        seq = Seq(query_frag)
        query_frag = seq.reverse_complement()

    #get to upper case
    ref_frag = ref_frag.upper()
    query_frag = query_frag.upper()
    ref_after_sv_frag = ref_after_sv_frag.upper()
    
    #test
#     print(len(ref_frag), len(ref_after_sv_frag), len(ref_after_sv_frag)-len(ref_frag))

    #TODO: find appropriate alignment parameters
    #paras: match, mismatch, open gap, extend gap
    alignment_beforeSV, alignment_afterSV = align_before_after_comb(query_frag, ref_frag, ref_after_sv_frag)    
    
    if hap == 1:
        comb_sv.analyzed_hap1 = True
        comb_sv.score_before_hap1 = alignment_beforeSV
        comb_sv.score_after_hap1 = alignment_afterSV
#         comb_sv.len_query_hap1 = query_end - query_start + 1
#         comb_sv.len_ref_hap1 = ref_end - ref_start + 1
    elif hap == 2:
        comb_sv.analyzed_hap2 = True
        comb_sv.score_before_hap2 = alignment_beforeSV
        comb_sv.score_after_hap2 = alignment_afterSV
#         comb_sv.len_query_hap2 = query_end - query_start + 1
#         comb_sv.len_ref_hap2 = ref_end - ref_start + 1
    
    return True
    
def get_comb_vali_info_len_only(comb_sv, hap, interval, contig_name_list, contig_pos_list, contig_name_dict, if_hg38, ref_rec, query_fasta_file, sv_idx_dict, region_len_m):
    
    if not get_comb_intervals_len_only(comb_sv, if_hg38, interval, contig_name_list, contig_pos_list, contig_name_dict, hap, ref_rec, region_len_m):
        return False
    
    if hap == 1:
#         query_rec = query_fasta_file.fetch(comb_sv.query_name_hap1)
        ref_start = comb_sv.ref_start_best_hap1
        ref_end = comb_sv.ref_end_best_hap1
        query_start = comb_sv.query_start_best_hap1
        query_end = comb_sv.query_end_best_hap1
        neg_strand = comb_sv.neg_strand_hap1
    elif hap == 2:
#         query_rec = query_fasta_file.fetch(comb_sv.query_name_hap2)
        ref_start = comb_sv.ref_start_best_hap2
        ref_end = comb_sv.ref_end_best_hap2
        query_start = comb_sv.query_start_best_hap2
        query_end = comb_sv.query_end_best_hap2
        neg_strand = comb_sv.neg_strand_hap2
    
#     if query_start >= len(query_rec) or query_end >= len(query_rec):
#         message = "bad_query_pos"
#         #write_err(output_file_name, message, g)
#         return False    
    
    #skip alignment in this mode (although alignment score is suggested for comb_sv)
    
#     query_frag = query_rec[query_start:query_end]
#     ref_frag = ref_rec[ref_start:ref_end]
    
    if hap == 1:
        comb_sv.analyzed_hap1 = True
#         comb_sv.score_before_hap1 = alignment_beforeSV
#         comb_sv.score_after_hap1 = alignment_afterSV
        comb_sv.len_query_hap1 = query_end - query_start + 1
        comb_sv.len_ref_hap1 = ref_end - ref_start + 1
    elif hap == 2:
        comb_sv.analyzed_hap2 = True
#         comb_sv.score_before_hap2 = alignment_beforeSV
#         comb_sv.score_after_hap2 = alignment_afterSV
        comb_sv.len_query_hap2 = query_end - query_start + 1
        comb_sv.len_ref_hap2 = ref_end - ref_start + 1
    
    return True
    
def get_comb_intervals_len_only(sv, if_hg38, interval, contig_name_list, contig_pos_list, contig_name_dict, hap, ref_rec, region_len_m):
    ref_rec_len = len(ref_rec)
#     region_len_m = 1000
    
    #test:
#     print(sv.sv_pos, sv.sv_stop)

    ref_start_start = max(get_align_info.getRefStart(sv.sv_pos, interval) - region_len_m, 0)
    ref_start_end = get_align_info.getRefStart(sv.sv_pos, interval)

    ref_end_start = get_align_info.getRefEnd(sv.sv_stop, interval)
    ref_end_end = min(get_align_info.getRefEnd(sv.sv_stop, interval) + region_len_m, get_align_info.getRefEnd(ref_rec_len, interval) - interval) 

    #first level key: chr index as an int
    int_ref_name = get_align_info.get_int_chr_name(sv.ref_name, if_hg38)
    lo_list_index = int_ref_name - 1
    first_key = lo_list_index
    
    #second level key: ref list_pos
    min_diff_len = float('inf')
    ref_start_best = ref_start_end
    ref_end_best = ref_end_start
    query_start_best = None
    query_end_best = None
    start_contig_name_ctr = -1
    end_contig_name_ctr = -1
    neg_strand = False
    
    #test
#     print(ref_end_end, ref_end_start, ref_end_start - ref_end_end)
#     print(ref_start_start, ref_start_end, ref_start_start - ref_start_end)

    for ref_end_cand in range(ref_end_end, ref_end_start, -interval):
        for ref_start_cand in range(ref_start_start, ref_start_end, interval):
            if ref_end_cand <= ref_start_cand:
                continue

            #second_key_start = str(ref_start_cand)
            #second_key_end = str(ref_end_cand)
            second_key_start = ref_start_cand//interval
            second_key_end = ref_end_cand//interval

            start_contig_name_ctr_cand = contig_name_list[first_key][second_key_start]
            end_contig_name_ctr_cand = contig_name_list[first_key][second_key_end]

            if start_contig_name_ctr_cand == -1 or end_contig_name_ctr_cand == -1:
                #print("wrong second key")
                message = "wrong_sec_key"
                #write_err(output_file_name, message, g)
                continue

            if start_contig_name_ctr_cand != end_contig_name_ctr_cand:
                #print("Not same contig")
                message = "not_same_contig"
                #write_err(output_file_name, message, g)
                continue

            query_start = contig_pos_list[first_key][second_key_start]
            query_end = contig_pos_list[first_key][second_key_end]

            neg_strand_tep = False
            #in case: negtive strand
            if query_end < query_start:
                tem = query_start
                query_start = query_end
                query_end = tem
                neg_strand_tep = True

            #take the best relative length to be the optimal interval
            if abs((query_end - query_start) - (ref_end_cand - ref_start_cand) - sv.length) < min_diff_len:
                min_diff_len = abs((query_end - query_start) - (ref_end_cand - ref_start_cand) - sv.length)
                ref_start_best = ref_start_cand
                ref_end_best = ref_end_cand
                query_start_best = query_start
                query_end_best = query_end
                start_contig_name_ctr = start_contig_name_ctr_cand
                end_contig_name_ctr = end_contig_name_ctr_cand

                if neg_strand_tep:
                    neg_strand = True
                else:
                    neg_strand = False

    if query_start_best == None or query_end_best == None:
        #print("Wrong query pos")
        message = "Wrong_query_pos"
        #write_err(output_file_name, message, g)
        return False

    if ref_start_best == sv.sv_pos:
        ref_start_best = ref_start_best - 1
    if ref_end_best == sv.sv_stop:
        ref_end_best = ref_end_best + 1
    if query_end_best == query_start_best:
        query_end_best += 1
    
    #test
#     print(ref_start_best, ref_end_best)
#     print(query_start_best, query_end_best)

    if hap == 1:
        sv.ref_start_best_hap1 = ref_start_best
        sv.ref_end_best_hap1 = ref_end_best
        sv.query_start_best_hap1 = query_start_best
        sv.query_end_best_hap1 = query_end_best
        sv.query_name_hap1 = contig_name_dict[start_contig_name_ctr]
        sv.neg_strand_hap1 = neg_strand
    elif hap == 2:
        sv.ref_start_best_hap2 = ref_start_best
        sv.ref_end_best_hap2 = ref_end_best
        sv.query_start_best_hap2 = query_start_best
        sv.query_end_best_hap2 = query_end_best
        sv.query_name_hap2 = contig_name_dict[start_contig_name_ctr]
        sv.neg_strand_hap2 = neg_strand

    return True

def align_before_after_comb(query_seq, ref_seq_1, ref_seq_2):
    #used third_fil when build groups
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    #aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    #aligner.score_only = True
    alignment_beforeSV = aligner.score(query_seq, ref_seq_1)
    alignment_afterSV = aligner.score(query_seq, ref_seq_2)
        
    return alignment_beforeSV, alignment_afterSV
    
def get_comb_vali_info(comb_sv, hap, interval, contig_name_list, contig_pos_list, contig_name_dict, if_hg38, ref_rec, query_fasta_file, sv_idx_dict, region_len_m):
    if not get_comb_intervals_len_only(comb_sv, if_hg38, interval, contig_name_list, contig_pos_list, contig_name_dict, hap, ref_rec, region_len_m):
        return False
    
    if hap == 1:
        query_rec = query_fasta_file.fetch(comb_sv.query_name_hap1)
        ref_start = comb_sv.ref_start_best_hap1
        ref_end = comb_sv.ref_end_best_hap1
        query_start = comb_sv.query_start_best_hap1
        query_end = comb_sv.query_end_best_hap1
        neg_strand = comb_sv.neg_strand_hap1
    elif hap == 2:
        query_rec = query_fasta_file.fetch(comb_sv.query_name_hap2)
        ref_start = comb_sv.ref_start_best_hap2
        ref_end = comb_sv.ref_end_best_hap2
        query_start = comb_sv.query_start_best_hap2
        query_end = comb_sv.query_end_best_hap2
        neg_strand = comb_sv.neg_strand_hap2
    
    if query_start >= len(query_rec) or query_end >= len(query_rec):
        message = "bad_query_pos"
        #write_err(output_file_name, message, g)
        return False    
    
    #definitely need alignment score for comb_sv
    
    query_frag = query_rec[query_start:query_end]
    ref_frag = ref_rec[ref_start:ref_end]
    ref_after_sv_frag = get_comb_ref_frag_after_sv(comb_sv, ref_rec, sv_idx_dict, ref_start, ref_end)
    
    #neg strand
    if neg_strand:
        seq = Seq(query_frag)
        query_frag = seq.reverse_complement()

    #get to upper case
    ref_frag = ref_frag.upper()
    query_frag = query_frag.upper()
    ref_after_sv_frag = ref_after_sv_frag.upper()
    
    #test
#     print(len(ref_frag), len(ref_after_sv_frag), len(ref_after_sv_frag)-len(ref_frag))

    #TODO: find appropriate alignment parameters
    #paras: match, mismatch, open gap, extend gap
    alignment_beforeSV, alignment_afterSV = align_before_after_comb(query_frag, ref_frag, ref_after_sv_frag)    
    
    if hap == 1:
        comb_sv.analyzed_hap1 = True
        comb_sv.score_before_hap1 = alignment_beforeSV
        comb_sv.score_after_hap1 = alignment_afterSV
        comb_sv.len_query_hap1 = query_end - query_start + 1
        comb_sv.len_ref_hap1 = ref_end - ref_start + 1
    elif hap == 2:
        comb_sv.analyzed_hap2 = True
        comb_sv.score_before_hap2 = alignment_beforeSV
        comb_sv.score_after_hap2 = alignment_afterSV
        comb_sv.len_query_hap2 = query_end - query_start + 1
        comb_sv.len_ref_hap2 = ref_end - ref_start + 1
    
    return True

def get_comb_ref_frag_after_sv(comb_sv, ref_rec, sv_idx_dict, ref_start, ref_end):
    cur_seq = ""
    cur_pos = ref_start
    
    #test
#     print("ref_start", ref_start)
    
    #note comb_sv sv are not overlapping
    sv_idx_list = comb_sv.idx
    for idx in sv_idx_list:
        sv = sv_idx_dict[idx]
        sv_pos = sv.sv_pos
        sv_stop = sv.sv_stop
        
        cur_seq += ref_rec[cur_pos:sv_pos]
        
        #test
#         print("cur_pos", cur_pos, "sv_pos", sv_pos, "cur_len", len(cur_seq), len(ref_rec))
        
        if sv.sv_type == "INS":
            cur_seq += sv.ins_seq
            #also skipped ref seq of INS, if any
            cur_pos = sv_stop
        elif sv.sv_type == "DEL":
            #skipped deleted seq
            cur_pos = sv_stop
            
    #test
#     print("cur_pos", cur_pos, "ref_end", ref_end)
    assert cur_pos <= ref_end
    assert len(ref_rec) >= ref_end
    
    cur_seq += ref_rec[cur_pos:ref_end]
    return cur_seq
    
    
    
##################################################################
##################################################################
#functions from agg_sv
