## %%
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
from Bio import Align
import mappy
import os

#function to get seq record from a fasta file
#seq_name: string, file_name:string with suffix
def getSeqRec(seq_name, file_name):
    #fasta_file = SeqIO.parse(file_name, "fasta")
    #for record in fasta_file:
    #    if seq_name == record.id.strip():
    #       return record
    fasta_file = pysam.FastaFile(file_name)
    seq = fasta_file.fetch(seq_name)
    return seq

def getQuerySeqRec(seq_name):
    seq = query_fasta_file.fetch(seq_name)
    return seq

#get the most close 'interval(500)' before sv_pos
def getRefStart(sv_pos, interval):
    start = math.floor(sv_pos/interval) * interval
    return start

#get the most close 'interval(500)' after sv_stop
def getRefEnd(sv_stop, interval):
    end = math.ceil(sv_stop/interval) * interval
    return end

#lookup in a dict given a key
def lookupDict(key, dictionary):
    result = dictionary.get(key)
    return result

#get depth of coverage given position of a chromosome and a target bam_file
#TODO: probably will not do read-depth filtering for giab callset:
#1. they should be reliable;
#2. they are a combination of multiple sources: don't know depth;
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

#write error cases with error messages
def write_err(output_file_name, message):
    f = open(output_file_name, "a")
    f.write(str(counter) + "\t")
    f.write(str(0) + "\t")
    f.write(str(0) + "\t")
    f.write(str(-1) + "\t")
    f.write(str(0) + "\t")
    f.write(str(0) + "\t")
    f.write(str(0) + "\t")
    f.write(str(0) + "\t")
    f.write(str(0) + "\t")
    f.write(str(0) + "\t")
    f.write(str(0) + "\t")
    f.write(str(0) + "\t")
    f.write(str(0) + "\t")
    f.write(str(0) + "\t")
    f.write(str(0) + "\t")
    f.write(message + "\t")
    f.write("\n")
    f.close()
    
#return chr name as an integer
#X = 23, Y = 24
def get_int_chr_name(name, if_hg38):
    #remove "chr"
    if if_hg38:
        chr_name = name[3:]
    else:
        chr_name = name
    
    if chr_name.upper() == 'X':
        int_name = 23
    elif chr_name.upper() == 'Y':
        int_name = 24
    else:
        int_name = int(chr_name)
        
    return int_name
    
#build map
def build_map(chr_len, interval, liftover_file, if_hg38):
    contig_name_list = []
    contig_pos_list = []
    for i in range(1, 25):
        lo_length = chr_len[i-1]
        #set the initial values to -1
        contig_name_list.append(np.zeros(int(lo_length/interval) + 1, dtype='int16') - 1)
        contig_pos_list.append(np.zeros(int(lo_length/interval) + 1, dtype='uint32'))

    #build a dictionary for contig names: to avoid store too many str
    contig_name_dict = dict()

    with open(liftover_file) as f:
        contig_name_ctr = -1
        pre_contig_name = ""
        for line in f:
            record = line.strip().split()
            
            int_ref_name = get_int_chr_name(record[4], if_hg38)
            ref_pos = int(record[5])
            contig_pos = int(record[1])

            #store contig name in a dict to save memory
            contig_name = record[0]
            if contig_name != pre_contig_name:
                #contig_name_ctr starting from 0
                contig_name_ctr += 1
                contig_name_dict[contig_name_ctr] = contig_name
                pre_contig_name = contig_name

            #chr index: 0-23
            lo_list_index = int_ref_name - 1
            #pos in the chr list
            lo_list_pos = int(ref_pos//interval)

            contig_name_list[lo_list_index][lo_list_pos] = contig_name_ctr
            contig_pos_list[lo_list_index][lo_list_pos] = contig_pos
    f.close()
    return contig_name_list, contig_pos_list, contig_name_dict

def build_map_compress(chr_len, interval, liftover_file, if_hg38):
    contig_name_list = []
    contig_pos_list = []
    for i in range(1, 25):
        lo_length = chr_len[i-1]
        #set the initial values to -1
        contig_name_list.append(np.zeros(int(lo_length/interval) + 1, dtype='int16') - 1)
        contig_pos_list.append(np.zeros(int(lo_length/interval) + 1, dtype='uint32'))

    #build a dictionary for contig names: to avoid store too many str
    contig_name_dict = dict()
    
    with open(liftover_file) as f:
        contig_name_ctr = -1
        pre_contig_name = ""
        for line in f:
            record = line.strip().split()
            int_ref_name = get_int_chr_name(record[0], if_hg38)
            #store contig name in a dict to save memory
            contig_name = record[1]
            if contig_name != pre_contig_name:
                #contig_name_ctr starting from 0
                contig_name_ctr += 1
                contig_name_dict[contig_name_ctr] = contig_name
                pre_contig_name = contig_name

            #chr index: 0-23
            lo_list_index = int_ref_name - 1

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
                    ref_pos = int(info_list[0]) + i*interval
                    if forward:
                        contig_pos = int(info_list[1]) + i*interval
                    else:
                        contig_pos = int(info_list[1]) - i*interval
                    #pos in the chr list
                    lo_list_pos = int(ref_pos//interval)

                    contig_name_list[lo_list_index][lo_list_pos] = contig_name_ctr
                    contig_pos_list[lo_list_index][lo_list_pos] = contig_pos        
    f.close()
    return contig_name_list, contig_pos_list, contig_name_dict

#get alignment score of given seq pair
def get_align_score(seq1, seq2):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    #aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    #aligner.score_only = True
    alignment_score = aligner.score(seq1, seq2)
    return alignment_score

#invert sequence
def inversion_seq(seq):
    inverted_seq = ''
    for i in seq:
        inverted_seq = i + inverted_seq
    return inverted_seq

#get different chromosomes' tandem regions start/end index
def get_chr_tandem_shart_end_list(tandem_info, if_hg38):
    start_list = [0] * 24
    end_list = [0] * 24
    
    cur_chr = tandem_info[0][0]
#     print(tandem_info[chr_index_ctr][0], " start ", chr_index_ctr)
    start_list[get_int_chr_name(cur_chr, if_hg38) - 1] = 0
    for i, rec in enumerate(tandem_info[1:]):
        if rec[0] != cur_chr:
#             print(cur_chr, " end ", i-1)
            end_list[get_int_chr_name(cur_chr, if_hg38) - 1] = i-1
            cur_chr = rec[0]
#             print(cur_chr, " start ", i)
            start_list[get_int_chr_name(cur_chr, if_hg38) - 1] = i
#     print(cur_chr, " end ", len(tandem_info)-1)
    end_list[get_int_chr_name(cur_chr, if_hg38) - 1] = len(tandem_info)-1
    return start_list, end_list

#check if in tandem repeats regions
def ol_tandem(ref_name, sv_pos, sv_stop, tandem_start_list, tandem_end_list, if_hg38, tandem_info):
    ref_idx = get_int_chr_name(ref_name, if_hg38)
    tandem_start_idx = tandem_start_list[ref_idx - 1]
    tandem_end_idx = tandem_end_list[ref_idx - 1]
    
    #index of tandem regions that contain sv start/end
    tandem_region_sv_start = tandem_region_sv_end = -1
    
    for i in range(tandem_start_idx, tandem_end_idx + 1):
        if int(tandem_info[i][1]) <= sv_pos and int(tandem_info[i][2]) >= sv_pos:
            tandem_region_sv_start = i
            break
            
    for i in range(tandem_start_idx, tandem_end_idx + 1):
        if int(tandem_info[i][1]) <= sv_stop and int(tandem_info[i][2]) >= sv_stop:
            tandem_region_sv_end = i
            break
    
    return [tandem_region_sv_start, tandem_region_sv_end]

def get_intervals(sv, cur_ref_name, ref_fasta_file, tandem_start_list, tandem_end_list, if_hg38, tandem_info, interval, memory_limit, contig_name_list, contig_pos_list, contig_name_dict, hap, ref_rec):
    if sv.sv_type in ['INS', 'DEL', 'DUP:TANDEM', 'DUP']:
        region_len_m = 500
        
        check_tandem = ol_tandem(sv.ref_name, sv.sv_pos, sv.sv_stop, tandem_start_list, tandem_end_list, if_hg38, tandem_info)
        if check_tandem[0] != -1:
            ref_start_start = max(getRefStart(int(tandem_info[check_tandem[0]][1]), interval) - region_len_m, 0)
            ref_start_end = getRefStart(sv.sv_pos, interval)
        else:
            ref_start_start = max(getRefStart(sv.sv_pos, interval) - region_len_m, 0)
            ref_start_end = getRefStart(sv.sv_pos, interval)
        if check_tandem[1] != -1:
            ref_end_start = getRefEnd(sv.sv_stop, interval)
            ref_end_end = min(getRefEnd(int(tandem_info[check_tandem[1]][2]), interval) + region_len_m, getRefEnd(len(ref_rec), interval) - interval)
        else:
            ref_end_start = getRefEnd(sv.sv_stop, interval)
            ref_end_end = min(getRefEnd(sv.sv_stop, interval) + region_len_m, getRefEnd(len(ref_rec), interval) - interval) 
        
        #if tandem repeat regions too large
        if check_tandem[0] != -1 and check_tandem[1] != -1:
            if int(tandem_info[check_tandem[1]][2]) - int(tandem_info[check_tandem[0]][1]) > memory_limit:
                ref_start_start = max(getRefStart(sv.sv_pos, interval) - region_len_m, 0)
                ref_start_end = getRefStart(sv.sv_pos, interval)
                ref_end_start = getRefEnd(sv.sv_stop, interval)
                ref_end_end = min(getRefEnd(sv.sv_stop, interval) + region_len_m, getRefEnd(len(ref_rec), interval) - interval)

        #first level key: chr index as an int
        int_ref_name = get_int_chr_name(sv.ref_name, if_hg38)
        lo_list_index = int_ref_name - 1
        first_key = lo_list_index

        #second level key: ref list_pos

        if sv.idx % 1000 == 0:
            print(sv.idx)   

        min_diff_len = 100000000
        ref_start_best = ref_start_start
        ref_end_best = ref_end_start
        query_start_best = None
        query_end_best = None
        start_contig_name_ctr = -1
        end_contig_name_ctr = -1

        neg_strand = False

        for ref_end_cand in range(ref_end_start, ref_end_end, interval):
            for ref_start_cand in range(ref_start_end, ref_start_start, -interval):
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
    #if inv, no need to use the flexible interval
    #bc/ flexi int will be found by relative length
    elif sv.sv_type == 'INV':
        region_len = 200

        check_tandem = ol_tandem(sv.ref_name, sv.sv_pos, sv.sv_stop, tandem_start_list, tandem_end_list, if_hg38, tandem_info)
        if check_tandem[0] != -1:
            ref_start = getRefStart(int(tandem_info[check_tandem[0]][1]), interval)
        else:
            ref_start = getRefStart(sv.sv_pos, interval)
        if check_tandem[1] != -1:
            ref_end = getRefEnd(int(tandem_info[check_tandem[1]][2]), interval)
        else:
            ref_end = getRefEnd(sv.sv_stop, interval)

#         ref_start = getRefStart(sv_pos, interval)
#         ref_end = getRefEnd(sv_stop, interval)

        #first level key: chr index as an int
        int_ref_name = get_int_chr_name(sv.ref_name, if_hg38)
        lo_list_index = int_ref_name - 1
        first_key = lo_list_index

        if sv.idx % 1000 == 0:
            print(sv.idx)  

        ref_start = ref_start-region_len
        ref_end = ref_end+region_len
            
        second_key_start = (ref_start)//interval
        second_key_end = (ref_end)//interval

        start_contig_name_ctr = contig_name_list[first_key][second_key_start]
        end_contig_name_ctr = contig_name_list[first_key][second_key_end]

        if start_contig_name_ctr == -1 or end_contig_name_ctr == -1:
            #print("wrong second key")
            message = "wrong_sec_key"
            #write_err(output_file_name, message, g)
            return False

        if start_contig_name_ctr != end_contig_name_ctr:
            #print("Not same contig")
            message = "not_same_contig"
            #write_err(output_file_name, message, g)
            return False

        query_start = contig_pos_list[first_key][second_key_start]
        query_end = contig_pos_list[first_key][second_key_end]

        neg_strand = False
        #in case: negtive strand
        if query_end < query_start:
            tem = query_start
            query_start = query_end
            query_end = tem
            neg_strand = True
            
        if ref_start == sv.sv_pos:
            ref_start = ref_start - 1
        if ref_end == sv.sv_stop:
            ref_end = ref_end + 1
        if query_end == query_start:
            query_end += 1
        
        if hap == 1:
            sv.ref_start_best_hap1 = ref_start
            sv.ref_end_best_hap1 = ref_end
            sv.query_start_best_hap1 = query_start
            sv.query_end_best_hap1 = query_end
            sv.query_name_hap1 = contig_name_dict[start_contig_name_ctr]
            sv.neg_strand_hap1 = neg_strand
        elif hap == 2:
            sv.ref_start_best_hap2 = ref_start
            sv.ref_end_best_hap2 = ref_end
            sv.query_start_best_hap2 = query_start
            sv.query_end_best_hap2 = query_end
            sv.query_name_hap2 = contig_name_dict[start_contig_name_ctr]
            sv.neg_strand_hap2 = neg_strand
        
        return True

#for calls large than memory_limit
def get_large_intervals(sv, cur_ref_name, ref_fasta_file, tandem_start_list, tandem_end_list, if_hg38, tandem_info, interval, memory_limit, contig_name_list, contig_pos_list, contig_name_dict, hap, ref_rec):
        region_len = int(0.5 * abs(sv.length))
        #region_len = 200

        check_tandem = ol_tandem(sv.ref_name, sv.sv_pos, sv.sv_stop, tandem_start_list, tandem_end_list, if_hg38, tandem_info)
        if check_tandem[0] != -1:
            ref_start = getRefStart(int(tandem_info[check_tandem[0]][1]), interval)
        else:
            ref_start = getRefStart(sv.sv_pos, interval)
        if check_tandem[1] != -1:
            ref_end = getRefEnd(int(tandem_info[check_tandem[1]][2]), interval)
        else:
            ref_end = getRefEnd(sv.sv_stop, interval)

#         ref_start = getRefStart(sv_pos, interval)
#         ref_end = getRefEnd(sv_stop, interval)

        #first level key: chr index as an int
        int_ref_name = get_int_chr_name(sv.ref_name, if_hg38)
        lo_list_index = int_ref_name - 1
        first_key = lo_list_index

        if sv.idx % 1000 == 0:
            print(sv.idx)  

        ref_start = max(1, ref_start-region_len)
        second_key_start = int((ref_start)//interval)
        
        ref_end = min((len(contig_name_list[first_key]) - 1) * interval, ref_end+region_len)
        second_key_end = int((ref_end)//interval)

        start_contig_name_ctr = contig_name_list[first_key][second_key_start]
        end_contig_name_ctr = contig_name_list[first_key][second_key_end]

        if start_contig_name_ctr == -1 or end_contig_name_ctr == -1:
            #print("wrong second key")
            message = "wrong_sec_key"
            #write_err(output_file_name, message, g)
            return False

        if start_contig_name_ctr != end_contig_name_ctr:
            #print("Not same contig")
            message = "not_same_contig"
            #write_err(output_file_name, message, g)
            return False

        query_start = contig_pos_list[first_key][second_key_start]
        query_end = contig_pos_list[first_key][second_key_end]

        neg_strand = False
        #in case: negtive strand
        if query_end < query_start:
            tem = query_start
            query_start = query_end
            query_end = tem
            neg_strand = True
            
        if ref_start == sv.sv_pos:
            ref_start = ref_start - 1
        if ref_end == sv.sv_stop:
            ref_end = ref_end + 1
        if query_end == query_start:
            query_end += 1
        
        if hap == 1:
            sv.ref_start_best_hap1 = ref_start
            sv.ref_end_best_hap1 = ref_end
            sv.query_start_best_hap1 = query_start
            sv.query_end_best_hap1 = query_end
            sv.query_name_hap1 = contig_name_dict[start_contig_name_ctr]
            sv.neg_strand_hap1 = neg_strand
        elif hap == 2:
            sv.ref_start_best_hap2 = ref_start
            sv.ref_end_best_hap2 = ref_end
            sv.query_start_best_hap2 = query_start
            sv.query_end_best_hap2 = query_end
            sv.query_name_hap2 = contig_name_dict[start_contig_name_ctr]
            sv.neg_strand_hap2 = neg_strand
        
        return True    

#write vali info
def write_vali_info(g, sv, hap):
    g.write(str(sv.idx) + "\t")
    
    if hap == 1:
        alignment_beforeSV = sv.score_before_hap1
        alignment_afterSV = sv.score_after_hap1
        query_length = sv.len_query_hap1
        ref_length = sv.len_ref_hap1
        ref_start = sv.ref_start_best_hap1
        ref_end = sv.ref_end_best_hap1
    elif hap == 2:
        alignment_beforeSV = sv.score_before_hap2
        alignment_afterSV = sv.score_after_hap2
        query_length = sv.len_query_hap2
        ref_length = sv.len_ref_hap2
        ref_start = sv.ref_start_best_hap2
        ref_end = sv.ref_end_best_hap2

    g.write(str(alignment_beforeSV) + "\t")
    #g.write(str(1) + "\t")
    g.write(str(alignment_afterSV) + "\t")
    #g.write(str(1) + "\t")

    #g.write(str(len(query_start_dic)) + "\t")
    g.write(str(1) + "\t")

    g.write(str(sv.sv_type) + "\t")
    g.write(str(hap) + "\t")
    g.write(str(query_length) + "\t")
    g.write(str(ref_length) + "\t")
    g.write(str(sv.length) + "\t")
    #g.write(str(sv_stop-sv_pos+1) + "\t")
    g.write(str(sv.ref_name) + "\t")
    g.write(str(sv.sv_pos) + "\t")
    g.write(str(sv.sv_stop) + "\t")
    g.write(str(ref_start) + "\t")
    g.write(str(ref_end) + "\t")
    g.write(str(sv.sv_type) + "\t")
    g.write("noerr" + "\t")
    g.write("\n")
    
def align_before_after(output_dir, sv, query_seq, ref_seq_1, ref_seq_2):
    #within length limit
    if not sv.is_third_fil:
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
    else:
        h = open(output_dir+"tmp_query.fasta", "w")
        h.write('>' + str(sv.idx) + "\n")
        h.write(query_seq + "\n")
        h.close()
#         aligner = mappy.Aligner(fn_idx_in=output_dir+"tmp_query.fasta", scoring=[1,1,2,1])
        aligner = mappy.Aligner(fn_idx_in=output_dir+"tmp_query.fasta")
        #if not alignment: raise Exception("ERROR: failed to load/build index")
        aligner_beforeSV = aligner.map(ref_seq_1, seq2=None, cs=False, MD=False)
        aligner_afterSV = aligner.map(ref_seq_2, seq2=None, cs=False, MD=False)

        #test
#         for agt in aligner_beforeSV:
#             alignment_beforeSV = len(query_seq) - (len(ref_seq_1) - agt.mlen)
#             break
#         for agt in aligner_afterSV:
#             alignment_afterSV = len(query_seq) - (len(ref_seq_2) - agt.mlen)
#             break        
    
        try:
            agt_before = next(aligner_beforeSV)
        except:
            os.remove(output_dir+"tmp_query.fasta")
            return None, None
        
        try:
            agt_after = next(aligner_afterSV)
        except:
            os.remove(output_dir+"tmp_query.fasta")
            return None, None
        
        alignment_beforeSV = len(query_seq) - (len(ref_seq_1) - agt_before.mlen)
        alignment_afterSV = len(query_seq) - (len(ref_seq_2) - agt_after.mlen)

        os.remove(output_dir+"tmp_query.fasta")
        
    return alignment_beforeSV, alignment_afterSV
    
#get vcf file and run score_callset on each SV record
def get_vali_info(output_dir, vcf_file, query_file, hap, ref_file, interval, 
                  contig_name_list, contig_pos_list, contig_name_dict, memory_limit, if_hg38, chr_list,
                  tandem_start_list, tandem_end_list, tandem_info, sv_list, seq_resolved):
    #query_file = query_file2
    #hap = 2
    name_str = "assem" + str(hap)
    chromosome = "all"
    output_file_name = output_dir + "align_info_" + name_str + "_chr" + chromosome + ".txt"
    
    
    query_fasta_file = pysam.FastaFile(query_file)
    
    ref_fasta_file = pysam.FastaFile(ref_file)
    cur_ref_name = ""

    g = open(output_file_name, "w")
    for sv in sv_list:
        
        #filtered by centromere, non-cov
        if sv.is_sec_fil:
            continue
            
        if cur_ref_name != sv.ref_name:
            cur_ref_name = sv.ref_name
            ref_rec = ref_fasta_file.fetch(cur_ref_name)
            
        #filtered by size
        if not sv.is_third_fil:
            #Search for the best second_key_start and second_key_end in the regions
            if not get_intervals(sv, cur_ref_name, ref_fasta_file, tandem_start_list, tandem_end_list, if_hg38, tandem_info, interval, memory_limit, contig_name_list, contig_pos_list, contig_name_dict, hap, ref_rec):
                continue
        else:
            #get the intervals, won't search 
            if not get_large_intervals(sv, cur_ref_name, ref_fasta_file, tandem_start_list, tandem_end_list, if_hg38, tandem_info, interval, memory_limit, contig_name_list, contig_pos_list, contig_name_dict, hap, ref_rec):
                continue
        
        if hap == 1:
            query_rec = query_fasta_file.fetch(sv.query_name_hap1)
            ref_start = sv.ref_start_best_hap1
            ref_end = sv.ref_end_best_hap1
            query_start = sv.query_start_best_hap1
            query_end = sv.query_end_best_hap1
            neg_strand = sv.neg_strand_hap1
        elif hap == 2:
            query_rec = query_fasta_file.fetch(sv.query_name_hap2)
            ref_start = sv.ref_start_best_hap2
            ref_end = sv.ref_end_best_hap2
            query_start = sv.query_start_best_hap2
            query_end = sv.query_end_best_hap2
            neg_strand = sv.neg_strand_hap2

        if query_start >= len(query_rec) or query_end >= len(query_rec):
            message = "bad_query_pos"
            #write_err(output_file_name, message, g)
            continue

        #case 1: DEL
        if sv.sv_type == "DEL":
            #query and ref seq fragment
            
            query_frag = query_rec[query_start:query_end]
            ref_frag = ref_rec[ref_start:ref_end]
            #TODO: this is for DEL
            #TODO: check +-1
            ref_afterSV_frag1 = ref_rec[ref_start:sv.sv_pos]
            ref_afterSV_frag2 = ref_rec[sv.sv_stop:ref_end]

            #alignment starts here
            if neg_strand:
                seq = Seq(query_frag)
                query_frag = seq.reverse_complement()
                #query_frag = query_frag.reverse_complement()

            #get to upper case
            ref_frag = ref_frag.upper()
            ref_afterSV_frag1 = ref_afterSV_frag1.upper()
            ref_afterSV_frag2 = ref_afterSV_frag2.upper()
            query_frag = query_frag.upper()

            #TODO: find appropriate alignment parameters
            #paras: match, mismatch, open gap, extend gap
            alignment_beforeSV, alignment_afterSV = align_before_after(output_dir, sv, str(query_frag), str(ref_frag), str(ref_afterSV_frag1) + str(ref_afterSV_frag2))

        #case 2: INS
        elif sv.sv_type in ["INS"]:
            #query and ref seq fragment
            query_frag = query_rec[query_start:query_end]
            ref_frag = ref_rec[ref_start:ref_end]
            
            if seq_resolved:
                ins_seq = sv.ins_seq
                ref_afterSV_frag1 = ref_rec[ref_start:sv.sv_pos]
                ref_afterSV_frag2 = ref_rec[sv.sv_stop:ref_end]

                
                if neg_strand:
                    seq = Seq(query_frag)
                    query_frag = seq.reverse_complement()
#                     seq = Seq(ins_seq)
#                     ins_seq = seq.reverse_complement()
                
                #get to upper case
                ref_frag = ref_frag.upper()
                ref_afterSV_frag1 = ref_afterSV_frag1.upper()
                ref_afterSV_frag2 = ref_afterSV_frag2.upper()
                query_frag = query_frag.upper()
                ins_seq = ins_seq.upper()
                
                alignment_beforeSV, alignment_afterSV = align_before_after(output_dir, sv, str(query_frag), str(ref_frag), str(ref_afterSV_frag1) + ins_seq + str(ref_afterSV_frag2))
            else:
                #for INS, not using relative score
                alignment_beforeSV = 1
                alignment_afterSV = 1

        #case 3: INV
        elif sv.sv_type == "INV":
            #query and ref seq fragment
            query_frag = query_rec[query_start:query_end]
            ref_frag = ref_rec[ref_start:ref_end]
            ref_afterSV_frag1 = ref_rec[ref_start:sv.sv_pos]
            ref_afterSV_frag2 = ref_rec[sv.sv_stop:ref_end]
            ref_inv_seq = ref_rec[sv.sv_pos:sv.sv_stop]

            #alignment starts here
            if neg_strand:
                seq = Seq(query_frag)
                query_frag = seq.reverse_complement()
                #query_frag = query_frag.reverse_complement()

            #reversse and complement
            inv_seq = Seq(ref_inv_seq)
            ref_inv_seq = inv_seq.reverse_complement()
            #ref_inv_seq = inversion_seq(str(ref_inv_seq))

            #get to upper case
            ref_frag = ref_frag.upper()
            query_frag = query_frag.upper()
            ref_inv_seq = ref_inv_seq.upper()
            ref_afterSV_frag1 = ref_afterSV_frag1.upper()
            ref_afterSV_frag2 = ref_afterSV_frag2.upper()
            
            alignment_beforeSV, alignment_afterSV = align_before_after(output_dir, sv, str(query_frag), str(ref_frag), str(ref_afterSV_frag1) + str(ref_inv_seq) + str(ref_afterSV_frag2))

#             #within length limit
#             if not sv.is_third_fil:
#                 #TODO: find a appropriate alignment parameters
#                 #paras: match, mismatch, open gap, extend gap
#                 aligner = Align.PairwiseAligner()
#                 aligner.mode = 'global'
#                 #aligner.mode = 'local'
#                 aligner.match_score = 1
#                 aligner.mismatch_score = -1
#                 aligner.open_gap_score = -1
#                 aligner.extend_gap_score = -0.5
#                 #aligner.score_only = True
#                 alignment_beforeSV = aligner.score(str(query_frag), str(ref_frag))
#                 alignment_afterSV = aligner.score(str(query_frag), str(ref_afterSV_frag1) + str(ref_inv_seq)
#                                         + str(ref_afterSV_frag2))
        
        #case 4: tandem dup
        elif sv.sv_type in ['DUP:TANDEM', 'DUP']:
            #query and ref seq fragment
            query_frag = query_rec[query_start:query_end]
            ref_frag = ref_rec[ref_start:ref_end]
            #TODO: this is for DEL
            #TODO: check +-1
            ref_afterSV_frag1 = ref_rec[ref_start:sv.sv_stop]
            ref_afterSV_frag2 = ref_rec[sv.sv_stop:ref_end]
            dup_seq = ref_rec[sv.sv_pos:sv.sv_stop]

            #alignment starts here

            if neg_strand:
                seq = Seq(query_frag)
                query_frag = seq.reverse_complement()
                #query_frag = query_frag.reverse_complement()

            #get to upper case
            ref_frag = ref_frag.upper()
            ref_afterSV_frag1 = ref_afterSV_frag1.upper()
            ref_afterSV_frag2 = ref_afterSV_frag2.upper()
            dup_seq = dup_seq.upper()
            query_frag = query_frag.upper()
            
            alignment_beforeSV, alignment_afterSV = align_before_after(output_dir, sv, str(query_frag), str(ref_frag), str(ref_afterSV_frag1) + str(dup_seq) + str(ref_afterSV_frag2))

#             #within length limit
#             if not sv.is_third_fil:
#                 #TODO: find a appropriate alignment parameters
#                 #paras: match, mismatch, open gap, extend gap
#                 aligner = Align.PairwiseAligner()
#                 aligner.mode = 'global'
#                 #aligner.mode = 'local'
#                 aligner.match_score = 1
#                 aligner.mismatch_score = -1
#                 aligner.open_gap_score = -1
#                 aligner.extend_gap_score = -0.5
#                 #aligner.score_only = True
#                 alignment_beforeSV = aligner.score(str(query_frag), str(ref_frag))
#                 alignment_afterSV = aligner.score(str(query_frag), str(ref_afterSV_frag1) + str(dup_seq) + str(ref_afterSV_frag2))

        if hap == 1:
            sv.analyzed_hap1 = True
            sv.score_before_hap1 = alignment_beforeSV
            sv.score_after_hap1 = alignment_afterSV
            sv.len_query_hap1 = len(query_frag)
            sv.len_ref_hap1 = len(ref_frag)
        elif hap == 2:
            sv.analyzed_hap2 = True
            sv.score_before_hap2 = alignment_beforeSV
            sv.score_after_hap2 = alignment_afterSV
            sv.len_query_hap2 = len(query_frag)
            sv.len_ref_hap2 = len(ref_frag)
            
        write_vali_info(g, sv, hap)
    g.close()


def main():
    #get command line input
    output_dir = sys.argv[1] + "/"
    vcf_file = sys.argv[2]
    #ref fasta file
    ref_file = sys.argv[3]
    #assembly fasta files
    query_file1 = sys.argv[4]
    query_file2 = sys.argv[5]
    liftover_file1 = sys.argv[6]
    liftover_file2 = sys.argv[7]
    if_hg38_input = sys.argv[8]
    tandem_file = sys.argv[9]

    #constants
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
                    "chr21", "chr22", "chrX"]
    else:
        chr_list = ["1", "2", "3", "4", "5",
                    "6", "7", "8", "9", "10",
                    "11", "12", "13", "14", "15",
                    "16", "17", "18", "19", "20",
                    "21", "22", "X"]
    #interval length
    interval = 20
    #approximate length of chromosomes
    chr_len = [250000000, 244000000, 199000000, 192000000, 182000000, 
                172000000, 160000000, 147000000, 142000000, 136000000, 
                136000000, 134000000, 116000000, 108000000, 103000000, 
                90400000, 83300000, 80400000, 59200000, 64500000, 
                48200000, 51400000, 157000000, 59400000]
    #max length of allowed alignment
    memory_limit = 50000

    #tandem repeats regions file
    with open(tandem_file) as f:
        reader = csv.reader(f, delimiter="\t")
        tandem_info = list(reader)
    f.close()    
    
    #get tandem start and end list
    tandem_start_list, tandem_end_list = get_chr_tandem_shart_end_list(tandem_info, if_hg38)
     
    #build map and get validation info haplotype 1
    contig_name_list, contig_pos_list, contig_name_dict = build_map(chr_len, interval, liftover_file1, if_hg38)
    get_vali_info(output_dir, vcf_file, query_file1, 1, ref_file, interval, 
                  contig_name_list, contig_pos_list, contig_name_dict, memory_limit, if_hg38, chr_list,
                  tandem_start_list, tandem_end_list, tandem_info)
    
    #build map and get validation info haplotype 2
    contig_name_list, contig_pos_list, contig_name_dict = build_map(chr_len, interval, liftover_file2, if_hg38)
    get_vali_info(output_dir, vcf_file, query_file2, 2, ref_file, interval, 
                  contig_name_list, contig_pos_list, contig_name_dict, memory_limit, if_hg38, chr_list,
                  tandem_start_list, tandem_end_list, tandem_info)
    

#main function and pack all the steps by both haplotypes
if __name__ == "__main__":
    main()
