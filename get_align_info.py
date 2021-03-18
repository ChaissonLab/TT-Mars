# %%
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

#get the most close 'interval(500)' after sv_end
def getRefEnd(sv_end, interval):
    end = math.ceil(sv_end/interval) * interval
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

#get alignment score of given seq pair
def get_align_score(seq1, seq2):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
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

#get vcf file and run score_callset on each SV record
def get_vali_info(output_dir, vcf_file, query_file, hap, ref_file, interval, 
                  contig_name_list, contig_pos_list, contig_name_dict, memory_limit, if_hg38, chr_list):
    f = pysam.VariantFile(vcf_file,'r')
    #query_file = query_file2
    #hap = 2
    name_str = "assem" + str(hap)
    chromosome = "all"
    output_file_name = output_dir + "align_info_" + name_str + "_chr" + chromosome + ".txt"

    query_fasta_file = pysam.FastaFile(query_file)
    ref_fasta_file = pysam.FastaFile(ref_file)
    ref_name = ""

    g = open(output_file_name, "w")
    for counter, rec in enumerate(f.fetch()):
        name = rec.chrom
        
        if name not in chr_list:
            continue

        #choose target chr
        #if str(name) != chromosome:
        #    continue

        sv_type = rec.info['SVTYPE']

        #for testing purpose, include dup as tandem dup
        if sv_type not in ['DEL', 'INS', 'INV', 'DUP:TANDEM', 'DUP']:
        #test INV DUP
        #if sv_type not in ['INV', 'DUP:TANDEM', 'DUP']:
        #if sv_type != 'DEL' and sv_type != 'INS':
            #print("Wrong type!")
            continue

        #sv_len = abs(rec.info["SVLEN"][0])
        #without abs value: rela length closer to 1 the better
        sv_len = rec.info['SVLEN'][0]

        sv_pos = rec.pos
        sv_end = rec.stop
        sv_ref_seq = rec.ref
        #TODO: update here when multi-calls in one records allowed
        #sv_query_seq = rec.alts[0]

        #find given fragment (on ref) to queries

        #get ref seq
        if ref_name != name:
            ref_name = name
            ref_rec = ref_fasta_file.fetch(ref_name)

        ref_start = getRefStart(sv_pos, interval)
        ref_end = getRefEnd(sv_end, interval)

        #first level key: chr index as an int
        int_ref_name = get_int_chr_name(ref_name, if_hg38)
        lo_list_index = int_ref_name - 1
        first_key = lo_list_index

        #second level key: ref list_pos

        if counter % 1000 == 0:
            print(counter)

        #if ins or del or tandem dup
        #Search for the best second_key_start and second_key_end in the regions
        if sv_type in ['INS', 'DEL', 'DUP:TANDEM', 'DUP']:
        #if sv_type == 'INS' or sv_type == 'DEL':
            region_len_m = 500
            min_diff_len = 100000000
            ref_start_best = ref_start
            ref_end_best = ref_end
            query_start_best = None
            query_end_best = None
            start_contig_name_ctr = -1
            end_contig_name_ctr = -1

            neg_strand = False

            for ref_end_cand in range(ref_end, ref_end+region_len_m, interval):
                for ref_start_cand in range(ref_start, ref_start-region_len_m, -interval):

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
                    if abs((query_end - query_start) - (ref_end_cand - ref_start_cand) - sv_len) < min_diff_len:
                        min_diff_len = abs((query_end - query_start) - (ref_end_cand - ref_start_cand) - sv_len)
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
                print("Wrong query pos")
                message = "Wrong_query_pos"
                #write_err(output_file_name, message, g)
                continue

            ref_start = ref_start_best
            ref_end = ref_end_best
            query_start = query_start_best
            query_end = query_end_best
        
                
#         elif sv_type == 'INV':        
#             region_len_m = 500
#             max_rela_score = -100000000
#             ref_start_best = ref_start
#             ref_end_best = ref_end
#             query_start_best = None
#             query_end_best = None
#             start_contig_name_ctr = -1
#             end_contig_name_ctr = -1
#             cur_contig_name_ctr = -100

#             neg_strand = False

#             for ref_end_cand in range(ref_end, ref_end+region_len_m, interval):
#                 for ref_start_cand in range(ref_start, ref_start-region_len_m, -interval):

#                     #second_key_start = str(ref_start_cand)
#                     #second_key_end = str(ref_end_cand)
#                     second_key_start = ref_start_cand//interval
#                     second_key_end = ref_end_cand//interval

#                     start_contig_name_ctr_cand = contig_name_list[first_key][second_key_start]
#                     end_contig_name_ctr_cand = contig_name_list[first_key][second_key_end]

#                     if start_contig_name_ctr_cand == -1 or end_contig_name_ctr_cand == -1:
#                         #print("wrong second key")
#                         message = "wrong_sec_key"
#                         #write_err(output_file_name, message, g)
#                         continue

#                     if start_contig_name_ctr_cand != end_contig_name_ctr_cand:
#                         #print("Not same contig")
#                         message = "not_same_contig"
#                         #write_err(output_file_name, message, g)
#                         continue

#                     query_start = contig_pos_list[first_key][second_key_start]
#                     query_end = contig_pos_list[first_key][second_key_end]

#                     neg_strand_tep = False
#                     #in case: negtive strand
#                     if query_end < query_start:
#                         tem = query_start
#                         query_start = query_end
#                         query_end = tem
#                         neg_strand_tep = True

#                     #take the best relative length to be the optimal interval
#                     if cur_contig_name_ctr != start_contig_name_ctr_cand:
#                         cur_contig_name_ctr == start_contig_name_ctr_cand
#                         query_name = contig_name_dict[start_contig_name_ctr_cand]
#                         query_rec = query_fasta_file.fetch(query_name)
        
#                     ref_frag = ref_rec[ref_start_cand:ref_end_cand]
#                     ref_afterSV_frag1 = ref_rec[ref_start_cand:sv_pos]
#                     ref_afterSV_frag2 = ref_rec[sv_end:ref_end_cand]
#                     ref_inv_seq = ref_rec[sv_pos:sv_end]     
#                     query_frag = query_rec[query_start:query_end]   
                    
#                     if neg_strand_tep:
#                         from Bio.Seq import Seq
#                         seq = Seq(query_frag)
#                         query_frag = seq.reverse_complement()

#                     #reversse and complement
#                     from Bio.Seq import Seq
#                     inv_seq = Seq(ref_inv_seq)
#                     ref_inv_seq = inv_seq.reverse_complement()

#                     #get to upper case
#                     ref_frag = ref_frag.upper()
#                     query_frag = query_frag.upper()
#                     ref_inv_seq = ref_inv_seq.upper()
#                     ref_afterSV_frag1 = ref_afterSV_frag1.upper()
#                     ref_afterSV_frag2 = ref_afterSV_frag2.upper()

#                     if len(str(query_frag)) > memory_limit or len(str(ref_frag)) > memory_limit:
#                         message = "memory_limit"
#                         #write_err(output_file_name, message)
#                         continue

#                     if len(str(query_frag)) == 0 or len(str(ref_frag)) == 0 or len(str(ref_afterSV_frag1)) + len(str(ref_inv_seq)) + len(str(ref_afterSV_frag2)) == 0:
#                         message = "Wrong seq!!!"
#                         #write_err(output_file_name, message)
#                         continue
                    
#                     cur_before_score = get_align_score(str(query_frag), str(ref_frag))
#                     cur_after_score = get_align_score(str(query_frag), str(ref_afterSV_frag1) + str(ref_inv_seq) + str(ref_afterSV_frag2))
                    
#                     if cur_before_score == 0:
#                         continue
#                     cur_rela_score = round((cur_after_score - cur_before_score)/abs(cur_before_score), 2)
                    
#                     if cur_rela_score > max_rela_score:
#                         max_rela_score = cur_rela_score
#                         ref_start_best = ref_start_cand
#                         ref_end_best = ref_end_cand
#                         query_start_best = query_start
#                         query_end_best = query_end
#                         start_contig_name_ctr = start_contig_name_ctr_cand
#                         end_contig_name_ctr = end_contig_name_ctr_cand

#                         if neg_strand_tep:
#                             neg_strand = True
#                         else:
#                             neg_strand = False

#             if query_start_best == None or query_end_best == None:
#                 print("Wrong query pos")
#                 message = "Wrong_query_pos"
#                 #write_err(output_file_name, message, g)
#                 continue

#             ref_start = ref_start_best
#             ref_end = ref_end_best
#             query_start = query_start_best
#             query_end = query_end_best
        
        #if inv, no need to use the flexible interval
        #bc/ flexi int will be found by relative length
        elif sv_type == 'INV':
            region_len = 200
            
            #second_key_start = str(ref_start_cand)
            #second_key_end = str(ref_end_cand)
            second_key_start = (ref_start-region_len)//interval
            second_key_end = (ref_end+region_len)//interval
            
            start_contig_name_ctr = contig_name_list[first_key][second_key_start]
            end_contig_name_ctr = contig_name_list[first_key][second_key_end]

            if start_contig_name_ctr == -1 or end_contig_name_ctr == -1:
                #print("wrong second key")
                message = "wrong_sec_key"
                #write_err(output_file_name, message, g)
                continue

            if start_contig_name_ctr != end_contig_name_ctr:
                #print("Not same contig")
                message = "not_same_contig"
                #write_err(output_file_name, message, g)
                continue

            query_start = contig_pos_list[first_key][second_key_start]
            query_end = contig_pos_list[first_key][second_key_end]

            neg_strand = False
            #in case: negtive strand
            if query_end < query_start:
                tem = query_start
                query_start = query_end
                query_end = tem
                neg_strand = True    
                
        #Get the sequences
        if ref_start == sv_pos:
            ref_start = ref_start - 1

        if ref_end == sv_end:
            ref_end = ref_end + 1

        if query_end == query_start:
            query_end += 1
            #continue

        query_name = contig_name_dict[start_contig_name_ctr]
        query_rec = query_fasta_file.fetch(query_name)

        if query_start >= len(query_rec) or query_end >= len(query_rec):
            message = "bad_query_pos"
            #write_err(output_file_name, message, g)
            continue

        #case 1: DEL
        if sv_type == "DEL":
            #query and ref seq fragment
            query_frag = query_rec[query_start:query_end]
            ref_frag = ref_rec[ref_start:ref_end]
            #TODO: this is for DEL
            #TODO: check +-1
            ref_afterSV_frag1 = ref_rec[ref_start:sv_pos]
            ref_afterSV_frag2 = ref_rec[sv_end:ref_end]

            #alignment starts here
            if neg_strand:
                from Bio.Seq import Seq
                seq = Seq(query_frag)
                query_frag = seq.reverse_complement()
                #query_frag = query_frag.reverse_complement()

            #get to upper case
            ref_frag = ref_frag.upper()
            ref_afterSV_frag1 = ref_afterSV_frag1.upper()
            ref_afterSV_frag2 = ref_afterSV_frag2.upper()
            query_frag = query_frag.upper()

            #TODO: fragments too long will cause memory problem

            if len(str(query_frag)) > memory_limit or len(str(ref_frag)) > memory_limit:
                message = "memory_limit"
                #write_err(output_file_name, message, g)
                continue
            if len(str(ref_afterSV_frag1)) + len(str(ref_afterSV_frag2)) > memory_limit:
                message = "memory_limit"
                #write_err(output_file_name, message, g)
                continue

            #TODO: find appropriate alignment parameters
            #paras: match, mismatch, open gap, extend gap

            aligner = Align.PairwiseAligner()
            aligner.mode = 'global'
            aligner.match_score = 1
            aligner.mismatch_score = -1
            aligner.open_gap_score = -1
            aligner.extend_gap_score = -0.5
            #aligner.score_only = True
            alignment_beforeSV = aligner.score(str(query_frag), str(ref_frag))
            alignment_afterSV = aligner.score(str(query_frag), str(ref_afterSV_frag1) + str(ref_afterSV_frag2))

            #alignment_beforeSV = pairwise2.align.globalms(str(query_frag), str(ref_frag), 1, -1, -1, -0.5, score_only = True)
            #alignment_afterSV = pairwise2.align.globalms(str(query_frag), str(ref_afterSV_frag1) 
            #                        + str(ref_afterSV_frag2), 1, -1, -1, -0.5, score_only = True)
            #get correct query info format

        #case 2: INS
        elif sv_type in ["INS"]:
            #query and ref seq fragment
            query_frag = query_rec[query_start:query_end]
            ref_frag = ref_rec[ref_start:ref_end]
            #TODO: this is for DEL
            #TODO: check +-1
#             ref_afterSV_frag1 = ref_rec[ref_start:sv_pos]
#             ref_afterSV_frag2 = ref_rec[sv_end:ref_end]
#             ins_seq = sv_query_seq

#             #alignment starts here

#             if neg_strand:
#                 from Bio.Seq import Seq
#                 seq = Seq(query_frag)
#                 query_frag = seq.reverse_complement()
#                 #query_frag = query_frag.reverse_complement()

#             #get to upper case
#             ref_frag = ref_frag.upper()
#             ref_afterSV_frag1 = ref_afterSV_frag1.upper()
#             ref_afterSV_frag2 = ref_afterSV_frag2.upper()
#             ins_seq = ins_seq.upper()
#             query_frag = query_frag.upper()

#             #TODO: fragments too long will cause memory problem

#             if len(str(query_frag)) > memory_limit or len(str(ref_frag)) > memory_limit:
#                 message = "memory_limit"
#                 #write_err(output_file_name, message, g)
#                 continue
#             if len(str(ref_afterSV_frag1)) + len(str(ins_seq)) + len(str(ref_afterSV_frag2)) > memory_limit:
#                 message = "memory_limit"
#                 #write_err(output_file_name, message, g)
#                 continue

#             if len(str(query_frag)) == 0 or len(str(query_frag)) == 0 or len(str(ref_afterSV_frag1)) + len(str(ins_seq)) + len(str(ref_afterSV_frag2)) == 0:
#                 message = "Wrong seq!!!"
#                 #write_err(output_file_name, message, g)
#                 continue

#             #TODO: find a appropriate alignment parameters
#             #paras: match, mismatch, open gap, extend gap
#             aligner = Align.PairwiseAligner()
#             aligner.mode = 'global'
#             aligner.match_score = 1
#             aligner.mismatch_score = -1
#             aligner.open_gap_score = -1
#             aligner.extend_gap_score = -0.5
#             #aligner.score_only = True
#             alignment_beforeSV = aligner.score(str(query_frag), str(ref_frag))
#             alignment_afterSV = aligner.score(str(query_frag), str(ref_afterSV_frag1) + str(ins_seq) + str(ref_afterSV_frag2))

            #for INS, not using relative score
            alignment_beforeSV = 1
            alignment_afterSV = 1

        #case 3: INV
        elif sv_type == "INV":
            #query and ref seq fragment
            query_frag = query_rec[query_start:query_end]
            ref_frag = ref_rec[ref_start-region_len:ref_end+region_len]
            ref_afterSV_frag1 = ref_rec[ref_start-region_len:sv_pos]
            ref_afterSV_frag2 = ref_rec[sv_end:ref_end+region_len]
            ref_inv_seq = ref_rec[sv_pos:sv_end]
            
#             query_frag = query_rec[query_start:query_end]
#             ref_frag = ref_rec[ref_start:ref_end]
#             ref_afterSV_frag1 = ref_rec[ref_start:sv_pos]
#             ref_afterSV_frag2 = ref_rec[sv_end:ref_end]
#             ref_inv_seq = ref_rec[sv_pos:sv_end]

            #alignment starts here

            if neg_strand:
                from Bio.Seq import Seq
                seq = Seq(query_frag)
                query_frag = seq.reverse_complement()
                #query_frag = query_frag.reverse_complement()

            #inversed query seq

            #reverse only no complement
            #from Bio.Seq import MutableSeq
            #inv_seq = MutableSeq(ref_inv_seq)
            #inv_seq.reverse()
            #ref_inv_seq = inv_seq.toseq()

            #reversse and complement
            from Bio.Seq import Seq
            inv_seq = Seq(ref_inv_seq)
            ref_inv_seq = inv_seq.reverse_complement()
            #ref_inv_seq = inversion_seq(str(ref_inv_seq))

            #get to upper case
            ref_frag = ref_frag.upper()
            query_frag = query_frag.upper()
            ref_inv_seq = ref_inv_seq.upper()
            ref_afterSV_frag1 = ref_afterSV_frag1.upper()
            ref_afterSV_frag2 = ref_afterSV_frag2.upper()

            if len(str(query_frag)) > memory_limit or len(str(ref_frag)) > memory_limit:
                message = "memory_limit"
                #write_err(output_file_name, message)
                continue

            if len(str(query_frag)) == 0 or len(str(ref_frag)) == 0 or len(str(ref_afterSV_frag1)) + len(str(ref_inv_seq)) + len(str(ref_afterSV_frag2)) == 0:
                message = "Wrong seq!!!"
                #write_err(output_file_name, message)
                continue

            #TODO: find a appropriate alignment parameters
            #paras: match, mismatch, open gap, extend gap
            aligner = Align.PairwiseAligner()
            aligner.mode = 'global'
            aligner.match_score = 1
            aligner.mismatch_score = -1
            aligner.open_gap_score = -1
            aligner.extend_gap_score = -0.5
            #aligner.score_only = True
            alignment_beforeSV = aligner.score(str(query_frag), str(ref_frag))
            alignment_afterSV = aligner.score(str(query_frag), str(ref_afterSV_frag1) + str(ref_inv_seq)
                                    + str(ref_afterSV_frag2))
        
        #case 4: tandem dup
        elif sv_type in ['DUP:TANDEM', 'DUP']:
            #query and ref seq fragment
            query_frag = query_rec[query_start:query_end]
            ref_frag = ref_rec[ref_start:ref_end]
            #TODO: this is for DEL
            #TODO: check +-1
            ref_afterSV_frag1 = ref_rec[ref_start:sv_end]
            ref_afterSV_frag2 = ref_rec[sv_end:ref_end]
            dup_seq = ref_rec[sv_pos:sv_end]

            #alignment starts here

            if neg_strand:
                from Bio.Seq import Seq
                seq = Seq(query_frag)
                query_frag = seq.reverse_complement()
                #query_frag = query_frag.reverse_complement()

            #get to upper case
            ref_frag = ref_frag.upper()
            ref_afterSV_frag1 = ref_afterSV_frag1.upper()
            ref_afterSV_frag2 = ref_afterSV_frag2.upper()
            dup_seq = dup_seq.upper()
            query_frag = query_frag.upper()

            #TODO: fragments too long will cause memory problem

            if len(str(query_frag)) > memory_limit or len(str(ref_frag)) > memory_limit:
                message = "memory_limit"
                #write_err(output_file_name, message, g)
                continue
            if len(str(ref_afterSV_frag1)) + len(str(dup_seq)) + len(str(ref_afterSV_frag2)) > memory_limit:
                message = "memory_limit"
                #write_err(output_file_name, message, g)
                continue

            if len(str(query_frag)) == 0 or len(str(query_frag)) == 0 or len(str(ref_afterSV_frag1)) + len(str(dup_seq)) + len(str(ref_afterSV_frag2)) == 0:
                message = "Wrong seq!!!"
                #write_err(output_file_name, message, g)
                continue

            #TODO: find a appropriate alignment parameters
            #paras: match, mismatch, open gap, extend gap
            aligner = Align.PairwiseAligner()
            aligner.mode = 'global'
            aligner.match_score = 1
            aligner.mismatch_score = -1
            aligner.open_gap_score = -1
            aligner.extend_gap_score = -0.5
            #aligner.score_only = True
            alignment_beforeSV = aligner.score(str(query_frag), str(ref_frag))
            alignment_afterSV = aligner.score(str(query_frag), str(ref_afterSV_frag1) + str(dup_seq) + str(ref_afterSV_frag2))

            #alignment_beforeSV = pairwise2.align.globalms(str(query_frag), str(ref_frag), 1, -1, -1, -0.5, score_only = True)
            #alignment_afterSV = pairwise2.align.globalms(str(query_frag), str(ref_afterSV_frag1) 
            #                        + str(ref_afterSV_frag2), 1, -1, -1, -0.5, score_only = True)
            #get correct query info format        
        
        
            
        #Need to store the following information in order:
        # SV_index score_before score_after no_found SV_type hap:1 or 2
        # (counter alignment_beforeSV[0][2] alignment_afterSV[0][2] len(query_start_dic) sv_type hap)
        # contig_seg_len ref_seg_len SV_len
        # (len(query_frag) len(ref_frag) sv_end-sv_pos+1)
        # chr_name SV_start SV_end interval_start interval_end
        # (ref_name sv_pos sv_end ref_start ref_end)
        # error message
        # (err_mes)
        # if_score_inc in one of a haplotypes (current)
        # 
        # if_align_failed failure_reasons (in current haplotype)
        # 
        # genotype from vcf: 0/1 or 1/1 

        g.write(str(counter) + "\t")

        g.write(str(alignment_beforeSV) + "\t")
        #g.write(str(1) + "\t")
        g.write(str(alignment_afterSV) + "\t")
        #g.write(str(1) + "\t")

        #g.write(str(len(query_start_dic)) + "\t")
        g.write(str(1) + "\t")

        g.write(str(sv_type) + "\t")
        g.write(str(hap) + "\t")
        g.write(str(len(query_frag)) + "\t")
        g.write(str(len(ref_frag)) + "\t")
        g.write(str(sv_len) + "\t")
        #g.write(str(sv_end-sv_pos+1) + "\t")
        g.write(str(ref_name) + "\t")
        g.write(str(sv_pos) + "\t")
        g.write(str(sv_end) + "\t")
        g.write(str(ref_start) + "\t")
        g.write(str(ref_end) + "\t")
        g.write(str(sv_type) + "\t")
        g.write("noerr" + "\t")
        g.write("\n")
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
                    "chr21", "chr22", "chrX", "chrY"]
    else:
        chr_list = ["1", "2", "3", "4", "5",
                    "6", "7", "8", "9", "10",
                    "11", "12", "13", "14", "15",
                    "16", "17", "18", "19", "20",
                    "21", "22", "X", "Y"]
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
     
    #build map and get validation info haplotype 1
    contig_name_list, contig_pos_list, contig_name_dict = build_map(chr_len, interval, liftover_file1, if_hg38)
    get_vali_info(output_dir, vcf_file, query_file1, 1, ref_file, interval, 
                  contig_name_list, contig_pos_list, contig_name_dict, memory_limit, if_hg38, chr_list)
    
    #build map and get validation info haplotype 2
    contig_name_list, contig_pos_list, contig_name_dict = build_map(chr_len, interval, liftover_file2, if_hg38)
    get_vali_info(output_dir, vcf_file, query_file2, 2, ref_file, interval, 
                  contig_name_list, contig_pos_list, contig_name_dict, memory_limit, if_hg38, chr_list)
    

#main function and pack all the steps by both haplotypes
if __name__ == "__main__":
    main()