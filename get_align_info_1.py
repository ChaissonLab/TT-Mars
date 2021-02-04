# %%
#Results: 
# SV_index SV_type
# score_before_sv_assem1 score_after_sv_assem1 score_before_sv_assem2 score_after_sv_assem2
# contig_seg_len ref_seg_len SV_len
# chr_name SV_start SV_end interval_start interval_end
# if_score_inc in one of a haplotypes
# if_align_failed failure_reasons
# genotype from vcf: 0/1 or 1/1

#not produce now:
# comb_score_before comb_score_after better_assem (1 or 2)


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
  
#get command line input
#n = len(sys.argv)
output_dir = sys.argv[1] + "/"
vcf_file = sys.argv[2]
bam_file1 = sys.argv[3]
bam_file2 = sys.argv[4]
ref_file = sys.argv[5]
query_file1 = sys.argv[6]
query_file2 = sys.argv[7]

# %%
#constants

interval = 20

# %%
#functions

#Will be modified

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
def getRefStart(sv_pos):
    start = math.floor(sv_pos/interval) * interval
    return start

#get the most close 'interval(500)' after sv_end
def getRefEnd(sv_end):
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

'''
#infer the corresponding pos on query of a pos on ref
def infer_sv_query_pos(read, pos):
    #query_sequence: read sequence bases, including soft clipped bases; use this one!!!
    #query_alignment_start: This the index of the first base in seq that is not soft-clipped
    #query_alignment_sequence: This is a substring of seq that excludes flanking bases that were soft clipped
    ref_align_start = read.reference_start
    ref_align_end = read.reference_end
    #pos on query after the SC
    cur_query_pos = read.query_alignment_start
    cur_ref_pos = ref_align_start
    
    if pos < ref_align_start or pos > ref_align_end:
        return -1
    
    #opeartion index that consume reference
    consume_ref = [0,2,3,7,8]
    #operation index that consume query
    consume_query = [0,1,4,7,8]
    
    cigar = read.cigartuples
    #test
    print(pos)
    for tup in cigar:
        #test
        print(tup)
        #cigar string can contain hardclip
        if tup[0] in [5, 6]:
            continue
        #if softclip
        if tup[0] == 4:
            cur_query_pos = tup[1]
            continue
            
        if ((tup[0] in [0,7,8]) and cur_ref_pos <= pos and cur_ref_pos + tup[1] > pos):
            return cur_query_pos + (pos - cur_ref_pos)
        elif ((tup[0] in [2, 3]) and cur_ref_pos <= pos and cur_ref_pos + tup[1] > pos):
            return cur_query_pos
        
        if tup[0] in consume_ref:
            cur_ref_pos += tup[1]
        if tup[0] in consume_query:
            cur_query_pos += tup[1]
    
    return -1

'''

#Modify the dictionay to reduce memory usage
dict1 = dict()
with open(output_dir + "lo_pos_assem1_withSed_result.bed") as f:
    for line in f:
        record = line.strip().split()
        #first level key: refname
        first_key = record[4]
        #int is used as str in dict
        #second level key: refpos
        second_key = record[5]
        #third level key: queryname
        third_key = record[0]
        if first_key in dict1:
            #if second_key in dict1[first_key]:
                #TODO: Note: do NOT allow one ref pos have multiple pos on query contigs now!!!
                #dict1[first_key][second_key][third_key] = record[1]
            #else:
            if second_key not in dict1[first_key]:
                #third_dict = {third_key: record[1]}
                third_dict = [third_key, record[1]]
                dict1[first_key][second_key] = third_dict
        else:
            #third_dict = {third_key: record[1]}
            third_dict = [third_key, record[1]]
            second_dict = {second_key: third_dict}
            dict1[first_key] = second_dict
f.close()

'''
dict2 = dict()
with open(output_dir + "lo_pos_assem2_withSed_result.bed") as f:
    for line in f:
        record = line.strip().split()
        #first level key: refname
        first_key = record[4]
        #int is used as str in dict
        #second level key: refpos
        second_key = record[5]
        #third level key: queryname
        third_key = record[0]
        if first_key in dict2:
            #if second_key in dict2[first_key]:
                #TODO: Note: do NOT allow one ref pos have multiple pos on query contigs now!!!
                #dict2[first_key][second_key][third_key] = record[1]
            #else:
            if second_key not in dict2[first_key]:
                #third_dict = {third_key: record[1]}
                third_dict = [third_key, record[1]]
                dict2[first_key][second_key] = third_dict
        else:
            #third_dict = {third_key: record[1]}
            third_dict = [third_key, record[1]]
            second_dict = {second_key: third_dict}
            dict2[first_key] = second_dict
f.close()
'''

# %%
#get vcf file and run score_callset on each SV record

#build a dictionary:
#key: position on hg37
#value: lifted over position on hg38
#files are resutls from lo_ref.py

#TODO: now the liftover result is for DEL only

#build a list to store positions of DEL cases failed to be lifted over
'''
DEL_fail_pos = []
with open("../DEL_calls_hg38_pos_from_hg37_fail.bed") as f:
    reader = csv.reader(f, delimiter="\t")
    DEL_fail_pos_raw = list(reader)
f.close()

for record in DEL_fail_pos_raw:
	#if the line is not failure information
	if len(record) > 1:
		#appened a list of str
		DEL_fail_pos.append(record)

#build a dictionary for the cases not in DEL_fail_pos

with open("../DEL_calls_hg37_pos.bed") as f:
    reader = csv.reader(f, delimiter="\t")
    DEL_calls_hg37_pos = list(reader)
f.close()

with open("../DEL_calls_hg38_pos_from_hg37.bed") as f:
    reader = csv.reader(f, delimiter="\t")
    DEL_calls_hg38_pos_from_hg37 = list(reader)
f.close()

dict_hg37_to_hg38 = dict()

counter_success = 0
for record in DEL_calls_hg37_pos:
	temp = [record[0], record[1], record[2]]
	if temp in DEL_fail_pos:
		continue
	else:
		#key: chrname_startpos_endpos
		key = str(record[0]) + str(record[1]) + str(record[2])
		dict_hg37_to_hg38[key] = DEL_calls_hg38_pos_from_hg37[counter_success]
		counter_success = counter_success + 1
'''


# %%
#get vcf file and run score_callset on each SV record

#TODO: filter out centromere cases

f = pysam.VariantFile(vcf_file,'r')

#dict1/dict2
dict_on = dict1
query_file = query_file1
hap = 1
name_str = "assem1"
sam_file = bam_file1

chromosome = "all"

output_file_name = output_dir + "align_info_" + name_str + "_chr" + chromosome + ".txt"


query_fasta_file = pysam.FastaFile(query_file)
ref_fasta_file = pysam.FastaFile(ref_file)
ref_name = ""

for counter, rec in enumerate(f.fetch()):
    #counter of all DEL

    #test 1129 3636
    #if counter < 3636:
    #    continue

    #if counter > 110:
    #    break
        
    name = rec.chrom

    #choose target chr
    #if str(name) != chromosome:
    #    continue

    sv_type = rec.info['SVTYPE']
    
    if sv_type != 'DEL' and sv_type != 'INS':
        print("Wrong type!")
        continue
    #sv_len = abs(rec.info["SVLEN"][0])
    #without abs value: rela length closer to 1 the better
    sv_len = rec.info['SVLEN'][0]
    #TODOL double check the start for different types
    sv_pos = rec.pos
    sv_end = rec.stop
    #if counter == 1271:
    #    break
    
    sv_ref_seq = rec.ref
    #TODO: update here when multi-calls in one records allowed
    sv_query_seq = rec.alts[0]

    #TODO: use switch instead of if/else
    #find given fragment (on ref) to queries
    
    #get ref seq
    if ref_name != name:
        ref_name = name
        ref_rec = ref_fasta_file.fetch(ref_name)

    ref_start = getRefStart(sv_pos)
    ref_end = getRefEnd(sv_end)
    #don't know query name
    #query_name

    #first level key: refname
    first_key = ref_name
    #int is used as str in dict
    #second level key: refpos

    
    if counter % 1000 == 0:
        print(counter)
    
    #Search for the best second_key_start and second_key_end in the regions
    region_len_m = 500
    min_diff_len = 100000000
    ref_start_best = ref_start
    ref_end_best = ref_end
    query_start_best = None
    query_end_best = None
        
    neg_strand = False
    
    for ref_end_cand in range(ref_end, ref_end+region_len_m, 20):
        for ref_start_cand in range(ref_start, ref_start-region_len_m, -20):
            second_key_start = str(ref_start_cand)
            second_key_end = str(ref_end_cand)
                
            first_res = lookupDict(first_key, dict_on)
            if first_res:
                second_start_res = lookupDict(second_key_start, first_res)
                if second_start_res:
                    query_start_dic = second_start_res
                else:
                    print("wrong second start key")
                    message = "wrong_sec_start_key"
                    #write_err(output_file_name, message)
                    continue

                second_end_res = lookupDict(second_key_end, first_res)
                if second_end_res:
                    query_end_dic = second_end_res
                else:
                    print("wrong second end key")
                    message = "wrong_sec_end_key"
                    #write_err(output_file_name, message)
                    continue
            else:
                print("wrong first key")
                message = "wrong_first_key"
                #write_err(output_file_name, message)
                continue 

            #print(len(query_start_dic))

            #query_name = next(iter(query_start_dic))
            query_name = query_start_dic[0]

            if query_name != query_end_dic[0]:
                print("Not same contig")
                message = "not_same_contig"
                #write_err(output_file_name, message)
                continue
            query_start = int(query_start_dic[1])
            query_end = int(query_end_dic[1])
            
            neg_strand_tep = False
            #in case: negtive strand
            if query_end < query_start:
                tem = query_start
                query_start = query_end
                query_end = tem
                neg_strand_tep = True
                
            if abs((query_end - query_start) - (ref_end_cand - ref_start_cand) - sv_len) < min_diff_len:
                min_diff_len = abs((query_end - query_start) - (ref_end_cand - ref_start_cand) - sv_len)
                ref_start_best = ref_start_cand
                ref_end_best = ref_end_cand
                query_start_best = query_start
                query_end_best = query_end
                if neg_strand_tep:
                    neg_strand = True
                else:
                    neg_strand = False
    
    if query_start_best == None or query_end_best == None:
        message = "Wrong_query_pos"
        write_err(output_file_name, message)
        continue
    
    ref_start = ref_start_best
    ref_end = ref_end_best
    query_start = query_start_best
    query_end = query_end_best
    
    if ref_start == sv_pos:
        ref_start = ref_start - 1
    
    if ref_end == sv_end:
        ref_end = ref_end + 1

    if query_end == query_start:
        query_end += 1
        #continue
        
    query_rec = query_fasta_file.fetch(query_name)
    
    if query_start >= len(query_rec) or query_end >= len(query_rec):
        message = "bad_query_pos"
        write_err(output_file_name, message)
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
        memory_limit = 20000
        if len(str(query_frag)) > memory_limit or len(str(ref_frag)) > memory_limit:
            message = "memory_limit"
            write_err(output_file_name, message)
            continue
        if len(str(ref_afterSV_frag1)) + len(str(ref_afterSV_frag2)) > memory_limit:
            message = "memory_limit"
            write_err(output_file_name, message)
            continue

        #TODO: find a appropriate alignment parameters
        #paras: match, mismatch, open gap, extend gap
        alignment_beforeSV = pairwise2.align.globalms(str(query_frag), str(ref_frag), 1, -1, -1, -0.5)
        alignment_afterSV = pairwise2.align.globalms(str(query_frag), str(ref_afterSV_frag1) 
                                + str(ref_afterSV_frag2), 1, -1, -1, -0.5)
	    #get correct query info format
        
    elif sv_type == "INS":
        #query and ref seq fragment
        query_frag = query_rec[query_start:query_end]
        ref_frag = ref_rec[ref_start:ref_end]
        #TODO: this is for DEL
        #TODO: check +-1
        ref_afterSV_frag1 = ref_rec[ref_start:sv_pos]
        ref_afterSV_frag2 = ref_rec[sv_end:ref_end]
        ins_seq = sv_query_seq

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
        ins_seq = ins_seq.upper()
        query_frag = query_frag.upper()

        #TODO: fragments too long will cause memory problem
        memory_limit = 20000
        if len(str(query_frag)) > memory_limit or len(str(ref_frag)) > memory_limit:
            message = "memory_limit"
            write_err(output_file_name, message)
            continue
        if len(str(ref_afterSV_frag1)) + len(str(ins_seq)) + len(str(ref_afterSV_frag2)) > memory_limit:
            message = "memory_limit"
            write_err(output_file_name, message)
            continue
            
        if len(str(query_frag)) == 0 or len(str(query_frag)) == 0 or len(str(ref_afterSV_frag1)) + len(str(ins_seq)) + len(str(ref_afterSV_frag2)) == 0:
            message = "Wrong seq!!!"
            write_err(output_file_name, message)
            continue

        #TODO: find a appropriate alignment parameters
        #paras: match, mismatch, open gap, extend gap
        alignment_beforeSV = pairwise2.align.globalms(str(query_frag), str(ref_frag), 1, -1, -1, -0.5)
        alignment_afterSV = pairwise2.align.globalms(str(query_frag), str(ref_afterSV_frag1) + str(ins_seq)
                                + str(ref_afterSV_frag2), 1, -1, -1, -0.5)
	    #get correct query info format


    #TODO: check the no of contig for given position!!!
    
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
    # 

    f = open(output_file_name, "a")
    f.write(str(counter) + "\t")
    f.write(str(alignment_beforeSV[0][2]) + "\t")
    #f.write(str(1) + "\t")
    f.write(str(alignment_afterSV[0][2]) + "\t")
    #f.write(str(1) + "\t")
    f.write(str(len(query_start_dic)) + "\t")
    f.write(str(sv_type) + "\t")
    f.write(str(hap) + "\t")
    f.write(str(len(query_frag)) + "\t")
    f.write(str(len(ref_frag)) + "\t")
    f.write(str(sv_len) + "\t")
    #f.write(str(sv_end-sv_pos+1) + "\t")
    f.write(str(ref_name) + "\t")
    f.write(str(sv_pos) + "\t")
    f.write(str(sv_end) + "\t")
    f.write(str(ref_start) + "\t")
    f.write(str(ref_end) + "\t")
    f.write(str(sv_type) + "\t")
    f.write("noerr" + "\t")
    f.write("\n")
    f.close()

# %%
