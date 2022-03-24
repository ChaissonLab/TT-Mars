#Get confidence scores of assembly intervals
import pysam
import sys 
import os
import csv
from Bio import SeqIO

class indexed_read:
    # constructor 
    def __init__(self, read):
        self.ori_read = read
        self.read_cur_ref_pos = read.reference_start
        self.cur_cigar_tuple = 0

def get_interval_cigar(cigar_tuples, start_pos, end_pos, read, consume_ref):
    #t
    #read_ref_start <= start_pos (100)
    #read_ref_end >= end_pos (200) - 1
    #cigar covering the interval
    read_ref_start = read.read_cur_ref_pos
    interval_cigar = []
    #first cigar tup entering the interval
    first_tup = read.cur_cigar_tuple
    #cur_ref_pos: the end of previous stop +1
    #print(read_ref_start, first_tup)
    cur_ref_pos = read_ref_start
    for i in range(read.cur_cigar_tuple, len(cigar_tuples)):
        tup = cigar_tuples[i]
        if tup[0] in consume_ref:
            cur_ref_pos += tup[1]
            #if entering the interval
            if cur_ref_pos > start_pos:
                first_tup = i
                break
    #get the interval_cigar
    if cur_ref_pos >= end_pos:
        interval_cigar.append((tup[0], end_pos - start_pos))
        read.read_cur_ref_pos = cur_ref_pos - tup[1]
        read.cur_cigar_tuple = first_tup
    else:
        interval_cigar.append((tup[0], cur_ref_pos - start_pos))
        for i in range(first_tup + 1, len(cigar_tuples)):
            tup = cigar_tuples[i]
            if tup[0] not in consume_ref:
                interval_cigar.append(tup)
            else:
                cur_ref_pos += tup[1]
                if cur_ref_pos < end_pos:
                    interval_cigar.append(tup)
                else:
                    interval_cigar.append((tup[0], tup[1] - (cur_ref_pos - end_pos)))
                    read.read_cur_ref_pos = cur_ref_pos - tup[1]
                    read.cur_cigar_tuple = i                    
                    break
    #note interval_cigar length may not sum to the interval length: including consuming query only
    return interval_cigar

#return percentage of bad bases
def bad_base_prct(interval_cigar, interval_len):
    bad_bases = 0
    for tup in interval_cigar:
        #if tup[0] not in [0, 7]:
        if tup[0] != 7:
            bad_bases += tup[1]
    return bad_bases/interval_len

#get contig seq
def getSeqRec(seq_name, file_name):
    #fasta_file = SeqIO.parse(file_name, "fasta")
    #for record in fasta_file:
    #    if seq_name == record.id.strip():
    #       return record
    fasta_file = pysam.FastaFile(file_name)
    seq = fasta_file.fetch(seq_name)
    return seq

#write score
def write_score(contig_name, start_pos, end_pos, score, g):
    g.write(str(contig_name) + "\t")
    g.write(str(start_pos) + "\t")
    g.write(str(end_pos) + "\t")
    g.write(str(score))   
    g.write("\n") 
    
#parse reads
def parse_reads(output_file_name, contig_list, interval_len, consume_ref, cut_off, samfile, min_depth, max_depth):
    # #test
    # contig_ctr = 0
    g = open(output_file_name, "w")
    for contig in contig_list:
        chr_name = contig[0]
    #     #test
    #     contig_ctr += 1
    #     print(contig_ctr)
        print(chr_name)
        contig_len = contig[1]
        #test
        #stop when parsed no_of_reads
        #no_of_reads = 400000
        #counter_read = 0
        #interval is left-closed, starting from 0
        cur_start_pos = 0
        cur_end_pos = cur_start_pos + interval_len

        ol_reads = []
        read_interval_cigars = []

        #TODO: check empty iteration?
        #print(contig_len)
        #loop through iter, index will not be reset
        #iter = samfile.fetch(chr_name, 187700, 188900)
    #     iter = samfile.fetch(chr_name, 1000, 1100)
    #     0 23 11 0.48
        iter = samfile.fetch(chr_name)

    #     for rec in iter:
    #         ol_reads.append(rec)
    #         break

        while cur_end_pos < contig_len:
            #add reads to ol_reads
            #the last one added can be out of range

            if len(ol_reads) == 0:
                try:
                    # get the next item
                    rec = next(iter)
                    ol_reads.append(indexed_read(rec))
                    # do something with element
                except StopIteration:
                    # if StopIteration is raised, break from loop
                    break     
    #             for rec in iter:
    #                 ol_reads.append(rec)
    #                 break            

            if ol_reads[len(ol_reads)-1].ori_read.reference_start < cur_end_pos:
                for rec in iter:
                    ol_reads.append(indexed_read(rec))
                    # if not overlapping, stop adding
                    if ol_reads[len(ol_reads)-1].ori_read.reference_start >= cur_end_pos:
                        break;

    #         test_ind_read = indexed_read(ol_reads[0].ori_read)
    #         print(test_ind_read.cur_cigar_tuple)
    #         test_ind_read.cur_cigar_tuple = 5
    #         print(test_ind_read.cur_cigar_tuple)
    #         test_change_tuple(test_ind_read)
    #         print(test_ind_read.cur_cigar_tuple)
    #         break

    #         class indexed_read:
    #             # constructor 
    #             def __init__(self, read):
    #                 self.ori_read = read
    #                 self.read_cur_ref_pos = read.reference_start
    #                 self.cur_cigar_tuple = 0


            #delete bottom reads that not overlapping with current read
            #note that the non-overlapping reads may not be at the bottom: lengths are different
            #thus, get rid of the non-overlapping reads at the bottom first
            #then check overlapping (with indels) every time before indels detection
            while ol_reads[0].ori_read.reference_end < cur_start_pos:
                ol_reads.pop(0)
                if len(ol_reads) == 0:
                    #print(chr_name, cur_start_pos)
                    break

            #skip if depth 0 or too large
            if len(ol_reads) == 0 or len(ol_reads) > max_depth:
                #TODO: output bad score 0 for these intervals
                write_score(chr_name, cur_start_pos, cur_end_pos, 0, g)
                cur_start_pos += interval_len
                cur_end_pos = cur_start_pos + interval_len    
                continue

            #get the reads in iter
    #         for rec in iter:
    #             if rec.reference_start > cur_start_pos or rec.reference_end < cur_end_pos:
    #                 continue
    #             ol_reads.append(rec)
            valid_ol_reads_ctr = 0

            for cur_read in ol_reads:
                #test
                #print(cur_read.reference_start, cur_read.reference_end)
                #if cur_read.reference_end < cur_start_pos or cur_read.reference_start >= cur_end_pos:
                if cur_read.ori_read.reference_end < (cur_end_pos - 1) or cur_read.ori_read.reference_start > cur_start_pos:
                    continue
                valid_ol_reads_ctr += 1
                cur_cigar_tuples = cur_read.ori_read.cigartuples
                #test
                #print(cur_cigar_tuples)
                interval_cigar = get_interval_cigar(cur_cigar_tuples, cur_start_pos, cur_end_pos, cur_read, consume_ref)
                read_interval_cigars.append(interval_cigar)
                #test
                #print(interval_cigar)

            bad_read_ctr = 0
            #Question: more elegent way to check empty in this case?
            #Coverage below threshold will be counted as 0 score intervals
            if len(read_interval_cigars) > min_depth:
                for interval_cigar in read_interval_cigars:
                    #calculate percentage of matched bases
                    #print(bad_base_prct(interval_cigar, interval_len))
                    if bad_base_prct(interval_cigar, interval_len) >= cut_off:
                        #test
                        #print(cur_start_pos)
                        bad_read_ctr += 1
                #write score
                write_score(chr_name, cur_start_pos, cur_end_pos, 1 - round(bad_read_ctr/valid_ol_reads_ctr, 2), g)
            else:
                write_score(chr_name, cur_start_pos, cur_end_pos, 0, g)

            #test
            #if bad_read_ctr > 0:
            #if round(bad_read_ctr/len(ol_reads), 2) > 0.3:    
            #    print(cur_start_pos, valid_ol_reads_ctr, bad_read_ctr, round(bad_read_ctr/valid_ol_reads_ctr, 2))
    #         if cur_start_pos % 100000 == 0:
    #             print(cur_start_pos)
    #         print(cur_start_pos, valid_ol_reads_ctr, bad_read_ctr, round(bad_read_ctr/valid_ol_reads_ctr, 2))
    #         break

            cur_start_pos += interval_len
            cur_end_pos = cur_start_pos + interval_len

            #ol_reads = []
            read_interval_cigars = []
            #if cur_end_pos >= 2000:
            #    break
    g.close()
    
# new parse reads
def parse_reads_1(output_file_name, contig_list, interval_len, consume_ref, cut_off, samfile, min_depth, max_depth, safe_len):
    # #test
    # contig_ctr = 0
    g = open(output_file_name, "w")
    for contig in contig_list:
        chr_name = contig[0]
    #     #test
    #     contig_ctr += 1
    #     print(contig_ctr)
        print(chr_name)
        contig_len = contig[1]
        #test
        #stop when parsed no_of_reads
        #no_of_reads = 400000
        #counter_read = 0
        #interval is left-closed, starting from 0
        cur_start_pos = 0
        cur_end_pos = cur_start_pos + interval_len

        ol_reads = []
        read_interval_cigars = []

        #TODO: check empty iteration?
        #print(contig_len)
        #loop through iter, index will not be reset
        #iter = samfile.fetch(chr_name, 187700, 188900)
    #     iter = samfile.fetch(chr_name, 1000, 1100)
    #     0 23 11 0.48
        iter = samfile.fetch(chr_name)

    #     for rec in iter:
    #         ol_reads.append(rec)
    #         break

        while cur_end_pos < contig_len:
            #add reads to ol_reads
            #the last one added can be out of range

            if len(ol_reads) == 0:
                try:
                    # get the next item
                    rec = next(iter)
                    ol_reads.append(indexed_read(rec))
                    # do something with element
                except StopIteration:
                    # if StopIteration is raised, break from loop
                    break     
    #             for rec in iter:
    #                 ol_reads.append(rec)
    #                 break            

            if ol_reads[len(ol_reads)-1].ori_read.reference_start < cur_end_pos:
                for rec in iter:
                    ol_reads.append(indexed_read(rec))
                    # if not overlapping, stop adding
                    if ol_reads[len(ol_reads)-1].ori_read.reference_start >= cur_end_pos:
                        break

            #delete bottom reads that not overlapping with current read
            #note that the non-overlapping reads may not be at the bottom: lengths are different
            #thus, get rid of the non-overlapping reads at the bottom first
            #then check overlapping (with indels) every time before indels detection
            while ol_reads[0].ori_read.reference_end < cur_start_pos:
                ol_reads.pop(0)
                if len(ol_reads) == 0:
                    #print(chr_name, cur_start_pos)
                    break

            #skip if depth 0 or too large
            if len(ol_reads) == 0 or len(ol_reads) > max_depth:
                #TODO: output bad score 0 for these intervals
                write_score(chr_name, cur_start_pos, cur_end_pos, 0, g)
                cur_start_pos += interval_len
                cur_end_pos = cur_start_pos + interval_len    
                continue

            #get the reads in iter
    #         for rec in iter:
    #             if rec.reference_start > cur_start_pos or rec.reference_end < cur_end_pos:
    #                 continue
    #             ol_reads.append(rec)
            valid_ol_reads_ctr = 0
            good_reads_ctr = 0

            for cur_read in ol_reads:
                #test
                #print(cur_read.reference_start, cur_read.reference_end)
                #if cur_read.reference_end < cur_start_pos or cur_read.reference_start >= cur_end_pos:
                if cur_read.ori_read.reference_end < (cur_end_pos - 1) or cur_read.ori_read.reference_start > cur_start_pos:
                    continue
                valid_ol_reads_ctr += 1
                #count reads extend the interval by at least save_len as good
                if cur_read.ori_read.reference_end >= min((cur_end_pos + safe_len), contig_len-1) and \
                    cur_read.ori_read.reference_start <= max((cur_start_pos - safe_len), 0):
                    good_reads_ctr += 1

    #         bad_read_ctr = 0
    #         #Question: more elegent way to check empty in this case?
    #         #Coverage below threshold will be counted as 0 score intervals
            if valid_ol_reads_ctr > min_depth:
                #write score
                write_score(chr_name, cur_start_pos, cur_end_pos, round(good_reads_ctr/valid_ol_reads_ctr, 2), g)
            else:
                write_score(chr_name, cur_start_pos, cur_end_pos, 0, g)

            cur_start_pos += interval_len
            cur_end_pos = cur_start_pos + interval_len

            #ol_reads = []
            read_interval_cigars = []
    g.close()
    
# new parse reads
def parse_reads_0831(output_file_name, contig_list, interval_len, consume_ref, cut_off, samfile, min_depth, max_depth, safe_len):
    # #test
    # contig_ctr = 0
    g = open(output_file_name, "w")
    for contig in contig_list:
        contig_name = contig[0]
        print(contig_name)
        contig_len = contig[1]
        
        #interval is left-closed, starting from 0: eg: [0, 100)
        cur_start_pos = 0
        cur_end_pos = cur_start_pos + interval_len

        ol_reads = []

        #loop through iter, index will not be reset
        #iter = samfile.fetch(contig_name, 187700, 188900)
    #     iter = samfile.fetch(contig_name, 1000, 1100)
        iter = samfile.fetch(contig_name)

        while cur_end_pos < contig_len:
            #add reads to ol_reads
            #the last one added can be out of current interval
            if len(ol_reads) == 0:
                try:
                    # get the next item
                    rec = next(iter)
                    ol_reads.append(indexed_read(rec))
                    # do something with element
                except StopIteration:
                    # if StopIteration is raised, break from loop
                    break     
                    
    #             for rec in iter:
    #                 ol_reads.append(rec)
    #                 break            

            if ol_reads[len(ol_reads)-1].ori_read.reference_start < cur_end_pos:
                for rec in iter:
                    ol_reads.append(indexed_read(rec))
                    # if not overlapping, stop adding
                    if ol_reads[len(ol_reads)-1].ori_read.reference_start >= cur_end_pos:
                        break

            #delete bottom reads that not overlapping with current read
            #note that the non-overlapping reads may not be at the bottom: lengths are different
            #thus, get rid of the non-overlapping reads at the bottom first
            #then check overlapping (with indels) every time before indels detection
            while ol_reads[0].ori_read.reference_end < cur_start_pos:
                ol_reads.pop(0)
                if len(ol_reads) == 0:
                    #print(contig_name, cur_start_pos)
                    break

            #skip if depth 0
            if len(ol_reads) == 0:
                write_score(contig_name, cur_start_pos, cur_end_pos, 0, g)
                cur_start_pos += interval_len
                cur_end_pos = cur_start_pos + interval_len    
                continue
                
            valid_ol_reads_ctr = 0
            good_reads_ctr = 0

            for cur_read in ol_reads:
                #test
                #print(cur_read.reference_start, cur_read.reference_end)
                #if cur_read.reference_end < cur_start_pos or cur_read.reference_start >= cur_end_pos:
                #make sure read has to cover the interval
                if cur_read.ori_read.reference_end < (cur_end_pos - 1) or cur_read.ori_read.reference_start > cur_start_pos:
                    continue
                valid_ol_reads_ctr += 1
                #count reads extend the interval by at least save_len as good
                if cur_read.ori_read.reference_end >= min((cur_end_pos + safe_len), contig_len-1) and \
                    cur_read.ori_read.reference_start <= max((cur_start_pos - safe_len), 0):
                    good_reads_ctr += 1

            write_score(contig_name, cur_start_pos, cur_end_pos, good_reads_ctr, g)

            cur_start_pos += interval_len
            cur_end_pos = cur_start_pos + interval_len
    g.close()
    
#main function
def main():
    #get command line input
    #n = len(sys.argv)
    output_dir = sys.argv[1] + "/"
    #read to assembly bam files
    read2assem1_bamfile = sys.argv[2]
    read2assem2_bamfile = sys.argv[3]
    centromere_file = sys.argv[4]
    #assembly fasta files
    assem1_fasta = sys.argv[5]
    assem2_fasta = sys.argv[6]
    #output names
    output_name1 = sys.argv[7]
    output_name2 = sys.argv[8]
    
    interval_len = 100
    safe_len = 1000
    #if not matches bases less than cut_off, the read is good in currect interval
    cut_off = 0.1
    #opeartion index that consume reference
    consume_ref = [0,2,3,7,8]
    #operation index that consume query
    consume_query = [0,1,4,7,8]    
    #min valid depth for an interval, below will have score 0
    min_depth = 15
    max_depth = 2000
    
    #hap1
    output_file = output_dir + output_name1
    samfile = pysam.AlignmentFile(read2assem1_bamfile, "rb")
    contig_list = []
    for seq_record in SeqIO.parse(assem1_fasta, "fasta"):
        #print(seq_record.id)
        #print(repr(seq_record.seq))
        #print(len(seq_record))
        contig_list.append([seq_record.id, len(seq_record)])
    #parse reads and get scores of assembly intervals and write results
    parse_reads_0831(output_file, contig_list, interval_len, consume_ref, cut_off, samfile, min_depth, max_depth, safe_len)
    
    #hap2
    output_file = output_dir + output_name2
    samfile = pysam.AlignmentFile(read2assem2_bamfile, "rb")
    contig_list = []
    for seq_record in SeqIO.parse(assem2_fasta, "fasta"):
        #print(seq_record.id)
        #print(repr(seq_record.seq))
        #print(len(seq_record))
        contig_list.append([seq_record.id, len(seq_record)])
    #parse reads and get scores of assembly intervals and write results
    parse_reads_0831(output_file, contig_list, interval_len, consume_ref, cut_off, samfile, min_depth, max_depth, safe_len)    
    
if __name__ == "__main__":
    main()