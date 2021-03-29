#Get confidence scores of assembly intervals
import pysam
import sys 
import os
import csv
from Bio import SeqIO

def get_interval_cigar(cigar_tuples, start_pos, end_pos, read_ref_start, consume_ref):
    #read_ref_start <= start_pos (100)
    #read_ref_end >= end_pos (200) - 1
    #cigar covering the interval
    interval_cigar = []
    #first cigar tup entering the interval
    first_tup = -1
    #cur_ref_pos: the end of previous stop +1
    cur_ref_pos = read_ref_start
    for i in range(0, len(cigar_tuples)):
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
def parse_reads(output_file_name, contig_list, interval_len, consume_ref, cut_off, samfile):
    # #test
    # contig_ctr = 0
    g = open(output_file_name, "w")
    for contig in contig_list:
        chr_name = contig[0]
    #     #test
    #     contig_ctr += 1
    #     print(contig_ctr)
    #     print(chr_name)
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
                    ol_reads.append(rec)
                    # do something with element
                except StopIteration:
                    # if StopIteration is raised, break from loop
                    break     
    #             for rec in iter:
    #                 ol_reads.append(rec)
    #                 break            

            if ol_reads[len(ol_reads)-1].reference_start < cur_end_pos:
                for rec in iter:
                    ol_reads.append(rec)
                    # if not overlapping, stop adding
                    if ol_reads[len(ol_reads)-1].reference_start >= cur_end_pos:
                        break;

            #delete bottom reads that not overlapping with current read
            #note that the non-overlapping reads may not be at the bottom: lengths are different
            #thus, get rid of the non-overlapping reads at the bottom first
            #then check overlapping (with indels) every time before indels detection
            while ol_reads[0].reference_end < cur_start_pos:
                ol_reads.pop(0)
                if len(ol_reads) == 0:
                    #print(chr_name, cur_start_pos)
                    break

            if len(ol_reads) == 0:
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
                if cur_read.reference_end < cur_start_pos or cur_read.reference_start >= cur_end_pos:
                    continue
                valid_ol_reads_ctr += 1
                cur_cigar_tuples = cur_read.cigartuples
                #test
                #print(cur_cigar_tuples)
                interval_cigar = get_interval_cigar(cur_cigar_tuples, cur_start_pos, cur_end_pos, cur_read.reference_start, consume_ref)
                read_interval_cigars.append(interval_cigar)
                #test
                #print(interval_cigar)

            bad_read_ctr = 0
            #Question: more elegent way to check empty in this case?
            if len(read_interval_cigars) > 0:
                for interval_cigar in read_interval_cigars:
                    #calculate percentage of matched bases
                    #print(bad_base_prct(interval_cigar, interval_len))
                    if bad_base_prct(interval_cigar, interval_len) >= cut_off:
                        #test
                        #print(cur_start_pos)
                        bad_read_ctr += 1
            #write score
            write_score(chr_name, cur_start_pos, cur_end_pos, 1 - round(bad_read_ctr/valid_ol_reads_ctr, 2), g)

            #test
            #if bad_read_ctr > 0:
            #if round(bad_read_ctr/len(ol_reads), 2) > 0.3:    
            #    print(cur_start_pos, valid_ol_reads_ctr, bad_read_ctr, round(bad_read_ctr/valid_ol_reads_ctr, 2))
            #if cur_start_pos % 100000 == 0:
            #    print(cur_start_pos)

            cur_start_pos += interval_len
            cur_end_pos = cur_start_pos + interval_len

            #ol_reads = []
            read_interval_cigars = []
            #if cur_end_pos >= 2000:
            #    break
    g.close()
    
#main function
def main():
    #get command line input
    #n = len(sys.argv)
    output_dir = sys.argv[1] + "/"
    #read to assembly bam files
    #bam_file = "/panfs/qcb-panasas/jianzhiy/data/pacbio/HG002/h1_toassem_nosec_0326.bam"
    #bam_file = "/panfs/qcb-panasas/jianzhiy/data/pacbio/HG002/h2_toassem_nosec_0326.bam"
    read2assem1_bamfile = sys.argv[2]
    read2assem2_bamfile = sys.argv[3]
    centromere_file = sys.argv[4]
    #assembly fasta files
    #assem1_fasta = "/home/cmb-16/nobackups/mjc/data_download/HG002/CCS/assembly.1.consensus.fasta"
    assem1_fasta = sys.argv[5]
    assem2_fasta = sys.argv[6]
    #output names
    output_name1 = sys.argv[7]
    output_name2 = sys.argv[8]
    
    interval_len = 100
    #save_len = 2000
    #if ot matches bases less than cut_off, the read is good in currect interval
    cut_off = 0.1
    #opeartion index that consume reference
    consume_ref = [0,2,3,7,8]
    #operation index that consume query
    consume_query = [0,1,4,7,8]    
    
    #hap1
    output_file = output_dir + output_name1
    #centromere_file = "/home/cmb-16/mjc/quentin/illuVeri/use_mcutils/data_files/centromere_hg37.txt"
    samfile = pysam.AlignmentFile(read2assem1_bamfile, "rb")
    contig_list = []
    for seq_record in SeqIO.parse(assem1_fasta, "fasta"):
        #print(seq_record.id)
        #print(repr(seq_record.seq))
        #print(len(seq_record))
        contig_list.append([seq_record.id, len(seq_record)])
    #parse reads and get scores of assembly intervals and write results
    parse_reads(output_file, contig_list, interval_len, consume_ref, cut_off, samfile)
    
    #hap2
    output_file = output_dir + output_name2
    #centromere_file = "/home/cmb-16/mjc/quentin/illuVeri/use_mcutils/data_files/centromere_hg37.txt"
    samfile = pysam.AlignmentFile(read2assem2_bamfile, "rb")
    contig_list = []
    for seq_record in SeqIO.parse(assem2_fasta, "fasta"):
        #print(seq_record.id)
        #print(repr(seq_record.seq))
        #print(len(seq_record))
        contig_list.append([seq_record.id, len(seq_record)])
    #parse reads and get scores of assembly intervals and write results
    parse_reads(output_file, contig_list, interval_len, consume_ref, cut_off, samfile)    
    
if __name__ == "__main__":
    main()