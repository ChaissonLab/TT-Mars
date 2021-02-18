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
  
#get command line input
#n = len(sys.argv)
output_dir = sys.argv[1] + "/"
vcf_file = sys.argv[2]
bam_file1 = sys.argv[3]
bam_file2 = sys.argv[4]
ref_file = sys.argv[5]
query_file1 = sys.argv[6]
query_file2 = sys.argv[7]
liftover_file = sys.argv[8]

# %%
#constants

#interval length
interval = 20

#approximate length of chromosomes
chr_len = [250000000, 244000000, 199000000, 192000000, 182000000, 
            172000000, 160000000, 147000000, 142000000, 136000000, 
            136000000, 134000000, 116000000, 108000000, 103000000, 
            90400000, 83300000, 80400000, 59200000, 64500000, 
            48200000, 51400000, 157000000, 59400000]

# %%
#functions

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
    
#return chr name as an integer
#X = 23, Y = 24
def get_int_chr_name(name):
    #remove "chr"
    if if_hg38:
        chr_name = name[3:]
    
    if chr_name.upper() == 'X':
        int_name = 23
    elif chr_name.upper() == 'Y':
        int_name = 24
    else:
        int_name = int(chr_name)
        
    return int_name
    

# %%

#mapping
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
        int_ref_name = get_int_chr_name(record[4])
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


