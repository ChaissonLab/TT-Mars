#get chr length from ref fasta

import pysam
import sys 

#if_hg38_str = True
#ref_file = "./hg38.no_alts.fasta"
    
ref_file = sys.argv[1]
if_hg38_str = sys.argv[2]

if if_hg38_str == "True":
    if_hg38 = True
else:
    if_hg38 = False

chr_len = []
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

ref_fasta_file = pysam.FastaFile(ref_file)
for chr_name in chr_list:
    ref_rec = ref_fasta_file.fetch(chr_name)
    chr_len.append(len(ref_rec))
print(chr_len)