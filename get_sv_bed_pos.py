# %%
#bedtools filter out SVs that are too close to end of contigs

#0817
#TODO
#1. Change append to write;

#Generate plots
import csv
import math
#pysam: https://pysam.readthedocs.io/en/latest/usage.html
import pysam
#print(pysam.__version__)

import numpy as np 
import matplotlib.pyplot as plt  

import sys 
  
#get command line input
#n = len(sys.argv)
output_dir = sys.argv[1]
vcf_file = sys.argv[2]

#vcf_file = "../../../callset_files/naive_caller/0810_v1_result/lra_ccs_0810_v1_sorted.vcf.gz"
#test
#vcf_file = "../../callset_files/naive_caller/test_example.vcf"
#output_dir = "../../output/v4.0/naive_caller_0810_v1/"

# %%
#process liftover results
'''
with open(output_dir + "contig_positions_assem1_draft_result.bed") as f:
    reader = csv.reader(f, delimiter="\t")
    contig_positions_assem1_result = list(reader)
f.close()

for record in contig_positions_assem1_result:
    name = record[0]
    start = record[1]
    end = record[2]
    if int(start) >= int(end):
        continue
    g = open(output_dir + "contig_positions_assem1.bed", "a")
    g.write(str(name) + "\t")
    g.write(str(start) + "\t")
    g.write(str(end))
    g.write("\n")
    g.close()

with open(output_dir + "contig_positions_assem2_draft_result.bed") as f:
    reader = csv.reader(f, delimiter="\t")
    contig_positions_assem2_result = list(reader)
f.close()

for record in contig_positions_assem2_result:
    name = record[0]
    start = record[1]
    end = record[2]
    if int(start) >= int(end):
        continue
    g = open(output_dir + "contig_positions_assem2.bed", "a")
    g.write(str(name) + "\t")
    g.write(str(start) + "\t")
    g.write(str(end))
    g.write("\n")
    g.close()

'''

# %%
# %%
#bedtools filter out SVs that are too close to each other

#create bedfiles

f = pysam.VariantFile(vcf_file, 'r')

#create bedfile contains SVs' positions
g = open(output_dir + "SV_positions.bed", "w")
for counter, rec in enumerate(f.fetch()):
	ref_name = rec.chrom
	sv_type = rec.info['SVTYPE']
	sv_len = rec.rlen
	#TODOL double check the start for different types
	sv_pos = rec.pos
	sv_end = rec.stop

	g.write(str(ref_name) + "\t")
	g.write(str(sv_pos) + "\t")
	g.write(str(sv_end))
	g.write("\n")
g.close()



#using bedtools to filter SV
#bedtools intersect -a ../giab_data_output/SV_positions.bed -b ../../run_manta_data/manta_output/contig_positions_assem1_result_new.bed -u > ../giab_data_output/exclude__assem1_result.bed
#bedtools intersect -a ../giab_data_output/SV_positions.bed -b ../../run_manta_data/manta_output/contig_positions_assem2_result_new.bed -u > ../giab_data_output/exclude__assem2_result.bed