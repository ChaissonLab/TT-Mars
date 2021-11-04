# %%
#bedtools filter out SVs that are too close to end of contigs

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

