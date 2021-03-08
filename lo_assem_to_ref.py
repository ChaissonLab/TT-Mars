# %%
#generate bed file that contains positions on ref to be lifted to contigs
import csv
import math
#pysam: https://pysam.readthedocs.io/en/latest/usage.html
import pysam
#print(pysam.__version__)
import sys 
  
#get command line input
#n = len(sys.argv)
output_dir = sys.argv[1] + "/"
assembly_bam_file_hap1 = sys.argv[2]
assembly_bam_file_hap2 = sys.argv[3]
#print("\nName of Python script:", sys.argv[0])

interval = 20

#check if k-th bit of a given number is set or not 
def isKthBitSet(n, k): 
    if n & (1 << (k - 1)): 
        return True 
    else: 
        return False

#for assem1
#generate bed file liftover

#TODO: change output 
g = open(output_dir + "lo_pos_assem1_withSed.bed", "w")

samfile = pysam.AlignmentFile(assembly_bam_file_hap1, "rb")
for record in samfile:
    #TODO: solve the unmapped segment problem???
    if not isKthBitSet(record.flag, 3):
        #samLiftover's bed file needs ref name for both directions
        ref_name = record.reference_name
        query_name = record.query_name
        #TODO: check not matching with paf
        query_start = record.query_alignment_start
        query_end = record.query_alignment_end
        ref_start = record.reference_start
        ref_end = record.reference_end
        query_start = math.ceil(query_start/interval) * interval
        query_end = math.floor(query_end/interval) * interval
        ref_start = math.ceil(ref_start/interval) * interval
        ref_end = math.floor(ref_end/interval) * interval

        for i in range(ref_start, ref_end + interval, interval):
            g.write(ref_name + "\t")
            g.write(str(i) + "\t")
            g.write(str(i+1) + "\t")
            g.write(query_name + "\t")
            g.write("\n")

g.close()

#for assem2
g = open(output_dir + "lo_pos_assem2_withSed.bed", "w")

samfile = pysam.AlignmentFile(assembly_bam_file_hap2, "rb")
for record in samfile:
    #TODO: solve the unmapped segment problem???
    if not isKthBitSet(record.flag, 3):
        #samLiftover's bed file needs ref name for both directions
        ref_name = record.reference_name
        query_name = record.query_name
        #TODO: check not matching with paf
        query_start = record.query_alignment_start
        query_end = record.query_alignment_end
        ref_start = record.reference_start
        ref_end = record.reference_end
        query_start = math.ceil(query_start/interval) * interval
        query_end = math.floor(query_end/interval) * interval
        ref_start = math.ceil(ref_start/interval) * interval
        ref_end = math.floor(ref_end/interval) * interval

        for i in range(ref_start, ref_end + interval, interval):
            g.write(ref_name + "\t")
            g.write(str(i) + "\t")
            g.write(str(i+1) + "\t")
            g.write(query_name + "\t")
            g.write("\n")

g.close()
