
# %%
#In this file: Liftover hg38 reference to HG002 assembly contigs using samLiftover
#liftover results in ../../data_files/hg19/lo_pos_assem1_withSed_result.bed and ../../data_files/hg19/lo_pos_assem2_withSed_result.bed

# %%
#minimap2 results: no secondary alignments

#minimap2 -x asm20 -ac --secondary=no /panfs/qcb-panasas/jingwenr/reference/CLR/human_hs37d5.fasta /home/cmb-16/nobackups/mjc/data_download/HG002/CCS/assembly.1.consensus.fasta > ../../data_files/hg19/HG002_to_hg37_asm20_woSed_assem1.sam
#minimap2 -x asm20 -ac --secondary=no /panfs/qcb-panasas/jingwenr/reference/CLR/human_hs37d5.fasta /home/cmb-16/nobackups/mjc/data_download/HG002/CCS/assembly.2.consensus.fasta > ../../data_files/hg19/HG002_to_hg37_asm20_woSed_assem2.sam


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
output_dir = sys.argv[1]
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

#TODO: change output dir
#g = open("../../data_files/hg19/lo_pos_assem1_withSed.bed", "w")
g = open(output_dir + "lo_pos_assem1_withSed.bed", "w")

#samfile = pysam.AlignmentFile("../../data_files/hg19//HG002_to_hg37_asm20_woSed_assem1.bam", "rb")
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
#g = open("../../data_files/hg19/lo_pos_assem2_withSed.bed", "w")
g = open(output_dir + "lo_pos_assem2_withSed.bed", "w")

#samfile = pysam.AlignmentFile("../../data_files/hg19/HG002_to_hg37_asm20_woSed_assem2.bam", "rb")
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

# %%
#run samLiftover

#/home/cmb-16/mjc/quentin/illuVeri/use_mcutils/mcutils/bin/samLiftover ../../data_files/hg19/HG002_to_hg37_asm20_woSed_assem1.sam ../../data_files/hg19/lo_pos_assem1_withSed.bed ../../data_files/hg19/lo_pos_assem1_withSed_result.bed --dir 1
#mapdb 17971
#example results: chr1/contig_1   382     383     chr1/contig_1   1       51184500        51184501

#/home/cmb-16/mjc/quentin/illuVeri/use_mcutils/mcutils/bin/samLiftover ../../data_files/hg19/HG002_to_hg37_asm20_woSed_assem2.sam ../../data_files/hg19/lo_pos_assem2_withSed.bed ../../data_files/hg19/lo_pos_assem2_withSed_result.bed --dir 1
#mapdb 18803




# %%
#get file to double check in IGV

#whole alignment bam file
#minimap2 -x asm20 -ac --secondary=no /panfs/qcb-panasas/jingwenr/reference/CLR/human_hs37d5.fasta /home/cmb-16/nobackups/mjc/data_download/HG002/CCS/assembly.1.consensus.fasta | samtools sort -o HG002_to_hg37_asm20_woSed_assem1.bam
#samtools index -@4 HG002_to_hg37_asm20_woSed_assem1.bam

#minimap2 -x asm20 -ac --secondary=no /panfs/qcb-panasas/jingwenr/reference/CLR/human_hs37d5.fasta /home/cmb-16/nobackups/mjc/data_download/HG002/CCS/assembly.2.consensus.fasta | samtools sort -o HG002_to_hg37_asm20_woSed_assem2.bam
#samtools index -@4 HG002_to_hg37_asm20_woSed_assem2.bam
