#Get regions where there are reads < k bp
import pysam
import sys 
import os
import sys 
  
#get command line input
#n = len(sys.argv)
output_dir = sys.argv[1]
bam_file1 = sys.argv[2]
bam_file2 = sys.argv[3]

k = 250000

#bam_file = "/scratch2/jingwenr/assembly/HG00514/HG00514_CCS_SS_PG_PRR/h1.bam"
#bam_file1 = "../../data_files/hg19//HG002_to_hg37_asm20_woSed_assem1.bam"
#bam_file2 = "../../data_files/hg19//HG002_to_hg37_asm20_woSed_assem2.bam"

if_hg38 = True
#output_dir = "../../output/v4.0/HG00514_CCS_SS_PG_PRR/"
#output_dir = "../../output/v4.0/naive_caller_0810_v1/"

chr_list = ["1", "2", "3", "4", "5",
            "6", "7", "8", "9", "10",
            "11", "12", "13", "14", "15",
            "16", "17", "18", "19", "20",
            "21", "22", "X", "Y"]

#chr_list = ["1"]

samfile = pysam.AlignmentFile(bam_file1, "rb")
g = open(output_dir + "assem1_short_reads_regions_" + str(k) + ".bed", "w")
for chr_name in chr_list:
    #test
    #print(chr_name)
    #stop when parsed no_of_reads
    #no_of_reads = 400000
    counter_read = 0
    
    #loop through iter, index will not be reset
    if if_hg38:
        ref_name = "chr" + chr_name
    else:
        ref_name = chr_name
        
    iter = samfile.fetch(ref_name)
    
    for rec in iter:
        #print(rec.infer_query_length())
        if (rec.infer_query_length() < k):
            g.write(str(ref_name) + "\t")
            g.write(str(rec.reference_start) + "\t")
            g.write(str(rec.reference_end))
            g.write("\n")
g.close()

samfile = pysam.AlignmentFile(bam_file2, "rb")
g = open(output_dir + "assem2_short_reads_regions_" + str(k) + ".bed", "w")
for chr_name in chr_list:
    #test
    #print(chr_name)
    #stop when parsed no_of_reads
    #no_of_reads = 400000
    counter_read = 0
    
    #loop through iter, index will not be reset
    if if_hg38:
        ref_name = "chr" + chr_name
    else:
        ref_name = chr_name
        
    iter = samfile.fetch(ref_name)
    
    for rec in iter:
        #print(rec.infer_query_length())
        if (rec.infer_query_length() < k):
            g.write(str(ref_name) + "\t")
            g.write(str(rec.reference_start) + "\t")
            g.write(str(rec.reference_end))
            g.write("\n")
g.close()


######################################################################################
######################################################################################

#Get regions on ref where its not covered by at least one of the assembly

#output_dir = "../../output/v4.0/HG00514_CCS_SS_PG_PRR/"
#output_dir = "../../output/v4.0/naive_caller_0810_v1/"

chr_list = ["1", "2", "3", "4", "5",
            "6", "7", "8", "9", "10",
            "11", "12", "13", "14", "15",
            "16", "17", "18", "19", "20",
            "21", "22", "X", "Y"]

#chr_list = ["1"]

samfile = pysam.AlignmentFile(bam_file1, "rb")
g = open(output_dir + "assem1_non_cov_regions.bed", "w")
for chr_name in chr_list:
    #test
    #print(chr_name)
    #stop when parsed no_of_reads
    #no_of_reads = 400000
    cur_end = 0
    
    #loop through iter, index will not be reset
    if if_hg38:
        ref_name = "chr" + chr_name
    else:
        ref_name = chr_name
    iter = samfile.fetch(ref_name)
    
    for rec in iter:
        if rec.reference_start > cur_end:
            g.write(str(ref_name) + "\t")
            g.write(str(cur_end) + "\t")
            g.write(str(rec.reference_start))
            g.write("\n")
            cur_end = rec.reference_end
        else:
            if rec.reference_end > cur_end:
                cur_end = rec.reference_end

g.close()

samfile = pysam.AlignmentFile(bam_file2, "rb")
g = open(output_dir + "assem2_non_cov_regions.bed", "w")
for chr_name in chr_list:
    #test
    #print(chr_name)
    #stop when parsed no_of_reads
    #no_of_reads = 400000
    cur_end = 0
    
    #loop through iter, index will not be reset
    if if_hg38:
        ref_name = "chr" + chr_name
    else:
        ref_name = chr_name
    iter = samfile.fetch(ref_name)
    
    for rec in iter:
        if rec.reference_start > cur_end:
            g.write(str(ref_name) + "\t")
            g.write(str(cur_end) + "\t")
            g.write(str(rec.reference_start))
            g.write("\n")
            cur_end = rec.reference_end
        else:
            if rec.reference_end > cur_end:
                cur_end = rec.reference_end

g.close()