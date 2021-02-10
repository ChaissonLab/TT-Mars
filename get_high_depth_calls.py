#Get the bed file of 

import pysam
import sys 
import os
import numpy as np

#get calls which have too high avg depth, output an ecclude_bed file
import sys 
  
#get command line input
#n = len(sys.argv)
output_dir = sys.argv[1]
vcf_file = sys.argv[2]
bam_file = sys.argv[3]

#output_dir = "../../output/v4.0/naive_caller_0810_v1/"
#bam_file = "/panfs/qcb-panasas/jingwenr/result/HG002_CCS/HG002.lra.bam"
#avg_depth = 40.04
#bam_file = "/panfs/qcb-panasas/jianzhiy/data/illumina/HG00096/HG00096.bam"
#avg_depth = 31.9256
#vcf_file = "../../../callset_files/naive_caller/0810_v1_result/lra_ccs_0810_v1_sorted.vcf.gz"

avg_depth = 30
samfile = pysam.AlignmentFile(bam_file, "rb")

f = pysam.VariantFile(vcf_file,'r')
for counter, rec in enumerate(f.fetch()):
    g = open(output_dir + "exclude_high_depth.bed", "a")
    #get ref start and ref end
    name = rec.chrom
    sv_pos = rec.pos
    sv_end = rec.stop
    res = samfile.count_coverage(name, sv_pos, sv_end+1, quality_threshold = 0)
    #print(res)
    if round(np.sum(res)/(sv_end+1-sv_pos), 2) > 2*avg_depth:
        #test
        #print(counter, round(np.sum(res)/(sv_end+1-sv_pos), 2))
        g.write(str(name) + "\t")
        g.write(str(sv_pos) + "\t")
        g.write(str(sv_end))
        g.write("\n")
            
    g.close()
