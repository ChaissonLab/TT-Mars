#Get regions where there are reads < k bp
import pysam
import sys 
import os
import numpy as np

#Get regions on ref where its not covered by at least one of the assembly
def get_non_cover_regions(output_dir, assem_bam_file, hap, chr_list):
    #hap is an int = 1/2
    samfile = pysam.AlignmentFile(assem_bam_file, "rb")
    g = open(output_dir + "assem" + str(hap) + "_non_cov_regions.bed", "w")
    for chr_name in chr_list:
        #test
        #print(chr_name)
        #stop when parsed no_of_reads
        #no_of_reads = 400000
        cur_end = 0
        #loop through iter, index will not be reset
        #if if_hg38:
        #    ref_name = "chr" + chr_name
        #else:
        #    ref_name = chr_name
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

#Get regions where read depth > 2 * avg_read_depth
#For now, we filter calls by read depth
#In other words, here output calls having high read depth
#TODO: Too slow here
def get_high_depth_calls_info(output_dir, read_bam_file, vcf_file):
    samfile = pysam.AlignmentFile(read_bam_file, "rb")
    f = pysam.VariantFile(vcf_file,'r')
    #TODO: change open condition
    g = open(output_dir + "exclude_high_depth.bed", "a")
    for counter, rec in enumerate(f.fetch()):
        #get ref start and ref end
        name = rec.chrom
        sv_pos = rec.pos
        sv_end = rec.stop
        res = samfile.count_coverage(name, sv_pos, sv_end+1, quality_threshold = 0)
        #print(res)
        if round(np.sum(res)/(sv_end+1-sv_pos), 2) > 2*avg_read_depth:
            #test
            #print(counter, round(np.sum(res)/(sv_end+1-sv_pos), 2))
            g.write(str(name) + "\t")
            g.write(str(sv_pos) + "\t")
            g.write(str(sv_end))
            g.write("\n")
    g.close()

def main():
    #get command line input
    #n = len(sys.argv)
    output_dir = sys.argv[1] + "/"
    #assembly bam file
    bam_file1 = sys.argv[2]
    bam_file2 = sys.argv[3]
    if_hg38_str = sys.argv[4]
    #if_hg38 = True
    avg_read_depth = sys.argv[5]
    #reads bam file
    read_bam_file = sys.argv[6]
    #callset file
    vcf_file = sys.argv[7]
    
    #constants
    if if_hg38_str == "True":
        if_hg38 = True
    else:
        if_hg38 = False 
    
    chr_list = []
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
    
    #Output regions on ref where its not covered by at least one of the assembly
    get_non_cover_regions(output_dir, bam_file1, 1, chr_list)
    get_non_cover_regions(output_dir, bam_file2, 2, chr_list)
    
    #Get regions where read depth > 2 * avg_read_depth
    get_high_depth_calls_info(output_dir, read_bam_file, vcf_file)
    
if __name__ == "__main__":
    main()