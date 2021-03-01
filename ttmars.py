import sys

import get_conf_int.py as get_conf_int
import validate.py as validate
import get_align_info.py as get_align_info


#main function
def main():
    #input
    output_dir = sys.argv[1] + "/"
    if_hg38_input = sys.argv[2]
    centromere_file = sys.argv[3]
    #exclude_assem1_non_cover_file = sys.argv[4]
    #exclude_assem2_non_cover_file = sys.argv[5]
    #exclude_high_depth_file = sys.argv[6]
    #assembly bam files
    bam_file1 = sys.argv[4]
    bam_file2 = sys.argv[5]
    avg_read_depth = sys.argv[6]
    read_bam_file = sys.argv[7]
    vcf_file = sys.argv[8]
    
    #constants
    interval = 20
    if_hg38 = False
    if if_hg38_input == "True":
        if_hg38 = True
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
    
    #build centromere dictionary
    dict_centromere = validate.build_centro_dict(centromere_file)
    
    #build lists for excluded SV positions
    
    #Output regions on ref where its not covered by at least one of the assembly
    get_conf_int.get_non_cover_regions(output_dir, bam_file1, 1, chr_list)
    get_conf_int.get_non_cover_regions(output_dir, bam_file2, 2, chr_list)
    
    #Get regions where read depth > 2 * avg_read_depth
    get_conf_int.get_high_depth_calls_info(output_dir, read_bam_file, vcf_file)
    
    #get filtered sv info, using results from get_conf_int.py 
    exclude_assem1_non_cover, exclude_assem2_non_cover, exclude_high_depth = 
        validate.get_filtered_sv_pos(output_dir + "assem1_non_cov_regions.bed", 
                                     output_dir + "assem2_non_cov_regions.bed", 
                                     output_dir + "exclude_high_depth.bed")
    
    #get validation info files
    
    
    
    
    #validate by both haplotypes: return a dict containing validation info
    dict_comb = vali_info(output_dir, 
                          exclude_assem1_non_cover, 
                          exclude_assem2_non_cover, 
                          exclude_high_depth, 
                          "align_info_assem1_chrall.txt",
                          "align_info_assem2_chrall.txt")
    
    #write output
    validate.write_output(output_dir, dict_comb)

if __name__ == "__main__":
    main()