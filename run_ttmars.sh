#!/bin/bash

# conda create -n ttmars
# conda activate ttmars
# conda install -c bioconda pysam
# conda install -c anaconda numpy
# conda install -c bioconda mappy
# conda install -c conda-forge biopython
# conda install -c bioconda pybedtools

#sample: HG00096 HG00171 HG00513 HG00731 HG00732 HG00864 HG01114 HG01505 HG01596 HG03009
sample=HG00096
reference=path-to-reference_file/hg38.no_alts.fasta
vcf_file=path-to-target-vcf/callset.vcf
#assembly files, can be downloaded by download_asm.sh
asm_h1=path-to-assemblies/h1.fa
asm_h2=path-to-assemblies/h2.fa
output_dir=output-directory

#path to downloaded ttmars files: default ./ttmars_files/sample
files_dir=./ttmars_files/$sample

#provided centromere file
centro_file=centromere_hg38.txt
#provided tandem repeats file
tr_file=hg38_tandem_repeats.bed
#1: if male sample; 2: if female sample
num_X_chr=1

python ttmars.py "$output_dir" "$centro_file" "$files_dir"/assem1_non_cov_regions.bed "$files_dir"/assem2_non_cov_regions.bed "$vcf_file" "$reference" "$asm_h1" "$asm_h2" "$files_dir"/lo_pos_assem1_result_compressed.bed "$files_dir"/lo_pos_assem2_result_compressed.bed "$files_dir"/lo_pos_assem1_0_result_compressed.bed "$files_dir"/lo_pos_assem2_0_result_compressed.bed "$tr_file" 1000 "$num_X_chr" -s -g -w -d -i -v

# positional arguments:
#   output_dir            output directory
#   centromere_file       centromere file
#   assem1_non_cov_regions_file
#                         Regions that are not covered on hap1
#   assem2_non_cov_regions_file
#                         Regions that are not covered on hap2
#   vcf_file              input vcf file
#   ref_file              reference file
#   query_file1           assembly fasta file hap1
#   query_file2           assembly fasta file hap2
#   liftover_file1        liftover file hap1
#   liftover_file2        liftover file hap2
#   liftover_file1_0      liftover file hap1 asm to ref
#   liftover_file2_0      liftover file hap2 asm to ref
#   tandem_file           tandem repeats regions
#   region_len_m          region_len_m
#   {1,2}                 male sample 1, female sample 2

# optional arguments:
#   -h, --help            show this help message and exit
#   -n, --not_hg38        if reference is NOT hg38 (hg19)
#   -p, --passonly        if consider PASS calls only
#   -s, --seq_resolved    if consider sequence resolved calls (INS)
#   -w, --wrong_len       if count wrong length calls as True
#   -g, --gt_vali         conduct genotype validation
#   -i, --gt_info         index with GT info
#   -d, --phased          take phased information
#   -v, --vcf_out         output results as vcf files, must be used together with -f/--vcf_file
#   -f, --false_neg       output false negative, must be used together with -t/--truth_file and -f/--vcf_file
#   -t TRUTH_FILE, --truth_file TRUTH_FILE
#                         input truth vcf file, must be used together with -n/--false_neg