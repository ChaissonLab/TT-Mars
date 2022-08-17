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

python ttmars.py "$output_dir" "$centro_file" "$files_dir"/assem1_non_cov_regions.bed "$files_dir"/assem2_non_cov_regions.bed "$vcf_file" "$reference" "$asm_h1" "$asm_h2" "$files_dir"/lo_pos_assem1_result_compressed.bed "$files_dir"/lo_pos_assem2_result_compressed.bed "$files_dir"/lo_pos_assem1_0_result_compressed.bed "$files_dir"/lo_pos_assem2_0_result_compressed.bed "$tr_file"
# optional arguments: 
# -n/--not_hg38: if reference is NOT hg38 (hg19).
# -p/--passonly: if consider PASS calls only. 
# -s/--seq_resolved: if consider sequence resolved calls (INS). 
# -w/--wrong_len: if count wrong length calls as True.
# -g/--gt_vali: conduct genotype validation.

python combine.py "$output_dir" "$num_X_chr"
#optional arguments:
# -v/--vcf_out: output results as vcf files (tp (true positive), fp (false positive) and na), must be used together with -f/--vcf_file.  
# -f VCF_FILE/--vcf_file VCF_FILE: input vcf file, use as template.  
# -g/--gt_vali: conduct genotype validation.
# -n/--false_neg: output recall, must be used together with -t/--truth_file and -f/--vcf_file.
# -t/--truth_file: input truth vcf file, must be used together with -n/--false_neg.