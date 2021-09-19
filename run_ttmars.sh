#!/bin/bash
#SBATCH --account=mchaisso_100
#SBATCH --partition=qcb
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --mem=40gb
#SBATCH --time=2:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jianzhiy@usc.edu
source /project/mchaisso_100/cmb-16/quentin/software/miniconda3/bin/activate base

sample=HG00096
caller=wham
aligner=lra
reference=/panfs/qcb-panasas/jianzhiy/data/reference/hg38.no_alts.fasta
cd /project/mchaisso_100/cmb-16/quentin/ttmars/
output_dir=../test_output_0917/
asm_dir=/scratch2/jianzhiy/data/assemblies/"$sample"/
vcf_file=/scratch2/jianzhiy/ttmars/callset_files/$caller/$sample/"$sample"_filtered_sorted.vcf.gz
centro_file=centromere_hg38.txt
tr_file=hg38_tandem_repeats.bed
if_hg38=True
pass_only=True
seq_resolved=False  

python ttmars.py "$output_dir" "$if_hg38" "$centro_file" "$asm_dir"/"$aligner"/assem1_nool_sort.bam "$asm_dir"/"$aligner"/assem2_nool_sort.bam "$vcf_file" "$reference" "$asm_dir"/h1.fa "$asm_dir"/h2.fa "$asm_dir"/"$aligner"/lo_pos_assem1_result_compressed.bed "$asm_dir"/"$aligner"/lo_pos_assem2_result_compressed.bed "$tr_file" "$pass_only" "$seq_resolved"

python reg_dup.py "$output_dir" "$if_hg38" "$centro_file" "$asm_dir"/"$aligner"/assem1_nool_sort.bam "$asm_dir"/"$aligner"/assem2_nool_sort.bam "$vcf_file" "$reference" "$asm_dir"/h1.fa "$asm_dir"/h2.fa "$asm_dir"/"$aligner"/lo_pos_assem1_result_compressed.bed "$asm_dir"/"$aligner"/lo_pos_assem2_result_compressed.bed "$tr_file" "$asm_dir"/"$aligner"/lo_pos_assem1_0_result_compressed.bed "$asm_dir"/"$aligner"/lo_pos_assem2_0_result_compressed.bed "$pass_only"

python combine_dup.py "$output_dir"

# sample=HG00096
# caller=delly
# aligner=lra
# reference=/path-to-reference/hg38.no_alts.fasta
# output_dir=/path-to-output-dir/
# asm_dir=/assemblies/"$sample"/
# vcf_file=/path-to-vcf-file/"$sample".vcf
# centro_file=centromere_hg38.txt
# tr_file=hg38_tandem_repeats.bed
# if_hg38=True
# pass_only=True
# seq_resolved=False

# python ttmars.py "$output_dir" "$if_hg38" "$centro_file" "$asm_dir"/"$aligner"/assem1_nool_sort.bam "$asm_dir"/"$aligner"/assem2_nool_sort.bam "$vcf_file" "$reference" "$asm_dir"/h1.fa "$asm_dir"/h2.fa "$asm_dir"/"$aligner"/lo_pos_assem1_result_compressed.bed "$asm_dir"/"$aligner"/lo_pos_assem2_result_compressed.bed "$tr_file" "$pass_only" "$seq_resolved"

# python reg_dup.py "$output_dir" "$if_hg38" "$centro_file" "$asm_dir"/"$aligner"/assem1_nool_sort.bam "$asm_dir"/"$aligner"/assem2_nool_sort.bam "$vcf_file" "$reference" "$asm_dir"/h1.fa "$asm_dir"/h2.fa "$asm_dir"/"$aligner"/lo_pos_assem1_result_compressed.bed "$asm_dir"/"$aligner"/lo_pos_assem2_result_compressed.bed "$tr_file" "$asm_dir"/"$aligner"/lo_pos_assem1_0_result_compressed.bed "$asm_dir"/"$aligner"/lo_pos_assem2_0_result_compressed.bed "$pass_only"

# python combine_dup.py "$output_dir"