#######################################################
#######################################################
#1. Align assembly to reference

#use lra (https://github.com/ChaissonLab/LRA) to align asm to ref
lra index -CONTIG $reference
lra align -CONTIG $reference h1.fa -t 16 -p s | samtools sort -o assem1_sort.bam
lra align -CONTIG $reference h2.fa -t 16 -p s | samtools sort -o assem2_sort.bam

#index bam file
samtools index assem1_sort.bam
samtools index assem2_sort.bam

#######################################################
#######################################################
#2. Trim overlapping contigs

#trim overlapping contigs
python trim_overlapping_contigs.py assem1_sort.bam $output_dir $if_hg38
python trim_overlapping_contigs.py assem2_sort.bam $output_dir $if_hg38

#sort trimmed bam file
samtools sort $output_dir/assem1_sort_nool.bam -o $output_dir/assem1_nool_sort.bam
samtools sort $output_dir/assem2_sort_nool.bam -o $output_dir/assem2_nool_sort.bam

#index sorted trimmed file
samtools index $output_dir/assem1_nool_sort.bam
samtools index $output_dir/assem2_nool_sort.bam 

#######################################################
#######################################################
#3. Liftover

#convert to sam file
samtools view -h $output_dir/assem1_nool_sort.bam | samtools sort -O sam -o $output_dir/assem1_nool_sort.sam
samtools view -h $output_dir/assem2_nool_sort.bam | samtools sort -O sam -o $output_dir/assem2_nool_sort.sam

#liftover using samLiftover (https://github.com/mchaisso/mcutils): ref to asm lo
python lo_assem_to_ref.py $output_dir $output_dir/assem1_nool_sort.bam $output_dir/assem2_nool_sort.bam

samLiftover $output_dir/assem1_nool_sort.sam $output_dir/lo_pos_assem1.bed $output_dir/lo_pos_assem1_result.bed --dir 1
samLiftover $output_dir/assem2_nool_sort.sam $output_dir/lo_pos_assem2.bed $output_dir/lo_pos_assem2_result.bed --dir 1

#liftover using samLiftover (https://github.com/mchaisso/mcutils): asm to ref lo
python lo_assem_to_ref_0.py $output_dir $output_dir/assem1_nool_sort.bam $output_dir/assem2_nool_sort.bam

samLiftover $output_dir/assem1_nool_sort.sam $output_dir/lo_pos_assem1_0.bed $output_dir/lo_pos_assem1_0_result.bed --dir 0
samLiftover $output_dir/assem2_nool_sort.sam $output_dir/lo_pos_assem2_0.bed $output_dir/lo_pos_assem2_0_result.bed --dir 0    

#######################################################
#######################################################
#4. Compress liftover files

python compress_liftover.py $output_dir lo_pos_assem1_result.bed lo_pos_assem1_result_compressed.bed
python compress_liftover.py $output_dir lo_pos_assem2_result.bed lo_pos_assem2_result_compressed.bed
python compress_liftover.py $output_dir lo_pos_assem1_0_result.bed lo_pos_assem1_0_result_compressed.bed
python compress_liftover.py $output_dir lo_pos_assem2_0_result.bed lo_pos_assem2_0_result_compressed.bed

#######################################################
#######################################################
#5. Get non-covered regions

python get_conf_int.py $output_dir $output_dir/assem1_nool_sort.bam $output_dir/assem2_nool_sort.bam $if_hg38
