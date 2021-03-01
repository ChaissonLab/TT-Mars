# TT-Mars

TT-Mars: S**t**ructural Varian**t**s Assess**m**ent B**a**sed on Haplotype-**r**esolved A**s**semblies.  
A sequence content based structural variants (SVs) validation method which aims to efficiently and accurately assess SV calls through assembly contigs and established alignment information.  

## Requirement

Python 3 packages: pysam (https://github.com/pysam-developers/pysam).  
Softwares: samtools (https://github.com/samtools/samtools), samliftover (https://github.com/mchaisso/mcutils).

## Input

1. Output directory  
2. if_hg38: if reference is hg38 (True/False). If False, TT-Mars will use hs37d5 as reference  
3. centromere file: centromere_positions.bed  
4-5. Assembly alignment bam files: assembly1.bam and assembly2.bam (haplotype-resolved)  
6. Average read depth: the average read depth of the read file  
7. Read file: read.bam (for high-depth-read calls only, need to modify the method to get high-depth regions)  
8. A callset file: callset.vcf.gz  
9. Referemce genome sequence file: reference_genome.fasta  
10-11. Assembly sequence files: assembly1.fasta/.fa and assembly2.fasta/.fa (haplotype-resolved)
12-13. liftover file: assembly1_liftover.bed and assembly2_liftover.bed

## Usage

python 

## Output

ttmars_res.txt: chr start end type relative_length relative_score validation_result  
Example: chr1	893792	893827	DEL	1.03	3.18	True