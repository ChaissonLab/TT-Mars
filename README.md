# TT-Mars

TT-Mars: S**t**ructural Varian**t**s Assess**m**ent B**a**sed on Haplotype-**r**esolved A**s**semblies.  
A sequence content based structural variants (SVs) validation method which aims to efficiently and accurately assess SV calls through assembly contigs and established alignment information.  

## Requirement

Python 3 packages: pysam (https://github.com/pysam-developers/pysam).  
Softwares: samtools (https://github.com/samtools/samtools), samliftover (https://github.com/mchaisso/mcutils).

## Input

1. A callset file: callset.vcf.gz  
2. Read file: read.bam (for high-depth-read calls only, need to modify the method to get high-depth regions)  
2.1. The average read depth: the average read depth of the read file  
3. assembly1/2.bam (haplotype-resolved): [will go through the trimming process]
4. reference_genome.fasta
5. sample_assembly1/2.fasta/.fa (haplotype-resolved)
6. output_dir
7. if_hg38: if ref is hg38 (True/False)
8. centromere_hg37/38_.bed

## Usage

## Output