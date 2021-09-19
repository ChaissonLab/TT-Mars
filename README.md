# TT-Mars

TT-Mars: S**t**ructural Varian**t**s Assess**m**ent B**a**sed on Haplotype-**r**esolved A**s**semblies.  
A sequence content based structural variants (SVs) validation method which aims to efficiently and accurately assess SV calls through assembly contigs and established alignment information.  

## Requirement

Python 3 packages: pysam (https://github.com/pysam-developers/pysam).  
Softwares: samtools (https://github.com/samtools/samtools), samliftover (https://github.com/mchaisso/mcutils).

## Positional parameters

1. Output directory  
2. if hg38: if reference is hg38 (True/False). If False, TT-Mars will use hs37d5 as reference  
3. centromere file: centromere_positions.bed  
4. Assembly alignment bam file: assembly1.bam (haplotype-resolved)  
5. Assembly alignment bam file: assembly2.bam (haplotype-resolved)  
6. Average read depth: the average read depth of the read file  
7. Read file: read.bam (for high-depth-read calls only, need to modify the method to get high-depth regions)  
8. A callset file: callset.vcf.gz  
9. Referemce genome sequence file: reference_genome.fasta  
10. Assembly sequence files: assembly1.fasta/.fa (haplotype-resolved)  
11. Assembly sequence files: assembly2.fasta/.fa (haplotype-resolved)  
12. liftover file: assembly1_liftover.bed (provided)  
13. liftover file: assembly2_liftover.bed (provided)

## Usage

python ttmars.py output_dir if_hg38 centromere_positions.bed assembly1.bam assembly2.bam avg_read_depth read.bam callset.vcf.gz assembly1.fasta assembly2.fasta assembly1_liftover.bed assembly2_liftover.bed

## Example Output

ttmars_res.txt:  
(chr start end type relative_length relative_score validation_result)  
chr1	893792	893827	DEL	1.03	3.18	True

## Available Data

### Liftover files and assemblies  
| Samples      | Liftover Hap1 | Liftover Hap2     | Assembly Liftover Hap1 | Assembly Liftover Hap2     |
| :----:      |    :----:   |        :----: |    :----:   |        :----: |
| HG00096 | https://ndownloader.figshare.com/files/26663234 | https://ndownloader.figshare.com/files/26663231 |   |      |
| HG00171 | https://ndownloader.figshare.com/files/27072881 | https://ndownloader.figshare.com/files/27072878 |   |      |
| HG00513 | https://ndownloader.figshare.com/files/27073232 | https://ndownloader.figshare.com/files/27073241 |   |      |
| HG00731 | https://ndownloader.figshare.com/files/27073718 | https://ndownloader.figshare.com/files/27073721 |   |      |
| HG00732 | https://ndownloader.figshare.com/files/27074015 | https://ndownloader.figshare.com/files/27074018 |   |      |
| HG00864 | https://ndownloader.figshare.com/files/27076040 | https://ndownloader.figshare.com/files/27076085 |   |      |
| HG01596 | https://ndownloader.figshare.com/files/27076994 | https://ndownloader.figshare.com/files/27077000 |   |      |
| HG03009 |  |  |   |      |
