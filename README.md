# TT-Mars

TT-Mars: S**t**ructural Varian**t**s Assess**m**ent B**a**sed on Haplotype-**r**esolved A**s**semblies.

## Requirement

Python 3 packages: pysam (https://github.com/pysam-developers/pysam). 

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

## Accompanying Resources

### Liftover files  
| Samples      | Reference Liftover Hap1 | Reference Liftover Hap2 | Assembly Liftover Hap1 | Assembly Liftover Hap2 |
| :----:      |    :----:   |        :----: |    :----:   |        :----: |
| HG00096 | https://figshare.com/ndownloader/files/30817390 | https://figshare.com/ndownloader/files/30817384 | https://figshare.com/ndownloader/files/30817387  |  https://figshare.com/ndownloader/files/30817381   |
| HG00171 | https://figshare.com/ndownloader/files/30817402  | https://figshare.com/ndownloader/files/30817396 |  https://figshare.com/ndownloader/files/30817399 |  https://figshare.com/ndownloader/files/30817393    |
| HG00513 | https://figshare.com/ndownloader/files/30817411  | https://figshare.com/ndownloader/files/30817405 | https://figshare.com/ndownloader/files/30817408  |   https://figshare.com/ndownloader/files/30817414   |
| HG00731 |  https://figshare.com/ndownloader/files/30817426 | https://figshare.com/ndownloader/files/30817420 | https://figshare.com/ndownloader/files/30817423  |  https://figshare.com/ndownloader/files/30817417    |
| HG00732 |  https://figshare.com/ndownloader/files/30817435 | https://figshare.com/ndownloader/files/30817429 | https://figshare.com/ndownloader/files/30817432  |   https://figshare.com/ndownloader/files/30817438   |
| HG00864 |  https://figshare.com/ndownloader/files/30817450 | https://figshare.com/ndownloader/files/30817444 | https://figshare.com/ndownloader/files/30817447  |   https://figshare.com/ndownloader/files/30817441   |
| HG01114 |  https://figshare.com/ndownloader/files/30817459 | https://figshare.com/ndownloader/files/30817453 |  https://figshare.com/ndownloader/files/30817456 |   https://figshare.com/ndownloader/files/30817462   |
| HG01505 | https://figshare.com/ndownloader/files/30817471  | https://figshare.com/ndownloader/files/30817465 | https://figshare.com/ndownloader/files/30817468  |   https://figshare.com/ndownloader/files/30817474   |
| HG01596 |  https://figshare.com/ndownloader/files/30817486 | https://figshare.com/ndownloader/files/30817480 | https://figshare.com/ndownloader/files/30817483  |   https://figshare.com/ndownloader/files/30817477   |
| HG03009 | https://figshare.com/ndownloader/files/30817498 | https://figshare.com/ndownloader/files/30817492 |  https://figshare.com/ndownloader/files/30817495 |   https://figshare.com/ndownloader/files/30817489   |


### Genome coverage files  
| Samples      | Reference Liftover Hap1 | Reference Liftover Hap2 |
| :----:      |    :----:   |        :----: |
| HG00096 | https://figshare.com/ndownloader/files/30850246 | https://figshare.com/ndownloader/files/30850249 |
| HG00171 | https://figshare.com/ndownloader/files/30850258 | https://figshare.com/ndownloader/files/30850261 |
| HG00513 | https://figshare.com/ndownloader/files/30850639 | https://figshare.com/ndownloader/files/30850642 | 
| HG00731 |  https://figshare.com/ndownloader/files/30850663 | https://figshare.com/ndownloader/files/30850660 | 
| HG00732 | https://figshare.com/ndownloader/files/30850687 | https://figshare.com/ndownloader/files/30850681 |
| HG00864 | https://figshare.com/ndownloader/files/30850708 | https://figshare.com/ndownloader/files/30850711 | 
| HG01114 | https://figshare.com/ndownloader/files/30850726 | https://figshare.com/ndownloader/files/30850729 | 
| HG01505 | https://figshare.com/ndownloader/files/30850747  | https://figshare.com/ndownloader/files/30850744 | 
| HG01596 | https://figshare.com/ndownloader/files/30850768 | https://figshare.com/ndownloader/files/30850762 |
| HG03009 | https://figshare.com/ndownloader/files/30850777 | https://figshare.com/ndownloader/files/30850780 | 
