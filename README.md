# TT-Mars

TT-Mars: S**t**ructural Varian**t**s Assess**m**ent B**a**sed on Haplotype-**r**esolved A**s**semblies.

## Usage

0. Clone TT-Mars from github and `cd TT-Mars`.
1. Run dowaload_files.sh to download required files to `./ttmars_files`.
2. Create environment and activate: `conda create -n ttmars` and `conda activate ttmars_test`.
3. Install packages: `conda install -c bioconda pysam`, `conda install -c anaconda numpy`, `conda install -c bioconda mappy`, `conda install -c conda-forge biopython`, `conda install -c bioconda pybedtools`.
4. Run TT-Mars with following steps: more instructions in `run_ttmars.sh`.

`python ttmars.py output_dir if_hg38 centro_file files_dir/assem1_non_cov_regions.bed files_dir/assem2_non_cov_regions.bed vcf_file reference asm_h1 asm_h2 files_dir/lo_pos_assem1_result_compressed.bed files_dir/lo_pos_assem2_result_compressed.bed tr_file pass_only seq_resolved`

`python reg_dup.py output_dir if_hg38 centro_file files_dir/assem1_non_cov_regions.bed files_dir/assem2_non_cov_regions.bed vcf_file reference asm_h1 asm_h2 files_dir/lo_pos_assem1_result_compressed.bed files_dir/lo_pos_assem2_result_compressed.bed tr_file files_dir/lo_pos_assem1_0_result_compressed.bed files_dir/lo_pos_assem2_0_result_compressed.bed pass_only`

`python chrx.py output_dir if_hg38 centro_file files_dir/assem1_non_cov_regions.bed files_dir/assem2_non_cov_regions.bed vcf_file reference asm_h1 asm_h2 files_dir/lo_pos_assem1_result_compressed.bed files_dir/lo_pos_assem2_result_compressed.bed tr_file pass_only seq_resolved`

`python combine.py output_dir num_X_chr`

## Positional parameters

1. output_dir: Output directory.
2. if hg38: if reference is hg38 (True/False). 
3. centro_file: provided centromere file. 
4. tr_file: provided tandem repeats file.
5. vcf_file: callset file callset.vcf(.gz)  
6. reference: referemce file reference_genome.fasta  
7. asm_h1/2: assembly files assembly1/2.fa 
8. assem1_non_cov_regions.bed, assem2_non_cov_regions.bed, lo_pos_assem1_result_compressed.bed, lo_pos_assem2_result_compressed.bed, lo_pos_assem1_0_result_compressed.bed, lo_pos_assem2_0_result_compressed.bed: required files, downloaded to `./ttmars_files`

## Example Output

ttmars_res.txt:  
|chr| start| end| type| relative length| relative score| validation result|
|chr1|	893792|	893827|	DEL|	1.03|	3.18|	True|

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
| Samples      | Hap1 | Hap2 |
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
