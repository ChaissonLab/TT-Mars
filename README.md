# TT-Mars

TT-Mars: S**t**ructural Varian**t**s Assess**m**ent B**a**sed on Haplotype-**r**esolved A**s**semblies.

## Usage

0. Clone TT-Mars from github and `cd TT-Mars`. Python >= 3.8 is preferred.
1. Create environment and activate: `conda create -n ttmars` and `conda activate ttmars`.
2. Run `dowaload_files.sh` to download required files to `./ttmars_files`.
3. Run `download_asm.sh` to download assembly files of 10 samples from HGSVC.
4. Install packages: `conda install -c bioconda pysam`, `conda install -c anaconda numpy`, `conda install -c bioconda mappy`, `conda install -c conda-forge biopython`, `conda install -c bioconda pybedtools`.
5. Run TT-Mars with following steps: `run_ttmars.sh` includes more instructions. Users can run it to run TT-Mars after setting up.

The main program: run `python ttmars.py -h` for help.

`python ttmars.py output_dir files_dir centro_file vcf_file reference asm_h1 asm_h2 tr_file num_X_chr`

## Positional arguments

1. `output_dir`: Output directory.  
2. `files_dir`: Input files directory. `./ttmars_files/sample_name`. The directory where you store required files after running `dowaload_files.sh`.  
3. `centro_file`: provided centromere file.  
4. `vcf_file`: callset file callset.vcf(.gz).  
5. `reference`: referemce file reference_genome.fasta.  
6. `asm_h1`: assembly files assembly1.fa, which were downloaded after running `download_asm.sh`.  
7. `asm_h2`: assembly files assembly2.fa, which were downloaded after running `download_asm.sh`.  
8. `tr_file`: provided tandem repeats file. 
9. `num_X_chr`: if male sample: 1; if female sample: 2.

## Optional arguments

`-n/--not_hg38`: if reference is NOT hg38/chm13 (hg19).  
`-p/--passonly`: if consider PASS calls only.  
`-s/--seq_resolved`: if consider sequence resolved calls.  
`-w/--wrong_len`: if count wrong length calls as True.  
`-g/--gt_vali`: conduct genotype validation.  
`-i/--gt_info`: index with GT info. (For phased callsets)  
`-d/--phased `: take phased information. (For phased callsets)  
`-v/--vcf_out`: output results as vcf files (tp (true positive), fp (false positive) and na).  
`-f/--false_neg`: output recall, must be used together with `-t/--truth_file`.  
`-t/--truth_file`: input truth vcf file, must be used together with `-f/--false_neg`.  

## Example Output

ttmars_combined_res.txt:  
|SV index| relative length| relative score| validation result| chr| start| end| Type| Genotype Match|
| :----: | :----: |  :----: | :----: | :----: | :----: | :----: |:----: | :----: |
|0|	1.0|	3.48|	True|	chr1|	249912|	249912| INS| True|

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
| HG002 (hg19 ref) | https://figshare.com/ndownloader/files/31455682 | https://figshare.com/ndownloader/files/31455676 |  https://figshare.com/ndownloader/files/31455685 |   https://figshare.com/ndownloader/files/31455679   |


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
| HG002 (hg19 ref) | https://figshare.com/ndownloader/files/31455670 | https://figshare.com/ndownloader/files/31455673 | 
