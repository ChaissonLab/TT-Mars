from Bio import SeqIO
import csv
import math

from pysam import VariantFile

from Bio.Seq import Seq
from Bio import pairwise2
import sys    
    
def check_ol_coor(sv_coor_list, trimmed_coor):
    name = sv_coor_list[0]
    cur_start = sv_coor_list[1]
    cur_end = sv_coor_list[2]
    for rec in trimmed_coor:
        if str(name) == str(rec[0]):
            trim_start = int(rec[1])
            trim_end = int(rec[2])
            cur_start = int(cur_start)
            cur_end = int(cur_end)
            
            if cur_start <= trim_end and trim_start <= cur_start:
                return True
    return False

def filter_dipcall_for_fn_trimcoor(bcf_in_file, asm1_trimmed_file, asm2_trimmed_file, bcf_out_file):
#     print(sample_name)
    
    del_ctr = 0
    ins_ctr = 0
        
    bcf_in = VariantFile(bcf_in_file)

    ##############################################
    ##############################################
    #trim coordinate files
    with open(asm1_trimmed_file) as f:
        reader = csv.reader(f, delimiter="\t")
        asm1_trimmed_coor = list(reader)
    f.close()
    with open(asm2_trimmed_file) as f:
        reader = csv.reader(f, delimiter="\t")
        asm2_trimmed_coor = list(reader)
    f.close()
    ##############################################
    ##############################################
        
    vcfh = bcf_in.header
    bcf_out = VariantFile(bcf_out_file, 'w', header=vcfh)

    #g = open("./HG00514/HG00514_filtered.vcf", "a")
    for counter, rec in enumerate(bcf_in.fetch()):
        name = rec.chrom  
                
        cur_start = rec.pos
        cur_end = rec.stop
        cur_length = 1

        ##############################################
        ##############################################
        #filter dipcall sv if its overlapping with trimmed coordinates
        if check_ol_coor([name, cur_start, cur_end], asm1_trimmed_coor):
            continue
        if check_ol_coor([name, cur_start, cur_end], asm2_trimmed_coor):
            continue
        ##############################################
        ##############################################

        bcf_out.write(rec)
        
    #     if counter > 200:
    #         break    
        #bcf_out.write(new_rec)

    bcf_out.close()
    bcf_in.close()
    
#     print(del_ctr, ins_ctr)
    
def main():
    #hg38 samples
    sample_names = ["HG00096", "HG01505", "HG01596", "HG03009", "HG00731", "HG00171", "HG00864", "HG01114", "HG00513", "HG00732"]
    chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7',
                'chr8','chr9','chr10','chr11','chr12','chr13','chr14',
                'chr15','chr16','chr17','chr18','chr19','chr20','chr21',
                'chr22','chrX']

    truth_sv_size_lb = 50

    for sample_name in sample_names:
        bcf_in_file = sample_name + "_dip_sv_lenfilter_chrfilter_pass.vcf"
    
        asm1_trimmed_file = "assem1_sort_trimmed_coor.bed"
        asm2_trimmed_file = "assem2_sort_trimmed_coor.bed"
        
    #     filter_dipcall_for_fn(sample_name)
        filter_dipcall_for_fn_trimcoor(bcf_in_file, asm1_trimmed_file, asm2_trimmed_file, bcf_out_file)

if __name__ == "__main__":
    main()