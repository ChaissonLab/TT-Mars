
from Bio import SeqIO
import csv
import math
#pysam: https://pysam.readthedocs.io/en/latest/usage.html
import pysam
#print(pysam.__version__)

from Bio.Seq import Seq
from Bio import pairwise2

f = pysam.VariantFile('/home/cmb-16/mjc/quentin/illuVeri/callset_files/giab/HG002_GRCh37_1_22_v4.1_draft_benchmark_snpOnly_sorted.vcf.gz')

count_confident = 0

for info in dots_truvari_tp_fp:
    whole_info = dict_comb[str(info[2])]
    ref_name = whole_info[8]
    ref_start = int(whole_info[9])
    ref_end = int(whole_info[10])
    
    SNV_dist_start = 10000000000
    SNV_dist_end = 10000000000
    for counter, rec in enumerate(f.fetch()):
        name = rec.chrom
        if rec.samples[0]["GT"] == (0, 1) and name == ref_name:
            sv_pos = rec.pos
            if SNV_dist_start > abs(sv_pos - ref_start):
                SNV_dist_start = abs(sv_pos - ref_start)
            if SNV_dist_end > abs(sv_pos - ref_end):
                SNV_dist_end = abs(sv_pos - ref_end)
    smalleset_dist = 0
    if SNV_dist_start < SNV_dist_end:
        smalleset_dist = smalleset_dist
    else:
        smalleset_dist = SNV_dist_end
    
    print(smalleset_dist)

    count_confident += 1
    #if count_confident >10:
    #    break
f.close()