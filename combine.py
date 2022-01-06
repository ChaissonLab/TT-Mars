import csv
import sys

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("output_dir",
                    help="output directory")
parser.add_argument("no_X_chr",
                    choices=[1, 2],
                    help="male sample 1, female sample 2",
                    type=int)
parser.add_argument("-v",
                    "--vcf_out",
                    help="output results as vcf files, must be used together with -f/--vcf_file",
                    action="store_true")
parser.add_argument("-f",
                    "--vcf_file",
                    help="input vcf file using as template, must be used together with -v/--vcf_out")
parser.add_argument("-g",
                    "--gt_vali",
                    help="conduct genotype validation",
                    action="store_true")
args = parser.parse_args()

if bool(args.vcf_out) ^ bool(args.vcf_file):
    parser.error('-v/--vcf_out and -f/--vcf_file must be given together')

output_dir = args.output_dir + "/"

if int(args.no_X_chr) == 1:
    if_male = True
else:
    if_male = False
    
if args.vcf_out:
    if_vcf = True
    in_vcf_file = args.vcf_file
else:
    if_vcf = False
    
other_sv_res_file = output_dir+"ttmars_res.txt"
regdup_res_file = output_dir+"ttmars_regdup_res.txt"

with open(other_sv_res_file) as f:
    reader = csv.reader(f, delimiter="\t")
    other_sv_res = list(reader)
f.close()

with open(regdup_res_file) as f:
    reader = csv.reader(f, delimiter="\t")
    regdup_res = list(reader)
f.close()

sv_dict = {}
#results of SVs other than interspersed DUP
for rec in other_sv_res:
    ref_name = rec[0]
    sv_pos = int(rec[1])
    sv_end = int(rec[2])
    sv_type = rec[3]
    
    if not args.gt_vali:
        sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6]]
    else:
        sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6], rec[7]]
#interspersed DUP
for rec in regdup_res:
    ref_name = rec[0]
    sv_pos = int(rec[1])
    sv_end = int(rec[2])
    sv_type = rec[3]
    if (ref_name, int(sv_pos), int(sv_end), sv_type) in sv_dict:
        if not args.gt_vali:
            if rec[6] == 'True' and sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)][2] == 'False':
                sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6]]
        else:
            if rec[6] == 'True' and sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)][2] == 'False':
                sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6], rec[7]]
            elif rec[6] == 'True' and sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)][2] == 'True':
                if rec[7] == 'True' and sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)][3] == 'False':
                    sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6], rec[7]]
    else:
        if not args.gt_vali:
            sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6]]
        else:
            sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6], rec[7]]
        
if if_male:
    chrx_res_file = output_dir+"ttmars_chrx_res.txt"
    with open(chrx_res_file) as f:
        reader = csv.reader(f, delimiter="\t")
        chrx_res = list(reader)
    f.close()
    for rec in chrx_res:
        if len(rec) == 0:
            continue
        ref_name = rec[0]
        sv_pos = int(rec[1])
        sv_end = int(rec[2])
        sv_type = rec[3]
        
        if (ref_name, int(sv_pos), int(sv_end), sv_type) in sv_dict:
            if not args.gt_vali:
                if rec[6] == 'True' and sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)][2] == 'False':
                    sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6]]
            else:
                if rec[6] == 'True' and sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)][2] == 'False':
                    sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6], rec[7]]
                elif rec[6] == 'True' and sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)][2] == 'True':
                    if rec[7] == 'True' and sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)][3] == 'False':
                        sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6], rec[7]]
        else:
            if not args.gt_vali:
                sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6]]
            else:
                sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6], rec[7]]
        
#         if (ref_name, int(sv_pos), int(sv_end), sv_type) in sv_dict:
#             if rec[6] == 'True' and sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)][2] == 'False':
#                 sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6]]
#         else:
#             sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6]]    

g = open(output_dir + "/ttmars_combined_res.txt", "w")
for key in sv_dict:
    res = sv_dict[key]

    g.write(str(key[0]) + "\t")
    g.write(str(key[1]) + "\t")
    g.write(str(key[2]) + "\t")
    g.write(str(key[3]) + "\t")
    g.write(str(res[0]) + "\t")
    g.write(str(res[1]) + "\t")
    g.write(str(res[2]))
    
    if args.gt_vali:
        g.write("\t" + str(res[3]))
    
    g.write("\n")
g.close()

#if output vcf
if if_vcf:
    from pysam import VariantFile
    vcf_in = VariantFile(in_vcf_file)
    vcfh = vcf_in.header
    #vcfh.add_meta('INFO', items=[('ID',"TTMars"), ('Number',1), ('Type','String'),('Description','TT-Mars NA12878 results: TP, FA, NA or .')])
    vcfh.add_meta('INFO', items=[('ID',"GT_vali"), ('Number',1), ('Type','String'),('Description','TT-Mars GT validation (require flag -g): True, False or NA')])
    vcf_out_tp = VariantFile(output_dir+"ttmars_tp.vcf", 'w', header=vcfh)
    vcf_out_fp = VariantFile(output_dir+"ttmars_fp.vcf", 'w', header=vcfh)
    vcf_out_na = VariantFile(output_dir+"ttmars_na.vcf", 'w', header=vcfh)
    
    for rec in vcf_in.fetch():
        ref_name = rec.chrom
        sv_type = rec.info['SVTYPE']
        sv_pos = rec.pos
        sv_end = rec.stop
        
        try:
            validation_res = sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)][2]
            if validation_res == 'True':
                if args.gt_vali:
                    rec.info['GT_vali'] = sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)][3]
                vcf_out_tp.write(rec)
            elif validation_res == 'False':
                vcf_out_fp.write(rec)
        except:
            vcf_out_na.write(rec)

    vcf_out_tp.close()
    vcf_out_fp.close()
    vcf_out_na.close()
    vcf_in.close()    

#remove files
import os
for name in ['assem1_non_cov_regions.bed', 'assem2_non_cov_regions.bed',
             'exclude_assem1_non_cover.bed', 'exclude_assem2_non_cover.bed',
             'SV_positions.bed', 'ttmars_chrx_res.txt', 'ttmars_regdup_res.txt',
             'ttmars_res.txt', 'align_info_assem1_chrall.txt', 'align_info_assem2_chrall.txt',
             'all_reg_dup.fasta', 'all_reg_dup.fasta.fai']:
    if os.path.exists(output_dir + name):
        os.remove(output_dir + name)








