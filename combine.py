import csv
import sys

output_dir = sys.argv[1] + "/"
male = sys.argv[2]

if male == "True":
    if_male = True
else:
    if_male = False
    
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
for rec in other_sv_res:
    ref_name = rec[0]
    sv_pos = int(rec[1])
    sv_end = int(rec[2])
    sv_type = rec[3]
    sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6]]
for rec in regdup_res:
    ref_name = rec[0]
    sv_pos = int(rec[1])
    sv_end = int(rec[2])
    sv_type = rec[3]
    if (ref_name, int(sv_pos), int(sv_end), sv_type) in sv_dict:
        if rec[6] == 'True' and sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)][2] == 'False':
            sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6]]
    else:
        sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6]]
        
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
            if rec[6] == 'True' and sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)][2] == 'False':
                sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6]]
        else:
            sv_dict[(ref_name, int(sv_pos), int(sv_end), sv_type)] = [rec[4], rec[5], rec[6]]    

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
    g.write("\n")
g.close()

#remove files
# import os
# for name in ['assem1_non_cov_regions.bed', 'assem2_non_cov_regions.bed',
#              'exclude_assem1_non_cover.bed', 'exclude_assem2_non_cover.bed',
#              'SV_positions.bed', 'ttmars_chrx_res.txt', 'ttmars_regdup_res.txt',
#              'ttmars_res.txt', 'align_info_assem1_chrall.txt', 'align_info_assem2_chrall.txt']:
#     if os.path.exists(output_dir + name):
#         os.remove(output_dir + name)
    
    







