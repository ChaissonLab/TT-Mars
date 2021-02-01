# %%

#change truvari dir and output folder (wham/lumpy/manta)

#Generate plots
import csv
import math
#pysam: https://pysam.readthedocs.io/en/latest/usage.html
import pysam
#print(pysam.__version__)

import numpy as np 
import matplotlib.pyplot as plt  


# %%
#constants
#test stash

interval = 500
output_dir = "../output/naive_caller_0810_v1/"
bam_file = "/panfs/qcb-panasas/jingwenr/result/HG002_CCS/HG002.lra.bam"
avg_depth = 40.04
#bam_file = "/home/cmb-16/mjc/sample-datasets/giab/HG002/Illumina/HG002.hs37d5.ILL.bam"
#avg_depth = 33.3876


# %%
#build centromere position dictionary

#get centromere file
with open("../data_files/centromere_hg37.txt") as f:
    reader = csv.reader(f, delimiter="\t")
    centromere_raw = list(reader)
f.close()

#build dictionay
dict_centromere = dict()
pre_centro_chr = ""
start_centro_pos = ""
for record in centromere_raw:
    if record[0] == pre_centro_chr:
        end_centro_pos = str(record[2])
        dict_centromere[record[0]] = (start_centro_pos, end_centro_pos)
    else:
        start_centro_pos = str(record[1])
        pre_centro_chr = record[0]


# %%
#build lists for excluded SV positions

with open(output_dir + "exclude__assem1_result.bed") as f:
    reader = csv.reader(f, delimiter="\t")
    exclude__assem1_result = list(reader)
f.close()

with open(output_dir + "exclude__assem2_result.bed") as f:
    reader = csv.reader(f, delimiter="\t")
    exclude__assem2_result = list(reader)
f.close()


# %%
#update dictionary given a list of scores
#in dict str is stored
#TODO: modify here

#check if current cases is in the excluded list
def check_exclude(list_to_check):
	if (list_to_check in exclude__assem1_result) or (list_to_check in exclude__assem2_result):
		return True
	else:
		return False

#get depth of coverage given position of a chromosome
def get_depth(ref_name, ref_pos, bam_file):
    pos_str = ref_name + ':' + str(int(ref_pos) - 1) + '-' + str(ref_pos)
    res = pysam.depth("-r", pos_str, bam_file)
    if res=='':
        return 0
    start = 0
    end = len(res) - 1
    for i in range(len(res) - 1, -1, -1):
        if res[i] == '\t':
            start = i + 1
            break
    return int(res[start:end])

#test
#bam = pysam.AlignmentFile('/home/cmb-16/mjc/hgsvg/datasets/HG00514/Illumina/HG00514.bam', 'rb')
#bam.count_coverage('chr6', 1, 5, quality_threshold = 0)

#in align_info:
# SV_index score_before score_after no_found SV_type hap:1 or 2
# (counter alignment_beforeSV[0][2] alignment_afterSV[0][2] len(query_start_dic) sv_type hap)
# contig_seg_len ref_seg_len SV_len
# (len(query_frag) len(ref_frag) sv_end-sv_pos+1)
# chr_name SV_start SV_end interval_start interval_end
# (ref_name sv_pos sv_end ref_start ref_end)
# error message
# (err_mes)

def updateDict(dict_score, align_info):
	for record in align_info:
		#some cases are not aligned successfully due to function's (memory) issues: give -1 in that column
		if int(record[3]) > 0:
			#filter out centromere cases
			index = str(record[0])
			ref_name = str(record[9])
			ref_start = int(record[12])
			ref_end = int(record[13])
			centro_start = int(dict_centromere['chr'+ref_name][0])
			centro_end = int(dict_centromere['chr'+ref_name][1])

			#if SV in the exclude_list: start or end of a contig
			sv_pos = int(record[10])
			sv_end = int(record[11])
			list_to_check = [str(ref_name), str(sv_pos), str(sv_end)]
			if check_exclude(list_to_check):
				continue

			#if ref start or ref end in centromere, skip this case
			if (ref_start > centro_start and ref_start < centro_end) or (ref_end > centro_start and ref_end < centro_end):
				continue
            
			#filter out call at high coverage location
			ref_start_depth = get_depth(ref_name, ref_start, bam_file)
			ref_end_depth = get_depth(ref_name, ref_end, bam_file)
			if ((ref_start_depth > 2 * avg_depth) or (ref_end_depth > 2 * avg_depth)):
				continue

			if record[0] not in dict_score:
				dict_score[record[0]] = record[1:len(record)]
			else:
				#TODO: may change
				#choose the better relative score one:
				#if one of the before score is 0, choose the better relative length one:
				'''
				if float(dict_score[record[0]][0]) == 0 or float(record[1]) == 0:
					old_rela_len = (float(dict_score[record[0]][5]) - float(dict_score[record[0]][6]))/abs(float(dict_score[record[0]][7]))
					new_rela_len = (float(record[6]) - float(record[7]))/abs(float(record[8]))
					if new_rela_len < old_rela_len:
						dict_score.update({record[0]: record[1:len(record)]})
					continue
				old_rela_change = (float(dict_score[record[0]][1]) - float(dict_score[record[0]][0]))/abs(float(dict_score[record[0]][0]))
				new_rela_change = (float(record[2]) - float(record[1]))/abs(float(record[1]))
				if new_rela_change > old_rela_change:
					dict_score.update({record[0]: record[1:len(record)]})
				'''
				#choose the better relative length one:
				old_rela_len = (float(dict_score[record[0]][5]) - float(dict_score[record[0]][6]))/abs(float(dict_score[record[0]][7]))
				new_rela_len = (float(record[6]) - float(record[7]))/abs(float(record[8]))
				if abs(new_rela_len + 1) < abs(old_rela_len + 1):
					dict_score.update({record[0]: record[1:len(record)]})
	return dict_score


# %%

#scatter plot of relative length VS relative score

# %%
#plot relative length VS relative score

#plot as needed: when ass!_align_score.txt[i][3] != -1
#some in align_len_info_1, not in ass1_align_score.txt
#store plot info in align_len_plot_info_1 and align_len_plot_info_2

with open(output_dir + "align_info_assem1_chrall.txt") as f:
	reader = csv.reader(f, delimiter="\t")
	align_info_assem1 = list(reader)
f.close()

with open(output_dir + "align_info_assem2_chrall.txt") as f:
	reader = csv.reader(f, delimiter="\t")
	align_info_assem2 = list(reader)
f.close()

#dicts for single haplotypes and combined
dict_1 = dict()
dict_1 = updateDict(dict_1, align_info_assem1)
dict_2 = dict()
dict_2 = updateDict(dict_2, align_info_assem2)
dict_comb = dict()
dict_comb = updateDict(dict_comb, align_info_assem1)
dict_comb = updateDict(dict_comb, align_info_assem2)

align_plot_info_1 = []
for record in align_info_assem1:
	if record[0] in dict_1:
		#TODO: handle this case
		#if zero score_before
		if float(dict_1[record[0]][0]) == 0:
			continue
		#round to 2 decimal points
		rela_len = round((float(record[6]) - float(record[7]))/float(record[8]), 2)
		rela_score = round((float(dict_1[record[0]][1]) - float(dict_1[record[0]][0]))/abs(float(dict_1[record[0]][0])), 2)
		align_plot_info_1.append([rela_len, rela_score, int(record[0])])

align_plot_info_2 = []
for record in align_info_assem2:
	if record[0] in dict_2:
		#TODO: handle this case
		#if zero score_before
		if float(dict_2[record[0]][0]) == 0:
			continue
		#round to 2 decimal points
		rela_len = round((float(record[6]) - float(record[7]))/float(record[8]), 2)
		rela_score = round((float(dict_2[record[0]][1]) - float(dict_2[record[0]][0]))/abs(float(dict_2[record[0]][0])), 2)
		align_plot_info_2.append([rela_len, rela_score, int(record[0])])

#combine two haplotypes
#TODO: debug this part: check the results!!!
align_plot_info_comb = []
for record in dict_comb:
		#TODO: solve 0 score problem: should be included!!!
		#if zero score
		if float(dict_comb[record][0]) == 0:
			#rela_score = 0.01
			continue
		else:
			#round to 2 decimal points
			rela_score = round((float(dict_comb[record][1]) - float(dict_comb[record][0]))/abs(float(dict_comb[record][0])), 2)
		rela_len = round((float(dict_comb[record][5]) - float(dict_comb[record][6]))/float(dict_comb[record][7]), 2)

		align_plot_info_comb.append([rela_len, rela_score, int(record)])


# %%
#plot relative length VS relative score

#haplotype 1
fig = plt.figure(figsize = (10, 5)) 
ax=fig.add_axes([0,0,1,1])
ax.scatter([row[1] for row in align_plot_info_1], [row[0] for row in align_plot_info_1], color='b', s=5)
#ax.plot([-10000,2000], [-10000,2000])
#wrong length cases
#ax.scatter([row[1] for row in align_len_plot_wrong_len_hap1], [row[0] for row in align_len_plot_wrong_len_hap1], color='r', s=5)
ax.set_xlabel('relative score')
ax.set_ylabel('relative length')
ax.set_title('haplotype 1')
ax.set_xlim((-5,5))
ax.set_ylim((-2,2))
plt.show()

# %%
#plot relative length VS relative score

#haplotype 2
fig = plt.figure(figsize = (10, 5)) 
ax=fig.add_axes([0,0,1,1])
ax.scatter([row[1] for row in align_plot_info_2], [row[0] for row in align_plot_info_2], color='b', s=5)
#ax.plot([-10000,2000], [-10000,2000])
#wrong length cases
#ax.scatter([row[1] for row in align_len_plot_wrong_len_hap2], [row[0] for row in align_len_plot_wrong_len_hap2], color='r', s=5)
ax.set_xlabel('relative score')
ax.set_ylabel('relative length')
ax.set_title('haplotype 2')
ax.set_xlim((-5,5))
ax.set_ylim((-2,2))
plt.show()


# %%
#plot relative length VS relative score

#haplotypes combined
fig = plt.figure(figsize = (10, 5)) 
ax=fig.add_axes([0,0,1,1])
ax.scatter([row[1] for row in align_plot_info_comb], [row[0] for row in align_plot_info_comb], color='b', s=5)
#ax.plot([-10000,2000], [-10000,2000])
#wrong length cases
#ax.scatter([row[1] for row in align_len_plot_wrong_len_both], [row[0] for row in align_len_plot_wrong_len_both], color='r', s=5)
#plot reference lines
ax.plot([-0.0,-0.0], [-2,1])
ax.plot([2,-5], [-0.0,-0.0])
ax.set_xlabel('relative score')
ax.set_ylabel('relative length')
ax.set_title('haplotypes combined')
ax.set_xlim((-5,5))
ax.set_ylim((-3,3))
plt.show()

# %%
# add truvari's results to the plot

truvari_fp_index = []

f = pysam.VariantFile('/home/cmb-16/mjc/quentin/illuVeri/use_mcutils/run_truvari/truvari_12/fp.vcf','r')

for counter, rec in enumerate(f.fetch()):
	#test
	#if counter < 91:
	#    continue

	ref_name = rec.chrom
	sv_type = rec.info['SVTYPE']
	sv_len = rec.rlen
	#TODOL double check the start for different types
	sv_pos = rec.pos
	sv_end = rec.stop

	for key, record in dict_comb.items(): 
		if int(record[9]) == sv_pos and int(record[10]) == sv_end:
			truvari_fp_index.append(int(key))

align_plot_truvari_fp_both = []
for ind in truvari_fp_index:
    for record in align_plot_info_comb:
        if record[2] == ind:
            align_plot_truvari_fp_both.append(record)
            break


# %%
# add truvari's results to the plot
#from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('../../results/plot.pdf')

#haplotypes combined
fig = plt.figure(figsize = (8, 8)) 
ax=fig.add_axes([0,0,1,1])
ax.scatter([row[1] for row in align_plot_info_comb], [row[0] for row in align_plot_info_comb], color='b', s=5)
#ax.plot([-10000,2000], [-10000,2000])
#truvari fp cases
ax.scatter([row[1] for row in align_plot_truvari_fp_both], [row[0] for row in align_plot_truvari_fp_both], color='r', s=5)
#plot reference lines
ax.plot([-0.0,-0.0], [-5,5], linewidth=0.5)
ax.plot([5,-5], [-0.0,-0.0], linewidth=0.5)
#ax.set_xlabel('relative score')
#ax.set_ylabel('relative length')
#ax.set_title('Relative length VS Relative score for GIAB callset')
ax.set_xlim((-5,5))
ax.set_ylim((-5,5))
plt.xlabel('relative score')
plt.ylabel('relative length')
plt.title('Relative length VS Relative score for GIAB callset')
#plt.savefig('../../results/plot.pdf')
plt.show()




# %%
#get file to double check in IGV

#whole alignment bam file
#minimap2 -x asm20 -ac --secondary=no /panfs/qcb-panasas/jingwenr/reference/CLR/human_hs37d5.fasta /home/cmb-16/nobackups/mjc/data_download/HG002/CCS/assembly.1.consensus.fasta | samtools sort -o HG002_assem1_hs37d5.bam
#samtools index -@4 HG002_assem1_hs37d5.bam

#minimap2 -x asm20 -ac --secondary=no /panfs/qcb-panasas/jingwenr/reference/CLR/human_hs37d5.fasta /home/cmb-16/nobackups/mjc/data_download/HG002/CCS/assembly.2.consensus.fasta | samtools sort -o HG002_assem2_hs37d5.bam
#samtools index -@4 HG002_assem2_hs37d5.bam

# %%
#Analysis

# add truvari's results to the plot

#truvari's false positive
truvari_fp_index = []

no_tp_truvari = 0
no_fp_truvari = 0

f = pysam.VariantFile('/home/cmb-16/mjc/quentin/illuVeri/use_mcutils/run_truvari/truvari_12/fp.vcf','r')

for counter, rec in enumerate(f.fetch()):
    no_fp_truvari = no_fp_truvari + 1
    #test
    #if counter < 91:
    #    continue

    ref_name = rec.chrom
    sv_type = rec.info['SVTYPE']
    sv_len = rec.rlen
    #TODOL double check the start for different types
    sv_pos = rec.pos
    sv_end = rec.stop

    for key, record in dict_comb.items(): 
        if int(record[9]) == sv_pos and int(record[10]) == sv_end:
            if str(record[8]) == str(ref_name):
                truvari_fp_index.append(int(key))

align_plot_truvari_fp_both = []
for ind in truvari_fp_index:
    for record in align_plot_info_comb:
        if record[2] == ind:
            align_plot_truvari_fp_both.append(record)
            break

#truvari's true positive
truvari_tp_index = []

f = pysam.VariantFile('/home/cmb-16/mjc/quentin/illuVeri/use_mcutils/run_truvari/truvari_12/tp-call.vcf','r')

for counter, rec in enumerate(f.fetch()):
    no_tp_truvari = no_tp_truvari + 1
    #test
    #if counter < 91:
    #    continue

    ref_name = rec.chrom
    sv_type = rec.info['SVTYPE']
    sv_len = rec.rlen
    #TODOL double check the start for different types
    sv_pos = rec.pos
    sv_end = rec.stop

    for key, record in dict_comb.items(): 
        if int(record[9]) == sv_pos and int(record[10]) == sv_end:
            if str(record[8]) == str(ref_name):
                truvari_tp_index.append(int(key))

align_plot_truvari_tp_both = []
for ind in truvari_tp_index:
    for record in align_plot_info_comb:
        if record[2] == ind:
            align_plot_truvari_tp_both.append(record)
            break


# %%

# add truvari's results to the plot

#haplotypes combined
fig = plt.figure(figsize = (10, 5)) 
ax=fig.add_axes([0,0,1,1])
ax.scatter([row[1] for row in align_plot_info_comb], [row[0] for row in align_plot_info_comb], color='b', s=5)
#ax.plot([-10000,2000], [-10000,2000])
#truvari fp cases
ax.scatter([row[1] for row in align_plot_truvari_fp_both], [row[0] for row in align_plot_truvari_fp_both], color='r', s=5)
#truvari tp cases
ax.scatter([row[1] for row in align_plot_truvari_tp_both], [row[0] for row in align_plot_truvari_tp_both], color='g', s=5)
#plot reference lines
ax.plot([-0.0,-0.0], [-2,1])
ax.plot([2,-5], [-0.0,-0.0])
ax.set_xlabel('relative score')
ax.set_ylabel('relative length')
ax.set_title('haplotypes combined')
ax.set_xlim((-5,5))
ax.set_ylim((-5,5))
plt.show()


# %%

#get information to anlayze truvari's and my results

#total number of dots
no_dots = len(align_plot_info_comb)
#total number of truvari TP
#no_tp_truvari
#total number of truvari TP on the plot
no_dots_tp = len(align_plot_truvari_tp_both)
#total number of truvari FP
#no_fp_truvari
#total number of truvari FP on the plot
no_dots_fp = len(align_plot_truvari_fp_both)

#list of index of cases in the plot that have |rela length + 1| > 0.8, but TP by truvari
dots_truvari_tp_fp = []
dots_truvari_tp_wlen = []
dots_truvari_tp_tp = []
for record in align_plot_truvari_tp_both:
	if abs(float(record[0]) + 1) > 0.8:
		dots_truvari_tp_fp.append(record)
	elif abs(float(record[0]) + 1) < 0.1:
		dots_truvari_tp_tp.append(record)
	else:
		dots_truvari_tp_wlen.append(record)

#list of index of cases in the plot that have |rela length + 1| < 0.1, but FP by truvari
dots_truvari_fp_tp = []
dots_truvari_fp_wlen = []
dots_truvari_fp_fp = []
for record in align_plot_truvari_fp_both:
	if abs(float(record[0]) + 1) > 0.8:
		dots_truvari_fp_fp.append(record)
	elif abs(float(record[0]) + 1) < 0.1:
		dots_truvari_fp_tp.append(record)
	else:
		dots_truvari_fp_wlen.append(record)


# %%
print(no_dots, no_dots_tp, no_dots_fp)
#2853 2853 0
#%%
print(len(dots_truvari_tp_fp), len(dots_truvari_tp_wlen), len(dots_truvari_tp_tp))
#79 366 2408
#%%
print(len(dots_truvari_fp_tp), len(dots_truvari_fp_wlen), len(dots_truvari_fp_fp))
#0 0 0


# %%
dots_truvari_tp_fp

# %%
dots_truvari_tp_wlen

# %%
dots_truvari_fp_tp
# %%
dots_truvari_fp_wlen


# %%
dict_comb['2453']

# %%
dict_comb['18']

#%%
counter_my_fp = 0
counter_my_tp = 0
counter_my_wlen = 0
for record in align_plot_info_comb:
	if abs(float(record[0]) + 1) > 0.8:
		counter_my_fp = counter_my_fp + 1
	elif abs(float(record[0]) + 1) < 0.1:
		counter_my_tp = counter_my_tp + 1
	else:
		counter_my_wlen = counter_my_wlen + 1


#%%
print(counter_my_tp, counter_my_wlen, counter_my_fp)
		

# %%
dict_2['1116']

# %%
f = open("../../results/align_info_combined_chrall.txt", "a")
#counter = 1
counter_my_fp = 0
counter_my_tp = 0
counter_my_wlen = 0
for key, value in dict_comb.items():
	f.write(str(key) + "\t")
	for i in range(0, len(value) -1):
		f.write(str(value[i]) + "\t")
	rela_len = round((float(value[5]) - float(value[6]))/float(value[7]), 2)
	f.write(str(rela_len) + "\t")

	my_class = ""
	#distinguish FP, TP and wrong length by relative length
	if abs(rela_len + 1) > 0.8:
		counter_my_fp = counter_my_fp + 1
		my_class = "FP"
	elif abs(rela_len + 1) < 0.1:
		counter_my_tp = counter_my_tp + 1
		my_class = "TP"
	else:
		counter_my_wlen = counter_my_wlen + 1
		my_class = "WrongLen"

	f.write(str(my_class) + "\t")
	f.write("\n")
	#print(value)
	#counter = counter + 1
	#if counter > 2:
	#	break
f.close()


#%%
print(counter_my_tp, counter_my_wlen, counter_my_fp)

# %%
len(dict_comb)

# %%
