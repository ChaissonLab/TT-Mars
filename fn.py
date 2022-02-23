from pysam import VariantFile

#CONST
###########################################################
#fn const
max_dist_to_merge = 1500
max_dist_search = 1000
ratio_size_lb = 0.7

###########################################################
#ttmars const
chr_list = ["chr1", "chr2", "chr3", "chr4", "chr5",
            "chr6", "chr7", "chr8", "chr9", "chr10",
            "chr11", "chr12", "chr13", "chr14", "chr15",
            "chr16", "chr17", "chr18", "chr19", "chr20",
            "chr21", "chr22", "chrX"]

memory_limit = 100000
memory_min = 10

#valid types
valid_types = ['DEL', 'INS', 'INV', 'DUP:TANDEM', 'DUP']
#pass_only
if_pass_only = True
wrong_len = False
gt_vali = False
if_hg38 = True

###########################################################

#first_filter: type, PASS, chr_name
def first_filter(sv, sv_type):
    #type filter
    if sv_type not in valid_types:
        return True
    #PASS filter
    if if_pass_only:
        if 'PASS' not in sv.filter.keys():
            return True
    chr_name = sv.chrom
    #chr filter
    if chr_name not in chr_list:
        return True
    return False


#define class
class struc_var:
    def __init__(self, idx, ref_name, sv_type, sv_pos, sv_stop, length, gt):
        self.idx = idx
        self.ref_name = ref_name
        self.sv_pos = sv_pos
        self.sv_stop = sv_stop
        self.sv_type = sv_type
        self.length = length
        self.gt = gt
        #if the call is part of an aggregate SV
        self.is_agg = False
        #if second filtered out
        self.is_sec_fil = False
        self.is_third_fil = False
        
        self.query_name_hap1 = "NA"
        self.query_name_hap2 = "NA"
        
        self.ref_start_best_hap1 = -1
        self.ref_end_best_hap1 = -1
        self.query_start_best_hap1 = -1
        self.query_end_best_hap1 = -1
        
        self.ref_start_best_hap2 = -1
        self.ref_end_best_hap2 = -1
        self.query_start_best_hap2 = -1
        self.query_end_best_hap2 = -1
        
        self.analyzed_hap1 = False
        self.analyzed_hap2 = False
        
        self.len_query_hap1 = -1
        self.len_ref_hap1 = -1
        self.len_query_hap2 = -1
        self.len_ref_hap2 = -1
        
        self.score_before_hap1 = -1
        self.score_after_hap1 = -1
        self.score_before_hap2 = -1
        self.score_after_hap2 = -1
        
        self.neg_strand_hap1 = False
        self.neg_strand_hap2 = False
        
        self.ins_seq = ""
        self.if_seq_resolved = False
        
    def check_tp(self, rela_len, rela_score):
        result = True
        if self.sv_type in ['DEL', 'DUP', 'DUP:TANDEM']:
            if rela_score >= 0 and rela_score <= 2.5:
                if rela_len >= -0.05*rela_score + 0.8 and rela_len <= 0.05*rela_score + 1.2:
                    result = True
                else:
                    result = False
            elif rela_score > 2.5:
                if rela_len >= 0.675 and rela_len <= 1.325:
                    result = True
                else:
                    result = False
            else:
                result = False
        elif self.sv_type == 'INS':
            #not seq-resolved
            #if len(self.ins_seq) == 0:
            if not self.if_seq_resolved:
                if rela_len < 0.675 or rela_len > 1.325:
                    result = False
            #seq-resolved
            else:
                if rela_score >= 0 and rela_score <= 2.5:
                    if rela_len >= -0.05*rela_score + 0.8 and rela_len <= 0.05*rela_score + 1.2:
                        result = True
                    else:
                        result = False
                elif rela_score > 2.5:
                    if rela_len >= 0.675 and rela_len <= 1.325:
                        result = True
                    else:
                        result = False
                else:
                    result = False                
                
        elif self.sv_type == 'INV':
            if rela_score <= 0:
                result = False
        return result
    
    #TP when wrong length flag presents -- looser rules for TP
    def check_tp_wlen(self, rela_len, rela_score):
        result = True
        if self.sv_type in ['DEL', 'DUP', 'DUP:TANDEM']:
            if rela_score >= 0 and rela_score <= 2.5:
                if rela_len >= -0.05*rela_score + 0.6 and rela_len <= 0.05*rela_score + 1.4:
                    result = True
                else:
                    result = False
            elif rela_score > 2.5:
                if rela_len >= 0.475 and rela_len <= 1.525:
                    result = True
                else:
                    result = False
            else:
                result = False
        elif self.sv_type == 'INS':
            #not seq-resolved
            #if len(self.ins_seq) == 0:
            if not self.if_seq_resolved:
                if rela_len < 0.475 or rela_len > 1.525:
                    result = False
            #seq-resolved
            else:
                if rela_score >= 0 and rela_score <= 2.5:
                    if rela_len >= -0.05*rela_score + 0.6 and rela_len <= 0.05*rela_score + 1.4:
                        result = True
                    else:
                        result = False
                elif rela_score > 2.5:
                    if rela_len >= 0.475 and rela_len <= 1.525:
                        result = True
                    else:
                        result = False
                else:
                    result = False                
                
        elif self.sv_type == 'INV':
            if rela_score <= 0:
                result = False
        return result
        
    def print_info(self):
        print(self.idx, self.ref_name, self.sv_pos, self.sv_stop, self.sv_type, self.length, self.gt, self.is_agg, self.is_sec_fil, self.is_third_fil)
        
    def cal_rela_score(self, score_before, score_after):
        if score_before > -1 and score_before < 0:
            tmp_score_before = -1
            tmp_score_after = score_after + (tmp_score_before - score_before)
            return round((tmp_score_after - tmp_score_before) / abs(tmp_score_before), 2)
        
        elif score_before >= 0 and score_before < 1:
            tmp_score_before = 1
            tmp_score_after = score_after + (tmp_score_before - score_before)
            return round((tmp_score_after - tmp_score_before) / abs(tmp_score_before), 2)
        
        else:
            return round((score_after - score_before) / abs(score_before), 2)
        
    def cal_rela_len(self, query_len, ref_len):
        return round((query_len - ref_len) / self.length, 2)
        
    def get_vali_res(self):
        if (not self.analyzed_hap1) or (not self.analyzed_hap2):
            return -1
        
        if self.analyzed_hap1 and self.analyzed_hap2:
            rela_len_1 = self.cal_rela_len(self.len_query_hap1, self.len_ref_hap1)
            rela_len_2 = self.cal_rela_len(self.len_query_hap2, self.len_ref_hap2)
            
            rela_score_1 = self.cal_rela_score(self.score_before_hap1, self.score_after_hap1)
            rela_score_2 = self.cal_rela_score(self.score_before_hap2, self.score_after_hap2)
            
            if not wrong_len:
                res_hap1 = self.check_tp(rela_len_1, rela_score_1)
                res_hap2 = self.check_tp(rela_len_2, rela_score_2)
            else:
                res_hap1 = self.check_tp_wlen(rela_len_1, rela_score_1)
                res_hap2 = self.check_tp_wlen(rela_len_2, rela_score_2)
            
            gt_validate = False
#             if args.gt_vali:
            if gt_vali:
                if res_hap1 and res_hap2:
                    if self.gt == (1,1):
                        gt_validate = True
                elif res_hap1 or res_hap2:
                    if self.gt == (1,0) or self.gt == (0,1):
                        gt_validate = True
                
            if res_hap1 and res_hap2:
                if abs(rela_len_1 - 1) <= abs(rela_len_2 - 1):
                    return (res_hap1, rela_len_1, rela_score_1, gt_validate)
                else:
                    return (res_hap2, rela_len_2, rela_score_2, gt_validate)
            elif res_hap1:
                return (res_hap1, rela_len_1, rela_score_1, gt_validate)
            elif res_hap2:
                return (res_hap2, rela_len_2, rela_score_2, gt_validate)
            else:
                if abs(rela_len_1 - 1) <= abs(rela_len_2 - 1):
                    return (res_hap1, rela_len_1, rela_score_1, gt_validate)
                else:
                    return (res_hap2, rela_len_2, rela_score_2, gt_validate)
                
                
def idx_sv(input_vcf):
    f = VariantFile(input_vcf,'r')
    sv_list = []
    for count, rec in enumerate(f.fetch()):
        #get sv_type
        try:
            sv_type = rec.info['SVTYPE']
        except:
            print("invalid sv type info")
            continue

        if first_filter(rec, sv_type):
            continue

        #get sv length
        if sv_type == 'INV':
            sv_len = abs(rec.stop - rec.pos + 1)
        else:
            try:
                sv_len = rec.info['SVLEN'][0]
            except:
                try:
                    sv_len = rec.info['SVLEN']
                except:
                    sv_len = abs(rec.stop - rec.pos + 1)
                    #print("invalid sv length info")
    #         try:
    #             sv_len = rec.info['SVLEN'][0]
    #         except:
    #             sv_len = rec.info['SVLEN']
        #handle del length > 0:
        if sv_type == 'DEL':
            sv_len = -abs(sv_len)

        if abs(sv_len) < memory_min:
            continue

        #get gt
        #only taking the first sample genotype 
        if gt_vali:
            sv_gt = rec.samples[0]["GT"]
            #bad genotype
            if sv_gt not in [(1, 1), (1, 0), (0, 1)]:
                sv_gt = None
        else:
            sv_gt = None

    #     if len(rec.samples.values()) != 1:
    #         raise Exception("Wrong number of sample genotype(s)")
    #     gts = [s['GT'] for s in rec.samples.values()] 

        sv_list.append(struc_var(count, rec.chrom, sv_type, rec.pos, rec.stop, sv_len, sv_gt))   

    f.close()
    
    return sv_list


# merge close SVs
def merge_sv(sv_list):
    merged_sv_list = []
    cur_idx = sv_list[0].idx
    cur_chrom = sv_list[0].ref_name
    cur_start = sv_list[0].sv_pos
    cur_end = sv_list[0].sv_stop
    cur_length = sv_list[0].length
    cur_gt = None
    
    for sv in sv_list[1:]:
        if sv.ref_name != cur_chrom or sv.sv_stop - cur_end > max_dist_to_merge:
            if cur_length > 0:
                cur_type = 'INS'
            else:
                cur_type = 'DEL'
            merged_sv_list.append(struc_var(cur_idx, cur_chrom, cur_type, cur_start, cur_end, cur_length, cur_gt))

            cur_idx = sv.idx
            cur_chrom = sv.ref_name
            cur_start = sv.sv_pos
            cur_end = sv.sv_stop
            cur_length = sv.length
            cur_gt = None

            continue

        cur_end = sv.sv_stop
        cur_length += sv.length

    merged_sv_list.append(struc_var(cur_idx, cur_chrom, cur_type, cur_start, cur_end, cur_length, cur_gt))
    
    return merged_sv_list



#return T if two sv match
def match_sv(base_sv, cand_sv):
    if base_sv.ref_name != cand_sv.ref_name:
        #test
#         print("name")
        return False
    if base_sv.sv_type != cand_sv.sv_type:
        #test
#         print("type")
        return False
    if min(abs(base_sv.length), abs(cand_sv.length)) / max(abs(base_sv.length), abs(cand_sv.length)) < ratio_size_lb:
        #test
#         print("size")
        return False

    if base_sv.sv_pos - cand_sv.sv_stop > max_dist_search:
        return False
    
    if cand_sv.sv_pos - base_sv.sv_stop > max_dist_search:
        return False
    
    return True

#index chr name: chr1-chrX: 1-23
def idx_chr_name(chr_name):
    if not if_hg38:
        chr_name = 'chr' + chr_name
    return chr_list.index(chr_name) + 1

#compare two sv position
#return T if sv1 is before sv2
def compare_sv_positions(sv1, sv2):
    sv1_chr_idx = idx_chr_name(sv1.ref_name)
    sv2_chr_idx = idx_chr_name(sv2.ref_name)
    
    if sv1_chr_idx < sv2_chr_idx:
        return True
    elif sv1_chr_idx > sv2_chr_idx:
        return False
    
    if sv1.sv_pos < sv2.sv_pos:
        return True
    else:
        return False
    
#merge 2 sorted sv_list
def merge_sv_list(sv_list1, sv_list2):
    merged_list = []
    ptr1 = 0
    ptr2 = 0
    
    while(True):
        if compare_sv_positions(sv_list1[ptr1], sv_list2[ptr2]):
            merged_list.append(sv_list1[ptr1])
            ptr1 += 1
        else:
            merged_list.append(sv_list2[ptr2])
            ptr2 += 1
        
        if ptr1 == len(sv_list1) or ptr2 == len(sv_list2):
            break
    
    if ptr1 == len(sv_list1):
        for ptr in range(ptr2, len(sv_list2)):
            merged_list.append(sv_list2[ptr])
    elif ptr2 == len(sv_list2):
        for ptr in range(ptr1, len(sv_list1)):
            merged_list.append(sv_list1[ptr])
    
    return merged_list

#sort SV list by chr and start pos
def sort_sv_list(sv_list):
    if len(sv_list) == 0:
        return sv_list
    
    if len(sv_list) == 1:
        return sv_list
    
    if len(sv_list) == 2:
        if compare_sv_positions(sv_list[0], sv_list[1]):
            return sv_list
        else:
            return [sv_list[1], sv_list[0]]
    
    mid_pt = len(sv_list)//2
    sorted_sv_list = merge_sv_list(sort_sv_list(sv_list[:mid_pt]), sort_sv_list(sv_list[mid_pt:]))
    return sorted_sv_list
    
    
def count_tp_base(merged_sv_list_sorted, cand_sv_list_sorted):
    tp_base_ctr = 0
    cand_start_idx = 0
    for count, base_sv in enumerate(merged_sv_list_sorted):

        cur_chrom = base_sv.ref_name
        cur_start = base_sv.sv_pos
        cur_end = base_sv.sv_stop
        cur_length = base_sv.length

        cur_cand_start_idx = cand_start_idx
        new_cand_start_idx = cand_start_idx

        for count1, cand_sv in enumerate(cand_sv_list_sorted[cur_cand_start_idx:]):
            if (cur_start - cand_sv.sv_stop > max_dist_search and base_sv.ref_name == cand_sv.ref_name) or (idx_chr_name(base_sv.ref_name) > idx_chr_name(cand_sv.ref_name)):
                new_cand_start_idx += 1
                continue

            if match_sv(base_sv, cand_sv):
                tp_base_ctr += 1
                break

            if (cand_sv.sv_pos - cur_end > max_dist_search and base_sv.ref_name == cand_sv.ref_name) or (idx_chr_name(base_sv.ref_name) < idx_chr_name(cand_sv.ref_name)):
                break

        cand_start_idx = new_cand_start_idx
    return tp_base_ctr


def match_sv_dist_only(base_sv, cand_sv):
    if base_sv.ref_name != cand_sv.ref_name:
        #test
#         print("name")
        return False
#     if base_sv.sv_type != cand_sv.sv_type:
#         #test
# #         print("type")
#         return False
#     if min(abs(base_sv.length), abs(cand_sv.length)) / max(abs(base_sv.length), abs(cand_sv.length)) < ratio_size_lb:
#         #test
# #         print("size")
#         return False

    if base_sv.sv_pos - cand_sv.sv_stop > max_dist_search:
        return False
    
    if cand_sv.sv_pos - base_sv.sv_stop > max_dist_search:
        return False
    
    return True

def count_tp_base_dist_only(sv_list_sorted, cand_sv_list_sorted):
    tp_base_ctr = 0
    cand_start_idx = 0
    for count, base_sv in enumerate(sv_list_sorted):
        #test
    #     if count % 3000 == 0:
    #         print(count, tp_base_ctr)
    #         base_sv.print_info()

        cur_chrom = base_sv.ref_name
        cur_start = base_sv.sv_pos
        cur_end = base_sv.sv_stop
        cur_length = base_sv.length

        cur_cand_start_idx = cand_start_idx
        new_cand_start_idx = cand_start_idx

        for count1, cand_sv in enumerate(cand_sv_list_sorted[cur_cand_start_idx:]):
            if (cur_start - cand_sv.sv_stop > max_dist_search and base_sv.ref_name == cand_sv.ref_name) or (idx_chr_name(base_sv.ref_name) > idx_chr_name(cand_sv.ref_name)):
                new_cand_start_idx += 1
                continue

            if match_sv_dist_only(base_sv, cand_sv):
                tp_base_ctr += 1
                break

            if (cand_sv.sv_pos - cur_end > max_dist_search and base_sv.ref_name == cand_sv.ref_name) or (idx_chr_name(base_sv.ref_name) < idx_chr_name(cand_sv.ref_name)):
                break

        cand_start_idx = new_cand_start_idx
    return tp_base_ctr

###########################################################

def main():
    #hg38 samples
    chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7',
                'chr8','chr9','chr10','chr11','chr12','chr13','chr14',
                'chr15','chr16','chr17','chr18','chr19','chr20','chr21',
                'chr22','chrX']
    
    sv_list = idx_sv(truth_vcf)
    cand_sv_list = idx_sv(cand_vcf)

    sv_list_sorted = sort_sv_list(sv_list)
    cand_sv_list_sorted = sort_sv_list(cand_sv_list)

    tp_base_ctr = count_tp_base_dist_only(sv_list_sorted, cand_sv_list_sorted)
    recall = tp_base_ctr / len(sv_list_sorted)

    for sample_name in sample_names:
        bcf_in_file = sample_name + "_dip_sv_lenfilter_chrfilter_pass.vcf"
    
        asm1_trimmed_file = "assem1_sort_trimmed_coor.bed"
        asm2_trimmed_file = "assem2_sort_trimmed_coor.bed"
        
        filter_dipcall_for_fn_trimcoor(bcf_in_file, asm1_trimmed_file, asm2_trimmed_file, bcf_out_file)

if __name__ == "__main__":
    main()
