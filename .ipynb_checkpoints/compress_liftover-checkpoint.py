import csv
import sys
interval = 20

output_dir = sys.argv[1] + "/"
in_file_name = output_dir + sys.argv[2]
out_file_name = output_dir + sys.argv[3]

#consider reverse strand also
with open(in_file_name, "r") as file:
    reader = csv.reader(file, delimiter="\t")
    start_row = next(reader)
    chr_name = start_row[4]
    contig_name = start_row[0]
    chr_pos = int(start_row[5])
    contig_pos = int(start_row[1])
    pre_chr_pos = chr_pos
    pre_contig_pos = contig_pos
    ctr = 1
    forward_strand = True
    out_file = open(out_file_name, "w")    
    out_file.write(chr_name + "\t" + contig_name + "\t")
    out_file.write(str(chr_pos)+":"+str(contig_pos)+":")

    for row in reader:
        if (row[4] == chr_name) and (row[0] == contig_name) and (int(row[5]) - pre_chr_pos == interval) and \
            (int(row[1]) - pre_contig_pos == interval) and forward_strand:
                pre_chr_pos = int(row[5])
                pre_contig_pos = int(row[1])
                ctr += 1
        elif (row[4] == chr_name) and (row[0] == contig_name) and (int(row[5]) - pre_chr_pos == interval) and \
            (int(row[1]) - pre_contig_pos == -interval) and not forward_strand:
                pre_chr_pos = int(row[5])
                pre_contig_pos = int(row[1])
                ctr += 1             
        else:
            #if different chr or contig
            if row[4] != chr_name or row[0] != contig_name:
                #print(chr_name, contig_name)
                out_file.write(str(ctr)+":")
                if forward_strand:
                    out_file.write("+;")
                else:
                    out_file.write("-;")
                out_file.write("\n")
                chr_name = row[4]
                contig_name = row[0]
                chr_pos = int(row[5])
                contig_pos = int(row[1])
                pre_chr_pos = chr_pos
                pre_contig_pos = contig_pos
                ctr = 1
                forward_strand = True
                out_file.write(chr_name + "\t" + contig_name + "\t")
                out_file.write(str(chr_pos)+":"+str(contig_pos)+":")
            elif int(row[5]) - pre_chr_pos != interval:
                out_file.write(str(ctr)+":")
                if forward_strand:
                    out_file.write("+;")
                else:
                    out_file.write("-;")
                chr_pos = int(row[5])
                contig_pos = int(row[1])
                pre_chr_pos = chr_pos
                pre_contig_pos = contig_pos
                out_file.write(str(chr_pos)+":"+str(contig_pos)+":")
                ctr = 1
                forward_strand = True
            elif (int(row[1]) - pre_contig_pos != interval) and (int(row[1]) - pre_contig_pos != -interval):
                out_file.write(str(ctr)+":")
                if forward_strand:
                    out_file.write("+;")
                else:
                    out_file.write("-;")
                chr_pos = int(row[5])
                contig_pos = int(row[1])
                pre_chr_pos = chr_pos
                pre_contig_pos = contig_pos
                out_file.write(str(chr_pos)+":"+str(contig_pos)+":")
                ctr = 1
                forward_strand = True
            elif (int(row[1]) - pre_contig_pos == interval) and not forward_strand:
                out_file.write(str(ctr)+":")
                out_file.write("-;")
                chr_pos = int(row[5])
                contig_pos = int(row[1])
                pre_chr_pos = chr_pos
                pre_contig_pos = contig_pos
                out_file.write(str(chr_pos)+":"+str(contig_pos)+":")
                ctr = 1
                forward_strand = True
            elif (int(row[1]) - pre_contig_pos == -interval) and forward_strand:
                out_file.write(str(ctr)+":")
                out_file.write("+;")
                chr_pos = int(row[5])
                contig_pos = int(row[1])
                pre_chr_pos = chr_pos
                pre_contig_pos = contig_pos
                out_file.write(str(chr_pos)+":"+str(contig_pos)+":")
                ctr = 1
                forward_strand = False
    #write last lines
    out_file.close()