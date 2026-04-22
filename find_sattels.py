import sys
import getopt
import time
import pysam
from collections import defaultdict
import regex as re

def write_fastq(read, file_name):
    with open(file_name, 'a', encoding="utf-8") as file:
        file.write("@"+read.query_name+"\n")
        file.write(read.query_sequence+"\n")
        file.write("+"+"\n")
        file.write("F"*read.query_length+"\n")

def find_motif_sat(read : str, min_num : int, max_mismatches : int) -> bool:
    motif = "AAACTAGACAGAAGCAT"
    match = re.findall(re.compile(r'('+motif+'){s<='+str(max_mismatches)+'}'), read)
    if(len(match) >= min_num):
        return True
    else:
       return False
    
def find_motif_tel_fw(read : str, min_num : int) -> bool:
    motif = "TTAGGG"
    pattern = r'('+motif+'){'+str(min_num)+',100}'
    match = re.search(re.compile(pattern), read)
    if match:
        return True
    else:
        return False
    
def find_motif_tel_rev(read : str, min_num : int) -> bool:
    motif = "CCCTAA"
    pattern = r'('+motif+'){'+str(min_num)+',100}'
    match = re.search(re.compile(pattern), read)
    if match:
        return True
    else:
        return False


def parse_bam_test(file_in, file_out, file_tmp_bam) -> list:
    n_reads = 0
    n_mol = 0
    n_sat_mol = 0
    n_tel_mol = 0
    n_sattel_mol = 0

    size = 120000000
    count_sat = [0] * size
    print("list size: " + str(sys.getsizeof(count_sat) / 1024 / 1024))
    file_out = open("test_tmp.txt", 'w')
    for i in range(size):
       if count_sat[i] == 0:
           count_sat[i] += 1
       #file_out.write(f"{i}\t{count_sat[i]}\n")
    file_out.close()  

    return n_mol, n_reads, n_sat_mol, n_tel_mol, n_sattel_mol



#no tel mols x
#no sat mols x
#no sattel mols x
#sattel mols position
#no tel reads x
#no sat reads x
#no sattel reads x
#combination of both x
#leave out duplicates is_duplicate x

def parse_bam(file_in, file_out, file_sat_reads_bam, file_tel_reads_bam, min_num_tel, min_num_sat, max_mismatches_sat) -> list:

    min_reads_per_mol = 2
    min_tel_per_mol = 2
    min_sat_per_mol = 1
    min_sattel_per_mol = 1

    n_reads_total = 0
    n_reads_barcoded_total = 0
    n_reads_barcoded_unique = 0
    n_tel_reads = 0
    n_sat_reads = 0
    n_sattel_reads = 0
    n_mol = 0
    n_sat_mol = 0
    n_tel_mol = 0
    n_sattel_mol = 0
    n_sattel_reads_mol = 0
    max_n_mol = 120000000 # depends on maximum possible number of molecules (identifiers); ca. 920MB RAM needed per list below
    count_reads_per_mol = [0] * max_n_mol
    count_sat_per_mol = [0] * max_n_mol
    count_tel_per_mol = [0] * max_n_mol
    count_sattel_per_mol = [0] * max_n_mol
 
    bamfile = pysam.AlignmentFile(file_in, "rb")
    bamfile_tmp_sat = pysam.AlignmentFile(file_sat_reads_bam, "wb", template=bamfile)
    bamfile_tmp_tel = pysam.AlignmentFile(file_tel_reads_bam, "wb", template=bamfile)

    # scan read sequences for satellite and telomere repeat motifs
    for read in bamfile.fetch():
        n_reads_total += 1
        try:
            mi = read.get_tag("MI")
            n_reads_barcoded_total += 1
            if not read.is_duplicate:
                n_reads_barcoded_unique += 1
                seq = read.query_sequence
                #seq2 = "AAAAAAAAAATTAGGGTTAGGGTTAGGGTTAGGGTTTTTTTTTTAAACTAGACAGAAGCATAAACTAGACAGAAGCATAAAAAAAAAA"
                is_sat = find_motif_sat(seq, min_num_sat, max_mismatches_sat)
                is_tel = find_motif_tel_fw(seq, min_num_tel) or find_motif_tel_rev(seq, min_num_tel)
                n_sat_reads += is_sat
                n_tel_reads += is_tel
                n_sattel_reads += is_sat and is_tel
                count_reads_per_mol[mi] += 1
                count_sat_per_mol[mi] += is_sat
                count_tel_per_mol[mi] += is_tel
                count_sattel_per_mol[mi] += is_tel and is_sat
                # write telomeric and satellite reads into temporariy bam files, so we don't have to scan the whole input bam file again to add their position to sattel molecules (see below)
                if read.is_mapped:
                    if is_sat:
                        bamfile_tmp_sat.write(read)
                    if is_tel:
                        bamfile_tmp_tel.write(read)
        except KeyError:
            pass

    bamfile.close()
    bamfile_tmp_sat.close()
    bamfile_tmp_tel.close()

    pysam.index(file_sat_reads_bam)
    pysam.index(file_tel_reads_bam)

    # add read alignment positions to sattel molecules
    # satellite reads
    bamfile_tmp = pysam.AlignmentFile(file_sat_reads_bam, "rb")
    pos_sat_mol = {}
    for read in bamfile_tmp.fetch():
        if not read.is_duplicate:
            try:
                mi = read.get_tag("MI")
                if count_sat_per_mol[mi] >= min_sat_per_mol and count_tel_per_mol[mi] >= min_tel_per_mol:
                    if mi not in pos_sat_mol:
                        pos_sat_mol[mi] = []
                    pos_sat_mol[mi].append(f"{read.reference_name}:{read.reference_start}")        
            except KeyError:
                pass
    bamfile_tmp.close()
    #telomere reads
    bamfile_tmp = pysam.AlignmentFile(file_tel_reads_bam, "rb")
    pos_tel_mol = {}
    for read in bamfile_tmp.fetch():
        if not read.is_duplicate:
            try:
                mi = read.get_tag("MI")
                if count_sat_per_mol[mi] >= min_sat_per_mol and count_tel_per_mol[mi] >= min_tel_per_mol:
                    if mi not in pos_tel_mol:
                        pos_tel_mol[mi] = []
                    pos_tel_mol[mi].append(f"{read.reference_name}:{read.reference_start}")        
            except KeyError:
                pass
    bamfile_tmp.close()

    # summarize and write out results
    file = open(file_out, 'w')
    file.write("molecule_id\tn_reads\tn_tel\tn_sat\treads_sat_aln\treads_tel_aln\n")

    for i in range(max_n_mol):
        # only use and count molecules with at least two reads
        if count_reads_per_mol[i] >= min_reads_per_mol:
            n_mol += 1
            n_sat_mol += count_sat_per_mol[i] >= min_sat_per_mol
            n_tel_mol += count_tel_per_mol[i] >= min_tel_per_mol
            n_sattel_reads_mol += count_sattel_per_mol[i] >= min_sattel_per_mol and count_tel_per_mol[i] >= min_tel_per_mol

            # in case of a sattel molecule: write out detailed molecule information
            if count_sat_per_mol[i] >= min_sat_per_mol and count_tel_per_mol[i] >= min_tel_per_mol:
                n_sattel_mol += 1
                pos_sat = ",".join(pos_sat_mol[i])
                pos_tel = ",".join(pos_tel_mol[i])
                file.write(f"{i}\t{count_reads_per_mol[i]}\t{count_tel_per_mol[i]}\t{count_sat_per_mol[i]}\t{pos_sat}\t{pos_tel}\n")

    file.close()

    return n_mol, n_reads_total, n_reads_barcoded_total, n_reads_barcoded_unique, n_tel_reads, n_sat_reads, n_sattel_reads, n_sat_mol, n_tel_mol, n_sattel_mol, n_sattel_reads_mol


def usage():
    print("Here will soon be a help page")

def main():

    min_num_tel = 4 # default
    min_num_sat = 1 # default
    max_mm_sat = 2 # default
    
    script = sys.argv[0]
    argv = sys.argv[1:]
    short_args = "i:o:t:s:m:h"
    long_args = ["input_bam=", "out_prefix=", "min_num_tel=", "min_num_sat=", "max_mm_sat=", "help"]

    try:
        opts, args = getopt.getopt(argv, short_args, long_args)
    except getopt.GetoptError as err:
        print("Error: Invalid arguments")
        print(str(err))
        sys.exit(2)

    for name, value in opts:
        if name in ["-i", "--input_bam"]:
            file_in = value
        elif name in ["-o", "--out_prefix"]:
            file_out_prefix = value
        elif name in ["-t", "--min_num_tel"]:
            min_num_tel = int(value)
        elif name in ["-s", "--min_num_sat"]:
            min_num_sat = int(value)
        elif name in ["-m", "--max_mm_sat"]:
            max_mm_sat = int(value)
        elif name in ["-h", "--help"]:
            usage()
            sys.exit(0)


    #sample_name = sys.argv[1] #"test"
    #file_in = sys.argv[2] #"/Users/cbartenh/Programs/telomere_satellites/test.bam"
    #file_out_prefix = sys.argv[3] #"/Users/cbartenh/Programs/telomere_satellites/test"
    file_mol = file_out_prefix + "_molecule_counts.txt"     
    file_summary = file_out_prefix + "_summary.txt"
    file_sat_reads_bam = file_out_prefix + "_sat.bam"     
    file_tel_reads_bam = file_out_prefix + "_tel.bam"     

    time_start = time.time()
    n_mol, n_reads_total, n_reads_barcoded_total, n_reads_barcoded_unique, n_tel_reads, n_sat_reads, n_sattel_reads, n_sat_mol, n_tel_mol, n_sattel_mol, n_sattel_reads_mol = parse_bam(file_in, file_mol, file_sat_reads_bam, file_tel_reads_bam, min_num_tel, min_num_sat, max_mm_sat)
    time_end = time.time()
    print(f"Done ({time_end-time_start:.2f}s)")

    print(f"Number of total reads: {n_reads_total}")
    print(f"Number of total barcoded reads: {n_reads_barcoded_total}")
    print(f"Number of unique barcoded reads: {n_reads_barcoded_unique}")
    print(f"Number of unique telomeric reads: {n_tel_reads}")
    print(f"Number of unique satellite reads: {n_sat_reads}")
    print(f"Number of unique satellite-telomere split-reads: {n_sattel_reads}")
    print(f"Number of molecules: {n_mol}")
    print(f"Number of telomere molecules: {n_tel_mol}")
    print(f"Number of satellite molecules: {n_sat_mol}")
    print(f"Number of satellite-telomere molecules: {n_sattel_mol}")
    print(f"Number of satellite-telomere molecules with satellite-telomere reads: {n_sattel_reads_mol}")

    with open(file_summary, 'w', encoding="utf-8") as file:
        file.write(f"call: {script} -i {file_in} -o {file_out_prefix} -t {min_num_tel} -s {min_num_sat} -m {max_mm_sat}\n")
        file.write("n_reads_total\tn_reads_barcoded_total\tn_reads_barcoded_unique\tn_reads_barcoded_unique_tel\tn_reads_barcoded_unique_sat\tn_reads_barcoded_unique_sat_tel\tn_molecules_total\tn_molecules_tel\tn_molecules_sat\tn_molecules_sat_tel\tn_molecules_sat_tel_reads\n")
        file.write(f"{n_reads_total}\t{n_reads_barcoded_total}\t{n_reads_barcoded_unique}\t{n_tel_reads}\t{n_sat_reads}\t{n_sattel_reads}\t{n_mol}\t{n_tel_mol}\t{n_sat_mol}\t{n_sattel_mol}\t{n_sattel_reads_mol}\n")

if __name__ == "__main__":
    main()

