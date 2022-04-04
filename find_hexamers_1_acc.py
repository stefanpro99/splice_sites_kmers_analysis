import sys
import os
import subprocess
import time

"""

This code takes SpliSER combine output to detect all specified k-mers on the intronic/exonic side of given splice sites
and outputs their frequency across all splice-sites and their average splice site strength.
It updates the k-mer reference sequences based on the reference genome fasta file and variants in the vcf file.  

"""

start_time = time.time() # start measuring run time

file_sites = sys.argv[1]  
# file containing list of splice-sites to be analysed.
# Format: [chromosome]_[*]/SpliSER.[gene]_[position].tsv
# can be obtained by using ls 

vcf_file = sys.argv[2]
# vcf file with genomic variants, positions should match the fasta file

sites_dir = sys.argv[3]
# directory in which the desired SpliSER combine outputs are located

FastaFile = sys.argv[4]
# fasta file (.fa) of reference genome with matching positions to the supplied VCF file

output_file = sys.argv[5]
# prefix of the output files

K = int(sys.argv[6])
# kmer length

output_file_2 = open(output_file + "_avged.tsv", "w+")
# open file to store all analysed sites info

MAX_DEL_SIZE = 20
RANGE = K + MAX_DEL_SIZE
REF_AL = "0/0"
ALT_AL = "1/1"
TARGET_ecoID = 100
GENDERS = ["male", "female"]



def acc_with_allele_set(allele, head, line):
    # returns the set containing accessions harbouring the given allele in a given VCF line, based on given VCF header 
    if type(line) != list:
        line = line.strip().split("\t")
    return {vcf_ecoID_to_num(head[index])
            for index, value in enumerate(line)
            if value.split(":")[0] == allele and index >= 9}


def vcf_ecoID_to_num(accession):
    # takes accession id from the VCF encoding and returns it as integer 
    number = accession.split("-")[1]
    if number[0] == "0":
        number = number[1:]
    return int(number)

def get_site(line, ls = True):
    # returns chromosome and site position from the input line
    # example input format: /home/.../SpliSER_output/3R_1_950/SpliSER.FBgn0000003_6822509.tsv
    # set ls = False if sites list is a tsv 
    if ls:
        chr = line.strip().split("/")[-2].split("_")[0]
        site = line.strip().split("_")[-1].split(".")[0]
    else:
        chr,site = line.strip().split("\t")
    return chr, int(site)

positions_dict = {}

with open(file_sites, "r") as f:  # read splice sites from the input file as keys in positions_dict 
    f.readline()
    for line in f:

        chr, site = get_site(line)
        for offset in range(RANGE * -1, RANGE + 1):
            positions_dict[chr + "_" + str(site + offset)] = []

with open(vcf_file, "r") as vcf_file:  # read the vcf file
    line = vcf_file.readline()
    while "#" in line:  ## skip the # lines and store header
        head = line
        line = vcf_file.readline()
    head = head.strip().split("\t")

    accession_index = head.index(
        "DGRP-" + str(TARGET_ecoID))  # get corresponding index of the target accession in vcf file

    for line in vcf_file:  # for variant in vcf file
        line = line.strip().split("\t")
        if line[accession_index].split(":")[0] == ALT_AL:  # if target accession has alternative allele
            chr, pos, id, ref, alt, qual, fliter, info, format = line[:9]

            if "DEL" in id:  # if it's a deletion                            # if variant is a deletion
                deletion_len = len(ref) - 1  # get size of the deletion
                if deletion_len > MAX_DEL_SIZE:  # if it's longer than 20bp
                    id = "/"  # / = mark for unusable all positions covering deleted space
                for pos_del in range(int(pos), int(pos) + deletion_len):
                    if "_".join([chr, str(pos_del)]) in positions_dict:  # add deletion info to all positions under it

                        positions_dict["_".join([chr, str(pos_del)])] = id  # example id: 2L_2262_DEL_TTC_T   or / 

            elif "_".join(
                    [chr, pos]) in positions_dict:  # if insertion/SNP is around some splice-site -> store variant info
                positions_dict["_".join([chr, pos])] = id  # example id: 2L_2262_INS_T_TTC  or 2L_2262_SNP_C_G 


# example dictionary in positions_dict
# { 3R_123 : [[variant_id, {ref}, {alt}], .. [all other variants at position]], ... for all positions in range around splice sites}


def get_sequence(fastaFile, chrom, position, site_seq_size=30, bp_offset=0):
    # Get the sequence in range (+- 30bp default) from the given from a fasta file with option to use offset.
    # requires loading/installing module samtools on the commandline

    bp_seq = None
    bamview = subprocess.Popen(['samtools', 'faidx', str(fastaFile),
                                str(chrom) + ':' + str(int(position) - site_seq_size + bp_offset) + '-' + str(
                                    position + site_seq_size + bp_offset + 1)], stdout=subprocess.PIPE)

    bamstream = bamview.stdout
    for i, line in enumerate(bamstream):
        dline = line.decode('ascii')
        if i > 1:  # get multiple lines
            bp_seq = bp_seq + dline.rstrip()
        else:
            bp_seq = dline.rstrip()
    return bp_seq


def substitute_allele(ref_seq, site_pos, alt_info):
    # update reference sequence with alternative variants (SNPs and INDELs)
    # can take string or list, returns a list
    if type(ref_seq) != list:
        ref_seq = [bp for bp in ref_seq]  # -30 to +30 idx(pos) = 30
    chr, var_pos, var_type, ref_allele, alt_allele = alt_info.split("_")
    var_idx = int(var_pos) - int(site_pos) + RANGE
    if var_type == "DEL":
        for i in range(max(0, var_idx), min(var_idx + len(ref_allele), len(ref_seq))):
            ref_seq[i] = ""
    ref_seq[var_idx] = alt_allele
    return ref_seq


def left_right(pos, partner):
    # evaluates if splice site is on left or right side of intron from forward strand perspective
    # returns side and updated position of splice site adjusted for SpliSER position detection of splice junctions.
    if int(pos) > int(partner):
        side = "right"
        offset = 0
    else:
        side = "left"
        offset = 1
    return side, str(int(pos) + offset)


def complemet(seq):
    # returns reverse complement of a DNA sequence
    out = ""
    dict = {"A": "T", "G": "C", "T": "A", "C": "G"}
    for base in seq:
        out = dict[base] + out
    return out


def donor_acceptor(side, strand):
    # function that evaluates if splice_site is a donor or acceptor based on input side and strand
    if side == "left":
        if strand == "+":
            return "donor"
        else:
            return "acceptor"
    elif side == "right":
        if strand == "+":
            return "acceptor"
        else:
            return "donor"


def get_kmer(ref_seq, Side, range=30, k=6):
    # takes the reference (str) or updated (list) sequence and returns the k-mer (6 default) on the intronic side from splice site
    kmer_seq = ""
    i = range
    if Side == "left":
        while len(kmer_seq) < k:  # keep merging bases (0-k) from ref_seq
            kmer_seq += ref_seq[i]
            i += 1
        if len(kmer_seq) > k:
            kmer_seq = kmer_seq[:k]  # trim if there are insertions
    elif Side == "right":  # if right side
        while len(kmer_seq) < k:  # keep merging bases (0-k) in reverse order from ref_seq
            kmer_seq = ref_seq[i] + kmer_seq
            i -= 1
        if len(kmer_seq) > k:
            kmer_seq = kmer_seq[:k]  # trim if there are insertions
    return kmer_seq


def get_exon_kmer(ref_seq, Side, range=30, k=6):
    # takes the reference or updated (list) sequence and returns the desired k-mer (6 default) on the exonic side from splice site
    kmer_seq = ""

    if Side == "right":
        i = range + 1
        while len(kmer_seq) < k:  # keep merging bases (0-k) from ref_seq
            kmer_seq += ref_seq[i]
            i += 1
        if len(kmer_seq) > k:
            kmer_seq = kmer_seq[:k]  # trim if there are insertions
    elif Side == "left":  # if right side
        i = range - 1
        while len(kmer_seq) < k:  # keep merging bases (0-k) in reverse order from ref_seq
            kmer_seq = ref_seq[i] + kmer_seq
            i -= 1
        if len(kmer_seq) > k:
            kmer_seq = kmer_seq[:k]  # trim if there are insertions
    return kmer_seq


def record_kmers(samples_dict,kmers_dict, GENDERS, output_file_2, sites_2bp ):
    # update kmer info within a dictionary
    # write into a file (output_file_2) extra data about the detected site and kmer
    for gender in GENDERS:
        counts, SSE_sum, prints = samples_dict[gender]
        if counts > 0:
            Region, Site, Strand, Gene, Sample, SSE, site_type, two_bp, kmer_seq, exon_hexa_seq, Partners, Competitors, variants_list = prints

            if str(two_bp) == sites_2bp[site_type]:
                avg_SSE = SSE_sum / counts
                output_file_2.write("\t".join([str(x) for x in (
                    Region, Site, Strand, Gene, Sample, avg_SSE, site_type, two_bp, kmer_seq,
                    exon_hexa_seq, Partners, Competitors, variants_list)]) + "\n")

                try:  # example { 'GTACGG' : [5.9, 6] ... }
                    kmers_dict[gender][site_type][kmer_seq][0] += float(avg_SSE)         # add SSE to the total sum
                    kmers_dict[gender][site_type][kmer_seq][1] += 1                      # increase count of sites in which the kmer is detected
                except KeyError:
                    kmers_dict[gender][site_type][kmer_seq] = [float(avg_SSE), 1]



count_outputs = 0
sites_2bp = {"donor": "GT", "acceptor": "AG"}
kmers_dict = {gender: {"donor": {}, "acceptor": {}} for gender in GENDERS}

for file in os.listdir(sites_dir):
    if ".tsv" in str(file):
        # print(str(file), time.time() - start_time)
        count_outputs += 1
        with open(str(sites_dir) + file, "r") as f:
            f.readline()
            current_site = ""                           # keep track of splice site as multiple are stored in one file
            count_samples = 0
            for line in f:
                # read information about each site detected by SpliSER
                Sample, Region, Site, Strand, Gene, SSE, alpha_count, beta1_count, beta2Simple_count, beta2Cryptic_count, beta2_weighted, Partners, Competitors = line.strip().split(
                    "\t")
                srr, ecoID, gender = Sample.split("_")
                if Site != current_site:           # if new site is seen, output info for previously seen replicates and reset any counts
                    if count_samples > 0 :         # if any samples passed the filter for the target accession

                        record_kmers(samples_dict,kmers_dict, GENDERS, output_file_2, sites_2bp )  # record data for the detected samples

                    samples_dict = {gender: [0, 0, []] for gender in GENDERS}          # reset samples info for new site
                    count_samples = 0
                    current_site = Site

                if ecoID == str(TARGET_ecoID):  # if target accession
                    if int(alpha_count) + int(beta1_count) + int(
                            beta2Simple_count) > 10:  # only consider splice-sites in samples with > 10 read coverage

                        Side, Site = left_right(Site, Partners.split(":")[0][
                                                      1:])  # compute the side of the splice site and update position accordingly
                        site_type = donor_acceptor(Side, Strand)  # compute if site is donor or acceptor
                        ref_seq = get_sequence(FastaFile, Region, int(Site), site_seq_size=RANGE,
                                               bp_offset=0)  # get reference sequence from a fasta file via samtools
                        Key = Region + "_" + Site
                        avoid = False
                        variants_list = []
                        for offset in range(RANGE * -1, RANGE):  # for positions -30 to +30 bp from site
                            Site_nearby = int(Site) + offset
                            Key = Region + "_" + str(Site_nearby)
                            if Key in positions_dict:  # if there is a variant (SNP/INDEL) at position
                                alt_info = positions_dict[Key]

                                if alt_info and alt_info != "/":  # if not deletion > 20bp

                                    chr, var_pos, var_type, ref_allele, alt_allele = alt_info.split("_")

                                    variants_list.append(alt_info)
                                    ref_seq = substitute_allele(ref_seq, Site,
                                                                alt_info)  # update ref seq for alt alleles
                                    if ref_seq[RANGE] == "":
                                        avoid = True
                                        break

                                elif alt_info == "/":  # if deletion > 20bp around site
                                    avoid = True
                                    break

                        if not avoid:  # if site is usable -> store its kmer data 
                            kmer_seq = get_kmer(ref_seq, Side, range= RANGE, k = K)           
                            exon_hexa_seq = get_exon_kmer(ref_seq, Side, range = RANGE, k= K)

                            if Strand == "-":                  # reverse complement kmers on - strand
                                kmer_seq = complemet(kmer_seq)
                                exon_hexa_seq = complemet(exon_hexa_seq)

                            if site_type == "donor":         # get the splice site dimer sequence (GT/AG if cannonical)
                                two_bp = kmer_seq[:2]
                            else:                            # if acceptor
                                two_bp = kmer_seq[-2:]

                            prints = (
                                Region, Site, Strand, Gene, Sample, SSE, site_type, two_bp, kmer_seq, exon_hexa_seq,
                                Partners, Competitors, variants_list)
                            #print("\t".join([str(x) for x in prints]))
                            count_samples += 1          # store data for replicates
                            samples_dict[gender][0] += 1
                            samples_dict[gender][1] += float(SSE)
                            samples_dict[gender][2] = prints

                elif int(ecoID) > TARGET_ecoID:     # if accession is after the target accession 
                                                    # store replicate info (if any)
                                                    
                    if count_samples > 0:  # if any samples with enough reads for the target accession

                        record_kmers(samples_dict,kmers_dict, GENDERS, output_file_2, sites_2bp)  # update info for previously detected samples

                        samples_dict = {gender: [0, 0, []] for gender in GENDERS}  # reset samples info and samples count
                        count_samples = 0

            if count_samples > 0 :
                record_kmers(samples_dict, kmers_dict, GENDERS, output_file_2, sites_2bp ) # update info for the last splice site in the file

for gender in kmers_dict.keys():  # write output files male/female if gendered and acceptor/donor
    for site_type in kmers_dict[gender].keys():
        with open(output_file + "_" + site_type + "_" + gender + ".txt", "w+") as f:
            for key in kmers_dict[gender][site_type].keys():
                SSE_total, count = kmers_dict[gender][site_type][key]
                SSE_avg = SSE_total / count
                f.write("\t".join([key, str(SSE_avg), str(count)]) + "\n")

print("done in ", (time.time() - start_time) / 60, "minutes")    # output run time

