import re
import os
import glob

#quality control function
def checkSeq(s):
    import re
    s = s.strip()
    s = s.upper()
    global clean
    clean = 0

    #check if multiple of 3
    if len(s) % 3 != 0:
        clean = 1

    #check no non ATGC
    if clean == 0:
        notATCG = re.sub("[ATCG]","",s)
        if len(notATCG) > 0:
            clean = 1

    #check starts with ATG
    #don't do this for bacterial genomes as many start NTG
    #if clean == 0:
        #first_cdn = s[0:3]
        #if first_cdn != "ATG":
            #clean = 1

    #check if ends in a stop
    stops = ["TAA", "TAG", "TGA"]
    if clean == 0:
        last_cdn = s[-3:]
        if last_cdn not in stops:
            clean = 1

    #check for internal stops
    stops = ["TAA", "TAG", "TGA"]
    if clean == 0:
        cdns = re.findall('[A-Za-z]{3}', s[:-3])
        for stp in stops:
            if stp in cdns:
                clean = 1

    return clean


#get GC information function
def GCcount(sequence):
    nucs = "ACTG"
    GC = "GC"

    count_nucs = {}
    for n in nucs:
        count_nucs[n] = sequence.count(n)

    total = 0
    for i in nucs:
        total = total + count_nucs[i]

    tot_GC = 0
    for j in GC:
        tot_GC = tot_GC + count_nucs[j]

    GC_count = round(tot_GC / total, 2)
    GC_perc = GC_count * 100
    return GC_perc
    #print(GC_perc)

def GC3count(sequence):
    if len(sequence) % 3 != 0:
        raise ValueError("The length of the sequence must be a multiple of 3.")

    third_positions = sequence[2::3]
    GC3_count = third_positions.count('G') + third_positions.count('C')
    #print(GC3_count)

    GC3_perc = round((GC3_count / len(third_positions)) * 100,2) if len(third_positions) > 0 else 0
    return GC3_perc


cwd = os.getcwd()

os.chdir(f"{cwd}/all_genomes")

allfiles = glob.glob("*.fna")

for f in allfiles:
    genus = re.search("(.*?).fna", f)
    genus = genus.group()
    genus = re.sub(".fna", "", genus)

    os.chdir(cwd)

    if not os.path.exists(f"{cwd}/all_seqs"):
        os.makedirs(f"{cwd}/all_seqs")
    os.chdir(f"{cwd}/all_seqs")

    outfile_five = open(f"{genus}_five.csv", "x")
    outfile_five.write("gene_name,prot_id,number_of_exons,five_prime_seq\n")

    outfile_core = open(f"{genus}_core.csv", "x")
    outfile_core.write("gene_name,prot_id,number_of_exons,gen_core_seq\n")

    outfile_gc = open(f"{genus}_GC.csv", "x")
    outfile_gc.write("gene_name,prot_id,number_of_exons,GC,GC3\n")

    os.chdir(f"{cwd}/all_genomes")

    infile = open(f, "r")
    all_genes = infile.read()
    infile.close()


    genelist = all_genes.split(">l") #some entries contain "<" or ">" between exons
    genelist = genelist[1:]

    genelist_len = len(genelist)

    gene_name_list = []

    cnt_gn = 0
    cnt_gn_good = 0

    for gn in genelist:
        cnt_gn = cnt_gn + 1
        #if cnt_gn % 100 == 0:
            #print(f"Analysing gene {cnt_gn}")

        # split entry name and sequence into separate lines
        g_all = gn.splitlines(True)
        top_line = g_all[0]

        # each gene sequence should be a single line
        seq = g_all[1:]
        seq = re.sub("[^A-Za-z]", "", "".join(seq))
        seq = seq.upper()

        # quality control for sequence
        checkSeq(seq)
        if clean == 0:

            # get gene name from entry name (top line)
            if re.search("\[gene=.*?\]", top_line):
                gene = re.search("\[gene=(.*?)\]", top_line)
                gene_name = gene.group(1)
            else:
                gene_name = "NA"

            # get protein id from entry name (top line)
            if re.search("protein_id=.*?]", top_line):
                prot = re.search("protein_id=(.*?)]", top_line)
                prot_id = prot.group(1)
            else:
                prot_id = "NA"


            # get exon number:
            # get exon locations from entry name (top line)
            location = re.search("location=(.*?)\]", top_line)
            locs = location.group(1)

            if re.search("[<>]", locs):
                print(f'Rubbish annotation: {locs}')
                continue

            num_exons = locs.count(",") + 1
            num_dots = locs.count("..")


            if num_exons != num_dots:
                print(f'Gene = {gene_name}, number {cnt_gn}: exon structure is odd, {locs}')
                continue


            if gene_name != "NA":
                if gene_name in gene_name_list:
                    #print(f"Already analysed {gene_name}")
                    continue
                else:
                    gene_name_list.append(gene_name)
            #cnt_gn_good = cnt_gn_good + 1


            # get seq at 5' end and genetic core
            seq_split = []

            if len(seq) < 70: #check if entry seq is long enough - go higher? the literature uses 50-codons long seq
                print(f'Gene = {gene_name}, number {cnt_gn}, seq length {len(seq)}: sequence is too short')
                continue
            else:
                for i in range(0, len(seq), 3):
                    seq_split.append(seq[i:i+3])

                seq_five = seq_split[0:11]
                seq_five = "".join(seq_five)

                g_core = (len(seq_split)-1)//2
                seq_core = seq_split[g_core:(g_core+11)]
                seq_core = "".join(seq_core)

                outfile_five.write(f'{gene_name},{prot_id},{num_exons},{seq_five}\n')
                outfile_core.write(f'{gene_name},{prot_id},{num_exons},{seq_core}\n')

            # get GC information
            CDS = seq[30:]
            GC_perc = GCcount(CDS)
            GC3_perc = GC3count(CDS)

            outfile_gc.write(f'{gene_name},{prot_id},{num_exons},{GC_perc},{GC3_perc}\n')

            cnt_gn_good = cnt_gn_good + 1

    outfile_five.close()
    outfile_core.close()
    outfile_gc.close()

    print(f'For {genus}: of {genelist_len}, {cnt_gn_good} are decent')



