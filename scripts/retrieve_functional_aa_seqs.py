#!/usr/bin/python3

# from a protein sequence alignment, retrieve functional amino acid sequences (stop codon
# within 3' end of 5% of median alignment length and < 20% of deletions of median length
# and output aa alignment without defectives.
# Wenjie Deng
# 2021-10-28

import re
import argparse

def main(ntfile, aafile, outfile, gene):
    geneProtein = {
        "gag": "Gag",
        "pol": "Pol",
        "env": "Env"
    }
    functionalaafile = aafile.replace(".fasta", "_functional.fasta")
    defectiveaafile = aafile.replace(".fasta", "_defective.fasta")
    functionalntfile = ntfile.replace(".fasta", "_functional.fasta")
    defevtiventfile = ntfile.replace(".fasta", "_defective.fasta")
    name = ""
    names, alignlens, lens = [], [], []
    count, seqcount, hxb2alignlen, maxalignlen, missing5count, notMcount, functionalcount, deletioncount, prematurecount = 0, 0, 0, 0, 0, 0, 0, 0, 0
    nameAaseq, nameNtseq, nameStatus = ({} for i in range(3))
    with open(aafile, "r") as ifp:
        for line in ifp:
            line = line.strip()
            linematch = re.search("^>(\S+)", line)
            if linematch:
                name = linematch.group(1)
                names.append(name)
                nameAaseq[name] = ""
                count += 1
            else:
                nameAaseq[name] += line

    with open(ntfile, "r") as nfp:
        for line in nfp:
            line = line.strip()
            linematch = re.search("^>(\S+)", line)
            if linematch:
                name = linematch.group(1)
                pname = name.replace("_"+gene+"_", "_"+geneProtein[gene]+"_")
                nameNtseq[pname] = ""
            else:
                nameNtseq[pname] += line

    for name in names:
        aaseq = nameAaseq[name]
        if re.search("HXB2", name):
            hxb2alignlen = len(aaseq)
        else:
            seqcount += 1
            alignlen = len(aaseq)
            if alignlen > maxalignlen:
                maxalignlen = alignlen
            alignlens.append(alignlen)
            seq = aaseq.replace("-", "")
            lens.append(len(seq))
    if hxb2alignlen > maxalignlen:
        maxalignlen = hxb2alignlen

    sortedlens = sorted(lens)
    sortedalignlens = sorted(alignlens)
    mediancount = int(float(seqcount) / 2)
    medianlen = sortedlens[mediancount]
    medianalignlen = sortedalignlens[mediancount]

    for name in names:
        seq = nameAaseq[name]
        seqlen = len(seq)
        gapstripseq = seq.replace("-", "")
        if re.search("HXB2", name):
            nameStatus[name] = 1
            for i in range(seqlen, maxalignlen):
                nameAaseq[name] += "-"
        else:
            if re.search("^-", seq):
                missing5count += 1
            elif (gene != "pol" and re.search("^M", seq) is None):
                notMcount += 1
            elif seqlen >= medianalignlen * 0.95:
                if len(gapstripseq) >= medianlen * 0.8:
                    nameStatus[name] = 1
                    functionalcount += 1
                    for i in range(seqlen, maxalignlen):
                        nameAaseq[name] += "-"
                else:
                    deletioncount += 1
            else:
                prematurecount += 1

    with open(functionalaafile, "w") as fafp:
        with open(defectiveaafile, "w") as dafp:
            with open(functionalntfile, "w") as fnfp:
                with open(defevtiventfile, "w") as dnfp:
                    for name in names:
                        aaseq = nameAaseq[name]
                        ntseq = nameNtseq[name]
                        gname = name.replace("_"+geneProtein[gene]+"_", "_"+gene+"_")
                        if nameStatus.get(name):
                            fafp.write(">" + name + "\n" + aaseq + "\n")
                            fnfp.write(">" + gname + "\n" + ntseq + "\n")
                        else:
                            dafp.write(">" + name + "\n" + aaseq + "\n")
                            dnfp.write(">" + gname + "\n" + ntseq + "\n")

    with open(outfile, "w") as ofp:
        ofp.write(aafile + "," + geneProtein[gene] + "," + str(seqcount) + "," + str(functionalcount) + "," + str(
            prematurecount) + "," + str(deletioncount)
                  + "," + str(missing5count) + "," + str(notMcount) + "," + str(medianlen) + "\n")

    log = "Processing " + str(seqcount) + " sequences, " + str(functionalcount) + " functional, " + str(
        missing5count) + " missing 5' end, " + str(notMcount) + " not start with 'M', " + str(deletioncount) + " with big deletion (> 20%), " + str(
        prematurecount) + " premature stop (< 95%), median length: " + str(medianlen)
    return log

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("ntfile", help="input NT fasta file")
    parser.add_argument("aafile", help="input AA fasta file")
    parser.add_argument("outfile", help="output summary csv file")
    parser.add_argument("gene", help="gene name")
    args = parser.parse_args()
    ntfile = args.ntfile
    aafile = args.aafile
    outfile = args.outfile
    gene = args.gene

    main(ntfile, aafile, outfile, gene)
