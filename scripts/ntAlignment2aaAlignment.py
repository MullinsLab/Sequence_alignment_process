#!/usr/bin/python3

# translates any nucleotide alignment into amino acid alignment
# xxxx------xxxxx --> X--XX
# xxxxx------xxxx --> XX--X
# Wenjie Deng
# 2021-10-29

import re
import argparse

def main(infile, outfile, gene, protein):
    name, refname, refseq = '', '', ''
    nameSeq = {}
    names = []
    count = 0
    codon2aa = {
        "---": "-",
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        "TTA": "L",
        "TTG": "L",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "TTT": "F",
        "TTC": "F",
        "ATG": "M",
        "TGT": "C",
        "TGC": "C",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "AGT": "S",
        "AGC": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TGG": "W",
        "CAA": "Q",
        "CAG": "Q",
        "AAT": "N",
        "AAC": "N",
        "CAT": "H",
        "CAC": "H",
        "GAA": "E",
        "GAG": "E",
        "GAT": "D",
        "GAC": "D",
        "AAA": "K",
        "AAG": "K",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "AGA": "R",
        "AGG": "R",
        "TAA": "*",
        "TAG": "*",
        "TGA": "*",
    }

    with open(infile, "r") as ifp:
        for line in ifp:
            line = line.strip()
            linematch = re.search(">(\S+)", line)
            if linematch:
                name = linematch.group(1)
                aaname = name.replace("_"+gene+"_", "_"+protein+"_")
                names.append(aaname)
                nameSeq[aaname] = ""
                count += 1
            else:
                nameSeq[aaname] += line.upper()

    with open(outfile, "w") as ofp:
        for name in names:
            seq = nameSeq[name]
            seqlen = len(seq)
            nts = list(seq)
            aaseq, partialaaseq, gaps, gapflag = '', '', '', ''
            codonnts, codongaps = [], []
            for i in range(seqlen):
                if nts[i] == "-":
                    codongaps.append("-")
                elif re.match("[ACGT]", nts[i]):
                    codonnts.append(nts[i])

                if len(codongaps) == 3:
                    codongaps = []
                    gaps += "-"
                    if len(codonnts) == 0:
                        aaseq += gaps
                        gaps, gapflag, partialaaseq = "", "", ""
                    elif len(codonnts) == 1:
                        gapflag = "head"
                    elif len(codonnts) == 2:
                        gapflag = "tail"

                if len(codonnts) == 3:
                    codon = codonnts.pop(0)
                    codon += codonnts.pop(0)
                    codon += codonnts.pop(0)
                    aa = codon2aa[codon]
                    if gapflag == "head":
                        partialaaseq = gaps + aa
                    elif gapflag == "tail":
                        partialaaseq = aa + gaps
                    else:
                        partialaaseq = aa
                    aaseq += partialaaseq
                    gaps, gapflag, partialaaseq = "", "", ""
                    if aa == "*":
                        break
            ofp.write(">" + name + "\n" + aaseq + "\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input nuclotide sequence alignment fasta file")
    parser.add_argument("outfile", help="output amino acid sequence alignment fasta file")
    parser.add_argument("gene", help="gene name")
    parser.add_argument("protein", help="protein name")
    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile
    gene = args.gene
    protein = args.protein

    main(infile, outfile, gene, protein)

