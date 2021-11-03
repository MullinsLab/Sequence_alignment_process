#!/usr/bin/python3

# extract alignment portion based on the coordinates of reference sequence
# Wenjie Deng
# 2021-10-29

import re
import argparse

def main(infile, outfile, region, gene, start, end):
    name, refname, refseq = '', '', ''
    nameSeq = {}
    names = []
    with open(infile, "r") as ifp:
        for line in ifp:
            line = line.strip()
            linematch = re.search(">(\S+)", line)
            if linematch:
                name = linematch.group(1)
                name = name.replace("_" + region, "_" + gene)
                namematch = re.search("HXB2", name)
                if namematch:
                    refname = name
                else:
                    names.append(name)
                nameSeq[name] = ""
            else:
                nameSeq[name] += line.upper()

    refseq = nameSeq[refname]
    refnts = list(refseq)
    alignlen = len(refnts)
    refseqnogaps = refseq.replace("-", "")
    reflen = len(refseqnogaps)
    if end == 0:
        end = reflen
    idx, startidx, endidx = 0, 0, 0
    for i in range(alignlen):
        if re.match("[A-Z]", refnts[i]):
            idx += 1
            if idx == start:
                startidx = i
            elif idx == end:
                endidx = i
                break

    with open(outfile, "w") as ofp:
        refextract = refseq[startidx:endidx + 1]
        ofp.write(">" + refname + "\n" + refextract + "\n")
        for name in names:
            seq = nameSeq[name]
            seqextract = seq[startidx:endidx + 1]
            ofp.write(">" + name + "\n" + seqextract + "\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infasta", help="input sequence alignment fasta file")
    parser.add_argument("outfile", help="output extracted sequence alignment fasta file")
    parser.add_argument("region", help="amplicon region")
    parser.add_argument("gene", help="extracted gene name")
    parser.add_argument("-s", "--start", help="extract start position in HXB2 reference sequence (default: 1, beginning of HXB2)", default=1, nargs='?', const=1, type=int)
    parser.add_argument("-e", "--end", help="extract end position in HXB2 reference sequence (default: 0, end of HXB2)", default=0, nargs='?', const=1, type=int)
    args = parser.parse_args()
    infile = args.infasta
    outfile = args.outfile
    region = args.region
    gene = args.gene
    start = args.start
    end = args.end

    main(infile, outfile, region, gene, start, end)


