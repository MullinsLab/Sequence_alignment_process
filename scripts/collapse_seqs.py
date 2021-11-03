#!/usr/bin/python3

# collapse sequences into unique sequences and output a name file showing the relationship between unique sequences and original sequences
# Wenjie Deng
# 2021-10-25

import sys
import re
import argparse

def main(infile, outfile):
    names = []
    nametag, seqname, seq = "", "", ""
    count = 0
    nameseq, seqcount, seqnames = ({} for i in range(3))
    fields = infile.split("/")
    filename = fields[-1]
    match = re.search(r'(.*?).fasta$', filename)
    if match:
        nametag = match.group(1)
    else:
        sys.exit("Not correct fasta file extension, must be '.fasta'")

    with open(infile, "r") as fp:
        for line in fp:
            line = line.strip()
            linematch = re.match(r'^>(\S+)', line)
            if linematch:
                count += 1
                seqname = linematch.group(1)
                names.append(seqname)
                nameseq[seqname] = ""
            else:
                nameseq[seqname] += line.upper()

    for name in names:
        seq = nameseq[name]
        if seqcount.get(seq) is None:
            seqcount[seq] = 0
        seqcount[seq] += 1

        if seqnames.get(seq) is None:
            seqnames[seq] = []
        seqnames[seq].append(name)

    uniqcount = 0
    namefile = outfile.replace("fasta", "name")

    with open(outfile, "w") as out:
        with open(namefile, "w") as fw:
            for seq in sorted(seqcount, key=seqcount.get, reverse=True):
                uniqcount += 1
                name = nametag + "_" + str(uniqcount) + "_" + str(seqcount[seq])
                out.write(">" + name + "\n")
                out.write(seq + "\n")
                fw.write(name + "\t" + ','.join(seqnames[seq]) + "\n")

    log = "processed " + str(count) + " sequences, " + str(uniqcount) + " unique sequences"
    return log

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input sequence fasta file")
    parser.add_argument("outfile", help="output collapsed sequence fasta file")
    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile

    main(infile, outfile)