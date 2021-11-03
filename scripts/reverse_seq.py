#!/usr/bin/python3

# reverse sequences in a sequence fasta file
# Wenjie Deng
# 2021-10-29

import re
import argparse

def main(infile, outfile):
    name = ""
    names = []
    nameSeq = {}
    with open(infile, "r") as ifp:
        for line in ifp:
            line = line.strip()
            if re.search("^>", line):
                name = line
                names.append(name)
                nameSeq[name] = ""
            else:
                nameSeq[name] += line

    with open(outfile, "w") as ofp:
        for name in names:
            ofp.write(name + "\n" + nameSeq[name][::-1] + "\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input sequence fasta file")
    parser.add_argument("outfile", help="output reversed sequence fasta file")
    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile

    main(infile, outfile)