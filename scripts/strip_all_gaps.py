#!/usr/bin/python3

# Purpose: parse the fasta alignment file, remove the columns which all gaps
# Input: sequence alignment fasta file
# Output: cleaned sequence alignment file
# Author: Wenjie Deng
# Date: 2021-10-29

import sys
import re
import argparse

def main(infile, outfile):
    names = []
    name = ""
    count, rmcount = 0, 0
    nameSeq, nameNts, removeSite = ({} for i in range(3))
    with open(infile, "r") as ifp:
        for line in ifp:
            line = line.strip()
            linematch = re.search("^>(\S+)", line)
            if linematch:
                name = linematch.group(1)
                nameSeq[name] = ""
                namematch = re.search("HXB2", name)
                if namematch:
                    refname = name
                else:
                    count += 1
                    names.append(name)
            else:
                nameSeq[name] += line

    alignlen = len(nameSeq[refname])
    for name in names:
        seq = nameSeq[name]
        if len(seq) != alignlen:
            sys.exit("sequence not aligned")
        nameNts[name] = list(seq)

    for i in range(alignlen):
        gapcount = 0
        for name in names:
            if nameNts[name][i] == "-":
                gapcount += 1
        if gapcount == count:
            removeSite[i] = 1
            rmcount += 1

    with open(outfile, "w") as ofp:
        for name in names:
            ofp.write(">" + name + "\n")
            for i in range(alignlen):
                if removeSite.get(i) is None:
                    ofp.write(nameNts[name][i])
            ofp.write("\n")
    log = "total " + str(count) + " sequences, alignment length " + str(alignlen) + ", " + str(
        rmcount) + " all gap columns removed"
    return log

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input sequence alignment fasta file")
    parser.add_argument("outfile", help="output gap stripped sequence alignment fasta file")
    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile

    main(infile, outfile)
