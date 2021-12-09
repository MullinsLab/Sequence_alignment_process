#!/usr/bin/python

# calculate consensus sequence from an AA sequence alignment (0 majority and ignore gaps)
# majority character (excluding "-") gets consensus.
# If two or more AA share the majority or "-" is majority, "X" will be consensus.
# "-" gets consensus if all gaps at the position

import sys
import re
import os
import argparse
from collections import OrderedDict

def consensus(infile, outfile):
    outconsfile, outaafile, subject, name = "", "", "", ""
    names, consAAs = [], []
    nameSeq, nameAAs = {}, {}
    count, processed, flag, alignlen = 0, 0, 0, 0
    filename = infile
    if re.search("/", infile):
        filefields = infile.split("/")
        filename = filefields[-1]
    infilematch = re.search('^(.*?)_(\d+)_(.*?)_(.*?)_', filename)
    if infilematch:
        subject = infilematch.group(1) + "_" + infilematch.group(2) + "_" + infilematch.group(
            3) + "_" + infilematch.group(4)
    else:
        sys.exit("file name is not formatted: " + infile)

    with open(infile, "r") as ifp:
        for line in ifp:
            line = line.strip()
            linematch = re.search('^>(\S+)', line)
            if linematch:
                name = linematch.group(1)
                namematch = re.search("HXB2", name)
                if namematch:
                    flag = 0
                else:
                    names.append(name)
                    count += 1
                    flag = 1
                    nameSeq[name] = ""
            else:
                if flag == 1:
                    nameSeq[name] += line

    for name in names:
        processed += 1
        seq = nameSeq[name]
        seqlen = len(seq)
        if processed == 1:
            alignlen = seqlen
        else:
            if seqlen != alignlen:
                sys.exit("sequence not aligned: " + name)

        nameAAs[name] = list(seq)

    for i in range(alignlen):
        posAAcount = {}
        for name in names:
            aa = nameAAs[name][i]
            if posAAcount.get(aa) is None:
                posAAcount[aa] = 0
            namematch = re.search('_(\d+)$', name)
            if namematch:
                duplicates = namematch.group(1)
                posAAcount[aa] += int(duplicates)
            else:
                posAAcount[aa] += 1

        if len(posAAcount) == 1:
            consAAs.append(list(posAAcount.keys())[0])
        else:
            sortedPosAAcount = OrderedDict(sorted(posAAcount.items(), key=lambda x: x[1], reverse=True))
            if list(sortedPosAAcount.values())[0] > list(sortedPosAAcount.values())[1]:
                if list(sortedPosAAcount.keys())[0] == '-':
                    consAAs.append('X')
                else:
                    consAAs.append(list(sortedPosAAcount.keys())[0])
            elif list(sortedPosAAcount.values())[0] == list(sortedPosAAcount.values())[1]:
                consAAs.append("X")
            else:
                sys.exit("impossible, something wrong!")

    consname = subject + "_AA_functional_consensus"
    consseq = ''.join(consAAs)

    with open(outfile, "w") as afp:
        afp.write(">" + consname + "\n")
        afp.write(consseq + "\n")
        for name in names:
            afp.write(">" + name + "\n")
            afp.write(''.join(nameAAs[name]) + "\n")

    return consname, consseq

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input AA sequence alignment fasta file")
    parser.add_argument("outfile", help="input AA sequences and consensus alignment fasta file")
    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile

    consensus(infile, outfile)

