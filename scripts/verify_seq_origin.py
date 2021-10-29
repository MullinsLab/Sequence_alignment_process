#!/usr/bin/python

# verify the sequences between two files
# Wenjie Deng
# 2021-10-29

import re
import argparse

def main(infile, originfile):
    print("input1: " + infile)
    print("input2: " + originfile)
    name, originname = "", ""
    count, origincount = 0, 0
    nameseq, originnameseq = ({} for i in range(2))

    with open(originfile, "r") as ofp:
        for line in ofp:
            line = line.strip()
            linematch = re.search('^>(\S+)', line)
            if linematch:
                origincount += 1
                name = linematch.group(1)
                originnameseq[name] = ""
            else:
                line = line.upper()
                line = line.replace("-", "")
                originnameseq[name] += line

    with open(infile, "r") as ifp:
        for line in ifp:
            line = line.strip()
            linematch = re.search('^>(\S+)', line)
            if linematch:
                name = linematch.group(1)
                if re.search("HXB2", name) is None:
                    count += 1
                nameseq[name] = ""
            else:
                line = line.upper()
                line = line.replace("-", "")
                nameseq[name] += line

    if (origincount != count):
        print("number of sequences are not same. count " + str(count) + ", original count " + str(origincount))

    exactcount, containcount, origincontaincount, mismatchcount, missingcount = 0, 0, 0, 0, 0
    for originname in originnameseq:
        originalseq = originnameseq[originname]
        if nameseq.get(originname) is None:
            missingcount += 1
            print(originname + " missing in " + infile + " file")
        else:
            seq = nameseq[originname]
            if (seq == originalseq):
                exactcount += 1
            else:
                originsearch = re.search(seq, originalseq)
                if originsearch:
                    origincontaincount += 1
                else:
                    seqsearch = re.search(originalseq, seq)
                    if seqsearch:
                        containcount += 1
                    else:
                        mismatchcount += 1
                        print("sequence " + originname + " are not same, original: " + originalseq + ", seq: " + seq)

    print("total " + str(count) + " sequences, " + str(origincount) + " original sequences, " + str(
        exactcount) + " exact match, " + str(origincontaincount) + " contained in " + originfile + ", " + str(
        containcount) + " contained in " + infile + ", " + str(missingcount) + " missed in " + infile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input sequence fasta file")
    parser.add_argument("originfile", help="input sequence fasta file to compare")
    args = parser.parse_args()
    infile = args.infile
    originfile = args.originfile

    main(infile, originfile)
