#!/usr/bin/python

# append sequences
# Wenjie Deng
# 2021-10-29

import argparse

def main(infile1, infile2, outfile):
    print("input1: " + infile1)
    print("input2: " + infile2)
    print("output: " + outfile)
    with open(outfile, "w") as ofp:
        with open(infile1, "r") as i1fp:
            for line in i1fp:
                line = line.strip()
                ofp.write(line+"\n")
        with open(infile2, "r") as i2fp:
            for line in i2fp:
                line = line.strip()
                ofp.write(line+"\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile1", help="input sequence fasta file1")
    parser.add_argument("infile2", help="input sequence fasta file2")
    parser.add_argument("outfile", help="output combined sequence fasta file")
    args = parser.parse_args()
    infile1 = args.infile1
    infile2 = args.infile2
    outfile = args.outfile

    main(infile1, infile2, outfile)