#!/usr/bin/python3

# profile-profile alignment via Muscle
# Wenjie Deng
# 2021-11-10

import os
import argparse

def profile_align(infile1, infile2, outfile):
    os.system("muscle -quiet -profile -in1 "+infile1+" -in2 "+infile2+" -out "+outfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile1", help="input reference sequence alignment fasta file")
    parser.add_argument("infile2", help="input sequence alignment fasta file")
    parser.add_argument("outfile", help="output aligned sequence fasta file")
    args = parser.parse_args()
    infile1 = args.infile1
    infile2 = args.infile2
    outfile = args.outfile

    profile_align(infile1, infile2, outfile)