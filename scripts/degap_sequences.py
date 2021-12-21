#!/usr/bin/python3

# remove gaps in sequences. output degapsed sequence file with the same name of input file in the user named directory (provided by "-d" or "--dir")
# if "-d" or "--dir" is not defined, the output file will be stored in the directory of "degap_sequences" in the working directory by default

import sys, re, os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="input sequence fasta file")
parser.add_argument("-d", "--dir", help="output directory to hold output degaps sequence fasta file", nargs="?", const=1, type=str, default="degap_sequences")
args = parser.parse_args()
infile = args.infile
outdir = args.dir
count = 0
if os.path.isdir(outdir) is False:
    os.mkdir(outdir)
if re.search("/", infile):
    fields = infile.split("/")
    outfile = outdir+"/"+fields[-1]
else:
    outfile = outdir+"/"+infile
print("=== processing "+infile+" ===")
with open(outfile, "w") as ofp:
    with open(infile, "r") as ifp:
        for line in ifp:
            line = line.strip()
            if re.search("^>", line):
                count += 1
                if count > 1:
                    ofp.write("\n")
                ofp.write(line+"\n")
            else:
                line = line.replace("-", "")
                line = line.upper()
                ofp.write(line)
    ofp.write("\n")

print("total "+str(count)+" sequences.")