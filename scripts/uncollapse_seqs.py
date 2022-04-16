#!/usr/bin/python3

# uncollapse unique sequences into original sequences based on name file
# Wenjie Deng
# 2021-10-29

import re
import sys
import argparse

def main(infile, namefile, outfile):
    names = []
    nameseq, uniqnamenames = ({} for i in range(2))
    seqname = ""
    uniqcount, seqcount = 0, 0
    with open(infile, "r") as fp:
        for line in fp:
            line = line.strip()
            linematch = re.match(r'^>(\S+)', line)
            if linematch:
                uniqcount += 1
                seqname = linematch.group(1)
                names.append(seqname)
                nameseq[seqname] = ""
            else:
                nameseq[seqname] += line.upper()

    with open(namefile, "r") as namefp:
        for line in namefp:
            line = line.strip()
            [uniqname, allname] = line.split("\t")
            uniqnamenames[uniqname] = allname.split(",")

    with open(outfile, "w") as outfp:
        for name in names:
            if uniqnamenames.get(name):
                namematch = re.search('_(\d+)$', name)
                if namematch:
                    duplicates = namematch.group(1)
                    arraylen = len(uniqnamenames[name])
                    if (int(duplicates) == arraylen):
                        for originalname in uniqnamenames[name]:
                            seqcount += 1
                            outfp.write(">" + originalname + "\n" + nameseq[name] + "\n")
                    else:
                        sys.exit("duplicates: " + str(duplicates) + " does not equal to array length:" + str(arraylen))
            else:
                seqcount += 1
                outfp.write(">" + name + "\n" + nameseq[name] + "\n")
    log = "processed " + str(uniqcount) + " unique sequences, " + str(seqcount) + " uncollapsed sequences"
    print(log)
    return log

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input collapsed sequence fasta file")
    parser.add_argument("namefile", help="input name file")
    parser.add_argument("outfile", help="output uncollapsed sequence fasta file")
    args = parser.parse_args()
    infile = args.infile
    namefile = args.namefile
    outfile = args.outfile

    main(infile, namefile, outfile)