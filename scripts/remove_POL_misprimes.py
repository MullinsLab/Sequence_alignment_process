#!/usr/bin/python3

# remove misprimed POL sequences
# Wenjie Deng
# 2021-11-12

import sys
import re
import argparse

def main(infile, outfile, summaryfile):
    count, fullcount, count1, count2, count3, count4 = 0, 0, 0, 0, 0, 0
    names = []
    nameSeq = {}
    refname = ""
    with open(infile, "r") as ifp:
        for line in ifp:
            line = line.strip()
            linematch = re.search(">(\S+)", line)
            if linematch:
                name = linematch.group(1)
                nameSeq[name] = ""
                namematch = re.search("HXB2", name)
                if namematch:
                    refname = name
                else:
                    names.append(name)
                    count += 1
            else:
                nameSeq[name] += line

    ntidx, idx1, idx2, idx3 = 0, 0 ,0, 0
    refnts = list(nameSeq[refname])
    for i in range(len(refnts)):
        refnt = refnts[i]
        ntmatch = re.match("[A-Z]", refnt)
        if ntmatch:
            ntidx += 1
            if ntidx == 2908:
                idx1 = i
            elif ntidx == 2993:
                idx2 = i
            elif ntidx == 3008:
                idx3 = i
                break

    with open(outfile, "a") as ofp:
        for name in names:
            seq = nameSeq[name]
            seqmatch = re.search("[A-Z]$", seq)
            if seqmatch:
                fullcount += 1
                seq = seq.replace("-", "")
                ofp.write(">"+name+"\n"+seq+"\n")
            else:
                nts = list(seq)
                if (re.search("[A-Z]", nts[idx1]) and re.search("^-+$", seq[idx1+1:])):
                    count1 += 1
                elif (re.search("[A-Z]", nts[idx2]) and re.search("^-+$", seq[idx2+1:])):
                    count2 += 1
                elif (re.search("[A-Z]", nts[idx3]) and re.search("^-+$", seq[idx3+1:])):
                    count3 += 1
                else:
                    count4 += 1

    with open(summaryfile, "w") as sfp:
        sfp.write(infile+","+str(count)+","+str(fullcount)+","+str(count1)+","+str(count2)+","+str(count3)+","+str(count4)+"\n")

    log = "Processed " + str(count) + " sequences, " + str(fullcount) + " full length of POL, " + str(
        count1) + " trimmed at 2908, " + str(count2) + " trimmed at 2993, " + str(
        count3) + " trimmed at 3008, " + str(count4) + " trimmed at other sites"
    return log


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input sequence alignment with HXB2 POL fasta file")
    parser.add_argument("outfile", help="output POL misprimes removed fasta file")
    parser.add_argument("summaryfile", help="output summary csv file")
    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile
    summaryfile = args.summaryfile

    main(infile, outfile, summaryfile)