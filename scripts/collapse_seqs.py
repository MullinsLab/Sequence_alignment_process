#!/usr/bin/python3

# collapse sequences into unique sequences and output a name file showing the relationship between unique sequences and original sequences
# Wenjie Deng
# 2021-10-25

import sys
import re
import argparse

def main(infile, outfile):
    seqname, namefile, id = '', '', ''
    count = 0
    idseqnames, nameseq, seqCount, idseqcount, uniqDup = ({} for i in range(5))
    match = re.match(r'(.*).fasta$', infile)
    if match:
        id = match.group(1)
        namefile = outfile.replace("fasta", "name")
    else:
        sys.exit("Not correct fasta file extension, must be '.fasta'")
    
    with open(infile, "r") as ifp:
        for line in ifp:
            line = line.strip()
            linematch = re.search('^>(\S+)', line)
            if linematch:
                seqname = linematch.group(1)
                if id not in idseqnames:
                    idseqnames[id] = []
                idseqnames[id].append(seqname)
                count += 1
                nameseq[seqname] = ""
            else:
                line = line.upper()
                nameseq[seqname] += line

    for id in idseqnames:
        for name in idseqnames[id]:
            seq = nameseq[name]
            if id not in seqCount:
                seqCount[id] = {}
            if seq not in seqCount[id]:
                seqCount[id][seq] = 0
            seqCount[id][seq] += 1
            if id not in idseqcount:
                idseqcount[id] = 0
            idseqcount[id] += 1
            if id not in uniqDup:
                uniqDup[id] = {}
            if seq not in uniqDup[id]:
                uniqDup[id][seq] = []
            uniqDup[id][seq].append(name)

    with open(outfile, "w") as ofp:
        with open(namefile, "w") as nfp:
            for id in sorted(seqCount):
                uniqcount = 0
                countseqs = {}
                for seq in seqCount[id]:
                    cnt = seqCount[id][seq]
                    if cnt not in countseqs:
                        countseqs[cnt] = []
                    countseqs[cnt].append(seq)

                counts = list(countseqs.keys())
                counts.sort(reverse=True)

                for cnt in (counts):
                    for seq in sorted(countseqs[cnt]):
                        uniqcount += 1
                        name = id+"_"+str(uniqcount)+"_"+str(seqCount[id][seq])
                        ofp.write(">"+name+"\n"+seq+"\n")
                        nfp.write(name + "\t" + ','.join(uniqDup[id][seq]) + "\n")
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