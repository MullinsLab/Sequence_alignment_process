#!/usr/bin/python3

# rearrange homopolymer Ts at the beginning of pol to avoid the frame shift caused by the homopolymer T insertions due
# to the left aligning
# Wenjie Deng
# 2022-01-14

import re
import argparse

def main(infile, outfile):
    name, refname, refseq = '', '', ''
    nameSeq, homopolymerTinsStatus = {}, {}
    names, gapidxs = [], []
    count, tlen, gapflag, homopolymerTinsseqcount = 0, 0, 0, 0
    idx = -1
    with open(infile, "r") as ifp:
        for line in ifp:
            line = line.strip()
            linematch = re.search(">(\S+)", line)
            if linematch:
                name = linematch.group(1)
                names.append(name)
                nameSeq[name] = ""
                count += 1
                if re.search("HXB2", name):
                    refname = name
            else:
                nameSeq[name] += line.upper()

    hxb2nts = list(nameSeq[refname])
    for nt in hxb2nts:
        idx += 1
        if nt == "A":
            break
        if gapflag:
            gapidxs.append(idx)
        else:
            if nt == "T":
                tlen += 1
                if tlen == 6:
                    gapflag = 1

    if gapidxs:
        nameNts = {}
        alen = gapidxs[-1] + 1
        for name in names:
            nameNts[name] = list(nameSeq[name])
            for idx in gapidxs:
                if nameNts[name][idx] != "-":
                    seq = ""
                    for i in range(alen):
                        seq += nameNts[name][i]
                    seq = seq.replace("-", "")
                    if re.search("^T+$", seq) and len(seq) > 6:
                        nameNts[name][idx] = "-"
                        if homopolymerTinsStatus.get(name) is None:
                            homopolymerTinsStatus[name] = 1
                            homopolymerTinsseqcount += 1

    if homopolymerTinsseqcount:
        with open(outfile, "w") as ofp:
            for name in names:
                seq = "".join(nameNts[name])
                ofp.write(">"+name+"\n"+seq+"\n")

    return homopolymerTinsseqcount

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input original pol sequence alignment fasta file")
    parser.add_argument("outfile", help="output refined pol sequence alignment fasta file")
    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile

    main(infile, outfile)

