#!/usr/bin/python3

#########################################################################################################
# Program: group_seqs_into_samples.py
# Purpose: group sequences based on sample ID (V705_xxxx_xxx_GP|POL|REN_[rpt]) and output each samples's
# sequence fasta file
# Author: Wenjie Deng
# Date: 2022-01-11
#########################################################################################################

import re
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input sequence fasta file")
    args = parser.parse_args()
    infile = args.infile
    name, subject, umi = "", "", ""
    readcount, writecount = 0, 0
    subjectStatus, subjectNames, nameSeq = {}, {}, {}
    subjects = []
    with open(infile, "r") as ifp:
        for line in ifp:
            line = line.strip()
            if line:
                linematch = re.search("^>(\S+)", line)
                if linematch:
                    readcount += 1
                    name = linematch.group(1)
                    if re.search("\|", name):
                        fields = name.split("|")
                        name = fields[-1]
                    namefields = name.split("_")
                    umi = namefields.pop()
                    subject = "_".join(namefields)
                    if subjectStatus.get(subject) is None:
                        subjectNames[subject] = []
                        subjectStatus[subject] = 1
                        subjects.append(subject)
                    subjectNames[subject].append(name)
                    nameSeq[name] = ""
                else:
                    nameSeq[name] += line

    for sbjct in subjects:
        outfile = sbjct + "_pblib.fasta"
        count = 0
        print("\n=== Process sample "+sbjct+" ===")
        with open(outfile, "w") as ofp:
            for name in subjectNames[sbjct]:
                count += 1
                writecount += 1
                ofp.write(">"+name+"\n"+nameSeq[name]+"\n")
        print("** write "+str(count)+" sequences in "+outfile+" **")

    print("\nProcessed total "+str(readcount)+" sequences, write total "+str(writecount)+" sequences into sample files.")
