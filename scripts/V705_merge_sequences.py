#!/usr/bin/python3

#########################################################################################################
# Program: V705_merge_sequences.py
# Purpose: In a directory with sequence fasta files:
# merge sequences in different time points or repeated sequencing runs based on the sample id (V705_xxxx)
# Author: Wenjie Deng
# Date: 2021-12-15
#########################################################################################################

import sys, re, os
import argparse
import glob
from collections import defaultdict

def nested_dict(n, type):
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", help="directory to hold input sequence fasta file", nargs="?", const=1, type=str, default=".")
    args = parser.parse_args()
    dir = args.dir

    sampleRegionTpStatus = nested_dict(3, int)
    sampleRegionTP = nested_dict(2, list)
    sampleRegionIdStatus = nested_dict(3, int)
    sampleRegionNameSeq = nested_dict(3, str)
    outdir = dir + "/merge_sequences_outputs"
    if os.path.isdir(outdir) is False:
        os.mkdir(outdir)
    for file in glob.glob(os.path.join(dir, '*.fasta')):
        filefields = file.split("/")
        filename = filefields[-1]
        filenamematch = re.search("^(V\d+)_(\d+)_(.*?)_([A-Z]+)_", filename)
        if filenamematch:
            sample = filenamematch.group(1)+"_"+filenamematch.group(2)
            tp = filenamematch.group(3)
            region = filenamematch.group(4)
            if sampleRegionTpStatus[sample][region][tp] == 0:
                sampleRegionTpStatus[sample][region][tp] = 1
                sampleRegionTP[sample][region].append(tp)
            name = ""
            with open(file, "r") as fp:
                for line in fp:
                    line = line.strip()
                    if line:
                        linematch = re.search("^>(\S+)(.*)", line)
                        if linematch:
                            id = linematch.group(1)
                            description = linematch.group(2)
                            # remove viroverse ID if there is in the sequence name
                            if re.search("\|", id):
                                fields = id.split("|")
                                id = fields[1]
                            if re.search("_rpt\d?_", id):
                                id = re.sub("_rpt\d?_", "_", id)
                            if sampleRegionIdStatus[sample][region][id] == 0:
                                sampleRegionIdStatus[sample][region][id] = 1
                                name = id + description
                            else:
                                sampleRegionIdStatus[sample][region][id] += 1
                                name = id + "-" + str(sampleRegionIdStatus[sample][region][id]) + description
                        else:
                            line = line.replace("-", "")
                            sampleRegionNameSeq[sample][region][name] += line
        else:
            sys.exit("Not a right file format: "+filename)

for sample in sorted(sampleRegionTP):
    for rg in sorted(sampleRegionTP[sample]):
        tps = sorted(sampleRegionTP[sample][rg])
        tp = "-".join(tps)
        outfile = outdir + "/" + sample + "_" + tp + "_" + rg +"_pblib.fasta"
        count = 0
        with open(outfile, "w") as ofp:
            for name in sorted(sampleRegionNameSeq[sample][rg]):
                ofp.write(">"+name+"\n"+sampleRegionNameSeq[sample][rg][name]+"\n")
                count += 1
        print("Write merged file "+outfile+", total "+str(count)+" sequences")