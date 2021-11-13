#!/usr/bin/python3

##########################################################################################
# Program: V705_POL_misprimes_removal.py
# Purpose: In a POL directory with sequence fasta files,
# adds reference HXB2 POL sequence,
# align sequences and HXB2 via Muscle
# retrieve full length of sequences
# multiprocessing
# Author: Wenjie Deng
# Date: 2021-11-12
##########################################################################################

import sys, re, os
import argparse
import glob
import append_seqs
import muscle_align
import remove_POL_misprimes
from multiprocessing import Pool


def worker(file, outdir, logdir, reffile):
    fields = file.split("/")
    filename = fields[-1]
    filenamefields = filename.split(".")
    sampleid = filenamefields[0]
    sampleoutdir = outdir + "/" + sampleid
    if os.path.isdir(sampleoutdir) is False:
        os.mkdir(sampleoutdir)

    logfilename = filename.replace(".fasta", ".log")
    logfile = logdir + "/" + logfilename
    summaryfile = logfile.replace(".log", "_summary.csv")

    with open(logfile, "w") as lfp:
        print("\n" + "=== Processing file " + file + " ===")
        lfp.write("=== Processing file " + file + " ===" + "\n")
        # append reference (HXB2) sequence
        withreffile = sampleoutdir + "/" + filename.replace(".fasta", "_withHXB2.fasta")
        append_seqs.main(reffile, file, withreffile)
        lfp.write("** Append HXB2 POL sequence to " + file + " **" + "\n")
        lfp.write("input1: " + reffile + "\n")
        lfp.write("input2: " + file + "\n")
        lfp.write("output: " + withreffile + "\n")

        # align sequences
        alignfile = withreffile.replace(".fasta", "_align.fasta")
        muscle_align.main(withreffile, alignfile)
        lfp.write("** Align " + withreffile + " **\n")
        lfp.write("input: " + withreffile + "\n")
        lfp.write("output: " + alignfile + "\n")

        # retrieve full length POL
        fulllenpolfile = sampleoutdir + "/" + filename
        log = remove_POL_misprimes.main(alignfile, fulllenpolfile, summaryfile)
        lfp.write("** Remove POL misprimes in " + alignfile + " **\n")
        lfp.write("input: " + alignfile + "\n")
        lfp.write("output: " + fulllenpolfile + "\n")
        lfp.write("summary: " + summaryfile + "\n")
        lfp.write(log + "\n")

		# output success info
        lfp.write("*** Succeed ***\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", help="directory to hold input sequence fasta file", nargs="?", const=1, type=str, default=".")
    parser.add_argument("-p", "--processes", help="number of processes for multiprocessing", nargs="?", const=1, type=int,
                        default="1")
    args = parser.parse_args()
    dir = args.dir
    proc = args.processes

    outdir = dir+"/POL_misprimes_removal_outputs"
    logdir = dir+"/POL_misprimes_removal_logs"
    if os.path.isdir(outdir) is False:
        os.mkdir(outdir)
    if os.path.isdir(logdir) is False:
        os.mkdir(logdir)
    scriptpath = os.path.dirname(__file__)
    reffile = ''
    if scriptpath == '.':
        reffile = "../references/HXB2_POL.fasta"
    elif re.search("scripts", scriptpath):
        refpath = scriptpath.replace("scripts", "references")
        reffile = refpath + "/HXB2_POL.fasta"
    else:
        reffile = "/opt/V705_alignment_process/references/HXB2_POL.fasta"
    
    files = []
    for file in glob.glob(os.path.join(dir, '*_POL_*.fasta')):
        files.append(file)

    pool = Pool(proc)
    pool.starmap(worker, [(file, outdir, logdir, reffile) for file in files])

    pool.close()
    pool.join()

    alllogfile = logdir + "/POL_misprimes_removal.log"
    allsummaryfile = logdir + "/POL_misprimes_removal_summary.csv"
    with open(alllogfile, "w") as afp:
        with open(allsummaryfile, "w") as sfp:
            sfp.write(
                "file,total,full_length,trim_2908,trim_2993,trim_3008,others" + "\n")
            for file in files:
                fields = file.split("/")
                filename = fields[-1]
                logfilename = filename.replace(".fasta", ".log")
                logfile = logdir + "/" + logfilename
                summaryfilename = filename.replace(".fasta", "_summary.csv")
                summaryfile = logdir + "/" + summaryfilename
                with open(logfile, "r") as lfp:
                    for line in lfp:
                        afp.write(line)
                    afp.write("\n")
                with open(summaryfile, "r") as tfp:
                    for line in tfp:
                        sfp.write(line)


