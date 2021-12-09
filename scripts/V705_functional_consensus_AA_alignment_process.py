#!/usr/bin/python3

##########################################################################################
# Program: V705_functional_consensus_AA_alignment_process.py
# Purpose: In a Gag/Pol/Env directory with AA sequence alignment fasta files:
# calculate AA consensus via 0 majority and ignoring gaps,
# align AA consensus sequences via Muscle,
# profile-profile align reference AA alignment with consensus alignment
# multiprocessing
# Author: Wenjie Deng
# Date: 2021-12-08
##########################################################################################

import sys, re, os
import argparse
import glob
import calc_AA_consensus
import muscle_align
import profile_align
from multiprocessing import Pool, Manager

def worker(file, outdir, logdir, nameSeq):
    fields = file.split("/")
    filename = fields[-1]
    logfilename = filename.replace(".fasta", ".log")
    logfile = logdir + "/" + logfilename

    with open(logfile, "w") as lfp:
        # calculate functional AA consensus
        sampleAAwithConsfile = outdir + "/" + filename.replace(".fasta", "_withCons.fasta")
        print("=== Calculate consensus sequences in " + file + " ===\n")
        name, seq = calc_AA_consensus.consensus(file, sampleAAwithConsfile)
        nameSeq[name] = seq
        lfp.write("=== Calculate consensus sequences in " + file + " ===" + "\n")
        lfp.write("input: " + file + "\n")
        lfp.write("output: " + sampleAAwithConsfile + "\n")

        # output success info
        lfp.write("*** Succeed ***\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("refalignfile", help="reference sequence alignment fasta file")
    parser.add_argument("-d", "--dir", help="directory to hold input sequence fasta file", nargs="?", const=1, type=str, default=".")
    parser.add_argument("-p", "--processes", help="number of processes for multiprocessing", nargs="?", const=1, type=int,
                        default="1")
    args = parser.parse_args()
    reffile = args.refalignfile
    dir = args.dir
    proc = args.processes

    outdir = dir + "/functional_consensus_AA_alignment_process_outputs"
    logdir = dir + "/functional_consensus_AA_alignment_process_logs"
    if os.path.isdir(outdir) is False:
        os.mkdir(outdir)
    if os.path.isdir(logdir) is False:
        os.mkdir(logdir)
    files = []
    for file in glob.glob(os.path.join(dir, '*.fasta')):
        files.append(file)

    manager = Manager()
    consnameSeq = manager.dict()
    pool = Pool(proc)
    pool.starmap(worker, [(file, outdir, logdir, consnameSeq) for file in files])

    pool.close()
    pool.join()

    filefields = files[0].split("/")
    filename = filefields[-1]
    filenamefields = filename.split("_")
    id = filenamefields[0] + "_" + filenamefields[3]
    alllogfile = logdir + "/" + id +"_functional_consensus_AA_alignment_process.log"
    consdir = outdir + "/functional_consensus"
    consensusfile = consdir + "/" + id +"_AA_functional_consensus.fasta"
    print("=== Write consensus sequences into " + consensusfile + " ===\n")
    if os.path.isdir(consdir) is False:
        os.mkdir(consdir)
    with open(consensusfile, "w") as cfp:
        for name in consnameSeq:
            cfp.write(">"+name+"\n"+consnameSeq[name]+"\n")
    with open(alllogfile, "w") as afp:
        for file in files:
            fields = file.split("/")
            filename = fields[-1]
            logfilename = filename.replace(".fasta", ".log")
            logfile = logdir + "/" + logfilename
            with open(logfile, "r") as lfp:
                for line in lfp:
                    afp.write(line)
                afp.write("\n")
        afp.write("=== Write consensus sequences into " + consensusfile + " ===\n\n")

        # align consensus sequences
        print("=== Align consensus sequences in "+consensusfile+" ===\n")
        consensusalignfile = consensusfile.replace(".fasta", "_align.fasta")
        muscle_align.main(consensusfile, consensusalignfile)
        afp.write("=== Align consensus sequences in "+consensusfile+" ===\n")
        afp.write("input: " + consensusfile + "\n")
        afp.write("output: " + consensusalignfile + "\n\n")

        # align reference alignment with consensus alignment via profile-profile align
        print("=== Align " + reffile + " and " + consensusfile + " ===\n")
        outfile = consensusalignfile.replace("_align.fasta", "_withRef_align.fasta")
        profile_align.profile_align(reffile, consensusalignfile, outfile)
        afp.write("=== Align " + reffile + " and " + consensusfile + " ===\n")
        afp.write("input1: " + reffile + "\n")
        afp.write("input2: " + consensusalignfile + "\n")
        afp.write("output: " + outfile + "\n")
