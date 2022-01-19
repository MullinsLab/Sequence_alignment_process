#!/usr/bin/python3

##########################################################################################
# Program: V705_alignment_process.py
# Purpose: In a GP/POL/REN directory with sequence fasta files,
# collapses sequences into unique sequences,
# adds reference HXB2 GP/POL/REN sequence,
# left-align unique sequences and HXB2 via Muscle
# uncollapse and verify sequences
# extract alignment into region of gag/pol/env
# translate NA sequence alignment into AA sequence alignment
# retrieve functional protein sequence alignment
# collapse functional protein sequence alignment
# Author: Wenjie Deng
# Date: 2021-10-29
# Modified: multiprocessing
# Date: 2021-11-03
# Modified: change the structure of output files based on sample id
# Date: 2021-11-05
##########################################################################################

import sys, re, os
import argparse
import glob
from argparse import ArgumentParser

import collapse_seqs
import append_seqs
import reverse_seq
import muscle_align
import uncollapse_seqs
import verify_seq_origin
import extract_alignment_portion
import refine_pol_start_homopolymerT
import ntAlignment2aaAlignment
import retrieve_functional_aa_seqs
import strip_all_gaps
from multiprocessing import Pool
import shutil


def worker(file, outdir, logdir, refpath):
    regionGene = {
        "GP": "Gag",
        "POL": "Pol",
        "REN": "Env"
    }
    geneStart = {
        "Gag": 86,
        "Pol": 247,
        "Env": 243
    }
    geneEnd = {
        "Gag": 1588,
        "Pol": 3258,
        "Env": 2813
    }
    fields = file.split("/")
    filename = fields[-1]
    filenamefields = filename.split(".")
    sampleid = filenamefields[0]
    sampleoutdir = outdir + "/" + sampleid
    if os.path.isdir(sampleoutdir) is False:
        os.mkdir(sampleoutdir)
    region = ''
    if re.search("_GP_", file):
        region = "GP"
    elif re.search("_POL_", file):
        region = "POL"
    elif re.search("_REN_", file):
        region = "REN"
    else:
        sys.exit("Error: there is no region of GP/POL/REN in sequence file name: " + file)
    logfilename = filename.replace(".fasta", ".log")
    logfile = logdir + "/" + logfilename
    tallyfilename = filename.replace(".fasta", "_functional_tally.csv")
    tallyfile = logdir + "/" + tallyfilename

    with open(logfile, "w") as lfp:
        # collapse sequences
        collapsedfile = sampleoutdir + "/" + filename.replace(".fasta", "_collapsed.fasta")
        print("\n" + "=== Processing file " + file + " ===")
        collapselog = collapse_seqs.main(file, collapsedfile)
        lfp.write("=== Processing file " + file + " ===" + "\n")
        lfp.write("** Collapse sequences in " + file + " **" + "\n")
        lfp.write("input: " + file + "\n")
        lfp.write("output: " + collapsedfile + "\n")
        lfp.write(collapselog + "\n")

        # append reference (HXB2) sequence
        reffile = refpath + "/HXB2_" + region + ".fasta"
        withreffile = collapsedfile.replace(".fasta", "_withRef.fasta")
        append_seqs.main(reffile, collapsedfile, withreffile)
        lfp.write("** Append reference sequence to " + collapsedfile + " **" + "\n")
        lfp.write("input1: " + reffile + "\n")
        lfp.write("input2: " + collapsedfile + "\n")
        lfp.write("output: " + withreffile + "\n")

        # reverse sequences
        reversedfile = withreffile.replace(".fasta", "_rvs.fasta")
        reverse_seq.main(withreffile, reversedfile)
        lfp.write("** reverse sequences in " + withreffile + " **" + "\n")
        lfp.write("input: " + withreffile + "\n")
        lfp.write("output: " + reversedfile + "\n")

        # align reversed sequences
        reversedalignfile = reversedfile.replace(".fasta", "_align.fasta")
        muscle_align.main(reversedfile, reversedalignfile)
        lfp.write("** Align collapsed file " + reversedfile + " **\n")
        lfp.write("input: " + reversedfile + "\n")
        lfp.write("output: " + reversedalignfile + "\n")

        # reverse back to make a left-aligned alignment
        withrefalignfile = collapsedfile.replace(".fasta", "_withRef_align.fasta")
        reverse_seq.main(reversedalignfile, withrefalignfile)
        lfp.write("** reverse sequences in " + reversedalignfile + " **" + "\n")
        lfp.write("input: " + reversedalignfile + "\n")
        lfp.write("output: " + withrefalignfile + "\n")

        # uncollapse alignment
        uncollapsealignfile = withrefalignfile.replace("_collapsed", "")
        namefile = collapsedfile.replace(".fasta", ".name")
        uncollapselog = uncollapse_seqs.main(withrefalignfile, namefile, uncollapsealignfile)
        lfp.write("** Uncollapse alignment " + withrefalignfile + " **" + "\n")
        lfp.write("input: " + withrefalignfile + "\n")
        lfp.write("namefile: " + namefile + "\n")
        lfp.write("output: " + uncollapsealignfile + "\n")
        lfp.write(uncollapselog + "\n")

        # verify sequences
        verifylog = verify_seq_origin.main(uncollapsealignfile, file)
        lfp.write("** Verify sequences between " + uncollapsealignfile + " **" + "\n")
        lfp.write("input1: " + uncollapsealignfile + "\n")
        lfp.write("input2: " + file + "\n")
        lfp.write(verifylog + "\n")

        # retrieve gag/pol/env
        gene = regionGene[region]
        genedir = sampleoutdir + "/" + gene
        if os.path.isdir(genedir) is False:
            os.mkdir(genedir)
        sgene = geneStart[gene]
        egene = geneEnd[gene]
        uncollapsealignfilefields = uncollapsealignfile.split("/")
        uncollapsealignfilename = uncollapsealignfilefields[-1]
        genefilename = uncollapsealignfilename.replace("_" + region + "_", "_" + gene + "_NT_")
        genefilename = genefilename.replace("_align.", ".")
        genefile = genedir + "/" + genefilename
        extract_alignment_portion.main(uncollapsealignfile, genefile, region, gene, sgene, egene)
        lfp.write("** Extract " + gene + " at HXB2 position of " + str(sgene) + " to " + str(egene) + " **" + "\n")
        lfp.write("input: " + uncollapsealignfile + "\n")
        lfp.write("output: " + genefile + "\n")
        lfp.write("region: " + region + "\n")
        lfp.write("gene: " + gene + "\n")
        lfp.write("gene start: " + str(sgene) + "\n")
        lfp.write("gene end: " + str(egene) + "\n")

        if re.search("_Pol_", genefile):
            # refine Pol sequences at the beginning of 6 homopolymer Ts
            originalpolfile = genefile.replace(".fasta", "_origin.fasta")
            refinedpolfile = genefile.replace(".fasta", "_refined.fasta")
            shutil.copyfile(genefile, originalpolfile)
            lfp.write("** Refine Pol alignment for the beginning homopolymer T insertion in " + originalpolfile + " **" + "\n")
            homopolymerTinscount = refine_pol_start_homopolymerT.main(originalpolfile, refinedpolfile)
            if homopolymerTinscount:
                lfp.write("input: " + originalpolfile + "\n")
                lfp.write("output: " + refinedpolfile + "\n")
                lfp.write(str(homopolymerTinscount)+" sequences with beginning homopolymer T insertions\n")

                # strip all gap columns
                striplog = strip_all_gaps.main(refinedpolfile, genefile, False)
                lfp.write("** Strip all gap columns in " + refinedpolfile + " **" + "\n")
                lfp.write("input: " + refinedpolfile + "\n")
                lfp.write("output: " + genefile + "\n")
                lfp.write(striplog + "\n")
            else:
                lfp.write("No sequence with beginning homopolymer T insertions\n")

        # translation
        geneaafile = genefile.replace("_NT_", "_AA_")
        ntAlignment2aaAlignment.main(genefile, geneaafile)
        lfp.write("** Translate nucleotide to amino acid sequences for " + genefile + " **" + "\n")
        lfp.write("input: " + genefile + "\n")
        lfp.write("output: " + geneaafile + "\n")

        #retrieve functional protein sequences and write summary
        funclog = retrieve_functional_aa_seqs.main(geneaafile, tallyfile, gene)
        lfp.write("** Retrieve functional amino acid sequences in " + geneaafile + " **" + "\n")
        lfp.write("input: " + geneaafile + "\n")
        lfp.write("output: " + tallyfile + "\n")
        lfp.write("gene: " + gene + "\n")
        lfp.write(funclog + "\n")

        # remove HXB2 and strip all gap columns
        functionalaafile = geneaafile.replace(".fasta", "_functional.fasta")
        gapstripfile = functionalaafile.replace("_withRef_", "_")
        striplog = strip_all_gaps.main(functionalaafile, gapstripfile, True)
        lfp.write("** Remove HXB2 and strip all gap columns in " + functionalaafile + " **" + "\n")
        lfp.write("input: " + functionalaafile + "\n")
        lfp.write("output: " + gapstripfile + "\n")
        lfp.write(striplog + "\n")

        #collapse functional protein sequence alignment
        collapsedgapstripfile = gapstripfile.replace(".fasta", "_collapsed.fasta")
        colpslog = collapse_seqs.main(gapstripfile, collapsedgapstripfile)
        lfp.write("** Collapse sequences in " + gapstripfile + " **" + "\n")
        lfp.write("input: " + gapstripfile + "\n")
        lfp.write("output: " + collapsedgapstripfile + "\n")
        lfp.write(colpslog + "\n")

		# output success info
        lfp.write("*** Succeed ***\n")

if __name__ == '__main__':
    parser: ArgumentParser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", help="directory to hold input sequence fasta file", nargs="?", const=1, type=str, default=".")
    parser.add_argument("-p", "--processes", help="number of processes for multiprocessing", nargs="?", const=1, type=int,
                        default="1")
    args = parser.parse_args()
    dir = args.dir
    proc = args.processes

    outdir = dir+"/alignment_process_outputs"
    logdir = dir+"/alignment_process_logs"
    if os.path.isdir(outdir) is False:
        os.mkdir(outdir)
    if os.path.isdir(logdir) is False:
        os.mkdir(logdir)
    scriptpath = os.path.dirname(__file__)
    refpath = ''
    if scriptpath == '.':
    	refpath = "../references"
    elif re.search("scripts", scriptpath):
    	refpath = scriptpath.replace("scripts", "references")
    else:
        refpath = "/opt/V705_alignment_process/references"

    files = []
    for file in glob.glob(os.path.join(dir, '*.fasta')):
        files.append(file)

    pool = Pool(proc)
    pool.starmap(worker, [(file, outdir, logdir, refpath) for file in files])

    pool.close()
    pool.join()

    alllogfile = logdir + "/alignment_process.log"
    summaryfile = logdir + "/functional_protein_summary.csv"
    with open(alllogfile, "w") as afp:
        with open(summaryfile, "w") as sfp:
            sfp.write(
                "file,gene,sequences,functional,premature_stop,big_deletion,missing_5',not_start_M,median_sequence_length" + "\n")
            for file in files:
                fields = file.split("/")
                filename = fields[-1]
                logfilename = filename.replace(".fasta", ".log")
                logfile = logdir + "/" + logfilename
                tallyfilename = filename.replace(".fasta", "_functional_tally.csv")
                tallyfile = logdir + "/" + tallyfilename
                with open(logfile, "r") as lfp:
                    for line in lfp:
                        afp.write(line)
                    afp.write("\n")
                with open(tallyfile, "r") as tfp:
                    for line in tfp:
                        sfp.write(line)


