#!/usr/bin/python

##########################################################################################
# Program: V705_alignment_process.pl
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
##########################################################################################

import sys, re, os
import argparse
import glob
import collapse_seqs
import append_seqs
import reverse_seq
import muscle_align
import uncollapse_seqs
import verify_seq_origin
import extract_alignment_portion
import ntAlignment2aaAlignment
import retrieve_functional_aa_seqs
import strip_all_gaps

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", help="directory to hold input sequence fasta file", nargs="?", const=1, type=str, default=".")
args = parser.parse_args()
dir = args.dir
outdir = "alignment_process_output"
if os.path.isdir(outdir) is False:
    os.mkdir(outdir)
scriptpath = os.path.dirname(__file__)
refpath = scriptpath.replace("scripts", "references")
workdir = os.getcwd()
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
filecount = 0
for file in glob.glob(os.path.join(dir, '*.fasta')):
    fields = file.split("/")
    filename = fields[-1]
    filecount += 1
    region = ""
    if re.search("_GP_", file):
        region = "GP"
    elif re.search("_POL_", file):
        region = "POL"
    elif re.search("_REN_", file):
        region = "REN"
    else:
        sys.exit("Error: there is no region of GP/POL/REN in sequence file name: " + file)

    # collapse sequences
    collapsedfile = outdir + "/" + filename.replace(".fasta", "_collapsed.fasta")
    print("\n"+"=== Processing file "+file+" ===")
    print("== Collapse sequences in "+file+" ==")
    collapse_seqs.main(file, collapsedfile)

    # append reference (HXB2) sequence
    print("== Append reference sequence to "+collapsedfile+" ==")
    reffile = refpath + "/HXB2_" + region + ".fasta"
    withreffile = collapsedfile.replace(".fasta", "_withRef.fasta")
    append_seqs.main(reffile, collapsedfile, withreffile)

    # reverse sequences
    print("== reverse sequences in "+withreffile+" ==")
    reversedfile = withreffile.replace(".fasta", "_rvs.fasta")
    reverse_seq.main(withreffile, reversedfile)

    # align reversed sequences
    print("== Align collapsed file " + reversedfile + " ==")
    reversedalignfile = reversedfile.replace(".fasta", "_align.fasta")
    muscle_align.main(reversedfile, reversedalignfile)

    # reverse back to make a left-aligned alignment
    print("== reverse sequences in " + reversedalignfile + " ==")
    alignmentdir = outdir + "/alignments"
    if os.path.isdir(alignmentdir) is False:
        os.mkdir(alignmentdir)
    withrefalignfile = alignmentdir + "/" + filename.replace(".fasta", "_collapsed_withRef_align.fasta")
    reverse_seq.main(reversedalignfile, withrefalignfile)

    # uncollapse alignment
    print("== Uncollapse alignment " + withrefalignfile + " ==")
    uncollapsealignfile = withrefalignfile.replace("_collapsed", "")
    namefile = collapsedfile.replace(".fasta", ".name")
    uncollapse_seqs.main(withrefalignfile, namefile, uncollapsealignfile)

    # verify sequences
    print("== Verify sequences between " + uncollapsealignfile + " and " + file + " ==")
    verify_seq_origin.main(uncollapsealignfile, file)

    # retrieve gag/pol/env
    gene = regionGene[region]
    genedir = alignmentdir + "/" + gene
    if os.path.isdir(genedir) is False:
        os.mkdir(genedir)
    sgene = geneStart[gene]
    egene = geneEnd[gene]
    print("== Extract "+gene+" at HXB2 position of "+str(sgene)+" to "+str(egene)+ " ==")
    uncollapsealignfilefields = uncollapsealignfile.split("/")
    uncollapsealignfilename = uncollapsealignfilefields[-1]
    genefilename = uncollapsealignfilename.replace("_"+region+"_", "_"+gene+"_NT_")
    genefilename = genefilename.replace("_align.", ".")
    genefile = genedir + "/" + genefilename
    extract_alignment_portion.main(uncollapsealignfile, genefile, region, gene, sgene, egene)

    # translation
    print("== Translate nucleotide to amino acid sequences for "+genefile+" ==")
    geneaafile = genefile.replace("_NT_", "_AA_")
    ntAlignment2aaAlignment.main(genefile, geneaafile)

    # retrieve functional protein sequences and write summary
    print("== Retrieve functional amino acid sequences in "+geneaafile+" ==")
    summaryfile = "functional_summary.csv"
    if filecount == 1:
        with open(summaryfile, "w") as sfp:
            sfp.write("file,gene,sequences,functional,premature_stop,big_deletion,missing_5',not_start_M,median_sequence_length"+"\n")
    retrieve_functional_aa_seqs.main(geneaafile, summaryfile, gene)

    # remove HXB2 and strip all gap columns
    functionalaafile = geneaafile.replace(".fasta", "_functional.fasta")
    gapstripfile = functionalaafile.replace("_withRef_", "_")
    print("== Remove HXB2 and strip all gap columns in "+functionalaafile+" ==")
    strip_all_gaps.main(functionalaafile, gapstripfile)

    # collapse functional protein sequence alignment
    print("== Collapse sequences in "+gapstripfile+" ==")
    collapsedgapstripfile = gapstripfile.replace(".fasta", "_collapsed.fasta")
    collapse_seqs.main(gapstripfile, collapsedgapstripfile)
