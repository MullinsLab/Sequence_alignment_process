#!/usr/bin/python

########################################################
# In a directory with GP alignments and a directory with POL alignments,
# extract the overlap region of GP and POL from GP and POL alignments of the same subject id,
# combine overlapped alignments,
# align combined sequences,
# run fasttree on the alignment
# input: GP directory path and POL directory path
# output: combined alignment fasta file and tree file
# Wenjie Deng
# 2022-02-04
########################################################
import os
import sys
import re
import argparse
from collections import defaultdict
import extract_alignment_portion
import append_seqs
import muscle_align

def nested_dict(n, type):
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-gd", "--gpdir", help="directory to hold GP alignment fasta files", nargs="?", const=1, type=str,
                        default=".")
    parser.add_argument("-pd", "--poldir", help="directory to hold input POL alignment fasta files", nargs="?", const=1, type=str,
                        default="POL")
    parser.add_argument("-od", "--outdir", help="directory to hold input sequence fasta file", nargs="?", const=1,
                        type=str,
                        default="Output")
    args = parser.parse_args()
    gpdir = args.gpdir
    poldir = args.poldir
    outdir = args.outdir

    id = ""
    lineageStatus, idStatus = {}, {}
    lineages, ids = [], []
    idLineageCount = nested_dict(2, int)

    for filename in os.listdir(gpdir):
        gpfile = os.path.join(gpdir, filename)
        if os.path.isfile(gpfile):
            if re.search("\.fasta$", gpfile):
                print("=== Processing "+gpfile+" ===")
                polfile = ""
                filematch = re.search('^(.*?)_(.*?)_(.*).fasta$', filename)
                if filematch:
                    fields = filename.split("_")
                    id = fields[0] + "_" + fields[1] + "_"
                    for fname in os.listdir(poldir):
                        file = os.path.join(poldir, fname)
                        if os.path.isfile(file):
                            if re.search(id, file):
                                polfile = file
                                break
                    if polfile:
                        if os.path.isdir(outdir) is False:
                            os.mkdir(outdir)
                        gpovlpfilename = id + "GP_overlap.fasta"
                        polovlpfilename = id + "POL_overlap.fasta"
                        gpovlpfile = os.path.join(outdir, gpovlpfilename)
                        polovlpfile = os.path.join(outdir, polovlpfilename)
                        print("== Extract GP-POL overlap region from "+gpfile+"==")
                        extract_alignment_portion.main(gpfile, gpovlpfile, "GP", "GP", 1135, 0, True)
                        print("== Extract GP-POL overlap region from " + polfile + "==")
                        extract_alignment_portion.main(polfile, polovlpfile, "POL", "POL", 1, 1653, True)
                        print("== Combine "+gpovlpfile+" and "+polovlpfile+" files ==")
                        gppolovlpfilename = id +"GP-POL_overlap.fasta"
                        gppolovlpfile = os.path.join(outdir, gppolovlpfilename)
                        append_seqs.main(gpovlpfile, polovlpfile, gppolovlpfile)
                        print("== Align "+gppolovlpfile+" via Muscle ==")
                        gppolovlpalignfile = gppolovlpfile.replace(".fasta", "_align.fasta")
                        muscle_align.main(gppolovlpfile, gppolovlpalignfile)
                        print("== Run FastTree on "+gppolovlpalignfile+" ==")
                        treefile = gppolovlpalignfile.replace(".fasta", ".tre")
                        os.system("FastTree -nt -gtr < "+gppolovlpalignfile+" > "+treefile)
                    else:
                        print("No corresponding POL alignments in "+poldir+"/ for "+gpfile+"\n")
                        continue
                else:
                    sys.exit("file name not formatted: " + gpfile)

    print("== All done! ==\n")

