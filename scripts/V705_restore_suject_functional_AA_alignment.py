#!/usr/bin/python3

##########################################################################################
# Program: V705_restore_subject_functional_AA_alignment.py
# Purpose: From a cross-subjects functional consensus AA alignment:
# restore each subject's functional AA alignment based on subject consensus and collapsed sequence alignment,
# uncollapse to subject's uncollapsed functional AA alignment based on name file
# Author: Wenjie Deng
# Date: 2022-11-21
##########################################################################################

import sys
import re
import os
import argparse
import verify_seq_origin
import uncollapse_seqs

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("consensusalignfile", help="cross-subjects consensus sequence alignment fasta file")
    parser.add_argument("-ad", "--alignmentdir",
                        help="directory to hold subject consensus and collapsed sequences alignment fasta files",
                        nargs="?", const=1, type=str)
    parser.add_argument("-nd", "--namefiledir",
                        help="directory to hold subject collapsed sequences name files",
                        nargs="?", const=1, type=str)
    parser.add_argument("-sd", "--seqdir",
                        help="directory to hold subject uncollapsed sequence fasta files",
                        nargs="?", const=1, type=str)
    args = parser.parse_args()
    consalignfile = args.consensusalignfile
    aligndir = args.alignmentdir
    namedir = args.namefiledir
    seqdir = args.seqdir

    cwd = os.getcwd()
    collapsedoutdir = cwd + "/restore_subject_functional_collapsed_AA_alignment_outputs"
    uncollapsedoutdir = cwd + "/restore_subject_functional_uncollapsed_AA_alignment_outputs"
    if os.path.isdir(collapsedoutdir) is False:
        os.mkdir(collapsedoutdir)
    if os.path.isdir(uncollapsedoutdir) is False:
        os.mkdir(uncollapsedoutdir)

    consnames, collapsenames = [], []
    consnameSeq = {}
    name, refname, refseq = '', '', ''
    with open(consalignfile, "r") as cfp:
        for line in cfp:
            line = line.strip()
            if line:
                linematch = re.search(r"^>(\S+)", line)
                if linematch:
                    name = linematch.group(1)
                    if re.search("HXB2", name):
                        refname = name
                    else:
                        consnames.append(name)
                    consnameSeq[name] = ''
                else:
                    consnameSeq[name] += line
    refseq = consnameSeq[refname].replace("Z", "*")

    for consname in consnames:
        print("=== processing "+consname+" ===")
        sid, alignfile, namefile, seqfile = '', '', '', ''
        consseq = consnameSeq[consname].replace("Z", "*")
        namematch = re.search(r"^(.*?)_functional_consensus", consname)
        if namematch:
            sid = namematch.group(1)
        else:
            sys.exit("name not formatted: "+consname)

        for filename in os.listdir(aligndir):
            file = os.path.join(aligndir, filename)
            if os.path.isfile(file):
                if re.search(sid, filename):
                    alignfile = file
                    break
        if alignfile == '':
            sys.exit("No alignment file for sid " + sid)

        print("== processing "+alignfile+" ==")
        collapsenames = []
        nameAAs = {}
        collapsename = ''
        with open(alignfile, "r") as afp:
            for line in afp:
                line = line.strip()
                if line:
                    lmatch = re.search(r"^>(\S+)", line)
                    if lmatch:
                        collapsename = lmatch.group(1)
                        if re.search("_functional_consensus", collapsename) is None:
                            collapsenames.append(collapsename)
                            nameAAs[collapsename] = []
                    else:
                        if re.search("_functional_consensus", collapsename):
                            consseqnogap = consseq.replace("-", "")
                            seqnogap = line.replace("-", "")
                            if consseqnogap != seqnogap:
                                print("* consensus sequences do not match: "+consseqnogap+" vs. "+seqnogap+" *")
                                if re.search(consseqnogap, seqnogap):
                                    print("reviewed consensus was trimmed compared to the origianl consensus")
                                else:
                                    sys.exit("reviewed and original consensus sequences are totally different")
                        else:
                            nameAAs[collapsename] = list(line)

        outcollapsefile = os.path.join(collapsedoutdir, sid+"_functional_collapsed.fasta")
        print("= restoring collapsed functional alignment for "+sid+" =")
        with open(outcollapsefile, "w") as ofp:
            ofp.write(">"+refname+"\n"+refseq+"\n")
            collapseseqcount = 0
            consaas = list(consseq)
            for collapsename in collapsenames:
                collapseseqcount += 1
                idx = -1
                aas = []
                for i in range(len(consaas)):
                    if consaas[i] == '-':
                        aas.append('-')
                    else:
                        idx += 1
                        aas.append(nameAAs[collapsename][idx])
                ofp.write(">"+collapsename+"\n")
                ofp.write(''.join(aas) + "\n")
            print("restored total "+str(collapseseqcount)+" collapsed sequences")

        for filename in os.listdir(namedir):
            file = os.path.join(namedir, filename)
            if os.path.isfile(file):
                if re.search(sid, filename):
                    namefile = file
                    break
        if namefile == '':
            sys.exit("No name file for sid " + sid)

        print("== processing " + namefile + " ==")
        print("= uncollapsing functional sequences for "+sid+" =")
        outuncollapsefile = os.path.join(uncollapsedoutdir, sid + "_functional_uncollapsed.fasta")
        uncollapse_seqs.main(outcollapsefile, namefile, outuncollapsefile)

        for filename in os.listdir(seqdir):
            file = os.path.join(seqdir, filename)
            if os.path.isfile(file):
                if re.search(sid, filename):
                    seqfile = file
                    break
        if seqfile == '':
            sys.exit("No uncollapsed sequence file for sid " + sid)

        print("== processing " + seqfile + " ==")
        print("= verifying uncollapsed sequences for " + sid + " =")
        verify_seq_origin.main(outuncollapsefile, seqfile)
        print()
