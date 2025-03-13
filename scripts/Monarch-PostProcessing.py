#! /usr/bin/python3
import sys
import os
from Bio import SeqIO
from Monarch_Modules import *




# ------------------------------------------------------------------------------
#                          Functions for part 2
# ------------------------------------------------------------------------------

def ReadEnsemblesFull(ensembles):
    EnsDict = dict()
    with open(ensembles,'r') as file:
        for line in file.readlines():
            row = line.strip().split("\t")
            EnsDict[row[4]] = [str(row[0]),str(row[1]),int(row[2]),int(row[3]),
                               str(row[4]),int(row[5]),str(row[6])]
    return EnsDict

def GetReadSeqs(EnsDict, reads, JuncSeqPath):
    with open(reads) as tmp:
        for read in SeqIO.parse(tmp, "fasta"):
            if read.id not in EnsDict.keys():
                pass
            else:
                # 1) Append read seq
                EnsDict[read.id].append(read.seq)
                # 2) Read read tmp file and...
                AllSeqs = (open(os.path.join(JuncSeqPath, read.id.replace("/", "_")), 'r').
                           readline().
                           strip().
                           split('\t'))
                # 3) Append all seqs contained in AllSeqs file
                for seq in AllSeqs:
                    EnsDict[read.id].append(seq)

def PrintResults(EnsDict, ens_full_seqs):
    with open(ens_full_seqs, 'a') as out:
        for l in EnsDict.values():
            for v in l:
                out.write(str(v)+'\t')
            out.write('\n')

# ------------------------------------------------------------------------------

def main():
    # part 1
    # construct ensembles with read names
    input = sys.argv[1]
    output_dir = str(sys.argv[2])
    prefix = str(sys.argv[3])
    # -------------------------------------
    ens_full = os.path.join(output_dir, prefix + ".ensembles.full.readnames.tmp")
    ens_simple = os.path.join(output_dir, prefix + ".ensembles.simple")

    ensembles = GetEnsembles(input, ens_full, transcript_seqs = False,
                             header = False)

    PrintEnsemblesSimple(ensembles, ens_simple)

    # part 2
    # get read sequences for ensembles.full
    # also get remaining cols inside the tmp file named with the ReadName

    reads = sys.argv[4]
    ens_full_seqs = os.path.join(output_dir, prefix + ".ensembles.full.readseqs.tmp")
    EnsDict = ReadEnsemblesFull(ens_full)
    JuncSeqPath= os.path.join(output_dir, "tmp")
    GetReadSeqs(EnsDict, reads, JuncSeqPath)
    PrintResults(EnsDict, ens_full_seqs)

if __name__ == '__main__':
    main()
