#!/usr/bin/env python

'''
Purpose:
This script, using default values, determines and plots the CpG islands in 
relation to a given feature "type" (e.g. "gene" or "mRNA") from a GFF file 
which corresponds to the user-provided fasta file.

Note:
CpG Islands are determined by ObEx = (Observed CpG) / (Expected CpG) ,
default threshold > 1. 

Where Expected CpG = (count(C) * count(G)) / WindowSize

Usage:
python cpg_gene.py FastaFile Gff_File OutFile.png

Default optional parameters:
    -s, Step Size, default = 50
    -w, Window Size, default = 200
    -oe, Minimum Observed Expected CpG, default = 1
    -gc, Minimum GC, default = .5
    -r Range from ATG, or provided feature, default = 5000
    -f, GFF Feature, default = "gene"
    -i, Gene ID from GFF, default = ""
'''

import sys
import os
import argparse

from collections import Counter
from Bio import SeqIO
import cpgmod
import gffutils

import pandas as pd
import numpy as np
from ggplot import *


# Capture command line args, with or without defaults
if __name__ == '__main__':
    # Parse the arguments
    LineArgs = cpgmod.parseArguments()

# Populate vars with args
FastaFile = LineArgs.FastaFile
GffFile = LineArgs.GffFile
OutFile = LineArgs.FileOut
Step = LineArgs.s
WinSize = LineArgs.w
ObExthresh = LineArgs.oe
GCthresh = LineArgs.gc
StartRange = LineArgs.r
FeatGFF = LineArgs.f
ID_Feat = LineArgs.i

# Gather all possible CpG islands
MergedRecs = []
print "Parsing sequences...\n"
for SeqRecord in SeqIO.parse(FastaFile, "fasta"):
    print SeqRecord.id
    # Determine if sequences and args are acceptable
    cpgmod.arg_seqcheck(SeqRecord, WinSize, Step)
    # Pre-determine number of islands
    NumOfChunks = cpgmod.chunks(SeqRecord, WinSize, Step)
    # Return array of SeqRec class (potential CpG island) instances
    SeqRecList = cpgmod.compute(SeqRecord, Step, NumOfChunks, WinSize)
    MergedRecs = MergedRecs + SeqRecList

# Create GFF DB
GffDb = gffutils.create_db(GffFile, dbfn='GFF.db', force=True, keep_order=True, 
                            merge_strategy='merge', sort_attribute_values=True,
                            disable_infer_transcripts=True,
                            disable_infer_genes=True)

print "\nGFF Database Created...\n"

# Filter out SeqRec below threshold 
DistArr = []
for Rec in MergedRecs:
    Cond1 = Rec.expect() > 0     
    if Cond1 == True:
        ObEx = (Rec.observ() / Rec.expect())
        Cond2 = ObEx > ObExthresh 
        Cond3 = Rec.gc_cont() > GCthresh
        if Cond2 and Cond3:
            # Query GFF DB for closest gene feature *or provided feature*
            Arr = cpgmod.get_closest(Rec, GffDb, StartRange, FeatGFF, ID_Feat)
            if Arr <> False:
                Arr.append(ObEx)
                DistArr.append(Arr)

print "CpG Islands predicted...\n"
print "Generating Figure...\n"

# Releasing SeqRecs
MergedRecs = None
SeqRecList = None

# Pre-check DistArr Results
if len(DistArr) < 2:
    print "WARNING, "+ str(len(DistArr)) + " sites were found."
    print "Consider changing parameters.\n"

# Generate Figure:
ObExRes = pd.DataFrame({
                    'gene' : [],
                    'xval': [],
                    'yval': []})

try:
    Cnt = 0
    for Dist in DistArr:
        Cnt += 1
        print "PROGRESS: "+str(Cnt) +" of "+ str(len(DistArr))
        ObExdf = pd.DataFrame({
                        'gene': [Dist[2]],
                        'xval': [Dist[1]],
                        'yval': [Dist[3]]})
        ObExFram = [ObExRes, ObExdf]
        ObExRes = pd.concat(ObExFram, ignore_index=True)
    p = ggplot(aes(x='xval', y='yval'), data=ObExRes) \
        + geom_point() \
        + ylab("Observed/Expected CpG") \
        + xlab("Position (bp) Relative to (ATG = 0)") \
        + ggtitle("Predicted CpG Island Position Relative to ATG")
    p.save(OutFile)
except IndexError as e:
    print 'Error: '+ str(e)
    sys.exit('Exiting script...')
print p

# Remove GFF DB
os.remove('GFF.db')
