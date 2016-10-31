#!/usr/bin/env python
'''
module for cpg_gene script
'''
import Bio.SeqUtils
import gffutils
import argparse

# Retrieve args from command line, usage help
def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("FastaFile", help="Fasta File.", type=str)
    parser.add_argument("GffFile", help="Matching GFF File.", type=str)
    parser.add_argument("FileOut", help="OutFile.png Name.", type=str)

    # Optional arguments
    parser.add_argument("-s", help="Step Size.", type=int, default=50)
    parser.add_argument("-w", help="Window Size.", type=int, default=200)
    parser.add_argument("-oe", help="Minimum Observed Expected CpG.", 
                        type=float, default=1.)
    parser.add_argument("-gc", help="Minimum GC.", type=float, default=.5)
    parser.add_argument("-r", help="Range from ATG.", type=int, default=5000)
    parser.add_argument("-f", help="GFF Feature", type=str, default="gene")
    parser.add_argument("-i", help="GFF Attribute Gene ID", type=str, 
                        default="")

    # Parse arguments
    args = parser.parse_args()

    return args

# Insure WinSize, SeqRecord, and Step are correct lengths
def arg_seqcheck(SeqRecord, WinSize, Step):
    if not ((type(WinSize) == type(0)) and (type(Step) == type(0))):
        raise Exception("**NOTE type(WinSize) and type(Step) must be int.")
    if Step > WinSize:
        raise Exception("**NOTE Step must not be larger than WinSize.")
    if WinSize > len(SeqRecord.seq):
        raise Exception('''**NOTE WinSize must not 
                        be larger than Sequence length.''')

# Class contains GC%, expected CpG, and Observed CpG
class SeqRec:
    def __init__(self, SeqId, Seq, WinSize, WinStart):
        self.SeqId = SeqId
        self.Seq = Seq.upper()
        self.WinSize = WinSize
        self.WinStart = WinStart

    def observ(self):
        CG = self.Seq.count('CG')
        return CG

    def expect(self):
        G = self.Seq.count('G')
        C = self.Seq.count('C')
        ExpGC = (G * C) / float(self.WinSize)
        return ExpGC

    def gc_cont(self):
        G = self.Seq.count('G')
        perGC = Bio.SeqUtils.GC(self.Seq) / float(100)
        return perGC

# Pre-Determine the number of CpG (potential) islands
def chunks(SeqRecord, WinSize, Step):
    NumOfChunks = ((len(SeqRecord)-WinSize)/Step)+1
    return NumOfChunks

# Push chunked seq data through SeqRec class
# Return array of SeqRec instances
def compute(SeqRecord, Step, NumOfChunks, WinSize):
    SeqRecList = []
    GC = 0
    for i in range(0, NumOfChunks*Step, Step):
        WinStart = i
        # Creating SeqRec instance, appending to array
        SeqRecList.append(SeqRec(SeqRecord.id, SeqRecord.seq[i:i + WinSize], 
                                    WinSize, WinStart))
    return SeqRecList

# Compare SeqRec instances to desired feature in GFF DB
# Returns distance to closest feature within provided ranges def= 5000bp
def get_closest(Rec, GffDb, StartRange, FeatGFF, ID_Feat):
    DistList = []
    Region = gffutils.FeatureDB.region(GffDb, 
                                        seqid=Rec.SeqId, featuretype=FeatGFF)    
    for R in Region:
        #print R.id
        DistVal = R.start - Rec.WinStart
        Cond1 = abs(DistVal) < StartRange
        if len(ID_Feat) > 0:
            Cond2 = R.id == ID_Feat
        else:
            Cond2 = True
        if Cond1 and Cond2:
            # Array of abs dist value, true value, and GFF feature name (ID=)
            DistList.append([abs(DistVal), DistVal, R.id])
    SortedList = sorted(DistList, key=lambda l:l[0])
    if len(SortedList) == 0:
        return False
    elif len(SortedList) > 0:
        return SortedList[0]
