#!/bin/python

import sys
import argparse
from string import *
from os import path

parser = argparse.ArgumentParser(description='parse out statistics from VCF ')

parser.add_argument('-f','--file', action="append", dest='files', type=str, help='freec file1', nargs='?')
#parser.add_argument('-mean', action="store_true", dest='mean', help='print mean distance & copy number',)
#parser.add_argument('-min', action="store_true", dest='min', help='print min overlap & mean copy number',)
#parser.add_argument('-max', action="store_true", dest='max', help='print combined distance & mean copy number',)

args = parser.parse_args()

c=dict()
c[None] = ""
c["gain"]="+"
c["loss"]="-"

def readCNV(cnvfile):
    line = cnvfile.readline()
#    print >>sys.stderr, cnvfile,"LINE:"+line+":",line is None,len(line),len(line)==0
    #if len(line)==0: line = None
#    if line is not None:
    if len(line)==0:
        (chrom,st,en,cn,cnv) = (None,None,None,None,None)
    else:
        F = line.split()
        (chrom,st,en,cn,cnv) = F[0:5]
        st = int(st); en=int(en); cn=float(cn)
    cnv = c[cnv]
    return (chrom,st,en,cn,cnv)

CNVs = dict()
fileI=-1
#print args.files
for cnvfile in args.files:
    fileI+=1
    sampCNV = dict()
    cnvs = open(cnvfile,'r')
    (chrom,st,en,cn,cnv) = readCNV(cnvs)
#    print "\t".join(map(str,[chrom,st,en,cn,cnv]))
    while chrom is not None:
#        print "\t".join([chrom,str(st)])
        for i in range(st,en):
#            print chrom,i
            if (chrom,i) not in CNVs:
                CNVs[(chrom,i)] = [""]*len(args.files)
            if CNVs[(chrom,i)][fileI] != "" and CNVs[(chrom,i)][fileI] != cnv:
                print >>sys.stderr, chrom+" "+str(i)+" "+cnvfile+" found ("+CNVs[(chrom,i)][fileI]+"), overwriting "+cnv
#            sampCNV[(chr,i)] = cnv
                CNVs[(chrom,i)][fileI] = "!"
            else:
                CNVs[(chrom,i)][fileI] = cnv
        #    sampCNV
        (chrom,st,en,cn,cnv) = readCNV(cnvs)
        
print "\t".join(["CHR","ST","EN"]+args.files)

lastChr = None
lastI = -1
lastCalls = []
startChr = None
startI = -1
startCalls = []
for (chrom,i) in sorted(CNVs):
    calls = CNVs[(chrom,i)]
    if chrom != lastChr or i != lastI+1 or calls != lastCalls:
        if lastChr is not None:
            #print "\t".join([lastChr,str(lastI)]+map(str,lastCalls))
            print "\t".join(map(str,[lastChr,startI,lastI]+lastCalls))
        startI=i; startChr=chrom; startCalls=calls
    lastI=i; lastChr=chrom; lastCalls=calls
print "\t".join(map(str,[lastChr,startI,lastI]+lastCalls))
#print "\t".join([chrom,str(i)]+map(str,CNVs[(chrom,i)]))
