#!/bin/python

import sys
import argparse
from string import *
from os import path

parser = argparse.ArgumentParser(description='parse out statistics from VCF ')

parser.add_argument('-f1','--file1', action="store", dest='file1', type=str, help='freec file1', nargs='?', default=None)
parser.add_argument('-f2','--file2', action="store", dest='file2', type=str, help='freec file2', nargs='?', default=None)
parser.add_argument('-mean', action="store_true", dest='mean', help='print mean distance & copy number',)
parser.add_argument('-min', action="store_true", dest='min', help='print min overlap & mean copy number',)
parser.add_argument('-max', action="store_true", dest='max', help='print combined distance & mean copy number',)

args = parser.parse_args()
file1 = open(args.file1,'r')
file2 = open(args.file2,'r')

def readCNV(cnvfile):
    line = cnvfile.readline()
#    print >>sys.stderr, cnvfile,"LINE:"+line+":",line is None,len(line),len(line)==0
    #if len(line)==0: line = None
#    if line is not None:
    if len(line)==0:
        (chr,st,en,cn,cnv) = (None,None,None,None,None)
    else:
        (chr,st,en,cn,cnv) = line.split()
        st = int(st); en=int(en); cn=int(cn)
    return (chr,st,en,cn,cnv)

#line1 = file1.readline()
#line2 = file2.readline()
#(chr1,st1,en1,cn1,cnv1) = line1.split()
#(chr2,st2,en2,cn2,cnv2) = line2.split()
(chr2,st2,en2,cn2,cnv2) = readCNV(file2)
(chr1,st1,en1,cn1,cnv1) = readCNV(file1)

while (chr1 is not None and chr2 is not None):
#    print "\t".join(map(str,[chr1,st1,en1,chr2,st2,en2,cn1,cn2,cnv1,cnv2])),

    if chr1 == chr2 and (en1 > st2 and st1 < en2):
        pcOl = round((int(en1)-int(st1))/(float(en2)-int(st2)),2)        
#        print pcOl,
        minst = min(st1,st2)
        maxst = max(st1,st2)
        minen = min(en1,en2)
        maxen = max(en1,en2)
        meanst=(st1+st2)/2
        meanen=(en1+en2)/2
        meancn=(cn1+cn2)/2
        if args.mean:
            print "\t".join(map(str,[chr1,meanst,meanen,meancn,cnv1,pcOl]))
        elif args.max:
            print "\t".join(map(str,[chr1,minst,maxen,meancn,cnv1,pcOl]))
        elif args.min: 
            print "\t".join(map(str,[chr1,maxst,minen,meancn,cnv1,pcOl]))
        else:
            print "\t".join(map(str,[chr1,maxst,minen,minst,maxen,cn1,cn2,cnv1,pcOl]))
    else:
#        print ".",
        pass
    
    if chr1 > chr2:
#        print "1>2",
        (chr2,st2,en2,cn2,cnv2) = readCNV(file2)
    elif chr2 > chr1:
#        print "2>1",
        (chr1,st1,en1,cn1,cnv1) = readCNV(file1)
    elif chr1 == chr2 and (st1 > st2 or (st1==st2 and en1>en2)):
#        print "1>2",
        (chr2,st2,en2,cn2,cnv2) = readCNV(file2)
    elif chr2 == chr1 and (st2 > st1 or (st1==st2 and en2>en1)):
#        print "2>1",
        (chr1,st1,en1,cn1,cnv1) = readCNV(file1)
    elif chr1 == chr2 and st1 == st2 and en1==en2:
#        print "2=1",
        (chr1,st1,en1,cn1,cnv1) = readCNV(file1)
        (chr2,st2,en2,cn2,cnv2) = readCNV(file2)
#    print ""
