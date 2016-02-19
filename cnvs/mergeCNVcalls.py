#!/bin/python

import sys
import argparse
from string import *
from os import path
import re 
parser = argparse.ArgumentParser(description='parse out statistics from VCF ')

parser.add_argument('-f','--freec', action="store", dest='freec', type=str, help='freec file1', nargs='?', default=None)
parser.add_argument('-pD','--pindelD', action="store", dest='pindelD', type=str, help='pindel deletions file', nargs='?', default=None)
parser.add_argument('-pTD','--pindelTD', action="store", dest='pindelTD', type=str, help='pindel tandem duplications file', nargs='?', default=None)
parser.add_argument('-c','--cnvnator', action="store", dest='cnvnator', type=str, help='pindel tandem duplications file', nargs='?', default=None)
parser.add_argument('--mean', action="store_true", dest='mean', help='print mean distance & copy number',)
parser.add_argument('--min', action="store_true", dest='min', help='print min overlap & mean copy number',)
parser.add_argument('--max', action="store_true", dest='max', help='print combined distance & mean copy number',)

args = parser.parse_args()

chrP,stP,enP,cnP,cnvP = None,None,None,None,None
chrC,stC,enC,cnC,cnvC = None,None,None,None,None
chrD,stD,enD,cnD,cnvD = None,None,None,None,None
chrTD,stTD,enTD,cnTD,cnvTD = None,None,None,None,None

def readFREEC(cnvfile):
    line = cnvfile.readline()
#    print >>sys.stderr, cnvfile,"LINE:"+line+":",line is None,len(line),len(line)==0
    #if len(line)==0: line = None
#    if line is not None:
    if len(line)==0:
        (chr,st,en,cn,cnv) = (None,None,None,None,None)
    else:
        F = line.split()
        (chr,st,en,cn,cnv) = F[0:5]
        st = int(st); en=int(en); cn=int(cn)
    return (chr,st,en,cn,cnv)

def readCNVNATOR(cnvfile):
    line = cnvfile.readline()
#    print >>sys.stderr, cnvfile,"LINE:"+line+":",line is None,len(line),len(line)==0
    #if len(line)==0: line = None
#    if line is not None:
    if len(line)==0:
        (chr,st,en,cn,cnv) = (None,None,None,None,None)
    else:
        F = line.split()
        cn=float(F[3]) #copy numbers normalized to 1 in cnvnator
        if F[0] == "deletion": cnv="loss"
        elif F[0] == "duplication": cnv="gain"
        chr,st,en = re.split("[\:\-]",F[1])
        st = int(st); en=int(en); #cn=int(cn)
    return (chr,st,en,cn,cnv)

def readPINDELS(fileD,fileTD):
    global chrP,stP,enP,cnP,cnvP
    global chrD,stD,enD,cnD,cnvD
    global chrTD,stTD,enTD,cnTD,cnvTD
    #if an init state, get both
    if (chrTD,stTD,enTD,cnTD,cnvTD) == (None,None,None,None,None) and (chrD,stD,enD,cnD,cnvD) == (None,None,None,None,None):
        (chrTD,stTD,enTD,cnTD,cnvTD) = readPINDEL(fileTD)
        (chrD,stD,enD,cnD,cnvD) = readPINDEL(fileD)
    #if last compared was TD, get new TD
    if ((chrP,stP,enP,cnP,cnvP) == (chrTD,stTD,enTD,cnTD,cnvTD)):
        (chrTD,stTD,enTD,cnTD,cnvTD) = readPINDEL(fileTD)
    #if last compared was D, get new D
    elif ((chrP,stP,enP,cnP,cnvP) == (chrD,stD,enD,cnD,cnvD)):
        (chrD,stD,enD,cnD,cnvD) = readPINDEL(fileD)

    #if one file has ended, return the other, else return the min
    if chrD is None:
        (chrP,stP,enP,cnP,cnvP) = (chrTD,stTD,enTD,cnTD,cnvTD)
    elif chrTD is None:
        (chrP,stP,enP,cnP,cnvP) = (chrD,stD,enD,cnD,cnvD)
    else:
        (chrP,stP,enP,cnP,cnvP) = min((chrD,stD,enD,cnD,cnvD),(chrTD,stTD,enTD,cnTD,cnvTD))
    return (chrP,stP,enP,cnP,cnvP)
    
def readPINDEL(pinfile):
#    print >>sys.stderr, cnvfile,"LINE:"+line+":",line is None,len(line),len(line)==0
    #if len(line)==0: line = None
#    if line is not None:
    line="INIT"; cn=1;cnv="."
    header=False

    line = pinfile.readline()
    while line:
        isheader=re.match('#', line)
#        print line, len(line),isheader
        if isheader:
#            print "header found"
            header=True
            line = pinfile.readline()
#            print line
            F = line.split()
#            print "\t".join(F[0:15])
            if F[1]=="D": cnv="loss"; cn=0
            elif F[1]=="TD": cnv="gain"; cn=2
#            (chr,st,en) = (F[7],F[9],F[10])
            (chr,st,en) = (F[7],F[12],F[13])   #use largest pindel range
            st = int(st); en=int(en); 
#            print "\t".join(map(str,[chr,st,en,(en-st),cnv,cn]))
            return (chr,st,en,cn,cnv)
        else:
            pass
#            print "no header"
        line = pinfile.readline()
    return (None,None,None,None,None)


#line1 = file1.readline()
#line2 = file2.readline()
#(chr1,st1,en1,cn1,cnv1) = line1.split()
#(chr2,st2,en2,cn2,cnv2) = line2.split()
#get tandem dupications and deletions, put earlies into pindel file 
#(chrD,stD,enD,cnD,cnvD) = readPINDEL(pindelD)
#(chrTD,stTD,enTD,cnTD,cnvTD) = readPINDEL(pindelTD)
#(chrP,stP,enP,cnP,cnvP) = min((chrD,stD,enD,cnD,cnvD),(chrTD,stTD,enTD,cnTD,cnvTD))

if args.pindelD and args.freec and args.cnvnator:
    print >>sys.stderr, "only making pairwise comparisons for the moment"
    exit(1)
elif args.pindelD and args.freec:
    freec = open(args.freec,'r')
    pindelD = open(args.pindelD,'r')
    pindelTD = open(args.pindelTD,'r')
    (chrP,stP,enP,cnP,cnvP) = readPINDELS(pindelD,pindelTD)
    (chrF,stF,enF,cnF,cnvF) = readFREEC(freec)
    readCNV1 = lambda : readPINDELS(pindelD,pindelTD)
    readCNV2 = lambda : readFREEC(freec)

elif args.pindelD and args.cnvnator:
    pindelD = open(args.pindelD,'r')
    pindelTD = open(args.pindelTD,'r')
    cnvnator = open(args.cnvnator,'r')
    (chrP,stP,enP,cnP,cnvP) = readPINDELS(pindelD,pindelTD)
    (chrC,stC,enC,cnC,cnvC) = readCNVNATOR(cnvnator)
    readCNV1 = lambda : readPINDELS(pindelD,pindelTD)
    readCNV2 = lambda : readCNVNATOR(cnvnator)

elif args.freec and args.cnvnator:
    freec = open(args.freec,'r')
    cnvnator = open(args.cnvnator,'r')
    (chrF,stF,enF,cnF,cnvF) = readFREEC(freec)
    (chrC,stC,enC,cnC,cnvC) = readCNVNATOR(cnvnator)
    readCNV1 = lambda : readFREEC(freec)
    readCNV2 = lambda : readCNVNATOR(cnvnator)

last1=(None,None,None,None,None)
last2=(None,None,None,None,None)

(chr1,st1,en1,cn1,cnv1) = readCNV1()
(chr2,st2,en2,cn2,cnv2) = readCNV2()

while (chr1 is not None and chr2 is not None):

    if 1==2:
        if (chr2,st2,en2,cn2,cnv2) == last2:
            print "\t".join(map(str,[chr1,st1,en1,"   .   ","   .   ","   .   ",cn1,cn2,cnv1,cnv2]))
        elif (chr1,st1,en1,cn1,cnv1) == last1:
            print "\t".join(map(str,["   .   ","   .   ","   .   ",chr2,st2,en2,cn1,cn2,cnv1,cnv2]))
        else:
            print "\t".join(map(str,[chr1,st1,en1,chr2,st2,en2,cn1,cn2,cnv1,cnv2]))

    if chr1 == chr2 and (en1 > st2 and st1 < en2) and (cnv1==cnv2):
        pcOl = round((int(en1)-int(st1))/(float(en2)-int(st2)),2)        
#        print pcOl,
        minst = min(st1,st2)
        maxst = max(st1,st2)
        minen = min(en1,en2)
        maxen = max(en1,en2)
        meanst=(st1+st2)/2
        meanen=(en1+en2)/2
        meancn=(cn1+cn2)/2.0
        
        if args.mean:
            print "\t".join(map(str,[chr1,meanst,meanen,meancn,cnv1,pcOl]))
        elif args.max:
            print "\t".join(map(str,[chr1,minst,maxen,meancn,cnv1,pcOl]))
        elif args.min: 
            print "\t".join(map(str,[chr1,maxst,minen,meancn,cnv1,pcOl]))
        else:
            print "\t".join(map(str,[chr1,maxst,minen,minst,maxen,cn1,cn2,cnv1,cnv2,pcOl]))
#            print "\t".join([chr1,str(st1)+"-"+str(en1),str(en1-st1),str(st2)+"-"+str(en2),str(en2-st2),str(pcOl)])
    else:
#        print ".",
        pass

    last2=(chr2,st2,en2,cn2,cnv2)
    last1=(chr1,st1,en1,cn1,cnv1)
    
    if chr1 > chr2:
        (chr2,st2,en2,cn2,cnv2) = readCNV2()
    elif chr2 > chr1:
        (chr1,st1,en1,cn1,cnv1) = readCNV1()
    elif chr1 == chr2 and (st1 > st2 or (st1==st2 and en1>en2)):
        (chr2,st2,en2,cn2,cnv2) = readCNV2()
    elif chr2 == chr1 and (st2 > st1 or (st1==st2 and en2>en1)):
        (chr1,st1,en1,cn1,cnv1) = readCNV1()
    elif chr1 == chr2 and st1 == st2 and en1==en2:
        (chr1,st1,en1,cn1,cnv1) = readCNV1()
        (chr2,st2,en2,cn2,cnv2) = readCNV2()
