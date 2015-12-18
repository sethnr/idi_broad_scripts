#!/usr/bin/python

import vcf
import sys
import argparse
from string import *

parser = argparse.ArgumentParser(description='from CHR : POS : 1/0 file, add true/false to VCF')

parser.add_argument('-v','-vcf', action="store", dest='vcfFile', type=str, help='vcf file', nargs='?', default=None)
parser.add_argument('-f','-file', action="store", dest='valueFile', type=str, help='valueFile', nargs='?', default=None)

parser.add_argument('-i','--infoname', action="store", dest='valname', type=str, help='make combined outfile named <outfile> instead of separate files', nargs='?', default="VALUE")
#parser.add_argument('-i','--defaultSize', action="store", dest='defSize', type=int, help='default size for STRs', nargs='?', default=0)


args = parser.parse_args()

vcfFile1 = open(args.vcfFile,'r')
vcfReader=vcf.Reader(vcfFile1)

vcfoutF1 = replace(args.vcfFile,'.vcf','.'+args.valname+'.vcf')
vcfoutF1 = replace(vcfoutF1,'.vcf.gz','.'+args.valname+'.vcf')
#print >>sys.stdout, "VCFOUT1",vcfoutF1
vcfoutF1 = open(vcfoutF1,'w')
vcfout1=vcf.Writer(vcfoutF1,vcfReader)

valfile = open(args.valueFile,'r')

vals = dict()
for vline in valfile:
    c,p,v = vline.split()
    #if p == 451480: print c,p,v,bool(v) 
    vals[(c,int(p))] = bool(int(v))

for rec in vcfReader:
    #if rec.is_indel:
    if (rec.CHROM,rec.POS) in vals:
        rec.INFO[args.valname]=vals[(rec.CHROM,rec.POS)]
    else:
        print >>sys.stderr, "warning: ",rec.CHROM+":"+str(rec.POS)," not found in valfile"
    vcfout1.write_record(rec)

