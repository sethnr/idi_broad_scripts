#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
import collections
import copy

parser = argparse.ArgumentParser(description='get allele numbers table')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile1', type=str, help='vcfFile1', nargs='?', default=None)
parser.add_argument('-i','--info', action="store", dest='info', type=str, help='vcfFile1', nargs='?', default=None)

args = parser.parse_args()

vcfFile1 = open(args.vcfFile1,'r')
reader1=vcf.Reader(vcfFile1)

if args.info is None:
    print "\t".join(["chr","pos","type","subtype","alleles"]+reader1.samples)
else:
    print "\t".join(["chr","pos","type","subtype","alleles"]+reader1.samples+[args.info])

    
for rec in reader1:
    if rec.num_unknown > 2: continue

    alleles = "/".join([rec.REF] + map(str,rec.ALT))
    vartype="COMPLEX";
    if rec.is_indel: vartype="INDEL"
    if rec.is_snp: vartype="SNP"
    print "\t".join(map(str,[rec.CHROM,rec.POS,vartype,rec.var_subtype, alleles])),
    GTindex=rec.FORMAT.split(":").index("GT")
#    print GTindex
    for call in rec.samples:
        callData=call.data[GTindex]
        if callData is not None:
            print "\t"+callData,
        else:
            print "\t.",
    if args.info is not None:
        print "\t",
        if args.info in rec.INFO:
            print rec.INFO[args.info][0],
    print ""
