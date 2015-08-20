#!/usr/bin/python

import sys
import argparse
from string import *
from subprocess import call, check_output
import os.path as path
from math import floor,ceil


import numpy as np
import vcf
#from pysam import TabixFile
import tabix

if len(sys.argv)==1:
    print >>sys.stdout, sys.argv[0] + """::
  This script will split a vcf into blocks of Nkb and 
  summarize variables from the INFO field of a VCF - 
  selecting VAR:Count will show a count of the no of variants 
  and the proportion of the total and VAR:Mean will show the count 
  and the mean of the value for that INFO field.
    eg: in 10kb blocks across 300kb region, get mean of all mafs and 
  count all vars with annotation field:
  
    python blockSummariseVCF.py \\
      -v myVCFfile.vcf -b 10  -r Pf3d7_07_v3:450000-750000 \\
      -i MAF:Mean -i ANO:Count
"""
    exit(0)

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile', type=str, help='vcf file', nargs='?', default=None)
parser.add_argument('-o','--outname', action="store", dest='valname', type=str, help='make outfile named <vcfTo>.<infoname>.vcf. ', nargs='+', default=None)

parser.add_argument('-b','--blockSize', action="store", dest='blocksize', type=int, help='default blocksize for summary (kb)', nargs='?', default=1)

parser.add_argument('-r','--region', action="store", dest='regions', type=str, help='regions over which to summarize', nargs='+', default=None)

parser.add_argument('-i','--infofield', action="append", dest='info', type=str, help='invo values which should be summarized', default=None)

args = parser.parse_args()
outfile = sys.stdout

if not path.isfile(args.vcfFile+".tbi"):
    #make tabix index for vcf if not tabixed
    if path.splitext(args.vcfFile) != '.gz':
        call(["bgzip",args.vcfFile])
        args.vcfFile += '.gz'
    call(["tabix","-p","vcf",args.vcfFile])
vcfFile = open(args.vcfFile,'r')
reader=vcf.Reader(vcfFile)

infotypes=['Count','Mean']
infos=[]
infotype=dict()
print args.info
for info in args.info:
    print info
    i, t = info.split(':')
    infos += [i]
    infotype[i]=t


#blocksize in kb => block in bases
block=args.blocksize * 1000


#P:RINT HEADER:

print >>outfile, "\t".join(['CHROM','ST','EN']),
for vartype in ['SNP','INDEL']:
    #PRINT TOTAL NO OF SNPs/INDELs IN BLOCK
    print >>outfile, "\t"+vartype,
    for i in infos:
        #IF COUNT, PRINT OUT COUNT AND PC OF TOTAL SNPs/INDELs
        if infotype[i]=='Count':
            print >>outfile,"\t"+i+"\t"+i+'.pc',                
        #IF MEAN, PRINT OUT COUNT AND MEAN
        elif infotype[i]=='Mean':
            print >>outfile,"\t"+i+"\t"+i+'.mn',
print ""


for region in args.regions:
    (chrom, locus) = region.split(":")
    (regst, regen) = locus.split("-")
    regst = int(floor(int(regst) / block)*block)
    regen = int(ceil(int(regen) / block)*block)
    for st in range(regst,regen,block):
        en = st+block

        ##################
        # PROCESS GENOTYPES IN BLOCK
        #################
        blockreader = reader.fetch(chrom, st, en)
        vars = 0
        #setup null count matrix

        counts=dict()
        counts[('SNP','ALL')]=0
        counts[('INDEL','ALL')]=0


        for vartype in ['SNP','INDEL']:
            for i in infos:
                counts[(vartype,i)] = 0
                if infotype[i] == 'Mean':
                    counts[(vartype,i+'sum')] = 0
                    
        for rec in blockreader:
            if rec.is_indel:
                vartype='INDEL'
            elif rec.is_snp:
                vartype='SNP'
            else:
                print >>sys.stderr, rec.CHROM+":"+str(rec.POS)+"-"+rec.REF+'/'.join(rec.ALT)+" is neither SNP or INDEL"
            counts[(vartype,'ALL')] +=1
            for i in infos:
                if i in rec.INFO:
                    counts[(vartype,i)] += 1
                    if infotype[i] == 'Mean':
                        counts[(vartype,i+'sum')] += float(rec.INFO[i])
                                
        
        #PRINT BLOCK POSITIONS
        print "\t".join(map(str,[chrom,st,en])),
        for vartype in ['SNP','INDEL']:
            #PRINT TOTAL NO OF SNPs/INDELs IN BLOCK
            print >>outfile, "\t"+str(counts[(vartype,'ALL')]),
            for i in infos:
                #IF COUNT, PRINT OUT COUNT AND PC OF TOTAL SNPs/INDELs
                if infotype[i]=='Count':
                    count = counts[(vartype,i)]
                    if counts[(vartype,'ALL')] > 0:
                        pc = count/float(counts[(vartype,'ALL')])
                    else:
                        pc=float(0)
                    print >>outfile,"\t"+str(count)+"\t"+str(round(pc,3)),
                    
                #IF MEAN, PRINT OUT COUNT AND MEAN
                elif infotype[i]=='Mean':
                    count = counts[(vartype,i)]
                    if count > 0:
                        mean = counts[(vartype,i+'sum')]/float(count)
                    else:
                        mean=float(0)
                    print >>outfile,"\t"+str(count)+"\t"+str(round(mean,3)),
                    

                elif infotype[i]=='Enumerate':
                    pass
                    #not written this bit yet
        print ""
