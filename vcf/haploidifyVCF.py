#!/usr/bin/python

import vcf
import sys
import argparse
from string import *
import collections
import copy

parser = argparse.ArgumentParser(description='get STR length differences from disco VCF')

parser.add_argument('-v','--vcf', action="store", dest='vcfFile1', type=str, help='vcfFile1', nargs='?', default=None)

args = parser.parse_args()

vcfFile1 = open(args.vcfFile1,'r')
reader1=vcf.Reader(vcfFile1)

vcfoutF1 = replace(args.vcfFile1,'.vcf','.HAPLOID.vcf')
vcfoutF1 = replace(vcfoutF1,'.vcf.gz','.vcf')
print >>sys.stdout, "VCFOUT1",vcfoutF1
vcfoutF1 = open(vcfoutF1,'w')

reader2 = copy.copy(reader1)

vcfout1=vcf.Writer(vcfoutF1,reader2)

makeCallData = vcf.model.make_calldata_tuple(("GT","ALTP","REFP","GP"))

print >>sys.stderr, reader1.formats


GTindex=None
for rec in reader1:
#    calls = dict()
#    rec.samples = None
    newSamples = list()
    GTindex=rec.FORMAT.split(":").index("GT")
    GPindex=rec.FORMAT.split(":").index("GP")
    format = rec.FORMAT.split(":")
    makeCallData = vcf.model.make_calldata_tuple(tuple(rec.FORMAT.split(":")))
    
    for call in rec.samples:
#        print call.data
#        print rec.FORMAT
        callData= list(call.data)
        dipGT = callData[GTindex]
        
        if dipGT is not None:
            dipAlleles = dipGT.split('/')
            callData[GTindex]
            callData[GTindex]=str(dipAlleles[0])
        
        dipGP=callData[GPindex]
        if dipGP is not None:
            if len(dipGP) > 3:
                print >>sys.stderr, rec.CHROM,rec.POS,len(rec.ALT),dipGP 
            newGP = list()
            for i in range(0,len(rec.ALT)+1):
                iiIndex=(i*(i+1)/2)+i
                newGP += [dipGP[iiIndex]]
            callData[GPindex]=newGP

        newCallString = ""

        for i in range(0,len(callData)):
            if callData[i] is None: callData[i]='.'
            if callData[i] == 'None': callData[i]='.'
            if type(callData[i]) is list: 
                for j in range(0,len(callData[i])):
                    if callData[i][j] is None:  callData[i][j]='.'
                callData[i]=",".join(map(str,callData[i]))
            newCallString += format[i]+'="'+str(callData[i])+'", '
        newCallData = eval("makeCallData("+newCallString+")")
        
#        print call.data
#        print callData
#        print newCallData
        newCall = vcf.model._Call(rec,call.sample,newCallData) 
        newSamples += [newCall]
#    print rec.POS, rec.samples
    rec.samples=newSamples
#    print rec.POS, rec.samples

    vcfout1.write_record(rec)
    
#makeCallData(namedtuple(callData))
        

        #print >> sys.stderr,  call.sample
#        if sampleKey[call.sample] not in newSampleNames:
#            newSampleNames += [sampleKey[call.sample]]
#            call.sample=sampleKey[call.sample]
#            newSamples += [call]
    #print >>sys.stderr, len(rec.samples), len(newSamples)
#    rec.samples = newSamples
