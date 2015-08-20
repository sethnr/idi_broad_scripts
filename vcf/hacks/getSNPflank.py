#!/usr/bin/python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
import os.path as path
import gzip
import argparse

########
# get some vars
########
parser = argparse.ArgumentParser(description='build index of genome matches in 3k blocks')
parser.add_argument('-f','--fasta', action="store", dest='fasta', type=str, help='target genome to index (fasta)', nargs='?')
parser.add_argument('-n','--SNP', action="store", dest='snpLocs', type=str, help='bed file for SNP locations', nargs='?', default=None)
parser.add_argument('-d','--flank','--flank_distance', action="store", dest='flank', type=int, help='region around SNP to parse', nargs='?')

args = parser.parse_args()

SNPs=[]
if args.snpLocs is not None:
    for line in open(args.snpLocs,'r'):
        F = line.strip().split("\t")
        [chrom, pos] = F[:2]
        pos = int(pos)
        SNPs += [(chrom,pos)]


if fasta[-2:] == "gz": 
  seqfile = gzip.open(fasta)
else:
  seqfile = open(fasta)

seq = SeqIO.parse(seqfile,'fasta')

for chrom in seq:
  seqlen = len(chrom.seq)
  cSNPs = [(c,p) for (c,p) in SNPs where c==chr.name]
  
  print chr.name, seqlen
  for (c,p) in cSNPs:
    start = p-args.flank
    end = p+args.flank
    title = chr.name+":"+
    seq = chr.seq[start:end]
    print >>sys.stdout, ">".title."\n".seq

