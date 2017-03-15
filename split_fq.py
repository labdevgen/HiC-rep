# -*- coding: utf-8 -*-
"""
Created on Wed Feb 08 19:07:50 2017

@author: FishmanVS
"""
from Bio import SeqIO
import sys

fname=sys.argv[1]
read_1 = []
read_2 = []

for read in SeqIO.parse(fname,"fastq"):
    read_1.append(read[:50])
    read_2.append(read[50:])

SeqIO.write(read_1,open(fname+"_1.fq","w"),"fastq")
SeqIO.write(read_2,open(fname+"_2.fq","w"),"fastq")
print "Done"