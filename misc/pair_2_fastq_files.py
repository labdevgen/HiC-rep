from Bio import SeqIO
import sys
import numpy as np

R1="../../out/Battulin2015/2500k_bowtie2_R1.sam.unique.fq.trim30.fq"
#R2="../../out/Battulin2015/2500k_bowtie2_R2.sam.unique.fq.trim30.fq"
R2="../../out/Battulin2015/2500k_bowtie2_R2.sam.unique.fq"
in1 = SeqIO.parse(open(R1),"fastq")
in2 = SeqIO.parse(open(R2),"fastq")



in1 = SeqIO.parse(open(R1),"fastq")
reads1=dict([(r1.id,r1) for r1 in in1])
in2 = SeqIO.parse(open(R2),"fastq")
reads2=dict([(r2.id,r2) for r2 in in2])

shared=np.intersect1d(reads1.keys(),reads2.keys())

SeqIO.write([reads1[i] for i in shared],open(R1+".paired","w"),"fastq")
SeqIO.write([reads2[i] for i in shared],open(R2+".paired","w"),"fastq")