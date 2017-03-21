import numpy as np
import sys
sys.path.append("/mnt/storage/home/vsfishman/tmp/HiC_repeats/scripts")
from miscellaneous import startPmDebug,startLog
from NGSseq_utils import *
from Hi_C_mapper import *
startPmDebug()

basedir = "/mnt/storage/home/vsfishman/tmp/HiC_repeats/out/Battulin2015/"

uniqueSam_R1 = basedir+"2500k_bowtie2_R1.sam"
uniqueSam_R2 = basedir+"2500k_bowtie2_R2.sam"
guesResults = basedir+"2500k.unique.fq.trim30.fq.realigned.coords"
guesResults = basedir + "2500k"+".unique.fq"+".trim30.fq"+"R2Untrimmed.paired.coords"

startLog(filename=None)
mapper = Cbowtie2_Hi_C_mapper("","","HindIII",logFileName=None)

reads_pairs = {}

for line in open(guesResults):
	line = line.strip().split()
	reads_pairs[line[0]] = map(int,line[1:])

dist1 = []
	
r1notFound = 0
r2notFound = 0
	
for r1,status in sam_file_analyzer_generator(uniqueSam_R1,check_multimapFunction=mapper.hic_check_multimap):
	if status == 1:
		try:
			if reads_pairs[r1.query_name][0] != r1.reference_id:
				dist1.append(-1)
			else:
				dist1.append(abs(r1.reference_start-reads_pairs[r1.query_name][1]))
		except:
			r1notFound += 1

dist2 = []
for r2,status in sam_file_analyzer_generator(uniqueSam_R2,check_multimapFunction=mapper.hic_check_multimap):
	if status == 1:
		try:
			if reads_pairs[r2.query_name][2] != r2.reference_id:
				dist2.append(-1)
			else:
				dist2.append(abs(r2.reference_start-reads_pairs[r2.query_name][3]))
		except:
			r2notFound += 1
		
np.savetxt(basedir+"dist1_R2Untrimmed.txt",dist1)
np.savetxt(basedir+"dist2_R2Untrimmed.txt",dist2)
print dist1,dist2
print r1notFound,r2notFound
