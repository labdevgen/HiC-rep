import sys
sys.path.append("/mnt/storage/home/vsfishman/tmp/HiC_repeats/scripts")
from miscellaneous import startPmDebug
startPmDebug()



fname = sys.argv[1]

with open(fname) as fin:
	fin.readline()
	fin.readline()
	fin.readline()
	reps = []
	for line in fin:
		line = line.split()
		line = [l for l in line if l.strip() != "" ]
		reps.append(line[9].upper())
	reps=set(reps)
	
print len(reps)

fname_repMasker = sys.argv[2]
with open(fname_repMasker) as fin:
	repMaskerReps=[]
	for line in fin:
		if line.startswith(">"):
			repMaskerReps.append(line.split()[0][1:].upper())
			
print len(reps),sum([i in repMaskerReps for i in reps])
for i in reps:
	if not i in repMaskerReps:
		print i