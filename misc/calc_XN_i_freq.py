import sys
import numpy as np
f = open(sys.argv[1])

nmbrs = []

count = 0
print "Analysing file ",sys.argv[1]
for line in f:
	if "NH:i:" in line:
		line = line.strip().split("NH:i:")
		nmbrs.append(int(line[-1].split()[0]))
	count += 1
#	if count > 100000:
#		break
print "Counting"
nmbrs = np.bincount(nmbrs)

print "Done"
if nmbrs[0] != 0:
	raise
nmbrs = np.array([val/(ind+1.) for ind,val in enumerate(nmbrs[1:])])
np.savetxt(sys.argv[1]+".xinumb",nmbrs)