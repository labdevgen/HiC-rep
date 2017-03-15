import pysam # to deal with sam files
import logging
from miscellaneous import startLog

def open_sam_or_bam(filename,mode="r",template=None):
	if filename.endswith(".bam") and mode.find("b")<0:
		if mode == "wh":
			logging.warning("'wh' mode cannot be used for bam file. Changing to 'wb' mode")
		mode[1] = "b"
	
	if template!=None:
		if mode[0]=="w":
			return pysam.AlignmentFile(filename,mode,template)
		else:
			logging.warn("Template provided for read mode, ignoring")
			return pysam.AlignmentFile(filename,mode)
	else:
		return pysam.AlignmentFile(filename,mode)


def default_check_multimap(read):
		#What should be considered as multi-map is a question
		#defaul is to use mapq field, with values 0 or 1.
		#however this function shell be override for each aligner
		return read.mapping_quality >= 2


def sam_file_analyzer_generator(inSam,check_multimapFunction=default_check_multimap):
		#Analyze inSam file to classify reads depending on alignment results.

		#What should be considered as multi-map is a question
		#defaul is to use default_check_multimap function, see above
		#however this function shell be overrided for each aligner

		#returns generator of tuples:
		#{read:status}
		#where statuses are 
		#1 = "unique"
		#2 = "repeats"
		#0 = "unmapped"
		#and read in Pysam AlignmentObject
		
		#inSam - filename of input sam file

		alfile = open_sam_or_bam(inSam,"r")
		
		for read in alfile.fetch(until_eof=True):
			if read.is_unmapped:
				yield (read,0)
			elif check_multimapFunction(read):
				yield (read,2)
			else:
				yield (read,1)

		
def analyze_sam_file(inSam,check_multimapFunction=default_check_multimap):
		#Analyze inSam file to classify reads ids depending on alignment results.

		#What should be considered as multi-map is a question
		#defaul is to use default_check_multimap function, see above
		#however this function shell be overrided for each aligner

		#returns dictionary:
		#{read:status}
		#where statuses are 
		#1 = "unique"
		#2 = "repeats"
		#0 = "unmapped"
		
		#inSam - filename of input sam file
		result = {}
		statistics = {}
		statistics[inSam]={"unmapped":0,"repeat":0,"unique":0}
		categories = {0:"unmapped",1:"unique",2:"repeat"}
		alfile = open_sam_or_bam(inSam,"r")
		
		for read,category in sam_file_analyzer_generator(inSam,check_multimapFunction):
			result[read.query_name] = category
			statistics[inSam][categories[category]] += 1

		total = float(sum(statistics[inSam].values()))
		if total == 0:
			logging.info(inSam+": Sam file with 0 reads")
		else:
			logging.info("Analysis of file "+inSam+"\n"+
				"\n".join([i+"\t"+str(statistics[inSam][i])+"("+str(statistics[inSam][i]*100/total)+\
				"%) reads" for i in statistics[inSam].keys()]))
		return result

		


def splif_fastq_depending_on_alignment_results(inFastq,inSam,outputs,check_multimapFunction=default_check_multimap):
	#Function creates new Fastq file from inFastq where kept only reads that are
	#in accordance with alignment class (based on inSam)
	
	#inFastq - input Fastq filename
	#inSam - referemce Sam filename
	#outputs - dict
	#alignment_class : handler
	#where alignment_classes are 
	#1 = "unique"
	#2 = "repeats"
	#0 = "unmapped"
	#and handler is either file handler or
	#None to throw away reads
	#
	#see analyze_sam_file for check_multimapFunction parameter
		
		startLog(None)
		alignment_results = analyze_sam_file(inSam,check_multimapFunction=check_multimapFunction)
		stats={}
		for i in [0,1,2]:
			if not i in outputs.keys():
					outputs[i] = None
			else:
				stats[i] = 0
		
		with open(inFastq,"r") as fin:
			for line in fin:
				alignment_class = alignment_results[line[1:].strip().split()[0]]
				if outputs[alignment_class]!=None:
						stats[alignment_class] += 1
						outputs[alignment_class].write(line)
						outputs[alignment_class].write(fin.next())
						outputs[alignment_class].write(fin.next())
						outputs[alignment_class].write(fin.next())
				else:
					fin.next()
					fin.next()
					fin.next()
		
		logging.info("Fastq "+inFastq+" filtered:\n"+
		"\n".join([str(i)+":\t"+str(j) for i,j in stats.iteritems()]))

def generate_fastq_depending_on_alignment_results(inSam,outputs,check_multimapFunction=default_check_multimap):
	#Function creates new Fastq file from inSam where kept only reads that are
	#in accordance with alignment class
	
	#inSam - referemce Sam filename
	#outputs - dict {alignment_class : handler}
	#where alignment_classes are 
	#1 = "unique"
	#2 = "repeats"
	#0 = "unmapped"
	#and handler is either file handler or
	#None to throw away reads
	#
	#see analyze_sam_file for check_multimapFunction parameter
		
		startLog(None)
		stats = dict([(i,0) for i in [0,1,2]])
		for read,category in sam_file_analyzer_generator(inSam,check_multimapFunction):
			stats[category] += 1
			if category in outputs.keys():
				outputs[category].write("@"+read.query_name+"\n"+
										read.query_sequence+
										"\n+\n")
				for q in read.query_qualities:
					outputs[category].write(chr(q+33))
				outputs[category].write("\n")
		logging.info("Sam file "+inSam+" filtered:\n"+
		"\n".join([str(i)+":\t"+str(j) for i,j in stats.iteritems()]))


def trim_reads_in_fq(inFastq,Nbases,outFastq,minlen=10):
#Trim several bases from each read in inFastq
#inFastq - filename of input fastq
#Nbases - number of basese to trim
#outFastq - output fastq filename

	with open(inFastq) as fin, open(outFastq,"w") as fout:
		name = fin.readline()
		while name:
			if name == "\n" or name == "":
				logging.warning("Emty line in file "+inFastq+" ==>stop iteration")
				break
			bases=fin.readline()
			if len(bases)>Nbases+minlen:
				fout.write(name)
				fout.write(bases[Nbases+1:])
				fout.write(fin.readline())
				fout.write(fin.readline()[Nbases+1:])
			else:
				fin.readline()
				fin.readline()
			name = fin.readline()
			
def compare_sam_headers(header1,header2):
	#compare headers of two sam files
	#returns True if references names (starting with @SQ) are similar in both headers
	
	return sum([i in header2['SQ'] for i in header1['SQ']])==len(header2['SQ'])==len(header1['SQ'])