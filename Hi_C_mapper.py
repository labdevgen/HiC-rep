#Hi-C mapping class by Minja
#

import Bio.Restriction #for RE site definition
from Bio import SeqIO #to dead with fastq files
import subprocess #to call aligners
import pysam # to deal with sam files
import logging 
import os
from miscellaneous import ensure_path_replacment,startLog
from NGSseq_utils import analyze_sam_file,open_sam_or_bam,default_check_multimap,compare_sam_headers
import random
import numpy as np


#Hi-C-pro based mapping pepline
#Step 1. Map original fastq file
#Step 2. Fileter unique hits
#Step 3. Collect repetitive and not aligned reads
#Step 4. Find RE site
#Step 5. Genereate new file with trimmed reads
#Step 6. Map new reads file

	

class CHi_C_mapper():

	def hic_check_multimap(self,read):
		#should be overrided for each aligner
		return default_check_multimap(read)

	def analyze_sam_file_old(self,inSam): #TODO - remove this func
	#Analyze sam file to return reads ids depending on alignment results.
	#returns dictionary:
	#{"unique":set_of_read_idxs,"repeats":set_of_read_idxs,"unmapped":set_of_read_idxs}
	
	#inSam - filename of input sam file
		raise
		result = {"unique":[],"repeats":[],"unmapped":[]}
		
		alfile = open_sam_or_bam(inSam,"r")
		
		for read in alfile.fetch(until_eof=True):
			if read.is_unmapped:
				result["unmapped"].append(read.query_name)
			elif self.hic_check_multimap(read):
				result["repeats"].append(read.query_name)
			else:
				result["unique"].append(read.query_name)
		logging.info("Analysis of file "+inSam+"\n"+
			"\n".join([i+"\t"+str(len(result[i]))+" reads" for i in result.keys()]))
		return result
	

	def merge_sam_files(self,mainSamFname,secondarySamFname,outSamFname,keepMultiHits=True,keepUnAligned=True,multiHitFile=None,UnAlignedFile=None):
	#Function merges 2sam files
	#It considers that only reads that were not aligned in mainSam are present in secondarySam
	#If these reads are aligned in secondarySam they are taken from this file, else they are removed
	#outSam - file name of resulting file with unique reads
	#If keepMultiHits multiple alignments are kept.
	#When multiHitFile provided multihit kept there, otherwize in outSam
	#If keepUnAligned unaligned reads are kept in outSam.unAligned.sam
	#When UnAlignedFile provided unaligned kept there, otherwize in outSam
		
		#open input files
		mainSam = open_sam_or_bam(mainSamFname,"r")
		secondarySam = open_sam_or_bam(secondarySamFname,"r")
		#chek that fileheaders of main and secondary sam files are similar
		try: 
			header1 = mainSam.header
			header2 = secondarySam.header
			assert compare_sam_headers(header1,header2)
		except:
			logging.error("Sam files do not match")
			raise Exception("Sam files do not math")
		#end chek
		
		#open output files
		outSam = open_sam_or_bam(outSamFname,"wh",template=mainSam)
		
		ext = os.path.splitext(outSamFname)[-1]
		ext = ext if ext != "" else ".sam"
		if keepMultiHits:
			if multiHitFile != None:
				outSamMulti = open_sam_or_bam(multiHitFile,
											"wh",template=mainSam)
			else:
				outSamMulti = outSam
		if keepUnAligned:
			if UnAlignedFile != None:
				outSamUnaligned = open_sam_or_bam(UnAlignedFile,
											"wh",template=mainSam)
			else:
				outSamUnaligned = outSam

			
		#analyze main file
		for read in mainSam.fetch(until_eof=True):
			if not read.is_unmapped:
				if self.hic_check_multimap(read):
					if keepMultiHits:
						outSamMulti.write(read)
				else:
					outSam.write(read)
				
		#analyze secondary file
		for read in secondarySam.fetch(until_eof=True):
			if read.is_unmapped:
				if keepUnAligned:
					outSamUnaligned.write(read)
			elif self.hic_check_multimap(read):
				if keepMultiHits:
					outSamMulti.write(read)
			else:
				outSam.write(read)
		
		#close files
		mainSam.close()
		secondarySam.close()
		outSam.close()
		if keepMultiHits: outSamMulti.close()
		if keepUnAligned: outSamUnaligned.close()


	def generate_trimmed_reads_fastq(self,inFastq,inSam,outFastq,minTrimLen=10,searchFilledInSeq=True):
	#Function creates new Fastq file from inFastq where kept only reads that are
	#not-aligned or not-unique according to inSam
	#reads that are kept searched for enzyme site 
	#and all 5' seq starting from +1 bp after RE site trimmed
	#if enzyme site not found or treamed read length < minTrimLen
	#read won't be kept
	
	#inFastq - input Fastq filename
	#searchFilledInSeq - if true, search for seq originating from RE-site fill-in
	#					 if false, search original RE site
	#inSam - referemce Sam filename
	#outFastq - filename for output Fastq

		alignment_results = analyze_sam_file(inSam,self.hic_check_multimap)
		
		NtreamedReads = 0
		if searchFilledInSeq:
			searchedSeq = self.enzymeFilledInSeq
		else:
			searchedSeq = self.enzymeSeq
		with open(inFastq,"r") as fin, open(outFastq,"w") as fout:
			for line in fin:
				if alignment_results[line[1:].strip().split()[0]]==0: #0 for unmapped read
					seq = fin.next()
					
					trimpos = seq.upper().find(searchedSeq)

					if trimpos<=0 or trimpos+len(searchedSeq) < minTrimLen:
						fin.next()
						fin.next()
					else:
						trimpos += len(searchedSeq)
						fout.write(line)
						fout.write(seq[:trimpos]+"\n")
						fout.write(fin.next())
						fout.write(fin.next()[:trimpos]+"\n")
						NtreamedReads += 1
				else:
					fin.next()
					fin.next()
					fin.next()
		
		logging.info("Number of treamed reads: "+str(NtreamedReads))
			
	def perform_alignmen(self,inFastq,outSam):
	#Perform alignment of inFastq and save info to outSam
			raise Exception("Aligner function is not defined")

	def __init__(self,inFastq,outSam,enzymeName,replaceFiles=False,
					aligner_function=None,aligner_parameters=None,
					logFileName="HiCmapper.log"):
		#inFastq - input fastq file with reads file name
		#outSam - file name of the sam file for alignment
		#aligner_function - function which call aligning software. Function must accepts agrs inFastq and outSam. Optional
		#aligner_parameters - list of additional parameters to be sent to aligner. Paramteres to pass to alignment software
		#example: "--very-sencitive -x path-to-file -p 8"
		#enzymeName - name of RE enzyme used for Hi-C data
		#replaceFiles - replace output and temp files
		#logFileName - file name for logs
		
		startLog(logFileName,loglevel=logging.DEBUG)
		self.statistics = {} #TODO intrpduce special statistics class
		if aligner_function != None:
			self.perform_alignmen==aligner_function
		
		self.aligner_parameters = aligner_parameters
		
		assert enzymeName in Bio.Restriction.Restriction_Dictionary.rest_dict.keys()
		self.enzyme = eval("Bio.Restriction."+enzymeName) #eval is dangerous. Always use check of str that you pass
		self.enzymeSeq = self.enzyme.site.upper()
		st = self.enzyme.elucidate().find("^") #example: EcoRI.elucidate() returns 'G^AATT_C'
		end = self.enzyme.elucidate().find("_")
		self.enzymeFilledInSeq = (self.enzymeSeq[:end-1] + self.enzymeSeq[st:]).upper()
		logging.debug("Uing filled-in enzyme seq: "+self.enzymeFilledInSeq)
		
		
		self.inFastq = inFastq
		self.outSam = outSam
		self.replace_files = replaceFiles
		
		#hic_check_multimap is a function to check wheather read has multiple alignments
		#this function should be defined separately for each aligner
		if not hasattr(self, 'hic_check_multimap'): 
			logging.warning("Using default function to check multiple alignments")
			self.hic_check_multimap = default_check_multimap
	def pepline(self):
		logging.debug("Starting mapping pepline")
		#check that we can (re)write all files
		ensure_path_replacment(self.outSam+".all.sam",self.replace_files) #sam file for initial alignment
		ensure_path_replacment(self.outSam+".trimmed.fq",self.replace_files) #fq for unaligned reads
		ensure_path_replacment(self.outSam+".trimmed.sam",self.replace_files) #sam file for treamed reads alignment
		ensure_path_replacment(self.outSam,self.replace_files) #resulting merged sam file
		
		#align reads
		logging.info("Performing alignment")
		self.perform_alignmen(self.inFastq,self.outSam+".all.sam")
		
		#extract unaligned reads and try to tream them
		logging.info("Generating fastq with treamed reads")
		self.generate_trimmed_reads_fastq(self.inFastq,self.outSam+".all.sam",self.outSam+".trimmed.fq")
		
		#align treamed reads
		logging.info("Aligning treamed reads")
		self.perform_alignmen(self.outSam+".trimmed.fq",self.outSam+".trimmed.sam")
		analyze_sam_file(self.outSam+".trimmed.sam",self.hic_check_multimap)
		
		#merge both sam files and keep multihits
		logging.info("Merging sam files")
		self.merge_sam_files(self.outSam+".all.sam",self.outSam+".trimmed.sam",self.outSam,keepMultiHits=True,keepUnAligned=True)


class Cbowtie2_Hi_C_mapper(CHi_C_mapper):
	def hic_check_multimap(self,read):
		#for bowtie the best way to find multiple alignments is to check AS:Xi optional flag
		#which provides a score for the best-scoring alignment found other than the alignment reported
		return read.has_tag("XS:i")
		
	def perform_alignmen(self,inFastq,outSam):
		if self.aligner_parameters != None:
			call_args = ["bowtie2","-U",inFastq,"-S",outSam]+self.aligner_parameters.split()
		else:
			raise Exception("Please provide bowtie2 reference as -x parameter")
		logging.info("Starting bowtie2 with command:\n"+" ".join(call_args))
		out = subprocess.check_output(call_args,stderr=subprocess.STDOUT)
		logging.info("Bowtie2 output:\n"+out)

class Cbwa_Hi_C_mapper(CHi_C_mapper):
	def hic_check_multimap(self,read):
		if	read.has_tag("XT"):
			return read.get_tag("XT")=="R"
		else:
			return False

	def perform_alignmen(self,inFastq,outSam):
		ensure_path_replacment(outSam+".sai",self.replace_files) #for bwa we additionally need sai files
		if self.aligner_parameters != None:
			if not isinstance(self.aligner_parameters, basestring):
				call_args_aln = ["bwa","aln"]+self.aligner_parameters[0].split()+[inFastq,outSam+".sai"]
				call_args_samse = ["bwa","samse"]+self.aligner_parameters[1].split()+[outSam+".sai",inFastq]
			else:
				logging.warning("bwa mapper uses 2 commands for alignment: aln and samse \
				you only provided sigle args for alignment. \
				It will be considered that both aln and samse shoud revieve these args")
		else:
			raise Exception("Please provide bwa reference.fasta as additional parameter")
		logging.info("Starting bwa with command:\n"+" ".join(call_args_aln))
		with open(outSam+".sai","w") as out_file:
			out = subprocess.check_call(call_args_aln,stdout=out_file)
		with open(outSam,"w") as out_file:
			out = subprocess.check_call(call_args_samse,stdout=out_file)
		
		logging.info("bwa executed successfully")
		
class CRepeats_Localizer():
	#the main purpos of this class is to localize repeatitve reads
	#Instance of the class got 2 main variables: genome and 2 sam files with (repetitive) reads
	#For repeatitve reads, sam file contains multiple (all possible) alignments
	#Using the second sam file, best alignment in first sam file will be found
	#(best means with minimal distance between reads)
	
	def __init__(self,R1sam,R2sam,replaceFiles = True,logFileName=None):
		#R1sam, R2sam - sam files with left and right reads
		#genome - the mirnylab genome instance
		#replaceFiles - ensure that output files can be rewritten if exist
		self.R1sam = R1sam
		self.R2sam = R2sam
		self.replaceFiles = replaceFiles
		startLog(filename=logFileName)
		logging.info("Starting repeat localizer")
		random.seed()
	
	def check_header(self,header,genome):
		#check that reference names in header are same as chromosome labels in genome
		return sum([i in header['SQ'] for i in genome.label2idx.keys()])==len(header2['SQ'])==len(genome.label2idx)
	
	def findBestAlignments(self,fout):
		#fout - name of the output file
		
		maxdist = 100000000000
		
		ensure_path_replacment(fout,self.replaceFiles)
		
		R1 = open_sam_or_bam(self.R1sam,"r")
		R2 = open_sam_or_bam(self.R2sam,"r")
		
		#chek that fileheaders of main and secondary sam files are similar
		#and in accordance with the genome object
		try: 
			header1 = R1.header
			header2 = R2.header
			assert compare_sam_headers(header1,header2)
		except:
			logging.error("Sam files do not match")
			raise Exception("Sam files do not math")
		#end chek
		
		logging.debug("Reading "+self.R1sam)
		##############Process R1 sam file and keep aligning info#############
		self.reads1 = {} #dict of tuples:
		#read_name: (r1_ref,r1_start)
		
		self.read_pairs={} #dict of tuples:
		#read_name: (r1_ref,r1_start,r2_ref,r2_start,dist)

		self.R1_unmapped = 0
		for r1 in R1.fetch(until_eof=True):
			if not r1.is_unmapped:
				try:
					self.reads1[r1.query_name].append((r1.reference_id,r1.reference_start))
				except:
					self.reads1[r1.query_name]=[(r1.reference_id,r1.reference_start)]
					self.read_pairs[r1.query_name]=(-1,-1,-1,-1,maxdist)
			else:
				self.R1_unmapped += 1

		logging.debug("Reading "+self.R2sam)
		##############Process R2 sam file and find best alignments################
		self.SE = 0
		self.R2_unmapped = 0
		for r2 in R2.fetch(until_eof=True):
			if r2.is_unmapped:
				self.R2_unmapped += 1
				continue
			
			try:	 
				current_dist = self.read_pairs[r2.query_name][4]
			except:
				self.SE += 1
				continue

			for r1 in self.reads1[r2.query_name]: #go through all r1 locations
				if r2.reference_id != r1[0]: #trans read
					if current_dist == maxdist: #if we don't have any pare yet
						self.read_pairs[r2.query_name] = r1+(r2.reference_id,r2.reference_start,-1) #set current
					elif current_dist == -1 and random.randint(0,1) == 0: #if we have trans pair, with 50% chance
						self.read_pairs[r2.query_name] = r1+(r2.reference_id,r2.reference_start,-1) #set current

					current_dist = -1 #set current distance to trans

				else: #cis read
					dist = abs(r1[1]-r2.reference_start) #calculate distance in nt
					if dist < current_dist or current_dist == -1: #if it's closer than we had
						self.read_pairs[r2.query_name] = r1+(r2.reference_id,r2.reference_start,-1) #set current
						current_dist = dist #update current distance

		##############Write results to the output################
		
		logging.debug("Writing results to "+fout)
		
		with open(fout,"w") as out:
			for i in self.read_pairs.keys():
				if self.read_pairs[i][-1] != maxdist:
					out.write(i+"\t"+"\t".join(map(str,self.read_pairs[i]))+"\n")
		
		logging.debug("Done")
		
	def print_stats(self):
		try:
			self.SE
		except:
			logging.error("Please run findBestAlignments first. \
			This call to print_stats will have no effect")
		else:
			logging.info("In R1: "+str(self.R1_unmapped)+" unmapped")
			logging.info("In R2: "+str(self.R2_unmapped)+" unmapped")
			logging.info("In R1: "+str(self.SE)+" single end")
			logging.info("In R2: "+str(len(self.reads1)-len(self.read_pairs))+" single end")
			logging.info("Destibution of alignments number per read in R1")
			c=np.bincount([len(i) for i in self.reads1])
			logging.info("\n".join([str(ind)+":"+str(val) \
						for ind,val in enumerate(c) if val != 0]))