from Hi_C_mapper import *
from miscellaneous import startPmDebug
from NGSseq_utils import *
startPmDebug()

basedir = "/mnt/storage/home/vsfishman/tmp/HiC_repeats/"

inFastq = basedir + "input/Battulin2015/2500k/Fib_full2.2500k.fq_R1.fastq"
#inFastq = basedir + "input/Battulin2015/Fib_full2.5k.fq_1.fq"
#inFastq = basedir + "input/Battulin2015/mmSample.fq"

#outSam = basedir + "out/Battulin2015/test_bwa.sam"
outSam = basedir + "out/Battulin2015/2500k_bowtie2_R1.sam"
enzymeName="HindIII"

mapper = Cbowtie2_Hi_C_mapper(inFastq,
					outSam,
					enzymeName,
					replaceFiles=True,
					aligner_function=None,
					aligner_parameters="--very-sensitive -p 10 -x "+basedir+"fasta/mm10",
					logFileName=basedir+"/logs/HiCmapper.log")

# mapper = Cbwa_Hi_C_mapper(inFastq,
					# outSam,
					# enzymeName,
					# replaceFiles=True,
					# aligner_function=None,
					# aligner_parameters=["-t 7 "+basedir+"fasta/mm10.fa",basedir+"fasta/mm10.fa"],
					# logFileName=basedir+"/logs/HiCmapper.log")
#mapper.pepline()

#now let's trim some reads and re-align it
from NGSseq_utils import splif_fastq_depending_on_alignment_results

unaligned = open(basedir+"out/Battulin2015/test.sam.all.unaligned.fq","w")
#unaligned = open(outSam+".unaligned.fq","w")
unique =  open(outSam+".unique.fq","w")
generate_fastq_depending_on_alignment_results(outSam,
						{1:unique},
						check_multimapFunction=mapper.hic_check_multimap)
#unaligned.close()
unique.close()

trim_reads_in_fq(outSam+".unique.fq",30,outSam+".unique.fq"+".trim30.fq")
mapper.perform_alignmen(outSam+".unique.fq"+".trim30.fq",outSam+".unique.fq"+".trim30.fq"+".realigned.fq")