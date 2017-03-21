from Hi_C_mapper import *
from miscellaneous import startPmDebug
from NGSseq_utils import *
startPmDebug()

basedir = "/mnt/storage/home/vsfishman/tmp/HiC_repeats/"

inFastq = basedir + "input/Battulin2015/2500k/Fib_full2.2500k.fq_R2.fastq"
#inFastq = basedir + "input/Battulin2015/Fib_full2.5k.fq_1.fq"
#inFastq = basedir + "input/Battulin2015/mmSample.fq"

#outSam = basedir + "out/Battulin2015/test_bwa.sam"
outSam = basedir + "out/Battulin2015/2500k_bowtie2_R2.sam"
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

# now let's trim some reads and re-align it
from NGSseq_utils import splif_fastq_depending_on_alignment_results

#unique =  open(outSam+".unique.fq","w")
repeats =  open(outSam+".repeat.fq","w")
generate_fastq_depending_on_alignment_results(outSam,
						{2:repeats},
						check_multimapFunction=mapper.hic_check_multimap)
repeats.close()

#trim_reads_in_fq(outSam+".unique.fq",30,outSam+".unique.fq"+".trim30.fq")
#mapper.aligner_parameters += " -a"
#mapper.perform_alignmen(outSam+".unique.fq"+".trim30.fq",outSam+".unique.fq"+".trim30.fq"+".realigned.sam")