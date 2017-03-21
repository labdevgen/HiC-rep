from Hi_C_mapper import *
from miscellaneous import startPmDebug
from NGSseq_utils import *
startPmDebug()

basedir = "/mnt/storage/home/vsfishman/tmp/HiC_repeats/"

R1sam = basedir + "out/Battulin2015/2500k_bowtie2_R1.sam"+".unique.fq"+".trim30.fq.paired"+".sam"
#R2sam = basedir + "out/Battulin2015/2500k_bowtie2_R2.sam"+".unique.fq"+".trim30.fq.paired"+".sam"
R2sam = basedir + "out/Battulin2015/2500k_bowtie2_R2.sam.unique.paired.realigned.sam"

RepeatLocalizer = CRepeats_Localizer(R1sam,R2sam,logFileName=basedir+"logs/HiC_rep_localizer.log")
RepeatLocalizer.findBestAlignments(basedir + "out/Battulin2015/2500k"+".unique.fq"+".trim30.fq"+"R2Untrimmed.paired.coords")
RepeatLocalizer.print_stats()

