#!/bin/bash

################################################################################
## mothurGreatLakes.sh
##  A shell script with SLURM directives to run mothur 16S microbiome analysis
##  jobs on the UM ARC-TS Great Lakes cluster with minimal user intervention.
##  This script is specifically tailored for running by Dickson Lab members but
##  can easily be tweaked for alternative uses. This script follows Pat Schloss'
##  Mothur MiSeq_SOP (https://www.mothur.org/wiki/MiSeq_SOP) very closely.
##
##  If you are seeing this section for the first time, I will strongly suggest
##  looking at the GitHub Wiki page (link below) for instructions to best
##  use this script.
##
## Author: piyuranjan\@gmail.com
## Source: https://github.com/piyuranjan/DicksonLabScripts/blob/master/mothurGreatLakes.sh
## Wiki: https://github.com/piyuranjan/DicksonLabScripts/wiki/mothurGreatLakes.sh
################################################################################

## SBatch preamble ##

#SBATCH --job-name=mothurGL
#SBATCH --account=rodickso
#SBATCH --mail-type=ALL
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
##SBATCH --mem=20G
#SBATCH --mem-per-cpu=2G
#SBATCH --time=1-00:00:00
#SBATCH --export=ALL

## End SBatch Preamble ##

### Listing compute nodes allocated to the job
if [[ $SLURM_JOB_NODELIST ]] ; then
	echo "Running on"
	scontrol show hostnames $SLURM_JOB_NODELIST
fi
### Adding Conda init to subshell
if [ -f "$(conda info --base)/etc/profile.d/conda.sh" ]; then
	. "$(conda info --base)/etc/profile.d/conda.sh"
fi

## Job Commands ##

### Variables that can be overridden on the command line
if [ -z $FQDIR ]; then #set fastq directory for finding MiSeq files
	FQDIR='.'
else
	FQDIR=${FQDIR%/}
fi
if [ -z $PROC ] ; then #Set procs/threads for the following commands to use; default: Job-CPUs
	PROC=$SLURM_CPUS_PER_TASK
fi
if [ -z $MOTHUROUTPREFIX ]; then #Set Mothur output file prefix; default: Job-Name
	MOTHUROUTPREFIX=$SLURM_JOB_NAME
fi
if [ -z $SILVAPATH ]; then #Set Silva database path
	SILVAPATH='/nfs/turbo/umms-rodickso/Databases/MothurDB/SilvaV132/silva.nr_v132.regionV4.align'
fi
if [ -z $RDPREFPATH ]; then #Set RDP database fasta path
	RDPREFPATH='/nfs/turbo/umms-rodickso/Databases/MothurDB/RDP16/trainset16_022016.rdp.fasta'
fi
if [ -z $RDPTAXPATH ]; then #Set RDP database taxonomy path
	RDPTAXPATH='/nfs/turbo/umms-rodickso/Databases/MothurDB/RDP16/trainset16_022016.rdp.tax'
fi

### Activate Mothur and check the package version
conda activate envMothur #Deactivate this line if you have Mothur installed differently and don't like it throwing an error statement in your logs. It will not impede the job though.
which mothur
if [[ $? -ne 0 ]]; then #Kill job if mothur cannot be seen on PATH.
	>&2 echo "[ERROR] Problem finding mothur."
	exit 1
fi
mothur -v|perl -pe 's/^.*Linux/Linux/;'

### Extract files and prepare them for use
parallel -j $PROC 'pigz -dc {} >{/.}' ::: $FQDIR/*.fastq.gz #Will extract .gz files in parallel to PWD
ls -1 *.fastq|perl -ne 'chomp;$n=$_;$n=~s/-//g;print `mv $_ $n\n`;' #Rename all fastq to remove '-' from name

### Execute the full Mothur workload
/usr/bin/time -f "\nCommand stats:\nProc:\tElapsed Time = %E,\tPerc CPU = %P,\nMem:\tAvg Total Mem = %KKB,\tPeak Mem = %MKB,\nExit Status: %x" \
mothur "#\
set.current(inputdir=.,outputdir=.,processors=$PROC);\
make.file(type=fastq, prefix=$MOTHUROUTPREFIX, delim=*);\
make.contigs(file=$MOTHUROUTPREFIX.files);\
summary.seqs();\
screen.seqs(fasta=current, maxambig=0, maxlength=275);\
unique.seqs();\
count.seqs(name=current, group=current);\
align.seqs(fasta=current, reference=$SILVAPATH);\
screen.seqs(fasta=current, count=current, start=1968, end=11550, maxhomop=8);\
summary.seqs(fasta=current, count=current);\
filter.seqs(fasta=current, vertical=T, trump=.);\
unique.seqs(fasta=current, count=current);\
pre.cluster(fasta=current, count=current, diffs=2);\
chimera.vsearch(fasta=current, count=current, dereplicate=t);\
remove.seqs(fasta=current, accnos=current, name=current);\
summary.seqs(fasta=current, count=current);\
classify.seqs(fasta=current, count=current, reference=$RDPREFPATH, taxonomy=$RDPTAXPATH, cutoff=80);\
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);\
summary.tax(taxonomy=current, count=current);\
dist.seqs(fasta=current, cutoff=0.03);\
cluster(column=current, count=current, cutoff=0.03);\
make.shared(list=current, count=current, label=0.03);\
classify.otu(list=current, count=current, taxonomy=current, label=0.03);\
get.oturep(column=current, count=current, fasta=current, list=current);\
get.current();\
"|perl -pe 's/^.*Linux/Linux/;'

### Copy important files with more legible names
OUTDIR="ResultFiles"
mkdir $OUTDIR
SHARED="$OUTDIR/$MOTHUROUTPREFIX.shared"
cp $MOTHUROUTPREFIX.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared $SHARED
TAXONOMY="$OUTDIR/$MOTHUROUTPREFIX.cons.taxonomy"
cp $MOTHUROUTPREFIX.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy $TAXONOMY
REP="$OUTDIR/$MOTHUROUTPREFIX.rep.fasta"
cp $MOTHUROUTPREFIX.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.fasta $REP
MOTHURLOG=$(ls -1 mothur.*.logfile|tail -1)
LOG="$OUTDIR/$MOTHUROUTPREFIX.log"
cp $MOTHURLOG $LOG
if [ -f "$SHARED" -a -f "$TAXONOMY" -a -f "$REP" -a -f "$LOG" ]; then
	FINISH="\n------\
	\n\nThis job has finished successfully: $SLURM_JOB_NAME\
	\n\nDirect Mothur log for this run can be found in: $MOTHURLOG\
	\n\nJob logs, after this job finishes can be found in:\
	\n- $SLURM_JOB_NAME-$SLURM_JOB_ID.out\
	\n- $SLURM_JOB_NAME-$SLURM_JOB_ID.err\
	\n\n---\
	\n\nThis Mothur run was processed with the following database paths:\
	\n- $SILVAPATH \n- $RDPREFPATH \n- $RDPTAXPATH\
	\n\nFiles after mothur processing are:\
	\n- $SHARED \n- $TAXONOMY \n- $REP \n- $LOG"
	echo -e $FINISH >>$LOG
	echo -e $FINISH
else
	ERR="\n------\
	\n\nERROR: Somthing has gone wrong with the mothur execution! Please check:\
	\n- $SLURM_JOB_NAME-$SLURM_JOB_ID.out \n- $SLURM_JOB_NAME-$SLURM_JOB_ID.err \n- $MOTHURLOG"
	echo -e $ERR >>$LOG
	echo -e $ERR
fi


## End Job Commands ##
