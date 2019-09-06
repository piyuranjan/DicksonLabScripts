#!/bin/bash

################################################################################
## mothurGreatLakes.sh
##  A shell script with SLURM directives to run mothur microbial analysis jobs
##  on the UM ARC-TS Great Lakes cluster with minimal user intervention.
##  This script is specifically tailored for running by Dickson Lab members but
##  can easily be tweaked for alternative uses. This script follows Pat Schloss'
##  Mothur MiSeq_SOP (https://www.mothur.org/wiki/MiSeq_SOP) very closely.
##
## Usage:
##  Newbie: Submit a batch of fastq.gz from 16S Illumina run; recommended
##    only for a small run
##    sbatch mothurGreatLakes.sh
##  Techie: Submit a batch of fastq.gz from a large 16S Illumina run
##    sbatch --export=ALL,FQDIR="../MiSeq_M02127_2019/" --job-name=ardsLungs mothurGreatLakes.sh
##  Expert: Use a full node (36 cores) for the run
##    sbatch --export=ALL,PROC=36,FQDIR="../MiSeq_M02127_2019/" --cpus-per-task=36 --job-name=ardsLungs mothurGreatLakes.sh
##  Wizard: 
##    sbatch --export=ALL,PROC=36,FQDIR="../MiSeq_M02127_2019/" -c 18 -J ardsLungs --mem 20G -a ghuff -p largemem -t 6:00:00 mothurGreatLakes.sh
##
## Input:
##  The script is pre-configured to pull-in all information that it needs but
##    users can make the following choices or custom overrides if they wish.
##    After these modifications, users are free to override any SLURM/SBATCH
##    directive on the command line while submission.
## --job-name=<prefix>
##     Default: mothurGL; This prefix will be the name of the Great Lakes job
##       and the prefix of the output files. User override strongly suggested.
## --export=ALL,FQDIR=<Path/To/Gzipped/Fastq>
##     Default: PWD; Path to gzipped fastq files. User override strongly
##       suggested. Folder for fastq.gz should not be where you run mothur.
## --export=ALL,PROC=<#Processors>
##     Default: 10; Processor threads to run the job with. Change this
##       parameter if desired. Please be sure to pair it with the SLURM
##       override for processors --cpus-per-task=<#Processors>.
## --cpus-per-task=<#Processors>
##     Default: 10; Processors for SLURM to allocate to the job. Please be
##       sure to pair it with the PROC variable for the mothur command with
##       --export=ALL,PROC=<#Processors>.
## 
## Output:
##  Once the cluster job finishes, it will prepare the following three files
##    as output along with all the intermediate and original files generated
##    by the mothur program during execution in the WorkDirectory (PWD)
##    - <prefix>.shared
##    - <prefix>.cons.taxonomy
##    - <prefix>.rep.fasta
##
## Dependencies: Mothur package, Silva and RDP datasets
##  Mothur installation: This script assumes the latest version of mothur is
##    installed via the conda package manager in the environment envMothur.
##    If you are not following the same technique, please adjust the mothur
##    installation accordingly. This script requires that the mothur command
##    is available in environment PATH on execution.
##  Silva database: The location for the most up to date silva database for
##    the Dickson lab is in the Turbo drive at
##    /nfs/turbo/umms-rodickso/Databases/MothurDB/SilvaV132
##  RDP database: The location for the most up to date RDP database for the
##    Dickson lab is in the Turbo Drive at
##    /nfs/turbo/umms-rodickso/Databases/MothurDB/RDP16
##
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

## Listing compute nodes allocated to the job
if [[ $SLURM_JOB_NODELIST ]] ; then
	echo "Running on"
	scontrol show hostnames $SLURM_JOB_NODELIST
fi

## Adding Conda init to subshell
if [ -f "$CONDA_PREFIX/etc/profile.d/conda.sh" ]; then
	. "$CONDA_PREFIX/etc/profile.d/conda.sh"
fi

## Job Commands ##

### Activate Mothur and check the package version
conda activate envMothur #deactivate this line if you have Mothur installed differently.
which mothur
mothur -v|perl -pe 's/^.*Linux/Linux/;'

### Set procs/threads for the following commands to use
if [ -z $PROC ] ; then
	PROC=10
fi
### Set fastq directory for the following commands to use
if [ -z $FQDIR ]; then
	# FQDIR='.'
	FQDIR='../MiSeq_M02127_2019_R409-Woods_July2019'
else
	FQDIR=${FQDIR%/}
fi

### Extract files and prepare them for use
parallel -j $PROC 'pigz -dc {} >{/.}' ::: $FQDIR/*.fastq.gz #will extract .gz files in parallel to PWD
ls -1 *.fastq|perl -ne 'chomp;$n=$_;$n=~s/-//g;print `mv $_ $n\n`;' #rename all fastq to remove - from name

### Execute the full Mothur workload
/usr/bin/time -f "\nCommand stats:\nProc:\tElapsed Time = %E,\tPerc CPU = %P,\nMem:\tAvg Total Mem = %KKB,\tPeak Mem = %MKB,\nExit Status: %x" \
mothur "#\
set.current(processors=$PROC);\
make.file(inputdir=., type=fastq, prefix=stability);\
make.contigs(file=stability.files);\
summary.seqs(fasta=stability.trim.contigs.fasta);\
screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxambig=0, maxlength=275, maxhomop=8);\
summary.seqs(fasta=stability.trim.contigs.good.fasta);\
unique.seqs(fasta=stability.trim.contigs.good.fasta);\
count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups);\
summary.seqs(count=stability.trim.contigs.good.count_table);\
align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=/nfs/turbo/umms-rodickso/Databases/MothurDB/SilvaV132/silva.nr_v132.regionV4.align);\
summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table);\
screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary=stability.trim.contigs.good.unique.summary, start=1968, end=11550);\
summary.seqs(fasta=stability.trim.contigs.good.unique.good.align, count=stability.trim.contigs.good.good.count_table);\
filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.);\
unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table);\
pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=2);\
chimera.vsearch(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t);\
remove.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos);\
summary.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table);\
classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=/nfs/turbo/umms-rodickso/Databases/MothurDB/RDP16/trainset16_022016.rdp.fasta, taxonomy=/nfs/turbo/umms-rodickso/Databases/MothurDB/RDP16/trainset16_022016.rdp.tax, cutoff=80);\
remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);\
summary.tax(taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.taxonomy, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table);\
dist.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.03);\
cluster(column=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table);\
make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03);\
classify.otu(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.pick.taxonomy, label=0.03);\
get.oturep(column=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list);\
"|perl -pe 's/^.*Linux/Linux/;'

### Copy important files with more useful names
SHARED="$SLURM_JOB_NAME.shared"
cp stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared $SHARED
TAXONOMY="$SLURM_JOB_NAME.cons.taxonomy"
cp stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy $TAXONOMY
REP="$SLURM_JOB_NAME.rep.fasta"
cp stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.fasta $REP
MOTHURLOG=$(ls -1 mothur.*.logfile)
if [ -f "$SHARED" && -f "$TAXONOMY" && -f "$REP" ]; then
	echo "Files after mothur processing are:\n- $SHARED \n- $TAXONOMY \n- $REP"
	echo "Job logs, after this job finishes can be found in \n- $SLURM_JOB_NAME-$SLURM_JOB_ID.out \n- $SLURM_JOB_NAME-$SLURM_JOB_ID.err"
	echo "Mothur log for this run can be found in $MOTHURLOG"
else
	echo "Somthing has gone wrong with the mothur execution. Please check \n- $SLURM_JOB_NAME-$SLURM_JOB_ID.out \n- $SLURM_JOB_NAME-$SLURM_JOB_ID.err \n- $MOTHURLOG"
fi

## End Job Commands ##
