# DicksonLabScripts
Scripts for everyday work in the Dickson lab

---

[`mothurGreatLakes.sh`](./mothurGreatLakes.sh) A SLURM shell script for running mothur on Great Lakes.

A shell script with SLURM directives to run mothur 16S microbiome analysis jobs on the UM ARC-TS Great Lakes cluster with minimal user intervention.

This script is specifically tailored for running by members of the Dickson Lab but can easily be tweaked for alternative uses. This script follows Pat Schloss' Mothur [MiSeq_SOP](https://www.mothur.org/wiki/MiSeq_SOP) very closely.

- Usage and configuration can be found on the [Wiki page](https://github.com/piyuranjan/DicksonLabScripts/wiki/mothurGreatLakes.sh)
- Make a copy of this script in an empty work directory on the scratch filesystem on the Great Lakes and execute with `sbatch`.
- A quick use case: `sbatch --export=ALL,FQDIR="../MiSeq_M02127_2019/" --job-name=ardsLungs mothurGreatLakes.sh`

---

[`summarizeFastq.pl`](./summarizeFastq.pl) Summarizes a fastq file uisng seqtk.

- Usage with no arguments and full help with `-h`.
- Will process gzipped fastq (.fastq.gz) the same way.
- Multiple files can be submitted as positional arguments.
- Code is meant for serial execution but can be parallelized on multiple files using GNU Parallel.
- Examples include all use cases.

---
