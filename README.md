# DicksonLabScripts
Scripts for everyday work in the Dickson lab

---

[`mothurGreatLakes.sh`](./mothurGreatLakes.sh) A SLURM shell script for running mothur on Great Lakes.

A shell script with SLURM directives to run mothur microbial analysis jobs on the UM ARC-TS Great Lakes cluster with minimal user intervention.
This script is specifically tailored for running by Dickson Lab members but can easily be tweaked for alternative uses. This script follows Pat Schloss' Mothur [MiSeq_SOP](https://www.mothur.org/wiki/MiSeq_SOP) very closely.
- Usage and configuration can be found on the [Wiki page](https://github.com/piyuranjan/DicksonLabScripts/wiki/mothurGreatLakes.sh)
- Make a copy of this script in an empty work directory on the scratch filesystem on the Great Lakes and execute with `sbatch`.
- A quick use case: `sbatch --export=ALL,FQDIR="../MiSeq_M02127_2019/" --job-name=ardsLungs mothurGreatLakes.sh`

---
