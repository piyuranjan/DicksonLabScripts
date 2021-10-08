# DicksonLabScripts

Scripts for everyday work in the Dickson lab.

## Setup

### By cloning

Repository can be cloned and put on the login profile if it is frequently used. Users are recommended to use [GitHub Desktop](https://desktop.github.com/) to clone and keep the repository updated for its ease. However, this can also be achieved with replacing the `<Path/To/Folder>` with the location of the repository in the following commands.

Clone and update.

```
$ git clone https://github.com/piyuranjan/DicksonLabScripts.git
$ cd <Path/To/Folder>/DicksonLabScripts/ && git pull
```

To make these scripts easily available from any location (great when using often), they can be added to the system `PATH` variable. To do this, add the following code to a login profile file like the `~/.profile` or `~/.bashrc`.

```
## DicksonLabScripts GitHub repo addition
export PATH=<Path/To/Folder>/DicksonLabScripts:$PATH
```

### By pulling scripts on-demand

A single script can be pulled for use from this repository by replacing `<scriptName>` with name of one of the scripts available in this repository with the following command.

```
$ wget https://raw.githubusercontent.com/piyuranjan/DicksonLabScripts/main/<scriptName>
```

After pulling on-demand, scripts will need to be addressed with their path for execution.

---

### [`mothurGreatLakes.sh`](./mothurGreatLakes.sh)

A SLURM shell script for running `mothur` on Great Lakes.

This shell script includes SLURM directives to run `mothur` 16S microbiome analysis jobs on the UM ARC-TS Great Lakes cluster with minimal user intervention. This script is specifically tailored for running by the members of the Dickson Lab but can easily be tweaked for alternative uses. This script follows Pat Schloss' Mothur [MiSeq_SOP](https://www.mothur.org/wiki/MiSeq_SOP) very closely.

Go to the [Wiki page](https://github.com/piyuranjan/DicksonLabScripts/wiki/mothurGreatLakes.sh) for full instructions and configuration. Quick steps for execution are below.

- Upload the gzipped fastq for the 16S data in a directory on the scratch filesystem on the Great Lakes cluster.
- No need to rename fastq files. It will be automatically taken care of during execution.
- Copy the script in an empty work directory (preferably in a separate location than the 16S fastq) on the scratch filesystem on the Great Lakes cluster.
- Execute with sbatch as: `sbatch --export=ALL,FQDIR="<path/to/fastq/>" -J <jobName> mothurGreatLakes.sh`

[`mothurGreatLakes_ghuff.sh`](./mothurGreatLakes_ghuff.sh) is the Huffnagle lab version of the same script. Directions to use are going to be the same over on the [Wiki page](https://github.com/piyuranjan/DicksonLabScripts/wiki/mothurGreatLakes.sh) except for modification of the filename.

---

### [`summarizeFastq.pl`](./summarizeFastq.pl)

Summarize fastq file[s] using `seqtk`.

Features:

- Will process gzipped fastq (*.fastq.gz).
- Multiple fastq can be submitted as positional arguments.
- Can be parallelized on multiple files using GNU `parallel`.
- Examples are included that demonstrate use cases.

---

