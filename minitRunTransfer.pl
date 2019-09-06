#!/usr/bin/env perl
use strict;

die "## This script transfers the passed fastq files from the MinIT device via rsync. ##
#### Makes a run folder *-fastqPass inside the current directory.
#### Performs repetitive transfer with a sleep cycle unless killed (use Ctrl+C to kill).
#### Please put this script in your PATH variable to execute it easily.

Usage: $0 [ExperimentID] [SampleID-if-used]\n" unless defined $ARGV[0];

##Assign variables
chomp(@ARGV);
my $expId=$ARGV[0];
my $sampleId=$expId;
$sampleId=$ARGV[1] if(defined $ARGV[1]);
my $sleepSec=30; #Sleeping seconds, can be changed here for desired effect.

##prepare local folder
my $runFolder=$expId;
$runFolder.=$sampleId if(defined $ARGV[1]);
$runFolder.="-fastqPass";
mkdir $runFolder unless(-d $runFolder);

while(1)
	{
	print "\n### Looking for new files in the MinIT for run: $expId\n";
	print `rsync -vrthP minit\@10.42.0.1::MinitData/$expId/$sampleId/*/fastq_pass/*.fastq $runFolder`;
	
	print "### Sleeping for $sleepSec seconds now\n";
	sleep($sleepSec);
	}
