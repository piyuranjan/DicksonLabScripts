#!/usr/bin/env perl

use Getopt::Long;
use POSIX qw(strftime);
use Time::HiRes qw(time);
use Time::Seconds;
use strict;
use warnings;
use v5.10;

use diagnostics; #only while in development

#########################
######### Usage #########
#########################

sub Info
	{
	say "
################################################################################
## Summarize fastq[.gz] file(s) with basic statistics obtained
## via seqtk.
##
## Input:
##  - .fastq or .fastq.gz filename(s) as positional argument(s) or/and
##  - pattern(s) that match all .fastq[.gz] files in a path.
## 
## Output:
##  A tab delimited file (or on STDOUT) with the fields: #Dataset,
##  #Reads, #Bases, #MinLen, #MaxLen, #AvgLen, #AvgQ, #AvgA, #AvgA,
##  #AvgC, #AvgG, #AvgT, #AvgN
##
## Dependencies: seqtk >= 1.2-r94
##
## Author: piyuranjan\@gmail.com
## Source:
##  https://github.com/piyuranjan/PerlScripts/blob/master/summarizeFastq.pl
################################################################################
	";
	}
sub Usage
	{
	say "
Usage:\n$0 [options] <*.fastq[.gz]|fastqFile1..n> >[outputFile]

Arguments: at least one of the following is required
 *.fastq[.gz]   Feel free to use a pattern for .fastq or .fastq.gz files.
 fastqFile1..n  Specify files as positional arguments.

Options:
 -c|colNames    Only show column names and exit 0.
 -f|force       Force overwrite outFile if it exists. Depends on -o|outFile.
 -n|noHeader    Do not include header line in the output.
 -o|outFile     [str:optional] Send output to a file. STDOUT otherwise.
 -h|help        Show more help and exit 0.
 -v|verbose     Use multiple times to increase verbosity.
 -q|quiet       Silent execution except for significant/fatal errors (ERR).
                 No warnings (WARN), statistics (STAT) and information (INFO)
                 when -q engaged. Make sure you know what you're doing.
 --debug        Set max verbosity for debugging.
	";
	}
sub Examples
	{
	say "
Examples:

### Newbie: Find all .fastq.gz in PWD and summarize them.
$0 *.fastq.gz

### Techie: Summarize a list of .fastq[.gz] file names stored
### in another file fastqFiles.txt. Store results in summary.txt.
cat fastqFiles.txt|xargs $0 >summary.txt

### Expert: Find all *.fastq.gz in Path/to/Fastq/ and all *.fastq in PWD
### and summarize them together in the file summary.txt.
$0 -v -v -o summary.txt Path/To/Fastq/*.fastq.gz *.fastq

### Wizard: Pair it with GNU Parallel to summarize all *.fastq.gz files in Path/to/Fastq/ in parallel on 6 threads
parallel -j 6 $0 -v {} ::: Path/to/Fastq/*.fastq.gz | sort | uniq >summary.txt

### Ghost: Parallel wizardry with absolute silence (except for errors)
parallel -j 6 $0 -q {} ::: Path/to/Fastq/*.fastq.gz | sort | uniq >summary.txt
	";
	}
sub StatusInfo
	{
	say "
Status message categories:
 - ERR: Fatal errors!
 - INFO: Information about status of the program or time step.
 - STAT: Statistics about time or other parameters.
 - WARN: Non-fatal warnings!

Exit codes:
 0  Successful
 1  Unsuccessful or need arguments
 2  Problem with a dependency
 3  Overwrite permission
	";
	}


#########################
######### Main ##########
#########################

my $startTime=time();

### Preset column names for output
my @colNames=qw(Dataset Reads Bases MinLen MaxLen AvgLen AvgQ AvgA AvgC AvgG AvgT AvgN);

## Option reading and configuration
my($colNames,$force,$noHeader,$outFile,$threads,$help,$verbose,$quiet,$debug);
$colNames=0; #set default show only colNames to false
$noHeader=0; #set default no-header to false
$threads=1; #set default threads to 1
$verbose=1; #set default verbose to 1
unless(GetOptions(
				'c|colNames' =>\$colNames,
				'f|force' => \$force,
				'n|noHeader' => \$noHeader,
				'o|outFile=s' => \$outFile,
				't|threads=i' => \$threads,
				'h|help' => \$help,
				'v|verbose+' => \$verbose,
				'q|quiet' => \$quiet,
				'debug' => \$debug))
	{Usage;exit 1;} #quit with error code
if($help) #quit with full help
	{Info;Usage;Examples;StatusInfo;exit 0;}
if($colNames) #show only the header line and exit
	{PrintHeader(\@colNames);exit 0;}
### Set necessary arguments
unless(defined $ARGV[0])
	{warn TimeStamp($verbose)."ERR: Need arguments.\n";Usage;exit 1;}
### Set option behaviors
$verbose=0 if($quiet);
if($debug) #override verbosity to highest value for debugging code
	{$verbose=100;warn TimeStamp($verbose)."INFO: Debug mode enabled. Setting verbosity to max.\n";}
if($force && !$outFile) #warn that option force is dependent on option outFile
	{warn TimeStamp($verbose)."WARN: Unnecessary use of force. Can't see option -o|outfile specified.\n";}

## Check dependencies
my $seqtkVersion="1.2-r94";
unless(CheckSeqtk($verbose,$seqtkVersion))
	{warn TimeStamp($verbose)."ERR: Exiting\n";exit 2;}

## Output redirection if outFile is given
my $OUT;
if(defined $outFile)
	{
	if(-e $outFile) #outFile overwrite check
		{
		warn TimeStamp($verbose)."WARN: Output file exists: $outFile\n" if($verbose);
		if($force) #check if forced to overwrite
			{warn TimeStamp($verbose)."WARN: Overwriting with -f|force\n" if($verbose);}
		else #stop if force option is not used
			{warn TimeStamp($verbose)."ERR: Need elevation to overwrite. Use -f|force\n";exit 3;}
		}
	open($OUT,">",$outFile) or die $!;
	warn TimeStamp($verbose)."INFO: Output will be written to: $outFile\n" if($verbose>1);
	select $OUT; #redirect STDOUT to $OUT for writing in file
	}
### For completeness, remember to revert redirection when its use has finished.
### Not necessary if you are not planning to write anything else to STDOUT or another file

## Parse and prepare input files
my @fqFiles;
foreach my $arg(@ARGV)
	{
	my @files=glob("$arg");
	push(@fqFiles,@files);
	}
warn TimeStamp($verbose)."STAT: Found ".scalar(@fqFiles)." file(s) to summarize\n" if($verbose);
if($verbose>1) #list all files found
	{warn TimeStamp($verbose)."INFO: List of files:\n";warn "$_\n" foreach @fqFiles;}

## Execute seqtk on all fastq files and display output
warn TimeStamp($verbose)."INFO: The header line has ".scalar(@colNames)." columns.\n" if($verbose>3);
PrintHeader(\@colNames) unless($noHeader); #print the header line
my $processedFiles=0;
foreach my $fqFile(@fqFiles)
	{
	my $recordTime=time();
	unless(-r $fqFile) #skip file if not found/readable
		{warn TimeStamp($verbose)."WARN: Skipping fastq; not found/readable: $fqFile\n" if($verbose);next;}
	
	### Prepare and run seqtk
	my $seqtkCommand="seqtk fqchk $fqFile|head -4";
	warn TimeStamp($verbose)."INFO: Running seqtk command: `$seqtkCommand`\n" if($verbose>3);
	my @seqtkOut=`$seqtkCommand`;
	if($verbose>3)
		{warn TimeStamp($verbose)."INFO: Finished seqtk command:\n";warn $_ foreach @seqtkOut;}
	unless($seqtkOut[3]) #warn user if the file format has not been proper fastq
		{
		warn TimeStamp($verbose)."WARN: Sequences are not in fastq format; will skip processing: $fqFile\n" if($verbose);
		warn TimeStamp($verbose)."WARN: Consider checking the fastq file or run with -v -v (extra verbose) to investigate.\n" if($verbose);
		next;
		}
	
	### Parse seqtk output to create summary elements
	my %summary;
	$summary{$_}='' foreach(@colNames); #initialize all elements of summary
	$summary{Dataset}=$fqFile;
	($summary{MinLen},$summary{MaxLen},$summary{AvgLen})=$seqtkOut[0]=~/^min_len:\h+(\d*?);\h+max_len:\h+(\d*?);\h+avg_len:\h+(\d*?\.?\d*?);/;
	($summary{Bases},$summary{AvgA},$summary{AvgC},$summary{AvgG},$summary{AvgT},$summary{AvgN},$summary{AvgQ})=$seqtkOut[2]=~/^ALL\h+(\d*?)\h+(\d*?\.?\d*?)\h+(\d*?\.?\d*?)\h+(\d*?\.?\d*?)\h+(\d*?\.?\d*?)\h+(\d*?\.?\d*?)\h+(\d*?\.?\d*?)\h+/;
	($summary{Reads})=$seqtkOut[3]=~/^1\h+(\d*?)\h+/;
	### If some of these elements are undefined or zero, means a problem.
	### Not adding exceptions/error handlers here because seqtk seems to handle them well.
	
	### Print and count successful summary
	print "$summary{$_}\t" foreach(@colNames);say "";
	$processedFiles++;
	
	### Print time taken to process this fastq
	my $diffTime=Time::Seconds->new(int((time()-$recordTime)+0.5));
	my $spentTime=$diffTime->pretty;
	warn TimeStamp($verbose)."STAT: Spent $spentTime for $fqFile\n" if($verbose>1);
	}

### Revert output redirection when finished
if(defined $outFile)
	{close($OUT);select STDOUT;}

### Remove file (if) created, if errors have happend to a point that no fastqFile has been processed.
###  Don't remove if verbose is at level debug.
if($outFile && $processedFiles<1 && $verbose<100)
	{
	warn TimeStamp($verbose)."WARN: Seems no records to write, will delete $outFile.\n" if($verbose);
	unlink $outFile;
	}

### Print finishing status with total time taken
my $diffTime=Time::Seconds->new(int((time()-$startTime)+0.5));
my $spentTime=$diffTime->pretty;
warn TimeStamp($verbose)."STAT: Finished summarizing $processedFiles file(s) in $spentTime\n" if($verbose);


#########################
###### Subroutines ######
#########################

sub PrintHeader #Prints \t spaced header column names to STDOUT (or a fileHandle) with a '#' apended before names
	{
	### Order of Input: \@columnNames (a reference to array); optional $fileHandle
	### Order of Output: none
	my @colNames=@{$_[0]};
	my $FH=$_[1];
	
	if(defined $FH) #print to filehandle if provided
		{print $FH "#$_\t" foreach @colNames;say "";}
	else #default to STDOUT
		{print "#$_\t" foreach @colNames;say "";}
	}

sub CheckSeqtk #Checks if seqtk is available on PATH, returns good/bad status
	{
	### Order of Input: verbosityLevel-optional; versionOfSeqtkToCheck-optional
	### Order of Output: SeqtkStatus:0-notFound,1-good,-1-lower
	warn TimeStamp($verbose)."INFO: Starting function CheckSeqtk\n" if($verbose>2);
	
	## Read input and set defaults
	my $verbose=$_[0];
	$verbose//=1;
	my $needVersion=$_[1];
	$needVersion//="1.2-r94";
	
	## Check seqtk and extract version
	warn TimeStamp($verbose)."INFO: Running seqtk command: `seqtk 2>&1`\n" if($verbose>3);
	my @seqtkOut=`seqtk 2>&1`;
	if($verbose>3)
		{warn TimeStamp($verbose)."INFO: Finished seqtk command:\n";warn $_ foreach @seqtkOut;}
	### If seqtk is not found in PATH, suggest easiest installation and return with 0
	unless(@seqtkOut)
		{
		warn TimeStamp($verbose)."ERR: seqtk not found in PATH\n";
		warn TimeStamp($verbose)."INFO: Please install seqtk via `conda install seqtk` or `sudo apt install seqtk`\n";
		warn TimeStamp($verbose)."INFO: Or please check out https://github.com/lh3/seqtk \n";
		return(0);
		}
	my $version=(split(/\h+/,$seqtkOut[2]))[-1]; chomp($version);
	warn TimeStamp($verbose)."INFO: Found seqtk version: $version!\n" if($verbose>2);
	
	## Compare seqtk version
	if($version ge $needVersion)
		{
		warn TimeStamp($verbose)."INFO: seqtk $version found is same/higher than $needVersion needed. Good to go!\n" if($verbose>2);
		return(1);
		}
	else
		{
		warn TimeStamp($verbose)."WARN: seqtk $version found is lower than $needVersion.\n" if($verbose);
		warn TimeStamp($verbose)."WARN: Program may not work as intended.\n" if($verbose);
		return(-1);
		}
	}

sub TimeStamp #Provides the current time in various formats
	{
	### Order of Input: verbosityLevel-optional
	### Order of Output: "[formattedCurrentTimeAsString] "
	
	## Read input and set defaults
	my $verbose=$_[0]; #Read in the verbosity level
	$verbose//=1; #Set verbose to 1 if nothing was passed
	
	## Find current time, format according to verbosity level
	my $currentTime=time();
	my $formattedTime=strftime("%Y%m%d %H:%M:%S", localtime $currentTime);
	if($verbose>=3) #add microseconds if verbose 3
		{$formattedTime.=sprintf(".%06d", ($currentTime-int($currentTime))*1000000);}
	elsif($verbose>=2) #add miliseconds if verbose 2
		{$formattedTime.=sprintf(".%03d", ($currentTime-int($currentTime))*1000);}
	$formattedTime="[$formattedTime] ";
	return($formattedTime);
	}
