#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2010 Malachi Griffith

#Process data lane-by-lane (one unit of data is: flowcell -> lane -> R1 or R2

#Steps.
#1.) Get an export file path and other parameters from the command line
#2.) Create a temp working dir for output and move to this dir
#3.) Create a fastq file from the export file
#4.) Run hmmSplicer using this fastq file
#    - Source genome db is here (hg18)
#    /projects/malachig/sequence_databases/hg18_genome
#5.) Rename the output files according to: flowcell -> lane -> R1 or R2
#    - Each run produces 3 output files: log.txt, junction.final.bed, junction.nonCanonical.bed
#6.) Move these output files to a results dir
#    /projects/alexa2/hmmSplicer
#    Place each file in a subdirectory according to its library ID
#7.) Delete all unneeded temporary files

#Example hmmSplicer command:
#/home/malachig/tools/python/python_install/bin/python /home/malachig/tools/hmmSplicer/hmmSplicer-0.9.0/runHMM.py -p 4 -o /projects/malachig/sequence_databases/hg18_genome/hmmSplicerResults/ -i 30CGLAAXX_Lane7_R1.fastq -g ./all_human.fa -p 4 -t /nobackup/tmp_malachig/ -q p -r True -k 100000

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;
use Benchmark;

my $read_file = '';
my $read_file_type = '';
my $library_id = '';
my $flowcell_lane_read = '';
my $analysis_dir = '';
my $ref_genome_db = '';
my $python_path = '';
my $hmm_splicer_path = '';
my $maq_dir = '';
my $cpus = '';
my $quality_flag = '';
my $test = '';
my $nobackup_dir = '';

GetOptions ('read_file=s'=>\$read_file, 'read_file_type=s'=>\$read_file_type, 'library_id=s'=>\$library_id, 'flowcell_lane_read=s'=>\$flowcell_lane_read, 'analysis_dir=s'=>\$analysis_dir, 
            'ref_genome_db=s'=>\$ref_genome_db, 'python_path=s'=>\$python_path, 'hmm_splicer_path=s'=>\$hmm_splicer_path, 'maq_dir=s'=>\$maq_dir,
            'cpus=i'=>\$cpus, 'quality_flag=s'=>\$quality_flag, 'test=i'=>\$test, 'nobackup_dir=s'=>\$nobackup_dir);

if ($read_file && $read_file_type && $library_id && $flowcell_lane_read && $analysis_dir && $ref_genome_db && $python_path && $hmm_splicer_path && $cpus && $quality_flag){
  #Parameters ok
}else{
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  print GREEN, "\n\nExample: hmmSplicerRun.pl  --read_file=/projects/rgoyaprj2/projects/breast_cancer/illumina/SA063/sequence/RNA/illumina/42HWLAAXX_1_1_export.txt.bz2  --read_file_type=export  --library_id=SA063  --flowcell_lane_read=42HWLAAXX_1_1  --analysis_dir=/projects/alexa2/hmmSplicer/  --ref_genome_db=/projects/malachig/sequence_databases/hg18_genome/all_human.fa  --python_path=/home/malachig/tools/python/python_install/bin/python  --hmm_splicer_path=/home/malachig/tools/hmmSplicer/hmmSplicer-0.9.0/runHMM.py  --maq_dir=/home/malachig/tools/MAQ/maq-0.7.1/  --cpus=2  --quality_flag=[p|s|t]  [ --test=1  --nobackup_dir=/nobackup/ ]\n\n", RESET;
  exit();
}

#Make sure read file type is allowed
unless ($read_file_type =~ /^export$|^fastq$|^qseq$/i){
  print RED, "\n\n--read_file_type must be one of: export | fastq | qseq | ...", RESET;
  exit();
}
chomp($quality_flag);
unless($quality_flag =~ /^p$|^s$|^t$/){
  print RED, "\n\n--quality_flag must be one of [p | s | t]: 'p' = Phred format.  's' = old solexa format (pre1.3).  't' = new solexa format (pipeline 1.3+)", RESET;
  exit();
}

#Check input files and dirs
unless (-e $analysis_dir && -d $analysis_dir){
  print RED, "\n\nAnalysis dir ($analysis_dir) is not valid\n\n", RESET;
  exit();
}
my $results_dir = "$analysis_dir"."results/";
my $log_dir = "$analysis_dir"."logs/";
my $logfile = "$log_dir"."$flowcell_lane_read".".log.txt";

unless (-e $read_file){
  print RED, "\n\nExport file ($read_file) is not valid\n\n", RESET;
  exit();
}
unless (-e $maq_dir && -d $maq_dir){
  print RED, "\n\nMaq script dir ($maq_dir) is not valid\n\n", RESET;
  exit();
}
unless (-e $ref_genome_db){
  print RED, "\n\nReference genome ($ref_genome_db) is not valid\n\n", RESET;
  exit();
}
unless (-e $python_path){
  print RED, "\n\nPython path ($python_path) is not valid\n\n", RESET;
  exit();
}
unless (-e $hmm_splicer_path){
  print RED, "\n\nHmmSplicer path ($hmm_splicer_path) is not valid\n\n", RESET;
  exit();
}
unless ($analysis_dir =~ /.*\/$/){
  $analysis_dir .= "/";
}
unless ($maq_dir =~ /.*\/$/){
  $maq_dir .= "/";
}
my $maq_script_dir = "$maq_dir"."scripts/";

unless(-e $results_dir && -d $results_dir){
  mkdir($results_dir);
}
unless(-e $log_dir && -d $log_dir){
  mkdir($log_dir);
}

#Global time variable
my $t1 = new Benchmark;

#0.) Define file paths and check to see if the job was already done...
my $lib_dir = "$results_dir"."$library_id/";
my $lib_fc_dir = "$lib_dir"."$flowcell_lane_read/";
my $tmp_dir = "$lib_fc_dir"."temp/";
my $jfb_file = "$lib_fc_dir"."junction.final.bed";
my $jnb_file = "$lib_fc_dir"."junction.nonCanonical.bed";
my $hmm_log_file = "$lib_fc_dir"."log.txt";

#If the user specified to use 'nobackup' space, create a temp dir there instead
if ($nobackup_dir){
  unless ($nobackup_dir =~ /.*\/$/){
    $nobackup_dir .= "/";
  }  
  unless (-e $nobackup_dir && -d $nobackup_dir){
    print RED, "\n\nNobackup space specified does not appear valid: $nobackup_dir\n\n", RESET;
    exit();
  }
  $tmp_dir = "$nobackup_dir"."$flowcell_lane_read/";
}

if (-e $jfb_file && -e $jnb_file && -e $hmm_log_file){
  print RED, "\n\nResults files already present, aborting...\n\n", RESET;
  exit();
}


open (LOG, ">$logfile") || die "\n\nCould not open log file: $logfile\n\n";
print LOG "Begin processing of a single lane of data by hmmSplicerRun.pl";
print LOG "\n\nUser specified the following options:\nexport_file = $read_file\nlibrary_id = $library_id\nflowcell_lane_read = $flowcell_lane_read\nanalysis_dir = $analysis_dir\nref_genome_db = $ref_genome_db\npython_path = $python_path\nhmm_splicer_path = $hmm_splicer_path\nmaq_dir = $maq_dir\ncpus = $cpus";


#1.) Create a temp working dir for output and move to this dir
unless (-e $lib_dir && -d $lib_dir){
  mkdir($lib_dir);
}
if (-e $lib_fc_dir && -d $lib_fc_dir){
  system("rm -fr $lib_fc_dir");
  mkdir($lib_fc_dir);
}else{
  mkdir($lib_fc_dir);
}

#Create a temp sub dir to be used for the hmmSplicer run
if (-e $tmp_dir){
  system("rm -fr $tmp_dir");
  mkdir($tmp_dir);
}else{
  mkdir($tmp_dir);
}

#Deal with compression
my $fastq_file = "$tmp_dir"."$flowcell_lane_read".".fq";
my $read_fh_cat;
if ($read_file =~ /\.txt$/){
  $read_fh_cat = "cat $read_file";
}elsif($read_file =~ /\.gz$/){
  $read_fh_cat = "zcat $read_file";
}elsif($read_file =~ /\.bz2$/){
  $read_fh_cat = "bzcat $read_file";
}else{
  print RED, "\n\nExtension on read file not recognized (.txt | .txt.gz | .txt.bz): $read_file\n\n", RESET;
  print LOG "\n\nExtension on read file not recognized (.txt | .txt.gz | .txt.bz): $read_file\n\n";
  exit();
}

print BLUE, "\n\nDealing with input read file format...", RESET;
if ($read_file_type =~ /^export$/i){
  #EXPORT FILES
  #2.) Get an export file path and other parameters from the command line
  #    - Note export file could be .bz, .gz, or .txt

  #3.) Create a fastq file from the export file
  #Use MAQ conversion script for this conversion
  my $export_to_fq_cmd = "$read_fh_cat | $maq_script_dir"."fq_all2std.pl export2std > $fastq_file";

  #TESTING
  if ($test){
    $export_to_fq_cmd = "$read_fh_cat | $maq_script_dir"."fq_all2std.pl export2std | head -n 400000 > $fastq_file";
  }
  print BLUE, "\n\nConverting export file to standard fastq format...", RESET;
  print BLUE, "\n$export_to_fq_cmd", RESET;
  print LOG "\n\nConverting export file to standard fastq format...";
  print LOG "\n$export_to_fq_cmd";
  system($export_to_fq_cmd);
}elsif ($read_file_type =~ /^fastq$/i){
  #FASTQ FILES

  if ($quality_flag =~ /^s$/){
    #If a quality flag of 's' was specified, use MAQ to convert qualities and update the quality flag to 'p'
    my $convert_fq_to_fq_cmd = "$read_fh_cat | $maq_dir"."maq sol2sanger - $fastq_file";
    if ($test){
      $convert_fq_to_fq_cmd = "$read_fh_cat | head -n 400000 | $maq_dir"."maq sol2sanger - $fastq_file";
    }
    $quality_flag = "p";
    print BLUE, "\n\nConverting solexa fastq file to standard fastq format...", RESET;
    print BLUE, "\n$convert_fq_to_fq_cmd", RESET;
    print LOG "\n\nConverting solexa fastq file to standard fastq format...";
    print LOG "\n$convert_fq_to_fq_cmd";
    system($convert_fq_to_fq_cmd);
  }else{
    #Otherwise just make a copy of the fastq
    my $fq_to_fq_cmd = "$read_fh_cat > $fastq_file";
    if ($test){
      $fq_to_fq_cmd = "$read_fh_cat | head -n 400000 > $fastq_file";
    }
    print BLUE, "\n\nMaking simple copy of file already in fastq format...", RESET;
    print BLUE, "\n$fq_to_fq_cmd", RESET;
    print LOG "\n\nMaking simple copy of file already in fastq format...";
    print LOG "\n$fq_to_fq_cmd";
    system($fq_to_fq_cmd);
  }

}elsif ($read_file_type =~ /^qseq$/i){
  #QSEQ FILES  
  #Use MAQ conversion script for this conversion
  my $qseq_to_fq_cmd = "$read_fh_cat | $maq_script_dir"."fq_all2std.pl qseq2std > $fastq_file";

  #TESTING
  if ($test){
    $qseq_to_fq_cmd = "$read_fh_cat | $maq_script_dir"."fq_all2std.pl qseq2std | head -n 400000 > $fastq_file";
  }
  print BLUE, "\n\nConverting qseq file to standard fastq format...", RESET;
  print BLUE, "\n$qseq_to_fq_cmd", RESET;
  print LOG "\n\nConverting qseq file to standard fastq format...";
  print LOG "\n$qseq_to_fq_cmd";
  system($qseq_to_fq_cmd);

}else{  
  print RED, "\n\nFile type not understood!\n\n", RESET;
  exit();
}



#4.) Run hmmSplicer using this fastq file
#    - Source genome db is here (hg18)
#    /projects/malachig/sequence_databases/hg18_genome
my $stdout = "$tmp_dir"."hmm_stdout.txt";
my $stderr = "$tmp_dir"."hmm_stderr.txt";
my $hmm_splicer_cmd = "$python_path -E $hmm_splicer_path -p $cpus -o $lib_fc_dir -i $fastq_file -g $ref_genome_db -t $tmp_dir -q $quality_flag -r True -k 100000 1>$stdout 2>$stderr";

if ($test){
  $hmm_splicer_cmd = "$python_path -E $hmm_splicer_path -p $cpus -o $lib_fc_dir -i $fastq_file -g $ref_genome_db -t $tmp_dir -q $quality_flag -r True -k 100000";
}

print BLUE, "\n\nRunning hmmSplicer...", RESET;
print BLUE, "\n$hmm_splicer_cmd\n", RESET;
print LOG "\n\nRunning hmmSplicer...";
print LOG "\n$hmm_splicer_cmd\n";
system($hmm_splicer_cmd);

#5.) Rename the output files according to: flowcell -> lane -> R1 or R2
#    - Each run produces 3 output files: log.txt, junction.final.bed, junction.nonCanonical.bed


#6.) Move these output files to a results dir
#    /projects/alexa2/hmmSplicer
#    Place each file in a subdirectory according to its library ID


#7.) Delete all unneeded temporary files

#8.) Perform some basic quality checks to determine if the run was successful
#    - Did the hmmSplicer log file get to the last step?
#    - Are the results files present
#    - How many canonical junctions were found
my $pass = 1;
if (-e $jfb_file && -e $jnb_file && -e $hmm_log_file){
  print BLUE, "\n\nFound all three results files", RESET;
  print LOG "\n\nFound all three results files";
}else{
  print RED, "\n\nDid not find all three results files", RESET;
  print LOG "\n\nDid not find all three results files";
  $pass = 0;
}

#Get time since last check
my $t2 = new Benchmark;
my $td = timediff($t2, $t1);
my $time_string = timestr($td);
my $seconds = "N";
my $minutes = "N";
my $hours = "N";
if ($time_string =~ /(\d+)\s+wallclock/){
  $seconds = $1;
  $minutes = sprintf("%.2f", ($seconds/60));
  $hours = sprintf("%.2f", ($minutes/60));
}

print BLUE, "\n\nProcessing time = $seconds seconds | $minutes minutes | $hours hours", RESET;
print LOG "\n\nProcessing time = $seconds seconds | $minutes minutes | $hours hours";

if ($pass){
  print "\n\nSCRIPT COMPLETE\n\n";
  print LOG "\n\nSCRIPT COMPLETE\n\n";
  my $rm_cmd = "rm -fr $tmp_dir";
  system($rm_cmd);
}else{
  print "\n\nJOB FAILED\n\n";
  print LOG "\n\nJOB FAILED\n\n";
}
close(LOG);

print "\n\n";

exit();


