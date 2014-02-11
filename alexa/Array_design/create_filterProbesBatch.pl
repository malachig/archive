#!/usr/bin/perl -w

#Written by Malachi Griffith
#Take a master file of probe info and divide into into smaller blocks and places them in the specified folder:
#Then create a command for the perl script filterProbes.pl specifying this file, neccessary paramters and output files desired

#The resulting batch job can be submitted to a cluster of computers or simply run on a single one as: 'bash batch_file.sh'

#***NOTE the master probe file used here must be ordered according to probe and probeset ids!!!  i.e. so that all the probes of a probeset are grouped sequentially

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

#ALEXA libraries
#When a script is initiated, use the full path of the script location at execution to add the perl module libraries to @INC
#This should allow this scripts to work regardless of the current working directory or the script location (where it was unpacked).
#The /utilities directory must remain in the same directory as this script but the entire code directory can be moved around
BEGIN {
  my $script_dir = &File::Basename::dirname($0);
  push (@INC, $script_dir);
}
use utilities::utility qw(:all);

my $script_dir = &File::Basename::dirname($0);
my $probe_filter_bin = "$script_dir"."/filterProbes.pl";

#Initialize command line options
my $master_probe_file = '';
my $probesets_bin_size = '';
my $working_dir = '';
my $target_length = '';
my $target_tm = '';
my $tm_range = '';
my $complexity_region_limit = '';
my $simfold_score_limit = '';
my $pairfold_score_limit = '';
my $mrna_hit_percent_limit = '';
my $est_hit_percent_limit = '';
my $enst_hit_percent_limit = '';
my $probe_hit_percent_limit = '';
my $genomic_hit_percent_limit = '';
my $withinGene_hit_percent_limit = '';
my $probeset_size = '';
my $verbose = '';
my $print_probes = '';

GetOptions ('master_probe_file=s'=>\$master_probe_file, 'probesets_bin_size=i'=>\$probesets_bin_size, 'working_dir=s'=>\$working_dir,
	    'target_length=i'=>\$target_length, 'target_tm=f'=>\$target_tm, 'tm_range=f'=>\$tm_range, 'complexity_region_limit=i'=>\$complexity_region_limit,
	    'simfold_score_limit=f'=>\$simfold_score_limit, 'pairfold_score_limit=f'=>\$pairfold_score_limit,
	    'mrna_hit_percent_limit=f'=>\$mrna_hit_percent_limit, 'est_hit_percent_limit=f'=>\$est_hit_percent_limit,
	    'enst_hit_percent_limit=f'=>\$enst_hit_percent_limit, 'probe_hit_percent_limit=f'=>\$probe_hit_percent_limit,
	    'genomic_hit_percent_limit=f'=>\$genomic_hit_percent_limit, 'withinGene_hit_percent_limit=f'=>\$withinGene_hit_percent_limit,
	    'probeset_size=i'=>\$probeset_size, 'verbose=s'=>\$verbose, 'print_probes=s'=>\$print_probes);

#Provide instruction to the user
print GREEN, "\n\nSpecify the MASTER input probe file using: --master_probe_file", RESET;
print GREEN, "\nSpecify the full path to a working directory for temp files using: --working_dir", RESET;
print GREEN, "\nSpecify the number of probesets to place in each temp probe file using: --probesets_bin_size", RESET;
print GREEN, "\n\tIf you have memory issues running this script, try reducing this bin size", RESET;
print GREEN, "\n\tFor negative control probes (which have very large probesets) you may have to use a bin size of 100 or smaller", RESET;
print GREEN, "\nThis script creates a batch file for filtering a large probe file by processing pieces of it", RESET;
print GREEN, "\nYou must specify the following parameters for filtering the probes:", RESET;
print GREEN, "\n\tTarget length using: --target_length (say 36 bp)", RESET;
print GREEN, "\n\tTarget Tm using: --target_tm (say 67.0 degrees)", RESET;
print GREEN, "\n\tAcceptable range of Tm around this target using: --tm_range (say 3.0 degrees)", RESET;
print GREEN, "\n\tThe minimum sized region of low complexity that will cause a probe to fail using: --complexity_region_limit (say 6 bp)", RESET;
print GREEN, "\n\tThe cutoff SimFold score (within-probe/internal folding) that will cause a probe to fail using: --simfold_score_limit (say -7.0)", RESET;
print GREEN, "\n\tThe cutoff PairFold score (self-self/between probe folding) that will cause a probe to fail using: --pairfold_score_limit (say -19.5)", RESET;
print GREEN, "\n\tThe mRNA specifity hit percent length cutoff using: --mrna_hit_percent_limit (say 50.0%)", RESET;
print GREEN, "\n\tThe EST specifity hit percent length cutoff using: --est_hit_percent_limit (say 50.0%)", RESET;
print GREEN, "\n\tThe Ensembl specifity hit percent length cutoff using: --enst_hit_percent_limit (say 50.0%)", RESET;
print GREEN, "\n\tThe Probe specifity hit percent length cutoff using: --probe_hit_percent_limit  (say 75.0%)", RESET;
print GREEN, "\n\tThe Genomic specificity hit percent length cutoff using: --genomic_hit_percent_limit (say 50.0%)", RESET;
print GREEN, "\n\tThe Within gene specificity hit percent length cutoff using:  --withinGene_hit_percent_limit (say 50.0%)", RESET;
print GREEN, "\n\tThe desired probeset size (number of probes per single exon, junction, etc.) using: --probeset_size (say 2)", RESET;
print GREEN, "\n\tFor verbose output, use: verbose=yes\n\n", RESET;
print GREEN, "\n\nIf you just want to generate an updated batch_file, use: --print_probes=no", RESET;
print GREEN, "\n\nExample create_filterProbesBatch.pl  --master_probe_file=exonJunctionProbes.txt  --probesets_bin_size=75000  --working_dir=/home/user/alexa/ALEXA_version/filter   --target_length=36  --target_tm=67.0  --tm_range=3.0  --complexity_region_limit=6  --simfold_score_limit=-7.0  --pairfold_score_limit=-19.5  --mrna_hit_percent_limit=65.0  --est_hit_percent_limit=65.0  --enst_hit_percent_limit=65.0  --probe_hit_percent_limit=65.0  --genomic_hit_percent_limit=65.0  --withinGene_hit_percent_limit=65.0  --probeset_size=3  --verbose=no  --print_probes=yes\n\n", RESET;

unless ($master_probe_file && $probesets_bin_size && $working_dir && $target_length && $target_tm && $tm_range && $complexity_region_limit && $simfold_score_limit && $pairfold_score_limit && $mrna_hit_percent_limit && $est_hit_percent_limit && $enst_hit_percent_limit && $probe_hit_percent_limit && $genomic_hit_percent_limit && $withinGene_hit_percent_limit && $probeset_size && $verbose && $print_probes){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Create temp directories in the working directory: 
#  - $temp_probe_dir ('probe'), $temp_filter_dir (filtered) and $temp_log_dir (log)
unless ($working_dir =~ /.*\/$/){
  $working_dir = "$working_dir"."/";
}
#Make sure the specified directory is valid
unless (-e $working_dir && -d $working_dir){
  print RED, "\nSpecified working directory does not appear to be valid: $working_dir\n\n", RESET;
  exit();
}

#Get the base name of the probe supplied probe file
my $probe_file_name = &File::Basename::basename($master_probe_file);
my $base_name;
if ($probe_file_name =~ /^(.*)\./){
  $base_name = $1;
}else{
  $base_name = $probe_file_name;
}

#Create temp directories
my $sub_dir = &createNewDir('-path'=>$working_dir, '-new_dir_name'=>$base_name);
my $temp_probe_dir = &createNewDir('-path'=>$sub_dir, '-new_dir_name'=>"probe");
my $temp_filter_dir = &createNewDir('-path'=>$sub_dir, '-new_dir_name'=>"filtered");
my $temp_log_dir = &createNewDir('-path'=>$sub_dir, '-new_dir_name'=>"log");

my $batch_file = "$sub_dir"."filter_"."$base_name".".sh";
print BLUE, "\nWriting batch file commands to: $batch_file\n\n", RESET;

#Open the output batch file
open (BATCHFILE, ">$batch_file") || die "\nCould not open batch file: $batch_file\n\n";

my $probe_count = 0;
my $total_probe_count = 0;
my $probe_files_printed = 0;
my $block_num = 1;
my ($probe_file,$out_file,$temp_log_file,$header_line);

open (PROBEFILE, "$master_probe_file") || die "\nCould not open probe file: $master_probe_file";

#Find out the number of probesets in the input file
my %probesets;
my $num_probes = 0;
while (<PROBEFILE>){
  my @line = split("\t", $_);

  if ($line[1] =~ /^\d+/){
    $probesets{$line[1]}{t}='';
    $num_probes++;
  }
}
my $probesets_count = keys %probesets;
close PROBEFILE;

print BLUE, "\nFound a total of $probesets_count probesets in the master probe file\n\n", RESET;

open (PROBEFILE, "$master_probe_file") || die "\nCould not open probe file: $master_probe_file";

#open the first output file
$probe_file = "$temp_probe_dir"."probes"."$block_num";
$out_file = "$temp_filter_dir"."filtered_probes$block_num";
$temp_log_file = "$temp_log_dir"."filter_LOG$block_num";

#Print 1st batch file line
print BATCHFILE "$probe_filter_bin --infile $probe_file  --target_length $target_length  --target_tm $target_tm  --tm_range $tm_range  --complexity_region_limit $complexity_region_limit  --simfold_score_limit $simfold_score_limit  --pairfold_score_limit $pairfold_score_limit  --mrna_hit_percent_limit $mrna_hit_percent_limit  --est_hit_percent_limit $est_hit_percent_limit  --enst_hit_percent_limit $enst_hit_percent_limit  --probe_hit_percent_limit $probe_hit_percent_limit  --genomic_hit_percent_limit $genomic_hit_percent_limit  --withinGene_hit_percent_limit $withinGene_hit_percent_limit  --probeset_size $probeset_size  --outfile $out_file --logfile $temp_log_file  --verbose $verbose\n";
$probe_files_printed++;

unless ($print_probes eq "no"){
  open (TEMP_FILE, ">$probe_file") || die "\nCould not open temp probe file: $probe_file\n\n";
}

my $probesets_processed_count = 0;
my $probes_printed = 0;
my $first_line = 1;

my %probesets_processed;

while (<PROBEFILE>){
  my $line_record = $_;

  #Grab the header line
  if ($first_line == 1){
    $header_line = $_;
    unless ($print_probes eq "no"){
      print TEMP_FILE "$header_line";
    }
   $first_line = 0;
    next();
  }

  #Watch for end of file
  if ($total_probe_count == $num_probes){
    next();
  }

  $total_probe_count++;

  #Get the probeset ID of this line
  my @line = split("\t", $line_record);
  my $probeset_id;
  if ($line[1] =~ /^\d+/){
    $probeset_id = $line[1];
  }else{
    print RED, "\nProbeset ID missing or not understood\n\n", RESET;
    exit();
  }

  #If this probeset has already been started or the target number of probesets in the bin has not been reached, print it out to the current file
  #otherwise start the next probe file
  if ($probesets_processed{$probeset_id} || $probesets_processed_count < $probesets_bin_size){
    $probesets_processed{$probeset_id}{t} = '';
    $probesets_processed_count = keys %probesets_processed;
    unless ($print_probes eq "no"){
      print TEMP_FILE "$line_record";
      $probes_printed++;
    }
  }else{
    #Once a full block of probesets is reached, and as long as we are not in the middle of a probeset, start on the next file
    %probesets_processed = ();
    $probesets_processed{$probeset_id}{t} = '';
    $probesets_processed_count = keys %probesets_processed;

    #Close the current output file
    unless ($print_probes eq "no"){
       close TEMP_FILE;
    }
    $block_num++;
    $probe_file = "$temp_probe_dir"."probes"."$block_num";
    $out_file = "$temp_filter_dir"."filtered_probes$block_num";
    $temp_log_file = "$temp_log_dir"."filter_LOG$block_num";

    #Open the new output file
    unless ($print_probes eq "no"){
      open (TEMP_FILE, ">$probe_file") || die "\nCould not open temp probe file: $probe_file\n\n";

      #Print the header line to this file
      print TEMP_FILE "$header_line";

      #print the current probe to this new file
      print TEMP_FILE "$line_record";
      $probes_printed++;
    }

    print BATCHFILE "$probe_filter_bin --infile $probe_file  --target_length $target_length  --target_tm $target_tm  --tm_range $tm_range  --complexity_region_limit $complexity_region_limit  --simfold_score_limit $simfold_score_limit  --pairfold_score_limit $pairfold_score_limit  --mrna_hit_percent_limit $mrna_hit_percent_limit  --est_hit_percent_limit $est_hit_percent_limit  --enst_hit_percent_limit $enst_hit_percent_limit  --probe_hit_percent_limit $probe_hit_percent_limit  --genomic_hit_percent_limit $genomic_hit_percent_limit  --withinGene_hit_percent_limit $withinGene_hit_percent_limit  --probeset_size $probeset_size  --outfile $out_file --logfile $temp_log_file  --verbose $verbose\n";

    $probe_files_printed++;
  }
}

unless ($print_probes eq "no"){
   close TEMP_FILE;
}

#Create commands to concatenate filtered probe and log files and clean up directories

#Create master log file and clean-up log dir
my $master_log_file = "$working_dir"."filter_"."$base_name"."_LOG.txt";
print BATCHFILE "\ncat $temp_log_dir* > $master_log_file";
print BATCHFILE "\nrm -fr $temp_log_dir";

#Clean up temp probe dir
print BATCHFILE "\nrm -fr $temp_probe_dir";

#Join probe files.  First need to find the total number of probes in them
my $cat_path = "$temp_filter_dir"."filtered*";
my $probe_join_bin = "$script_dir"."/joinProbeFiles.pl";
my $joined_probe_file = "$working_dir"."$base_name"."_filtered.txt";
print BATCHFILE "\nFILE_COUNT=`ls $cat_path | wc -l`";
print BATCHFILE "\nTOTAL_LINES=`cat $temp_filter_dir* | wc -l`";
print BATCHFILE "\nPROBE_COUNT=\$((\$TOTAL_LINES - \$FILE_COUNT))";
print BATCHFILE "\n$probe_join_bin  --source_dir=$temp_filter_dir  --probes_expected=`echo \$PROBE_COUNT`  --output_file=$joined_probe_file";

#Clean up filtered probe dir
print BATCHFILE "\nrm -fr $temp_filter_dir";

close PROBEFILE;
close BATCHFILE;

print BLUE, "\nFound a total of $num_probes probes in the input file", RESET;
print BLUE, "\nPrinted a total of $probes_printed probes", RESET;
print BLUE, "\nThese were printed in blocks of $probesets_bin_size probesets resulting in a total of $probe_files_printed files\n\n", RESET;

exit();


