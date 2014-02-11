#!/usr/bin/perl -w

#Written by Malachi Griffith
#Generate simple statistics from probes files containing probe characteristics
#For example: average Tm, simFold score, pairFold score, etc.

#This script examines a directory of probe files, and generates statistics for several columns of these files using R
#1.) Create a working directory as specified by the user
#2.) Parse the probe directory and create temp tiles
#    - Print temp files containing only the columns to be summarized for each probe file
#    - Also print a seperate file containing all experimental probes together (all but negative control)
#    - Also make note of the location of a probe_coverage file specified by the user
#3.) Test the output directory specified by the user
#4.) Create an R script to process each group of data (each temp file created)
#    - create commands to summarize data columns and generate figures in the output directory
#5.) Run the R script and redirect the output to a log file of the R command output

#NOTE: If you run this on a 64-bit machine make sure you specify a 64-bit version of R
#NOTE2: If you run this on a remote host via an SSH connection you need to import your display settings when you connect
#       - Use ssh -X hostname
#       - Otherwise jpeg files will not be created
#       - If you have trouble generating jpeg files you can change the R commands to use a postscript device instead

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

my $r_bin = '';               #Full path to R binary
my $probe_dir = '';           #Full path to directory containing probe files
my $temp_dir = '';            #Full path working directory, sub-directory will be created within this
my $results_dir = '';         #Full path to desired location of output files
my $logfile = '';

GetOptions ('r_bin=s'=>\$r_bin,'probe_dir=s'=>\$probe_dir, 'temp_dir=s'=>\$temp_dir,
	    'results_dir=s'=>\$results_dir, 'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nThis generates statistics for probe files in a specified directory", RESET;
print GREEN, "\nUsage:", RESET;
print GREEN, "\n\tSpecify the full path to a 'R' binary using: --r_bin", RESET;
print GREEN, "\n\tSpecify the full path to a directory containing probe files using: --probe_dir", RESET;
print GREEN, "\n\tSpecify the full path to a working directory using: --temp_dir\n", RESET;
print GREEN, "\n\tSpecify the full path to a directory for output of files: --results_dir", RESET;
print GREEN, "\n\tSpecify a logfile for output using: --logfile", RESET;
print GREEN, "\n\nExample: getProbeStats.pl  --r_bin=/usr/bin/R  --probe_dir=/home/user/alexa/ALEXA_version/unfiltered_probes/  --temp_dir=/home/user/alexa/ALEXA_version/stats/  --results_dir=/home/user/alexa/ALEXA_version/stats/unfiltered_probes/  --logfile=/home/user/alexa/ALEXA_version/logs/getProbeStats_unfiltered_LOG.txt\n\n", RESET;

unless ($r_bin && $probe_dir && $temp_dir && $results_dir && $logfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";
print LOG "\nUser Specified the following options:\nr_bin = $r_bin\nprobe_dir = $probe_dir\ntemp_dir = $temp_dir\nresults_dir = $results_dir\nlogfile = $logfile\n\n";

#Make sure the R binary exists
unless (-e $r_bin){
  print RED, "\nR binary: $r_bin does not appear valid!\n\n", RESET;
  exit();
}

#1.) Create a working directory as specified by the user
my $working_dir = &createNewDir('-path'=>$temp_dir, '-new_dir_name'=>"probe_stats_temp");

#2.) Test the output directory specified by the user
unless ($results_dir =~ /.*\/$/){
  $results_dir = "$results_dir"."/";
}
#First make sure the specified base path exists and is a directory
unless (-e $results_dir && -d $results_dir){
  print RED, "\nSpecified results directory: $results_dir does not appear valid!\n\n", RESET;
  close (LOG);
  exit();
}

#Make a figures directory within the specified results directory
my $figures_dir = &createNewDir('-path'=>$results_dir, '-new_dir_name'=>"figures");


#3.) Parse the probe directory and create temp tiles
#    - Print temp files containing only the columns to be summarized for each probe file
#    - Also print a seperate file containing all experimental probes together (all but negative control)
#    - Also make note of the location of a probe_coverage file specified by the user

#Create list of all desired columns to be summarized
my %desired_columns;
$desired_columns{ProbeSet_ID}{name} = "Probes Per Probeset";
$desired_columns{ProbeSet_ID}{display_order} = 1;
$desired_columns{Probe_length}{name} = "Probe Length";
$desired_columns{Probe_length}{display_order} = 2;
$desired_columns{Probe_Tm}{name} = "Probe Tm";
$desired_columns{Probe_Tm}{display_order} = 3;
$desired_columns{masked_bases}{name} = "Masked Bases"; #All 'NA' for NC probes
$desired_columns{masked_bases}{display_order} = 4;
$desired_columns{SimFold_score}{name} = "Hairpin Folding Energy";
$desired_columns{SimFold_score}{display_order} = 5;
$desired_columns{PairFold_score}{name} = "Self-self Folding Energy";
$desired_columns{PairFold_score}{display_order} = 6;
$desired_columns{MdustBases_cutoff_11}{name} = "Low Complexity Bases";
$desired_columns{MdustBases_cutoff_11}{display_order} = 7;
$desired_columns{genomic_TargetHits}{name} = "BLAST Hits to Target Region"; #All 0 for NC probes
$desired_columns{genomic_TargetHits}{display_order} = 8;
$desired_columns{'genomic_Non-TargetHits'}{name} = "Genomic Non-specific Hits";
$desired_columns{'genomic_Non-TargetHits'}{display_order} = 9;
$desired_columns{'genomic_Non-TargetHits_over50'}{name} = "Genomic Non-specific Hits >50% probe length";
$desired_columns{'genomic_Non-TargetHits_over50'}{display_order} = 10;
$desired_columns{'genomic_Non-TargetHits_over75'}{name} = "Genomic Non-specific Hits >75% probe length";
$desired_columns{'genomic_Non-TargetHits_over75'}{display_order} = 11;
$desired_columns{'genomic_largestNon-TargetHitLength'}{name} = "Size of Genomic Non-specific Hit";
$desired_columns{'genomic_largestNon-TargetHitLength'}{display_order} = 12;
$desired_columns{'genomic_largestNon-TargetPercentIdentity'}{name} = "Percent Identity of Non-specific Hits";
$desired_columns{'genomic_largestNon-TargetPercentIdentity'}{display_order} = 13;
$desired_columns{CycleCount}{name} = "NimbleGen Cycle Count";
$desired_columns{CycleCount}{display_order} = 14;

my @desired_columns = keys %desired_columns;

my $temp_files_ref = &parseProbeFiles('-probe_dir'=>$probe_dir, '-working_dir'=>$working_dir);

#4.) Create an R script to process each group of data (each temp file created)
#    - create commands to summarize data columns and generate figures in the output directory
#    - Save a log of the R commands run
my $r_script_path = "$working_dir"."probe_stats.R";
&createRScript('-r_script'=>$r_script_path, '-working_dir'=>$working_dir, '-results_dir'=>$results_dir, '-figures_dir'=>$figures_dir,
	       '-files'=>$temp_files_ref);

#5.) Run the R script and redirect the output to a log file of the R command output
my $r_output = "$results_dir"."probe_stats_R_output_LOG.txt";
my $command = "$r_bin --no-save < $r_script_path > $r_output";

print BLUE, "\n\nExecuting command: $command", RESET;
print BLUE, "\nRedirecting R output to: $r_output", RESET;
print BLUE, "\nBecause R loads files so slowly, this may take a while\n\n", RESET;
print LOG "\n\nExecuting command: $command";
print LOG "\nRedirecting R output to: $r_output";
print LOG "\nBecause R loads files so slowly, this may take a while\n\n";

system ("$command");

#Clean up the temp directory unless the user says otherwise
print YELLOW, "\nRemove temp working dir: $working_dir (y/n)? ", RESET;
my $answer = <>;
chomp($answer);

if ($answer eq "y"){
  my $command = "rm -r $working_dir";
  system ($command);
}else{
  print YELLOW, "\nLeaving it alone\n", RESET;
}
print "\n\n";
close (LOG);

exit();



#######################################################################################
#Parse the input probe files and generate temp output files                           #
#######################################################################################
sub parseProbeFiles{
  my %args = @_;
  my $probe_dir = $args{'-probe_dir'};
  my $working_dir = $args{'-working_dir'};

  #1.) Get files from this directory
  my %files;

  unless ($probe_dir =~ /.*\/$/){
    $probe_dir = "$probe_dir"."/";
  }
  #First make sure the specified base path exists and is a directory
  unless (-e $probe_dir && -d $probe_dir){
    print RED, "\nSpecified directory: $probe_dir does not appear valid!\n\n", RESET;
    close (LOG);
    exit();
  }

  #Get all the input files in the repeat masked result directory
  print BLUE, "\nSearching $probe_dir for files\n", RESET;
  print LOG "\nSearching $probe_dir for files\n";

  my $file_count = 0;
  opendir(DIRHANDLE, "$probe_dir") || die "\nCannot open directory: $probe_dir\n\n";
  my @test_files = readdir(DIRHANDLE);

  #Create a list of all observed columns
  my %master_column_list;
  my $master_column_count = 0;

 FILE:foreach my $test_file (@test_files){
    my $file_path = "$probe_dir"."$test_file";

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      print BLUE, "\n\t$file_path  is a directory - skipping", RESET;
      print LOG "\n\t$file_path  is a directory - skipping";
      next();
    }

    #Only process the first line of each file for now
    my %columns;
    open (PROBE_FILE, "$file_path") || die "\nCould not open probe file: $file_path\n\n";
    while (<PROBE_FILE>){
      my $line_record = $_;
      chomp ($line_record);

      my @line = split ("\t", $line_record);

      #Watch for the header line which is assumed to contain neccessary columns
      my $column_count = 0;
      foreach my $column (@line){
	$columns{$column}{column_pos} = $column_count;
	$column_count++;

	unless ($master_column_list{$column}){
	  $master_column_count++;
	  $master_column_list{$column}{column_pos} = $master_column_count;
	}
      }

      #If the file has the neccessary columns, add it to the list of files to be processed
      if ($columns{'Probe_Count'} && $columns{'Probe_Tm'}){
	$file_count++;
	$files{$file_count}{file_path} = $file_path;
	$files{$file_count}{file_name} = $test_file;
	$files{$file_count}{source_dir} = $probe_dir;
	$files{$file_count}{columns} = \%columns;

	my $new_file_name = "$test_file"."_simple";
	$files{$file_count}{new_file_name} = $new_file_name;
	$files{$file_count}{new_file_path} = "$working_dir"."$new_file_name";

      }else{
	print YELLOW, "\nFile: $file_path does not appear to be a valid probe file (expected column missing) - skipping\n\n", RESET;
	print LOG "\nFile: $file_path does not appear to be a valid probe file (expected column missing) - skipping\n\n";
      }
      close (PROBE_FILE);
      next FILE;  #skip to next file after processing only the first line
    }
  }

  #2.) Summarize the files found - list old and new new names
  print BLUE, "\n\nFile summary:", RESET;
  print LOG "\n\nFile summary:";
  foreach my $file_count (sort {$a <=> $b} keys %files){
    print BLUE, "\n\tFile: $file_count\n\t\tINFILE: $files{$file_count}{file_path}\n\t\tTEMP_FILE: $files{$file_count}{new_file_path}", RESET;
    print LOG "\n\tFile: $file_count\n\t\tINFILE: $files{$file_count}{file_path}\n\t\tTEMP_FILE: $files{$file_count}{new_file_path}";
  }
  print BLUE, "\n\n", RESET;
  print LOG "\n\n";


  #3.) Go through each file, get probe data and write it to a temp output file which includes data for all probe types (except NC)
  my $all_probes_file_name = "allProbes.txt_simple";
  my $all_probes_file_path = "$working_dir"."$all_probes_file_name";
  my $combined_nmax = 0;


  open (ALL_FILE, ">$all_probes_file_path") || die "\nCould not open all probes file: $all_probes_file_path\n\n";
  foreach my $column (sort {$master_column_list{$a}->{column_pos} <=> $master_column_list{$b}->{column_pos}} keys %master_column_list){
    unless ($desired_columns{$column}){
      next();
    }
    print ALL_FILE "$column\t";
  }
  print ALL_FILE "\n";

  my $global_max_id = 0;
  foreach my $file_count (sort {$a <=> $b} keys %files){
    my $nmax = 0;

    print BLUE, "\nBegin processing the input probe file: $files{$file_count}{file_path}\n\n", RESET;
    print LOG "\nBegin processing the input probe file: $files{$file_count}{file_path}\n\n";

    my %columns = %{$files{$file_count}{columns}};

    #Use linux 'cut' commands to dump desired columns to create a new simplified temp file for each input file
    my $column_numbers = '';
    my $column_count = 0;
    foreach my $desired_column (@desired_columns){
      if ($columns{$desired_column}){
	$column_count++;
	my $current_col = $columns{$desired_column}{column_pos} + 1;
	$column_numbers = "$column_numbers"."$current_col".",";
      }
    }
    $files{$file_count}{column_count} = $column_count;

    my $command = "/usr/bin/cut -f $column_numbers $files{$file_count}{file_path} > $files{$file_count}{new_file_path}";
    system ("$command");
    print BLUE, "\nCreate new file with command: $command\n\n", RESET;

    #If the current probe is a negative control probe, skip it
    #Go through all of the possible columns in the master column list
    #If it is not a desired column, skip it
    #Find the corresponding value for this column in the current file
    #If the value is not defined in the current file, set the value to na
    #print the value

    open (PROBE_FILE, "$files{$file_count}{file_path}") || die "\nCould not open probe file: $files{$file_count}{file_path}\n\n";
    my $header = 1;
    while (<PROBE_FILE>){

      #skip the first line
      if ($header == 1){
	$files{$file_count}{max_id} = 0;
	$header = 0;
	next();
      }
      $nmax++;

      my $line_record = $_;
      chomp ($line_record);
      my @line = split ("\t", $line_record);

      if ($line[0] > $files{$file_count}{max_id}){
	$files{$file_count}{max_id} = $line[0];
      }
      if ($line[0] > $global_max_id){
	$global_max_id = $line[0];
      }

      if ($columns{'Probe_Type'}){
	my $probe_type = $line[$columns{'Probe_Type'}{column_pos}];

	if ($probe_type eq "Control-Negative"){
	  next();
	}
      }

      $combined_nmax++;

      foreach my $column (sort {$master_column_list{$a}->{column_pos} <=> $master_column_list{$b}->{column_pos}} keys %master_column_list){
	unless ($desired_columns{$column}){
	  next();
	}
	#Is this value in the current file? If not set to 'na'
	my $data_value = 'na';

	if ($columns{$column}){
	  $data_value = $line[$columns{$column}{column_pos}];
	}
	print ALL_FILE "$data_value\t";
      }
      print ALL_FILE "\n";
    }
    $files{$file_count}{nmax} = $nmax;
    close (PROBE_FILE);
  }
  close (ALL_FILE);

  #Count the number of columns in the combined file
  my $combined_column_count = 0;
  foreach my $column (sort {$master_column_list{$a}->{column_pos} <=> $master_column_list{$b}->{column_pos}} keys %master_column_list){
    unless ($desired_columns{$column}){
      next();
    }
    $combined_column_count++;
  }
  #Add the all probes file to the list of files to be processed
  $file_count++;
  $files{$file_count}{new_file_name} = $all_probes_file_name;
  $files{$file_count}{new_file_path} = $all_probes_file_path;
  $files{$file_count}{columns} = \%master_column_list;
  $files{$file_count}{column_count} = $combined_column_count;
  $files{$file_count}{nmax} = $combined_nmax;
  $files{$file_count}{max_id} = $global_max_id + 1;

  return(\%files);
}


#####################################################################################################
#Create R script to generate stats                                                                  #
#####################################################################################################
sub createRScript{
  my %args = @_;
  my $r_script_path = $args{'-r_script'};
  my $working_dir = $args{'-working_dir'};
  my $results_dir = $args{'-results_dir'};
  my $figures_dir = $args{'-figures_dir'};
  my $files_ref = $args{'-files'};

  my $figure_count = 0;

  print BLUE, "\nCreating R script: $r_script_path\n\n", RESET;
  print LOG "\nCreating R script: $r_script_path\n\n";

  #open R script
  open (R_OUT, ">$r_script_path") || die "\nCould not open output R script: $r_script_path\n\n";

  print R_OUT "#Begin R script for summarizing probe stats for files located in: $working_dir\n\n";
  #Generic tasks that can be completed for each input file
  foreach my $file_count (sort {$files_ref->{$a}->{max_id} <=> $files_ref->{$b}->{max_id}} keys %{$files_ref}){

    #Open input file and create data frame
    my $file_name = $files_ref->{$file_count}->{new_file_name};
    print R_OUT "\n\n#$file_count.) Process input file:  $file_name\n";
    print R_OUT "datadir = \"$working_dir\"\n";
    print R_OUT "setwd(datadir)\n";

    #Use 'read.table' method
    my $nmax = ($files_ref->{$file_count}->{nmax})+1;
    print R_OUT "\nprobe_data = read.table(\"$file_name\", header=T, quote=\"\", sep=\"\\t\", comment.char=\"\", na.strings=\"na\", nrows=$nmax, colClasses=\"numeric\")\n";

    #Determine the number of columns in this file and use this to generate a list of quotes for the scan function
    my $column_count = $files_ref->{$file_count}->{column_count};
    my $quotes = '';
    for (my $i = 0; $i < $column_count-1; $i++){
      $quotes = "$quotes"."\"\",";
    }
    $quotes = "$quotes"."\"\"";

    #Get the columns found in this file
    my $columns_ref = $files_ref->{$file_count}->{columns};

    #Change to an output directory
    print R_OUT "\n#Set output directory for results to: $figures_dir\n";
    print R_OUT "datadir = \"$figures_dir\"\n";
    print R_OUT "setwd(datadir)\n";

    #Go through each desired column, if it is defined in the current file, process it as follows:
    #Create number datapoints, number non-zero datapoints, summaries, boxplots and histograms for each data column
    foreach my $desired_column (sort {$desired_columns{$a}->{display_order} <=> $desired_columns{$b}->{display_order}} keys %desired_columns){

      unless ($columns_ref->{$desired_column}){
	next();
      }

      my $data_name = $desired_columns{$desired_column}{name};

      #Replace '-' with '.' in column names (because R does this)
      $desired_column =~ tr/\-/\./;
      print R_OUT "\n##Data column: $desired_column ($data_name)\n";

      print R_OUT "length(probe_data[,\"$desired_column\"])\n";

      print R_OUT "length(which(probe_data[,\"$desired_column\"] > 0))\n";

      print R_OUT "summary(probe_data[,\"$desired_column\"])\n\n";

      my $dataset;
      if ($file_name =~ /^(.*)\./){
	$dataset = $1;
      }else{
	$dataset = $file_name;
      }

      #Special cases for certain columns of data
      #Gene_ID
      if ($desired_column eq "ProbeSet_ID"){

	$figure_count++;
	my $figure_number = &getFigureNumber('-current_figure_number'=>$figure_count);
	my $figure_name = "f"."$figure_number"."_ProbesPerProbeset"."_$dataset"."_hist".".jpg";

	print R_OUT "probes_per_probeset = table(probe_data[,\"$desired_column\"])\n";
	print R_OUT "jpeg(file = \"$figure_name\", quality=100, res=600)\n";
	print R_OUT "hist(probes_per_probeset, main=\"Probes per ProbeSet (Target Region)\", xlab=\"ProbeSet Size\", ylab=\"ProbeSet Count ($dataset)\", col=\"blue\")\n";
	print R_OUT "dev.off()\n";

	$figure_count++;
	$figure_number = &getFigureNumber('-current_figure_number'=>$figure_count);
	$figure_name = "f"."$figure_number"."_ProbesPerProbeset"."_$dataset"."_bplot".".jpg";

	print R_OUT "jpeg(file = \"$figure_name\", quality=100, res=600)\n";
	print R_OUT "boxplot(probes_per_probeset, main=\"Probes per ProbeSet (Target Region)\", xlab = \"ProbeSet Size\", ylab =\"ProbeSet Count ($dataset)\",  col=\"red\")\n";
	print R_OUT "dev.off()\n";

	next();
      }

      #Masked_bases - deal with cases where all probes are 'NA' for this column
      if ($desired_column eq "masked_bases"){

	$figure_count++;
	my $figure_number = &getFigureNumber('-current_figure_number'=>$figure_count);
	my $figure_name = "f"."$figure_number"."_ProbesPerGene"."_$dataset"."_bplot".".jpg";

	print R_OUT "na_values = length(which(is.na(probe_data[,\"$desired_column\"])))\n";
	print R_OUT "if (na_values == (length(probe_data[,\"$desired_column\"]))){print(\"Column: $desired_column is all na's\")}else{\n";
	print R_OUT "jpeg(file = \"$figure_name\", quality=100, res=600)\n";
	print R_OUT "boxplot(probe_data[,\"$desired_column\"], main=\"Distribution of $data_name\", xlab = \"($dataset)\", ylab =\"$data_name\", col=\"red\")\n";
	print R_OUT "dev.off()\n";

	$figure_count++;
	$figure_number = &getFigureNumber('-current_figure_number'=>$figure_count);
	$figure_name = "f"."$figure_number"."_$desired_column"."_$dataset"."_hist".".jpg";
	print R_OUT "jpeg(file = \"$figure_name\", quality=100, res=600)\n";
	print R_OUT "hist(probe_data[,\"$desired_column\"], main=\"Distribution of $data_name\", xlab =\"$data_name\", ylab=\"Probe Count ($dataset)\", col=\"blue\")\n";
	print R_OUT "dev.off()\n";
	print R_OUT "}";

	next();
      }


      #Generic figures for most columns
      $figure_count++;
      my $figure_number = &getFigureNumber('-current_figure_number'=>$figure_count);

      my $figure_name = "f"."$figure_number"."_$desired_column"."_$dataset"."_bplot".".jpg";
      print R_OUT "jpeg(file = \"$figure_name\", quality=100, res=600)\n";

      print R_OUT "boxplot(probe_data[,\"$desired_column\"], main=\"Distribution of $data_name\", xlab = \"($dataset)\", ylab =\"$data_name\", col=\"red\")\n";

      print R_OUT "dev.off()\n";

      $figure_count++;
      $figure_number = &getFigureNumber('-current_figure_number'=>$figure_count);

      $figure_name = "f"."$figure_number"."_$desired_column"."_$dataset"."_hist".".jpg";
      print R_OUT "jpeg(file = \"$figure_name\", quality=100, res=600)\n";

      print R_OUT "hist(probe_data[,\"$desired_column\"], main=\"Distribution of $data_name\", xlab =\"$data_name\", ylab=\"Probe Count ($dataset)\", col=\"blue\")\n";

      print R_OUT "dev.off()\n";

    }
  }

  close (R_OUT);
  return();
}


########################################################################################################
#Provide formatted figure number
########################################################################################################
sub getFigureNumber{
  my %args = @_;
  my $current_figure = $args{'-current_figure_number'};

  my $figure_number;
  if ($current_figure =~ /^\d{1}$/){
    $figure_number = "00"."$current_figure";
  }elsif ($current_figure =~ /^\d{2}$/){
    $figure_number = "0"."$current_figure";
  }elsif ($current_figure =~ /^\d{3}$/){
    $figure_number = "$current_figure";
  }else{
    print RED, "\nCurrent figure number: $current_figure not understood or out of range\n\n";
    exit();
  }

  return ($figure_number);
}


