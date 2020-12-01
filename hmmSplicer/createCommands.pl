#!/usr/bin/perl -w
# This script will merge a series of bed files from hmmSplicer into a single BED file
# Use this script to create commands for hmmSplicers runs 
# The main task here is to find the raw data path, and figure out the location and format of raw read data...
# The input for this job will be a table of LIMS information for all libraries, and a list of library IDs and library names for the project

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;
use Benchmark;

my $library_list = '';
my $lims_data = '';
my $analysis_dir = '';
my $project = '';
my $ref_genome = '';

GetOptions ('library_list=s'=>\$library_list, 'lims_data=s'=>\$lims_data, 'analysis_dir=s'=>\$analysis_dir, 'project=s'=>\$project, 'ref_genome=s'=>\$ref_genome);

if ($library_list && $lims_data && $analysis_dir && $project && $ref_genome){
  #Parameters ok
}else{
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  print GREEN, "\nSpecify a file containing library ID to library name mappings using: --library_list", RESET;
  print GREEN, "\nSpecify a file containing LIMS data for all libraries of all types using: --lims_data", RESET;
  print GREEN, "\n\nExample: createCommands.pl  --library_list=/projects/alexa2/hmmSplicer/SA_TN_Breast/jobs/LibraryList.txt  --lims_data=/archive/solexa1_4/analysis/projects/all_library_status.out  --analysis_dir=/projects/alexa2/hmmSplicer/  --project=SA_TN_Breast  --ref_genome=/projects/malachig/sequence_databases/hg18_genome/all_human.fa\n\n", RESET;
  exit();
}
unless($analysis_dir =~ /\/$/){
  $analysis_dir .= "/";
}
unless(-e $analysis_dir && -d $analysis_dir){
  print RED, "\n\nAnalysis dir ($analysis_dir) does not appear valid\n\n", RESET;
  exit();
}
unless (-e $ref_genome){
  print RED, "\n\nCould not find reference genome file: $ref_genome\n\n", RESET;
  exit();
}


#Load the library list
my %libs = %{&getLibraryList()};
#print Dumper %libs;

#Load a list of all lanes of WTSS data for each library
my %lanes = %{&getLaneData()};

#Get the read1 and read2 files 
#Determine their type (export?, fastq?, merged qseq?, qseq tiles?, merged seq?, seq tiles?)
#Determine their compression status
print BLUE, "\n\nExamining data paths, searching for read files and determining format and compression status", RESET;
my $files_found = 0;
my $broken_dirs = 0;
my @suffix = qw(.txt .txt.gz .txt.bz2);
foreach my $lc (sort {$lanes{$a}{lib_id} cmp $lanes{$b}{lib_id}} keys %lanes){
  my $lib_id = $lanes{$lc}{lib_id};
  my $data_path = $lanes{$lc}{data_path};
  my $flowcell = $lanes{$lc}{flowcell};
  my $lane = $lanes{$lc}{lane};

  #qseq summary file names
  my $r1_concat_qseq_a = "$data_path"."basecalls/qseq_summary/s_"."$lane"."_1_concat_qseq";
  my $r2_concat_qseq_a = "$data_path"."basecalls/qseq_summary/s_"."$lane"."_2_concat_qseq";
  my $r1_concat_qseq_b = "$data_path"."bustard/qseq_summary/s_"."$lane"."_1_concat_qseq";
  my $r2_concat_qseq_b = "$data_path"."bustard/qseq_summary/s_"."$lane"."_2_concat_qseq";

  #export files (e.g. 61F92AAXX_2_1_export.txt.bz2)
  my $r1_export = "$data_path"."$flowcell"."_"."$lane"."_1_export";
  my $r2_export = "$data_path"."$flowcell"."_"."$lane"."_2_export";

  #Export files, shadow filtered (e.g. 436CFAAXX_8_1_export_f1b6.txt.gz)
  my $r1_export_f1b6 = "$data_path"."$flowcell"."_"."$lane"."_1_export_f1b6";
  my $r2_export_f1b6 = "$data_path"."$flowcell"."_"."$lane"."_2_export_f1b6";

  my %test_files;
  $test_files{'1'}{r1_path} = $r1_concat_qseq_a;
  $test_files{'1'}{r2_path} = $r2_concat_qseq_a;
  $test_files{'1'}{type} = "qseq";
  $test_files{'2'}{r1_path} = $r1_concat_qseq_b;
  $test_files{'2'}{r2_path} = $r2_concat_qseq_b;
  $test_files{'2'}{type} = "qseq";
  $test_files{'3'}{r1_path} = $r1_export;
  $test_files{'3'}{r2_path} = $r2_export;
  $test_files{'3'}{type} = "export";
  $test_files{'4'}{r1_path} = $r1_export_f1b6;
  $test_files{'4'}{r2_path} = $r2_export_f1b6;
  $test_files{'4'}{type} = "export";


  print YELLOW, "\n\n\tCHECKING PATH: '$data_path' (lib:$lib_id  flowcell:$flowcell  lane:$lane)", RESET;

  #Now check for each file
  my $r1_file_found = 0;
  my $r2_file_found = 0;
  foreach my $tf (sort {$a <=> $b} keys %test_files){
    my $r1_base_path = $test_files{$tf}{r1_path};
    my $r2_base_path = $test_files{$tf}{r2_path};

    my $type = $test_files{$tf}{type};
    #Add each suffix
    foreach my $suffix (@suffix){
      my $r1_full_path = "$r1_base_path"."$suffix";
      my $r2_full_path = "$r2_base_path"."$suffix";

      #Look for the read1 and read2 files AS A PAIR - it could cause problems to have different file formats for R1 and R2 downstream
      #print CYAN, "\n\t\t$r1_full_path\n\t\t$r2_full_path", RESET;
      if (-e $r1_full_path && -s $r1_full_path && $r1_file_found == 0 && -e $r2_full_path && -s $r2_full_path && $r2_file_found == 0){
        $r1_file_found = 1;
        $r2_file_found = 1;
        print YELLOW, "\n\t\tFound R1 file: $r1_full_path (type: $type $suffix)", RESET;
        print YELLOW, "\n\t\tFound R2 file: $r2_full_path (type: $type $suffix)", RESET;
        $files_found+=2;

        #Store the lane read files for R1 and R2
        $lanes{$lc}{r1_path} = $r1_full_path;
        $lanes{$lc}{r2_path} = $r2_full_path;
        $lanes{$lc}{file_format} = $type;
        $lanes{$lc}{suffix} = $suffix;
      }
    }
  }

  #Make note of files not found
  unless ($r1_file_found == 1 && $r2_file_found){
    print RED, "\n\t\tCould not find R1 and/or R2 file in: '$data_path'", RESET;
    $broken_dirs++;
  }
}

my $lib_count = keys %libs;
my $lane_count = keys %lanes;
my $expected_r12_files = $lane_count*2;

print BLUE, "\n\nSUMMARY STATS:", RESET;
print BLUE, "\nTotal libraries considered = $lib_count", RESET;
print BLUE, "\nTotal passing lanes found = $lane_count", RESET;
print BLUE, "\nTotal files expected = $expected_r12_files", RESET;
print BLUE, "\nTotal files found = $files_found", RESET;
print BLUE, "\nBroken dirs = $broken_dirs", RESET;

#Print output files. (1) Data files.  (2) Broken dirs (3) hmmSplicer commands 
my $data_files_file = "$analysis_dir"."$project/jobs/DataFilesList.txt";
my $broken_dirs_file = "$analysis_dir"."$project/jobs/BrokenDirs.txt";
my $hmmSplicer_jobs_file = "$analysis_dir"."$project/jobs/hmmSplicerJobs.sh";
my $mergeJunctions_file = "$analysis_dir"."$project/jobs/mergeJunctionBed.sh";
my $mergeJunctionsNonCan_file = "$analysis_dir"."$project/jobs/mergeNonCanJunctionBed.sh";
my $annotateJunctionFile = "$analysis_dir"."$project/jobs/annotateObservedJunctions.sh";

open(DATAFILES, ">$data_files_file") || die "\n\nCould not open data files for output: $data_files_file\n\n";
open(BROKENDIRS, ">$broken_dirs_file") || die "\n\nCould not open broken dirs for output: $broken_dirs_file\n\n";
open(HMMSPLICER, ">$hmmSplicer_jobs_file") || die "\n\nCould not open commands for output: $hmmSplicer_jobs_file\n\n";
open(MERGE, ">$mergeJunctions_file") || die "\n\nCould not open commands for output: $mergeJunctions_file\n\n";
open(MERGE2, ">$mergeJunctionsNonCan_file") || die "\n\nCould not open commands for output: $mergeJunctionsNonCan_file\n\n";
open(ANNOTATE, ">$annotateJunctionFile") || die "\n\nCould not open commands for output: $annotateJunctionFile\n\n";

print DATAFILES "LibraryID\tLibraryName\tFlowcell\tLane\tRead\tReadLength\tFileFormat\tSuffix\tPath\n";
foreach my $lc (sort {$lanes{$a}{sort_key} cmp $lanes{$b}{sort_key}} keys %lanes){
  my $lib_id = $lanes{$lc}{lib_id};
  my $lib_name = $lanes{$lc}{lib_name};
  my $data_path = $lanes{$lc}{data_path};
  my $flowcell = $lanes{$lc}{flowcell};
  my $lane = $lanes{$lc}{lane};
  my $read_length = $lanes{$lc}{read_length};
  my $r1_full_path = $lanes{$lc}{r1_path};
  my $r2_full_path = $lanes{$lc}{r2_path};
  my $file_format = $lanes{$lc}{file_format};
  my $suffix = $lanes{$lc}{suffix};
  my $flowcell_lane_read1 = "$flowcell"."_"."$lane"."_"."1";
  my $flowcell_lane_read2 = "$flowcell"."_"."$lane"."_"."2";

  #Store the lib name associated with each lib id
  $libs{$lib_id}{name} = $lib_name;

  if (defined($r1_full_path) && defined($r2_full_path) && defined($file_format) && defined($suffix)){
    #List of data files
    print DATAFILES "$lib_id\t$lib_name\t$flowcell\t$lane\t1\t$read_length\t$file_format\t$suffix\t$r1_full_path\n";
    print DATAFILES "$lib_id\t$lib_name\t$flowcell\t$lane\t2\t$read_length\t$file_format\t$suffix\t$r2_full_path\n";

    #hmmSplicer command
    print HMMSPLICER "source /home/malachig/.bashrc; /home/malachig/svn/hmmSplicer/hmmSplicerRun.pl  --read_file=$r1_full_path  --read_file_type=$file_format  --library_id=$lib_id  --flowcell_lane_read=$flowcell_lane_read1  --analysis_dir=$analysis_dir$project/  --ref_genome_db=$ref_genome  --python_path=/home/malachig/tools/python3/python_install/bin/python  --hmm_splicer_path=/home/malachig/tools/hmmSplicer/hmmSplicer-0.9.0/runHMM.py  --maq_dir=/home/malachig/tools/MAQ/maq-0.7.1/  --cpus=1  --quality_flag=p  --nobackup_dir=/tmp/\n";
    print HMMSPLICER "source /home/malachig/.bashrc; /home/malachig/svn/hmmSplicer/hmmSplicerRun.pl  --read_file=$r2_full_path  --read_file_type=$file_format  --library_id=$lib_id  --flowcell_lane_read=$flowcell_lane_read2  --analysis_dir=$analysis_dir$project/  --ref_genome_db=$ref_genome  --python_path=/home/malachig/tools/python3/python_install/bin/python  --hmm_splicer_path=/home/malachig/tools/hmmSplicer/hmmSplicer-0.9.0/runHMM.py  --maq_dir=/home/malachig/tools/MAQ/maq-0.7.1/  --cpus=1  --quality_flag=p  --nobackup_dir=/tmp/\n";

  }else{
    #List of broken dirs
    print BROKENDIRS "$lib_id\t$lib_name\t$flowcell\t$lane\t$data_path\n";
  }
}

#Create a merge bed command to merge the results for all lanes of each library
foreach my $lib (sort {$a cmp $b} keys %libs){
  my $name = $libs{$lib}{name};
  print MERGE "/home/malachig/svn/hmmSplicer/mergeJunctionBed.pl  --results_dir=$analysis_dir$project/results/$lib/  --dir_depth=1  --summary_name=$lib  --track_name=$name  --ref_genome=$ref_genome\n";
  print MERGE2 "/home/malachig/svn/hmmSplicer/mergeNonCanJunctionBed.pl  --results_dir=$analysis_dir$project/results/$lib/  --dir_depth=1  --summary_name=$lib  --track_name=$name  --ref_genome=$ref_genome  --min_read_count=2\n";
}

#Create a merge bed command to merge the results for the entire project
print MERGE "/home/malachig/svn/hmmSplicer/mergeJunctionBed.pl  --results_dir=$analysis_dir$project/results/  --dir_depth=2  --summary_name=$project  --track_name=$project  --ref_genome=$ref_genome\n";
print MERGE2 "/home/malachig/svn/hmmSplicer/mergeNonCanJunctionBed.pl  --results_dir=$analysis_dir$project/results/  --dir_depth=2  --summary_name=$project  --track_name=$project  --ref_genome=$ref_genome  --min_read_count=2\n";


#Create an annotate observed junctions command for the entire project
print ANNOTATE "/home/malachig/svn/hmmSplicer/annotateObservedJunctions.pl  --bedtools_bin_dir=/home/malachig/tools/BEDTools-Version-2.10.1/bin/  --analysis_dir=$analysis_dir  --project=$project  --ref_junction_file=$analysis_dir/ReferenceAnnotations/hg18/ALL.junc  --ref_ec_file=$analysis_dir/ReferenceAnnotations/hg18/ALL.ExonContent";

close(DATAFILES);
close(BROKENDIRS);
close(HMMSPLICER);
close(MERGE);
close(MERGE2);
close(ANNOTATE);


print "\n\n";

exit();


#################################################################################################################
#Get library list                                                                                               #
#################################################################################################################
sub getLibraryList{
  my %libs;
  open (LIBS, "$library_list") || die "\n\nCould not open library list: $library_list\n\n";
  while(<LIBS>){
    chomp($_);
    #my @line = split(" +", $_);
    my @line = split ' ',$_,2;
    $libs{$line[0]}{name}="unknown";
  }
  close(LIBS);
  return(\%libs);
}


#################################################################################################################
#Load a list of all lanes of WTSS data for each library                                                         #
#################################################################################################################
sub getLaneData{
  my %lanes;
  open (LANES, "$lims_data") || die "\n\nCould not open lims data: $lims_data\n\n";
  my $lc = 0;
  while(<LANES>){
    chomp($_);

    #Store useful data:
    #[0]=Name
    #[1]=Library Id
    #[2]=Flowcell
    #[3]=Lane
    #[4]=data_path
    #[6]=read_length
    #[7]=total_reads
    #[13]=billable
    #[16]=lims_run_validation
    #[18]=upper_protocol
    #[19]=lower_protocol
    #[21]=project

    my @line = split("\\|", $_);

    #Remove extra white space
    foreach my $l (@line){
      $l =~ s/ //g;
      $line[$#line]=$l;
    }

    my $lib_name = $line[0];
    my $lib_id = $line[1];
    my $flowcell = $line[2];
    my $lane = $line[3];
    my $data_path = $line[4];
    my $read_length = $line[6];
    my $billable = $line[13];
    my $lims_run_validation = $line[16];
    my $upper_protocol = $line[18];
    my $lower_protocol = $line[19];
    my $project = $line[21];

    #Unless the library is in the target list, skip it
    if (defined($lib_id)){
      unless ($libs{$lib_id}){
        next();
      }
    }else{
      next();
    }

    #Skip all lanes except for WTSS lanes ('SLX-Transcriptome')
    if(defined($upper_protocol)){
      unless($upper_protocol =~ /SLX\-Transcriptome/i){
        next();
      }
    }else{
      next();
    }

    #Skip lanes that are not billable
    unless($billable =~ /yes/i){
      next();
    }
    #Skip lanes that are not SLX-PET
    unless($lower_protocol =~ /SLX\-PET/i){
      next();
    }

    #Skip lanes that are not approved in the LIMS
    unless($lims_run_validation =~ /Approved/i){
      next();
    }

    #Skip lanes that do not have a data path defined yet
    if(defined($data_path)){
      unless($data_path =~ /\S+/){
        next();
      }
    }else{
      next();
    }

    my $sort_key = "$lib_id"."_"."$flowcell"."_"."$lane";

    #Store useful values
    $lc++;
    $lanes{$lc}{sort_key} = $sort_key;
    $lanes{$lc}{lib_name} = $lib_name;
    $lanes{$lc}{lib_id} = $lib_id;
    $lanes{$lc}{flowcell} = $flowcell;
    $lanes{$lc}{lane} = $lane;
    $lanes{$lc}{data_path} = $data_path;
    $lanes{$lc}{read_length} = $read_length;
    $lanes{$lc}{billable} = $billable;
    $lanes{$lc}{lims_run_validation} = $lims_run_validation;
    $lanes{$lc}{upper_protocol} = $upper_protocol ;
    $lanes{$lc}{lower_protocol} = $lower_protocol;
    $lanes{$lc}{project} = $project;
  }
  close(LANES);
  return(\%lanes);
}

