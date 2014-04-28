#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#1.) Import a genome regions file to determine the number of results files expected for each library/feature combination
#2.) For each library specified and feature type, go to the log directory and count the number of log files that contain "SCRIPT COMPLETE"
#3.) If the correct number of files is found, go to the results files and join all results files into one master file
#4.) Delete the log files
#5.) Delete the temporary results files

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $source_dir = '';
my $batch_dir = '';
my $libraries = '';
my $ensembl_version = '';
my $annotation_dir = '';
my $junction_seq_size = '';
my $boundary_seq_size = '';
my $regions_file_version = '';
my $no_delete = '';

GetOptions ('source_dir=s'=>\$source_dir, 'batch_dir=s'=>\$batch_dir, 'libraries=s'=>\$libraries, 'ensembl_version=i'=>\$ensembl_version, 'annotation_dir=s'=>\$annotation_dir, 'junction_seq_size=i'=>\$junction_seq_size, 'boundary_seq_size=i'=>\$boundary_seq_size, 'regions_file_version=i'=>\$regions_file_version, 'no_delete=s'=>\$no_delete);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the source dir for the analysis using:  --source_dir", RESET;
print GREEN, "\n\tSpecify the directory containing batch expression jobs using: -batch_dir", RESET;
print GREEN, "\n\tSpecify the libraries to be checked as a comma seperated list using:  --libraries", RESET;
print GREEN, "\n\tSpecify the ensembl version used for this analysis using: --ensembl_version", RESET;
print GREEN, "\n\tSpecify the path to the annotation files using: --annotation_dir", RESET;
print GREEN, "\n\tSpecify the junction database sequence length using:  --junction_seq_size", RESET;
print GREEN, "\n\tSpecify the boundary database sequence length using:  --boundary_seq_size", RESET;
print GREEN, "\n\tSpecify the regions file version using: --regions_file_version (e.g. 50, 250)", RESET;
print GREEN, "\n\tFor testing purposes you can prevent log and temp files from being deleted by using: --no_delete=1", RESET;
print GREEN, "\n\nExample:  checkExpressionJobs.pl  --source_dir=/projects/malachig/solexa/  --batch_dir=/projects/malachig/solexa/batch_jobs/Neuroblastoma/  --libraries='HS0499,HS0502'  --ensembl_version=53  --annotation_dir=/projects/malachig/sequence_databases/hs_53_36o/  --junction_seq_size=62  --boundary_seq_size=62  --regions_file_version=50  --no_delete=1\n\n", RESET;

unless ($source_dir && $batch_dir && $libraries && $ensembl_version && $regions_file_version && $annotation_dir && $junction_seq_size && $boundary_seq_size){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}
$source_dir = &checkDir('-dir'=>$source_dir, '-clear'=>"no");
$batch_dir = &checkDir('-dir'=>$batch_dir, '-clear'=>"no");
$annotation_dir = &checkDir('-dir'=>$annotation_dir, '-clear'=>"no");

my $test1 = 1;
my $test2 = 1;

#1.) Import a genome regions file to determine the number of results files expected for each library/feature combination
my $regions_file = "$annotation_dir"."Regions_"."$regions_file_version"."_Genes.txt";
print BLUE, "\n\nGetting regions ...", RESET;
open(REGIONS, "$regions_file") || die "\nCould not open regions file: $regions_file\n\n";
my $r = 0;
my %regions;
my $header = 1;
while(<REGIONS>){
  if ($header == 1){$header = 0; next();}
  $r++;
  chomp($_);
  my @data = split("\t", $_);
  $regions{$r}{chromosome} = $data[0];
  $regions{$r}{region} = $data[1];
  $regions{$r}{chr_start} = $data[2];
  $regions{$r}{chr_end} = $data[3];
  $regions{$r}{size} = $data[4];
  $regions{$r}{gene_count} = $data[5];
}
close(REGIONS);
my $jobs_expected = keys %regions;

print BLUE, "\n\nBased on the specified regions file we are expecting $jobs_expected jobs for each library/feature combination", RESET;


my @libs = split(",", $libraries);


my @feature_dirs = qw(ENST Junctions Boundaries Introns Intergenics);
my %repeat_jobs;
$r = 0;

foreach my $lib (@libs){
  print BLUE, "\n\nProcessing library: $lib", RESET;
  foreach my $feature_dir (@feature_dirs){

    print BLUE, "\n\n\tProcessing feature dir: $feature_dir", RESET;
    my $feature_dir_full = "$feature_dir"."_v"."$ensembl_version";
    my $log_dir = "$source_dir"."logs/$lib/generateExpressionValues/$feature_dir_full/";
    my $cmd = "tail $log_dir"."*.txt | grep \"SCRIPT COMPLETE\" | wc -l";
    my $results_dir = "$source_dir"."read_records/$lib/$feature_dir_full/Summary/results/";

    opendir(DIRHANDLE, "$log_dir") || die "\nCannot open directory: $log_dir\n\n";
    my @tmp_files = readdir(DIRHANDLE);
    closedir(DIRHANDLE);
    unless (scalar(@tmp_files) > 2){
      print MAGENTA, "\n\tNo log files found, nothing to do...", RESET;
      next();
    }

    print MAGENTA, "\n\t$cmd = ", RESET;
    my $result = `$cmd`;
    chomp($result);
    print MAGENTA, "$result", RESET;

    if ($result == 0){
      #If the directory is empty, there is no point in proceeding
      print MAGENTA, "\n\tNo log files found, nothing to do...", RESET;
      next();
      
    }elsif ($result == $jobs_expected){
      #If all the log files were found and indicated job completion, proceed to the next steps
      print MAGENTA, "\n\tFound the expected number of jobs ... proceeding to joining and clean-up step...", RESET;
      
      #a.) get list of files to be joined, divided by file types
      opendir(DIRHANDLE, "$results_dir") || die "\nCannot open directory: $results_dir\n\n";
      my @files = readdir(DIRHANDLE);
      closedir(DIRHANDLE);  
      my %file_types;
      foreach my $file(@files){
        my $file_path = "$results_dir"."$file";
        if ($file =~ /\S+\_\S+\_(\w+)\.txt/){
          if ($file_types{$1}){
            $file_types{$1}{count}++;
            push(@{$file_types{$1}{file_list}}, $file_path);
          }else{
            $file_types{$1}{count} = 1;
            my @tmp;
            push(@tmp, $file_path);
            $file_types{$1}{file_list} = \@tmp;
            my $cmd = "head -n 1 $file_path";
            my $header = `$cmd`;
            chomp($header);
            $file_types{$1}{header} = $header;
          }
        }
      }
      #b.) join files of each file type.  Skip headers
      foreach my $file_type (keys %file_types){
        print MAGENTA, "\n\t\tJoining files of type: $file_type", RESET;
        my @file_list = @{$file_types{$file_type}{file_list}};
        my $tmp_file = "$source_dir"."read_records/$lib/$feature_dir_full/Summary/$lib"."_"."$file_type"."_v"."$ensembl_version"."_TEMP.txt";
        my $summary_file = "$source_dir"."read_records/$lib/$feature_dir_full/Summary/$lib"."_"."$file_type"."_v"."$ensembl_version".".txt";
        open (TMP, ">$tmp_file") || die "\nCould not open temp output file: $tmp_file\n\n";
        print MAGENTA, "\n\t\t\t$tmp_file", RESET;
        foreach my $file (@file_list){
          my $header = 1;
          open(FILE, "$file") || die "\nCould not open file: $file\n\n";
          while(<FILE>){
            if ($header == 1){$header = 0; next();}
            print TMP "$_";
          }
          close(FILE);
        }
        close(TMP);

        #Sort the resulting file and add the header to the top
        my $sort_file = "$source_dir"."read_records/$lib/$feature_dir_full/Summary/$lib"."_"."$file_type"."_v"."$ensembl_version"."_SORT.txt";
        my $sort_cmd = "sort $tmp_file > $sort_file";
        system($sort_cmd);
        open (OUT, ">$summary_file") || die "\nCould not open summary output file: $summary_file\n\n";
        open (SORT, "$sort_file") || die "\nCould not open sort output file: $summary_file\n\n";
        print OUT "$file_types{$file_type}{header}\n";
        while(<SORT>){
          print OUT "$_";
        }
        close(OUT);
        close(SORT);
        my $rm_cmd1 = "rm -f $tmp_file";
        my $rm_cmd2 = "rm -f $sort_file";
        system($rm_cmd1);
        system($rm_cmd2);
      }

      #Clean up tmp results files and log files
      if ($no_delete){
        print MAGENTA "\n\t\tNo-delete option specified - leaving temp results files and log files in place\n\n", RESET;
      }else{
        print MAGENTA "\n\t\tCleaning up temp results files and log files", RESET;
        my $log_dir = "$source_dir"."logs/$lib/generateExpressionValues/$feature_dir_full/";
        #system ("rm -fr $results_dir");
        #system ("mkdir $results_dir");
        system ("rm -fr $log_dir");
        system ("mkdir $log_dir");
        print "\n\n";
      }

    }else{
      #If some log files were missing or indicated the job was not complete, we need to find those jobs
      #Note that we are processing one library at a time here, but the batch files contained jobs for both libraries!

      my $missing_jobs = $jobs_expected-$result;
      print MAGENTA, "\n\t$missing_jobs jobs were missing or incomplete.  Adding incomplete jobs to repair list ...", RESET;

      my $batch_file = "$batch_dir"."generateExpressionValues_"."$feature_dir"."_v"."$ensembl_version".".sh";
      unless (-e $batch_file){
        print RED, "\nCould not find batch file: $batch_file\n\n", RESET;
        exit();
      }
      my %jobs;
      my $c = 0;
      open(BATCH, "$batch_file") || die "\nCould not open batch_file: $batch_file";
      while(<BATCH>){

        #Unless the current job corresponds to the current library skip it
        unless ($_ =~ /$lib/){
          next();
        }

        chomp($_);
        $c++;
        $jobs{$c}{command} = $_;
        if ($_ =~ /\-\-log\_file\=(\S+)/){
          $jobs{$c}{logfile} = $1;
        }else{
          print RED, "\nCould not find log file name in command: $_\n\n", RESET;
          exit();
        }
      }
      close(BATCH);
      foreach my $job (sort {$a <=> $b} keys %jobs){
        if (-e $jobs{$job}{logfile}){
          my $cmd = "tail $jobs{$job}{logfile} | grep \"SCRIPT COMPLETE\" | wc -l";
          my $result = `$cmd`;
          unless ($result >= 1){
            $r++;
            $repeat_jobs{$r}{command} = $jobs{$job}{command};
          }
        }else{
          $r++;
          $repeat_jobs{$r}{command} = $jobs{$job}{command};
        }
      }
    }
  }
}

my $repeat_job_count = keys %repeat_jobs;
my $repeat_job_file = "$batch_dir"."generateExpressionValues_REPEAT.sh";
if ($repeat_job_count >= 1){
  print BLUE, "\n\nPrinting a total of $repeat_job_count jobs to the REPEAT batch file: $repeat_job_file\n\n", RESET;
  open (REPEAT, ">$repeat_job_file") || die "\nCould not open REPEAT job file: $repeat_job_file\n\n";
  foreach my $job (sort {$a <=> $b} keys %repeat_jobs){
    print REPEAT "$repeat_jobs{$job}{command}\n";
  }
  close(REPEAT);
  $test1 = 0;
}else{
  print BLUE, "\n\nThere were no failed jobs!\n\n", RESET;
}


#Finally check the joined expression results files that have been created to see if they have the correct number of rows
print BLUE, "\n\nChecking number of records in joined files versus annotation files", RESET;
my %results_files;
$results_files{1}{feature_dir} = "ENST";
$results_files{1}{data_type} = "GeneExpression";
$results_files{1}{annotation_file} = "genes/genes_annotated.txt.gz";
$results_files{2}{feature_dir} = "ENST";
$results_files{2}{data_type} = "ExonRegionExpression";
$results_files{2}{annotation_file} = "exonRegions/exonRegions_annotated.txt.gz";
$results_files{3}{feature_dir} = "Junctions";
$results_files{3}{data_type} = "JunctionExpression";
$results_files{3}{annotation_file} = "exonJunctions/exonJunctions_"."$junction_seq_size"."mers_annotated.txt.gz";
$results_files{4}{feature_dir} = "Boundaries";
$results_files{4}{data_type} = "BoundaryExpression";
$results_files{4}{annotation_file} = "exonBoundaries/exonBoundaries_"."$boundary_seq_size"."mers_annotated.txt.gz";
$results_files{5}{feature_dir} = "Introns";
$results_files{5}{data_type} = "IntronExpression";
$results_files{5}{annotation_file} = "introns/introns_annotated.txt.gz";
$results_files{6}{feature_dir} = "Introns";
$results_files{6}{data_type} = "ActiveIntronRegionExpression";
$results_files{6}{annotation_file} = "introns/activeIntronRegions.txt.gz";
$results_files{7}{feature_dir} = "Introns";
$results_files{7}{data_type} = "SilentIntronRegionExpression";
$results_files{7}{annotation_file} = "introns/silentIntronRegions.txt.gz";
$results_files{8}{feature_dir} = "Intergenics";
$results_files{8}{data_type} = "IntergenicExpression";
$results_files{8}{annotation_file} = "intergenics/intergenics_annotated.txt.gz";
$results_files{9}{feature_dir} = "Intergenics";
$results_files{9}{data_type} = "ActiveIntergenicRegionExpression";
$results_files{9}{annotation_file} = "intergenics/activeIntergenicRegions.txt.gz";
$results_files{10}{feature_dir} = "Intergenics";
$results_files{10}{data_type} = "SilentIntergenicRegionExpression";
$results_files{10}{annotation_file} = "intergenics/silentIntergenicRegions.txt.gz";

foreach my $result (sort {$a <=> $b} keys %results_files){

  #Get the expected number of records from the annotation file
  my $feature_dir = $results_files{$result}{feature_dir};
  my $data_type = $results_files{$result}{data_type};
  my $annotation_file = $results_files{$result}{annotation_file};
  my $annotation_path = "$annotation_dir"."$annotation_file";
  print YELLOW, "\n\n$feature_dir - $data_type ($annotation_path) = ", RESET;
  my $line_count = `zcat $annotation_path | wc -l`;
  chomp($line_count);
  print YELLOW, "$line_count", RESET;

  #Now check the number of records found in the result file for each library
  my @counts;
  foreach my $lib (@libs){
    my $feature_dir_full = "$feature_dir"."_v"."$ensembl_version";
    my $results_file = "$source_dir"."read_records/"."$lib"."/$feature_dir_full"."/Summary/"."$lib"."_"."$data_type"."_v"."$ensembl_version".".txt";
    print YELLOW, "\n\t\t$lib\t$results_file = ", RESET;

    if (-e $results_file){
      my $line_count = `cat $results_file | wc -l`;
      chomp($line_count);
      print YELLOW, "$line_count", RESET;
      push(@counts, $line_count);
    }else{
      print YELLOW, "0", RESET;
      push(@counts, 0);
    }
  }
  my $test_count = $counts[0];
  my $test = 1;
  foreach my $count (@counts){
    unless ($count == $test_count && $count > 1){
      $test = 0;
      $test2 = 0;
    }
  }
  unless($test){
    print RED, "\n\t\tCounts are not equal across libraries: (@counts) - Do not proceed until this is resolved", RESET;
  }
}

if ($test1 && $test2){
  print BLUE, "\n\nEverything looks good - should be safe to proceed", RESET;
}else{
  print RED, "\n\nData is missing - do not proceed - check for running or failed jobs", RESET;
}

print "\n\n";

exit();

