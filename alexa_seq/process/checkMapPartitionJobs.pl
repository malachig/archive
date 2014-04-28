#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to check that partitioning of map files into genome regions has completed successfully
#Success will be determined by checking for completion status in log files and also by checking for the correct number of region sub directories

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use BerkeleyDB;
use IO::File;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $analysis_dir = '';
my $ensembl_version = '';
my $project_name = '';
my $partition_file = '';

GetOptions ('analysis_dir=s'=>\$analysis_dir, 'ensembl_version=i'=>\$ensembl_version, 'project_name=s'=>\$project_name, 'partition_file=s'=>\$partition_file);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the analysis dir using: --analysis_dir", RESET;
print GREEN, "\n\tSpecify the ensembl version using:  --ensembl_version", RESET;
print GREEN, "\n\tSpecify the project_name using: --project_name", RESET;
print GREEN, "\n\tSpecify the partition_file using: --partition_file", RESET;

print GREEN, "\n\nExample: checkMapPartitionJobs.pl  --analysis_dir=/projects/malachig/alexa_seq/  --ensembl_version=53  --project_name=Breast  --partition_file=/projects/alexa/sequence_databases/hs_53_36o/Regions_250_Genes.txt\n\n", RESET;

unless ($analysis_dir && $ensembl_version && $project_name && $partition_file){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Check input files and dirs
$analysis_dir = &checkDir('-dir'=>$analysis_dir, '-clear'=>"no");
unless (-e $partition_file){
  print RED, "\n\nCould not find partition file: $partition_file", RESET;
  exit();
}
my $library_list_file = "$analysis_dir"."batch_jobs/$project_name/$project_name"."_Lib_Names.txt";
unless (-e $library_list_file){
  print RED, "\n\nCould not find library list file: $library_list_file", RESET;
  exit();
}

#Define list of classes to check
my %class;
$class{'ENST'}{dir_name} = "ENST_v$ensembl_version";
$class{'ENST'}{order} = 1;
$class{'NOVEL_JUNCTION'}{dir_name} = "Junctions_v$ensembl_version";
$class{'NOVEL_JUNCTION'}{order} = 0;
$class{'NOVEL_BOUNDARY'}{dir_name} = "Boundaries_v$ensembl_version";
$class{'NOVEL_BOUNDARY'}{order} = 0;
$class{'INTRON'}{dir_name} = "Introns_v$ensembl_version";
$class{'INTRON'}{order} = 0;
$class{'INTERGENIC'}{dir_name} = "Intergenics_v$ensembl_version";
$class{'INTERGENIC'}{order} = 0;


#Get the libraries to be checked for this library
open (LIBS, "$library_list_file") || die "\n\nCould not open library list file: $library_list_file\n\n";
my %libs;
my $lc = 0;
while(<LIBS>){
  $lc++;
  chomp($_);
  my @line = split ("\t", $_);
  if ($_ =~ /\w+/){
    $libs{$lc}{library_id} = $line[0];
  }
}
close(LIBS);

#Get the region list
open (REGIONS, "$partition_file") || die "\n\nCould not open partitions file: $partition_file\n\n";
my $header = 1;
my %regions;
my $r = 0;
while(<REGIONS>){
  chomp($_);
  my @line = split ("\t", $_);
  if ($header == 1){
    $header = 0;
    next();
  }
  $r++;
  $regions{$r}{chromosome} = $line[0];
  $regions{$r}{region} = $line[1];
}
close(REGIONS);
my $target_dirs = keys %regions;


#For each library, check for completion status in the log files
#For each mapping class check for existence of each region dir for each mapping type
print MAGENTA, "\n\nChecking for completion status in log files:", RESET;
my $grand_test = 1;
foreach my $lc (sort {$a <=> $b} keys %libs){
  my $library_id = $libs{$lc}{library_id};
  print BLUE, "\n\nChecking library: $library_id", RESET;

  foreach my $class (sort {$class{$a}->{order} <=> $class{$b}->{order}} keys %class){
    print BLUE, "\n\t$class", RESET;

    #Check log file
    my $log_file = "$analysis_dir"."logs/$library_id/partitionMapFiles/partitionMapFiles_"."$class"."_LOG.txt";
    my $test = 0;
    if (-e $log_file){
      open(LOG, "tail $log_file |") || die "\n\nCould not open log file: $log_file\n\n";
      while(<LOG>){
        if ($_ =~ /SCRIPT\sCOMPLETE/){
          $test = 1;
        }
      }
      close(LOG);
    }

    if ($test){
      print BLUE "\n\t\tJob complete according to log file", RESET;
    }else{
      print YELLOW, "\n\t\tLog file missing or incomplete", RESET;
      $grand_test = 0;
    }

    #Check partition dirs
    my $dirs_found = 0;
    my $class_dir_name = $class{$class}{dir_name};
    my $base_dir = "$analysis_dir"."read_records/$library_id/$class_dir_name/part/";
    foreach my $r (sort {$a <=> $b} keys %regions){
      my $sub_dir = "$base_dir"."$regions{$r}{chromosome}"."_"."$regions{$r}{region}";
      if (-e $sub_dir && -d $sub_dir){
        $dirs_found++;
      }else{
        print YELLOW, "\n\tCould not find dir: $sub_dir", RESET
      }
    }

    if ($dirs_found == $target_dirs){
      print BLUE "\n\t\tFound expected number of partition dirs", RESET;
    }else{
      print YELLOW, "\n\t\tPartition dirs missing (expecting $target_dirs but found $dirs_found)", RESET;
      $grand_test = 0;
    }

  }
}




if ($grand_test){
  print BLUE, "\n\n\nEverything looks good - should be safe to proceed", RESET;
}else{
  print RED, "\n\n\nFiles are missing - do not proceed - check for running or failed jobs", RESET;
}



print "\n\n";
exit();








