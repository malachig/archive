#!/usr/bin/perl -w
#Written by Malachi Griffith and Obi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Given a directory of mapping results, produce boxplots for position bias for data binned according to gene size
#i.e. boxplot of percent position of read for genes that are: 0-500, 500-1000, 1000-1500, ..., > 10000

#Also plot distribution of fragment sizes for paired-end reads where both reads were successfully mapped to the same gene

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use IO::File; 

#Load the ALEXA libraries
my $script_dir;
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
    $script_dir = $1;
  }
}
use utilities::utility qw(:all);

#Initialize command line options
my $data_dir = '';
my $working_dir = '';
my $results_dir = '';
my $record_limit = '';

GetOptions('data_dir=s'=>\$data_dir, 'working_dir=s'=>\$working_dir, 'results_dir=s'=>\$results_dir, 
           'record_limit=i'=>\$record_limit);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the complete path to a data directory containing summarized quality stats using: --data_dir", RESET;
print GREEN, "\n\tSpecify the complete path to a temporary working directory using: --working_dir", RESET;
print GREEN, "\n\tSpecify the complete path to a results directory using: --results_dir", RESET;
print GREEN, "\n\t\tMake sure that the two data types are the same if you do this!", RESET;
print GREEN, "\n\tSpecify the max number of records to process using: --record_limit [optional]", RESET;
print GREEN, "\n\nUsage: positionBiasVsGeneSize.pl  --data_dir=/projects/malachig/alexa_seq/read_records/HS04391/ENST_v53/  --working_dir=/projects/malachig/alexa_seq/figures_and_stats/HS04391/temp/  --results_dir=/projects/malachig/alexa_seq/figures_and_stats/HS04391/ENST_v53/\n\n", RESET;

my $column1 = 'R1_RelativePosition';
my $column2 = 'R2_RelativePosition';
my $column3 = 'R1_TranscriptSize';
my $column4 = 'R2_TranscriptSize';
my $column5 = 'DistanceBetweenReads_Transcript';

if($record_limit){
  print YELLOW, "\n\nLimiting analysis to the first $record_limit records processed\n", RESET;
}

#Check user supplied options
unless ($data_dir && $working_dir && $results_dir){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}

$working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"yes");
$results_dir = &checkDir('-dir'=>$results_dir, '-clear'=>"no");

#For debugging limit to only this number of files - set to 0 for no limit
my $test_limit = 0;

#Get list of input files containing data from the specified directory
my %files;
my $count = 0;
my @data_dirs = split (" ", $data_dir);
foreach my $dir (@data_dirs){
  $dir = &checkDir('-dir'=>$dir, '-clear'=>"no");
  &getDataFiles('-input_dir'=>$dir);
}

#Grab data from the specified column1 for all files in the input directory
#If specified by the user, also grab column2 and join these data to a single column
my $temp_file = $working_dir."temp_data.txt";
my $temp_file2 = $working_dir."temp_fragsize_data.txt";

#Create a hash to store counts of fragment sizes, keep track of total number of fragment sizes
my %frag_size_counts;
my $total_frag_count = 0;

#Use an array to set the transcript size ranges for binning data
#Open a temp file for each of these bins and store the file handles in a hash (keyed on the bin value)
my @sizes = qw(500 1000 2000 3000 4000 5000 10000 15000 20000);
my %counts;
foreach my $bin (@sizes){
  my %tmp;
  my $i_size = 0.1;
  my $max = 100.1;
  for (my $i = 0; $i < $max; $i += $i_size){ 
    my $is = sprintf("%.1f", $i);
    $tmp{$is}=0;
  }
  $counts{$bin}{positions} = \%tmp;
}
&parseDataColumn('-files_ref'=>\%files);

#Now write out all data to a temp file
#One row for each position value, one column for each transcript size bin
open(TEMP, ">$temp_file") || die "\n\nCould not open temp file: $temp_file\n\n";
print TEMP "Position @sizes\n";
my $max = 100.1;
my $i_size = 0.1;
for (my $i = 0; $i < $max; $i += $i_size){ 
  my $is = sprintf("%.1f", $i);
  print TEMP "$is";
  foreach my $bin_size (sort {$a <=> $b} keys %counts){
    my $pos_ref = $counts{$bin_size}{positions};
    my $count = $pos_ref->{$is};
    print TEMP " $count";
  }
  print TEMP "\n";
}
close(TEMP);

#Run R scripts to create positionBias plots
my $r_script = "$script_dir/R_bin/positionBias.R";
my $r_cmd = "$r_script $temp_file $results_dir";
print BLUE, "\n\nExecuting the R command:\n\n$r_cmd\n\n", RESET;
system($r_cmd);

#Write fragment sizes to temp file
#One row for each fragment size and count
if ($total_frag_count>0){
  open(TEMP2, ">$temp_file2") || die "\n\nCould not open temp file: $temp_file2\n\n";
  print TEMP2 "fragment_size\tcount\n";
  foreach my $frag_size (sort{$a<=>$b} keys %frag_size_counts){
    print TEMP2 "$frag_size\t$frag_size_counts{$frag_size}\n";
  }
  close(TEMP2);

  #Run R scripts to create fragment size distribution plots
  my $r_script2 = "$script_dir/R_bin/FragmentSizeDistribution.R";
  my $r_cmd2 = "$r_script2 $temp_file2 $results_dir";
  print BLUE, "\n\nExecuting the R command:\n\n$r_cmd2\n\n", RESET;
  system($r_cmd2);
}else{
  print RED, "\n\nNo fragment sizes found, skipping plot creation.\nCheck files in: $data_dir\n\n", RESET;
}

#Clean up temp files
my $rm_cmd = "rm -f $temp_file";
my $rm_cmd2 = "rm -f $temp_file2";
system($rm_cmd);
system($rm_cmd2);

exit();


###########################################################################################################
#Get data files and the columns of each                                                                   #
###########################################################################################################
sub getDataFiles{
  my %args = @_;
  my $dir = $args{'-input_dir'};

  my $dh = opendir(DIR, $dir) || die "\nCould not open directory: $dir\n\n";

  my @files = readdir(DIR);
  
  foreach my $file (@files){
    my %columns;
    my $header = 1;
    
    chomp($file);
    unless ($file =~ /\.txt\.gz$/){
      next();
    }
    if (-d $file){
      next();
    }
    $count++;

    $files{$count}{file_name} = $file;
    $files{$count}{file_path} = $dir.$file;

    #Get the header values for this file
    open (FILE, "zcat $dir$file |") || die "\nCould not open file: $dir$file";

    while(<FILE>){
      if ($header == 1){
        $header = 0;
        chomp($_);
        my @header = split("\t", $_);
        my $pos = 0;
        foreach my $head (@header){
          $columns{$head}{column_position} = $pos;
          $pos++;
        }
        last();
      }
    }
    close(FILE);
    $files{$count}{columns} = \%columns;

    if ($test_limit > 0 && $count == $test_limit){
      print RED, "\nLimiting this analysis to $test_limit files for debugging purposes\n\n", RESET;
      last();
    }
  }

  closedir(DIR);

  my $files_count = keys %files;
  print BLUE, "\n\nFound $files_count files to be processed (all .txt.gz files in the specified directory)", RESET;

  return();
}


###########################################################################################################
#Parse a specific data column1 (and column2 if specified) and produce a temp file with this data          #
#- Temp file will be created with only a single column of data                                            #
###########################################################################################################
sub parseDataColumn{
  my %args = @_;
  my $files_ref = $args{'-files_ref'};

  print BLUE, "\n\nProcessing files and parsing data: $column1 and $column2\n", RESET;

  my $counter = 0;
  my $fc = 0;

  foreach my $file (sort {$files_ref->{$a}->{file_name} cmp $files_ref->{$b}->{file_name}} keys %{$files_ref}){

    my $columns_ref = $files_ref->{$file}->{columns};

    unless ($columns_ref->{$column1} && $columns_ref->{$column2} && $columns_ref->{$column3} && $columns_ref->{$column4} && $columns_ref->{$column5}){
      print YELLOW, "\n\tFile: $files_ref->{$file}->{file_name} does not have the neccessary columns: $column1 $column2 $column3 $column4 $column5\n\n", RESET; 
      next();
    }
    $fc++;

    my $file_path = $files_ref->{$file}->{file_path};
    my $column1_pos = $columns_ref->{$column1}->{column_position};
    my $column2_pos = $columns_ref->{$column2}->{column_position};
    my $column3_pos = $columns_ref->{$column3}->{column_position};
    my $column4_pos = $columns_ref->{$column4}->{column_position};
    my $column5_pos = $columns_ref->{$column5}->{column_position};
   
    print BLUE, "\n\tProcessing: $files_ref->{$file}->{file_name}", RESET;
    open (FILE, "zcat $file_path |") || die "\nCould not open file: $file_path";

    my $first_line = 1;
    while(<FILE>){
      #Skip the header line
      if ($first_line == 1){
        $first_line = 0;
        next();
      }
      $counter++;

      #If a record limit was specified, watch out for it
      if ($record_limit){
        if($counter > $record_limit){
          close(FILE);
          last();
        }
      }

      chomp($_);
      my @line = split("\t", $_);
      my $pos1 = $line[$column1_pos];
      my $pos2 = $line[$column2_pos];
      my $size1 = $line[$column3_pos];
      my $size2 = $line[$column4_pos];
      my $frag_size = $line[$column5_pos];

      #Process Read 1 alignment - skip non-numerical values such as 'NA'
      if (($pos1 =~ /^\d+$|^\d+\.\d+$/) && ($size1 =~ /^\d+$/)){
 
        #Using the transcript size for this read alignment, place it in the correct bin and write it to the appropriate file
        my $bin = $sizes[scalar(@sizes)-1]; #Initialize bin to largest bin defined (will be used if a value is larger than the largest bin)
        foreach my $size_bin (@sizes){
          if ($size1 < $size_bin){
            $bin = $size_bin;
            my $pos_ref = $counts{$bin}{positions};
            if (defined($pos_ref->{$pos1})){
              $pos_ref->{$pos1}++;
            }
            last();
          }
        }
      }

      #Process Read 2 alignment
      if (($pos2 =~ /^\d+$|^\d+\.\d+$/) && ($size2 =~ /^\d+$/)){

        #Using the transcript size for this read alignment, place it in the correct bin and write it to the appropriate file
        my $bin = $sizes[scalar(@sizes)-1]; #Initialize bin to largest bin defined (will be used if a value is larger than the largest bin)
        foreach my $size_bin (@sizes){
          if ($size2 < $size_bin){
            $bin = $size_bin;
            my $pos_ref = $counts{$bin}{positions};
            if (defined($pos_ref->{$pos2})){
              $pos_ref->{$pos2}++;
            }
            last();
          }
        }
      }

      #Process fragment size from transcript alignment
      if ($frag_size =~ /^\d+$/){
        $frag_size_counts{$frag_size}++;
        $total_frag_count++;
      }
    }
    close(FILE);
  }

  return();
}









