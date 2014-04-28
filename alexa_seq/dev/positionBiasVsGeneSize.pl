#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Given a directory of mapping results, produce boxplots for position bias for data binned according to gene size
#i.e. boxplot of percent position of read for genes that are: 0-500, 500-1000, 1000-1500, ..., > 10000

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
my $data_type = '';
my $data_name = '';
my $units = '';
my $graph_type = '';
my $record_limit = '';

GetOptions('data_dir=s'=>\$data_dir, 'working_dir=s'=>\$working_dir, 'results_dir=s'=>\$results_dir, 
           'data_type=s'=>\$data_type, 'data_name=s'=>\$data_name, 'units=s'=>\$units, 'graph_type=s'=>\$graph_type, 'record_limit=i'=>\$record_limit);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the complete path to a data directory containing summarized quality stats using: --data_dir", RESET;
print GREEN, "\n\tSpecify the complete path to a temporary working directory using: --working_dir", RESET;
print GREEN, "\n\tSpecify the complete path to a results directory using: --results_dir", RESET;
print GREEN, "\n\t\tMake sure that the two data types are the same if you do this!", RESET;
print GREEN, "\n\tSpecify the data type using: --data_type (int/integer, float, string)", RESET;
print GREEN, "\n\tSpecify a text friendly version of this name using: --data_name (e.g. 'Fragment Size')", RESET;
print GREEN, "\n\tSpecify the name of the units for this data using: --units (e.g. 'bp') [optional]", RESET;
print GREEN, "\n\tSpecify the max number of records to process using: --record_limit [optional]", RESET;
print GREEN, "\n\nUsage: positionBiasVsGeneSize.pl  --data_dir=/projects/malachig/solexa/read_records/HS04391/ENST_v53/  --working_dir=/projects/malachig/solexa/figures_and_stats/HS04391/temp/  --results_dir=/projects/malachig/solexa/figures_and_stats/HS04391/ENST_v53/  --data_type=FLOAT  --data_name='MIP101 Library - Position Bias'  --units='Percent Position'\n\n", RESET;

my $column1 = 'R1_RelativePosition';
my $column2 = 'R2_RelativePosition';
my $column3 = 'R1_TranscriptSize';
my $column4 = 'R2_TranscriptSize';

if($record_limit){
  print YELLOW, "\n\nLimiting analysis to the first $record_limit records processed\n", RESET;
}

#Check user supplied options
unless ($data_dir && $working_dir && $results_dir && $data_type && $units && $data_name){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}
chomp($data_type);
unless($data_type =~ /int|integer|float|string/i){
  print RED, "\n\n--data_type specified is not a supported option!  Use one of : int|integer|float|string\n\n", RESET;
  exit();
}
chomp($units);

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

#Use an array to set the transcript size ranges for binning data
#Open a temp file for each of these bins and store the file handles in a hash (keyed on the bin value)
my @sizes = qw(500 1000 2000 3000 4000 5000 10000 15000 20000);
my %outfiles;
foreach my $size (@sizes){
  my $temp_file = "$working_dir"."Temp_$size"."bp.txt";
  my $fh = IO::File->new(">$temp_file") || die "\nCould not open temp file: $temp_file\n\n";
  $outfiles{$size}{file_path} = $temp_file;
  $outfiles{$size}{fh} = $fh;
  $outfiles{$size}{record_count} = 0;
}

#Write the header value for each of these temp files
foreach my $size (sort {$a <=> $b} keys %outfiles){
  my $fh = $outfiles{$size}{fh};
  print $fh "$size\n";
}


my @return = @{&parseDataColumn('-files_ref'=>\%files, '-data_type'=>$data_type)};
my $column_count = $return[0];
my $record_count = $return[1];

my $r_script = "$script_dir/R_bin/positionBias.R";
my $r_cmd = "$r_script $temp_file $column_count $record_count \'$data_name\' \'$units\' $results_dir";

print BLUE, "Executing the R command:\n\n$r_cmd\n\n", RESET;
system($r_cmd);


#Clean up temp file
my $rm_cmd = "rm -f $temp_file";
system($rm_cmd);

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
  my $data_type = $args{'-data_type'};

  print BLUE, "\n\nProcessing files and parsing data: $column1 and $column2 to temp file: $temp_file\n", RESET;

  my $counter = 0;
  foreach my $file (sort {$files_ref->{$a}->{file_name} cmp $files_ref->{$b}->{file_name}} keys %{$files_ref}){

    my $columns_ref = $files_ref->{$file}->{columns};

    unless ($columns_ref->{$column1} && $columns_ref->{$column2} && $columns_ref->{$column3} && $columns_ref->{$column4}){
      print YELLOW, "\n\tFile: $files_ref->{$file}->{file_name} does not have the neccessary columns: $column1 $column2 $column3 $column4\n\n", RESET; 
      next();
    }

    my $file_path = $files_ref->{$file}->{file_path};
    my $column1_pos = $columns_ref->{$column1}->{column_position};
    my $column2_pos = $columns_ref->{$column2}->{column_position};
    my $column3_pos = $columns_ref->{$column3}->{column_position};
    my $column4_pos = $columns_ref->{$column4}->{column_position};
   
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

      #DEBUG: test on single gene with known 3' end bias
      #unless ($line[5] eq 10108 && $line[18] eq 10108){
      #   next();
      #}

      #Process Read 1 alignment - skip non-numerical values such as 'NA'
      if (($pos1 =~ /^\d+$|^\d+\.\d+$/) && ($size1 =~ /^\d+$/)){
 
        #Using the transcript size for this read alignment, place it in the correct bin and write it to the appropriate file
        my $bin = $sizes[scalar(@sizes)-1]; #Initialize bin to largest bin defined (will be used if a value is larger than the largest bin)
        foreach my $size_bin (@sizes){
          if ($size1 < $size_bin){
            $bin = $size_bin;
            last();
          }
        }
        my $fh = $outfiles{$bin}{fh};
        $outfiles{$bin}{record_count}++;
        print $fh "$pos1\n";
      }

      #Process Read 2 alignment
      if (($pos2 =~ /^\d+$|^\d+\.\d+$/) && ($size2 =~ /^\d+$/)){

        #Using the transcript size for this read alignment, place it in the correct bin and write it to the appropriate file
        my $bin = $sizes[scalar(@sizes)-1]; #Initialize bin to largest bin defined (will be used if a value is larger than the largest bin)
        foreach my $size_bin (@sizes){
          if ($size2 < $size_bin){
            $bin = $size_bin;
            last();
          }
        }
        my $fh = $outfiles{$bin}{fh};
        $outfiles{$bin}{record_count}++;
        print $fh "$pos2\n";
      }

    }
    close(FILE);
  }

  #Determine the max record count for all bins
  my $max_record_count = 0;
  foreach my $size (sort {$a <=> $b} keys %outfiles){
    if ($outfiles{$size}{record_count} > $max_record_count){
      $max_record_count = $outfiles{$size}{record_count};
    }
  }

  print BLUE, "\nFound and wrote a total of $max_record_count valid records to the largest temp file\n\n", RESET;

  #create a string listing all files in order
  my $file_string = '';
  foreach my $size (sort {$a <=> $b} keys %outfiles){
    $file_string = "$file_string"."$outfiles{$size}{file_path} ";
  }

  #Paste all temp files together in order as a single command
  my $cmd_paste = "paste -d \',\' $file_string > $temp_file";
  system($cmd_paste);
  print YELLOW, "\n\n$cmd_paste\n\n", RESET;

  #Delete all temp files
  my $cmd_rm = "rm -f $file_string";
  system($cmd_rm);
  print YELLOW, "\n\n$cmd_rm\n\n", RESET;

  my $column_count = scalar(@sizes);
  print BLUE, "\n\nCreated a multi-column temp file with $column_count columns and $max_record_count lines\n\n", RESET;

  my @return = ($column_count, $max_record_count);
  return(\@return);
}









