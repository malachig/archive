#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Given a directory of mapping results, produce basic statistical summaries and graphs

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

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
my $column1 = '';
my $column2 = '';
my $data_type = '';
my $data_name = '';
my $units = '';
my $graph_type = '';
my $breaks = '';

GetOptions('data_dir=s'=>\$data_dir, 'working_dir=s'=>\$working_dir, 'results_dir=s'=>\$results_dir, 
           'column1=s'=>\$column1, 'column2=s'=>\$column2, 
           'data_type=s'=>\$data_type, 'data_name=s'=>\$data_name, 'units=s'=>\$units, 'graph_type=s'=>\$graph_type, 'breaks=i'=>\$breaks);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the complete path to a data directory containing summarized quality stats using: --data_dir", RESET;
print GREEN, "\n\tTo include multiple data directories, specify as a space seperated list", RESET;
print GREEN, "\n\tSpecify the complete path to a temporary working directory using: --working_dir", RESET;
print GREEN, "\n\tSpecify the complete path to a results directory using: --results_dir", RESET;
print GREEN, "\n\tSpecify the name of the data field (column1) you wish to summarize using:  --column1 (e.g. DistanceBetweenReads_Transcript)", RESET;
print GREEN, "\n\t\tYou can also join this column of data to a second column by specifying: --column2 (e.g. to combine R1 and R2 values)", RESET;
print GREEN, "\n\t\tMake sure that the two data types are the same if you do this!", RESET;
print GREEN, "\n\tSpecify the data type using: --data_type (int/integer, float, string)", RESET;
print GREEN, "\n\tSpecify a text friendly version of this name using: --data_name (e.g. 'Fragment Size')", RESET;
print GREEN, "\n\tSpecify the name of the units for this data using: --units (e.g. 'bp') [optional]", RESET;
print GREEN, "\n\tSpecify the type of graph you wish to generate using: --graph_type (Hist, Bplot, Pie, Bar, Dot, Stats)", RESET;
print GREEN, "\n\tIf you select 'Hist' you can also specify how many bins you want the data divided into using: --breaks (default is 100)", RESET;
print GREEN, "\n\nUsage: summarizeMapping.pl  --data_dir=/projects/malachig/solexa/read_records/HS04391/ENST_v49/  --working_dir=/projects/malachig/solexa/figures_and_stats/HS04391/temp/  --results_dir=/projects/malachig/solexa/figures_and_stats/HS04391/ENST_v49/  --column1='DistanceBetweenReads_Transcript'  --data_type=INT  --data_name='MIP101 Library - Fragment Size'  --units='bp' --graph_type='Hist'\n\n", RESET;

#Check user supplied options
unless ($data_dir && $working_dir && $results_dir && $column1 && $data_type && $units && $data_name && $graph_type){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}
unless($column2){
  $column2 = "NA";
}
chomp($data_type);
unless($data_type =~ /int|integer|float|string/i){
  print RED, "\n\n--data_type specified is not a supported option!  Use one of : int|integer|float|string\n\n", RESET;
  exit();
}
chomp($graph_type);
unless($graph_type =~ /hist|histogram|bplot|boxplot|pie|bar|dot|stats/i){
  print RED, "\n\n--graph_type specified is not a supported option!  Use one of : hist|histogram|bplot|boxplot|pie|stats\n\n", RESET;
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
my @data_dirs = split(" ", $data_dir);
foreach my $data_dir (@data_dirs){
  $data_dir = &checkDir('-dir'=>$data_dir, '-clear'=>"no");
  &getDataFiles('-input_dir'=>$data_dir);
}

#Grab data from the specified column1 for all files in the input directory
#If specified by the user, also grab column2 and join these data to a single column
my $temp_file = $working_dir."temp_data.txt";

#Generate figures according to user specified option

#HISTOGRAMS - applied to one massive column of data (all lanes combined)
if ($graph_type =~ /hist|histogram/i){
  
  #If the user specified a breaks option - use this, otherwise default to 100
  unless ($breaks){
    $breaks = 100;
  }

  my $data_count = &parseDataColumn_single('-files_ref'=>\%files, '-column1'=>$column1, '-column2'=>$column2, '-data_type'=>$data_type, '-temp_file'=>$temp_file);

  my $r_script = "$script_dir/R_bin/hist.R";
  my $r_cmd = "$r_script $temp_file $data_count \'$data_name\' \'$units\' $results_dir $breaks";

  print BLUE, "Executing the R command:\n\n$r_cmd\n\n", RESET;
  system($r_cmd);
}

#PIE/DOT/BAR CHARTS - applied to one massive column of data (all lanes combined)
if ($graph_type =~ /pie|bar|dot/i){
  my $data_count = &parseDataColumn_single('-files_ref'=>\%files, '-column1'=>$column1, '-column2'=>$column2, '-data_type'=>$data_type, '-temp_file'=>$temp_file);

  my $r_script = "$script_dir/R_bin/pie.R";
  my $r_cmd = "$r_script $temp_file $data_count \'$data_name\' \'$units\' $results_dir";

  print BLUE, "Executing the R command:\n\n$r_cmd\n\n", RESET;
  system($r_cmd);
}

#BOXPLOTS - applied lane-wise (i.e. one column per lane)
if ($graph_type =~ /bplot|boxplot/i){
  my @counts = @{&parseDataColumn_multi('-files_ref'=>\%files, '-column1'=>$column1, '-column2'=>$column2, '-data_type'=>$data_type, '-temp_file'=>$temp_file)};

  my $r_script = "$script_dir/R_bin/bplot.R";
  my $r_cmd = "$r_script $temp_file $counts[0] $counts[1] \'$data_name\' \'$units\' $results_dir";

  print BLUE, "Executing the R command:\n\n$r_cmd\n\n", RESET;
  system($r_cmd);

  #Print a text file to act as a key for lane names versus the simple number used in graphs
  if ($column2 eq 'NA'){
    my $key_file = $results_dir."singleRead_LaneKey.txt";
    open(KEY, ">$key_file") || die "\nCould not open key file: $key_file\n\n";

    my $count = 0;
    foreach my $file (sort {$files{$a}->{file_name} cmp $files{$b}->{file_name}} keys %files){
      $count++;
      print KEY "$count\t$files{$file}{file_name}\n";
    }
    close(KEY);
  }else{
    my $key_file = $results_dir."doubleRead_LaneKey.txt";
    open(KEY, ">$key_file") || die "\nCould not open key file: $key_file\n\n";
    
    my $count = 0;
    foreach my $file (sort {$files{$a}->{file_name} cmp $files{$b}->{file_name}} keys %files){
      $count++;
      print KEY "$count\t$files{$file}{file_name} - R1\n";
      $count++;
      print KEY "$count\t$files{$file}{file_name} - R2\n";
    }
    close(KEY);
  }
}

#BASIC STATS - Calculated within Perl script
#Number and percent of read pairs which were mapped AS a pair to the same gene (i.e. Read1 and Read2 are 'Top Hits' and map to the same gene)
#Number and percent of read pairs which are mapped to the same gene and have the following strand orientations: 
#FIRST: +/+, -/-, +/-, -/+
#SECOND: right (+/- OR -/+) vs. wrong (+/+, -/-)
if ($graph_type =~ /stats/i){
  &parseData_stats('-files_ref'=>\%files, '-results_dir'=>$results_dir);
}


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
sub parseDataColumn_single{
  my %args = @_;
  my $files_ref = $args{'-files_ref'};
  my $column1 = $args{'-column1'};
  my $column2 = $args{'-column2'};
  my $data_type = $args{'-data_type'};
  my $temp_file = $args{'-temp_file'};

  my $valid_record_count = 0;

  print BLUE, "\n\nProcessing files and parsing data: $column1 to temp file: $temp_file\n", RESET;

  unless($column2 eq "NA"){
    print BLUE, "\tAlso parsing and joining data: $column2 to temp file: $temp_file\n", RESET;
  }

  open (TEMP, ">$temp_file") || die "\nCould not open temp file: $temp_file\n\n";
  
  foreach my $file (sort {$files_ref->{$a}->{file_name} cmp $files_ref->{$b}->{file_name}} keys %{$files_ref}){

    my $columns_ref = $files_ref->{$file}->{columns};

    unless ($columns_ref->{$column1}){
      print YELLOW, "\n\tFile: $files_ref->{$file}->{file_name} does not have the specified column1: $column1\n\n", RESET; 
      next();
    }
    unless ($column2 eq "NA"){
      unless ($columns_ref->{$column2}){
        print YELLOW, "\n\tFile: $files_ref->{$file}->{file_name} does not have the specified column2: $column2\n\n", RESET; 
        next();
      }
    }

    my $file_path = $files_ref->{$file}->{file_path};
    my $column1_pos = $columns_ref->{$column1}->{column_position};

    my $column2_pos;
    unless ($column2 eq "NA"){
      $column2_pos = $columns_ref->{$column2}->{column_position};
    }

    print BLUE, "\n\tProcessing: $files_ref->{$file}->{file_name}", RESET;
    open (FILE, "zcat $file_path |") || die "\nCould not open file: $file_path";

    my $first_line = 1;
    while(<FILE>){
      #Skip the header line
      if ($first_line == 1){
        $first_line = 0;
        next();
      }
      chomp($_);
      my @line = split("\t", $_);
      my $data1 = $line[$column1_pos];
      chomp($data1);

      my $data2;
      unless ($column2 eq "NA"){
        $data2 = $line[$column2_pos];
        chomp($data2);
      }


      #Process INT data type columns
      if ($data_type =~ /int|integer/i){
        if ($data1 =~ /^\d+$/){
          $valid_record_count++;
          print TEMP "$data1\n";
        }
        unless($column2 eq "NA"){
          if ($data2 =~ /^\d+$/){
            $valid_record_count++;
            print TEMP "$data2\n";
          }
        }
      }


      #Process FLOAT data type columns
      if ($data_type =~ /flt|float/i){
        if ($data1 =~ /^\d+$|^\d+\.\d+$/){
          $valid_record_count++;
          print TEMP "$data1\n";
        }
        unless($column2 eq "NA"){
          if ($data2 =~ /^\d+$|^\d+\.\d+$/){
            $valid_record_count++;
            print TEMP "$data2\n";
          }
        }
      }


      #Process STRING data type columns (leave as is except remove NA's)
      if ($data_type =~ /string/i){
        unless ($data1 =~ /^NA$|$column1/){
          $valid_record_count++;

          #Special case for chromosome names
          if (($column1 =~ /chromosome/i) && ($data1 =~ /chrNT|chrc/)){
            $data1 = "chrNT-Hap";
          }
          if ($data1 =~ /^chr(\d)$/){
            $data1 = "chr0$1";
          }

          print TEMP "$data1\n";
        }

        unless($column2 eq "NA"){
          unless ($data2 =~ /^NA$|$column2/){
            $valid_record_count++;

            #Special case for chromosome names
            if (($column2 =~ /chromosome/i) && ($data2 =~ /chrNT|chrc/)){
              $data2 = "chrNT-Hap";
            }
            if ($data2 =~ /^chr(\d)$/){
              $data2 = "chr0$1";
            }

            print TEMP "$data2\n";
          }
        }
      }

    }
    close(FILE);
  }
  close(TEMP);
  print BLUE, "\nFound and wrote a total of $valid_record_count valid records to the temp file\n\n", RESET;
  return($valid_record_count);
}


###########################################################################################################
#Parse a specific data column1 (and column2 if specified) and produce a temp file with this data          #
#- Temp file will be created with only a DIFFERENT column for each set of data                            #
###########################################################################################################
sub parseDataColumn_multi{
  my %args = @_;
  my $files_ref = $args{'-files_ref'};
  my $column1 = $args{'-column1'};
  my $column2 = $args{'-column2'};
  my $data_type = $args{'-data_type'};
  my $temp_file1 = $args{'-temp_file'};

  my $max_line_count = 0;
  my $column_count = 0;
  my $file_count = 0;

  print BLUE, "\n\nProcessing files and parsing data: $column1 to temp file: $temp_file1\n", RESET;

  unless($column2 eq "NA"){
    print BLUE, "\tAlso parsing and joining data: $column2 to temp file: $temp_file1\n", RESET;
  }


  foreach my $file (sort {$files_ref->{$a}->{file_name} cmp $files_ref->{$b}->{file_name}} keys %{$files_ref}){
    $file_count++;
    my $current_line_count_a = 0;
    my $current_line_count_b = 0;

    my $columns_ref = $files_ref->{$file}->{columns};

    unless ($columns_ref->{$column1}){
      print YELLOW, "\n\tFile: $files_ref->{$file}->{file_name} does not have the specified column1: $column1\n\n", RESET; 
      next();
    }
    $column_count++;
    my $temp_file2a = "$temp_file"."$column_count";
    open (TEMP2a, ">$temp_file2a") || die "\nCould not open temp file2a: $temp_file2a\n\n";

    my $temp_file2b;
    unless ($column2 eq "NA"){
      unless ($columns_ref->{$column2}){
        print YELLOW, "\n\tFile: $files_ref->{$file}->{file_name} does not have the specified column2: $column2\n\n", RESET; 
        next();
      }
      $column_count++;
      $temp_file2b = "$temp_file"."$column_count";
      open (TEMP2b, ">$temp_file2b") || die "\nCould not open temp file2b: $temp_file2b\n\n";
    }

    my $file_path = $files_ref->{$file}->{file_path};
    my $column1_pos = $columns_ref->{$column1}->{column_position};
    my $column2_pos;
    unless ($column2 eq "NA"){
      $column2_pos = $columns_ref->{$column2}->{column_position};
    }

    print BLUE, "\n\tProcessing: $files_ref->{$file}->{file_name}", RESET;
    open (FILE, "zcat $file_path |") || die "\nCould not open file: $file_path";

    my $first_line = 1;
    while(<FILE>){
      #Skip the header line
      if ($first_line == 1){
        $first_line = 0;
        next();
      }

      chomp($_);
      my @line = split("\t", $_);
      my $data1 = $line[$column1_pos];
      chomp($data1);

      my $data2;
      unless ($column2 eq "NA"){
        $data2 = $line[$column2_pos];
        chomp($data2);
      }

      #Process FLOAT data type columns
      if ($data_type =~ /flt|float/i){
        if ($data1 =~ /^\d+$|^\d+\.\d+$/){
          $current_line_count_a++;
          print TEMP2a "$data1\n";
        }
        unless($column2 eq "NA"){
          if ($data2 =~ /^\d+$|^\d+\.\d+$/){
            $current_line_count_b++;
            print TEMP2b "$data2\n";
          }
        }
      }else{
        print RED, "\n\nparseDataColumn_multi() only works with --data_type=FLOAT\n\n", RESET;
        exit();
      }

    }
    close(FILE);
    close(TEMP2a);
    unless($column2 eq "NA"){
      close(TEMP2b);
    }

    #Keep track of the maximum lines for any column written to the output file (this file will have have columns of different lengths)
    if ($current_line_count_a > $max_line_count){
      $max_line_count = $current_line_count_a;
    }
    if ($current_line_count_b > $max_line_count){
      $max_line_count = $current_line_count_b;
    }


  }

  #create a string listing all files in order
  my $file_string = '';
  for (my $i = 1; $i <= $column_count; $i++){
    $file_string = "$file_string"."$temp_file"."$i ";
  }

  #Paste all temp files together in order as a single command
  my $cmd_paste = "paste -d \',\' $file_string > $temp_file";
  system($cmd_paste);
  #print YELLOW, "\n\n$cmd_paste\n\n", RESET;

  #Delete all temp file
  my $cmd_rm = "rm -f $file_string";
  system($cmd_rm);
  #print YELLOW, "\n\n$cmd_rm\n\n", RESET;

  print BLUE, "\n\nCreated a multi-column temp file $column_count columns and $max_line_count lines\n\n", RESET;

  my @return = ($column_count, $max_line_count);
  return(\@return);
}


###########################################################################################################
#Parse data files and generate basic statistics                                                           #
#Summarize read pairing and strand pairing statistics                                                     #
###########################################################################################################
sub parseData_stats{
  my %args = @_;
  my $files_ref = $args{'-files_ref'};
  my $results_dir = $args{'-results_dir'};

  print BLUE, "\n\nProcessing files and generating basic summary stats\n", RESET;

  #Lane Counters
  my $total_read_pairs_lane = 0;
  my $same_gene_read_pairs_lane = 0;
  my $ambiguous_read_pairs_lane = 0;

  #Grand Counters
  my $total_read_pairs = 0;
  my $mapped_read_pairs = 0;
  my $same_gene_read_pairs = 0;
  my $ambiguous_read_pairs = 0;
  my $strand_plus_plus = 0;
  my $strand_minus_minus = 0;
  my $strand_plus_minus = 0;
  my $strand_minus_plus = 0;

  #target columns
  my $column1 = 'R1_HitType'; 
  my $column2 = 'R2_HitType'; 
  my $column3 = 'R1_GeneID'; 
  my $column4 = 'R2_GeneID'; 
  my $column5 = 'R1_Strand'; 
  my $column6 = 'R2_Strand'; 

  foreach my $file (sort {$files_ref->{$a}->{file_name} cmp $files_ref->{$b}->{file_name}} keys %{$files_ref}){

    my $columns_ref = $files_ref->{$file}->{columns};
    my $file_path = $files_ref->{$file}->{file_path};

    unless($columns_ref->{$column1}->{column_position} && $columns_ref->{$column2}->{column_position} && $columns_ref->{$column3}->{column_position} && $columns_ref->{$column4}->{column_position} && $columns_ref->{$column5}->{column_position} && $columns_ref->{$column6}->{column_position}){
      print RED, "\nNeccessary column for parseData_stats() is missing in file: $file_path\n\n", RESET;
      exit();
    }

    my $column1_pos = $columns_ref->{$column1}->{column_position};
    my $column2_pos = $columns_ref->{$column2}->{column_position};
    my $column3_pos = $columns_ref->{$column3}->{column_position};
    my $column4_pos = $columns_ref->{$column4}->{column_position};
    my $column5_pos = $columns_ref->{$column5}->{column_position};
    my $column6_pos = $columns_ref->{$column6}->{column_position};

    print BLUE, "\n\tProcessing: $files_ref->{$file}->{file_name}", RESET;
    open (FILE, "zcat $file_path |") || die "\nCould not open file: $file_path";

    my $first_line = 1;
    while(<FILE>){
      #Skip the header line
      if ($first_line == 1){
        $first_line = 0;
        next();
      }

      chomp($_);
      my @line = split("\t", $_);
      my $r1_hit_type = $line[$column1_pos];
      my $r2_hit_type = $line[$column2_pos];
      my $r1_gene_id = $line[$column3_pos];
      my $r2_gene_id = $line[$column4_pos];
      my $r1_strand = $line[$column5_pos];
      my $r2_strand = $line[$column6_pos];

      #print YELLOW, "\n$r1_hit_type\t$r2_hit_type\t$r1_gene_id\t$r2_gene_id\t$r1_strand\t$r2_strand", RESET;

      #Summarize types of read pair mappings
      $total_read_pairs++;
      $total_read_pairs_lane++;

      #Count the subset of total read pairs where both reads are MAPPED as an unambiguous hit
      if ($r1_hit_type eq "Top_Hit" && $r2_hit_type eq "Top_Hit"){
        $mapped_read_pairs++;
        
        #Count the subset of mapped reads where both reads map to the SAME gene
        if ($r1_gene_id == $r2_gene_id){
          $same_gene_read_pairs++;
          $same_gene_read_pairs_lane++;
          

          #Count the read pairs with strand values of: +/+, -/-. +/-, and -/+
          if ($r1_strand eq "+" && $r2_strand eq "+"){
            $strand_plus_plus++;
          }
          if ($r1_strand eq "-" && $r2_strand eq "-"){
            $strand_minus_minus++;
          }
          if ($r1_strand eq "+" && $r2_strand eq "-"){
            $strand_plus_minus++;
          }
          if ($r1_strand eq "-" && $r2_strand eq "+"){
            $strand_minus_plus++;
          }
        }
      }

      #Compare read pairs that are both Top_Hits to the same gene ($same_gene_read_pairs) versus read pairs that both map ambiguously ($ambiguous_read_pairs)
      if ($r1_hit_type eq "Ambiguous" && $r2_hit_type eq "Ambiguous"){
        $ambiguous_read_pairs++;
        $ambiguous_read_pairs_lane++;
      }
    }

    #For each lane summarize the total_read_pairs, same_gene_read_pairs, and ambiguous_read_pairs 
    my $total_read_pairs_lane_p = sprintf("%.2f", (($total_read_pairs_lane/$total_read_pairs_lane)*100));
    my $ambiguous_read_pairs_lane_p = sprintf("%.2f", (($ambiguous_read_pairs_lane/$total_read_pairs_lane)*100));
    my $same_gene_read_pairs_lane_p = sprintf("%.2f", (($same_gene_read_pairs_lane/$total_read_pairs_lane)*100));
    my $read_pairs_combined_lane =  $same_gene_read_pairs_lane + $ambiguous_read_pairs_lane;
    my $ambiguous_read_pairs_lane_pp = sprintf("%.2f", (($ambiguous_read_pairs_lane/$read_pairs_combined_lane)*100));
    my $same_gene_read_pairs_lane_pp = sprintf("%.2f", (($same_gene_read_pairs_lane/$read_pairs_combined_lane)*100));
   
    print BLUE, "\n\t\tTotal read pairs:\t$total_read_pairs_lane\t$total_read_pairs_lane_p%", RESET;
    print BLUE, "\n\t\tSame gene read pairs:\t$same_gene_read_pairs_lane\t$same_gene_read_pairs_lane_p% (% of total)\t$same_gene_read_pairs_lane_pp% (% of mapped read pairs)", RESET;
    print BLUE, "\n\t\tAmbiguous read pairs:\t$ambiguous_read_pairs_lane\t$ambiguous_read_pairs_lane_p% (% of total)\t$ambiguous_read_pairs_lane_pp% (% of mapped read pairs)\n", RESET;

    $total_read_pairs_lane = 0;
    $ambiguous_read_pairs_lane = 0;
    $same_gene_read_pairs_lane = 0;
    close(FILE);
  }

  my $strand_incorrect = $strand_plus_plus + $strand_minus_minus;
  my $strand_correct = $strand_plus_minus + $strand_minus_plus;

  #Calculate percents for each value
  my $mapped_read_pairs_p = sprintf("%.2f", (($mapped_read_pairs/$total_read_pairs)*100));
  my $same_gene_read_pairs_p = sprintf("%.2f", (($same_gene_read_pairs/$mapped_read_pairs)*100));
  my $ambiguous_read_pairs_p =  sprintf("%.2f", (($ambiguous_read_pairs/$total_read_pairs)*100));
  my $strand_plus_plus_p =  sprintf("%.2f", (($strand_plus_plus/$same_gene_read_pairs)*100));
  my $strand_minus_minus_p =  sprintf("%.2f", (($strand_minus_minus/$same_gene_read_pairs)*100));
  my $strand_plus_minus_p =  sprintf("%.2f", (($strand_plus_minus/$same_gene_read_pairs)*100));
  my $strand_minus_plus_p =  sprintf("%.2f", (($strand_minus_plus/$same_gene_read_pairs)*100));
  my $strand_correct_p = sprintf("%.2f", (($strand_correct/$same_gene_read_pairs)*100));
  my $strand_incorrect_p = sprintf("%.2f", (($strand_incorrect/$same_gene_read_pairs)*100));

  my $read_pairs_combined =  $same_gene_read_pairs + $ambiguous_read_pairs;
  my $ambiguous_read_pairs_pp = sprintf("%.2f", (($ambiguous_read_pairs/$read_pairs_combined)*100));
  my $same_gene_read_pairs_pp = sprintf("%.2f", (($same_gene_read_pairs/$read_pairs_combined)*100));
   

  #Print out the resulting stats to a stats file
  my $stats_file = "$results_dir"."$data_name".".txt";
  print BLUE, "\nPrinting results to output file: $stats_file\n\n", RESET;
  open(STATS, ">$stats_file") || die "\nCould not open output stats file: $stats_file\n\n";

  print BLUE, "MAPPING STATS:\n", RESET;
  print BLUE, "Total read pairs with any hit to the transcriptome (one read or both):\t$total_read_pairs\n", RESET;
  print BLUE, "Total of these read pairs where both reads map unambiguously to a single gene:\t$mapped_read_pairs\t$mapped_read_pairs_p%\n", RESET;
  print BLUE, "Total of these read pairs where both reads map to the SAME gene:\t$same_gene_read_pairs\t$same_gene_read_pairs_p% (% of total)\t$same_gene_read_pairs_pp% (% of mapped read pairs)\n", RESET;
  print BLUE, "Total of these read pairs where both reads map ambiguously:\t$ambiguous_read_pairs\t$ambiguous_read_pairs_p% (% of total)\t$ambiguous_read_pairs_pp% (% of mapped read pairs)\n", RESET;
  print BLUE, "\nSTRAND STATS (for SAME gene mapping reads only):\n", RESET;
  print BLUE, "+/+:\t$strand_plus_plus\t$strand_plus_plus_p%\n", RESET;
  print BLUE, "-/-:\t$strand_minus_minus\t$strand_minus_minus_p%\n", RESET;
  print BLUE, "+/-:\t$strand_plus_minus\t$strand_plus_minus_p%\n", RESET;
  print BLUE, "-/+:\t$strand_minus_plus\t$strand_minus_plus_p%\n", RESET;
  print BLUE, "\nCorrect:\t$strand_correct\t$strand_correct_p%\n", RESET;
  print BLUE, "Incorrect:\t$strand_incorrect\t$strand_incorrect_p%\n", RESET;

  print STATS "MAPPING STATS:\n";
  print STATS "Total read pairs with any hit to the transcriptome (one read or both):\t$total_read_pairs\n";
  print STATS "Total of these read pairs where both reads map unambiguously to a single gene:\t$mapped_read_pairs\t$mapped_read_pairs_p%\n";
  print STATS "Total of these read pairs where both reads map to the SAME gene:\t$same_gene_read_pairs\t$same_gene_read_pairs_p% (% of total)\t$same_gene_read_pairs_pp% (% of mapped read pairs)\n";
  print STATS "Total of these read pairs where both reads map ambiguously:\t$ambiguous_read_pairs\t$ambiguous_read_pairs_p% (% of total)\t$ambiguous_read_pairs_pp% (% of mapped read pairs)\n";
  print STATS "\nSTRAND STATS (for SAME gene mapping reads only):\n";
  print STATS "+/+:\t$strand_plus_plus\t$strand_plus_plus_p%\n";
  print STATS "-/-:\t$strand_minus_minus\t$strand_minus_minus_p%\n";
  print STATS "+/-:\t$strand_plus_minus\t$strand_plus_minus_p%\n";
  print STATS "-/+:\t$strand_minus_plus\t$strand_minus_plus_p%\n";
  print STATS "\nCorrect:\t$strand_correct\t$strand_correct_p%\n";
  print STATS "Incorrect:\t$strand_incorrect\t$strand_incorrect_p%\n";

  close(STATS);

  return();
}




