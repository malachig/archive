#!/usr/bin/perl -w
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to determine the library size of each library and store this information in a simple text file
#'Library Size' is the total number of unambiguously mapping reads (bases actually) for a library
#It is calculated by counting reads with one of the following classes: ENST_U INTRON_U INTERGENIC_U NOVEL_JUNCTION_U NOVEL_BOUNDARY_U
#These values are determined for each lane of each library and then a grand total is calculated

#The output file will contain a simple matrix
#Columns: LANE, ENST_Reads, INTRON_Reads, INTERGENIC_Reads, NOVEL_JUNCTION_Reads, NOVEL_BOUNDARY_Reads, REPEAT_Reads, UNAMBIGUOUS_Reads, AMBIGUOUS_Reads, TOTAL_Reads, LIBRARY_SIZE_Reads, ENST_Bases, INTRON_Bases, INTERGENIC_Reads, NOVEL_JUNCTION_Reads, NOVEL_BOUNDARY_Reads, REPEAT_Bases, UNAMBIGUOUS_Bases, AMBIGUOUS_Bases, TOTAL_Bases, LIBRARY_SIZE_Bases
#Rows: Flowcell_Lane1, Flowcell_Lane2, ..., TOTAL

#The value that will be used to normalize expression values to account for library sizes is: Column->LIBRARY_SIZE_Bases  AND   Row->TOTAL

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use BerkeleyDB;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $analysis_dir = '';
my $library_id = '';

GetOptions ('library_id=s'=>\$library_id,'analysis_dir=s'=>\$analysis_dir);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the analysis dir using:  --analysis_dir", RESET;
print GREEN, "\n\tSpecify the library ID using: --library_id", RESET;
print GREEN, "\n\nExample: calculateLibrarySize.pl  --analysis_dir=/projects/malachig/alexa_seq/  --library_id=HS1441\n\n", RESET;

unless ($analysis_dir && $library_id){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Form and check all neccessary file paths and directories
$analysis_dir = &checkDir('-dir'=>$analysis_dir, '-clear'=>"no");
my $read_records_dir = "$analysis_dir"."read_records/"."$library_id";
$read_records_dir = &checkDir('-dir'=>$read_records_dir, '-clear'=>"no");
my $summary_dir = "$read_records_dir"."Summary";
$summary_dir = &checkDir('-dir'=>$summary_dir, '-clear'=>"no");
my $outfile = "$summary_dir"."MappedReads.txt";


my $lib_size_class_list = "ENST_U INTRON_U INTERGENIC_U NOVEL_JUNCTION_U NOVEL_BOUNDARY_U";
my $u_class_list = "ENST_U INTRON_U INTERGENIC_U NOVEL_JUNCTION_U NOVEL_BOUNDARY_U Repeat_U";
my $a_class_list = "ENST_A INTRON_A INTERGENIC_A NOVEL_JUNCTION_A NOVEL_BOUNDARY_A Repeat_A";
my @summary_list = qw(ENST_U INTRON_U INTERGENIC_U NOVEL_JUNCTION_U NOVEL_BOUNDARY_U Repeat_U UNAMBIGUOUS AMBIGUOUS TOTAL LIBRARY_SIZE);


#Get files from this directory
print BLUE, "\n\nSearching $read_records_dir for read records files", RESET;

my %files;
opendir(DIRHANDLE, "$read_records_dir") || die "\nCannot open directory: $read_records_dir\n\n";
my @test_files = readdir(DIRHANDLE);
my $file_count = 0;

foreach my $test_file (sort @test_files){
  my $file_path = "$read_records_dir"."$test_file";
  my $first_read_pair_id;
  my $test = 1;

  #Skip directories within the specified directory
  if (-e $file_path && -d $file_path){
    #print YELLOW, "\n\t$file_path  is a directory - skipping", RESET;
    next();
  }

  #If the results file is compressed uncompress it
  unless ($file_path =~ /(.*)\.gz$/){
    print RED, "\nFound an uncompressed file: $file_path - make sure no files in progress are in this directory!!\n\n", RESET;
    exit();
  }
  $file_count++;
  print BLUE, "\n\t$file_path was added to the list of files to be processed", RESET;
  $files{$file_count}{name} = $test_file;
  $files{$file_count}{path} = $file_path;
}

#Get count for each lane
print BLUE, "\n\nProcessing these files and counting read mapping classes for each lanes", RESET;
foreach my $fc (sort {$a <=> $b} keys %files){
  my %counts;
  my $col_count = 0;
  foreach my $col (@summary_list){
    $col_count++;
    $counts{$col}{order} = $col_count;
    $counts{$col}{reads} = 0;
    $counts{$col}{bases} = 0;
  }
  #print Dumper %counts;

  my $file_name = $files{$fc}{name};
  my $file_path = $files{$fc}{path};
  my %columns;
  print BLUE, "\n\nProcessing: $file_name\n", RESET;

  open (READ, "zcat $file_path |") || die "\nCould not open file: $file_path\n\n";
  my $header = 1;
  my $counter = 0;
  while(<READ>){
    $counter++;
    chomp($_);
    my @line = split("\t", $_);

    if ($header == 1){
      my $column_count = 0;
      foreach my $column (@line){
        $columns{$column}{column_pos} = $column_count;
        $column_count++;
      }
      $header = 0;
      next();
    }

    if ($counter == 100000){
      $counter = 0;
      $| = 1; print BLUE, ".", RESET; $| = 0;
    }
    my $r1_id = $line[$columns{'Read1_ID'}{column_pos}];
    my $r2_id = $line[$columns{'Read2_ID'}{column_pos}];
    my $r1_status = $line[$columns{'Read1_Status'}{column_pos}];
    my $r2_status = $line[$columns{'Read2_Status'}{column_pos}];
    my $r1_seq_length = length($line[$columns{'Read1_Seq'}{column_pos}]);
    my $r2_seq_length = length($line[$columns{'Read2_seq'}{column_pos}]);

    #READ1
    if ($a_class_list =~ /($r1_status)/){
      #If read was mapped ambiguously to any class, add a count to the Ambiguous total
      $counts{AMBIGUOUS}{reads}++;
      $counts{AMBIGUOUS}{bases}+=$r1_seq_length;

      #Also add a count for the total mapped
      $counts{TOTAL}{reads}++;
      $counts{TOTAL}{bases}+=$r1_seq_length;
    }elsif ($u_class_list =~ /($r1_status)/){
      #If read was mapped unambiguously to any class, add a count to the Unambiguous total
      $counts{UNAMBIGUOUS}{reads}++;
      $counts{UNAMBIGUOUS}{bases}+=$r1_seq_length;

      #Also add a count for this specific unambiguous class
      $counts{$r1_status}{reads}++;
      $counts{$r1_status}{bases}+=$r1_seq_length;

      #Also add a count for the total mapped
      $counts{TOTAL}{reads}++;
      $counts{TOTAL}{bases}+=$r1_seq_length;
    }
    if ($lib_size_class_list =~ /($r1_status)/){
      #If read was mapped to a class counting toward the library size value, add it to the library size total
      $counts{LIBRARY_SIZE}{reads}++;
      $counts{LIBRARY_SIZE}{bases}+=$r1_seq_length;
    }

    #READ2
    if ($a_class_list =~ /($r2_status)/){
      #If read was mapped ambiguously to any class, add a count to the Ambiguous total
      $counts{AMBIGUOUS}{reads}++;
      $counts{AMBIGUOUS}{bases}+=$r2_seq_length;

      #Also add a count for the total mapped
      $counts{TOTAL}{reads}++;
      $counts{TOTAL}{bases}+=$r2_seq_length;
    }elsif ($u_class_list =~ /($r2_status)/){
      #If read was mapped unambiguously to any class, add a count to the Unambiguous total
      $counts{UNAMBIGUOUS}{reads}++;
      $counts{UNAMBIGUOUS}{bases}+=$r2_seq_length;

      #Also add a count for this specific unambiguous class
      $counts{$r2_status}{reads}++;
      $counts{$r2_status}{bases}+=$r2_seq_length;

      #Also add a count for the total mapped
      $counts{TOTAL}{reads}++;
      $counts{TOTAL}{bases}+=$r2_seq_length;
    }
    if ($lib_size_class_list =~ /($r2_status)/){
      #If read was mapped to a class counting toward the library size value, add it to the library size total
      $counts{LIBRARY_SIZE}{reads}++;
      $counts{LIBRARY_SIZE}{bases}+=$r2_seq_length;
    }
  }
  close(READ);

  #Store the counts for this lane
  #print Dumper %counts;
  $files{$fc}{counts} = \%counts;
}
  

#Now calculate the grand totals for all lanes
print BLUE, "\n\nCalculating GRAND totals", RESET;

my %total_counts;
my $col_count = 0;
foreach my $col (@summary_list){
  $col_count++;
  $total_counts{$col}{order} = $col_count;
  $total_counts{$col}{reads} = 0;
  $total_counts{$col}{bases} = 0;
}
foreach my $fc (sort {$a <=> $b} keys %files){
  my %counts = %{$files{$fc}{counts}};
  foreach my $count (keys %counts){
    $total_counts{$count}{order} = $counts{$count}{order};
    $total_counts{$count}{reads} += $counts{$count}{reads};
    $total_counts{$count}{bases} += $counts{$count}{bases};
  }
}

#Now print out the result to the output file:
open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
print OUT "LANE\tENST_Reads\tINTRON_Reads\tINTERGENIC_Reads\tNOVEL_JUNCTION_Reads\tNOVEL_BOUNDARY_Reads\tREPEAT_Reads\tUNAMBIGUOUS_Reads\tAMBIGUOUS_Reads\tTOTAL_Reads\tLIBRARY_SIZE_Reads\tENST_Bases\tINTRON_Bases\tINTERGENIC_Bases\tNOVEL_JUNCTION_Bases\tNOVEL_BOUNDARY_Bases\tREPEAT_Bases\tUNAMBIGUOUS_Bases\tAMBIGUOUS_Bases\tTOTAL_Bases\tLIBRARY_SIZE_Bases\n";
foreach my $fc (sort {$files{$a}->{name} cmp $files{$b}{name}} keys %files){
  my $file_name = $files{$fc}{name};
  my $name = $file_name;
  if ($file_name =~ /(.*)\.txt\.gz/){
    $name = $1;
  }
  print OUT "$name";

  my %counts = %{$files{$fc}{counts}};
  my $string = '';
  foreach my $count (sort {$counts{$a}{order} <=> $counts{$b}{order}} keys %counts){
    $string .= "\t$counts{$count}{reads}";
  }
  foreach my $count (sort {$counts{$a}{order} <=> $counts{$b}{order}} keys %counts){
    $string .= "\t$counts{$count}{bases}";
  }
  print OUT "$string\n";
}

print OUT "TOTAL";
my $string = '';
foreach my $count (sort {$total_counts{$a}{order} <=> $total_counts{$b}{order}} keys %total_counts){
  $string .= "\t$total_counts{$count}{reads}";
}
foreach my $count (sort {$total_counts{$a}{order} <=> $total_counts{$b}{order}} keys %total_counts){
  $string .= "\t$total_counts{$count}{bases}";
}
print OUT "$string\n";

print "\n\n";

exit();










