#!/usr/bin/perl -w
#Written by Malachi Griffith and Obi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to create expression matrix files.
#For each feature type, two files will be created. One with the normalized expression value and another with either 1 or 0 indicating whether this is above background.

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;
use Tie::File;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $project_name = '';
my $analysis_dir = '';
my $ensembl_version = '';
my $web_dir = '';

GetOptions ('project_name=s'=>\$project_name, 'analysis_dir=s'=>\$analysis_dir, 'ensembl_version=i'=>\$ensembl_version, 'web_dir=s'=>\$web_dir);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the main analysis directory using: --analysis_dir", RESET;
print GREEN, "\n\tSpecify the project name for this project using:  --project_name", RESET;
print GREEN, "\n\tSpecify the ensembl version used for this analysis using: --ensembl_version", RESET;
print GREEN, "\n\nExample:  createMatrixFiles.pl  --project_name=FL_Trans  --analysis_dir=/projects/malachig/alexa_seq/analysis  --ensembl_version=53  --web_dir=/projects/malachig/alexa_seq/www/htdocs/alexa-seq/FL_Trans/data\n\n", RESET;


unless ($project_name && $analysis_dir && $ensembl_version && $web_dir){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}
$analysis_dir = &checkDir('-dir'=>$analysis_dir, '-clear'=>"no");
$web_dir = &checkDir('-dir'=>$web_dir, '-clear'=>"no");


#Define feature types for expression data matrix files
my %feature_types;
$feature_types{"Boundaries"}{"Boundary"}="BoundaryExpression_v"."$ensembl_version";
$feature_types{"Boundaries"}{"KnownBoundary"}="KnownBoundaryExpression_v"."$ensembl_version";
$feature_types{"Boundaries"}{"NovelBoundary"}="NovelBoundaryExpression_v"."$ensembl_version";
$feature_types{"ENST"}{"ExonRegion"}="ExonRegionExpression_v"."$ensembl_version";
$feature_types{"ENST"}{"GeneExpression"}="GeneExpression_v"."$ensembl_version";
$feature_types{"Intergenics"}{"Intergenic"}="IntergenicExpression_v"."$ensembl_version";
$feature_types{"Intergenics"}{"ActiveIntergenicRegion"}="ActiveIntergenicRegionExpression_v"."$ensembl_version";
$feature_types{"Intergenics"}{"SilentIntergenicRegion"}="SilentIntergenicRegionExpression_v"."$ensembl_version";
$feature_types{"Introns"}{"Intron"}="IntronExpression_v"."$ensembl_version";
$feature_types{"Introns"}{"ActiveIntronRegion"}="ActiveIntronRegionExpression_v"."$ensembl_version";
$feature_types{"Introns"}{"SilentIntronRegion"}="SilentIntronRegionExpression_v"."$ensembl_version";
$feature_types{"Junctions"}{"Junction"}="JunctionExpression_v"."$ensembl_version";
$feature_types{"Junctions"}{"KnownJunction"}="KnownJunctionExpression_v"."$ensembl_version";
$feature_types{"Junctions"}{"NovelJunction"}="NovelJunctionExpression_v"."$ensembl_version";
$feature_types{"Transcripts"}{"Transcript"}="TranscriptExpression_v"."$ensembl_version";

#Get the library IDs for this project
my $lib_names_file = "$analysis_dir/batch_jobs/$project_name/$project_name"."_Lib_Names.txt";
my %lib_names;
open(LIBS, "$lib_names_file") || die "\n\nCould not open lib names file: $lib_names_file\n\n";
while(<LIBS>){
  chomp($_);
  my @line = split("\t", $_);
  $lib_names{$line[0]}{name} = $line[1];
}
close(LIBS);

my $dirs_created = 0;
#Make target dir for matrix files that website needs access to
my $matrix_dir = "$web_dir"."matrix/"; $dirs_created += &make_dir('-dir'=>$matrix_dir);
$matrix_dir = &checkDir('-dir'=>$matrix_dir, '-clear'=>"yes", '-recursive'=>"yes");

my @lib_list;
foreach my $lib (sort{$a cmp $b} keys %lib_names){
  push (@lib_list, $lib); #Important that this list is sorted the same way as when data are being printed to matrix files
}

#Create matrix files for expression data (feature x library) for each feature type (one feature type at a time)
foreach my $type (sort {$a cmp $b} keys %feature_types){
  foreach my $subtype (sort {$a cmp $b} keys %{$feature_types{$type}}){
    $| = 1; print MAGENTA, "\n\nCreating matrix file for $subtype", RESET; $| = 0; 

    #Get the header/column info for each file of each library for the current type
    $| = 1; print BLUE, "\n\tExamine the input files", RESET; $| = 0;
    my %files;
    my $master_line_count = 0;
    my $id_name;
    foreach my $lib (sort {$a cmp $b} keys %lib_names){
      my $lib_name = $lib_names{$lib}{name};
      my $expression_file = "$analysis_dir"."read_records/$lib/$type"."_v$ensembl_version/Summary/$lib"."_$feature_types{$type}{$subtype}".".txt";

      #Get the columns for each file and store a file handle
      #my $fh = IO::File->new("$expression_file") || die "\nCould not open output file handle: $expression_file\n\n";
      my @fa;
      tie @fa, 'Tie::File', "$expression_file", mode => O_RDONLY;

      #Determine data columns of interest
      #my %col_names; my $header = <$fh>; chomp $header; my @cols = split("\t",$header); my $i=0; foreach my $col (@cols){$col_names{$col}=$i; $i++;}
      my %col_names; my $header = $fa[0]; chomp $header; my @cols = split("\t",$header); my $i=0; foreach my $col (@cols){$col_names{$col}=$i; $i++;}
      $id_name=$cols[0];
      my $line_count = $#fa;

      #$files{$lib}{fh} = $fh;
      $files{$lib}{fn} = $expression_file;
      $files{$lib}{fa} = \@fa;
      $files{$lib}{cols} = \%col_names;
      $files{$lib}{line_count} = $line_count;
      $master_line_count = $line_count;
      $| = 1; print YELLOW, "\n\t\t$files{$lib}{fn}\t$files{$lib}{line_count}", RESET; $| = 0;
    }

    #Make sure the line count was consistent for all libraries
    $| = 1; print BLUE, "\n\tCheck for consistency of record counts", RESET;
    foreach my $lib (sort {$a cmp $b} keys %files){
      my $line_count = $files{$lib}{line_count};
      unless ($line_count == $master_line_count){
        print RED, "\n\nLine count of library expression files was not consistent!! - Aborting\n\n", RESET;
        exit();
      }
    }

    #Open the output matrix data files and write their headers
    my $matrix_file = "$matrix_dir"."Matrix_$feature_types{$type}{$subtype}".".txt";
    my $expressed_file = "$matrix_dir"."Expressed_$feature_types{$type}{$subtype}".".txt";

    open(MATRIX, ">$matrix_file") or die "can't open $matrix_file for write\n";
    open(EXPRESSED, ">$expressed_file") or die "can't open $expressed_file for write\n";

    print MATRIX "$id_name\tFID\tSeq_Name\tEnsEMBL_Gene_ID\t", join ("\t", @lib_list), "\n";
    print EXPRESSED "$id_name\tFID\tSeq_Name\tEnsEMBL_Gene_ID\t", join ("\t", @lib_list), "\n";


    #Write the matrix file
    my $c = 0;
    $| = 1; print BLUE, "\n\tWrite expression data to output matrix files ", RESET; $| = 0;
    for (my $i = 1; $i <= $master_line_count; $i++){
      $c++;
      if ($c == 100000){
        $| = 1; print BLUE, ".", RESET; $| = 0;
        $c = 0;
      }
      #First gather the data from line $i of each library file
      my %matrix_data;
      my %fids;
      my $master_fid;
      foreach my $lib (sort {$a cmp $b} keys %files){
        my %col_names = %{$files{$lib}{cols}};
        my $fa_ref = $files{$lib}{fa};
        my $line = $fa_ref->[$i];
        chomp $line;
        my @data = split ("\t", $line);
        my $fid=$data[$col_names{'FID'}];
        $master_fid=$fid;
        $fids{$fid}++;
        $matrix_data{$fid}{'ID'}=$data[0]; #Always grab primary ID from first column
        $matrix_data{$fid}{'Seq_Name'}=$data[$col_names{'Seq_Name'}];
        #For intergenics features, a single ENSG is not possible, instead concat the two flanking IDs or set to NA depending on what is available
        if ($subtype eq "Intergenic"){
          $matrix_data{$fid}{'EnsEMBL_Gene_ID'}="$data[$col_names{'Upstream_Gene_ID'}]"."_"."$data[$col_names{'Downstream_Gene_ID'}]";
        }elsif($subtype eq "ActiveIntergenicRegion" || $subtype eq "SilentIntergenicRegion"){
          $matrix_data{$fid}{'EnsEMBL_Gene_ID'}="NA";
        }else{
          $matrix_data{$fid}{'EnsEMBL_Gene_ID'}=$data[$col_names{'EnsEMBL_Gene_ID'}];
        }
        $matrix_data{$fid}{'Average_Coverage_NORM1'}{$lib}=$data[$col_names{'Average_Coverage_NORM1'}];
        $matrix_data{$fid}{'Expressed'}{$lib}=$data[$col_names{'Expressed'}];
      }

      #Make sure the same feature was being processed from each library file...
      my $lib_count = keys %files;
      my $fid_count = $fids{$master_fid};
      unless ($lib_count == $fid_count){
        print RED, "\n\nDid not find the expected number of FIDS ($lib_count libraries, but $fid_count FIDs)\n\n", RESET;
        exit();
      }

      #Now print out the matrix data to the output file
      foreach my $fid (sort{$a cmp $b} keys %matrix_data){
        print MATRIX "$matrix_data{$fid}{'ID'}\t$fid\t$matrix_data{$fid}{'Seq_Name'}\t$matrix_data{$fid}{'EnsEMBL_Gene_ID'}\t";
        print EXPRESSED "$matrix_data{$fid}{'ID'}\t$fid\t$matrix_data{$fid}{'Seq_Name'}\t$matrix_data{$fid}{'EnsEMBL_Gene_ID'}\t";
        my @data;
        foreach my $lib (sort{$a cmp $b} keys %{$matrix_data{$fid}{'Average_Coverage_NORM1'}}){
          push (@data, $matrix_data{$fid}{'Average_Coverage_NORM1'}{$lib});
        }
        my @data2;
        foreach my $lib (sort{$a cmp $b} keys %{$matrix_data{$fid}{'Expressed'}}){
          push (@data2, $matrix_data{$fid}{'Expressed'}{$lib});
        }
        print MATRIX join("\t", @data), "\n";
        print EXPRESSED join("\t", @data2), "\n";
      } 
    }

    #Close the output file handles
    close(MATRIX);
    close(EXPRESSED);

    #Close the file handles
    $| = 1; print BLUE, "\n\tUntie file array handles", RESET; $| = 0;
    foreach my $lib (sort {$a cmp $b} keys %files){
      #my $fh = $files{$lib}{fh};
      #$fh->close;
      my $fa_ref = $files{$lib}{fa};
      untie @{$fa_ref};
    }

    #gzip matrix file for website use
    $| = 1; print BLUE, "\n\tCompress the results files", RESET; $| = 0;
    my $gz_cmd = "gzip -f --best $matrix_file";
    system($gz_cmd);
    my $gz_cmd2 = "gzip -f --best $expressed_file";
    system($gz_cmd2);


  }
}

my $mem_message = &memoryUsage();
print YELLOW, "\n\n$mem_message\n\n", RESET;


exit();


####################################################################################################################
#Make dir                                                                                                          #
####################################################################################################################

sub make_dir{
  my %args = @_;
  my $dir = $args{'-dir'};
  unless (-e $dir && -d $dir){
    mkdir($dir);
    print MAGENTA, "\nCreating $dir", RESET;
    return(1);
  }
  print YELLOW, "\nAlready exists: mkdir $dir", RESET;
  return(0);
}

