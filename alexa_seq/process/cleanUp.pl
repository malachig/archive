#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Once an ALEXA-Seq analysis is finalized, use this script to clean-up all temp files that are no longer needed

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

my $analysis_dir = '';
my $project = '';
my $library_conf_file = '';
my $ensembl_version = '';

GetOptions ('analysis_dir=s'=>\$analysis_dir, 'project=s'=>\$project, 'library_conf_file=s'=>\$library_conf_file, 'ensembl_version=i'=>\$ensembl_version);

if ($analysis_dir && $project && $library_conf_file && $ensembl_version){
  #Parameters ok
}else{
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  print GREEN, "\n\nExample: cleanUp.pl  --analysis_dir=/projects/alexa2/alexa_seq/  --project=AllenBrain  --library_conf_file=/projects/alexa2/alexa_seq/config_files/project_config_files/AllenBrain/ALEXA_Seq_AllenBrain.conf  --ensembl_version=57\n\n", RESET;
  exit();
}
$analysis_dir = &checkDir('-dir'=>$analysis_dir, '-clear'=>"no");

#Make sure the user is actually ready to delete these temp files
print BLUE, "\n\nThis script is going to clean up many temp files.  DO NOT this unless you are sure that the ALEXA-Seq analysis has completed successfully!!!\n\n", RESET;
print YELLOW, "Proceed? (yes/no): ", RESET;
my $answer = <>;
chomp($answer);
if ($answer =~ /yes/i){
  print BLUE, "\n\nProceeding then ...\n\n", RESET;
}else{
  print RED, "\n\nAborting then ...\n\n", RESET;
  exit();
}


#Get library IDs and flowcell/Lane names
open(LIBS, "$library_conf_file") || die "\nCould not open lib file: $library_conf_file\n\n";
my %libs;
my %comps;
while(<LIBS>){
  #Skip empty lines
  unless ($_ =~ /\w+|\d+/){
    next();
  }
  chomp($_);
  my @data = split(/ +/, $_);
  if ($_ =~ /^LANE/){
    my $lib = $data[1];
    my $flowcell_lane = "$data[2]"."_Lane"."$data[3]";
    if ($libs{$lib}){
      push(@{$libs{$lib}{flowcell_lanes}}, $flowcell_lane);
    }else{
      my @tmp;
      push(@tmp, $flowcell_lane);
      $libs{$lib}{flowcell_lanes} = \@tmp;
    }
  }
  if ($_ =~ /^COMPARISON/){
    my $comp_name = $data[5];
    my $comp = "$data[1]"."_vs_"."$data[3]";
    $comps{$comp}{libA} = $data[1];
    $comps{$comp}{libB} = $data[3];
  }
}
close(LIBS);

#Define types
my %types1;
$types1{1}{name} = "Repeats";
$types1{2}{name} = "ENST_v"."$ensembl_version";
$types1{3}{name} = "Junctions_v"."$ensembl_version";
$types1{4}{name} = "Boundaries_v"."$ensembl_version";
$types1{5}{name} = "Introns_v"."$ensembl_version";
$types1{6}{name} = "Intergenics_v"."$ensembl_version";
$types1{7}{name} = "Transcripts_v"."$ensembl_version";

my %types2;
$types2{1}{name} = "repeats";
$types2{1}{abr1} = "R";
$types2{1}{abr2} = "R";
$types2{1}{abr3} = "r";
$types2{2}{name} = "transcripts";
$types2{2}{abr1} = "T";
$types2{2}{abr2} = "GE";
$types2{2}{abr3} = "t";
$types2{3}{name} = "exonJunctions";
$types2{3}{abr1} = "EJ";
$types2{3}{abr2} = "J";
$types2{3}{abr3} = "j";
$types2{4}{name} = "exonBoundaries";
$types2{4}{abr1} = "EB";
$types2{4}{abr2} = "B";
$types2{4}{abr3} = "b";
$types2{5}{name} = "introns";
$types2{5}{abr1} = "I";
$types2{5}{abr2} = "I";
$types2{5}{abr3} = "i";
$types2{6}{name} = "intergenics";
$types2{6}{abr1} = "IG";
$types2{6}{abr2} = "IG";
$types2{6}{abr3} = "ig";

my $dir = "";
my $dels = 0;
my $ls_cmd = "";
my $rm_cmd = "";

#1.) Temp files and dirs in 'batch_jobs'
my $batch_dir = "$analysis_dir"."batch_jobs/$project/";
print BLUE, "\n\nCleaning temp files and dirs in $batch_dir", RESET;
foreach my $lib (keys %libs){
  foreach my $t (sort {$a <=> $b} keys %types2){
    my $type_dir_name = $types2{$t}{name};
    my $abr1 = $types2{$t}{abr1};
    my $dir = "$batch_dir"."$lib/blast_vs_"."$type_dir_name/B$abr1"."_"."$lib";
    #$ls_cmd = "ls -d $dir";
    #system($ls_cmd);
    if (-e $dir && -d $dir){
      $rm_cmd = "rm -fr $dir";
      print MAGENTA, "\n\t$rm_cmd", RESET;
      system($rm_cmd);
      $dels++;
    }
  }
}
foreach my $t (sort {$a <=> $b} keys %types2){
  my $abr2 = $types2{$t}{abr2};
  my $dir = "$batch_dir"."expression/GEV_"."$abr2";
  #$ls_cmd = "ls -d $dir";
  #system($ls_cmd);
  if (-e $dir && -d $dir){
    $rm_cmd = "rm -fr $dir";
    print MAGENTA, "\n\t$rm_cmd", RESET;
    system($rm_cmd);
    $dels++;
  }
}
$dir = "$batch_dir"."expression/pmf/";
if (-e $dir && -d $dir){
  $rm_cmd = "rm -fr $dir";
  print MAGENTA, "\n\t$rm_cmd", RESET;
  system($rm_cmd);
  $dels++;
}
foreach my $t (sort {$a <=> $b} keys %types2){
  my $abr3 = $types2{$t}{abr3};
  my $dir = "$batch_dir"."process/p$abr3";
  #$ls_cmd = "ls -d $dir";
  #system($ls_cmd);
  if (-e $dir && -d $dir){
    $rm_cmd = "rm -fr $dir";
    print MAGENTA, "\n\t$rm_cmd", RESET;
    system($rm_cmd);
    $dels++;
  }
}
$dir = "$batch_dir"."website/GGH/";
if (-e $dir && -d $dir){
  $rm_cmd = "rm -fr $dir";
  print MAGENTA, "\n\t$rm_cmd", RESET;
  system($rm_cmd);
  $dels++;
}

#2.) Temp files and dirs in 'fasta_seq_data'
my $fasta_dir = "$analysis_dir"."fasta_seq_data/";
print BLUE, "\n\nCleaning temp files and dirs in $fasta_dir", RESET;
foreach my $lib (keys %libs){
  my @lanes = @{$libs{$lib}{flowcell_lanes}};
  foreach my $lane (@lanes){
    my $fgz = "$fasta_dir"."$lib/$lane/*.fa.gz";
    if (-e $fgz){
      $rm_cmd = "rm -f $fgz";
      print MAGENTA, "\n\t$rm_cmd", RESET;
      system($rm_cmd);
      $dels++;
    }
    $dir = "$fasta_dir"."$lib/$lane/fasta_blocks";
    if (-e $dir && -d $dir){
      $rm_cmd = "rm -fr $dir";
      print MAGENTA, "\n\t$rm_cmd", RESET;
      system($rm_cmd);
      $dels++;
      &make_dir('-dir'=>$dir);
    }
  }
}

#3.) Temp files and dirs in 'raw_seq_data'
my $raw_seq_dir = "$analysis_dir"."raw_seq_data/";
print BLUE, "\n\nCleaning temp files and dirs in $raw_seq_dir", RESET;
foreach my $lib (keys %libs){
  $dir = "$raw_seq_dir"."$lib";
  if (-e $dir && -d $dir){
    $rm_cmd = "rm -fr $dir/*";
    print MAGENTA, "\n\t$rm_cmd", RESET;
    system($rm_cmd);
    $dels++;
  }
}


#4.) Temp files and dirs in 'read_records'
my $read_records_dir = "$analysis_dir"."read_records/";
print BLUE, "\n\nCleaning temp files and dirs in $read_records_dir", RESET;
foreach my $lib (keys %libs){
  foreach my $t (sort {$a <=> $b} keys %types1){
    my $type_dir_name = $types1{$t}{name};
    $dir = "$read_records_dir"."$lib/"."$type_dir_name/part";
    if (-e $dir && -d $dir){
      $rm_cmd = "rm -fr $dir";
      print MAGENTA, "\n\t$rm_cmd", RESET;
      system($rm_cmd);
      $dels++;
    }
  }
}



print BLUE "\n\nPerformed $dels total deletion commands\n\n", RESET;


exit();


########################################################################################
#Make a dir                                                                            #
########################################################################################
sub make_dir{
  my %args = @_;
  my $dir = $args{'-dir'};
  unless (-e $dir && -d $dir){
    mkdir($dir);
    print MAGENTA, "\n\tCreating $dir", RESET;
    return(1);
  }
  #print YELLOW, "\nAlready exists: mkdir $dir", RESET;
  return(0);
}

