#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Create all neccessary analysis directories needed to run the ALEXA-Seq pipeline
#Start with a root dir and a file containing library and flowcell lane IDs

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
  print GREEN, "\n\nExample: createAnalysisDirs.pl  --analysis_dir=/projects/malachig/alexa_seq/  --project=Neuroblastoma  --library_conf_file=/projects/malachig/alexa_seq/config_files/project_config_files/Neuroblastoma/ALEXA_Seq_Neuroblastoma.conf  --ensembl_version=53\n\n", RESET;
  exit();
}
$analysis_dir = &checkDir('-dir'=>$analysis_dir, '-clear'=>"no");

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
my $repeat = "Repeats";
my $enst = "ENST_v"."$ensembl_version";
my $junction = "Junctions_v"."$ensembl_version";
my $boundary = "Boundaries_v"."$ensembl_version";
my $intron = "Introns_v"."$ensembl_version";
my $intergenic = "Intergenics_v"."$ensembl_version";
my $transcript = "Transcripts_v"."$ensembl_version";
my @types = ($repeat,$enst,$junction,$boundary,$intron,$intergenic,$transcript);

#Create dirs
my $dir = '';
my $dirs_created = 0;

#Create base dirs and project level dirs
$dir = "$analysis_dir"."temp/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."temp/complexity/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."temp/mdust/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."temp/website/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."temp/website/images/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."temp/website/$project/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."temp/website/$project/genes/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."temp/website/$project/images/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."temp/website/$project/data/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."temp/website/$project/ucsc/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."temp/website/$project/ucsc/DE/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."batch_jobs/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."batch_jobs/$project/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."batch_jobs/$project/process/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."batch_jobs/$project/expression/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."batch_jobs/$project/stats/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."batch_jobs/$project/website/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."blast_results/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."fasta_seq_data/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."figures_and_stats/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."figures_and_stats/DE/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."figures_and_stats/SI/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."figures_and_stats/DE/$project/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."figures_and_stats/SI/$project/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."logs/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."logs/website/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."logs/website/$project/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."logs/website/$project/genes/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."raw_seq_data/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$analysis_dir"."read_records/"; $dirs_created += &make_dir('-dir'=>$dir);

#Create comparison specific dirs
foreach my $lib (keys %libs){
  $dir = "$analysis_dir"."temp/website/$project/ucsc/$lib/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."temp/website/$project/ucsc/$lib/combined/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."temp/website/$project/ucsc/DE/$lib/"; $dirs_created += &make_dir('-dir'=>$dir);
  foreach my $comp (keys %comps){
    $dir = "$analysis_dir"."temp/website/$project/ucsc/DE/$lib/$comp/"; $dirs_created += &make_dir('-dir'=>$dir);
  }
}

#Create type specific dirs
foreach my $type (@types){
  $dir = "$analysis_dir"."figures_and_stats/DE/$project/$type/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."figures_and_stats/SI/$project/$type/"; $dirs_created += &make_dir('-dir'=>$dir);
}


#Create dirs specific to each library 
foreach my $lib (keys %libs){

  #Create base dirs for this library
  $dir = "$analysis_dir"."batch_jobs/$project/$lib/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."batch_jobs/$project/$lib/blast_vs_repeats/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."batch_jobs/$project/$lib/blast_vs_transcripts/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."batch_jobs/$project/$lib/blast_vs_exonJunctions/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."batch_jobs/$project/$lib/blast_vs_exonBoundaries/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."batch_jobs/$project/$lib/blast_vs_introns/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."batch_jobs/$project/$lib/blast_vs_intergenics/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."raw_seq_data/$lib/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."read_records/$lib/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."read_records/$lib/Summary/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."logs/$lib/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."logs/$lib/generateExpressionValues/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."logs/$lib/partitionMapFiles/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."fasta_seq_data/$lib/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."blast_results/$lib/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."figures_and_stats/$lib/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."figures_and_stats/$lib/LibraryQuality/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."figures_and_stats/$lib/Expression_v"."$ensembl_version/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."figures_and_stats/$lib/Generic/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."figures_and_stats/$lib/ENST_v"."$ensembl_version/"; $dirs_created += &make_dir('-dir'=>$dir);
  $dir = "$analysis_dir"."figures_and_stats/$lib/temp/"; $dirs_created += &make_dir('-dir'=>$dir);

  #Create type specific dirs for this library
  foreach my $type (@types){
    $dir = "$analysis_dir"."read_records/$lib/$type/"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."read_records/$lib/$type/temp/"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."read_records/$lib/$type/Summary/"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."read_records/$lib/$type/Summary/results"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."read_records/$lib/$type/Summary/temp"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."logs/$lib/generateExpressionValues/$type/"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."figures_and_stats/DE/$project/$type/"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."figures_and_stats/SI/$project/$type/"; $dirs_created += &make_dir('-dir'=>$dir);
  }

  #Create flowcell lane specific dirs for this library
  my @lanes = @{$libs{$lib}{flowcell_lanes}};
  foreach my $lane (@lanes){
    $dir = "$analysis_dir"."fasta_seq_data/$lib/$lane"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."fasta_seq_data/$lib/$lane/fasta_blocks/"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."logs/$lib/$lane/"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."blast_results/$lib/$lane/"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."blast_results/$lib/$lane/repeats/"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."blast_results/$lib/$lane/transcripts/"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."blast_results/$lib/$lane/exonJunctions/"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."blast_results/$lib/$lane/exonBoundaries/"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."blast_results/$lib/$lane/introns/"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."blast_results/$lib/$lane/intergenics/"; $dirs_created += &make_dir('-dir'=>$dir);
    $dir = "$analysis_dir"."blast_results/$lib/$lane/temp/"; $dirs_created += &make_dir('-dir'=>$dir);
  }
}

print BLUE, "\n\nCreated $dirs_created directories\n\n", RESET;


exit();

sub make_dir{
  my %args = @_;
  my $dir = $args{'-dir'};
  unless (-e $dir && -d $dir){
    mkdir($dir);
    print MAGENTA, "\nCreating $dir", RESET;
    return(1);
  }
  #print YELLOW, "\nAlready exists: mkdir $dir", RESET;
  return(0);
}

