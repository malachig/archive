#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#This script acts as a wrapper for the creation of other script jobs.
#Specifically it created a batch file containing jobs that calculate differential expression and splicing for each feature type
#This file must then be executed or farmed out to a cluster for execution

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $comparison_file = '';
my $ensembl_version = '';
my $script_dir = '';
my $analysis_dir = '';
my $project_name = '';
my $outfile = '';
my $comparison_type = '';
my $silent = '';

GetOptions ('comparison_file=s'=>\$comparison_file, 'ensembl_version=i'=>\$ensembl_version, 'script_dir=s'=>\$script_dir, 'analysis_dir=s'=>\$analysis_dir, 'project_name=s'=>\$project_name, 'outfile=s'=>\$outfile, 'comparison_type=s'=>\$comparison_type, 'silent=i'=>\$silent);

unless($silent){
  print GREEN, "\n\nSpecify a file describing the desired comparisons using: --comparison_file", RESET;
  print GREEN, "\n\tThis file should have the following format:", RESET;
  print GREEN, "\n\tLibraryA\tLibraryA_Name\tLibraryB\tLibraryB_Name\tComparison_Description\tComparison_name", RESET;
  print GREEN, "\nSpecify the ensembl version used for this analysis using: --ensembl_version", RESET;
  print GREEN, "\nSpecify the root script dir using:  --script_dir", RESET;
  print GREEN, "\nSpecify the root analysis dir using:  --analysis_dir", RESET;
  print GREEN, "\nSpecify the project name using: --project_name", RESET;
  print GREEN, "\nSpecify the name of the output job file using: --outfile", RESET;
  print GREEN, "\nSpecify the comparison type (DE or SI) using: --comparison_type", RESET;
  print GREEN, "\nUse --silent=1 to prevent routine messages from being printed", RESET;
  print GREEN, "\n\ne.g. createComparisonJobs.pl  --comparison_file=Lymphoma_Lib_Comparisons.txt  --ensembl_version=53  --script_dir=/home/malachig/svn/solexa_analysis/  --analysis_dir=/projects/malachig/solexa/  --project_name=Lymphoma  --outfile=calculateDifferentialExpression.sh  --comparison_type=DE\n\n", RESET;
}

unless ($comparison_file && $ensembl_version && $script_dir && $analysis_dir && $project_name && $outfile && $comparison_type){
  print RED, "\nRequired input parameter(s) missing or incorrect format!\n\n", RESET;
  exit();
}
unless ($comparison_type =~ /de|si/i){
  print RED, "\nComparison type must be DE or SI\n\n", RESET;
  exit();
}

#check dirs
$script_dir = &checkDir('-dir'=>$script_dir, '-clear'=>"no");
$analysis_dir = &checkDir('-dir'=>$analysis_dir, '-clear'=>"no");

my %types;
$types{1}{type} = "Gene";
$types{1}{group} = "ENST";
$types{1}{si} = 0;
$types{2}{type} = "Transcript";
$types{2}{group} = "Transcripts";
$types{2}{si} = 1;
$types{3}{type} = "ExonRegion";
$types{3}{group} = "ENST";
$types{3}{si} = 1;
$types{4}{type} = "Junction";
$types{4}{group} = "Junctions";
$types{4}{si} = 1;
$types{5}{type} = "KnownJunction";
$types{5}{group} = "Junctions";
$types{5}{si} = 1;
$types{6}{type} = "NovelJunction";
$types{6}{group} = "Junctions";
$types{6}{si} = 1;
$types{7}{type} = "Boundary";
$types{7}{group} = "Boundaries";
$types{7}{si} = 1;
$types{8}{type} = "KnownBoundary";
$types{8}{group} = "Boundaries";
$types{8}{si} = 1;
$types{9}{type} = "NovelBoundary";
$types{9}{group} = "Boundaries";
$types{9}{si} = 1;
$types{10}{type} = "Intron";
$types{10}{group} = "Introns";
$types{10}{si} = 1;
$types{11}{type} = "ActiveIntronRegion";
$types{11}{group} = "Introns";
$types{11}{si} = 1;
$types{12}{type} = "SilentIntronRegion";
$types{12}{group} = "Introns";
$types{12}{si} = 1;
$types{13}{type} = "Intergenic";
$types{13}{group} = "Intergenics";
$types{13}{si} = 0;
$types{14}{type} = "ActiveIntergenicRegion";
$types{14}{group} = "Intergenics";
$types{14}{si} = 0;
$types{15}{type} = "SilentIntergenicRegion";
$types{15}{group} = "Intergenics";
$types{15}{si} = 0;

#Import comparisons info from user specified file
my %comps;
open (IN, "$comparison_file") || die "\n\nCould not open input file: $comparison_file\n\n";
my $c = 0;
while(<IN>){
  $c++;
  chomp($_);
  my @line = split("\t", $_);
  $comps{$c}{libA} = $line[0];
  $comps{$c}{libA_name} = $line[1];
  $comps{$c}{libB} = $line[2];
  $comps{$c}{libB_name} = $line[3];
  $comps{$c}{description} = $line[4];
  $comps{$c}{name} = $line[5];
}
close (IN);

open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
foreach my $comp (sort {$a <=> $b} keys %comps){
  unless($silent){
    print MAGENTA, "\n\nProcessing comparison: $comp ($comps{$comp}{name})", RESET;
  }

  foreach my $type (sort {$a <=> $b} keys %types){
  unless($silent){
    print MAGENTA, "\n\tProcessing type: $type ($types{$type}{type})", RESET;
  }

  #Example command for differential expression calculations
  #/home/malachig/svn/solexa_analysis/calculateDifferentialExpression.pl  --dataname='Pre_vs_Post_Gene'  --fileA=/projects/malachig/solexa/read_records/HS0499/ENST_v53/Summary/HS0499_GeneExpression_v53.txt  --fileB=/projects/malachig/solexa/read_records/HS0502/ENST_v53/Summary/HS0502_GeneExpression_v53.txt  --columnA='Average_Coverage_RAW'   --columnB='Average_Coverage_RAW' --results_dir=/projects/malachig/solexa/figures_and_stats/DE/Neuroblastoma/ENST_v53/  --short_name='Gene'  --image_type='svg'  --libraryA_name='Post'  --libraryB_name='Pre'

  #Example command for differential splicing calculations
  #/home/malachig/svn/solexa_analysis/calculateDifferentialSplicing.pl  --dataname='Pre_vs_Post_Transcript'  --gene_fileA=/projects/malachig/solexa/read_records/HS0499/ENST_v53/Summary/HS0499_GeneExpression_v53.txt  --gene_fileB=/projects/malachig/solexa/read_records/HS0502/ENST_v53/Summary/HS0502_GeneExpression_v53.txt  --seq_fileA=/projects/malachig/solexa/read_records/HS0499/Transcripts_v53/Summary/HS0499_TranscriptExpression_v53.txt   --seq_fileB=/projects/malachig/solexa/read_records/HS0502/Transcripts_v53/Summary/HS0502_TranscriptExpression_v53.txt  --data_column='Average_Coverage_RAW'   --results_dir=/projects/malachig/solexa/figures_and_stats/SI/Neuroblastoma/Transcripts_v53/  --short_name='Transcript'  --image_type='svg'

    my $comp_name = $comps{$comp}{name};
    my $type_name = $types{$type}{type};
    my $group_name = $types{$type}{group};
    my $dataname = "$comp_name"."_"."$type_name";
    my $libA = $comps{$comp}{libA};
    my $libA_name = $comps{$comp}{libA_name};
    my $libB = $comps{$comp}{libB};
    my $libB_name = $comps{$comp}{libB_name};
    my $comp_id = "$libA"."_"."$libB";
    my $group_name_v = "$group_name"."_v"."$ensembl_version";
    my $gene_name_v = "ENST"."_v"."$ensembl_version";
    my $fileA_name = "$libA"."_"."$type_name"."Expression_v"."$ensembl_version".".txt";
    my $fileA_path = "$analysis_dir"."read_records/$libA/$group_name_v"."/Summary/$fileA_name";
    my $geneFileA_name = "$libA"."_"."Gene"."Expression_v"."$ensembl_version".".txt";
    my $geneFileA_path = "$analysis_dir"."read_records/$libA/$gene_name_v"."/Summary/$geneFileA_name";
    my $geneFileB_name = "$libB"."_"."Gene"."Expression_v"."$ensembl_version".".txt";
    my $geneFileB_path = "$analysis_dir"."read_records/$libB/$gene_name_v"."/Summary/$geneFileB_name";
    my $exonFileA_name = "$libA"."_"."ExonRegion"."Expression_v"."$ensembl_version".".txt";
    my $exonFileA_path = "$analysis_dir"."read_records/$libA/$gene_name_v"."/Summary/$exonFileA_name";
    my $exonFileB_name = "$libB"."_"."ExonRegion"."Expression_v"."$ensembl_version".".txt";
    my $exonFileB_path = "$analysis_dir"."read_records/$libB/$gene_name_v"."/Summary/$exonFileB_name";
    my $fileB_name = "$libB"."_"."$type_name"."Expression_v"."$ensembl_version".".txt";
    my $fileB_path = "$analysis_dir"."read_records/$libB/$group_name_v"."/Summary/$fileB_name";
    my $results_dir = "$analysis_dir"."figures_and_stats/DE/$project_name/$group_name_v";
    my $results_dir_si = "$analysis_dir"."figures_and_stats/SI/$project_name/$group_name_v";
    my $si_applicable = $types{$type}{si};
    $results_dir = &checkDir('-dir'=>$results_dir, '-clear'=>"no");

    if ($comparison_type =~ /de/i){
      #Create DE command
      print OUT "$script_dir"."calculateDifferentialExpression.pl  --dataname=\"$dataname\"  --fileA=$fileA_path  --fileB=$fileB_path  --columnA=Average_Coverage_RAW   --columnB=Average_Coverage_RAW  --results_dir=$results_dir  --short_name=$type_name  --image_type=svg  --libraryA_name=$libA_name  --libraryB_name=$libB_name  --comp_id=$comp_id\n";
    }elsif(($comparison_type =~ /si/i) && ($si_applicable)){
      #Create SI command
      print OUT "$script_dir"."calculateDifferentialSplicing.pl  --dataname=\"$dataname\"  --gene_fileA=$geneFileA_path  --gene_fileB=$geneFileB_path  --exon_fileA=$exonFileA_path  --exon_fileB=$exonFileB_path  --seq_fileA=$fileA_path   --seq_fileB=$fileB_path  --data_column=Average_Coverage_RAW  --results_dir=$results_dir_si  --short_name=$type_name  --image_type=svg  --comp_id=$comp_id\n"

    }
  }
}
close(OUT);

exit();
