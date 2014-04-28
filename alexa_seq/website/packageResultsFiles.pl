#!/usr/bin/perl -w
#Written by Malachi Griffith and Obi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to package up results files for a project to be made available on the website
#Files will be made available in a single large tarball file as well as individual (linked) files on the individual website summary pages.
#Target files: expression values for all feature types, differential expression values for all feature types, splicing index values for all feature types


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

my $project_name = '';
my $analysis_dir = '';
my $ensembl_version = '';
my $web_dir = '';

GetOptions ('project_name=s'=>\$project_name, 'analysis_dir=s'=>\$analysis_dir, 'ensembl_version=i'=>\$ensembl_version, 'web_dir=s'=>\$web_dir);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the main analysis directory using: --analysis_dir", RESET;
print GREEN, "\n\tSpecify the project name for this project using:  --project_name", RESET;
print GREEN, "\n\tSpecify the ensembl version used for this analysis using: --ensembl_version", RESET;
print GREEN, "\n\nExample:  packageResultsFiles.pl  --project_name=FL_Trans  --analysis_dir=/projects/malachig/alexa_seq/analysis  --ensembl_version=53  --web_dir=/projects/malachig/alexa_seq/www/htdocs/alexa-seq/FL_Trans/data\n\n", RESET;

unless ($project_name && $analysis_dir && $ensembl_version && $web_dir){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}
$analysis_dir = &checkDir('-dir'=>$analysis_dir, '-clear'=>"no");
$web_dir = &checkDir('-dir'=>$web_dir, '-clear'=>"no");

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
#Make target dir for package that will be deleted once tar created
my $target_dir = "$web_dir"."$project_name"."_ResultFiles/"; $dirs_created += &make_dir('-dir'=>$target_dir);
$target_dir = &checkDir('-dir'=>$target_dir, '-clear'=>"yes", '-recursive'=>"yes");

my $exp_dir = "$target_dir"."EXP/"; $dirs_created += &make_dir('-dir'=>$exp_dir);
my $de_dir = "$target_dir"."DE/"; $dirs_created += &make_dir('-dir'=>$de_dir);
my $ae_dir = "$target_dir"."AE/"; $dirs_created += &make_dir('-dir'=>$ae_dir);
my $matrix_dir = "$target_dir"."EXP/Matrix/"; $dirs_created += &make_dir('-dir'=>$matrix_dir);

my @lib_list;
foreach my $lib (sort{$a cmp $b} keys %lib_names){
  #Make directories for this library
  my $lib_dir = "$exp_dir"."$lib/"; $dirs_created += &make_dir('-dir'=>$lib_dir);

  #Copy expression files to expression dir for this library
  my $cp_cmd1 = "cp $analysis_dir"."read_records/$lib/*/Summary/*Expression*.txt $lib_dir";
  my $ls_cmd1 = "ls $analysis_dir"."read_records/$lib/*/Summary/*Expression*.txt";
  print BLUE, "\n\n$cp_cmd1", RESET;
  system($cp_cmd1);
}

#Copy differential expression files to DE dir
my $cp_cmd2 = "cp $analysis_dir"."figures_and_stats/DE"."/$project_name/*/*DE_Values_Sorted.txt $de_dir";
print BLUE, "\n\n$cp_cmd2", RESET;
system($cp_cmd2);

#Copy alternative expression files to AE dir
my $cp_cmd3 = "cp $analysis_dir"."figures_and_stats/SI/$project_name/*/*SI_Values_Sorted.txt $ae_dir";
print BLUE, "\n\n$cp_cmd3", RESET;
system($cp_cmd3);

#Copy matrix files to matrix dir
my $cp_cmd4 = "cp $web_dir"."matrix/* $matrix_dir";
print BLUE, "\n\n$cp_cmd4", RESET;
system($cp_cmd4);


#Make a tarball out the whole thing
chdir($web_dir);
my $tar_cmd = "tar -cf $project_name"."_ResultFiles.tar $project_name"."_ResultFiles/";
print BLUE, "\n\n$tar_cmd", RESET;
system($tar_cmd);

my $gz_cmd = "gzip -f --best $project_name"."_ResultFiles.tar";
system($gz_cmd);

my $rm_cmd = "rm -rf $project_name"."_ResultFiles/";
system($rm_cmd);

exit();


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

