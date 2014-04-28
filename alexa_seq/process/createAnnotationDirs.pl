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

my $annotation_dir = '';
my $species_build = '';

GetOptions ('annotation_dir=s'=>\$annotation_dir, 'species_build=s'=>\$species_build);

if ($annotation_dir && $species_build){
  #Parameters ok
}else{
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  print GREEN, "\n\nExample: createAnnotationDirs.pl  --annotation_dir=/projects/malachig/sequence_databases/  --species_build=hs_53_36o\n\n", RESET;
  exit();
}

$annotation_dir = &checkDir('-dir'=>$annotation_dir, '-clear'=>"no");


#Create dirs
my $dir = '';
my $dirs_created = 0;

#Create dirs
$dir = "$annotation_dir"."$species_build/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/alexa_db/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/genes/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/genes/temp/"; $dirs_created += &make_dir('-dir'=>$dir);

$dir = "$annotation_dir"."$species_build/exonBoundaries/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/exonBoundaries/temp/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/exonBoundaries/blastdb/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/exonBoundaries/bwadb/"; $dirs_created += &make_dir('-dir'=>$dir);

$dir = "$annotation_dir"."$species_build/exonJunctions/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/exonJunctions/temp/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/exonJunctions/blastdb/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/exonJunctions/bwadb/"; $dirs_created += &make_dir('-dir'=>$dir);
  
$dir = "$annotation_dir"."$species_build/exonRegions/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/exonRegions/temp/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/exonRegions/blastdb/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/exonRegions/blastdb/temp/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/exonRegions/bwadb/"; $dirs_created += &make_dir('-dir'=>$dir);

$dir = "$annotation_dir"."$species_build/intergenics/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/intergenics/temp/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/intergenics/blastdb/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/intergenics/bwadb/"; $dirs_created += &make_dir('-dir'=>$dir);
  
$dir = "$annotation_dir"."$species_build/introns/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/introns/temp/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/introns/blastdb/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/introns/bwadb/"; $dirs_created += &make_dir('-dir'=>$dir);

$dir = "$annotation_dir"."$species_build/transcripts/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/transcripts/temp/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/transcripts/blastdb/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/transcripts/bwadb/"; $dirs_created += &make_dir('-dir'=>$dir);

$dir = "$annotation_dir"."$species_build/logs/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/logs/createExonRegionDatabase/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/logs/createExonJunctionDatabase/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/logs/createExonBoundaryDatabase/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/logs/annotateExonJunctions/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/logs/annotateExonBoundaries/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/logs/createIntronDatabase/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/logs/createIntergenicDatabase/"; $dirs_created += &make_dir('-dir'=>$dir);

$dir = "$annotation_dir"."$species_build/mrna_est/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/mrna_est/partitions/"; $dirs_created += &make_dir('-dir'=>$dir);

$dir = "$annotation_dir"."$species_build/jobs/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/jobs/temp/"; $dirs_created += &make_dir('-dir'=>$dir);

$dir = "$annotation_dir"."$species_build/repeats/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/repeats/blastdb/"; $dirs_created += &make_dir('-dir'=>$dir);
$dir = "$annotation_dir"."$species_build/repeats/bwadb/"; $dirs_created += &make_dir('-dir'=>$dir);

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

