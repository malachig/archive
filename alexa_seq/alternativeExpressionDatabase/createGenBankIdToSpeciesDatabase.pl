#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Create a berkeley DB to store mappings of GenBank IDs to organism IDs

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;
use BerkeleyDB;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $ucsc_align_dir = '';
my $out_dir = '';

GetOptions ('ucsc_align_dir=s'=>\$ucsc_align_dir,'out_dir=s'=>\$out_dir);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a directory containing UCSC gbCdnaInfo.txt.gz files using: --ucsc_align_dir", RESET;
print GREEN, "\n\tSpecify an output directory for partitioned files using:  --out_dir", RESET;
print GREEN, "\n\nExample: createGenBankIdToSpeciesDatabase.pl  --ucsc_align_dir=/projects/malachig/sequence_databases/hs_53_36o/mrna_est/   --out_dir=/projects/malachig/sequence_databases/hs_53_36o/mrna_est/\n\n", RESET;

unless ($ucsc_align_dir && $out_dir){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Check the temp dir before proceeding.  Make sure it is empty
$ucsc_align_dir = &checkDir('-dir'=>$ucsc_align_dir, '-clear'=>"no");
$out_dir = &checkDir('-dir'=>$out_dir, '-clear'=>"no");

my $infile = "$ucsc_align_dir"."gbCdnaInfo.txt.gz";
my $db_file = "$out_dir"."GenBankToOrganism.btree";

print YELLOW, "\n\tCreating binary tree file: $db_file", RESET;
my %h;
tie(%h, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $db_file , -Flags => DB_CREATE) or die "can't open file $db_file: $! $BerkeleyDB::Error\n";

print BLUE, "\n\nCreate Berkeley DB : $db_file\n", RESET;
open (MAP, "zcat $infile |") || die "\nCould not open UCSC genbank-to-organism file: $db_file\n\n";
my $c = 0;
while(<MAP>){
  chomp($_);
  $c++;
  if ($c == 10000){
    $| = 1; print "."; $| = 0;
    $c = 0;
  }
  my @line = split("\t", $_);
  my $acc = $line[1];
  my $organism_id = $line[7]; 
  $h{$acc} = $organism_id;
  #print YELLOW, "\n$acc\t$organism_id", RESET;
}
close (MAP);
untie(%h);
print "\n\n";
exit();






