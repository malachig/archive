#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use website::WEB qw(:all);
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $annotation_dir = '';
my $junction_seq_size = '';
my $type = '';
my $id_list = '';
my $outfile = '';
my $min_skipped_bases = '';
my $max_skipped_bases = '';

GetOptions('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'annotation_dir=s'=>\$annotation_dir, 'junction_seq_size=i'=>\$junction_seq_size, 
           'type=s'=>\$type, 'id_list=s'=>\$id_list,
           'outfile=s'=>\$outfile, 'min_skipped_bases=i'=>\$min_skipped_bases, 'max_skipped_bases=i'=>\$max_skipped_bases);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script generates HTML pages to summarize expression, DE and SI results", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the base path to annotation files using:  --annotation_dir", RESET;
print GREEN, "\n\tSpecify the junction database sequence length using:  --junction_seq_size", RESET;
print GREEN, "\n\tSpecify the sequence type to be verified using:  --type", RESET;
print GREEN, "\n\tSpecify a file containing a list of ids of features to be verified using:  --id_list (first column should contain IDs)", RESET;
print GREEN, "\n\tSpecify an output file for design sequences (will be written in fasta format) using: --outfile", RESET;
print GREEN, "\n\tIf desired, specify a min and max allowable skipped bases using:  --min_skipped_bases (e.g. 100) and --max_skipped_bases (e.g. 500)", RESET;
print GREEN, "\n\nUsage:  getSequenceForPcrDesign.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --junction_seq_size=62  --annotation_dir=/projects/malachig/sequence_databases/hs_53_36o/  --type=junction  --id_list=ids.txt  --outfile=design_sequences.txt\n\n", RESET;

unless ($database && $server && $user && $password && $annotation_dir && $junction_seq_size && $type && $id_list && $outfile){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}
$annotation_dir = &checkDir('-dir'=>$annotation_dir, '-clear'=>"no");

#Get exon content for all genes
print BLUE, "\n\nGet exon content", RESET;
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};
my $g_storable_name = "$database"."_AllGenes_GeneInfo_WithSeq.storable";
my $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"yes", '-storable'=>$g_storable_name, '-silent'=>"yes");
my $ec_storable = "$database"."_AllGenes_ExonContent.storable";
my $gene_exon_content_ref = &getExonContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-storable'=>$ec_storable, '-silent'=>"yes");
$alexa_dbh->disconnect();

#Get the coordinates for this ID from the annotation file
my $annot_file;
if ($type =~ /junction/i){
  $annot_file = "$annotation_dir"."exonJunctions/exonJunctions_"."$junction_seq_size"."mers_annotated.txt.gz";
}

#Get IDs from list
open (IDS, "$id_list") || die "\n\nCould not open id list: $id_list\n\n";
my $header = 1;
my %ids;
my $order = 0;
while(<IDS>){
  chomp($_);
  my @line = split("\t", $_);
  if ($header == 1){
    $header = 0;
    next();
  }
  $order++;
  $ids{$line[0]}{order}=$order;
}
close(IDS);

open (IN, "zcat $annot_file |") || die "\nCould not open input file\n\n";
$header = 1;
my %columns = ();
while(<IN>){
  chomp($_);
  my @line = split ("\t", $_);
  if($header == 1){
    my $column_count = 0;
    foreach my $column (@line){
      $columns{$column}{column_pos} = $column_count;
      $column_count++;
    }
    $header = 0;
    next();
  }

  my $id = $line[0];

  if ($ids{$id}){
    #Create an exon content sequence for this gene with primer3 bracket flanking the target exon skip event
    #Also get the base count for the skipped portion to help determine product sizes later
    #print out in fasta format
    my $gene_id = $line[$columns{'Gene_ID'}{column_pos}];
    my $exons_skipped = $line[$columns{'Exons_Skipped'}{column_pos}];
    my $u1_start =  $line[$columns{'Unit1_start'}{column_pos}];
    my $u1_end =  $line[$columns{'Unit1_end'}{column_pos}];
    my $u2_start =  $line[$columns{'Unit2_start'}{column_pos}];
    my $u2_end =  $line[$columns{'Unit2_end'}{column_pos}];
    print BLUE, "\nFound feature: $id\tgene: $gene_id\tu1_start: $u1_start\tu2_end: $u2_end\n", RESET;

    my $gene_seq = $genes_ref->{$gene_id}->{sequence};

    #Get the exon content for this gene
    my $ec_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};

    my $skipped_bases = 0;
    my $design_seq = '';

    #Get all exon content blocks that have overlap with the coordinates identified
    foreach my $ec (sort {$ec_ref->{$a}->{start} <=> $ec_ref->{$b}->{start}} keys %{$ec_ref}){
      my $ec_start = $ec_ref->{$ec}->{start};
      my $ec_end = $ec_ref->{$ec}->{end};
      my $size = ($ec_end-$ec_start)+1;

      my $exon_seq = substr ($gene_seq, $ec_start, $size);

      #Does the current feature overlap the current exon content block - find the left and right exon
      if ($u1_start >= $ec_start && $u1_start <= $ec_end){
        #Found left exon
        $design_seq .= "$exon_seq\[";
        #print YELLOW, "$exon_seq\[", RESET;
      }elsif($u2_end >= $ec_start && $u2_end <= $ec_end){
        #Found right exon
        $design_seq .= "\]$exon_seq";
        #print YELLOW, "\]$exon_seq", RESET;
      }else{
        $design_seq .= "$exon_seq";
        #print CYAN, "$exon_seq", RESET;
      }

      #If this region is completely within the inner coordinates of the junction, count the size towards the skipped bases
      if ($ec_start > $u1_end && $ec_end < $u2_start){
        $skipped_bases+=$size;
      }
    }

    #If min and max skipped bases were specified, skip junctions that fail the criteria
    if ($min_skipped_bases && $max_skipped_bases){
      if ($skipped_bases < $min_skipped_bases || $skipped_bases > $max_skipped_bases){
        $ids{$id}{status} = 0;
        print YELLOW, "\n\tNumber of skipped bases ($skipped_bases) did not meet criteria ($min_skipped_bases - $max_skipped_bases)", RESET;
        next();
      }
      $ids{$id}{status} = 1;
      $ids{$id}{record} = ">EJ$id $genes_ref->{$gene_id}->{gene_name} S$exons_skipped SB$skipped_bases\n$design_seq\n";
    }else{
      $ids{$id}{record} = ">EJ$id $genes_ref->{$gene_id}->{gene_name} S$exons_skipped SB$skipped_bases\n$design_seq\n";
      $ids{$id}{status} = 1;
    }

    #Header line J$id Gene_Name Exons_skipped Bases_Skipped
    print ">EJ$id $genes_ref->{$gene_id}->{gene_name} S$exons_skipped SB$skipped_bases\n$design_seq\n";
  }
}
close(IN);

#Now printout the records to the output file in the same order as the input file
open (OUT, ">$outfile") || die "\nCould not open output file for writing\n\n";
foreach my $id (sort {$ids{$a}->{order} <=> $ids{$b}->{order}} keys %ids){
  if ($ids{$id}{status}){
    print OUT "$ids{$id}{record}";
  }
}
close(OUT);


exit();








