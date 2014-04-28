#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#This script takes an input file containing exons predicted by cufflinks or scripture (and then processed with gtf2exon.pl or bed2exon.pl)
#1.) Import a list of reference exons to search for (exon regions from ALEXA-seq) - store these by chromosome and strand
#ER100256        RTN2    19      -1      50684444        50684651

#2.) Go through a list of exons identified by Cufflinks or Scripture
#Scripture
#chr1:4220-4696  .       se      1(chr1:4219-4696)       0.0
#Cufflinks
#chrX:149405927-149406376        +       se      CUFF.3075476.1  0.1923128897

#Look for exons that overlap a target exon on the same strand
#Store the 'expression' value for the exon (cumulative score and number of matching exons)

#Overlap types
#1 = perfect match only
#2 = perfect match or contained within or flanking target exons
#3 = an overlap at all

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

my $target_regions_file = '';
my $comparison_regions_file = '';
my $overlap_type = '';
my $outfile = '';
my $min_percent_overlap = '';

GetOptions ('target_regions_file=s'=>\$target_regions_file, 'comparison_regions_file=s'=>\$comparison_regions_file, 'overlap_type=i'=>\$overlap_type, 'outfile=s'=>\$outfile, 'min_percent_overlap=f'=>\$min_percent_overlap);

unless ($target_regions_file && $comparison_regions_file && $overlap_type && $outfile){
  print GREEN, "\n\nExample: mergeCufflinksScriptureExons.pl  --target_regions_file=/projects/alexa2/alexa_seq/figures_and_stats/PlatformComparisons/qPCR/ValidationExonRegions.txt  --comparison_regions_file=/projects/alexa2/alexa_seq/figures_and_stats/PlatformComparisons/lane_by_lane_tophat/HS04391/cufflink_exons.csv  --overlap_type=2  --outfile=test.txt  [--min_percent_overlap=50]\n\n", RESET;
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Import the target regions
my %targets;
open (TR, $target_regions_file) || die "\n\nCould not open input file: $target_regions_file\n\n";
while(<TR>){
  chomp($_);
  my @line = split("\t", $_);
  my $id = $line[0];
  my $name = $line[1];
  my $chr = $line[2];
  my $strand = $line[3];
  my $start = $line[4];
  my $end = $line[5];
  $chr = "chr$chr";
  my $coord_string = "$chr:$start-$end";

  if ($strand eq "1"){
    $strand = "+";
  }
  if ($strand eq "-1"){
    $strand = "-";
  }
  my $chr_strand = "$chr"."$strand";
  if ($targets{$chr_strand}){
    my $ref = $targets{$chr_strand};
    $ref->{$coord_string}->{chromosome} = $chr;
    $ref->{$coord_string}->{start} = $start;
    $ref->{$coord_string}->{end} = $end;
    $ref->{$coord_string}->{strand} = $strand;
    $ref->{$coord_string}->{id} = $id;
    $ref->{$coord_string}->{name} = $name;
    $ref->{$coord_string}->{cum_exp} = 0;
    $ref->{$coord_string}->{matches} = 0;
  }else{
    my %tmp;
    $tmp{$coord_string}{chromosome} = $chr;
    $tmp{$coord_string}{start} = $start;
    $tmp{$coord_string}{end} = $end;
    $tmp{$coord_string}{strand} = $strand;
    $tmp{$coord_string}{id} = $id;
    $tmp{$coord_string}{name} = $name;
    $tmp{$coord_string}{cum_exp} = 0;
    $tmp{$coord_string}{matches} = 0;
    $targets{$chr_strand} = \%tmp;
  }
}
close(TR);

#Now go through each line in the Cufflinks/Scripture files and look for matching exons
open (CR, $comparison_regions_file) || die "\n\nCould not open infile: $comparison_regions_file\n\n";
while(<CR>){
  chomp($_);
  my @line = split("\t", $_);
  my $coord_string = $line[0];
  my $chr = '';
  my $start = '';
  my $end = '';
  if ($coord_string =~ /(.*)\:(\d+)\-(\d+)/){
    $chr = $1;
    $start = $2;
    $end = $3;
  }else{
    print RED, "\n\nExpecting coord string of the format chr1:1-100\n\n", RESET;
    exit();
  }

  my $strand = $line[1];
  my $type = $line[2];
  my $name = $line[3];
  my $exp = $line[4];
  
  my $chr_strand = "$chr"."$strand";

  #If there was nothing stored for this chr_strand, no point in proceeding
  unless($targets{$chr_strand}){
    next();
  }

  my $ref = $targets{$chr_strand};
  my $match_coord_string;

  my $overlap = 0;
  my $otype = "na";
  my $percent_overlap_f = 0;
  my $gene_name = '';
  if ($overlap_type == 1){
    #Exact match
    if ($ref->{$coord_string}){
      $overlap = 1;
      $match_coord_string = $coord_string;
      $otype = "Perfect";
      $percent_overlap_f = 100.00;
      $gene_name = $ref->{$match_coord_string}->{name};
    }
  }elsif($overlap_type == 2){
    #Exact match
    if ($ref->{$coord_string}){
      $overlap = 1;
      $match_coord_string = $coord_string;
      $otype = "Perfect";
      $percent_overlap_f = 100.00;
      $gene_name = $ref->{$match_coord_string}->{name};
    }else{
      #Or completely contained within, or completely flanking
      foreach my $cs (keys %{$ref}){
        my $ref_start = $ref->{$cs}{start};
        my $ref_end = $ref->{$cs}{end};
        if ($start >= $ref_start && $start <= $ref_end && $end >= $ref_start && $end <= $ref_end){
          #Contained within ref exon region
          $overlap = 1;
          $match_coord_string = "$chr:$ref_start-$ref_end";
          my $percent_overlap = ((($end - $start)+1)/(($ref_end - $ref_start)+1))*100;
          $percent_overlap_f = sprintf("%.2f", $percent_overlap);
          $gene_name = $ref->{$match_coord_string}->{name};
          $otype = "Contained";
        }elsif ($start <= $ref_start && $end >= $ref_end){
          #Flanking ref exon region
          $overlap = 1;
          $match_coord_string = "$chr:$ref_start-$ref_end";
          my $percent_overlap = ((($ref_end - $ref_start)+1)/(($end - $start)+1))*100;
          $percent_overlap_f = sprintf("%.2f", $percent_overlap);
          $gene_name = $ref->{$match_coord_string}->{name};
          $otype="Flanking";
        }
      }
    }
  }elsif($overlap_type == 3){
    #Exact match
    if ($ref->{$coord_string}){
      $overlap = 1;
      $match_coord_string = $coord_string;
      $otype = "Perfect";
      $percent_overlap_f = 100.00;
      $gene_name = $ref->{$match_coord_string}->{name};
    }else{
      #Or completely contained within, or completely flanking
      foreach my $cs (keys %{$ref}){
        my $ref_start = $ref->{$cs}{start};
        my $ref_end = $ref->{$cs}{end};
        if ($start >= $ref_start && $start <= $ref_end && $end >= $ref_start && $end <= $ref_end){
          #Contained within ref exon region
          $overlap = 1;
          $match_coord_string = "$chr:$ref_start-$ref_end";
          my $percent_overlap = ((($end - $start)+1)/(($ref_end - $ref_start)+1))*100;
          $percent_overlap_f = sprintf("%.2f", $percent_overlap);
          $gene_name = $ref->{$match_coord_string}->{name};
          $otype = "Contained";
        }elsif ($start <= $ref_start && $end >= $ref_end){
          #Flanking ref exon region
          $overlap = 1;
          $match_coord_string = "$chr:$ref_start-$ref_end";
          my $percent_overlap = ((($ref_end - $ref_start)+1)/(($end - $start)+1))*100;
          $percent_overlap_f = sprintf("%.2f", $percent_overlap);
          $gene_name = $ref->{$match_coord_string}->{name};
          $otype="Flanking";
        }elsif($start >= $ref_start && $start <= $ref_end){
          #Left overlap
          $overlap = 1;
          $match_coord_string = "$chr:$ref_start-$ref_end";          
          my $percent_overlap = ((($ref_end - $start)+1)/(($ref_end - $ref_start)+1))*100;
          $percent_overlap_f = sprintf("%.2f", $percent_overlap);
          $gene_name = $ref->{$match_coord_string}->{name};
          $otype="Left_Overlap";
        }elsif($end >= $ref_start && $end <= $ref_end){
          #right overlap
          $overlap = 1;
          $match_coord_string = "$chr:$ref_start-$ref_end";          
          my $percent_overlap = ((($end - $ref_start)+1)/(($ref_end - $ref_start)+1))*100;
          $percent_overlap_f = sprintf("%.2f", $percent_overlap);
          $gene_name = $ref->{$match_coord_string}->{name};
          $otype="Right_Overlap";  
        }
      }
    }
  }else{
    print RED, "\n\nOverlap type: $overlap_type not understood\n\n", RESET;
    exit();
  }

  if ($overlap){
    if ($min_percent_overlap){
      if ($percent_overlap_f > $min_percent_overlap){
        print "$gene_name\t$_\t$otype\t$percent_overlap_f\n";
        #Store the exp value...
        $ref->{$match_coord_string}->{cum_exp}+=$exp;
        $ref->{$match_coord_string}->{matches}++;
      }
    }else{
      print "$gene_name\t$_\t$otype\t$percent_overlap_f\n";
      $ref->{$match_coord_string}->{cum_exp}+=$exp;
      $ref->{$match_coord_string}->{matches}++;
    }
  }
}
close(CR);


#Now go through all the target regions identified and print out the expression values found (and how many exons they are derived from)
open (OUT, ">$outfile") || die "\n\nCould not open outfile: $outfile\n\n";
foreach my $chr_strand (sort keys %targets){  
  my $ref = $targets{$chr_strand};
  foreach my $cs (sort keys %{$ref}){
    my $chr = $ref->{$cs}->{chromosome};
    my $start = $ref->{$cs}->{start};
    my $end = $ref->{$cs}->{end};
    my $strand = $ref->{$cs}->{strand};
    my $id = $ref->{$cs}->{id};
    my $name = $ref->{$cs}->{name};
    my $cum_exp = $ref->{$cs}->{cum_exp};
    my $matches = $ref->{$cs}->{matches};

    print OUT "$cs\t$name\t$chr\t$strand\t$start\t$end\t$matches\t$cum_exp\n";
  }
}
close(OUT);

exit();



















