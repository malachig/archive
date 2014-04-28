package utilities::ucsc;
require Exporter;

#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

@ISA = qw( Exporter );
@EXPORT = qw();

@EXPORT_OK = qw(&importUcscAlignments);

%EXPORT_TAGS = (
     all => [qw(&importUcscAlignments)]
);

use strict;
use Data::Dumper;
use Term::ANSIColor qw(:constants);


###############################################################################################################################
#Load in alignments from an a standard UCSC align file (ESTs or mRNA - either target species or all other species)            #
###############################################################################################################################
sub importUcscAlignments{
  my %args = @_;
  my $file = $args{'-align_file'};
  my $type = $args{'-seq_type'};
  my $filter = $args{'-filter'};
  my $chr_filter = $args{'-chr_filter'};
  my $start_filter = $args{'-start_filter'};
  my $end_filter = $args{'-end_filter'};
  my $gb_org_map_ref = $args{'-gb_org_map'};
  my $import_type = $args{'-import_type'};

  unless($import_type eq "exon" || $import_type eq "junction"){
    print RED, "\nMust specify import type for importUcscAlignments()\n\n", RESET;
    exit();
  }

  #If the user specifies to store coords '-store_coords'=>1, then store the strand, @starts and @exon_sizes

  my %seq_coords;
  my $seq_import_count = 0;
  my $nr_exon_count = 0;

  $| = 1; print BLUE, "\n\nImporting UCSC alignments of the type: $type (from $file)", RESET; $| = 0;

  unless(-e $file){
    print YELLOW, "\n\tFile not found - returning empty hash\n\n", RESET;
    return(\%seq_coords);
  }

  open (ALIGN, "zcat $file |") || die "\nCould not open UCSC alignment file: $file\n\n";

  while (<ALIGN>){
    chomp($_);
    my @line = split("\t", $_);
    my $matches = $line[1];
    my $mismatches = $line[2];
    my $n_count = $line[4];
    my $query_insertions = $line[5];

    my $strand = $line[9];
    my $acc_id = $line[10];
    my $chr = $line[14];
    my $target_size = $line[15];
    my @exon_sizes = split(",", $line[19]);
    my @starts = split(",", $line[21]);

    my $species_id = $gb_org_map_ref->{$acc_id};
    unless($species_id){
      next();
    }

    #Note that for xenoMRNA and xenoESTs, the strand field has two strand values (this is because translated BLAT was used to create the alignment)
    #The first strand corresponds to the strand of the seq translated.  The second strand refers to the alignment to the genome
    #Since ESTs are often not cloned directionally, you need to try the translation of both strands
    #mRNAs are usually directional so they always appear as ++ or +-
    #In any case, from an alignment perspective, it is the second strand record that we should be trying to match to our known gene models
    #Note that for xmrna and xest tables, if the strand is '--' or '+-' then the coordinates are relative to the end of the chromosome instead of the beginning and need to be converted
    my $convert = 0;
    if ($strand eq "++" || $strand eq "-+"){
      $strand = "+";
    }
    if ($strand eq "--" || $strand eq "+-"){
      $strand = "-";
      $convert = 1;
    }

    my @adjusted_starts;
    my @ends;
    my @tmp_sizes = @exon_sizes;
    foreach my $start (@starts){
      my $exon_size = shift (@tmp_sizes);
      my $end = $start + $exon_size;

      #If neccessary, convert these coordinates
      if ($convert == 1){
        $end = $target_size-$start;
        $start = $end-$exon_size;
      }
      push (@ends, $end);
      push (@adjusted_starts, $start+1);
    }
    my $exon_count = scalar(@adjusted_starts);

    #Skip single exon alignments
    if ($exon_count == 1){next();}

    #Make sure the alignment has some overlap with the current target region
    my @tmp = (@adjusted_starts, @ends);
    my @tmp_sort = sort {$a <=> $b} @tmp;
    my $lower = $tmp_sort[0];
    my $upper = $tmp_sort[scalar(@tmp_sort)-1];
    my $chr_tmp = "chr"."$chr_filter";
    unless (($chr eq $chr_tmp) && ($lower >= $start_filter && $lower <= $end_filter) && ($upper >= $start_filter && $upper <= $end_filter)){
      next();
    }

    #Quality filter
    #Use no quality filter for mRNAs
    #Use stringent filter for ESTs
    #Use relaxed filter for xmRNA and xESTs (which will have more mismatches, gaps, etc. due to divergence rather than sequence errors)
    if ($filter == 0){
      #Do nothing 
    }elsif($filter == 1){
      unless($matches > 300 && $mismatches < 2 && $n_count < 2 && $query_insertions < 2){
        next();
      }
    }elsif($filter == 2){
      unless($matches > 200 && $n_count < 10){
        next();
      }
    }else{
      print RED, "\nFilter type provided to importUcscAlignments() not recognized!\n\n", RESET;
      exit();
    }

    my @coord_list;
    if ($import_type eq "exon"){
      #GENERATE EXONS
      for (my $i = 0; $i <= $exon_count-1; $i++){
        my $start = $adjusted_starts[$i]; 
        my $end = $ends[$i];
        $seq_import_count++;
        my $seq_exon = "$start"."_"."$end";
        push(@coord_list, $seq_exon);
      }
    }elsif($import_type eq "junction"){
      #GENERATE JUNCTIONS
      for (my $i = 0; $i < $exon_count-1; $i++){
        $seq_import_count++;
        my $seq_junction = "$ends[$i]"."_"."$adjusted_starts[$i+1]";
        push(@coord_list, $seq_junction);
      }
    }else{
      print RED, "\nImport type not recognized in importUcscAlignments() - must be 'exon' or 'junction'\n\n", RESET;
      exit();
    }

    #Store this junction in the appropriate part of the hash
    my $chr_strand = "$chr"."_"."$strand";

    if ($seq_coords{$chr_strand}){
      #If the current chromosome_strand combo has been observed before get the appropriate hash reference

      my $h_ref = $seq_coords{$chr_strand};
      foreach my $exon (@coord_list){
        if ($h_ref->{$exon}){
	  $h_ref->{$exon}->{c}++;
          my $s_ref = $h_ref->{$exon}->{s};
          $s_ref->{$species_id} = '';
        }else{
          $h_ref->{$exon}->{c} = 1;
          my %tmp;
          $tmp{$species_id} = '';
          $h_ref->{$exon}->{s} = \%tmp;
          $nr_exon_count++;
        }
      }
    }else{
      #If the current chromosome_strand combo has NOT been observed before, create the new hash and store as a berkleyDB file
      my %h;

      foreach my $exon (@coord_list){
        if ($h{$exon}){
	  $h{$exon}{c}++;
          my $s_ref = $h{$exon}{s};
          $s_ref->{$species_id} = '';
        }else{
          $h{$exon}{c} = 1;
          my %tmp;
          $tmp{$species_id} = '';
          $h{$exon}{s} = \%tmp;
          $nr_exon_count++;
        }
      }
      $seq_coords{$chr_strand} = \%h;
    }
  }

  close(ALIGN);
  print BLUE, "\n\tImported $seq_import_count $type $import_type coords ($nr_exon_count non-redundant) from the UCSC table file\n", RESET;

  return(\%seq_coords);
}


1;

