=head1 NAME

WEB.pm - Basic methods for generating HTML to summarize Illumina expression data

=head1 SYNOPSIS

use WEB qw(:all);

=head2 NOTE

Currently located in '~/svn/solexa_analysis/website/'

=head2 RECENT CHANGES

Various modifications.  Last modified 27 May 2009

=head1 DESCRIPTION

Generic utility for generating HTML to summarize Illumina expression data files.

=head1 EXAMPLES

use lib '~/svn/solexa_analysis';

use website::WEB qw(:all);

=head1 SEE ALSO

None

=head1 BUGS

Contact author via email

=head1 AUTHOR

Written by Malachi Griffith (malachig@bcgsc.ca) 

=head1 ACKNOWLEDGEMENTS

University of British Columbia Graduate Studies

Michael Smith Foundation for Health Research

Natural Sciences and Engineering Research Council

Genome British Columbia

Terry Fox Foundation

National Cancer Institute of Canada

BC Cancer Foundation

=head1 AFFLIATIONS

Malachi Griffith is supervised by Marco A. Marra

Genome Sciences Centre, BC Cancer Research Centre, BC Cancer Agency, UBC Faculty of Medicine - Medical Genetics

=head1 SUBROUTINES

=cut

package website_medium::WEB;
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

@EXPORT_OK = qw(&getAnnotations &annotateDataTypes &getDataPaths &getPartitions &getExpressionData &getDifferentialExpressionData &getDifferentialSplicingData &rankData &assignPartitions &generateExpressionContent &generateDifferentialExpressionContent &generateDifferentialSplicingContent &writePage);

%EXPORT_TAGS = (
     all => [qw(&getAnnotations &annotateDataTypes &getDataPaths &getPartitions &getExpressionData &getDifferentialExpressionData &getDifferentialSplicingData &rankData &assignPartitions &generateExpressionContent &generateDifferentialExpressionContent &generateDifferentialSplicingContent &writePage)]
);

use strict;
use Data::Dumper;
use Term::ANSIColor qw(:constants);

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::ALEXA_DB qw(:all);


###########################################################################################################################################
#Build an object to contain all annotations                                                                                               #
###########################################################################################################################################
sub getAnnotations{
  my %args = @_;
  my $database = $args{'-database'};
  my $dbh = $args{'-dbh'};
  my $annotation_dir = $args{'-annotation_dir'};
  my $junction_seq_size = $args{'-junction_seq_size'};
  my $boundary_seq_size = $args{'-boundary_seq_size'};
  my $types_list_ref = $args{'-types_list'};
  my $gene_id_list_ref = $args{'-gene_id_list'}; #Only specify if you wish to limit analysis to a specific list of genes

  my $gene_filter = 0;
  if ($gene_id_list_ref){
    $gene_filter = 1;
    my $gene_count = scalar(@{$gene_id_list_ref});
    print BLUE, "\n\nImporting annotations for all possible sequence types - for $gene_count genes\n", RESET;
  }else{
    print BLUE, "\n\nImporting annotations for all possible sequence types - for all genes\n", RESET;
  }

  my $db_name;
  if ($database =~ /ALEXA\_(\w+)/){
    $db_name = $1;
  }

  #Create a hash for all annotation types - each links to a hash of annotation info keyed on unique id for each seq type
  #Annotation types needed:
  #Gene Transcript ExonRegion Junction KnownJunction NovelJunction Boundary KnownBoundary NovelBoundary Intron ActiveIntronRegion SilentIntronRegion Intergenic ActiveIntergenicRegion SilentIntergenicRegion

  #Values that must be defined for all sequence feature types:
  #id: unique id for the feature
  #seq_name: name of sequence (gene name, transcript name, exon name, junction name, etc.)
  #chromosome
  #chr_start: start of sequence on chromosome (lower number always)
  #chr_end: end of sequence on chromosome (higher number always)
  #known: 1 if corresponds to an EnsEMBL sequence, 0 for hypothetical junctions, boundaries, etc.
  #supporting_seqs: number of supporting mRNA or EST sequences
  #conserved_species_count:  number of species with mRNA/EST seqs that match the event...

  #Global vars
  my %a;
  my $file = '';
  my $header = '';
  my %columns = ();


  #######
  #'Gene' - For genes only data will be retrieved from both the ALEXA database as well as an annotation file

  my @gene_ids = @{&getAllGenes ('-dbh'=>$dbh, '-gene_type'=>'All')};
  my $gene_storable = "$database"."_AllGenes_GeneInfo_NoSeq.storable";
  my $genes_ref = &getGeneInfo ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no", '-storable'=>$gene_storable, '-silent'=>"yes");

  foreach my $gene_id (keys %{$genes_ref}){
    $genes_ref->{$gene_id}->{gene_id} = $gene_id;
    $genes_ref->{$gene_id}->{seq_name} = $genes_ref->{$gene_id}->{gene_name};
    $genes_ref->{$gene_id}->{known} = 1;
    $genes_ref->{$gene_id}->{supporting_seqs} = "N/A";
    $genes_ref->{$gene_id}->{conserved_species_count} = "N/A";
    my $tmp = $genes_ref->{$gene_id}->{chr_start};
    if ($tmp > $genes_ref->{$gene_id}->{chr_end}){
      $genes_ref->{$gene_id}->{chr_start} = $genes_ref->{$gene_id}->{chr_end};
      $genes_ref->{$gene_id}->{chr_end} = $tmp;
    }
    my $chromosome = $genes_ref->{$gene_id}->{chromosome};
    if ($chromosome eq "M"){$chromosome = "MT";}
  }

  #Additional data from the gene annotation file
  $file = "$annotation_dir"."genes/genes_annotated.txt.gz";
  open(IN, "zcat $file |") || die "\nCould not open file: $file\n\n";
  $header = 1;
  %columns = ();
  while(<IN>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header == 1){
      my $column_count = 0;
      foreach my $column (@line){
        $columns{$column}{column_pos} = $column_count;
        $column_count++;
      }
      $header = 0;
      next();
    }
    my $id = $line[0];
    my $gene_id = $line[$columns{'Gene_ID'}{column_pos}];
    $genes_ref->{$gene_id}->{base_count} = $line[$columns{'Base_Count'}{column_pos}];
    $genes_ref->{$gene_id}->{unmasked_base_count} = $line[$columns{'UnMasked_Base_Count'}{column_pos}];
    $genes_ref->{$gene_id}->{coding_base_count} = $line[$columns{'Coding_Base_Count'}{column_pos}];
    $genes_ref->{$gene_id}->{fid} = $line[$columns{'FID'}{column_pos}];
  }
  close(IN);

  $a{'Gene'} = $genes_ref;
  print BLUE, ".", RESET;


  #############
  #'Transcript'
  if ($types_list_ref->{'Transcript'}){
    my %transcripts;
    $file = "$annotation_dir"."transcripts/transcripts_annotated.txt.gz";
    open(IN, "zcat $file |") || die "\nCould not open file: $file\n\n";
    $header = 1;
    %columns = ();
    while(<IN>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        my $column_count = 0;
        foreach my $column (@line){
          $columns{$column}{column_pos} = $column_count;
          $column_count++;
        }
        $header = 0;
        next();
      }
      my $id = $line[0];
      my $gene_id = $line[$columns{'Gene_ID'}{column_pos}];
      unless (($gene_filter) && (&skipGene('-gene_list'=>$gene_id_list_ref, '-gene_id'=>$gene_id))){
        $transcripts{$id}{gene_id} = $gene_id;
        $transcripts{$id}{base_count} = $line[$columns{'Base_Count'}{column_pos}];
        $transcripts{$id}{unmasked_base_count} = $line[$columns{'UnMasked_Base_Count'}{column_pos}];
        $transcripts{$id}{coding_base_count} = $line[$columns{'Coding_Base_Count'}{column_pos}];
        $transcripts{$id}{seq_name} = $line[$columns{'Seq_Name'}{column_pos}];
        $transcripts{$id}{known} = 1;
        $transcripts{$id}{supporting_seqs} = "N/A";
        $transcripts{$id}{conserved_species_count} = "N/A";
        $transcripts{$id}{transcript_size} = $line[$columns{'Transcript_Size'}{column_pos}];
        $transcripts{$id}{feature_list} = $line[$columns{'Feature_List'}{column_pos}];
        my $chromosome = $line[$columns{'Chromosome'}{column_pos}];
        if ($chromosome eq "M"){$chromosome = "MT";}
        $transcripts{$id}{chromosome} = $chromosome;
        $transcripts{$id}{chr_start} = $line[$columns{'Unit1_start_chr'}{column_pos}];
        $transcripts{$id}{chr_end} = $line[$columns{'Unit1_end_chr'}{column_pos}];
        $transcripts{$id}->{fid} = $line[$columns{'FID'}{column_pos}];
      }
    }
    close(IN);
    $a{'Transcript'} = \%transcripts;
    print BLUE, ".", RESET;
  }


  #############
  #'ExonRegion'
  if ($types_list_ref->{'ExonRegion'}){
    my %exon_regions;
    $file = "$annotation_dir"."exonRegions/exonRegions_annotated.txt.gz";
    open(IN, "zcat $file |") || die "\nCould not open file: $file\n\n";
    $header = 1;
    %columns = ();
    while(<IN>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        my $column_count = 0;
        foreach my $column (@line){
          $columns{$column}{column_pos} = $column_count;
          $column_count++;
        }
        $header = 0;
        next();
      }
      my $id = $line[0];
      my $gene_id = $line[$columns{'Gene_ID'}{column_pos}];
      unless (($gene_filter) && (&skipGene('-gene_list'=>$gene_id_list_ref, '-gene_id'=>$gene_id))){
        $exon_regions{$id}{gene_id} = $gene_id;
        $exon_regions{$id}{base_count} = $line[$columns{'Base_Count'}{column_pos}];
        $exon_regions{$id}{unmasked_base_count} = $line[$columns{'UnMasked_Base_Count'}{column_pos}];
        $exon_regions{$id}{coding_base_count} = $line[$columns{'Coding_Base_Count'}{column_pos}];
        $exon_regions{$id}{seq_name} = $line[$columns{'Seq_Name'}{column_pos}];
        $exon_regions{$id}{known} = 1;
        $exon_regions{$id}{supporting_seqs} = 0;
        $exon_regions{$id}{conserved_species_count} = $line[$columns{'Conserved_Species_Count'}{column_pos}];
        $exon_regions{$id}{fid} = $line[$columns{'FID'}{column_pos}];
        my $chromosome = $line[$columns{'Chromosome'}{column_pos}];
        if ($chromosome eq "M"){$chromosome = "MT";}
        $exon_regions{$id}{chromosome} = $chromosome;
        $exon_regions{$id}{chr_start} = $line[$columns{'Unit1_start_chr'}{column_pos}];
        $exon_regions{$id}{chr_end} = $line[$columns{'Unit1_end_chr'}{column_pos}];
        $exon_regions{$id}{model_start} = $line[$columns{'Unit1_start'}{column_pos}];
        $exon_regions{$id}{model_end} = $line[$columns{'Unit1_end'}{column_pos}];
        my $mrna_count = $line[$columns{'Supporting_mRNA_Count'}{column_pos}];
        my $est_count = $line[$columns{'Supporting_EST_Count'}{column_pos}];
        if ($mrna_count eq "na"){$mrna_count = 0;}
        if ($est_count eq "na"){$est_count = 0;}
        $exon_regions{$id}{supporting_seqs} += $mrna_count;
        $exon_regions{$id}{supporting_seqs} += $est_count;
      }
    }
    close(IN);
    $a{'ExonRegion'} = \%exon_regions;
    print BLUE, ".", RESET;
  }

  ####################################################################################################################
  #'Junction' 'KnownJunction' 'NovelJunction'  #Create a single object, but link to all three names for convenience...
  if ($types_list_ref->{'Junction'}){
    my %junctions;

    $file = "$annotation_dir"."exonJunctions/exonJunctions_"."$junction_seq_size"."mers_annotated.txt.gz";
    open(IN, "zcat $file |") || die "\nCould not open file: $file\n\n";
    $header = 1;
    %columns = ();
    while(<IN>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        my $column_count = 0;
        foreach my $column (@line){
          $columns{$column}{column_pos} = $column_count;
          $column_count++;
        }
        $header = 0;
        next();
      }
      my $id = $line[0];
      my $gene_id = $line[$columns{'Gene_ID'}{column_pos}];
      unless (($gene_filter) && (&skipGene('-gene_list'=>$gene_id_list_ref, '-gene_id'=>$gene_id))){
        $junctions{$id}{gene_id} = $gene_id;
        $junctions{$id}{base_count} = $line[$columns{'Base_Count'}{column_pos}];
        $junctions{$id}{unmasked_base_count} = $line[$columns{'UnMasked_Base_Count'}{column_pos}];
        $junctions{$id}{coding_base_count} = $line[$columns{'Coding_Base_Count'}{column_pos}];
        $junctions{$id}{seq_name} = $line[$columns{'Seq_Name'}{column_pos}];
        $junctions{$id}{known} = 0;
        $junctions{$id}{supporting_seqs} = 0;
        $junctions{$id}{conserved_species_count} = $line[$columns{'Conserved_Species_Count'}{column_pos}];
        $junctions{$id}{fid} = $line[$columns{'FID'}{column_pos}];
        $junctions{$id}{exons_skipped} = $line[$columns{'Exons_Skipped'}{column_pos}];
        my $chromosome = $line[$columns{'Chromosome'}{column_pos}];
        if ($chromosome eq "M"){$chromosome = "MT";}
        $junctions{$id}{chromosome} = $chromosome;
        $junctions{$id}{chr_start} = $line[$columns{'Unit1_start_chr'}{column_pos}];
        $junctions{$id}{chr_end} = $line[$columns{'Unit2_end_chr'}{column_pos}];
        $junctions{$id}{model_start} = $line[$columns{'Unit1_end'}{column_pos}];
        $junctions{$id}{model_end} = $line[$columns{'Unit2_start'}{column_pos}];
        if ($line[$columns{'Supporting_EnsEMBL_Count'}{column_pos}] > 0){
          $junctions{$id}{known} = 1;
        }
        my $mrna_count = $line[$columns{'Supporting_mRNA_Count'}{column_pos}];
        my $est_count = $line[$columns{'Supporting_EST_Count'}{column_pos}];
        if ($mrna_count eq "na"){$mrna_count = 0;}
        if ($est_count eq "na"){$est_count = 0;}
        $junctions{$id}{supporting_seqs} += $mrna_count;
        $junctions{$id}{supporting_seqs} += $est_count;
      }
    }
    close(IN);
    $a{'Junction'} = \%junctions;
    if ($types_list_ref->{'KnownJunction'}){
      $a{'KnownJunction'} = \%junctions;
    }
    if ($types_list_ref->{'NovelJunction'}){
      $a{'NovelJunction'} = \%junctions;
    }
    print BLUE, ".", RESET;
  }

  ####################################################################################################################
  #'Boundary' 'KnownBoundary' 'NovelBoundary'  #Create a single object, but link to all three names for convenience...
  if ($types_list_ref->{'Boundary'}){
    my %boundaries;

    $file = "$annotation_dir"."exonBoundaries/exonBoundaries_"."$boundary_seq_size"."mers_annotated.txt.gz";
    open(IN, "zcat $file |") || die "\nCould not open file: $file\n\n";
    $header = 1;
    %columns = ();
    while(<IN>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        my $column_count = 0;
        foreach my $column (@line){
          $columns{$column}{column_pos} = $column_count;
          $column_count++;
        }
        $header = 0;
        next();
      }
      my $id = $line[0];
      my $gene_id = $line[$columns{'Gene_ID'}{column_pos}];
      unless (($gene_filter) && (&skipGene('-gene_list'=>$gene_id_list_ref, '-gene_id'=>$gene_id))){
        $boundaries{$id}{gene_id} = $gene_id;
        $boundaries{$id}{base_count} = $line[$columns{'Base_Count'}{column_pos}];
        $boundaries{$id}{unmasked_base_count} = $line[$columns{'UnMasked_Base_Count'}{column_pos}];
        $boundaries{$id}{coding_base_count} = $line[$columns{'Coding_Base_Count'}{column_pos}];
        $boundaries{$id}{seq_name} = $line[$columns{'Seq_Name'}{column_pos}];
        $boundaries{$id}{known} = 0;
        $boundaries{$id}{supporting_seqs} = 0;
        $boundaries{$id}{conserved_species_count} = $line[$columns{'Conserved_Species_Count'}{column_pos}];
        $boundaries{$id}{fid} = $line[$columns{'FID'}{column_pos}];
        my $chromosome =  $line[$columns{'Chromosome'}{column_pos}];
        if ($chromosome eq "M"){$chromosome = "MT";}
        $boundaries{$id}{chromosome} = $chromosome;
        $boundaries{$id}{chr_start} = $line[$columns{'Unit1_start_chr'}{column_pos}];
        $boundaries{$id}{chr_end} = $line[$columns{'Unit1_end_chr'}{column_pos}];
        $boundaries{$id}{model_start} = $line[$columns{'Unit1_start'}{column_pos}];
        $boundaries{$id}{model_end} = $line[$columns{'Unit1_end'}{column_pos}];
        if ($line[$columns{'Supporting_EnsEMBL_Count'}{column_pos}] > 0){
          $boundaries{$id}{known} = 1;
        }
        my $mrna_count = $line[$columns{'Supporting_mRNA_Count'}{column_pos}];
        my $est_count = $line[$columns{'Supporting_EST_Count'}{column_pos}];
        if ($mrna_count eq "na"){$mrna_count = 0;}
        if ($est_count eq "na"){$est_count = 0;}
        $boundaries{$id}{supporting_seqs} += $mrna_count;
        $boundaries{$id}{supporting_seqs} += $est_count;
      }
    }
    close(IN);
    $a{'Boundary'} = \%boundaries;
    if ($types_list_ref->{'KnownBoundary'}){
      $a{'KnownBoundary'} = \%boundaries;
    }
    if ($types_list_ref->{'NovelBoundary'}){
      $a{'NovelBoundary'} = \%boundaries;
    }
    print BLUE, ".", RESET;
  }

  #########
  #'Intron'
  if ($types_list_ref->{'Intron'}){
    my %introns;
    $file = "$annotation_dir"."introns/introns_annotated.txt.gz";
    open(IN, "zcat $file |") || die "\nCould not open file: $file\n\n";
    $header = 1;
    %columns = ();
    while(<IN>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        my $column_count = 0;
        foreach my $column (@line){
          $columns{$column}{column_pos} = $column_count;
          $column_count++;
        }
        $header = 0;
        next();
      }
      my $id = $line[0];
      my $gene_id = $line[$columns{'Gene_ID_List'}{column_pos}];
      unless (($gene_filter) && (&skipGene('-gene_list'=>$gene_id_list_ref, '-gene_id'=>$gene_id))){
        $introns{$id}{gene_id} = $gene_id;
        $introns{$id}{base_count} = $line[$columns{'Base_Count'}{column_pos}];
        $introns{$id}{unmasked_base_count} = $line[$columns{'UnMasked_Base_Count'}{column_pos}];
        $introns{$id}{coding_base_count} = $line[$columns{'Coding_Base_Count'}{column_pos}];
        $introns{$id}{seq_name} = $line[$columns{'Seq_Name'}{column_pos}];
        $introns{$id}{known} = 0;
        $introns{$id}{supporting_seqs} = 0;
        $introns{$id}{conserved_species_count} = $line[$columns{'Conserved_Species_Count'}{column_pos}];
        $introns{$id}{fid} = $line[$columns{'FID'}{column_pos}];
        my $chromosome = $line[$columns{'Chromosome'}{column_pos}];
        if ($chromosome eq "M"){$chromosome = "MT";}
        $introns{$id}{chromosome} = $chromosome;
        $introns{$id}{chr_start} = $line[$columns{'Unit1_start_chr'}{column_pos}];
        $introns{$id}{chr_end} = $line[$columns{'Unit1_end_chr'}{column_pos}];
        $introns{$id}{model_start} = $line[$columns{'Unit1_start'}{column_pos}];
        $introns{$id}{model_end} = $line[$columns{'Unit1_end'}{column_pos}];
        my $mrna_count = $line[$columns{'Supporting_mRNA_Count'}{column_pos}];
        my $est_count = $line[$columns{'Supporting_EST_Count'}{column_pos}];
        if ($mrna_count eq "na"){$mrna_count = 0;}
        if ($est_count eq "na"){$est_count = 0;}
        $introns{$id}{supporting_seqs} += $mrna_count;
        $introns{$id}{supporting_seqs} += $est_count;
      }
    }
    close(IN);
    $a{'Intron'} = \%introns;
    print BLUE, ".", RESET;
  }

  #####################
  #'ActiveIntronRegion'
  if ($types_list_ref->{'ActiveIntronRegion'}){
    my %active_intron_regions;
    $file = "$annotation_dir"."introns/activeIntronRegions.txt.gz";
    open(IN, "zcat $file |") || die "\nCould not open file: $file\n\n";
    $header = 1;
    %columns = ();
    while(<IN>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        my $column_count = 0;
        foreach my $column (@line){
          $columns{$column}{column_pos} = $column_count;
          $column_count++;
        }
        $header = 0;
        next();
      }
      my $id = $line[0];
      my $gene_id = $line[$columns{'Gene_ID'}{column_pos}];
      unless (($gene_filter) && (&skipGene('-gene_list'=>$gene_id_list_ref, '-gene_id'=>$gene_id))){
        $active_intron_regions{$id}{gene_id} = $gene_id;
        $active_intron_regions{$id}{base_count} = $line[$columns{'Base_Count'}{column_pos}];
        $active_intron_regions{$id}{unmasked_base_count} = $line[$columns{'UnMasked_Base_Count'}{column_pos}];
        $active_intron_regions{$id}{coding_base_count} = $line[$columns{'Coding_Base_Count'}{column_pos}];
        $active_intron_regions{$id}{seq_name} = $line[$columns{'Seq_Name'}{column_pos}];
        $active_intron_regions{$id}{known} = 0;
        $active_intron_regions{$id}{supporting_seqs} = 1;
        $active_intron_regions{$id}{conserved_species_count} = $line[$columns{'Conserved_Species_Count'}{column_pos}];
        $active_intron_regions{$id}{fid} = $line[$columns{'FID'}{column_pos}];
        my $chromosome = $line[$columns{'Chromosome'}{column_pos}];
        if ($chromosome eq "M"){$chromosome = "MT";}
        $active_intron_regions{$id}{chromosome} = $chromosome;
        $active_intron_regions{$id}{chr_start} = $line[$columns{'Unit1_start_chr'}{column_pos}];
        $active_intron_regions{$id}{chr_end} = $line[$columns{'Unit1_end_chr'}{column_pos}];
        $active_intron_regions{$id}{model_start} = $line[$columns{'Unit1_start'}{column_pos}];
        $active_intron_regions{$id}{model_end} = $line[$columns{'Unit1_end'}{column_pos}];
        my $mrna_count = $line[$columns{'Supporting_mRNA_Count'}{column_pos}];
        my $est_count = $line[$columns{'Supporting_EST_Count'}{column_pos}];
        if ($mrna_count eq "na"){$mrna_count = 0;}
        if ($est_count eq "na"){$est_count = 0;}
        $active_intron_regions{$id}{supporting_seqs} += $mrna_count;
        $active_intron_regions{$id}{supporting_seqs} += $est_count;
      }
    }
    close(IN);
    $a{'ActiveIntronRegion'} = \%active_intron_regions;
    print BLUE, ".", RESET;
  }

  #####################
  #'SilentIntronRegion'
  if ($types_list_ref->{'SilentIntronRegion'}){
    my %silent_intron_regions;
    $file = "$annotation_dir"."introns/silentIntronRegions.txt.gz";
    open(IN, "zcat $file |") || die "\nCould not open file: $file\n\n";
    $header = 1;
    %columns = ();
    while(<IN>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        my $column_count = 0;
        foreach my $column (@line){
          $columns{$column}{column_pos} = $column_count;
          $column_count++;
        }
        $header = 0;
        next();
      }
      my $id = $line[0];
      my $gene_id = $line[$columns{'Gene_ID'}{column_pos}];
      unless (($gene_filter) && (&skipGene('-gene_list'=>$gene_id_list_ref, '-gene_id'=>$gene_id))){
        $silent_intron_regions{$id}{gene_id} = $gene_id;
        $silent_intron_regions{$id}{base_count} = $line[$columns{'Base_Count'}{column_pos}];
        $silent_intron_regions{$id}{unmasked_base_count} = $line[$columns{'UnMasked_Base_Count'}{column_pos}];
        $silent_intron_regions{$id}{coding_base_count} = $line[$columns{'Coding_Base_Count'}{column_pos}];
        $silent_intron_regions{$id}{seq_name} = $line[$columns{'Seq_Name'}{column_pos}];
        $silent_intron_regions{$id}{known} = 0;
        $silent_intron_regions{$id}{supporting_seqs} = 0;
        $silent_intron_regions{$id}{conserved_species_count} = $line[$columns{'Conserved_Species_Count'}{column_pos}];
        $silent_intron_regions{$id}{fid} = $line[$columns{'FID'}{column_pos}];
        my $chromosome = $line[$columns{'Chromosome'}{column_pos}];
        if ($chromosome eq "M"){$chromosome = "MT";}
        $silent_intron_regions{$id}{chromosome} = $chromosome;
        $silent_intron_regions{$id}{chr_start} = $line[$columns{'Unit1_start_chr'}{column_pos}];
        $silent_intron_regions{$id}{chr_end} = $line[$columns{'Unit1_end_chr'}{column_pos}];
        $silent_intron_regions{$id}{model_start} = $line[$columns{'Unit1_start'}{column_pos}];
        $silent_intron_regions{$id}{model_end} = $line[$columns{'Unit1_end'}{column_pos}];
        my $mrna_count = $line[$columns{'Supporting_mRNA_Count'}{column_pos}];
        my $est_count = $line[$columns{'Supporting_EST_Count'}{column_pos}];
        if ($mrna_count eq "na"){$mrna_count = 0;}
        if ($est_count eq "na"){$est_count = 0;}
        $silent_intron_regions{$id}{supporting_seqs} += $mrna_count;
        $silent_intron_regions{$id}{supporting_seqs} += $est_count;
      }
    }
    close(IN);
    $a{'SilentIntronRegion'} = \%silent_intron_regions;
    print BLUE, ".", RESET;
  }

  #############
  #'Intergenic'
  my %intergenics;
  if ($types_list_ref->{'Intergenic'}){
    $file = "$annotation_dir"."intergenics/intergenics_annotated.txt.gz";
    open(IN, "zcat $file |") || die "\nCould not open file: $file\n\n";
    $header = 1;
    %columns = ();
    while(<IN>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        my $column_count = 0;
        foreach my $column (@line){
          $columns{$column}{column_pos} = $column_count;
          $column_count++;
        }
        $header = 0;
        next();
      }
      my $id = $line[0];
      #Get IDs of closest gene to intergenic region
      my @list;
      my $us_gene = $line[$columns{'Upstream_Gene_ID'}{column_pos}];
      my $ds_gene = $line[$columns{'Downstream_Gene_ID'}{column_pos}];
      unless ($us_gene eq "NA"){push(@list, $us_gene);}
      unless ($ds_gene eq "NA"){push(@list, $ds_gene);}
      my $gene_id = "@list";

      unless (($gene_filter) && (&skipGene('-gene_list'=>$gene_id_list_ref, '-gene_id'=>$gene_id))){
        $intergenics{$id}{gene_id} = $gene_id;
        $intergenics{$id}{base_count} = $line[$columns{'Base_Count'}{column_pos}];
        $intergenics{$id}{unmasked_base_count} = $line[$columns{'UnMasked_Base_Count'}{column_pos}];
        $intergenics{$id}{coding_base_count} = $line[$columns{'Coding_Base_Count'}{column_pos}];
        $intergenics{$id}{seq_name} = $line[$columns{'Seq_Name'}{column_pos}];
        $intergenics{$id}{known} = 0;
        $intergenics{$id}{supporting_seqs} = 0;
        $intergenics{$id}{conserved_species_count} = $line[$columns{'Conserved_Species_Count'}{column_pos}];
        $intergenics{$id}{fid} = $line[$columns{'FID'}{column_pos}];
        my $chromosome = $line[$columns{'Chromosome'}{column_pos}];
        if ($chromosome eq "M"){$chromosome = "MT";}
        $intergenics{$id}{chromosome} = $chromosome;
        $intergenics{$id}{chr_start} = $line[$columns{'Unit1_start_chr'}{column_pos}];
        $intergenics{$id}{chr_end} = $line[$columns{'Unit1_end_chr'}{column_pos}];
        my $mrna_count = $line[$columns{'Supporting_mRNA_Count'}{column_pos}];
        my $est_count = $line[$columns{'Supporting_EST_Count'}{column_pos}];
        if ($mrna_count eq "na"){$mrna_count = 0;}
        if ($est_count eq "na"){$est_count = 0;}
        $intergenics{$id}{supporting_seqs} += $mrna_count;
        $intergenics{$id}{supporting_seqs} += $est_count;
      }
    }
    close(IN);
    $a{'Intergenic'} = \%intergenics;
    print BLUE, ".", RESET;
  }

  #########################
  #'ActiveIntergenicRegion'
  if ($types_list_ref->{'ActiveIntergenicRegion'}){
    my %active_intergenic_regions;
    $file = "$annotation_dir"."intergenics/activeIntergenicRegions.txt.gz";
    open(IN, "zcat $file |") || die "\nCould not open file: $file\n\n";
    $header = 1;
    %columns = ();
    while(<IN>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        my $column_count = 0;
        foreach my $column (@line){
          $columns{$column}{column_pos} = $column_count;
          $column_count++;
        }
        $header = 0;
        next();
      }
      my $id = $line[0];
      my $ig_id = $line[$columns{'Intergenic_ID'}{column_pos}];
      if ($intergenics{$ig_id}){
        my $gene_id = $intergenics{$ig_id}{gene_id};
        unless (($gene_filter) && (&skipGene('-gene_list'=>$gene_id_list_ref, '-gene_id'=>$gene_id))){
          $active_intergenic_regions{$id}{gene_id} = $gene_id;
          $active_intergenic_regions{$id}{base_count} = $line[$columns{'Base_Count'}{column_pos}];
          $active_intergenic_regions{$id}{unmasked_base_count} = $line[$columns{'UnMasked_Base_Count'}{column_pos}];
          $active_intergenic_regions{$id}{coding_base_count} = $line[$columns{'Coding_Base_Count'}{column_pos}];
          $active_intergenic_regions{$id}{seq_name} = $line[$columns{'Seq_Name'}{column_pos}];
          $active_intergenic_regions{$id}{known} = 0;
          $active_intergenic_regions{$id}{supporting_seqs} = 1;
          $active_intergenic_regions{$id}{conserved_species_count} = $line[$columns{'Conserved_Species_Count'}{column_pos}];
          $active_intergenic_regions{$id}{fid} = $line[$columns{'FID'}{column_pos}];
          my $chromosome = $line[$columns{'Chromosome'}{column_pos}];
          if ($chromosome eq "M"){$chromosome = "MT";}
          $active_intergenic_regions{$id}{chromosome} = $chromosome;
          $active_intergenic_regions{$id}{chr_start} = $line[$columns{'Unit1_start_chr'}{column_pos}];
          $active_intergenic_regions{$id}{chr_end} = $line[$columns{'Unit1_end_chr'}{column_pos}];
          my $mrna_count = $line[$columns{'Supporting_mRNA_Count'}{column_pos}];
          my $est_count = $line[$columns{'Supporting_EST_Count'}{column_pos}];
          if ($mrna_count eq "na"){$mrna_count = 0;}
          if ($est_count eq "na"){$est_count = 0;}
          $active_intergenic_regions{$id}{supporting_seqs} += $mrna_count;
          $active_intergenic_regions{$id}{supporting_seqs} += $est_count;
        }
      }
    }
    close(IN);
    $a{'ActiveIntergenicRegion'} = \%active_intergenic_regions;
    print BLUE, ".", RESET;
  }

  #########################
  #'SilentIntergenicRegion'
  if ($types_list_ref->{'SilentIntergenicRegion'}){
    my %silent_intergenic_regions;
    $file = "$annotation_dir"."intergenics/silentIntergenicRegions.txt.gz";
    open(IN, "zcat $file |") || die "\nCould not open file: $file\n\n";
    $header = 1;
    %columns = ();
    while(<IN>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        my $column_count = 0;
        foreach my $column (@line){
          $columns{$column}{column_pos} = $column_count;
          $column_count++;
        }
        $header = 0;
        next();
      }
      my $id = $line[0];
      my $ig_id = $line[$columns{'Intergenic_ID'}{column_pos}];
      if ($intergenics{$ig_id}){
        my $gene_id = $intergenics{$ig_id}{gene_id};
        unless (($gene_filter) && (&skipGene('-gene_list'=>$gene_id_list_ref, '-gene_id'=>$gene_id))){
          $silent_intergenic_regions{$id}{gene_id} = $gene_id;
          $silent_intergenic_regions{$id}{base_count} = $line[$columns{'Base_Count'}{column_pos}];
          $silent_intergenic_regions{$id}{unmasked_base_count} = $line[$columns{'UnMasked_Base_Count'}{column_pos}];
          $silent_intergenic_regions{$id}{coding_base_count} = $line[$columns{'Coding_Base_Count'}{column_pos}];
          $silent_intergenic_regions{$id}{seq_name} = $line[$columns{'Seq_Name'}{column_pos}];
          $silent_intergenic_regions{$id}{known} = 0;
          $silent_intergenic_regions{$id}{supporting_seqs} = 0;
          $silent_intergenic_regions{$id}{conserved_species_count} = $line[$columns{'Conserved_Species_Count'}{column_pos}];
          $silent_intergenic_regions{$id}{fid} = $line[$columns{'FID'}{column_pos}];
          my $chromosome = $line[$columns{'Chromosome'}{column_pos}];
          if ($chromosome eq "M"){$chromosome = "MT";}
          $silent_intergenic_regions{$id}{chromosome} = $chromosome;
          $silent_intergenic_regions{$id}{chr_start} = $line[$columns{'Unit1_start_chr'}{column_pos}];
          $silent_intergenic_regions{$id}{chr_end} = $line[$columns{'Unit1_end_chr'}{column_pos}];
          my $mrna_count = $line[$columns{'Supporting_mRNA_Count'}{column_pos}];
          my $est_count = $line[$columns{'Supporting_EST_Count'}{column_pos}];
          if ($mrna_count eq "na"){$mrna_count = 0;}
          if ($est_count eq "na"){$est_count = 0;}
          $silent_intergenic_regions{$id}{supporting_seqs} += $mrna_count;
          $silent_intergenic_regions{$id}{supporting_seqs} += $est_count;
        }
      }
    }
    close(IN);
    $a{'SilentIntergenicRegion'} = \%silent_intergenic_regions;
    print BLUE, ".", RESET;
  }

  #Summarize the annotations found...
  foreach my $type (sort keys %a){
    my $a_ref = $a{$type};
    my $c = keys %{$a_ref};
    print BLUE, "\n\tImported $c annotation records of the type: $type", RESET;
  }
  
  return(\%a);
}


###########################################################################################################################################
#If a gene is not in the allowed list and should be skipped, return 1                                                                     #
###########################################################################################################################################
sub skipGene{
  my %args = @_;
  my @gene_list = @{$args{'-gene_list'}};
  my $test_gene_string = $args{'-gene_id'};
  my $skip = 1;

  unless ($test_gene_string){
    print RED, "\nInvalid gene string in skipGene()!", RESET;
    exit();
  }

  #Test gene string may be a list of space seperate gene IDs (for introns and intergenic regions) or a single gene id
  my @test_genes = split(" ", $test_gene_string);
  my %genes_a;
  my %genes_b;

  foreach my $gene_id (@test_genes){
    $genes_a{$gene_id}=1;
  }
  foreach my $gene_id (@gene_list){
    $genes_b{$gene_id}=1;
  }
  
  foreach my $gene_id (keys %genes_a){
    if ($genes_b{$gene_id}){
      $skip = 0;
    }
  }

  return($skip);
}

###########################################################################################################################################
#Hardcode display order, etc. for data types
###########################################################################################################################################
sub annotateDataTypes{
  my %args = @_;
  my $expression_data_ref = $args{'-data_ref'};

  print BLUE, "\n\nAnnotating data paths ...", RESET;

  #Initialize values in case a data type is not hard coded here
  foreach my $data_type (keys %{$expression_data_ref}){
    $expression_data_ref->{$data_type}->{order} = 1;
    $expression_data_ref->{$data_type}->{level} = 1;
  }

  if ($expression_data_ref->{'Gene'}){
    $expression_data_ref->{'Gene'}->{order} = 1;
    $expression_data_ref->{'Gene'}->{level} = 1; 
  }
  if ($expression_data_ref->{'Transcript'}){
    $expression_data_ref->{'Transcript'}->{order} = 2;
    $expression_data_ref->{'Transcript'}->{level} = 1; 
  }
  if ($expression_data_ref->{'ExonRegion'}){
    $expression_data_ref->{'ExonRegion'}->{order} = 3;
    $expression_data_ref->{'ExonRegion'}->{level} = 1;
  }
  if ($expression_data_ref->{'Junction'}){
    $expression_data_ref->{'Junction'}->{order} = 4;
    $expression_data_ref->{'Junction'}->{level} = 1; 
  }
  if ($expression_data_ref->{'KnownJunction'}){
    $expression_data_ref->{'KnownJunction'}->{order} = 5;
    $expression_data_ref->{'KnownJunction'}->{level} = 2; 
  }
  if ($expression_data_ref->{'NovelJunction'}){
    $expression_data_ref->{'NovelJunction'}->{order} = 6;
    $expression_data_ref->{'NovelJunction'}->{level} = 2; 
  }
  if ($expression_data_ref->{'Boundary'}){
    $expression_data_ref->{'Boundary'}->{order} = 7;
    $expression_data_ref->{'Boundary'}->{level} = 1; 
  }
  if ($expression_data_ref->{'KnownBoundary'}){
    $expression_data_ref->{'KnownBoundary'}->{order} = 8;
    $expression_data_ref->{'KnownBoundary'}->{level} = 2; 
  }
  if ($expression_data_ref->{'NovelBoundary'}){
    $expression_data_ref->{'NovelBoundary'}->{order} = 9;
    $expression_data_ref->{'NovelBoundary'}->{level} = 2; 
  }
  if ($expression_data_ref->{'Intron'}){
    $expression_data_ref->{'Intron'}->{order} = 10;
    $expression_data_ref->{'Intron'}->{level} = 1; 
  }
  if ($expression_data_ref->{'ActiveIntronRegion'}){
    $expression_data_ref->{'ActiveIntronRegion'}->{order} = 11;
    $expression_data_ref->{'ActiveIntronRegion'}->{level} = 2; 
  }
  if ($expression_data_ref->{'SilentIntronRegion'}){
    $expression_data_ref->{'SilentIntronRegion'}->{order} = 12;
    $expression_data_ref->{'SilentIntronRegion'}->{level} = 2; 
  }
  if ($expression_data_ref->{'Intergenic'}){
    $expression_data_ref->{'Intergenic'}->{order} = 13;
    $expression_data_ref->{'Intergenic'}->{level} = 1; 
  }
  if ($expression_data_ref->{'ActiveIntergenicRegion'}){
    $expression_data_ref->{'ActiveIntergenicRegion'}->{order} = 14;
    $expression_data_ref->{'ActiveIntergenicRegion'}->{level} = 2; 
  }
  if ($expression_data_ref->{'SilentIntergenicRegion'}){
    $expression_data_ref->{'SilentIntergenicRegion'}->{order} = 15;
    $expression_data_ref->{'SilentIntergenicRegion'}->{level} = 2; 
  }
  return();
}


#############################################################################################################################################
#Get data paths from an input file                                                                                                          #
#############################################################################################################################################
sub getDataPaths{
  my %args = @_;
  my $file = $args{'-file'};

  print BLUE, "\n\nRetrieving data paths ...", RESET;

  unless (-e $file){
    print RED, "\n\nPaths file ($file) specified by user does not exist!  (&getDataPaths)\n\n", RESET;
    exit();
  }

  my %paths;

  open(PATHS, "$file") || die "\nCould not open file: $file\n\n";

  my $header = 1;
  my $line_order = 0;
  while(<PATHS>){
    chomp($_);
    my @line = split (" ", $_);
    #Skip header line
    if ($header == 1){
      $header = 0;
      next();
    }
    #Skip empty lines
    unless ($_ =~ /\w+|\d+/){
      next();
    }
    #Skip comment lines
    if ($_ =~ /^\#/){
      next();
    }
    $line_order++;
    my $library = $line[0];
    my $sample_name = $line[1];
    my $data_type = $line[2];
    my $data_path = $line[3];
    my $stats_path = $line[4];

    unless ($data_path =~ /\/$/){
      $data_path = "$data_path"."/";
    }
    unless ($stats_path =~ /\/$/){
      $stats_path = "$stats_path"."/";
    }
    unless (-e $data_path && -d $data_path){
      print RED, "\nDirectory: $data_path does not appear to be valid! - check paths file\n\n", RESET;
      exit();
    }
    unless (-e $stats_path && -d $stats_path){
      print RED, "\nDirectory: $stats_path does not appear to be valid! - check paths file\n\n", RESET;
      exit();
    }

    if ($paths{$data_type}){
      my $paths_ref = $paths{$data_type}{paths};
      $paths_ref->{$library}->{data_path} = $data_path;
      $paths_ref->{$library}->{stats_path} = $stats_path;
      $paths_ref->{$library}->{name} = $sample_name;
      $paths_ref->{$library}->{line_order} = $line_order;
    }else{
      my %tmp;
      $tmp{$library}{data_path} = $data_path;
      $tmp{$library}{stats_path} = $stats_path;
      $tmp{$library}{name} = $sample_name;
      $tmp{$library}{line_order} = $line_order;
      $paths{$data_type}{paths} = \%tmp;
    }
  }

  close(PATHS);

  return(\%paths);
}


###########################################################################################################################################
#Import partition values                                                                                                                  #
###########################################################################################################################################
sub getPartitions{
  my %args = @_;
  my $file = $args{'-file'};

  print BLUE, "\n\nImporting genome partition info ...", RESET;

  unless (-e $file){
    print RED, "\nPartitions file does not exist\n\n", RESET;
    exit();
  }

  my %partitions;

  open (PART, "$file") || die "\nCould not open partitions file $file\n\n";
  my $header = 1;
  my $c = 0;
  while(<PART>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header == 1){
      $header = 0;
      next();
    }
    $c++;
    my $chromosome = $line[0];
    
    if ($partitions{$chromosome}){
      my $part_ref = $partitions{$chromosome}{partitions};
      $part_ref->{$c}->{chromosome} = $line[0];
      $part_ref->{$c}->{region} = $line[1];
      $part_ref->{$c}->{start} = $line[2];
      $part_ref->{$c}->{end} = $line[3];
      $part_ref->{$c}->{size} = $line[4];
      $part_ref->{$c}->{gene_count} = $line[5];
    }else{
      my %tmp;
      $tmp{$c}{chromosome} = $line[0];
      $tmp{$c}{region} = $line[1];
      $tmp{$c}{start} = $line[2];
      $tmp{$c}{end} = $line[3];
      $tmp{$c}{size} = $line[4];
      $tmp{$c}{gene_count} = $line[5];
      $partitions{$chromosome}{partitions} = \%tmp;
    }
  }
  close(PART);
  
  my $count = keys %partitions;
  print BLUE, "\n\nImported genome partitions for $count chromosomes from $file\n\n", RESET;

  return(\%partitions);
}


#############################################################################################################################################
#Get expression data for all data types for all libraries specified in the input source directories                                         #
#############################################################################################################################################
sub getExpressionData{
  my %args = @_;
  my $paths_ref = $args{'-paths'};
  my $a_ref = $args{'-a_ref'};
  my $version = $args{'-ensembl_version'};
  my $import_option = $args{'-import_option'};
  my $types_list_ref = $args{'-types_list'};

  my %expression_data;

  print MAGENTA, "\n\nGetting EXPRESSION DATA for ensembl v$version analysis", RESET;

  foreach my $library (sort {$a cmp $b} keys %{$paths_ref}){
    my %data_dirs;
    my %library_expression_data;
    my $data_dir = $paths_ref->{$library}->{data_path};
    my $name = $paths_ref->{$library}->{name};
    my $line_order = $paths_ref->{$library}->{line_order};

    print BLUE, "\n\n\tImporting DATA for library: $library ($name $data_dir)", RESET;

    #1.) Get all available data directories within the specified path (must be a directory with v$version in the name)
    opendir(DIRHANDLE, "$data_dir") || die "\nCannot open directory: $data_dir\n\n";
    my @temp = readdir(DIRHANDLE);
    closedir(DIRHANDLE);

    foreach my $file (@temp){
      my $path = "$data_dir"."$file"."/";
      unless (-d $path){
        next();
      }
      if ($file =~ /v$version/){
        $data_dirs{$file}{path} = "$path";
      }else{
        next();
      }
    }

    my %data_files;
    #2.) Now get the actual expression files
    foreach my $data_dir (sort {$a cmp $b} keys %data_dirs){
      my $path =  $data_dirs{$data_dir}{path};
      my $summary_dir = "$path"."Summary/";

      #Get all available data files (files must contain 'Expresion' in their name)
      opendir(DIRHANDLE, "$summary_dir") || die "\nCannot open directory: $summary_dir\n\n";
      my @temp = readdir(DIRHANDLE);
      closedir(DIRHANDLE);

      foreach my $file (@temp){
        my $path = "$summary_dir"."$file";
        if (-d $path){
          next();
        }
        if ($file =~ /.+\_(\w+)Expression\_/){
          $data_files{$1}{path} = $path;
        }else{
          #print YELLOW, "\n\tFile: $file does not appear to be an expression file", RESET;
          next();
        }
      }
    }

    #3.) Now go through each data file and store data for that data type
    foreach my $data_type (sort {$a cmp $b} keys %data_files){
      unless($types_list_ref->{$data_type}){
        print YELLOW, "\n\t\tExcluding data_type: $data_type from import", RESET;
        next();
      }

      my $seq_a_ref = $a_ref->{$data_type};

      print BLUE, "\n\t\tProcessing data type: $data_type", RESET;
      my $file = $data_files{$data_type}{path};
      my %temp;
      my %columns;
      open(DATA, "$file") || die "\nCould not open file: $file\n\n";

      my $header = 1;
      while(<DATA>){
        chomp($_);
        my @line = split("\t", $_);
        if ($header == 1){
          my $column_count = 0;
          foreach my $column (@line){
            $columns{$column}{column_pos} = $column_count;
            $column_count++;
          }
          $header = 0;
          next();
        }
        #Get data (column names should be common to all expression data files)
        my $id = $line[0];

        #If this feature is not defined in the annotation object - skip it
        unless ($seq_a_ref->{$id}){
          next();
        }

        my $norm_expression = $line[$columns{'Average_Coverage_NORM1'}{column_pos}];
        unless ($norm_expression eq "NA"){
          $norm_expression = sprintf("%.5f", $norm_expression);
        }
        my $expressed = $line[$columns{'Expressed'}{column_pos}];
        my $percent_coverage_1x = $line[$columns{'Percent_Coverage_1x'}{column_pos}];
        my $percent_coverage_10x = "NA";
        if ($columns{'Percent_Coverage_10x'}){
          $percent_coverage_10x = $line[$columns{'Percent_Coverage_10x'}{column_pos}];
        }
        my $percent_gene_expression = "NA";
        if ($columns{'Percent_Gene_Expression'}){
          $percent_gene_expression = $line[$columns{'Percent_Gene_Expression'}{column_pos}];
        }

        #Get read count (may only be defined for gene expression file)
        my $rc = 0;
        my $read_count;
        if ($columns{'Read_Count'}){
          $rc = 1;
          $read_count = $line[$columns{'Read_Count'}{column_pos}];
        }

        #Only import expressed objects
        if ($import_option eq "expressed"){
          if ($expressed eq "1"){
            $temp{$id}{expressed} = $expressed;
            $temp{$id}{norm_expression} = $norm_expression;
            $temp{$id}{percent_coverage_1x} = $percent_coverage_1x;
            $temp{$id}{percent_coverage_10x} = $percent_coverage_10x;
            $temp{$id}{percent_gene_expression} = $percent_gene_expression;
            if ($rc){$temp{$id}{read_count} = $read_count;}
          }
        }elsif ($import_option eq "expressed_known"){
          my $known = $seq_a_ref->{$id}->{known};
          if ($expressed eq "1" || $known eq "1"){
            $temp{$id}{expressed} = $expressed;
            $temp{$id}{norm_expression} = $norm_expression;
            $temp{$id}{percent_coverage_1x} = $percent_coverage_1x;
            $temp{$id}{percent_coverage_10x} = $percent_coverage_10x;
            $temp{$id}{percent_gene_expression} = $percent_gene_expression;
            if ($rc){$temp{$id}{read_count} = $read_count;}
          }
        }elsif($import_option eq "all"){
           $temp{$id}{expressed} = $expressed;
           $temp{$id}{norm_expression} = $norm_expression;
           $temp{$id}{percent_coverage_1x} = $percent_coverage_1x;
           $temp{$id}{percent_coverage_10x} = $percent_coverage_10x;
           $temp{$id}{percent_gene_expression} = $percent_gene_expression;
           if ($rc){$temp{$id}{read_count} = $read_count;}
        }else{
          print RED, "\n\nImport option not understood in getExpressionData()\n\n", RESET;
          exit();
        }
      }

      close(DATA);
      my $data_count = keys %temp;
      print BLUE, "\n\t\t\tFound $data_count expression records", RESET;

      if ($expression_data{$data_type}){
        my $libraries_ref = $expression_data{$data_type}{libraries};
        $libraries_ref->{$library}->{data} = \%temp;
        $libraries_ref->{$library}->{line_order} = $line_order;
      }else{
        my %libraries;
        $libraries{$library}{data} = \%temp;
        $libraries{$library}{line_order} = $line_order;
        $expression_data{$data_type}{libraries} = \%libraries;
      }
    }

    my $stats_dir = $paths_ref->{$library}->{stats_path};
    print BLUE, "\n\n\tImporting JPEGs for library: $library ($name $stats_dir)", RESET;
    #If available get an svg of the distribution of data for this data type, copy the file to an 'images' dir in the webroot and add it to the webpage
    foreach my $data_type (sort {$a cmp $b} keys %data_files){
      unless($types_list_ref->{$data_type}){
        print YELLOW, "\n\t\tExcluding data_type: $data_type from import", RESET;
        next();
      }

      opendir(DIRHANDLE, "$stats_dir") || die "\nCannot open directory: $stats_dir\n\n";
      my @temp = readdir(DIRHANDLE);
      closedir(DIRHANDLE);

      my $string = "hist\_"."$data_type".".svg";
      foreach my $file (@temp){
        if ($file =~ /$string/){
          my $libraries_ref = $expression_data{$data_type}{libraries};
          $libraries_ref->{$library}->{svg_name} = "$file";
          $libraries_ref->{$library}->{svg_path} = "$stats_dir"."$file";

          print BLUE, "\n\t\tFound expression svg: $file", RESET;
        }
      }
    }
  }

  #Determine the unique record count for each data type
  foreach my $data_type (keys %expression_data){
    my $libraries_ref = $expression_data{$data_type}{libraries};
    my %temp;
    foreach my $library (keys %{$libraries_ref}){
      my $data_ref = $libraries_ref->{$library}->{data};
      foreach my $id (keys %{$data_ref}){
        $temp{$id}=1;
      }
    }
    $expression_data{$data_type}{record_count} = keys %temp;
  }

  return(\%expression_data);
}


#############################################################################################################################################
#Get differential expression data for all data types for all comparisons specified in the input source directories                           #
#############################################################################################################################################
sub getDifferentialExpressionData{
  my %args = @_;
  my $paths_ref = $args{'-paths'};
  my $a_ref = $args{'-a_ref'};
  my $version = $args{'-ensembl_version'};
  my $import_option = $args{'-import_option'};
  my $types_list_ref = $args{'-types_list'};

  my %differential_expression_data;

  print MAGENTA, "\n\nGetting DIFFERENTIAL EXPRESSION DATA for ensembl v$version analysis", RESET;
  foreach my $comparison (sort {$a cmp $b} keys %{$paths_ref}){
    my %data_dirs;
    my $data_dir = $paths_ref->{$comparison}->{data_path};
    my $comparison_name = $paths_ref->{$comparison}->{name};
    my $line_order = $paths_ref->{$comparison}->{line_order};

    print BLUE, "\n\n\tImporting DATA for comparison: $comparison ($comparison_name $data_dir)", RESET;

    #1.) Get all available data directories within the specified path (must be a directory with v$version in the name)
    opendir(DIRHANDLE, "$data_dir") || die "\nCannot open directory: $data_dir\n\n";
    my @temp = readdir(DIRHANDLE);
    closedir(DIRHANDLE);

    foreach my $file (@temp){
      my $path = "$data_dir"."$file"."/";
      unless (-d $path){
        next();
      }
      if ($file =~ /v$version/){
        $data_dirs{$file}{path} = "$path";
      }else{
        next();
      }
    }

    my %data_files;
    #2.) Now get the actual DE files
    foreach my $data_dir (sort {$a cmp $b} keys %data_dirs){
      my $path =  $data_dirs{$data_dir}{path};
      my $summary_dir = "$path"."/";

      #Get all available data files (files must contain 'DE' or 'GDE' in their name)
      opendir(DIRHANDLE, "$summary_dir") || die "\nCannot open directory: $summary_dir\n\n";
      my @temp = readdir(DIRHANDLE);
      closedir(DIRHANDLE);

      foreach my $file (@temp){
        my $path = "$summary_dir"."$file";
        if (-d $path){
          next();
        }

        if ($import_option eq "filtered"){
          if ($file =~ /$comparison_name\_(\w+)\_DE\_Values\_Significant\_MTC\.txt/){
            $data_files{$1}{path} = $path;
            $data_files{$1}{DE_type} = "DE";
            $data_files{$1}{string} = "$comparison_name"."\_"."$1"."_DE_Values_Significant_Hist_MTC.svg";
          }elsif($file =~ /$comparison_name\_(\w+)\_GDE\_Values\_Significant\_MTC\.txt/){
            $data_files{$1}{path} = $path;
            $data_files{$1}{DE_type} = "GDE";
            $data_files{$1}{string} = "$comparison_name"."\_"."$1"."_GDE_Values_Significant_Hist_MTC.svg";
          }else{
            next();
          }
        }elsif($import_option eq "noMTC"){
          if ($file =~ /$comparison_name\_(\w+)\_DE\_Values\_Significant\_NoMTC\.txt/){
            $data_files{$1}{path} = $path;
            $data_files{$1}{DE_type} = "DE";
            $data_files{$1}{string} = "$comparison_name"."\_"."$1"."_DE_Values_Significant_Hist_NoMTC.svg";
          }elsif($file =~ /$comparison_name\_(\w+)\_GDE\_Values\_Significant\_NoMTC\.txt/){
            $data_files{$1}{path} = $path;
            $data_files{$1}{DE_type} = "GDE";
            $data_files{$1}{string} = "$comparison_name"."\_"."$1"."_GDE_Values_Significant_Hist_NoMTC.svg";
          }else{
            next();
          }
        }elsif($import_option eq "all"){
           if ($file =~ /$comparison_name\_(\w+)\_DE\_Values\.txt/){
             $data_files{$1}{path} = $path;
             $data_files{$1}{DE_type} = "DE";
             $data_files{$1}{string} = "$comparison_name"."\_"."$1"."_DE_Values_Hist.svg";
           }elsif($file =~ /$comparison_name\_(\w+)\_GDE\_Values\.txt/){
             $data_files{$1}{path} = $path;
             $data_files{$1}{DE_type} = "GDE";
             $data_files{$1}{string} = "$comparison_name"."\_"."$1"."_GDE_Values_Hist.svg";
           }else{
             next();
           }
        }else{
          print RED, "\n\nImport option not understood in getDifferentialExpressionData()\n\n", RESET;
          exit();
        }
      }
    }

    #3.) Now go through each data file and store data for that data type
    foreach my $data_type (sort {$a cmp $b} keys %data_files){
      unless($types_list_ref->{$data_type}){
        print YELLOW, "\n\t\tExcluding data_type: $data_type from import", RESET;
        next();
      }

      my $seq_a_ref = $a_ref->{$data_type};
      print BLUE, "\n\t\tProcessing data type: $data_type", RESET;
      my $file = $data_files{$data_type}{path};
      my $DE_type = $data_files{$data_type}{DE_type};
      my %temp;
      my %columns;
      open(DATA, "$file") || die "\nCould not open file: $file\n\n";

      my $header = 1;
      while(<DATA>){
        chomp($_);
        my @line = split("\t", $_);
        if ($header == 1){
          my $column_count = 0;
          foreach my $column (@line){
            $columns{$column}{column_pos} = $column_count;
            $column_count++;
          }
          $header = 0;
          next();
        }
        #Get data (column names should be common to all expression data files)
        my $id = $line[0];

        #If this feature is not defined in the annotation object - skip it
        unless ($seq_a_ref->{$id}){
          next();
        }
        my $log2_diff = sprintf("%.5f", $line[$columns{'Log2_Diff'}{column_pos}]);
        my $fold_change = sprintf("%.5f", $line[$columns{'Fold_Change'}{column_pos}]);
        my $q_value = sprintf("%.3g", $line[$columns{'BH'}{column_pos}]);

        #If 'All' DE values were imported there will be a lot of 'NA' values where DE could not be assessed (feature was not expressed in one of the two conditions) - skip these
        if($DE_type eq "DE"){
          unless($fold_change eq "NA"){
            $temp{$id}{log2_diff} = $log2_diff;
            $temp{$id}{fold_change} = $fold_change;
            $temp{$id}{q_value} = $q_value;
          }
        }
        #Don't skip NA values for GDE because in many cases all values will be NA
        if($DE_type eq "GDE"){
          $temp{$id}{log2_diff} = $log2_diff;
          $temp{$id}{fold_change} = $fold_change;
          $temp{$id}{q_value} = $q_value;
        }
      }
      close(DATA);
      my $data_count = keys %temp;
      print BLUE, "\n\t\t\tFound $data_count differential expression records", RESET;

      if ($differential_expression_data{$data_type}){
        my $comparisons_ref = $differential_expression_data{$data_type}{comparisons};
        $comparisons_ref->{$comparison}->{data} = \%temp;
        $comparisons_ref->{$comparison}->{line_order} = $line_order;
      }else{
        my %comparisons;
        $comparisons{$comparison}{data} = \%temp;
        $comparisons{$comparison}{line_order} = $line_order;
        $differential_expression_data{$data_type}{comparisons} = \%comparisons;
      }
    }
    my $stats_dir = "$paths_ref->{$comparison}->{stats_path}";
    print BLUE, "\n\n\tImporting JPEGs for comparison: $comparison ($comparison_name $stats_dir)", RESET;
    #If available get a jpg of the distribution of data for this data type, copy the file to an 'images' dir in the webroot and add it to the webpage
    #Also get a jpeg of the correlation between normalized expression values for the two libraries
    foreach my $data_type (sort {$a cmp $b} keys %data_files){
      unless($types_list_ref->{$data_type}){
        print YELLOW, "\n\t\tExcluding data_type: $data_type from import", RESET;
        next();
      }
      my $comparisons_ref = $differential_expression_data{$data_type}{comparisons};
      my $string = $data_files{$data_type}{string};
      my $scan = "ls $stats_dir"."*/*.svg";
      my @temp = `$scan`;
      foreach my $file (@temp){
        chomp($file);
        if ($file =~ /$string/){
          my $name = '';
          if ($file =~ /\/(\w+\.svg)$/){
            $name = $1;
          }
          $comparisons_ref->{$comparison}->{svg_name} = "$name";
          $comparisons_ref->{$comparison}->{svg_path} = "$file";
          print BLUE, "\n\t\tFound differential expression svg: $file", RESET;
        }
      }

      $scan = "ls $stats_dir"."*/*.jpeg";
      @temp = `$scan`;
      my $string1 = "$comparison_name"."\_"."$data_type"."_Expression_Norm_Scatter.jpeg";
      my $string2 = "$comparison_name"."\_"."$data_type"."_Expression_Norm_SmoothScatter.jpeg";
      foreach my $file (@temp){
        chomp($file);
        if ($file =~ /$string1/){
          my $name = '';
          if ($file =~ /\/(\w+\.jpeg)$/){
            $name = $1;
          }
          $comparisons_ref->{$comparison}->{jpeg1_name} = "$name";
          $comparisons_ref->{$comparison}->{jpeg1_path} = "$file";
          print BLUE, "\n\t\tFound expression correlation Scatter jpeg: $file", RESET;
        }
        if ($file =~ /$string2/){
          my $name = '';
          if ($file =~ /\/(\w+\.jpeg)$/){
            $name = $1;
          }
          $comparisons_ref->{$comparison}->{jpeg2_name} = "$name";
          $comparisons_ref->{$comparison}->{jpeg2_path} = "$file";
          print BLUE, "\n\t\tFound expression correlation SmoothScatter jpeg: $file", RESET;
        }
      }
    }
  }

  #Determine the unique record count for each data type
  foreach my $data_type (keys %differential_expression_data){
    my $comparisons_ref = $differential_expression_data{$data_type}{comparisons};
    my %temp;
    foreach my $comparison (keys %{$comparisons_ref}){
      my $data_ref = $comparisons_ref->{$comparison}->{data};
      foreach my $id (keys %{$data_ref}){
        $temp{$id}=1;
      }
    }
    $differential_expression_data{$data_type}{record_count} = keys %temp;
  }

  return(\%differential_expression_data);
}


#############################################################################################################################################
#Get differential splicing data for all data types for all comparisons specified in the input source directories                           #
#############################################################################################################################################
sub getDifferentialSplicingData{
  my %args = @_;
  my $paths_ref = $args{'-paths'};
  my $a_ref = $args{'-a_ref'};
  my $version = $args{'-ensembl_version'};
  my $import_option = $args{'-import_option'};
  my $types_list_ref = $args{'-types_list'};

  my %differential_splicing_data;

  print MAGENTA, "\n\nGetting DIFFERENTIAL SPLICING DATA for ensembl v$version analysis", RESET;

  foreach my $comparison (sort {$a cmp $b} keys %{$paths_ref}){
    my %data_dirs;
    my $data_dir = $paths_ref->{$comparison}->{data_path};
    my $comparison_name = $paths_ref->{$comparison}->{name};
    my $line_order = $paths_ref->{$comparison}->{line_order};

    print BLUE, "\n\n\tImporting DATA for comparison: $comparison ($comparison_name $data_dir)", RESET;

    #1.) Get all available data directories within the specified path (must be a directory with v$version in the name)
    opendir(DIRHANDLE, "$data_dir") || die "\nCannot open directory: $data_dir\n\n";
    my @temp = readdir(DIRHANDLE);
    closedir(DIRHANDLE);

    foreach my $file (@temp){
      my $path = "$data_dir"."$file"."/";
      unless (-d $path){
        next();
      }
      if ($file =~ /v$version/){
        $data_dirs{$file}{path} = "$path";
      }else{
        next();
      }
    }

    my %data_files;
    #2.) Now get the actual expression files
    foreach my $data_dir (sort {$a cmp $b} keys %data_dirs){
      my $path =  $data_dirs{$data_dir}{path};
      my $summary_dir = "$path"."/";

      #Get all available data files (files must contain 'Expresion' in their name)
      opendir(DIRHANDLE, "$summary_dir") || die "\nCannot open directory: $summary_dir\n\n";
      my @temp = readdir(DIRHANDLE);
      closedir(DIRHANDLE);

      foreach my $file (@temp){
        my $path = "$summary_dir"."$file";
        if (-d $path){
          next();
        }

        if ($import_option eq "filtered"){
          if ($file =~ /$comparison_name\_(\w+)\_SI\_Values\_Sorted\_Cutoff\.txt/){
            $data_files{$1}{path} = $path;
            $data_files{$1}{SI_type} = "SI";
          }elsif ($file =~ /$comparison_name\_(\w+)\_GSI\_Values\_Sorted\_Cutoff\.txt/){
            $data_files{$1}{path} = $path;
            $data_files{$1}{SI_type} = "GSI";
          }else{
            #print YELLOW, "\n\tFile: $file does not appear to be an differential expression file", RESET;
            next();
          }
        }elsif($import_option eq "all"){
          if ($file =~ /$comparison_name\_(\w+)\_SI\_Values\.txt/){
            $data_files{$1}{path} = $path;
            $data_files{$1}{SI_type} = "SI";
          }elsif ($file =~ /$comparison_name\_(\w+)\_GSI\_Values\.txt/){
            $data_files{$1}{path} = $path;
            $data_files{$1}{SI_type} = "GSI";
          }else{
            #print YELLOW, "\n\tFile: $file does not appear to be an differential expression file", RESET;
            next();
          }
        }else{
          print RED, "\n\nImport option not understood in getDifferentialSplicingData()\n\n", RESET;
          exit();
        }
      }
    }

    #3.) Now go through each data file and store data for that data type
    foreach my $data_type (sort {$a cmp $b} keys %data_files){
      unless($types_list_ref->{$data_type}){
        print YELLOW, "\n\t\tExcluding data_type: $data_type from import", RESET;
        next();
      }

      my $seq_a_ref = $a_ref->{$data_type};
      print BLUE, "\n\t\tProcessing data type: $data_type", RESET;
      my $file = $data_files{$data_type}{path};
      my $SI_type = $data_files{$data_type}{SI_type};
      my %temp;
      my %columns;
      open(DATA, "$file") || die "\nCould not open file: $file\n\n";

      my $header = 1;
      while(<DATA>){
        chomp($_);
        my @line = split("\t", $_);
        if ($header == 1){
          my $column_count = 0;
          foreach my $column (@line){
            $columns{$column}{column_pos} = $column_count;
            $column_count++;
          }
          $header = 0;
          next();
        }
        #Get data (column names should be common to all expression data files)
        my $id = $line[0];
        #If this feature is not defined in the annotation object - skip it
        unless ($seq_a_ref->{$id}){
          next();
        }

        my $gene_fold_change = sprintf("%.5f", $line[$columns{'GENE_Fold_Change'}{column_pos}]);
        my $seq_fold_change = sprintf("%.5f", $line[$columns{'SEQ_Fold_Change'}{column_pos}]);
        my $si = sprintf("%.5f", $line[$columns{'SI'}{column_pos}]);
        my $reciprocal = $line[$columns{'Reciprocal'}{column_pos}];
        my $reciprocity = sprintf("%.5f", $line[$columns{'Reciprocity'}{column_pos}]);
        unless($reciprocal eq "1"){
          $reciprocity = "NA";
        }
        my $percent_seq_log2_de = sprintf("%.5f", $line[$columns{'percent_SEQ_Log2_DE'}{column_pos}]);
        my $q_value = sprintf("%.5f", $line[$columns{'BH'}{column_pos}]);

        #If 'All' SI values were imported there will be a lot of 'NA' values where DE could not be assessed (feature was not expressed in one of the two conditions) - skip
        if($SI_type eq "SI"){
          unless($si eq "NA"){
            $temp{$id}{gene_fold_change} = $gene_fold_change;
            $temp{$id}{seq_fold_change} = $seq_fold_change;
            $temp{$id}{si} = $si;
            $temp{$id}{reciprocal} = $reciprocal;
            $temp{$id}{reciprocity} = $reciprocity;
            $temp{$id}{percent_seq_log2_de} = $percent_seq_log2_de;
            $temp{$id}{q_value} = $q_value;
          }
        }
        #Don't skip NA values for GSI because in many cases all values will be NA
        if($SI_type eq "GSI"){
            $temp{$id}{gene_fold_change} = $gene_fold_change;
            $temp{$id}{seq_fold_change} = $seq_fold_change;
            $temp{$id}{si} = $si;
            $temp{$id}{reciprocal} = $reciprocal;
            $temp{$id}{reciprocity} = $reciprocity;
            $temp{$id}{percent_seq_log2_de} = $percent_seq_log2_de;
            $temp{$id}{q_value} = $q_value;
        }
      }
      close(DATA);
      my $data_count = keys %temp;
      print BLUE, "\n\t\t\tFound $data_count differential splicing records", RESET;

      my $data_type_comp = "$data_type"."_"."$comparison";
      $differential_splicing_data{$data_type_comp}{data} = \%temp;
      $differential_splicing_data{$data_type_comp}{data_type} = $data_type;
      $differential_splicing_data{$data_type_comp}{comparison} = $comparison;
      $differential_splicing_data{$data_type_comp}{line_order} = $line_order;
      $differential_splicing_data{$data_type_comp}{record_count} = keys %temp;
    }
    my $stats_dir = "$paths_ref->{$comparison}->{stats_path}";
    print BLUE, "\n\n\tImporting JPEGs for comparison: $comparison ($comparison_name $stats_dir)", RESET;
    #If available get a jpg of the distribution of data for this data type, copy the file to an 'images' dir in the webroot and add it to the webpage
    foreach my $data_type (sort {$a cmp $b} keys %data_files){
      unless($types_list_ref->{$data_type}){
        print YELLOW, "\n\t\tExcluding data_type: $data_type from import", RESET;
        next();
      }

      my $data_type_comp = "$data_type"."_"."$comparison";

      my $scan = "ls $stats_dir"."*/*.svg";
      my @temp = `$scan`;

      my $string = "$comparison_name"."_"."$data_type"."_G\*SI_Values_Cutoff_Hist.svg";
      foreach my $file (@temp){
        chomp($file);
        if ($file =~ /$string/){
          my $name = '';
          if ($file =~ /\/(\w+\.svg)$/){
            $name = $1;
          }
          $differential_splicing_data{$data_type_comp}{svg_name} = "$name";
          $differential_splicing_data{$data_type_comp}{svg_path} = "$file";

          print BLUE, "\n\t\tFound differential splicing svg: $file", RESET;
        }
      }
    }
  }

  return(\%differential_splicing_data);
}


#############################################################################################################################################
#Rank data in a hash reference using the specified value                                                                                    #
#############################################################################################################################################
sub rankData{
  my %args = @_;
  my $data_ref = $args{'-data_ref'};
  my $rank_variable = $args{'-rank_variable'};
  my $descending = $args{'-descending'};
  my $absolute = $args{'-absolute'};

  #A.) Assign ranks descending (i.e. highest thing gets rank 1 etc.)
  my %temp;
  my $max = 0;
  my $min = 100000000000000000000000000000000000000000000000000;

  #Convert NA's to 0
  foreach my $id (keys %{$data_ref}){
    if ($data_ref->{$id}->{$rank_variable} eq "NA"){
      $temp{$id}{data} = 0;
    }else{

      if ($absolute){
        $temp{$id}{data} = abs($data_ref->{$id}->{$rank_variable});
      }else{
        $temp{$id}{data} = $data_ref->{$id}->{$rank_variable};
      }

      if ($data_ref->{$id}->{$rank_variable} > $max){
        $max = $data_ref->{$id}->{$rank_variable};
      }
      if ($data_ref->{$id}->{$rank_variable} < $min){
        $min = $data_ref->{$id}->{$rank_variable};
      }
    }
  }

  #Assign ranks, ties will have the same rank
  my $rank = 1;

  if ($descending){
    #DECENDING
    my $current_data = $max;
    foreach my $id (sort {$temp{$b}->{data} <=> $temp{$a}->{data}} keys %temp){
      unless($temp{$id}{data} == $current_data){
        $rank++;
      }
      $current_data = $temp{$id}{data}; 
      $temp{$id}{rank} = $rank;

    }
  }else{
    #ASCENDING
    my $current_data = $min;
    foreach my $id (sort {$temp{$a}->{data} <=> $temp{$b}->{data}} keys %temp){
      unless($temp{$id}{data} == $current_data){
        $rank++;
      }
      $current_data = $temp{$id}{data}; 
      $temp{$id}{rank} = $rank;
    }
  }

  my $assigned_ranks = keys %temp;

  #NA's will all have the same highest rank
  $rank++;
  foreach my $id (keys %{$data_ref}){
    if ($temp{$id}){
      $data_ref->{$id}->{rank} = $temp{$id}{rank};
    }else{
      $data_ref->{$id}->{rank} = $rank;
    }
  }

  #Calculate rank percents (% of rank versus max rank assigned)
  foreach my $id (keys %{$data_ref}){
    $data_ref->{$id}->{rank_percent} = sprintf("%.2f", (($data_ref->{$id}->{rank}/$rank)*100));
  }

  return($assigned_ranks);
}


#############################################################################################################################################
#Assign each sequence record to a genome partition by examining its coordinates                                                             #
#############################################################################################################################################
sub assignPartitions{
  my %args = @_;
  my $a_ref = $args{'-a_ref'};
  my $partitions_ref = $args{'-partitions_ref'};
  my $assigned_count = 0;

  print BLUE, "\n\nAssigning features to genome partitions", RESET;

  foreach my $data_type (sort {$a cmp $b} keys %{$a_ref}){
    my $seq_a_ref = $a_ref->{$data_type};
    print BLUE, "\n\tProcessing: $data_type", RESET;

    foreach my $id (keys %{$seq_a_ref}){
      my $chromosome = $seq_a_ref->{$id}->{chromosome};
      my $seq_start = $seq_a_ref->{$id}->{chr_start};
      my $seq_end = $seq_a_ref->{$id}->{chr_end};

      $seq_a_ref->{$id}->{partition} = "NA";

      if ($partitions_ref->{$chromosome}){
        my $part_ref = $partitions_ref->{$chromosome}->{partitions};
        foreach my $part (keys %{$part_ref}){
          my $part_start = $part_ref->{$part}->{start};
          my $part_end = $part_ref->{$part}->{end};

          if (($seq_start >= $part_start && $seq_start <= $part_end) && ($seq_end >= $part_start && $seq_end <= $part_end)){
            $seq_a_ref->{$id}->{partition} = "chr"."$chromosome"."_"."$part_ref->{$part}->{region}";
            $assigned_count++;
            last();
          }
        }
      }

    }
  }

  return($assigned_count);
}


#############################################################################################################################################
#Generate HTML content for expression data for the top N% of expressed features                                                             #
#############################################################################################################################################
sub generateExpressionContent{
  my %args = @_;
  my $data_type = $args{'-data_type'};
  my $libraries_ref = $args{'-data_ref'};  #Library data for a single data type
  my $a_ref = $args{'-a_ref'};
  my $top_n_percent = $args{'-top_n_percent'};
  my $paths_ref = $args{'-paths'};
  my $web_images_dir = $args{'-web_images_dir'};
  my $track_url = $args{'-track_url'};
  my $ucsc_build = $args{'-ucsc_build'};
  my $matrix_file = $args{'-matrix_file'};
  my $matrix_file_name = $args{'-matrix_file_name'};

  my %print_values;

  my $content = '';

  #Determine the cutoff point
  my $cutoff = 0;
  my $max_data_count = 0;
  foreach my $library (keys %{$libraries_ref}){
    my $grand_read_count = 0;
    my $data_ref = $libraries_ref->{$library}->{data};
    my $data_count = keys %{$data_ref};
    if ($data_count > $max_data_count){
      $max_data_count = $data_count;
    }
    my $test = sprintf("%.2f", ($data_count*($top_n_percent/100)));
    if ($test > $cutoff){
      $cutoff= $test;
    }
    if ($data_type eq "Gene"){
      foreach my $id (keys %{$data_ref}){
        $grand_read_count += $data_ref->{$id}->{read_count};
      }
      $libraries_ref->{$library}->{grand_read_count} = $grand_read_count;
    }
  }
  $cutoff = sprintf("%.0f", $cutoff);
  #Make sure at least 100 values are shown
  if ($cutoff < 100){
    if ($max_data_count < 100){
      $cutoff = $max_data_count;
    }else{
      $cutoff = 100;
    }
  }

  print BLUE, "\n\t\tGenerating HTML to display $cutoff values (top $top_n_percent%)", RESET;

  #Add link to matrix file (feature x library) download
  unless ($matrix_file eq "NA"){
    $content .= "\n<P CLASS=\"Indented12LR_s16_bold\">Download data matrix file ($data_type x library): <A HREF=\"$matrix_file\">$matrix_file_name</A></P><BR>\n";
  }

  #Get annotation refs
  my $gene_a_ref = $a_ref->{'Gene'};
  my $seq_a_ref = $a_ref->{$data_type};

  #Build rows of data for each expression record ordered by rank and for only the top N % of expression records
  #Generate links for each seq to its Gene record page and to data in the custom UCSC track displaying the position corresponding to the record
  foreach my $library (sort {$libraries_ref->{$a}->{line_order} <=> $libraries_ref->{$b}->{line_order}} keys %{$libraries_ref}){
    my $data_ref = $libraries_ref->{$library}->{data};

    my $r=0;
    foreach my $id (sort {$data_ref->{$a}->{rank} <=> $data_ref->{$b}->{rank}} keys %{$data_ref}){
      $r++;
      if ($r > $cutoff){
        last();
      }
      #If this feature is not known and not supported by EST/mRNA (i.e. it appears NOVEL) - the row will be highlighted
      my $class = "Data1";
      if ($seq_a_ref->{$id}->{supporting_seqs} eq "0" && $seq_a_ref->{$id}->{known} eq "0"){
        $class = "Data3";
      }
      my $file_name_prefix = $seq_a_ref->{$id}->{partition};
      my $gene_id_list = $seq_a_ref->{$id}->{gene_id};
      my @gene_ids = split(" ", $gene_id_list);
      my @gene_names;
      my @gene_links;
      foreach my $gene_id (@gene_ids){
        if ($gene_a_ref->{$gene_id}){
          my $gene_name = $gene_a_ref->{$gene_id}->{gene_name};
          my $gene_link = "genes/"."$file_name_prefix"."/"."$gene_a_ref->{$gene_id}->{ensembl_g_id}".".htm";
          push(@gene_names, $gene_name);
          push(@gene_links, $gene_link);
        }
      }
      my $gene_name_link = '';
      foreach my $gene_name (@gene_names){
        my $gene_link = shift(@gene_links);
        $gene_name_link = "$gene_name_link"."<A HREF=\"$gene_link\">$gene_name</A> ";
      }
      my $seq_name = $seq_a_ref->{$id}->{seq_name};
  
      #Convert to log2 scale (first add 1 - to prevent NaN and -ve values)
      my $norm_expression;
      unless ($data_ref->{$id}->{norm_expression} eq "NA"){
        $norm_expression = sprintf("%.2f", ((log(($data_ref->{$id}->{norm_expression}) + 1))/(log(2))));
      }
      my $percent_gene_expression = $data_ref->{$id}->{percent_gene_expression};

      #For Genes only, determine the percentage of all library reads utilized by this gene
      if ($data_type eq "Gene"){
        my $gene_read_count = $data_ref->{$id}->{read_count};
        my $grand_read_count = $libraries_ref->{$library}->{grand_read_count};
        my $percent = sprintf("%.2f", (($gene_read_count/$grand_read_count)*100));
        $percent_gene_expression = $percent;
      }

      my $display_start = ($seq_a_ref->{$id}->{chr_start})-100;
      my $display_end = ($seq_a_ref->{$id}->{chr_end})+100;
      my $chromosome = $seq_a_ref->{$id}->{chromosome};
      if ($chromosome eq "MT"){$chromosome = "M";}

      #Public UCSC
      my $ucsc_link = "http://genome.ucsc.edu/cgi-bin/hgTracks?db="."$ucsc_build"."&position=chr$chromosome:$display_start-$display_end&hgt.customText=$track_url"."$library/combined/"."$file_name_prefix"."_merged.txt.gz&ctfile_"."$ucsc_build"."=";

      my $seq_name_link = $seq_name;
      if ($chromosome =~ /^\d$|^\d{2}$|^X$|^Y$|^M$/){
        $seq_name_link = "<A HREF=\"$ucsc_link\" TARGET=\"_blank\">$seq_name</A> ";

        #If this is a junction include the exons-skipped with the seq name
        if ($data_type =~ /Junction/){
          my $exons_skipped = $seq_a_ref->{$id}->{exons_skipped};
          $seq_name_link = "<A HREF=\"$ucsc_link\" TARGET=\"_blank\">$seq_name</A> (S$exons_skipped) ";
        }
      } 

      my $record = "<TD CLASS=\"$class\">$gene_name_link | $seq_name_link | $norm_expression | $percent_gene_expression%</TD>";

      if ($print_values{$r}){
        my $tmp = "$print_values{$r}{record}"." $record";
        $print_values{$r}{record} = $tmp;
      }else{
        $print_values{$r}{record} = "<TD CLASS=\"Data2\">$r</TD>"."$record";
      }
    }
  }

  #Grab the jpg of the distribution of data for this data type, copy the file to an 'images' dir in the webroot and add it to the webpage
  my $div_count = 0;
  foreach my $library (sort {$libraries_ref->{$a}->{line_order} <=> $libraries_ref->{$b}->{line_order}} keys %{$libraries_ref}){
    $div_count++;
    my $legend_div_name = "leg"."$div_count"."h";
    my $figure_div_name = "box"."$div_count"."h";
    my $library_name = $paths_ref->{$library}->{name};
    my $svg_name = $libraries_ref->{$library}->{svg_name};
    my $svg_path = $libraries_ref->{$library}->{svg_path};
    if ($svg_name){
      print BLUE, "\n\t\tProcessing image: $svg_name", RESET;

      #Copy this file to the images folder in the webroot
      my $cmd_cp = "cp -f $svg_path $web_images_dir"."$library/";
      system($cmd_cp);
      my $image_path = "images/$library/"."$svg_name";
      my $image_desc = "Distribution of expression values for library: $library_name and data type: $data_type";

      my $image_width=700;
      my $image_height=700;

      #Figure title
      $content .= "\n<P CLASS=\"Indented12LR_s16_bold\">Distribution of log2 expression values for library: $library_name and data type: $data_type</P>\n";
      $content .= "\n<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[$figure_div_name]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
      $content .= "<div id=\"$figure_div_name\">";

      #Figure legend block
      $content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$legend_div_name]\" data-openimage=\"images/QuestionOff_icon.gif\" data-closedimage=\"images/QuestionOn_icon.gif\"><img src=\"images/QuestionOn_icon.gif\" border=\"0\" /></a></P>";
      $content .= "\n<div id=\"$legend_div_name\" style=\"width: $image_width"."px;\">";
      $content .= "\n<P CLASS=\"Indented24LR_s16\">Distribution of all expression values for the feature type: $data_type.  Note that unlike in the individual gene pages these log2 expression values were not adjusted by adding 1 before converting to log2 scale.  This results in -ve expression values if the average coverage of the feature was less than 1 before conversion to log2 scale.  The 95th percentile of Silent Intergenic regions is reported to give an indication of the level of background expression noise.  *If you can not see the figure below, click <a href=\"../svg_help.htm\">here</a></P>";
      $content .= "\n</div>";

      #Figure block
      $content .= "\n<P CLASS=\"Indented24LR\"><OBJECT DATA=\"$image_path\" TYPE=\"image/svg+xml\" WIDTH=\"$image_width\" HEIGHT=\"$image_height\" TITLE=\"$image_desc\" CLASS=\"Pic_bordered\" STANDBY=\"$image_desc\"><EMBED SRC=\"$image_path\" TYPE=\"image/svg+xml\" WIDTH=\"$image_width\" HEIGHT=\"$image_height\" /></OBJECT></P><BR>";
      $content .= "\n<BR></div><BR>";

    }
  }

  #Print table header html
  $div_count++;
  my $legend_div_name = "leg"."$div_count"."s";
  $content .= "\n<P CLASS=\"Indented12LR_s16_bold\">Top 1% of expressed $data_type features</P>";
  $content .= "\n<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[$legend_div_name]\" data-openimage=\"images/QuestionOff_icon.gif\" data-closedimage=\"images/QuestionOn_icon.gif\"><img src=\"images/QuestionOff_icon.gif\" border=\"0\" /></a></P>";
  $content .= "\n<div id=\"$legend_div_name\" style=\"width: 800px;\">";
  $content .= "\n<P CLASS=\"Indented12LR_s16\">The following table provides a ranked list of the most highly expressed $data_type features.  Each column corresponds to an expression library.  Each table cell contains the Gene Name (which links to the ALEXA-Seq gene record), followed by the Feature Name (which links to the feature's coordinates in the UCSC Genome Browser and displays expression data), the log2 expression value 'Exp' for the feature, and the percent expression value '% Exp' ([feature expression/gene expression]*100).  For Gene records the percent displayed here represents the percentage of all reads for the library instead.  Features that are NOT supported by an EST/mRNA sequence are indicated in bold.  For exon junction features the number of exons skipped by the junction is indicated as 'Sn', where n is the number of exons skipped (e.g. S0 means no exons skipped, S1 means one exon skipped, etc.).</P><BR>";
  $content .= "\n</div>";
  $content .= "\n<TABLE CLASS=\"Data\">\n  <TR>\n    <TD CLASS=\"Head3\">RANK</TD>";
  foreach my $library (sort {$libraries_ref->{$a}->{line_order} <=> $libraries_ref->{$b}->{line_order}} keys %{$libraries_ref}){
    my $library_name = $paths_ref->{$library}->{name};
    $content .= "<TD CLASS=\"Head1\">$library_name (Gene | Feature | Exp | % Exp)</TD>";
  }
  $content .= "</TR>\n";

  #Print table data (each column corresponds to a single library
  foreach my $pv (sort {$a <=> $b} keys %print_values){
    $content .= "      <TR>$print_values{$pv}{record}</TR>\n";
  }
  $content .= "\n    </TD>\n  </TR>\n</TABLE><BR>";

  return(\$content);
}


#############################################################################################################################################
#Generate HTML content for differential expression data for all significant DE features                                                     #
#############################################################################################################################################
sub generateDifferentialExpressionContent{
  my %args = @_;
  my $data_type = $args{'-data_type'};
  my $comparisons_ref = $args{'-data_ref'};  #DE comparison data for a single data type
  my $a_ref = $args{'-a_ref'};
  my $exp_paths_ref = $args{'-exp_paths'};
  my $de_paths_ref = $args{'-de_paths'};
  my $web_images_dir = $args{'-web_images_dir'};
  my $track_url = $args{'-track_url'};
  my $ucsc_build = $args{'-ucsc_build'};
  my $content_type = $args{'-content_type'};
  my $div_count_start = $args{'-div_count_start'};
  my %print_values;

  my $content = '';

  #Determine the largest number of values to display
  my $max_values = 0;
  foreach my $comparison (keys %{$comparisons_ref}){
    my $data_ref = $comparisons_ref->{$comparison}->{data};
    my $data_count = keys %{$data_ref};
    if ($data_count > $max_values){
      $max_values = $data_count;
    }
  }

  #For performance reasons, do not display more than $max_rows rows in a table
  my $max_rows = 2500;
  if ($max_values > $max_rows){
    $max_values = $max_rows;
  }

  print BLUE, "\n\t\tGenerating HTML to display $max_values values (all significant DE events) or the top $max_rows, whichever is smaller", RESET;

  #Generate the first column
  my $r = 0;
  while ($r < $max_values){
    $r++;
    $print_values{$r}{record} = "<TD CLASS=\"Data2\">$r</TD>";
  }


  #Get annotation refs
  my $gene_a_ref = $a_ref->{'Gene'};
  my $seq_a_ref = $a_ref->{$data_type};

  #Build rows of data for each expression record ordered by rank and for only the top N % of expression records
  #Generate links for each seq to its Gene record page and to data in the custom UCSC track displaying the position corresponding to the record
  foreach my $comparison (sort {$comparisons_ref->{$a}->{line_order} <=> $comparisons_ref->{$b}->{line_order}} keys %{$comparisons_ref}){
    my $data_ref = $comparisons_ref->{$comparison}->{data};
    my $de_paths_ref = $args{'-de_paths'};

    my $r=0;
    foreach my $id (sort {$data_ref->{$a}->{rank} <=> $data_ref->{$b}->{rank}} keys %{$data_ref}){
      $r++;

      if ($r > $max_values){
        last();
      }

      #If this feature is not known and not supported by EST/mRNA (i.e. it appears NOVEL) - the row will be highlighted
      my $class = "Data1";
      if ($seq_a_ref->{$id}->{supporting_seqs} eq "0" && $seq_a_ref->{$id}->{known} eq "0"){
        $class = "Data3";
      }
      my $file_name_prefix = $seq_a_ref->{$id}->{partition};
      my $gene_id_list = $seq_a_ref->{$id}->{gene_id};
      my @gene_ids = split(" ", $gene_id_list);
      my @gene_names;
      my @gene_links;
      foreach my $gene_id (@gene_ids){
        if ($gene_a_ref->{$gene_id}){
          my $gene_name = $gene_a_ref->{$gene_id}->{gene_name};
          my $gene_link = "genes/"."$file_name_prefix"."/"."$gene_a_ref->{$gene_id}->{ensembl_g_id}".".htm";
          push(@gene_names, $gene_name);
          push(@gene_links, $gene_link);
        }
      }
      my $gene_name_link = '';
      foreach my $gene_name (@gene_names){
        my $gene_link = shift(@gene_links);
        $gene_name_link = "$gene_name_link"."<A HREF=\"$gene_link\">$gene_name</A> ";
      }

      my $seq_name = $seq_a_ref->{$id}->{seq_name};
  
      my $fold_change = sprintf("%.2f", $data_ref->{$id}->{fold_change});
      my $q_value = sprintf("%.3g", $data_ref->{$id}->{q_value});

      my $display_start = ($seq_a_ref->{$id}->{chr_start})-100;
      my $display_end = ($seq_a_ref->{$id}->{chr_end})+100;
      my $chromosome = $seq_a_ref->{$id}->{chromosome};
      if ($chromosome eq "MT"){$chromosome = "M";}

      #Public UCSC
      #Create a link to the raw data for both libraries in the current comparison
      my $libA = '';
      my $libB = '';
      if ($comparison =~ /(\w+)\_vs\_(\w+)/){
        $libA = $1;
        $libB = $2;
      }

      my $ucsc_link_A = "http://genome.ucsc.edu/cgi-bin/hgTracks?db="."$ucsc_build"."&position=chr$chromosome:$display_start-$display_end&hgt.customText=$track_url"."DE/$libA/$comparison/"."$file_name_prefix"."_DE_merged.txt.gz&ctfile_"."$ucsc_build"."=";
      my $ucsc_link_B = "http://genome.ucsc.edu/cgi-bin/hgTracks?db="."$ucsc_build"."&position=chr$chromosome:$display_start-$display_end&hgt.customText=$track_url"."DE/$libB/$comparison/"."$file_name_prefix"."_DE_merged.txt.gz&ctfile_"."$ucsc_build"."=";

      my ($libA_name, $libB_name, $linkA, $linkB);
      if ($content_type eq "DE"){
      $libA_name = $exp_paths_ref->{$libA}->{name};
      $libB_name = $exp_paths_ref->{$libB}->{name};
      $linkA = "$libA_name";
      $linkB = "$libB_name";
        if ($chromosome =~ /^\d$|^\d{2}$|^X$|^Y$|^M$/){
          $linkA = "<A HREF=\"$ucsc_link_A\" TARGET=\"_blank\">($libA_name)</A>";
          $linkB = "<A HREF=\"$ucsc_link_B\" TARGET=\"_blank\">($libB_name)</A>";
        } 
      }

      if ($content_type eq "GDE"){
        $libA_name = $libA;
        $libB_name = $libB;
        $linkA = "$libA_name";
        $linkB = "$libB_name";
      }

      my $record = "<TD CLASS=\"$class\">$gene_name_link | $seq_name | $linkA | $linkB | FC = $fold_change | q.value = $q_value</TD>";

      #If this is a junction include the exons-skipped with the seq name
      if ($data_type =~ /Junction/){
        my $exons_skipped = $seq_a_ref->{$id}->{exons_skipped};
        $record = "<TD CLASS=\"$class\">$gene_name_link | $seq_name (S$exons_skipped) | $linkA | $linkB | FC = $fold_change | q.value = $q_value</TD>";
      }

      my $tmp = "$print_values{$r}{record}"." $record";
      $print_values{$r}{record} = $tmp;
    }
    #If the max number of rows was not reached add the remainder...
    while ($r < $max_values){
      $r++;
      my $record = "<TD CLASS=\"Data1\">EMPTY</TD>";
      my $tmp = "$print_values{$r}{record}"." $record";
      $print_values{$r}{record} = $tmp;
    }
  }

  #Grab the a jpg/svg images for this data type, copy the files to an 'images' dir in the webroot and add it to the webpage
  my $div_count = $div_count_start;
  foreach my $comparison (sort {$comparisons_ref->{$a}->{line_order} <=> $comparisons_ref->{$b}->{line_order}} keys %{$comparisons_ref}){
    my $comparison_name = $de_paths_ref->{$comparison}->{name};
    my $svg_name = $comparisons_ref->{$comparison}->{svg_name};
    my $svg_path = $comparisons_ref->{$comparison}->{svg_path};
    my $jpeg1_name = $comparisons_ref->{$comparison}->{jpeg1_name};
    my $jpeg1_path = $comparisons_ref->{$comparison}->{jpeg1_path};
    my $jpeg2_name = $comparisons_ref->{$comparison}->{jpeg2_name};
    my $jpeg2_path = $comparisons_ref->{$comparison}->{jpeg2_path};

    $div_count++;
    my $comp_div_name = "box"."$div_count"."h";
    $content .= "\n<P CLASS=\"Indented12LR_s16_bold\">Differential expression plots for comparison: $comparison_name and data type: $data_type</P>";
    $content .= "\n<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[$comp_div_name]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
    $content .= "\n<div id=\"$comp_div_name\">";

    #Content for SVG showing distribution of DE values
    if ($svg_name){
      print BLUE, "\n\t\tProcessing image: $svg_name", RESET;
      
      #Copy this file to the images folder in the webroot
      my $cmd_cp = "cp -f $svg_path $web_images_dir"."$comparison/";
      system($cmd_cp);
      my $image_path = "images/$comparison/"."$svg_name";
      my $image_desc = "Distribution of differential expression values for comparison: $comparison_name and data type: $data_type";

      my $image_width=750;
      my $image_height=750;

      #Figure title
      $content .= "\n<P CLASS=\"Indented12LR_s16_bold\">Distribution of log2 differential expression values for comparison: $comparison_name and data type: $data_type</P>";

      #Figure legend
      $content .= "\n<P CLASS=\"Indented24LR_s16\">Distribution of all differential expression values that meet the p-value and fold-change cutoff for the feature type: $data_type.  The total number of significant DE features, as well as the max and min log2 DE observed are noted in the legend.  *If you can not see the figure below, click <a href=\"../svg_help.htm\">here</a></P>";

      #Figure
      $content .= "\n<P CLASS=\"Indented24LR\"><OBJECT DATA=\"$image_path\" TYPE=\"image/svg+xml\" WIDTH=\"$image_width\" HEIGHT=\"$image_height\" TITLE=\"$image_desc\" CLASS=\"Pic_bordered\" STANDBY=\"$image_desc\"><EMBED SRC=\"$image_path\" TYPE=\"image/svg+xml\" WIDTH=\"$image_width\" HEIGHT=\"$image_height\" /></OBJECT></P><BR>";
    }

    #Content for JPEG showing Scatter plot of expression values
    if ($jpeg1_name){
      print BLUE, "\n\t\tProcessing image: $jpeg1_name", RESET;
      
      #Copy this file to the images folder in the webroot
      my $cmd_cp = "cp -f $jpeg1_path $web_images_dir"."$comparison/";
      system($cmd_cp);
      my $image_path = "images/$comparison/"."$jpeg1_name";
      my $image_desc = "Scatter plot of expression values for comparison: $comparison_name and data type: $data_type";
      my $image_width=750;
      my $image_height=750;
      #Figure title
      $content .= "\n<P CLASS=\"Indented12LR_s16_bold\">Scatter plot of log2 expression values for comparison: $comparison_name and data type: $data_type</P>";
      
      #Figure legend
      $content .= "\n<P CLASS=\"Indented24LR_s16\">Correlation between expression values for all features that are expressed above background in one or both libraries for the feature type: $data_type.  Features that are differentially expressed (meet the p-value and fold-change cutoff) are indicated in magenta.  A linear model is fit to the data and the correlation by Spearman method is reported.  A loess model is also fit to illustrate the trend of the data.</P>";

      #Figure
      $content .= "\n<P CLASS=\"Indented24LR\"><IMG SRC=\"$image_path\" WIDTH=\"$image_width\" HEIGHT=\"$image_height\" TITLE=\"$image_desc\" ALT=\"$image_desc\" CLASS=\"Pic_bordered\"></P><BR>";
    }

    #Content for JPEG showing SmoothScatter plot of expression values
    if ($jpeg2_name){
      print BLUE, "\n\t\tProcessing image: $jpeg2_name", RESET;
      
      #Copy this file to the images folder in the webroot
      my $cmd_cp = "cp -f $jpeg2_path $web_images_dir"."$comparison/";
      system($cmd_cp);
      my $image_path = "images/$comparison/"."$jpeg2_name";
      my $image_desc = "SmoothScatter plot of expression values for comparison: $comparison_name and data type: $data_type";
      my $image_width=750;
      my $image_height=750;
      #Figure title
      $content .= "\n<P CLASS=\"Indented12LR_s16_bold\">SmoothScatter plot of log2 expression values for comparison: $comparison_name and data type: $data_type</P>";
      
      #Figure legend
      $content .= "\n<P CLASS=\"Indented24LR_s16\">Correlation between expression values for all features that are expressed above background in one or both libraries for the feature type: $data_type.  A linear model is fit to the data and the correlation by Spearman method is reported.  A loess model is also fit to illustrate the trend of the data.</P>";
      
      #Figure
      $content .= "\n<P CLASS=\"Indented24LR\"><IMG SRC=\"$image_path\" WIDTH=\"$image_width\" HEIGHT=\"$image_height\" TITLE=\"$image_desc\" ALT=\"$image_desc\" CLASS=\"Pic_bordered\"></P><BR>";
    }
    $content .= "\n<BR></div><BR>";
  }


  #Print table header html
  $div_count++;
  my $table_div_name = "box"."$div_count"."h";
  my $legend_div_name = "leg"."$div_count"."h";

  my $content_type_name;
  if ($content_type eq "DE"){$content_type_name = "pairwise";}
  if ($content_type eq "GDE"){$content_type_name = "groupwise";}

  $content .= "\n<P CLASS=\"Indented12LR_s16_bold\">Significant differentially expressed $data_type features for $content_type_name comparisons</P>";
  $content .= "\n<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[$table_div_name]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $content .= "\n<div id=\"$table_div_name\">";

  $content .= "\n<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[$legend_div_name]\" data-openimage=\"images/QuestionOff_icon.gif\" data-closedimage=\"images/QuestionOn_icon.gif\"><img src=\"images/QuestionOff_icon.gif\" border=\"0\" /></a></P>";
  $content .= "\n<div id=\"$legend_div_name\" style=\"width: 800px;\">";

  $content .= "\n<P CLASS=\"Indented12LR_s16\">The following table provides a ranked list of all significant differentially expressed $data_type features.  Each column corresponds to a pair-wise comparison of two libraries.  Each table cell contains the Gene Name (which links to the ALEXA-Seq gene record), the Feature Name, the name of each library being compared (which links to the feature's coordinates in the UCSC Genome Browser and displays expression data), the fold-change (FC) of the differential expression event, and the multiple testing corrected p-value for the feature. A bold row indicates that the feature is not currently supported by EST or mRNA sequence alignments. For exon junction features the number of exons skipped by the junction is indicated as 'Sn' where n is the number of exons skipped (e.g. S0 means no exons skipped, S1 means one exon skipped, etc.).</P><BR>";
  $content .= "\n</div>";

  $content .= "\n<TABLE CLASS=\"Data\">\n  <TR>\n    <TD CLASS=\"Head3\">RANK</TD>";
  foreach my $comparison (sort {$comparisons_ref->{$a}->{line_order} <=> $comparisons_ref->{$b}->{line_order}} keys %{$comparisons_ref}){
    my $comparison_name = $de_paths_ref->{$comparison}->{name};
    $content .= "<TD CLASS=\"Head1\">$comparison_name (Gene | Feature | Links | Values)</TD>";
  }
  $content .= "</TR>\n";

  #Print table data (each column corresponds to a single library
  foreach my $pv (sort {$a <=> $b} keys %print_values){
    $content .= "      <TR>$print_values{$pv}{record}</TR>\n";
  }
  $content .= "\n    </TD>\n  </TR>\n</TABLE><BR>";
  $content .= "\n</div><BR>";
  return(\$content);
}


#############################################################################################################################################
#Generate HTML content for differential splicing data for all significant SI features                                                       #
#############################################################################################################################################
sub generateDifferentialSplicingContent{
  my %args = @_;
  my $data_type_comp = $args{'-data_type_comp'};
  my $si_data_ref = $args{'-si_data_ref'};  #DE comparison data for a single data type
  my $a_ref = $args{'-a_ref'};
  my $exp_paths_ref = $args{'-exp_paths'};
  my $si_paths_ref = $args{'-si_paths'};
  my $web_images_dir = $args{'-web_images_dir'};
  my $track_url = $args{'-track_url'};
  my @rank_terms = @{$args{'-rank_terms'}};
  my $ucsc_build = $args{'-ucsc_build'};
  my $content_type = $args{'-content_type'};

  my $data_ref = $si_data_ref->{$data_type_comp}->{data};
  my $comparison = $si_data_ref->{$data_type_comp}->{comparison};
  my $data_type = $si_data_ref->{$data_type_comp}->{data_type};

  my %print_values;

  my $content = '';

  #Determine the number of values that will be displayed
  my $max_values = keys %{$data_ref};

  #For performance reasons, do not display more than $max_rows rows in a table
  my $max_rows = 2500;
  if ($max_values > $max_rows){
    $max_values = $max_rows;
  }

  print BLUE, "\n\t\t\tGenerating HTML to display $max_values values (all significant SI events) (comparison: $comparison) (data_type: $data_type)", RESET;

  #Generate the first column
  my $r = 0;
  while ($r < $max_values){
    $r++;
    $print_values{$r}{record} = "<TD CLASS=\"Data2\">$r</TD>";
  }

  #Get annotation refs
  my $gene_a_ref = $a_ref->{'Gene'};
  my $seq_a_ref = $a_ref->{$data_type};

  #Build rows of data for each expression record ordered by rank and for only the top N % of expression records
  #Generate links for each seq to its Gene record page and to data in the custom UCSC track displaying the position corresponding to the record
  $r = 0;
  foreach my $id (sort {$data_ref->{$a}->{rank} <=> $data_ref->{$b}->{rank}} keys %{$data_ref}){
    $r++;

    if ($r > $max_values){
      last();
    }

    #If this feature is not known and not supported by EST/mRNA (i.e. it appears NOVEL) - the row will be highlighted
    my $class = "Data1";
    if ($seq_a_ref->{$id}->{supporting_seqs} eq "0" && $seq_a_ref->{$id}->{known} eq "0"){
      $class = "Data3";
    }

    my $file_name_prefix = $seq_a_ref->{$id}->{partition};
    my $gene_id_list = $seq_a_ref->{$id}->{gene_id};
    my @gene_ids = split(" ", $gene_id_list);
    my @gene_names;
    my @gene_links;
    foreach my $gene_id (@gene_ids){
      my $gene_name = $gene_a_ref->{$gene_id}->{gene_name};
      my $gene_link = "genes/"."$file_name_prefix"."/"."$gene_a_ref->{$gene_id}->{ensembl_g_id}".".htm";
      push(@gene_names, $gene_name);
      push(@gene_links, $gene_link);
    }
    my $gene_name_link = '';
    foreach my $gene_name (@gene_names){
      my $gene_link = shift(@gene_links);
      $gene_name_link = "$gene_name_link"."<A HREF=\"$gene_link\">$gene_name</A> ";
    }

    my $seq_name = $seq_a_ref->{$id}->{seq_name};
  
    my $si = sprintf("%.2f", $data_ref->{$id}->{si});
    my $gene_fold_change = sprintf("%.2f", $data_ref->{$id}->{gene_fold_change});
    my $seq_fold_change = sprintf("%.2f", $data_ref->{$id}->{seq_fold_change});
    my $reciprocity;
    if ($data_ref->{$id}->{reciprocity} eq "NA"){
      $reciprocity = "N/A";
    }else{
      $reciprocity =  sprintf("%.2f", $data_ref->{$id}->{reciprocity});
    }
    my $percent_seq_log2_de = sprintf("%.2f", $data_ref->{$id}->{percent_seq_log2_de});

    my $q_value;
    if ($content_type eq "GSI"){
      $q_value = sprintf("%.3g", $data_ref->{$id}->{q_value});
    }else{
      $q_value = "N/A";
    }

    my $display_start = ($seq_a_ref->{$id}->{chr_start})-100;
    my $display_end = ($seq_a_ref->{$id}->{chr_end})+100;
    my $chromosome = $seq_a_ref->{$id}->{chromosome};
    if ($chromosome eq "MT"){$chromosome = "M";}

    #Public UCSC
    #Create a link to the raw data for both libraries in the current comparison
    my $libA = '';
    my $libB = '';
    if ($comparison =~ /(\w+)\_vs\_(\w+)/){
      $libA = $1;
      $libB = $2;
    }

    my $ucsc_link_A = "http://genome.ucsc.edu/cgi-bin/hgTracks?db="."$ucsc_build"."&position=chr$chromosome:$display_start-$display_end&hgt.customText=$track_url"."DE/$libA/$comparison/"."$file_name_prefix"."_DE_merged.txt.gz&ctfile_"."$ucsc_build"."=";
    my $ucsc_link_B = "http://genome.ucsc.edu/cgi-bin/hgTracks?db="."$ucsc_build"."&position=chr$chromosome:$display_start-$display_end&hgt.customText=$track_url"."DE/$libB/$comparison/"."$file_name_prefix"."_DE_merged.txt.gz&ctfile_"."$ucsc_build"."=";

    my ($libA_name, $libB_name, $linkA, $linkB);
      if ($content_type eq "SI"){
        $libA_name = $exp_paths_ref->{$libA}->{name};
        $libB_name = $exp_paths_ref->{$libB}->{name};
        $linkA = "$libA_name";
        $linkB = "$libB_name";
        if ($chromosome =~ /^\d$|^\d{2}$|^X$|^Y$|^M$/){
          $linkA = "<A HREF=\"$ucsc_link_A\" TARGET=\"_blank\">($libA_name)</A>";
          $linkB = "<A HREF=\"$ucsc_link_B\" TARGET=\"_blank\">($libB_name)</A>";
       } 
    }

    if ($content_type eq "GSI"){
      $libA_name = $libA;
      $libB_name = $libB;
      $linkA = "$libA_name";
      $linkB = "$libB_name";
    }

    my $record = "<TD CLASS=\"$class\">$gene_name_link | $seq_name | $linkA | $linkB</TD><TD CLASS=\"Data2\">$si</TD><TD CLASS=\"Data2\">$gene_fold_change</TD><TD CLASS=\"Data2\">$seq_fold_change</TD><TD CLASS=\"Data2\">$reciprocity</TD><TD CLASS=\"Data2\">$percent_seq_log2_de</TD><TD CLASS=\"Data2\">$q_value</TD>";

    #If this is a junction include the exons-skipped with the seq name
    if ($data_type_comp =~ /Junction/){
       my $exons_skipped = $seq_a_ref->{$id}->{exons_skipped};
       $record = "<TD CLASS=\"$class\">$gene_name_link | $seq_name (S$exons_skipped) | $linkA | $linkB</TD><TD CLASS=\"Data2\">$si</TD><TD CLASS=\"Data2\">$gene_fold_change</TD><TD CLASS=\"Data2\">$seq_fold_change</TD><TD CLASS=\"Data2\">$reciprocity</TD><TD CLASS=\"Data2\">$percent_seq_log2_de</TD><TD CLASS=\"Data2\">$q_value</TD>";
    }

    my $tmp = "$print_values{$r}{record}"." $record";
    $print_values{$r}{record} = $tmp;
  }

  #Grab the a jpg of the distribution of data for this data type, copy the file to an 'images' dir in the webroot and add it to the webpage
  my $comparison_name = $si_paths_ref->{$comparison}->{name};
  my $svg_name = $si_data_ref->{$data_type_comp}->{svg_name};
  my $svg_path = $si_data_ref->{$data_type_comp}->{svg_path};
  if ($svg_name){
    print BLUE, "\n\t\t\tProcessing image: $svg_name", RESET;

    #Copy this file to the images folder in the webroot
    my $cmd_cp = "cp -f $svg_path $web_images_dir"."$comparison/";
    system($cmd_cp);
    my $image_path = "images/$comparison/"."$svg_name";
    my $image_desc = "Distribution of differential splicing values for comparison: $comparison_name and data type: $data_type";

    my $image_width=750;
    my $image_height=750;

    #Figure title
    my $legend_div_name = "leg1h";
    my $figure_div_name = "box1h";

    $content .= "\n<P CLASS=\"Indented12LR_s16_bold\">Distribution of log2 SI values for comparison: $comparison_name and data type: $data_type</P>";
    $content .= "\n<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[$figure_div_name]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
    $content .= "<div id=\"$figure_div_name\">";

    #Figure legend
    $content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$legend_div_name]\" data-openimage=\"images/QuestionOff_icon.gif\" data-closedimage=\"images/QuestionOn_icon.gif\"><img src=\"images/QuestionOn_icon.gif\" border=\"0\" /></a></P>";
    $content .= "\n<div id=\"$legend_div_name\" style=\"width: $image_width"."px;\">";
    $content .= "\n<P CLASS=\"Indented24LR_s16\">Distribution of all splicing index values that meet the SI cutoff (indicated by red dotted lines) for the feature type: $data_type.  The total number of features with SI values exceeding the cutoff, as well as the max and min log2 SI observed are noted in the legend.  *If you can not see the figure below, click <a href=\"../svg_help.htm\">here</a></P>";
    $content .= "\n</div>";

    #Figure
    $content .= "\n<P CLASS=\"Indented24LR\"><OBJECT DATA=\"$image_path\" TYPE=\"image/svg+xml\" WIDTH=\"$image_width\" HEIGHT=\"$image_height\" TITLE=\"$image_desc\" CLASS=\"Pic_bordered\" STANDBY=\"$image_desc\"><EMBED SRC=\"$image_path\" TYPE=\"image/svg+xml\" WIDTH=\"$image_width\" HEIGHT=\"$image_height\" /></OBJECT></P><BR>";
    $content .= "\n<BR></div><BR>";
  }

  #Print table header html
  my ($link1,$link2,$link3,$link4,$link5,$link6);
  if ($content_type eq "SI"){
    $link1 = "SI_"."$data_type_comp"."_"."$rank_terms[0]".".htm";
    $link2 = "SI_"."$data_type_comp"."_"."$rank_terms[1]".".htm";
    $link3 = "SI_"."$data_type_comp"."_"."$rank_terms[2]".".htm";
    $link4 = "SI_"."$data_type_comp"."_"."$rank_terms[3]".".htm";
    $link5 = "SI_"."$data_type_comp"."_"."$rank_terms[4]".".htm";
  }
  if ($content_type eq "GSI"){
    $link1 = "GSI_"."$data_type_comp"."_"."$rank_terms[0]".".htm";
    $link2 = "GSI_"."$data_type_comp"."_"."$rank_terms[1]".".htm";
    $link3 = "GSI_"."$data_type_comp"."_"."$rank_terms[2]".".htm";
    $link4 = "GSI_"."$data_type_comp"."_"."$rank_terms[3]".".htm";
    $link5 = "GSI_"."$data_type_comp"."_"."$rank_terms[4]".".htm";
    $link6 = "GSI_"."$data_type_comp"."_"."$rank_terms[5]".".htm";
  }
  my $col1 = "<TD CLASS=\"Head1\">$comparison_name (Gene | Feature | Links)</TD>";
  my $col2 = "<TD CLASS=\"Head3\"><A HREF=\"$link1\">SI</A></TD>";
  my $col3 = "<TD CLASS=\"Head3\"><A HREF=\"$link2\">Gene FC</A></TD>";
  my $col4 = "<TD CLASS=\"Head3\"><A HREF=\"$link3\">Seq FC</A></TD>";
  my $col5 = "<TD CLASS=\"Head3\"><A HREF=\"$link4\">Reciprocity<br>(RI)</A></TD>";
  my $col6 = "<TD CLASS=\"Head3\"><A HREF=\"$link5\">Percent Seq DE<br>(PFC)</A></TD>";

  my $col7;
  if ($content_type eq "SI"){
  $col7 = "<TD CLASS=\"Head3\">q-value</TD>";
  } 
  if ($content_type eq "GSI"){
  $col7 = "<TD CLASS=\"Head3\"><A HREF=\"$link6\">q-value</A></TD>";
  } 

  $content .= "\n<P CLASS=\"Indented12LR_s16_bold\">Significant alternatively expressed $data_type features</P>";

  $content .= "\n<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[leg2s]\" data-openimage=\"images/QuestionOff_icon.gif\" data-closedimage=\"images/QuestionOn_icon.gif\"><img src=\"images/QuestionOff_icon.gif\" border=\"0\" /></a></P>";
  $content .= "\n<div id=\"leg2s\" style=\"width: 1000px;\">";

  $content .= "\n<P CLASS=\"Indented12LR_s16\">The following table provides a ranked list of alternatively expressed $data_type features for a single pair-wise library comparison.  The first column contains the Gene Name (which links to the ALEXA-Seq gene record), the Feature Name, and the name of each library being compared (which links to the feature's coordinates in the UCSC Genome Browser and displays expression data).  The 'SI' column reports the Splicing Index calculated for the feature.  The 'Gene FC' column reports the Fold-Change calculated for the entire gene to which the feature belongs.  The 'Seq FC' column reports the Fold-Change of the feature itself.  The 'Reciprocity' column reports the Reciprocity Index (RI) for the feature.  The 'Percent Seq DE' column reports the Percent Feature Contribution (PFC) value for the feature.  A description of the purpose and calculation of SI, RI and PFC values can be found in our manuscript.  Briefly, the SI value is a measure of the degree of change in expression of a feature (e.g. an exon) between two conditions relative to the change in expression at the gene level.  The RI value is a measure of the degree to which the change of the feature is reciprocal in direction to that observed for the gene overall.  For example, if an exon is up-regulated but the gene overall is down-regulated, this will give a higher RI value.  The PFC value is a measure of the degree of differential expression of the feature compared to the entire gene.  For example, if a gene is not changed overall between two conditions, but the feature is highly differentially expressed, this will give a higher PFC value.  To sort this table by each of the data values, simply click the column header.  A bold row indicates that the feature is not currently supported by EST or mRNA sequence alignments. For exon junction features the number of exons skipped by the junction is indicated as 'Sn' where n is the number of exons skipped (e.g. S0 means no exons skipped, S1 means one exon skipped, etc.).</P><BR>";
  $content .= "\n</div>";

  $content .= "\n<TABLE CLASS=\"Data\">\n  <TR>\n    <TD CLASS=\"Head3\">RANK</TD>";
  $content .= "$col1$col2$col3$col4$col5$col6$col7";
  $content .= "</TR>\n";

  #Print table data (each column corresponds to a single library
  foreach my $pv (sort {$a <=> $b} keys %print_values){
    $content .= "      <TR>$print_values{$pv}{record}</TR>\n";
  }
  $content .= "\n    </TD>\n  </TR>\n</TABLE><BR>";

  return(\$content);
}


#############################################################################################################################################
#Write and HTML page for a particular expression content object                                                                             #
#############################################################################################################################################
sub writePage{
  my %args = @_;
  my $path = $args{'-path'};
  my $title = $args{'-title'};
  my $content = $args{'-content'};
  my $css_path = $args{'-css_path'};
  my $alexa_home_path = $args{'-alexa_home_path'};
  my $alexa_seq_home_path = $args{'-alexa_seq_home_path'};
  my $summary_path = $args{'-summary_path'};
  my $genes_path = $args{'-genes_path'};
  my $search_path = $args{'-search_path'};
  my $meta_description = $args{'-meta_description'};
  my $meta_keywords = $args{'-meta_keywords'};
  my $google_analytics = $args{'-google_analytics'};
  my $google_analytics_id = $args{'-google_analytics_id'};
  my $collapse_div_script = $args{'-collapse_div_script'};
  my $jquery_min_script = $args{'-jquery_min_script'};
  my $div_count = $args{'-div_count'};

  #If the user specified to use collapsible divs, include the neccessary HEAD code

  #Create a list of generic data and legend boxs - each with a default show/hide status
  my $divs = '';
  for (my $i = 1; $i <= $div_count; $i++){
    my $n1 = "box"."$i"."s";
    my $n2 = "box"."$i"."h";
    my $n3 = "leg"."$i"."s";
    my $n4 = "leg"."$i"."h";
    $divs .= "animatedcollapse.addDiv(\'$n1\', \'fade=1, speed=500, hide=0\')\n"; #Content box - Show
    $divs .= "animatedcollapse.addDiv(\'$n2\', \'fade=1, speed=500, hide=1\')\n"; #Content box - Hide
    $divs .= "animatedcollapse.addDiv(\'$n3\', \'fade=1, speed=500, hide=0\')\n"; #Legend text - Show
    $divs .= "animatedcollapse.addDiv(\'$n4\', \'fade=1, speed=500, hide=1\')\n"; #Legend text - Hide
  }

  my $collapse_div_content = '';
  if ($collapse_div_script && $jquery_min_script){
    $collapse_div_content = <<"EOF";
<script type="text/javascript" src="$jquery_min_script"></script>
<script type="text/javascript" src="$collapse_div_script"></script>
<script type="text/javascript">
$divs
animatedcollapse.ontoggle=function(\$, divobj, state){}
animatedcollapse.init()
</script>
EOF
  }

  #If the user specified the google analytics option, add code to allow tracking of this webpage
  my $ga_content_bottom = '';
  if ($google_analytics){
    if ($google_analytics =~ /UA\-\d+\-\d+/){
    $ga_content_bottom = <<"EOF";

<script type=\"text/javascript\">
var gaJsHost = ((\"https:\" == document.location.protocol) ? \"https://ssl.\" : \"http://www.\");
document.write(unescape(\"%3Cscript src='\" + gaJsHost + \"google-analytics.com/ga.js' type='text/javascript'%3E%3C/script%3E\"));
</script>
<script type=\"text/javascript\">
var pageTracker = _gat._getTracker(\"UA-3266645-1\");
pageTracker._initData();
pageTracker._trackPageview();
</script>

EOF
    }
  }

  open (OUT, ">$path") || die "\nCould not open output file to write: $path ($!)\n\n";

print OUT <<"EOF";
<HTML>
<HEAD>
<TITLE>$title</TITLE>
<META name=\"description\" content=\"$meta_description\">
<META name=\"keywords\" content=\"$meta_keywords\">
<link rel="stylesheet" type="text/css" href="$css_path">
$collapse_div_content
</HEAD>
<BODY>

<!-- Start top navigation section -->
<P CLASS=\"Navigation\"><A HREF=\"$alexa_home_path\">ALEXA</A> | <A HREF=\"$alexa_seq_home_path\">ALEXA-Seq</A> | <A HREF=\"$summary_path\">Summary</A> | <A HREF=\"$genes_path\">Genes</A> | <A HREF=\"$search_path\">SEARCH</A></P><BR>
<!-- End top navigation section -->


<!-- Start content section -->
<BR>
<P CLASS=\"Indented12LR_s19_bold\">$title</P>
<BR>

${$content}

$ga_content_bottom

</BODY>
</HMTL>
<!-- End content section -->


<!-- Start bottom navigation section -->
<P CLASS=\"Navigation\"><A HREF=\"$alexa_home_path\">ALEXA</A> | <A HREF=\"$alexa_seq_home_path\">ALEXA-Seq</A> | <A HREF=\"$summary_path\">Summary</A> | <A HREF=\"$genes_path\">Genes</A> | <A HREF=\"$search_path\">SEARCH</A></P><BR>
<!-- End bottom navigation section -->

EOF

  close(OUT) || die "\nCould not close output file: $path ($!)\n\n";

  return();
}




1;

