#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Create a gene annotation database file for Gene features

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
use utilities::ALEXA_DB qw(:all);
use utilities::utility qw(:all);

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $chr_filter = '';
my $outdir = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'chr_filter=s'=>\$chr_filter, 'outdir=s'=>\$outdir);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify a single chromosome region to be processed using: --chr_filter", RESET;
print GREEN, "\n\tSpecify the path to be used for output files using: --outdir", RESET;

print GREEN, "\n\nExample: createGeneDatabase.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --chr_filter='3:16:121020102-126260001'  --outdir=/projects/malachig/sequence_databases/hs_53_36o/ensembl_genes_hs_53_36o/\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $chr_filter && $outdir){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

print BLUE, "\n\nUser supplied args.  database: $database\tserver: $server\tchr_filter: $chr_filter\toutdir: $outdir", RESET;

my $region_number;
my $start_filter;
my $end_filter;
if ($chr_filter =~ /(.*)\:(\d+)\:(\d+)\-(\d+)/){
  $chr_filter = $1;
  $region_number = $2;
  $start_filter = $3;
  $end_filter = $4;
  unless ($end_filter > $start_filter){
    print RED, "\nStart of range must be smaller than end ($chr_filter)\n\n", RESET;
    exit();
  }
}else{
  print RED, "\nFormat of chr_filter not understood: $chr_filter (should be of the form:  Y:1:1-9999001)\n\n", RESET;
  exit();
}
$outdir = &checkDir('-dir'=>$outdir, '-clear'=>"no");

my $outfile = "$outdir"."Genes_chr$chr_filter"."_"."$region_number".".txt";
open (GENE, ">$outfile") || die "\nCould not open outfile: $outfile\n\n", RESET;

my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

#Get all genes for the specified region
$| = 1; print BLUE, "\nGetting genes for region: $chr_filter:$start_filter-$end_filter", RESET; $| = 0;
my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All', '-chromosome'=>$chr_filter, '-range'=>"$start_filter-$end_filter")};


#Get the gene info for all genes for which reads were found on the current chromosome
$| = 1; print BLUE, "\nGetting gene info", RESET; $| = 0;
my $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

#Get masked gene sequence info
$| = 1; print BLUE, "\nGetting masked gene sequences", RESET; $| = 0;
my $genes_ref_masked = &getMaskedGene ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);

#Get the transcript info for all transcripts of these genes
$| = 1; print BLUE, "\nGetting transcript data as well as exons for each transcript", RESET; $| = 0;
my $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

#print Dumper $gene_transcripts_ref;

#Get chromosome coordinates for all EnsEMBL transcripts and Build UCSC tracks for all ensembl transcripts
$| = 1; print BLUE, "\n\nCalculating chromosome coordinates for the EXONS of each gene (only for genes of the current chromosome)", RESET; $| = 0;

print GENE "Gene_ID\tEnsEMBL_Gene_ID\tGene_Name\tDescription\tGene_Type\tGene_Evidence\tChromosome\tStrand\tUnit1_start_chr\tUnit1_end_chr\tTranscript_Count\tExon_Count\tExon_Content_Count\tBase_Count\tUnMasked_Base_Count\tCoding_Base_Count\tSeq_Name\n";

foreach my $gene_id (keys %{$gene_transcripts_ref}){

  my $chromosome = $genes_ref->{$gene_id}->{chromosome};
  my $chr_strand = $genes_ref->{$gene_id}->{chr_strand};
  my $chr_start = $genes_ref->{$gene_id}->{chr_start};
  my $chr_end = $genes_ref->{$gene_id}->{chr_end};
  my $gene_start = $genes_ref->{$gene_id}->{gene_start};
  my $gene_end = $genes_ref->{$gene_id}->{gene_end};

  #Get protein coding coordinates for this gene
  my $protein_bases_ref = &getProteinBases ('-gene_id'=>$gene_id, '-genes_ref'=>$genes_ref, '-gene_transcripts_ref'=>$gene_transcripts_ref);
  $genes_ref->{$gene_id}->{coding_base_count} = &getProteinBaseCount ('-protein_bases_ref'=>$protein_bases_ref, '-chr_start'=>$chr_start, '-chr_end'=>$chr_end);

  my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};
  $genes_ref->{$gene_id}->{trans_count} = keys %{$transcripts_ref};

  foreach my $trans_id (keys %{$transcripts_ref}){
    my $exons_ref = $transcripts_ref->{$trans_id}->{exons};
    foreach my $exon_id (keys %{$exons_ref}){
      my $start = $exons_ref->{$exon_id}->{exon_start};
      my $end = $exons_ref->{$exon_id}->{exon_end};
      my $coords_ref = &convertGeneCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$start, '-end_pos'=>$end, '-ordered'=>"yes");
      $exons_ref->{$exon_id}->{chr_start} = $coords_ref->{$gene_id}->{chr_start};
      $exons_ref->{$exon_id}->{chr_end} = $coords_ref->{$gene_id}->{chr_end};
      $exons_ref->{$exon_id}->{strand} = $coords_ref->{$gene_id}->{strand};
    }
  }
}

#Get exon content for genes of the chr
$| = 1; print BLUE, "\n\nGetting EXON CONTENT of each gene", RESET; $| = 0;
my $gene_exon_content_ref = &getExonContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);

my $counter = 0;
my $ref_exon_counter = 0;

foreach my $gene_id (@gene_ids){

  my %exon_content_chr;

  my $chromosome = $genes_ref->{$gene_id}->{chromosome};
  my $chr_strand = $genes_ref->{$gene_id}->{chr_strand};
  my $chr_start = $genes_ref->{$gene_id}->{chr_start};
  my $chr_end = $genes_ref->{$gene_id}->{chr_end};
  my $gene_start = $genes_ref->{$gene_id}->{gene_start};
  my $gene_end = $genes_ref->{$gene_id}->{gene_end};

  $counter++;
  if ($counter == 100){
    $counter = 0;
    $| = 1; print BLUE, ".", RESET; $| = 0;
  }

  #Calculate the size of each transcript by adding up the size of its exons
  my $size = 0;
  my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};

  $genes_ref->{$gene_id}->{base_count} = 0;           #Total number of bases covered by exon content of this gene
  $genes_ref->{$gene_id}->{unmasked_base_count} = 0;  #Total number of unmasked bases covered by exon content of this gene

  my $masked_gene_seq = $genes_ref_masked->{$gene_id}->{sequence};

  #Get the chromosome coordinates for exon content blocks
  foreach my $exon_id (sort {$exon_content_ref->{$a}->{start} <=> $exon_content_ref->{$b}->{start}} keys %{$exon_content_ref}){

    my $start = $exon_content_ref->{$exon_id}->{start};
    my $end = $exon_content_ref->{$exon_id}->{end};

    my $coords_ref = &convertGeneCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$start, '-end_pos'=>$end, '-ordered'=>"yes");
    $exon_content_chr{$exon_id}{chr_start} = $coords_ref->{$gene_id}->{chr_start};
    $exon_content_chr{$exon_id}{chr_end} = $coords_ref->{$gene_id}->{chr_end};
    $exon_content_chr{$exon_id}{strand} = $coords_ref->{$gene_id}->{strand};
    my $size = ($coords_ref->{$gene_id}->{chr_end} - $coords_ref->{$gene_id}->{chr_start})+1;
    $exon_content_chr{$exon_id}{size} = $size;
    $genes_ref->{$gene_id}->{base_count} += $size;

    my $masked_region_seq = substr ($masked_gene_seq, $start, $size);
    my $masked_n_count = ($masked_region_seq =~ tr/N/N/);
    my $unmasked_count = $size - $masked_n_count;
    $genes_ref->{$gene_id}->{unmasked_base_count} += $unmasked_count;

    #Store the exon-content coordinates
    $exon_content_ref->{$exon_id}->{chr_start} = $coords_ref->{$gene_id}->{chr_start};
    $exon_content_ref->{$exon_id}->{chr_end} = $coords_ref->{$gene_id}->{chr_end};
  }

  my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};

  #Assemble a reference set of exons (a superset of all non-redundant exons)
  my $trans_count = keys %{$transcripts_ref};
  $gene_transcripts_ref->{$gene_id}->{trans_count} = $trans_count;

  my %reference_exons;
  foreach my $trans_id (keys %{$transcripts_ref}){
    my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

    foreach my $trans_exon_id (keys %{$exons_ref}){
      my $trans_exon_start = $exons_ref->{$trans_exon_id}->{exon_start};
      my $trans_exon_end = $exons_ref->{$trans_exon_id}->{exon_end};

      #Check each of the reference exons to see if one of these is the same, otherwise add it to the list
      my $redundant_exon = 0;
      foreach my $ref_exon_id (keys %reference_exons){
        if ($trans_exon_start == $reference_exons{$ref_exon_id}{exon_start} && $trans_exon_end == $reference_exons{$ref_exon_id}{exon_end}){
          $redundant_exon = 1;
	}
      }
      #Unless the current transcript exon was found to be redundant, add it to the list
      unless ($redundant_exon == 1){
        $reference_exons{$trans_exon_id}{exon_start} = $trans_exon_start;
	$reference_exons{$trans_exon_id}{exon_end} = $trans_exon_end;
      }
    }
  }

  my $exon_count = keys %reference_exons;
  my $exon_content_count = keys %{$exon_content_ref};

  print GENE "$gene_id\t$genes_ref->{$gene_id}->{ensembl_g_id}\t$genes_ref->{$gene_id}->{gene_name}\t$genes_ref->{$gene_id}->{description}\t$genes_ref->{$gene_id}->{gene_type}\t$genes_ref->{$gene_id}->{evidence}\t$genes_ref->{$gene_id}->{chromosome}\t$genes_ref->{$gene_id}->{chr_strand}\t$genes_ref->{$gene_id}->{chr_start}\t$genes_ref->{$gene_id}->{chr_end}\t$trans_count\t$exon_count\t$exon_content_count\t$genes_ref->{$gene_id}->{base_count}\t$genes_ref->{$gene_id}->{unmasked_base_count}\t$genes_ref->{$gene_id}->{coding_base_count}\t$genes_ref->{$gene_id}->{gene_name}\n";

}

$alexa_dbh->disconnect();

close (GENE);

print "\n\nSCRIPT COMPLETE\n\n";

exit();




