#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to define exon regions 
#In the case of a normal exon from a gene without overlapping exons, the exon region is the same as the exon
#In the case of a gene with multiple transcripts that have overlapping exons with different boundaries, the exon regions reflect these

#Basic info to be printed for all exon regions:
#ALEXA Gene ID
#EnsEMBL Gene ID
#Gene Name
#Chromosome
#Start coord (gene based, always 1, 5' to 3') 
#End Coord (gene based, always equals length of gene including introns, 5' to 3')
#Strand
#Chromosome Start Coord (Start position on the chromosome)
#Chromosome End Coord (End position on the chromosome)
#Length (Number of bases in exon region)
#Exon Name (i.e. 1, 2, ,3.   OR in the case of regions extracted from a cluster of overlapping exons, 1a, 1b, 1c, ...
#Number supporting EnsEMBL transcripts
#Number supporting mRNAs
#Number supporting ESTs
#Transcript specific ID (Transcript ID if this region corresponds to a single transcript only)
#- 'S' denotes a region specific to a single transcript where only one transcript exists for the gene anyway
#- 'M' denotes a region specific to a single transcript where multiple transcripts were possible for that gene 

#- For each gene, get all exons and their gene coordinates via ALEXA
#(1.) For exons that are completely isolated from other exons (no overlapping coordinates)
#     - The coordinates will be as is.  Simply convert to chromosome coordinates and print out basic info

#(2.) For all exons that have some overlap, create clusters of these exons.
#     - (b) Identify all boundary points for the exons in the cluster - This will be a non-redundant list of the start/end positions
#           - Treat the regions between these as 'exons' report the coordinates, etc

#Basically the main goal of this script is to identify exon regions that will interogate the expression of individual exons in the most informative way possible

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use BerkeleyDB;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);
use utilities::ucsc qw(:all);

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $chr_filter = '';
my $ucsc_align_dir = '';
my $genbank_mapfile = '';
my $wiggle = '';
my $outfile = '';
my $fasta_file = '';
my $logfile = '';
my $verbose = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
            'chr_filter=s'=>\$chr_filter, 'ucsc_align_dir=s'=>\$ucsc_align_dir, 'genbank_mapfile=s'=>\$genbank_mapfile, 'wiggle=i'=>\$wiggle,
	    'fasta_file=s'=>\$fasta_file, 'outfile=s'=>\$outfile, 'logfile=s'=>\$logfile, 'verbose=s'=>\$verbose);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify a single chromosome to be processed using: --chr_filter", RESET;
print GREEN, "\n\tSpecify the name of the output annotation file using: --outfile", RESET;
print GREEN, "\n\tSpecify the name of an outuput fasta file using: --fasta_file", RESET;
print GREEN, "\n\tSpecify a directory containing UCSC mRNA/EST/xmRNA and xEST files using: --ucsc_align_dir", RESET;
print GREEN, "\n\tSpecify a berkeley DB containing GenBank-To-Species mappings using:  --genbank_mapfile", RESET;
print GREEN, "\n\tJunctions for xmRNA and xEST will be allowed to 'wiggle' by the amount specied by:  --wiggle", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of the design process using: --logfile", RESET;
print GREEN, "\n\tIf you want verbose output, use: --verbose=yes", RESET;

print GREEN, "\n\nExample: createExonRegionDatabase.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --chr_filter='3:16:121020102-126260001'  --ucsc_align_dir=/projects/malachig/sequence_databases/hs_53_36o/mrna_est/partitions/  --genbank_mapfile=/projects/malachig/sequence_databases/hs_53_36o/mrna_est/GenBankToOrganism.btree  --wiggle=3  --outfile=/projects/malachig/sequence_databases/hs_53_36o/ensembl_exonRegions_hs_53_36o/temp/exonRegions_annotated_3_16.txt  --fasta_file=/projects/malachig/sequence_databases/hs_53_36o/ensembl_exonRegions_hs_53_36o/blastdb/temp/exonRegions_3_16.fa  --logfile=/projects/malachig/sequence_databases/hs_53_36o/logs/createExonRegionDatabase/createExonRegionDatabase_3_16_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $chr_filter && $ucsc_align_dir && $genbank_mapfile && ($wiggle =~ /\d+/) && $outfile && $fasta_file && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

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

#Check the temp dir before proceeding.  Make sure it is empty
$ucsc_align_dir = &checkDir('-dir'=>$ucsc_align_dir, '-clear'=>"no");

if ($chr_filter eq "MT"){$chr_filter = "M";}
my $tmp_chr = "chr$chr_filter";

my $ucsc_mrna_table = "$ucsc_align_dir"."$tmp_chr"."_mrna.txt.gz";
my $ucsc_est_table = "$ucsc_align_dir"."$tmp_chr"."_est.txt.gz";
my $ucsc_xmrna_table = "$ucsc_align_dir"."$tmp_chr"."_xmrna.txt.gz";
my $ucsc_xest_table = "$ucsc_align_dir"."$tmp_chr"."_xest.txt.gz";

#Load the GenBank-To-Species mapfile
my %gb_org_map;
tie(%gb_org_map, 'BerkeleyDB::Btree', -Cachesize =>256000000, -Filename=> $genbank_mapfile , -Flags => DB_RDONLY) or die "can't open file $genbank_mapfile: $! $BerkeleyDB::Error\n";

#Global variable storing letters for labeling of multiple exon regions within an exon content block
my @letters = ('a' .. 'z', 'aa' .. 'zz');

#Keep track of the number of exon regions for which probes are designed.  This script will design many probes spanning each of these regions.
#A second script will try to choose the 'best' single probe from each region
my $gene_count = 0;

#Open the exon region and fasta output files
print BLUE, "\nAll data will be written to $outfile\n\n", RESET;
print LOG "\nAll data will be written to $outfile\n\n";

open (OUTFILE, ">$outfile") || die "\nCould not open outfile: $outfile";
print OUTFILE "ExonRegion_ID\tGene_ID\tEnsEMBL_Gene_ID\tGene_Name\tChromosome\tStrand\tUnit1_start\tUnit1_end\tUnit1_start_chr\tUnit1_end_chr\tBase_Count\tUnMasked_Base_Count\tCoding_Base_Count\tSupporting_EnsEMBL_Count\tSupporting_mRNA_Count\tSupporting_EST_Count\tSupporting_xmRNA_Count\tSupporting_xEST_Count\tConserved_Species_Count\tSeq_Name\tSpecific_Trans_ID\n";

open (FASTA, ">$fasta_file") || die "\nCould not open output fasta file: $fasta_file";


#1.) Get basic gene info
#Establish connection with the Alternative Splicing Expression database and get genes for the specified region
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All', '-chromosome'=>$chr_filter, '-range'=>"$start_filter-$end_filter")};
my $genes_ref;
my $gene_transcripts_ref;
my $gene_exon_content_ref;
my $masked_gene_ref;
&getBasicGeneInfo ('-gene_ids'=> \@gene_ids);
$alexa_dbh->disconnect();


#2.) Define clusters of overlapping exons and use these clusters to define a global set of 'Exon Regions'
#    - An 'Exon Region' has the coordinates of a normal exon except in the case where there are overlaping exons from multiple transcripts of a gene with different boundaries
#    - Organize these regions with unique exon region IDs but keyed according to gene ID (hash of hashes)
my $exon_region_id = 0;
my %gene_exon_regions;
&getExonClusters();


#3.) Name the exon regions (1,2,3... from 5' to 3' AND a,b,c... for multiple exon regions mapping to a single exon content block)
#    - Compare each exon regions to the exon content of their gene
&nameExonRegions();


#4.) Calculate chromosome coordinates for all exon regions
print BLUE, "\n\nCalculating chromosome coordinates for each exon region", RESET;
print LOG "\n\nCalculating chromosome coordinates for each exon region";

foreach my $gene_id (sort {$a <=> $b} keys %gene_exon_regions){
  my $exon_region_ref = $gene_exon_regions{$gene_id}{exon_regions};

  foreach my $er_id (sort {$a <=> $b} keys %{$exon_region_ref}){

    my $er_start = $exon_region_ref->{$er_id}->{unit1_start};
    my $er_end = $exon_region_ref->{$er_id}->{unit1_end};

    my $coords_ref = &convertGeneCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$er_start, '-end_pos'=>$er_end, '-ordered'=>"yes");

    $exon_region_ref->{$er_id}->{unit1_start_chr} = $coords_ref->{$gene_id}->{chr_start};
    $exon_region_ref->{$er_id}->{unit1_end_chr} = $coords_ref->{$gene_id}->{chr_end};
  }
}


#5.) Sequence support - Note the number of SEQs which support each exon region (i.e. the exon region is entirely contained within an exon of a known SEQ)
#    - SEQ type: mrna, est, xmrna, xest
my $mrna_exons = &importUcscAlignments('-import_type'=>"exon", '-align_file'=>$ucsc_mrna_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"mrna", '-filter'=>"0", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);
my $est_exons = &importUcscAlignments('-import_type'=>"exon", '-align_file'=>$ucsc_est_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"est", '-filter'=>"1", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);
my $xmrna_exons = &importUcscAlignments('-import_type'=>"exon", '-align_file'=>$ucsc_xmrna_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"xmrna", '-filter'=>"2", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);
my $xest_exons = &importUcscAlignments('-import_type'=>"exon", '-align_file'=>$ucsc_xest_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"xest", '-filter'=>"2", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);

&determineSequenceSupport('-mrna_file'=>$ucsc_mrna_table, '-est_file'=>$ucsc_est_table);


#6.) Determine which exon regions are specific to a single EnsEMBL transcript
#    - Note the number of EnsEMBL transcripts which support each exon region (i.e. the exon region is entirely contained within an exon of the EnsEMBL transcript)
#    - i.e. they are entirely contained within one EnsEMBL transcript and have NO overlap with any other EnsEMBL transcript
#    - Also note the number of EnsEMBL transcripts supporting each exon region (i.e. have an exon that the exon region is completely within)
&determineEnsemblSupport();


#7.) Print out the exon region annotation file
&printExonRegionFile();

untie(%gb_org_map);

print "\n\nSCRIPT COMPLETE\n\n";
print LOG "\n\nSCRIPT COMPLETE\n\n";

#Close the output file
close (OUTFILE);
close(FASTA);
close (LOG);
exit();


################################################################################################
#Create gene object - return hash with gene sequence, masked sequence, exon info, etc.         #
#Identify overlapping exons and define exon clusters                                           #
################################################################################################
sub getBasicGeneInfo{
  my %args = @_;
  my @gene_ids = @{$args{'-gene_ids'}};

  my %gene_object;

  #Get the raw gene sequence
  print BLUE, "\nGetting basic gene info for chr: $chr_filter:$start_filter-$end_filter", RESET;
  print LOG "\nGetting basic gene info for chr: $chr_filter:$start_filter-$end_filter";

  $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids,'-sequence'=>"yes");

  #Get the transcripts for these genes to allow comparison to actual known transcripts
  print BLUE, "\nGetting basic trancript and exon info", RESET;
  print LOG "\nGetting basic transcript and exon info";

  $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

  #Get exon content for all genes
  $gene_exon_content_ref = &getExonContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);

  $masked_gene_ref = &getMaskedGene ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);


  #Process each gene
  foreach my $gene_id (@gene_ids){

    if ($genes_ref->{$gene_id}->{chromosome} eq "MT"){$genes_ref->{$gene_id}->{chromosome} = "M";}

    #Get a list of non-redundant exon for this gene
    my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};

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

    my $nr_exon_count = keys %reference_exons;

    #Store the reference exons object (non-redundant exons) for later use
    $genes_ref->{$gene_id}->{reference_exons} = \%reference_exons;
    $genes_ref->{$gene_id}->{nr_exon_count} = $nr_exon_count;

    #Go through each exon and create clusters of overlapping exons.  Any exon that shares at least one bp with another exon will be placed in the same cluster
    #An exon that is completely isolated from others will comprise a cluster of one
    my $cluster_count = 0;
    my %clusters;

  EXON: foreach my $exon_id (sort {$reference_exons{$a}->{exon_start} <=> $reference_exons{$b}->{exon_start}} keys %reference_exons){

      #Get the start and end positions of the current exon
      my $exon_start = $reference_exons{$exon_id}{exon_start};
      my $exon_end = $reference_exons{$exon_id}{exon_end};

      #Go through each cluster and see if the current exon overlaps with one of the exons already present in one of these clusters
      foreach my $cluster_id (sort keys %clusters){
	my @tmp_exons = @{$clusters{$cluster_id}{exons}};
	foreach my $exon_id_test (@tmp_exons){

	  my $exon_start_test = $reference_exons{$exon_id_test}{exon_start};
	  my $exon_end_test = $reference_exons{$exon_id_test}{exon_end};

	  #See if the start of the current exon is within the range of the test exon
	  if ($exon_start >= $exon_start_test && $exon_start <= $exon_end_test){
	    push (@tmp_exons, $exon_id);
	    $clusters{$cluster_id}{exons} = \@tmp_exons;
	    next(EXON);
	  }
	  #See if the end of the current exon is within the range of the test exon
	  if ($exon_end >= $exon_start_test && $exon_end <= $exon_end_test){
	    push (@tmp_exons, $exon_id);
	    $clusters{$cluster_id}{exons} = \@tmp_exons;
	    next(EXON);
	  }
	  #See if the current exon completely flanks the test exon - if so it should be added to the cluster
	  if ($exon_start <= $exon_start_test && $exon_end >= $exon_end_test){
	    push (@tmp_exons, $exon_id);
	    $clusters{$cluster_id}{exons} = \@tmp_exons;
	    next(EXON);
	  }
	}
      }
      #If the current exon was not added to any of the current clusters - create a new cluster
      $cluster_count++;
      my @tmp_exons;
      push (@tmp_exons, $exon_id);
      $clusters{$cluster_count}{exons} = \@tmp_exons;
    }

    $genes_ref->{$gene_id}->{clusters} = \%clusters;
  }

  return();
}


################################################################################################
#Define clusters of overlaping exons
################################################################################################
sub getExonClusters{

  print BLUE, "\n\nDefining clusters of overlaping exons", RESET;
  print LOG "\n\nDefining clusters of overlaping exons";

  foreach my $gene_id (@gene_ids){

    $gene_count++;
    
    #Go through each exon cluster - for clusters of a single exon get the probe sequences, for clusters of multiple exons, define regions
    my $exons_ref = $genes_ref->{$gene_id}->{reference_exons};
    my $number_exons = keys %{$exons_ref};

    if ($verbose =~ /yes|y/i){
      print CYAN, "\nGENE: $gene_count\tALEXA: $gene_id\tENSG: $genes_ref->{$gene_id}->{ensembl_g_id} has $number_exons nr exons. Defining exon regions", RESET;
    }
    print LOG "\nGENE: $gene_count\tALEXA: $gene_id\tENSG: $genes_ref->{$gene_id}->{ensembl_g_id} has $number_exons nr exons. Defining exon regions";
    
    my $exon_clusters_ref = $genes_ref->{$gene_id}->{clusters};

    foreach my $cluster (sort {$a <=> $b} keys %{$exon_clusters_ref}){
      my $cluster_exons_aref = $exon_clusters_ref->{$cluster}->{exons};
      my $cluster_size = @{$cluster_exons_aref};

      if ($verbose =~ /yes|y/i){
        print BLUE, "\n\n\tProcessing cluster: $cluster consisting of $cluster_size exons", RESET;
      }
      print LOG "\n\n\tProcessing cluster: $cluster consisting of $cluster_size exons";
      

      my $probes_ref;

      #1.) Single-exon clusters
      if ($cluster_size == 1){
        my $exon_id =@{$cluster_exons_aref}[0];
        my $region_start = $exons_ref->{$exon_id}->{exon_start};
        my $region_end =  $exons_ref->{$exon_id}->{exon_end};
        my $base_count = ($region_end - $region_start) +1;

        if ($verbose =~ /yes|y/i){
          print BLUE, "\n\t\tSingle Exon Region: ($exons_ref->{$exon_id}->{exon_start} - $exons_ref->{$exon_id}->{exon_end}) covers $base_count bases", RESET;
        }
        print LOG "\n\t\tSingle Exon Region: ($exons_ref->{$exon_id}->{exon_start} - $exons_ref->{$exon_id}->{exon_end}) covers $base_count bases";
        
        $exon_region_id++;

        if ($gene_exon_regions{$gene_id}){

          #INSTANTIATE the Exon Region Record for SINGLE-exon clusters
          my $er_ref = $gene_exon_regions{$gene_id}{exon_regions};
          $er_ref->{$exon_region_id}->{unit1_start} = $region_start;
          $er_ref->{$exon_region_id}->{unit1_end} = $region_end;
          $er_ref->{$exon_region_id}->{base_count} = $base_count;

        }else{
          my %er;
          $er{$exon_region_id}{unit1_start} = $region_start;
          $er{$exon_region_id}{unit1_end} = $region_end;
          $er{$exon_region_id}{base_count} = $base_count;
          $gene_exon_regions{$gene_id}{exon_regions} = \%er;
        }

      }else{
        if ($verbose =~ /yes|y/i){
          print BLUE, "\n\n\tProcessing cluster: $cluster consisting of $cluster_size exons", RESET;
        }
        print LOG "\n\n\tProcessing cluster: $cluster consisting of $cluster_size exons";
        

        #2.) Multi-exon clusters
        #Get the most informative regions from the cluster of overlapping exons (non-overlaping where possible)
        my $exon_regions_ref = &defineExonRegions('-exon_object'=>$exons_ref, '-exon_ids'=>$cluster_exons_aref);

        #For each region identified attempt to extract probe sequences as done for single-exon clusters
        foreach my $region_id (sort {$a <=> $b} keys %{$exon_regions_ref}){
	  my $region_start = $exon_regions_ref->{$region_id}->{region_start};
	  my $region_end = $exon_regions_ref->{$region_id}->{region_end};
          my $base_count = $exon_regions_ref->{$region_id}->{base_count};

          if ($verbose =~ /yes|y/i){
	    print BLUE, "\n\t\tExon Region: ($exon_regions_ref->{$region_id}->{region_start} - $exon_regions_ref->{$region_id}->{region_end}) covers $base_count bases", RESET;
          }
	  print LOG "\n\t\tExon Region: ($exon_regions_ref->{$region_id}->{region_start} - $exon_regions_ref->{$region_id}->{region_end}) covers $base_count bases";
          
          #INSTANTIATE the Exon Region Record for MULTI-exon clusters
	  $exon_region_id++;
          if ($gene_exon_regions{$gene_id}){
            my $er_ref = $gene_exon_regions{$gene_id}{exon_regions};
            $er_ref->{$exon_region_id}->{unit1_start} = $region_start;
            $er_ref->{$exon_region_id}->{unit1_end} = $region_end;
            $er_ref->{$exon_region_id}->{base_count} = $base_count;

          }else{
            my %er;
            $er{$exon_region_id}{unit1_start} = $region_start;
            $er{$exon_region_id}{unit1_end} = $region_end;
            $er{$exon_region_id}{base_count} = $base_count;
            $gene_exon_regions{$gene_id}{exon_regions} = \%er;
          }
        }
      }
    }
    if ($verbose =~ /yes|y/i){
      print CYAN, "\n\nExon regions found: $exon_region_id", RESET;
      print CYAN, "\n****************************************************************************************************************************************\n", RESET;
    }
    print LOG "\n\nExon regions found: $exon_region_id";
    print LOG "\n****************************************************************************************************************************************\n";
  }

  print CYAN, "\nTotal Exon Regions Found: $exon_region_id\n", RESET;
  print LOG "\nTotal Exon Regions Found: $exon_region_id\n";

  return();
}


################################################################################################
#For exon clusters that have been derived from overlapping exons, identify informative regions #
#For each of these regions, make note of which of the exons from the cluster of overlapping    #
#exons are involved in each of the defined regions                                             #
################################################################################################
sub defineExonRegions{
  my %args = @_;

  my $exon_object_ref = $args{'-exon_object'};
  my $exon_ids_aref = $args{'-exon_ids'};

  #The exons in @exon_ids comprise a single cluster of overlapping exons
  #First try to identify regions of these exons that are unique to each exon

  #Compile a non-redundant list of start/end positions
  my %boundaries;
  my @boundaries;
  foreach my $exon_id (@{$exon_ids_aref}){
    #Use a hash to compile a non-redundant list of start/end positions
    my $exon_start = $exon_object_ref->{$exon_id}->{exon_start};
    my $exon_end = $exon_object_ref->{$exon_id}->{exon_end};
    $boundaries{$exon_start}{tmp}='na';
    $boundaries{$exon_end}{tmp}='na';
 }
  #Create a sorted array of these boundaries
  foreach my $bound (sort {$a <=> $b} keys %boundaries){
    push (@boundaries, $bound);
  }

  #Now consider the regions between the boundaries
  my $number_boundaries = @boundaries;
  my %regions;
  my $region_count;

  for (my $i = 0; $i < $number_boundaries-1; $i++){
    my $region_start = ($boundaries[$i]);
    my $region_end = ($boundaries[$i+1]);

    #Confirm that the region selected is actually valid and larger than the required probe size
    if ($region_end <= $region_start){
      print RED, "\nInvalid coordinate??", RESET;
      exit();
    }

    #Avoid any overlap between adjacent regions
    if ($i > 0){
      $region_start++;
    }
    my $base_count = ($region_end - $region_start)+1;

    $region_count++;
    $regions{$region_count}{region_start} = $region_start;
    $regions{$region_count}{region_end} = $region_end;
    $regions{$region_count}{base_count} = $base_count;
  }

  return(\%regions);
}

################################################################################################
#Set names for exon regions                                                                    #
################################################################################################
sub nameExonRegions{

  print BLUE, "\n\nNaming exon regions\n", RESET;
  print LOG "\n\nNaming exon regions\n", RESET;

  foreach my $gene_id (sort {$a <=> $b} keys %gene_exon_regions){

    my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};
    my $exon_region_ref = $gene_exon_regions{$gene_id}{exon_regions};

    my $exon_content_blocks = keys %{$exon_content_ref};
    my $exon_regions = keys %{$exon_region_ref};

    #1.) Assign 'exon order' to each exon content block
    my $exon_order = 0;
    foreach my $ec_id (sort {$exon_content_ref->{$a}->{start} <=> $exon_content_ref->{$b}->{start}} keys %{$exon_content_ref}){
      #Note that since we are dealing with GENE coordinates which are always stored 5' to 3', the strand does no matter at this stage
      $exon_order++;
      $exon_content_ref->{$ec_id}->{exon_order} = $exon_order;

      my %tmp;
      $exon_content_ref->{$ec_id}->{exon_region_matches} = \%tmp;
    }


    #2.) Now go through each Exon Region defined for this gene
    #    - Compare the coordinates of the exon region to the exon content coordinates for the gene 
    #    - Each exon region should be contained within a single exon content block

    foreach my $er_id (sort {$a <=> $b} keys %{$exon_region_ref}){
      my $ec_match_count = 0;
      my $ec_match_id;

      my $er_start = $exon_region_ref->{$er_id}->{unit1_start};
      my $er_end = $exon_region_ref->{$er_id}->{unit1_end};

      foreach my $ec_id (sort {$a <=> $b} keys %{$exon_content_ref}){

        my $ec_start = $exon_content_ref->{$ec_id}->{start};
        my $ec_end = $exon_content_ref->{$ec_id}->{end};


        if (($er_start >= $ec_start) && ($er_start <= $ec_end) && ($er_end >= $ec_start) && ($er_end <= $ec_end)){
          $ec_match_count++;
          $ec_match_id = $ec_id;
        }
      }

      if ($ec_match_count == 1 && $ec_match_id =~ /\d+/){
        my $er_match_ref = $exon_content_ref->{$ec_match_id}->{exon_region_matches};
        $er_match_ref->{$er_id}->{region_start} = $er_start;
      }else{
        print RED, "\n\nFound no exon content matches for an exon region or found more than one! (gene_id: $gene_id) ec_match_count = $ec_match_count (er_start = $er_start\ter_end = $er_end)\n", RESET;
        print Dumper $exon_content_ref;
        exit();
      }
    }

    #Now go through all of the exon content blocks (which have been mapped to all the exon regions) and assign names to each exon region
    foreach my $ec_id (sort {$exon_content_ref->{$a}->{start} <=> $exon_content_ref->{$b}->{start}} keys %{$exon_content_ref}){
      my $exon_order = $exon_content_ref->{$ec_id}->{exon_order};

      my $er_match_ref = $exon_content_ref->{$ec_id}->{exon_region_matches};
      my @label = @letters;
      foreach my $er_match_id (sort {$er_match_ref->{$a}->{region_start} <=> $er_match_ref->{$b}->{region_start}} keys %{$er_match_ref}){
        my $label = shift @label;
        my $name =  "ER"."$exon_order"."$label";
        $exon_region_ref->{$er_match_id}->{er_name} = "$name";
      }
    }
  } 
  return();
}


################################################################################################
#Determine mRNA/EST sequence support for each exon region                                      #
################################################################################################
sub determineSequenceSupport{
  my %args = @_;

  #Go through each exon region and determine if it is completely contained within an mRNA, EST, xmRNA or xEST exon
  my $mrna_supported_exon_regions = 0;
  my $est_supported_exon_regions = 0;
  my $xmrna_supported_exon_regions = 0;
  my $xest_supported_exon_regions = 0;

  my $counter = 0;

  print BLUE, "\n\nDetermining sequence support\n", RESET;
  print LOG "\n\nDetermining sequence support\n";

  foreach my $gene_id (sort {$a <=> $b} keys %gene_exon_regions){

    my $exon_region_ref = $gene_exon_regions{$gene_id}{exon_regions};

    my $chromosome = $genes_ref->{$gene_id}->{chromosome};
    my $strand = $genes_ref->{$gene_id}->{chr_strand};
    if ($strand eq "1"){
      $strand = "+";
    }else{
      $strand = "-";
    }
    my $chr_strand = "chr$chromosome"."_"."$strand";

    foreach my $er_id (sort {$a <=> $b} keys %{$exon_region_ref}){
      $counter++;
      if ($counter == 100){
        $counter = 0;
        $| = 1; print BLUE, ".", RESET; $| = 0;
      }

      #Initialize variables
      $exon_region_ref->{$er_id}->{supporting_mrna_count} = 0;
      $exon_region_ref->{$er_id}->{supporting_est_count} = 0;
      $exon_region_ref->{$er_id}->{supporting_xmrna_count} = 0;
      $exon_region_ref->{$er_id}->{supporting_xest_count} = 0;
      $exon_region_ref->{$er_id}->{conserved_species_count} = 0;

      my $er_start_chr = $exon_region_ref->{$er_id}->{unit1_start_chr};
      my $er_end_chr = $exon_region_ref->{$er_id}->{unit1_end_chr};

      #mRNA
      if ($mrna_exons->{$chr_strand}){
        my $mrna_ref = $mrna_exons->{$chr_strand};
        foreach my $mrna_exon (keys %{$mrna_ref}){

          my $string = $mrna_exon;
          if ($string =~ /(\d+)\_(\d+)/){
            my $mrna_exon_start = $1;
            my $mrna_exon_end = $2;

            if (($er_start_chr >= $mrna_exon_start) && ($er_start_chr <= $mrna_exon_end) && ($er_end_chr >= $mrna_exon_start) && ($er_end_chr <= $mrna_exon_end)){
              $exon_region_ref->{$er_id}->{supporting_mrna_count} += ($mrna_ref->{$mrna_exon}->{c});
            }
          }
        }
      }else{
        $exon_region_ref->{$er_id}->{supporting_mrna_count} = "na";
      }

      #EST
      if ($est_exons->{$chr_strand}){
        my $est_ref = $est_exons->{$chr_strand};
        foreach my $est_exon (keys %{$est_ref}){

          my $string = $est_exon;
          if ($string =~ /(\d+)\_(\d+)/){
            my $est_exon_start = $1;
            my $est_exon_end = $2;

            if (($er_start_chr >= $est_exon_start) && ($er_start_chr <= $est_exon_end) && ($er_end_chr >= $est_exon_start) && ($er_end_chr <= $est_exon_end)){
              $exon_region_ref->{$er_id}->{supporting_est_count} += ($est_ref->{$est_exon}->{c});
            }
          }
        }
      }else{
        $exon_region_ref->{$er_id}->{supporting_est_count} = "na";
      }

      #For OTHER SPECIES - allow some wiggle room in the coordinates of the boundary by making it a bit smaller
      $er_start_chr += $wiggle;
      $er_end_chr -= $wiggle;

      #Keep a list of all species IDs encountered for the current junction_id
      my %species_list;

      #xmRNA
      if ($xmrna_exons->{$chr_strand}){
        my $xmrna_ref = $xmrna_exons->{$chr_strand};
        foreach my $xmrna_exon (keys %{$xmrna_ref}){

          my $string = $xmrna_exon;
          if ($string =~ /(\d+)\_(\d+)/){
            my $xmrna_exon_start = $1;
            my $xmrna_exon_end = $2;

            if (($er_start_chr >= $xmrna_exon_start) && ($er_start_chr <= $xmrna_exon_end) && ($er_end_chr >= $xmrna_exon_start) && ($er_end_chr <= $xmrna_exon_end)){
              $exon_region_ref->{$er_id}->{supporting_xmrna_count} += ($xmrna_ref->{$xmrna_exon}->{c});
              my $s_ref = $xmrna_ref->{$xmrna_exon}->{s};
              foreach my $s (keys %{$s_ref}){
                $species_list{$s}='';
              }
            }
          }
        }
      }else{
        $exon_region_ref->{$er_id}->{supporting_xmrna_count} = "na";
      }

      #EST
      if ($xest_exons->{$chr_strand}){
        my $xest_ref = $xest_exons->{$chr_strand};
        foreach my $xest_exon (keys %{$xest_ref}){

          my $string = $xest_exon;
          if ($string =~ /(\d+)\_(\d+)/){
            my $xest_exon_start = $1;
            my $xest_exon_end = $2;

            if (($er_start_chr >= $xest_exon_start) && ($er_start_chr <= $xest_exon_end) && ($er_end_chr >= $xest_exon_start) && ($er_end_chr <= $xest_exon_end)){
              $exon_region_ref->{$er_id}->{supporting_xest_count} += ($xest_ref->{$xest_exon}->{c});
              my $s_ref = $xest_ref->{$xest_exon}->{s};
              foreach my $s (keys %{$s_ref}){
                $species_list{$s}='';
              }
            }
          }
        }
      }else{
        $exon_region_ref->{$er_id}->{supporting_xest_count} = "na";
      }

      unless ($exon_region_ref->{$er_id}->{supporting_mrna_count} =~ /^na$|^0$/ ){$mrna_supported_exon_regions++;}
      unless ($exon_region_ref->{$er_id}->{supporting_est_count} =~ /^na$|^0$/){$est_supported_exon_regions++;}
      unless ($exon_region_ref->{$er_id}->{supporting_xmrna_count} =~ /^na$|^0$/ ){$xmrna_supported_exon_regions++;}
      unless ($exon_region_ref->{$er_id}->{supporting_xest_count} =~ /^na$|^0$/){$xest_supported_exon_regions++;}

      #Count up all the species encountered in boundary matches from xmrnas and xests
      my $conserved_species_count = keys %species_list;
      $exon_region_ref->{$er_id}->{conserved_species_count} = $conserved_species_count;

    }
  }

  print BLUE, "\n\nFound mRNA support for $mrna_supported_exon_regions exon regions", RESET;
  print BLUE, "\nFound EST support for $est_supported_exon_regions exon regions", RESET;
  print BLUE, "\nFound xmRNA support for $xmrna_supported_exon_regions exon regions", RESET;
  print BLUE, "\nFound xEST support for $xest_supported_exon_regions exon regions\n\n", RESET;

  print LOG "\n\nFound mRNA support for $mrna_supported_exon_regions exon regions";
  print LOG "\nFound EST support for $est_supported_exon_regions exon regions";
  print LOG "\nFound xmRNA support for $xmrna_supported_exon_regions exon regions";
  print LOG "\nFound xEST support for $xest_supported_exon_regions exon regions\n\n";

  return();
}


################################################################################################
#Identify EnsEMBL transcript supporting each Exon Region                                       #
#- also identify exon regions specific to a single EnsEMBL transcript                          #
################################################################################################
sub determineEnsemblSupport{

  print BLUE, "\n\nDetermining EnsEMBL transcript support for each exon region", RESET;
  print LOG "\n\nDetermining EnsEMBL transcript support for each exon region";

  foreach my $gene_id (sort {$a <=> $b} keys %gene_exon_regions){
    my $exon_region_ref = $gene_exon_regions{$gene_id}{exon_regions};

    my $masked_gene_seq = $masked_gene_ref->{$gene_id}->{sequence};
    my $protein_bases_ref = &getProteinBases ('-gene_id'=>$gene_id, '-genes_ref'=>$genes_ref, '-gene_transcripts_ref'=>$gene_transcripts_ref);

    #Get the transcripts for this gene
    my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};
    $gene_transcripts_ref->{$gene_id}->{trans_count} = keys %{$transcripts_ref};

    my $counter = 0;
    foreach my $er_id (sort {$a <=> $b} keys %{$exon_region_ref}){
      $counter++;
      if ($counter == 100){
        $counter = 0;
        $| = 1; print BLUE, ".", RESET; $| = 0;
      }

      $exon_region_ref->{$er_id}->{specific_trans_id} = "na";
      $exon_region_ref->{$er_id}->{supporting_ensembl_count} = 0;
      my @temp;
      $exon_region_ref->{$er_id}->{transcripts} = \@temp;

      my $er_start = $exon_region_ref->{$er_id}->{unit1_start};
      my $er_end = $exon_region_ref->{$er_id}->{unit1_end};

      #Go through each transcript for the gene
      foreach my $trans_id (keys %{$transcripts_ref}){
        my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

        #Go through each of the exons of this transcript and see if one of them overlaps the current exon region
        my $overlap = 0;
        foreach my $trans_exon_id (keys %{$exons_ref}){
          my $e_start = $exons_ref->{$trans_exon_id}->{exon_start};
          my $e_end = $exons_ref->{$trans_exon_id}->{exon_end};

          if (($er_start >= $e_start) && ($er_start <= $e_end) && ($er_end >= $e_start) && ($er_end <= $e_end)){
            $overlap = 1;
          }
        }
        if ($overlap > 0){
          $exon_region_ref->{$er_id}->{supporting_ensembl_count}++;
          push(@{$exon_region_ref->{$er_id}->{transcripts}}, $trans_id)
        }
      }

      my $trans_id_modifier;
      if ($gene_transcripts_ref->{$gene_id}->{trans_count} == 1){
        $trans_id_modifier = "S";
      }else{
        $trans_id_modifier = "M";
      }	

      if ($exon_region_ref->{$er_id}->{supporting_ensembl_count} == 1){
        my $trans_id = @{$exon_region_ref->{$er_id}->{transcripts}}[0];
        $exon_region_ref->{$er_id}->{specific_trans_id} = "$trans_id_modifier"."_"."$trans_id";
      }

      #Get the number of masked bases and coding bases for this exon region
      #Masked Seq
      my $start = $exon_region_ref->{$er_id}->{unit1_start};
      my $size = $exon_region_ref->{$er_id}->{base_count};
      my $masked_exon_region_seq = substr ($masked_gene_seq, $start, $size);
      my $masked_n_count = ($masked_exon_region_seq =~ tr/N/N/);
      my $unmasked_base_count = $size - $masked_n_count;

      #Protein coding bases
      my $chr_start = $exon_region_ref->{$er_id}->{unit1_start_chr};
      my $chr_end = $exon_region_ref->{$er_id}->{unit1_end_chr};
      my $coding_base_count = &getProteinBaseCount ('-protein_bases_ref'=>$protein_bases_ref, '-chr_start'=>$chr_start, '-chr_end'=>$chr_end);

      $exon_region_ref->{$er_id}->{unmasked_base_count} = $unmasked_base_count;
      $exon_region_ref->{$er_id}->{coding_base_count} = $coding_base_count;

    }
  }

  return();
}


################################################################################################
#Print out the exon region annotation file                                                     #
################################################################################################
sub printExonRegionFile{
  my %args = @_;

  print BLUE, "\n\nPrinting exon region file for chr: $chr_filter:$start_filter-$end_filter\n\n", RESET;
  print LOG "\n\nPrinting exon region file for chr: $chr_filter:$start_filter-$end_filter\n\n";

  #ALEXA Gene ID
  #EnsEMBL Gene ID
  #Gene Name
  #Chromosome
  #Start coord (gene based, always 1, 5' to 3') 
  #End Coord (gene based, always equals length of gene including introns, 5' to 3')
  #Strand
  #Chromosome Start Coord (Start position on the chromosome)
  #Chromosome End Coord (End position on the chromosome)
  #Length (Number of bases in exon region)
  #Exon Name (i.e. 1, 2, ,3.   OR in the case of regions extracted from a cluster of overlapping exons, 1a, 1b, 1c, ...
  #Number supporting EnsEMBL transcripts
  #Number supporting mRNAs
  #Number supporting ESTs
  #Transcript specific ID (Transcript ID if this region corresponds to a single transcript only)
  #- 'S' denotes a region specific to a single transcript where only one transcript exists for the gene anyway
  #- 'M' denotes a region specific to a single transcript where multiple transcripts were possible for that gene 

  foreach my $gene_id (sort {$a <=> $b} keys %gene_exon_regions){
    my $exon_region_ref = $gene_exon_regions{$gene_id}{exon_regions};

    my $counter = 0;
    foreach my $er_id (sort {$a <=> $b} keys %{$exon_region_ref}){
      $counter++;
      if ($counter == 100){
        $counter = 0;
        $| = 1; print BLUE, ".", RESET; $| = 0;
      }

      my $unique_id = "$chr_filter"."_"."$region_number"."_"."$er_id";

      #Print out data line
      #"Gene_ID\tEnsEMBL_Gene_ID\tGene_Name\tChromosome\tStrand\tUnit1_start\tUnit1_end\tUnit1_start_chr\tUnit1_end_chr\tExon_Region_Length\tSupporting_EnsEMBL_Count\tSupporting_mRNA_Count\tSupporting_EST_Count\tExon_Region_Name\tSpecific_Trans_ID\n";


      #UnMasked Seq
      my $gene_seq = $genes_ref->{$gene_id}->{sequence};

      my $start = $exon_region_ref->{$er_id}->{unit1_start};
      my $size = $exon_region_ref->{$er_id}->{base_count};
      my $exon_region_seq = substr ($gene_seq, $start, $size);

      #Masked Seq
      my $masked_gene_seq = $masked_gene_ref->{$gene_id}->{sequence};

      print OUTFILE "$unique_id\t$gene_id\t$genes_ref->{$gene_id}->{ensembl_g_id}\t$genes_ref->{$gene_id}->{gene_name}\t$genes_ref->{$gene_id}->{chromosome}\t$genes_ref->{$gene_id}->{chr_strand}\t$exon_region_ref->{$er_id}->{unit1_start}\t$exon_region_ref->{$er_id}->{unit1_end}\t$exon_region_ref->{$er_id}->{unit1_start_chr}\t$exon_region_ref->{$er_id}->{unit1_end_chr}\t$exon_region_ref->{$er_id}->{base_count}\t$exon_region_ref->{$er_id}->{unmasked_base_count}\t$exon_region_ref->{$er_id}->{coding_base_count}\t$exon_region_ref->{$er_id}->{supporting_ensembl_count}\t$exon_region_ref->{$er_id}->{supporting_mrna_count}\t$exon_region_ref->{$er_id}->{supporting_est_count}\t$exon_region_ref->{$er_id}->{supporting_xmrna_count}\t$exon_region_ref->{$er_id}->{supporting_xest_count}\t$exon_region_ref->{$er_id}->{conserved_species_count}\t$exon_region_ref->{$er_id}->{er_name}\t$exon_region_ref->{$er_id}->{specific_trans_id}\n";

      if ($size >= 36){
        print FASTA ">$unique_id\n$exon_region_seq\n";
      }
    }
  }
  return();
}






