#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to analyze a file of exon boundary sequences, name them and determine basic info about them
#e.g.EnseMBL transcript support (i.e KNOWN vs NOVEL), EST and/or mRNA support, specific transcript IDs 

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
use utilities::ALEXA_DB qw(:all);
use utilities::utility qw(:all);
use utilities::ucsc qw(:all);

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $chr_filter = '';
my $target_size = '';
my $boundary_file = '';
my $ucsc_align_dir = '';
my $genbank_mapfile = '';
my $wiggle = '';
my $outfile = '';
my $logfile = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'chr_filter=s'=>\$chr_filter, 'target_size=i'=>\$target_size, 'wiggle=i'=>\$wiggle,
	    'ucsc_align_dir=s'=>\$ucsc_align_dir, 'genbank_mapfile=s'=>\$genbank_mapfile,
	    'boundary_file=s'=>\$boundary_file, 'outfile=s'=>\$outfile, 'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script annotates exon boundary sequences with basic info", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password\n", RESET;
print GREEN, "\n\tSpecify a genomic region to analyze using: --chr_filter", RESET;
print GREEN, "\n\tSpecify a tab-delimited input exon boundary file using: --boundary_file", RESET;
print GREEN, "\n\t\tIt is assumed the the 1st column contains a unique ID", RESET;
print GREEN, "\n\tSpecify the target size used to generate the database of hypothetical exon boundary extensions using:  --target_size", RESET;
print GREEN, "\n\tSpecify a directory containing UCSC mRNA/EST/xmRNA and xEST files using: --ucsc_align_dir", RESET;
print GREEN, "\n\tSpecify a berkeley DB containing GenBank-To-Species mappings using:  --genbank_mapfile", RESET;
print GREEN, "\n\tJunctions for xmRNA and xEST will be allowed to 'wiggle' by the amount specied by:  --wiggle", RESET;
print GREEN, "\n\tSpecify the name of the resulting outfile using: --outfile", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of this script: --logfile", RESET;
print GREEN, "\n\nExample: annotateExonBoundaries.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --chr_filter='3:16:121020102-126260001'  --target_size=62  --boundary_file=/projects/malachig/sequence_databases/hs_53_36o/ensembl_exonBoundaries_hs_53_36o/exonBoundaries_62mers.txt  --ucsc_align_dir=/projects/malachig/sequence_databases/hs_53_36o/mrna_est/partitions/  --genbank_mapfile=/projects/malachig/sequence_databases/hs_53_36o/mrna_est/GenBankToOrganism.btree  --wiggle=3  --outfile=/projects/malachig/sequence_databases/hs_53_36o/ensembl_exonBoundaries_hs_53_36o/temp/exonBoundaries_annotated_3_16.txt  --logfile=/projects/malachig/sequence_databases/hs_53_36o/logs/annotateExonBoundaries/annotateExonBoundaries_3_16_LOG.txt\n\n", RESET;

unless ($database && $server && $user && $password && $chr_filter && $target_size && $boundary_file && $ucsc_align_dir && $genbank_mapfile && ($wiggle =~ /\d+/) && $outfile && $logfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

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

#Get the UCSC alignment files (one each for target species mRNA and EST and one each for all other species mRNA and EST)
my $tmp_chr = "chr$chr_filter";
if ($tmp_chr eq "chrMT"){$tmp_chr = "chrM";}

my $ucsc_mrna_table = "$ucsc_align_dir"."$tmp_chr"."_mrna.txt.gz";
my $ucsc_est_table = "$ucsc_align_dir"."$tmp_chr"."_est.txt.gz";
my $ucsc_xmrna_table = "$ucsc_align_dir"."$tmp_chr"."_xmrna.txt.gz";
my $ucsc_xest_table = "$ucsc_align_dir"."$tmp_chr"."_xest.txt.gz";

#Load the GenBank-To-Species mapfile
my %gb_org_map;
tie(%gb_org_map, 'BerkeleyDB::Btree', -Cachesize =>256000000, -Filename=> $genbank_mapfile , -Flags => DB_RDONLY) or die "can't open file $genbank_mapfile: $! $BerkeleyDB::Error\n";

#Global variable storing letters for labeling of Acceptors and Donors
#my @letters = ('aa' .. 'zz');
my @letters = ('a' .. 'z');

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#Print out the parameters supplied by the user to the logfile for future reference
print LOG "\nUser Specified the following options:\ndatabase = $database\nboundary_file = $boundary_file\noutfile = $outfile\nlogfile = $logfile\n\n";

#1.) First get all the boundary IDs from the input file, as well as their coordinates
my %boundaries;
my %master_gene_list;

print BLUE, "\n1.) Parsing exon boundary file: $boundary_file for coordinate data ($chr_filter:$region_number:$start_filter-$end_filter)\n\n", RESET;
print LOG "\n1.) Parsing exon boundary file: $boundary_file for coordinate data ($chr_filter:$region_number:$start_filter-$end_filter)\n\n";


open (BOUNDARIES, "$boundary_file") || die "\nCould not open input boundary file: $boundary_file\n\n";

my $header = 1;
my %columns;
while (<BOUNDARIES>){
  chomp($_);
  my @line = split("\t", $_);

  if ($header == 1){
    my $column_count = 0;

    foreach my $column (@line){
      $columns{$column}{column_pos} = $column_count;
      $column_count++;
    }

    #Check for critical columns and their names
    unless ($columns{'Gene_ID'} && $columns{'Unit1_start'} && $columns{'Unit1_end'} && $columns{'Unit1_start_chr'} && $columns{'Unit1_end_chr'}){
      print RED, "\nCritical column missing or named incorrectly in boundary file: $boundary_file, check input file\n\n", RESET;
      close (LOG);
      exit();
    }
    $header = 0;
    next();
  }

  #Skip this boundary unless it is within the target region
  my $chr = $line[$columns{'Chromosome'}{column_pos}];
  my @tmp = ($line[$columns{'Unit1_start_chr'}{column_pos}], $line[$columns{'Unit1_end_chr'}{column_pos}]);
  my @tmp_sort = sort {$a <=> $b} @tmp;
  my $lower = $tmp_sort[0];
  my $upper = $tmp_sort[1];
  unless ($chr eq $chr_filter && ($lower >= $start_filter && $lower <= $end_filter) && ($upper >= $start_filter && $upper <= $end_filter)){
    next();
  }

  my $boundary_id = $line[0];
  $master_gene_list{$line[$columns{'Gene_ID'}{column_pos}]}{tmp} = '';

  $boundaries{$boundary_id}{gene_id} = $line[$columns{'Gene_ID'}{column_pos}];
  $boundaries{$boundary_id}{chromosome} = $line[$columns{'Chromosome'}{column_pos}];
  $boundaries{$boundary_id}{strand} = $line[$columns{'Strand'}{column_pos}];
  $boundaries{$boundary_id}{type} = $line[$columns{'Boundary_Type'}{column_pos}];

  $boundaries{$boundary_id}{start} = $line[$columns{'Unit1_start'}{column_pos}];
  $boundaries{$boundary_id}{end} = $line[$columns{'Unit1_end'}{column_pos}];

  $boundaries{$boundary_id}{start_chr} = $line[$columns{'Unit1_start_chr'}{column_pos}];
  $boundaries{$boundary_id}{end_chr} = $line[$columns{'Unit1_end_chr'}{column_pos}];

  #Make sure the boundary sequence length matches the the target length supplied by the user
  my $size = $line[$columns{'Base_Count'}{column_pos}];
  unless ($size == $target_size){
    print RED, "Length of boundary sequence ($size does not match the target size specified by the user ($target_size)!\n\n", RESET;
    exit();
  }

  #Format corrections to deal with UCSC
  if ($boundaries{$boundary_id}{strand} eq "-1"){
    $boundaries{$boundary_id}{strand} = "-";
  }else{
    $boundaries{$boundary_id}{strand} = "+";
  }
  $boundaries{$boundary_id}{chromosome} = "chr"."$boundaries{$boundary_id}{chromosome}";

}
close (BOUNDARIES);


#2.) Get gene info from ALEXA database

#Only use those genes actually found in the input boundary file
my @gene_ids = keys %master_gene_list;

my $gene_count = @gene_ids;
my $total_boundary_count = keys %boundaries;

print BLUE, "\nFound $total_boundary_count exon boundaries for a total of $gene_count genes\n\n", RESET;
print LOG "\nFound $total_boundary_count exon boundaries for a total of $gene_count genes\n\n";

#Now for each of these genes, get the transcripts and exons
print BLUE, "\n2.) Extracting Gene and Exon data from ALEXA\n\n", RESET;
print LOG "\n2.) Extracting Gene and Exon data from ALEXA\n\n";

my $genes_ref;
my $gene_transcripts_ref;
my $gene_exon_content_ref;
if ($gene_count > 0){
  my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
  $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");
  $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");
  $gene_exon_content_ref = &getExonContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);
  $alexa_dbh->disconnect();
}else{
  #Some regions may not have any multi-exon genes and therefore no boundaries
  #If there are no genes found, it is pointless to continue.
  print "\nSince no boundaries were found in this region ...\nSCRIPT COMPLETE\n";
  print LOG "\nSince no boundaries were found in this region ...\nSCRIPT COMPLETE\n";
  close(LOG);
  exit();
}


print BLUE, "\n\n2-A.) Generating reference exon objects\n\n", RESET;
print LOG "\n\n2-A.) Generating reference exon objects\n\n";

my $ref_exon_counter = 0;
foreach my $gene_id (@gene_ids){

  my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};
  my $trans_count = keys %{$transcripts_ref};

  $gene_transcripts_ref->{$gene_id}->{trans_count} = $trans_count;

  #2-A.) Assemble a reference set of exons (a superset of all non-redundant exons)
  my %reference_exons;

  foreach my $trans_id (sort {$a <=> $b} keys %{$transcripts_ref}){
    my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

    foreach my $trans_exon_id (sort {$exons_ref->{$a}->{exon_start} <=> $exons_ref->{$b}->{exon_start}} keys %{$exons_ref}){
      my $trans_exon_start = $exons_ref->{$trans_exon_id}->{exon_start};
      my $trans_exon_end = $exons_ref->{$trans_exon_id}->{exon_end};

      #Check each of the reference exons to see if one of these is the same, otherwise add it to the list
      my $redundant_exon = 0;
      my $matching_exon_count;
      foreach my $ref_exon (sort {$reference_exons{$a}->{exon_start} <=> $reference_exons{$b}->{exon_start}} keys %reference_exons){
	if ($trans_exon_start == $reference_exons{$ref_exon}{exon_start} && $trans_exon_end == $reference_exons{$ref_exon}{exon_end}){
	  $redundant_exon = 1;
          $matching_exon_count = $ref_exon;
	}
        #print "DEBUG: G: $gene_id\tT: $trans_id\tTE: $trans_exon_id\tRE: $ref_exon\t$trans_exon_start == $reference_exons{$ref_exon}{exon_start} && $trans_exon_end == $reference_exons{$ref_exon}{exon_end}) $redundant_exon\n";
      }
      #Unless the current transcript exon was found to be redundant, add it to the list
      if ($redundant_exon == 0){
        $ref_exon_counter++;
	$reference_exons{$ref_exon_counter}{exon_start} = $trans_exon_start;
	$reference_exons{$ref_exon_counter}{exon_end} = $trans_exon_end;
        my @transcripts;
        push(@transcripts, $trans_id);
        $reference_exons{$ref_exon_counter}{transcripts} = \@transcripts;
      }else{
        push(@{$reference_exons{$matching_exon_count}{transcripts}}, $trans_id);
      }
    }
  }

  #2-B.) Get arrays to represent the reference exons
  my $ref_exon_count = keys %reference_exons;

  my @reference_exon_starts;
  my @reference_exon_ends;
  foreach my $ref_exon_id (sort {$reference_exons{$a}->{exon_start} <=> $reference_exons{$b}->{exon_start}} keys %reference_exons){
    push (@reference_exon_starts, $reference_exons{$ref_exon_id}{exon_start});
    push (@reference_exon_ends, $reference_exons{$ref_exon_id}{exon_end});
  }

  unless($ref_exon_count == 1){
    $gene_transcripts_ref->{$gene_id}->{ref_exon_starts} = \@reference_exon_starts;
    $gene_transcripts_ref->{$gene_id}->{ref_exon_ends} = \@reference_exon_ends;
    $gene_transcripts_ref->{$gene_id}->{ref_exon_count} = $ref_exon_count;
    $gene_transcripts_ref->{$gene_id}->{ref_exons} = \%reference_exons;
  }


  my %starts;    #Non-redundant list of known start positions for all transcripts of the current gene
  my %ends;      #Non-redundant list of known end positions for all transcripts of the current gene
  foreach my $trans_id (keys %{$transcripts_ref}){
    my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

    #Build arrays of the exon starts and ends for this transcript
    foreach my $trans_exon_id (sort {$exons_ref->{$a}->{exon_start} <=> $exons_ref->{$b}->{exon_start}} keys %{$exons_ref}){
      my $trans_exon_start = $exons_ref->{$trans_exon_id}->{exon_start};
      my $trans_exon_end = $exons_ref->{$trans_exon_id}->{exon_end};
      $starts{$trans_exon_start}{label} = "na";
      $ends{$trans_exon_end}{label} = "na";
    }
  }

  #2-C.) Now assign 'exon order' to each exon content block and label the acceptors and donors within each exon content block as 'a .. z'
  my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};
  my $exon_count = keys %{$exon_content_ref};

  my $exon_order = 0;

  foreach my $exon_id (sort {$exon_content_ref->{$a}->{start} <=> $exon_content_ref->{$b}->{start}} keys %{$exon_content_ref}){
    my $ec_start = $exon_content_ref->{$exon_id}->{start};
    my $ec_end = $exon_content_ref->{$exon_id}->{end};

    #Note that since we are dealing with GENE coordinates which are always stored 5' to 3', the strand does no matter at this stage
    $exon_order++;
    $exon_content_ref->{$exon_id}->{exon_order} = $exon_order;

    #For this exon content block, identify the exon starts and ends that are within it:
    #Then label the donors (exon ends) and acceptor (exon starts) from A to Z (5' to 3')

    my %ec_starts; #Starts that are actually within the current exon content block
    my %ec_ends;   #Ends that are actually within the current exon content block

    foreach my $start (keys %starts){
      if ($start >= $ec_start && $start <= $ec_end){
	$ec_starts{$start}{tmp} = "na";
      }
    }
    foreach my $end (keys %ends){
      if ($end >= $ec_start && $end <= $ec_end){
	$ec_ends{$end}{tmp} = "na";
      }
    }

    #Now apply the labels to each Acceptor or Donor
    my @labels = @letters;
    foreach my $ec_start (sort {$a <=> $b} keys %ec_starts){
      $starts{$ec_start}{label} = shift (@labels);
    }
    @labels = @letters;
    foreach my $ec_end (sort {$a <=> $b} keys %ec_ends){
      $ends{$ec_end}{label} = shift (@labels);
    }
  }

  $gene_exon_content_ref->{$gene_id}->{acceptors} = \%starts;
  $gene_exon_content_ref->{$gene_id}->{donors} = \%ends;
}

print BLUE, "\n2-B.) Getting EXON CONTENT of each gene\n\n", RESET;
print LOG "\n2-B.) Getting EXON CONTENT of each gene\n\n";

my $counter = 0;
foreach my $gene_id (@gene_ids){

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

  $genes_ref->{$gene_id}->{exon_content_size} = 0; #Total number of bases covered by exons of this gene

  foreach my $exon_id (keys %{$exon_content_ref}){

    my $start = $exon_content_ref->{$exon_id}->{start};
    my $end = $exon_content_ref->{$exon_id}->{end};

    my $coords_ref = &convertGeneCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$start, '-end_pos'=>$end, '-ordered'=>"yes");

    $exon_content_ref->{$exon_id}->{chr_start} = $coords_ref->{$gene_id}->{chr_start};
    $exon_content_ref->{$exon_id}->{chr_end} = $coords_ref->{$gene_id}->{chr_end};
    $exon_content_ref->{$exon_id}->{size} = ($coords_ref->{$gene_id}->{chr_end} - $coords_ref->{$gene_id}->{chr_start})+1;
  }
}


#3.) Import all mRNA and EST observed exons
#    - NOTE that these databases are massive
#    - Organize into hashes according to chromosome and strand
#    - Use berkley DB to store these data structures temporarily and clean-up at the end of the script
#    - If a boundary sequence is completely contained within an EST/mRNA exon according to their coordinates this will be considered evidence that this exon extension has been observed in the data
#    - This is an overly simplistic view and this evidence should not be considered conclusive
my $mrna_exons = &importUcscAlignments('-import_type'=>"exon", '-align_file'=>$ucsc_mrna_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"mrna", '-filter'=>"0", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);
my $est_exons = &importUcscAlignments('-import_type'=>"exon", '-align_file'=>$ucsc_est_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"est", '-filter'=>"1", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);
my $xmrna_exons = &importUcscAlignments('-import_type'=>"exon", '-align_file'=>$ucsc_xmrna_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"xmrna", '-filter'=>"2", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);
my $xest_exons = &importUcscAlignments('-import_type'=>"exon", '-align_file'=>$ucsc_xest_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"xest", '-filter'=>"2", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);



#4.) Go through each exon boundary sequence and:
#    - Determine whether it corresponds to a known EnsEMBL transcript
#    - Determine whether is corresponds to a known EST or mRNA mapped to the same strand by UCSC
my $total_known_boundary_count = 0;
my $unknown_chromosome_records = 0;
my $mrna_supported_count = 0;
my $est_supported_count = 0;
my $xmrna_supported_count = 0;
my $xest_supported_count = 0;

print BLUE, "\n\nDetermining the EnsEMBL, mRNA and EST support for each boundary", RESET;
print LOG "\n\nDetermining the EnsEMBL, mRNA and EST support for each boundary";
print BLUE, "\n\nAt this time also annotate each boundary according to the exons and acceptor/donor they involve (e.g. E2Aa, E2Ab, E2Da, etc)\n\n", RESET;
print LOG "\n\nAt this time also annotate each boundary according to the exons and acceptor/donor they involve (e.g. E2Aa, E2Ab, E2Da, etc)\n\n";

my %failed_chrs;

$counter = 0;
foreach my $boundary_id (sort {$a <=> $b} keys %boundaries){

  $counter++;
  if ($counter == 100){
    $| = 1; print BLUE, ".", RESET; $| = 0;
    $counter = 0;
  }

  my $gene_id = $boundaries{$boundary_id}{gene_id};
  my $boundary_type = $boundaries{$boundary_id}{type};
  my $boundary_start = $boundaries{$boundary_id}{start};
  my $boundary_end = $boundaries{$boundary_id}{end};

  my $boundary_start_chr = $boundaries{$boundary_id}{start_chr};
  my $boundary_end_chr = $boundaries{$boundary_id}{end_chr};

  #Sanity check 
  unless ($boundary_end > $boundary_start){
    print RED, "\nBoundary Start and End positions do not make sense!\n\n", RESET;
    close (LOG);
    exit();
  }
  
  #Check if the exon boundary coordinates of this boundary correspond to a known EnsEMBL transcript
  my $ref_exons_ref = $gene_transcripts_ref->{$gene_id}->{ref_exons};
  $boundaries{$boundary_id}{specific_trans_id} = "na";

  #Loop through the exons of the gene for this boundary and see if any of them completely cover this boundary sequence
  $boundaries{$boundary_id}{known_boundary} = 0;
  $boundaries{$boundary_id}{trans_count} = 0;
  $boundaries{$boundary_id}{mrna_supported_boundary} = 0;
  $boundaries{$boundary_id}{est_supported_boundary} = 0;
  $boundaries{$boundary_id}{xmrna_supported_boundary} = 0;
  $boundaries{$boundary_id}{xest_supported_boundary} = 0;
  $boundaries{$boundary_id}{conserved_species_count} = 0;

  foreach my $ref_exon_id (keys %{$ref_exons_ref}){
    my $ref_exon_start = $ref_exons_ref->{$ref_exon_id}->{exon_start};
    my $ref_exon_end = $ref_exons_ref->{$ref_exon_id}->{exon_end};

    if (($boundary_start >= $ref_exon_start) && ($boundary_start <= $ref_exon_end) && ($boundary_end >= $ref_exon_start) && ($boundary_end <= $ref_exon_end)){
      $boundaries{$boundary_id}{known_boundary}++;
      my $trans_count = @{$ref_exons_ref->{$ref_exon_id}->{transcripts}};
      $boundaries{$boundary_id}{trans_count} += $trans_count;
      $boundaries{$boundary_id}{last_trans_id} = @{$ref_exons_ref->{$ref_exon_id}->{transcripts}}[0];
    }
  }

  if ($boundaries{$boundary_id}{known_boundary} > 0){
    $total_known_boundary_count++;

    #Also note whether it is specific to a single transcript 
    # - Always the case when there is only one transcript, but rarely the case when multiple transcripts exist
    # - Record the transcript ID for these cases, 
    # - Denote each transcript ID as 'S' or 'M' according to whether they correspond to a Single or Multiple transcript gene

    my $trans_id_modifier;
    if ($gene_transcripts_ref->{$gene_id}->{trans_count} == 1){
      $trans_id_modifier = "S";
    }else{
      $trans_id_modifier = "M";
    }	    

    if ($boundaries{$boundary_id}{trans_count} == 1){
      my $trans_id = $boundaries{$boundary_id}{last_trans_id};
      $boundaries{$boundary_id}{specific_trans_id} = "$trans_id_modifier"."_"."$trans_id";
    }
  }


  #Test the current boundary against exons from mRNAs/ESTs/xmRNAs/xESTs
  my $current_chr_strand = "$boundaries{$boundary_id}{chromosome}"."_"."$boundaries{$boundary_id}{strand}";
  &testBoundarySupport('-boundary_id'=>$boundary_id, '-boundary_start'=>$boundary_start_chr, '-boundary_end'=>$boundary_end_chr, '-chr_strand'=>$current_chr_strand, '-wiggle'=>$wiggle);


  #Annotate each exon-exon junction to identify which exons are involved
  #Use the 'exon_order' values assigned to each exon content block
  #Label will be according to the block of overlapping exons (if applicable)
  #In the case of overlapping exons, alternate acceptor, or donor positions will be labeled

  my $strand = $genes_ref->{$gene_id}->{chr_strand};
  my $exon_order_value = "na";

  my $exon_found = 0;

  my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};

  foreach my $exon_id (sort {$exon_content_ref->{$a}->{start} <=> $exon_content_ref->{$b}->{start}} keys %{$exon_content_ref}){

    my $ref_exon_start = $exon_content_ref->{$exon_id}->{start};
    my $ref_exon_end = $exon_content_ref->{$exon_id}->{end};
    my $ref_exon_order = $exon_content_ref->{$exon_id}->{exon_order};

    #Is the current boundary (start OR end depending on whether it is a IE or EI boundary) contained within this exon content block?
    if ($boundary_type eq "EI" && $boundary_start >= $ref_exon_start-1 && $boundary_start <= $ref_exon_end){
      $exon_found++;
      $exon_order_value = $ref_exon_order;
    }
    if ($boundary_type eq "IE" && $boundary_end >= $ref_exon_start && $boundary_end <= $ref_exon_end+1){
      $exon_found++;
      $exon_order_value = $ref_exon_order;
    }
  }

  if ($exon_found > 1){
    #Since we are comparing an exon boundary coordinate to exon-content blocks, we should always find only a single target exon
    print RED,  "\n\tAmbigous boundary to exon mapping!  Boundary: $boundary_id\tExons Found = $exon_found", RESET;
    exit();
  }elsif ($exon_found == 0){
    #Every boundary coordinate should correspond to at least one exon content block
    print RED,  "\n\tUnassigned boundary to exon mapping!  Boundary: $boundary_id ($boundary_start - $boundary_end) $boundary_type", RESET;
    exit();
  }else{
    #print YELLOW, "\n\tBoundary: $boundary_id\tStrand: $strand\tWill be annotated as: $exon_order_value", RESET;
  }

  #Assign the exon_names value for this junction
  if ($boundary_type eq "EI"){
    #Donor site of Exon-Intron boundaries
    my $donors_ref = $gene_exon_content_ref->{$gene_id}->{donors};

    #Calculate the exon-END position 
    my $exon_end = $boundary_start + ($target_size/2);

    #print "\n\nExon End: $exon_end\n";
    #print Dumper $donors_ref;

    my $donor_label = $donors_ref->{$exon_end}->{label};
    unless(defined($exon_order_value)){$exon_order_value=0;}
    unless(defined($donor_label)){$donor_label="x";}
    $boundaries{$boundary_id}{name} = "E"."$exon_order_value"."_"."D"."$donor_label";
  }
  if ($boundary_type eq "IE"){
    #Acceptor site or Intron-Exon boundaries
    my $acceptors_ref = $gene_exon_content_ref->{$gene_id}->{acceptors};

    #Calculate the exon-START position
    my $exon_start = ($boundary_start + ($target_size/2))+1;

    #print "\n\nExon start: $exon_start\n";
    #print Dumper $acceptors_ref;

    my $acceptor_label = $acceptors_ref->{$exon_start}->{label};
    unless(defined($exon_order_value)){$exon_order_value=0;}
    unless(defined($acceptor_label)){$acceptor_label="x";}
    $boundaries{$boundary_id}{name} = "E"."$exon_order_value"."_"."A"."$acceptor_label";
  }

}

print BLUE, "\n\n$total_known_boundary_count boundaries of $total_boundary_count correspond to known EnsEMBL transcripts\n\n", RESET;
print LOG "\n\n$total_known_boundary_count boundaries of $total_boundary_count correspond to known EnsEMBL transcripts\n\n";
print BLUE, "\n\n$mrna_supported_count boundaries of $total_boundary_count are supported by one or more mRNAs\n\n", RESET;
print LOG "\n\n$mrna_supported_count boundaries of $total_boundary_count are supported by one or more mRNAs\n\n", RESET;
print BLUE, "\n\n$est_supported_count boundaries of $total_boundary_count are supported by one or more ESTs\n\n", RESET;
print LOG "\n\n$est_supported_count boundaries of $total_boundary_count are supported by one or more ESTs\n\n", RESET;
print BLUE, "\n\n$xmrna_supported_count boundaries of $total_boundary_count are supported by one or more xmRNAs\n\n", RESET;
print LOG "\n\n$xmrna_supported_count boundaries of $total_boundary_count are supported by one or more xmRNAs\n\n", RESET;
print BLUE, "\n\n$xest_supported_count boundaries of $total_boundary_count are supported by one or more xESTs\n\n", RESET;
print LOG "\n\n$xest_supported_count boundaries of $total_boundary_count are supported by one or more xESTs\n\n", RESET;
print BLUE, "\n\n$unknown_chromosome_records boundaries of $total_boundary_count had an unrecognized chromosome name\n\n", RESET;
print LOG "\n\n$unknown_chromosome_records boundaries of $total_boundary_count had an unrecognized chromosome name\n\n";


#6.) Finally print out the results
print BLUE, "\nPrinting output to file: $outfile\n\n", RESET;
print LOG "\nPrinting output to file: $outfile\n\n";

open (BOUNDARIES, "$boundary_file") || die "\nCould not open input boundary file: $boundary_file\n\n";
open (OUTFILE, ">$outfile") || die "\nCould not open output boundary file: $outfile\n\n";

$header = 1;
while (<BOUNDARIES>){
  chomp($_);
  if ($header == 1){
    $header = 0;
    print OUTFILE "$_\tSupporting_EnsEMBL_Count\tSupporting_mRNA_Count\tSupporting_EST_Count\tSupporting_xmRNA_Count\tSupporting_xEST_Count\tConserved_Species_Count\tSeq_Name\tSpecific_Trans_ID\n";
    next();
  }
  my @line = split("\t", $_);
  my $boundary_id = $line[0];

  unless ($boundaries{$boundary_id}){
    next();
  }

  chomp ($_);
  print OUTFILE "$_\t$boundaries{$boundary_id}{known_boundary}\t$boundaries{$boundary_id}{mrna_supported_boundary}\t$boundaries{$boundary_id}{est_supported_boundary}\t$boundaries{$boundary_id}{xmrna_supported_boundary}\t$boundaries{$boundary_id}{xest_supported_boundary}\t$boundaries{$boundary_id}{conserved_species_count}\t$boundaries{$boundary_id}{name}\t$boundaries{$boundary_id}{specific_trans_id}\n";
}

print "SCRIPT COMPLETE";
print LOG "SCRIPT COMPLETE";

close (BOUNDARIES);
close (OUTFILE);
close (LOG);

untie(%gb_org_map);

exit();


##############################################################################################################################3
#Test support for a boundary against mRNA/EST/xmRNA/xEST
##############################################################################################################################3
sub testBoundarySupport{
  my %args = @_;
  my $boundary_id = $args{'-boundary_id'};
  my $current_chr_strand = $args{'-chr_strand'};
  my $boundary_start_chr = $args{'-boundary_start'};
  my $boundary_end_chr = $args{'-boundary_end'};
  my $wiggle = $args{'-wiggle'};

  #Check if the exon boundary coordinates correspond to a known mRNA sequence
  if ($mrna_exons->{$current_chr_strand}){
    my $mrna_exons_ref = $mrna_exons->{$current_chr_strand};
    my $mrna_overlap_found = 0;

    foreach my $mrna_exon (keys %{$mrna_exons_ref}){
      my $mrna_exon_start;
      my $mrna_exon_end;
      if ($mrna_exon =~ /(\d+)\_(\d+)/){
        $mrna_exon_start = $1;
        $mrna_exon_end = $2;
      }

      if (($boundary_start_chr >= $mrna_exon_start) && ($boundary_start_chr <= $mrna_exon_end) && ($boundary_end_chr >= $mrna_exon_start) && ($boundary_end_chr <= $mrna_exon_end)){
        $boundaries{$boundary_id}{mrna_supported_boundary}+=$mrna_exons_ref->{$mrna_exon}->{c};
        $mrna_overlap_found = 1;
      }
    }
    if ($mrna_overlap_found == 1){
      $mrna_supported_count++;
    }
  }else{
    $unknown_chromosome_records++;
    $boundaries{$boundary_id}{mrna_supported_boundary} = "na";

    if ($failed_chrs{$boundaries{$boundary_id}{chromosome}}){
      $failed_chrs{$boundaries{$boundary_id}{chromosome}}{count}++;
    }else{
      $failed_chrs{$boundaries{$boundary_id}{chromosome}}{count} = 1;
    }
  }

  #Check if the exon boundary coordinates correspond to a known EST sequence
  if ($est_exons->{$current_chr_strand}){
    my $est_exons_ref = $est_exons->{$current_chr_strand};
    my $est_overlap_found = 0;

    foreach my $est_exon (keys %{$est_exons_ref}){
      my $est_exon_start;
      my $est_exon_end;
      if ($est_exon =~ /(\d+)\_(\d+)/){
        $est_exon_start = $1;
        $est_exon_end = $2;
      }

      if (($boundary_start_chr >= $est_exon_start) && ($boundary_start_chr <= $est_exon_end) && ($boundary_end_chr >= $est_exon_start) && ($boundary_end_chr <= $est_exon_end)){
        $boundaries{$boundary_id}{est_supported_boundary}+=$est_exons_ref->{$est_exon}->{c};
        $est_overlap_found = 1;
      }
    }
    if ($est_overlap_found == 1){
      $est_supported_count++;
    }
  }else{
    $unknown_chromosome_records++;
    $boundaries{$boundary_id}{est_supported_boundary} = "na";
  }

  #For OTHER SPECIES - allow some wiggle room in the coordinates of the boundary by making it a bit smaller
  $boundary_start_chr += $wiggle;
  $boundary_end_chr -= $wiggle;

  #Keep a list of all species IDs encountered for the current junction_id
  my %species_list;

  #Check if the exon boundary coordinates correspond to a known xmRNA sequence
  if ($xmrna_exons->{$current_chr_strand}){
    my $xmrna_exons_ref = $xmrna_exons->{$current_chr_strand};
    my $xmrna_overlap_found = 0;

    foreach my $xmrna_exon (keys %{$xmrna_exons_ref}){
      my $xmrna_exon_start;
      my $xmrna_exon_end;
      if ($xmrna_exon =~ /(\d+)\_(\d+)/){
        $xmrna_exon_start = $1;
        $xmrna_exon_end = $2;
      }

      if (($boundary_start_chr >= $xmrna_exon_start) && ($boundary_start_chr <= $xmrna_exon_end) && ($boundary_end_chr >= $xmrna_exon_start) && ($boundary_end_chr <= $xmrna_exon_end)){
        $boundaries{$boundary_id}{xmrna_supported_boundary}+=$xmrna_exons_ref->{$xmrna_exon}->{c};
        $xmrna_overlap_found = 1;
        my $s_ref = $xmrna_exons_ref->{$xmrna_exon}->{s};
        foreach my $s (keys %{$s_ref}){
          $species_list{$s}='';
        }
      }
    }
    if ($xmrna_overlap_found == 1){
      $xmrna_supported_count++;
    }
  }else{
    $unknown_chromosome_records++;
    $boundaries{$boundary_id}{xmrna_supported_boundary} = "na";
  }

  #Check if the exon boundary coordinates correspond to a known xEST sequence
  if ($xest_exons->{$current_chr_strand}){
    my $xest_exons_ref = $xest_exons->{$current_chr_strand};
    my $xest_overlap_found = 0;

    foreach my $xest_exon (keys %{$xest_exons_ref}){
      my $xest_exon_start;
      my $xest_exon_end;
      if ($xest_exon =~ /(\d+)\_(\d+)/){
        $xest_exon_start = $1;
        $xest_exon_end = $2;
      }

      if (($boundary_start_chr >= $xest_exon_start) && ($boundary_start_chr <= $xest_exon_end) && ($boundary_end_chr >= $xest_exon_start) && ($boundary_end_chr <= $xest_exon_end)){
        $boundaries{$boundary_id}{xest_supported_boundary}+=$xest_exons_ref->{$xest_exon}->{c};
        $xest_overlap_found = 1;
        my $s_ref = $xest_exons_ref->{$xest_exon}->{s};
        foreach my $s (keys %{$s_ref}){
          $species_list{$s}='';
        }
      }
    }
    if ($xest_overlap_found == 1){
      $xest_supported_count++;
    }
  }else{
    $unknown_chromosome_records++;
    $boundaries{$boundary_id}{xest_supported_boundary} = "na";
  }

  #Count up all the species encountered in boundary matches from xmrnas and xests
  my $conserved_species_count = keys %species_list;
  $boundaries{$boundary_id}{conserved_species_count} = $conserved_species_count;

  return();
}



