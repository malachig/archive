#!/usr/bin/perl -w

#Written by Malachi Griffith
# The purpose of this script is to identify exons and/or conserved regions within the Target Genomic regions
#  - In this case there is only one large genomic region being targeted 
#  - Focus on the exons to start with and see how it goes
#  - Compile a list of non-redundant non-overlapping exon regions (exon content)
#  - Go through all exons regions identified and ensure the minimum size specified by the user is observed.  If not, increase the size
#  - Write the resulting list out to file
#  - See below for details

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

#ALEXA libraries
#When a script is initiated, use the full path of the script location at execution to add the perl module libraries to @INC
#This should allow this scripts to work regardless of the current working directory or the script location (where it was unpacked).
#The /utilities directory must remain in the same directory as this script but the entire code directory can be moved around
BEGIN {
  my $script_dir = &File::Basename::dirname($0);
  push (@INC, $script_dir);
}
use utilities::ALEXA_DB qw(:all);
use utilities::utility qw(:all);
use utilities::Descriptive;

#Initialize command line options
my $region_file = '';
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $min_exon_region_size = '';
my $output_file = '';
my $log_file = '';

GetOptions ('region_file=s'=>\$region_file, 'database=s'=>\$database, 'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'min_exon_region_size=i'=>\$min_exon_region_size, 'output_file=s'=>\$output_file, 'log_file=s'=>\$log_file);


print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the input file containing target genomic regions using: --region_file", RESET;
print GREEN, "\n\tSpecify the source alexa database and server for exons using: --database and --server", RESET;
print GREEN, "\n\tSpecify the corresponding user and password using: --user and --password", RESET;
print GREEN, "\n\tSpecify the minimum size of an exon region allowed using: --min_exon_region_size", RESET;
print GREEN, "\n\t\tIf an exon is smaller than the minimum size it will be expanded to this size", RESET;
print GREEN, "\n\tSpecify the output file using: --output_file", RESET;
print GREEN, "\n\nExample: getExonRegionsFromTargetRegions.pl  --region_file=TargetRegions.txt  --database=ALEXA_hs_48_36j  --server=jango.bcgsc.ca  --user=malachig  --password=pwd  --min_exon_region_size=250  --output_file=TargetExonRegions.txt  --log_file=getExonRegionsFromTargetRegions_LOG.txt", RESET;


unless ($region_file && $database && $server && $user && $password && $min_exon_region_size && $output_file && $log_file){
  print RED, "\n\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

open (LOG, ">$log_file") || die "\nCould not open log file: $log_file\n\n", RESET;

#Get target regions from the region file
my %regions;
my %sequence_names;
my $input_header;
&importTargetRegions('-region_file'=>$region_file);

#Get all genes from ALEXA (except pseudogenes)
my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'Non-pseudo')};

#Get info for each of these genes
my $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids,'-sequence'=>"no");

#Figure out which genes overlap the target regions
my @target_gene_ids = &findOverlapingGenes();

#Get exon content for only those genes overlapping the target regions
my $exon_content_ref = &getExonContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@target_gene_ids);

#Close database connection
$alexa_dbh->disconnect();


#Now go through the individual exon regions defined.  Identify those that overlap one of the target genomic regions
#At this time, also make sure each exon region identified is at least as big as the min size indicated
# - If not, increase the size.
# - This may cause overlap between adjacent regions but it should not matter since a probe will be selected if it overlaps any 1 region
my %target_exon_regions;
&identifyOverlapingExonRegions('-min_exon_region_size'=>$min_exon_region_size);


#Now there may still be overlaps between regions of different genes!!!
#Go through all the regions identified for each chromosome and merge any overlaps that remain
#At this time I will also merge and new overlaps that occured as a result of increasing exon sizes to the min
my %final_exon_regions;
&identifyFinalExonRegions();


#Finally print out the list of target exon regions
open (OUT, ">$output_file") || die "\nCould not open output file: $output_file\n\n";
print OUT "ExonRegion_ID\tGene_ID(s)\tGeneName(s)\tChromosome\tStart\tEnd\tSize\n";

foreach my $exon_region_count (sort {$final_exon_regions{$a}{start} <=> $final_exon_regions{$b}{start}} keys %final_exon_regions){
  print OUT "$exon_region_count\t$final_exon_regions{$exon_region_count}{gene_id}\t$final_exon_regions{$exon_region_count}{gene_name}\t$final_exon_regions{$exon_region_count}{chromosome}\t$final_exon_regions{$exon_region_count}{start}\t$final_exon_regions{$exon_region_count}{end}\t$final_exon_regions{$exon_region_count}{size}\n";
}
close(OUT);
close(LOG);

exit();


######################################################################################################################################
#Import the regions of interest                                                                                                      #
#- perform basic sanity checks on the input file and organize regions accorinding their target sequence (chromosome)                 #
######################################################################################################################################
sub importTargetRegions{
  my %args = @_;
  my $region_file = $args{'-region_file'};

  print YELLOW, "\n\nImporting target regions from: $region_file", RESET;
  print LOG "\n\nImporting target regions from: $region_file";

  open (REGIONS, "$region_file") || die "\nCould not open target region file: $region_file\n\n";

  my $first_line = 1;
  my $record_count = 0;

  while(<REGIONS>){
    chomp ($_);
    if ($first_line == 1){
      $first_line = 0;
      $input_header = $_;
      next();
    }
    $record_count++;
    my @line = split ("\t", $_);

    my $line = $_;
    my $region_id = $line[0];
    my $chr = $line[1];
    my $start = $line[2];
    my $end = $line[3];

    my $chromosome;
    if ($chr =~ /chr(.*)/){
      $chromosome = $1;
    }else{
      print RED, "\nChromosome name format not correct: ($chr) - check input file (line: $record_count + 1)!\n\n", RESET;
      exit();
    }

    #Check data formats:
    unless ($region_id =~ /^\d+$/){print RED, "\nRegion ID format not correct: ($region_id) - check input file (line: $record_count + 1)!\n\n", RESET; exit();}
    unless ($start =~ /^\d+$/){print RED, "\nStart coordinate format not correct: ($start) - check input file (line: $record_count + 1)!\n\n", RESET; exit();}
    unless ($end =~ /^\d+$/){print RED, "\nEnd coordinate format not correct: ($end) - check input file (line: $record_count + 1)!\n\n", RESET; exit();}

    #Sanity checks
    if ($start >= $end){
      print RED, "\nStart coordinate ($start) must be smaller than end coordinate ($end) - check input file!\n\n", RESET;
      exit();
    }
    if ($regions{$region_id}){
      print RED, "\nRegion ID: $region_id appears to be a duplicate - check input file!\n\n", RESET;
      exit();
    }

    $regions{$region_id}{record_count} = $record_count;
    $regions{$region_id}{chromosome} = $chromosome;
    $regions{$region_id}{start} = $start;
    $regions{$region_id}{end} = $end;
    $regions{$region_id}{region_size} = ($end - $start)+1;
    $regions{$region_id}{line} = $line;

    if ($sequence_names{$chromosome}){
      $sequence_names{$chromosome}{count}++;
      push (@{$sequence_names{$chromosome}{regions}}, $region_id);
      my $region_list_ref = $sequence_names{$chromosome}{region_list};
      $region_list_ref->{$region_id}->{start} = $start;
    }else{
      my @regions;
      my %region_list;
      push (@regions, $region_id);
      $region_list{$region_id}{start} = $start;
      $sequence_names{$chromosome}{count} = 1;
      $sequence_names{$chromosome}{regions} = \@regions;
      $sequence_names{$chromosome}{region_list} = \%region_list;
    }

  }
  close (REGIONS);

  my $region_count = keys %regions;
  my $seq_count = keys %sequence_names;

  print BLUE, "\n\n\tFound $region_count regions corresponding to $seq_count sequences (chromosomes)\n\n", RESET;
  print LOG "\n\n\tFound $region_count regions corresponding to $seq_count sequences (chromosomes)\n\n";

  return();
}


######################################################################################################################################
#Identify genes that overlap one or more of the target genomic regions
######################################################################################################################################
sub findOverlapingGenes{
  my @gene_ids;

  print BLUE, "\n\nSearching for genes which overlap one or more of the target genomic regions", RESET;
  print LOG "\n\nSearching for genes which overlap one or more of the target genomic regions";

  foreach my $gene_id (sort keys %{$genes_ref}){
    my $chromosome = $genes_ref->{$gene_id}->{chromosome};

    #Make sure this gene is on one of chromosomes targeted in the target regions list
    unless ($sequence_names{$chromosome}){
      next();
    }

    #Get the chromosome start and end coordinates for this gene
    my $chr_start = $genes_ref->{$gene_id}->{chr_start};
    my $chr_end = $genes_ref->{$gene_id}->{chr_end};

    #Sanity check
    if ($chr_start > $chr_end){
      print RED, "\nChr start is larger than end for gene: $gene_id\n\n", RESET;
      exit();
    }

    #Get the target regions from the chromosome
    my $chr_region_list_ref = $sequence_names{$chromosome}{region_list};

    my $overlap = 0;

    #Now go through the target region coordinates for this chromosome and see if any overlap with the current gene
    foreach my $region_id (keys %{$chr_region_list_ref}){
      my $region_start = $regions{$region_id}{start};
      my $region_end = $regions{$region_id}{end};

      if ($region_start > $region_end){
	print RED, "\nRegion start is larger than end for region: $region_id\n\n", RESET;
	exit();
      }

      if ($chr_start >= $region_start && $chr_start <= $region_end){
	$overlap = 1;
      }
      if ($chr_end >= $region_start && $chr_end <= $region_end){
	$overlap = 1;
      }
      if ($chr_start <= $region_start && $chr_end >= $region_end){
	$overlap = 1;
      }
    }
    #If an overlap was found between the current gene and any target region, add the gene to the list
    if ($overlap == 1){
      push (@gene_ids, $gene_id);
    }
  }

  my $genes_found = scalar(@gene_ids);
  print BLUE, "\n\tFound $genes_found genes which overlap one or more of the target genomic regions\n\n", RESET;
  print LOG "\n\tFound $genes_found genes which overlap one or more of the target genomic regions\n\n";

  return (@gene_ids);
}


######################################################################################################################################
sub identifyOverlapingExonRegions{
  my %args = @_;
  my $min_exon_region_size = $args{'-min_exon_region_size'};

  #1.) First calculate chromosome coordinates for all exon regions identified
  print BLUE, "\nGetting chromosome coordinates for exon regions\n\n", RESET;
  print LOG "\nGetting chromosome coordinates for exon regions\n\n";

  foreach my $gene_id (sort keys %{$exon_content_ref}){

    my $chr_strand = $genes_ref->{$gene_id}->{chr_strand};
    my $chr_start = $genes_ref->{$gene_id}->{chr_start};
    my $chr_end = $genes_ref->{$gene_id}->{chr_end};

    my $exon_regions_ref = $exon_content_ref->{$gene_id}->{exon_content};

    foreach my $exon_region_id (sort keys %{$exon_regions_ref}){

      my $start = $exon_regions_ref->{$exon_region_id}->{start};
      my $end = $exon_regions_ref->{$exon_region_id}->{end};

      #Convert the gene coordinates for this exon region to chromosome coordinates
      if ($chr_strand == 1){
	my $query_chr_start = $chr_start + $start - 1;
	my $query_chr_end = $chr_start + $end - 1;

	#Make sure the start and end are reported such that start is always smaller than end
	my $temp;
	if ($query_chr_start > $query_chr_end){
	  $temp = $query_chr_start;
	  $query_chr_start = $query_chr_end;
	  $query_chr_end = $temp;
	}

	#print YELLOW, "\n$gene_id\t$chromosome\texon: $exon_id\tstart: $query_chr_start\tend: $query_chr_end\tstrand: +", RESET;

	$exon_regions_ref->{$exon_region_id}->{chr_start} = $query_chr_start;
	$exon_regions_ref->{$exon_region_id}->{chr_end} = $query_chr_end;
	$exon_regions_ref->{$exon_region_id}->{strand} = "+";
	$exon_regions_ref->{$exon_region_id}->{size} = ($query_chr_end - $query_chr_start)+1;

      }elsif ($chr_strand == -1){

	my $query_chr_start = $chr_end - $end + 1;
	my $query_chr_end = $chr_end - $start + 1;

	#Make sure the start and end are reported such that start is always smaller than end
	my $temp;
	if ($query_chr_start > $query_chr_end){
	  $temp = $query_chr_start;
	  $query_chr_start = $query_chr_end;
	  $query_chr_end = $temp;
	}

	#print YELLOW, "\n$gene_id\t$chromosome\texon: $exon_id\tstart: $query_chr_start\tend: $query_chr_end\tstrand: -", RESET;

	$exon_regions_ref->{$exon_region_id}->{chr_start} = $query_chr_start;
	$exon_regions_ref->{$exon_region_id}->{chr_end} = $query_chr_end;
	$exon_regions_ref->{$exon_region_id}->{strand} = "-";
	$exon_regions_ref->{$exon_region_id}->{size} = ($query_chr_end - $query_chr_start)+1;

	}else{
	  print RED, "\nStrand format: $chr_strand not understood!\n\n", RESET;
	  exit();
	}
    }
  }


  #2.) Now Go through each of the exon regions identified and confirm that it overlaps a target genomic region
  #    - Also confirm that it is larger than the target size
  print BLUE, "\nDetermining which exon regions are actually within targeted genomic regions\n\n", RESET;
  print LOG "\nDetermining which exon regions are actually within targeted genomic regions\n\n";

  my $overlapping_exon_regions = 0;
  my $exon_regions_found = 0;
  my $expanded_regions = 0;
  my $bases_covered = 0;

  foreach my $gene_id (sort keys %{$exon_content_ref}){

    my $chromosome = $genes_ref->{$gene_id}->{chromosome};
    my $gene_name = $genes_ref->{$gene_id}->{gene_name};
    my $exon_regions_ref = $exon_content_ref->{$gene_id}->{exon_content};

    #Get the target regions for this chromosome
    my $chr_region_list_ref = $sequence_names{$chromosome}{region_list};

    foreach my $exon_region_id (sort keys %{$exon_regions_ref}){

      my $overlap = 0;

      my $exon_region_start = $exon_regions_ref->{$exon_region_id}->{chr_start};
      my $exon_region_end = $exon_regions_ref->{$exon_region_id}->{chr_end};
      my $exon_region_size = $exon_regions_ref->{$exon_region_id}->{size};

      #Compare the current exon region to all target genomic regions on the same chromosome
      foreach my $region_id (keys %{$chr_region_list_ref}){
	my $region_start = $regions{$region_id}{start};
	my $region_end = $regions{$region_id}{end};

	if ($exon_region_start >= $region_start && $exon_region_start <= $region_end){
	  $overlap = 1;
	}
	if ($exon_region_end >= $region_start && $exon_region_end <= $region_end){
	  $overlap = 1;
	}
	if ($exon_region_start <= $region_start && $exon_region_end >= $region_end){
	  $overlap = 1;
	}
      }

      #If overlap was found for this exon region add it to the list of %target_exon_regions
      if ($overlap == 1){
	$overlapping_exon_regions++;

	#Check the size of this region and adjust if neccessary
	if ($exon_region_size < $min_exon_region_size){
	  $expanded_regions++;
	  my $diff = $min_exon_region_size - $exon_region_size;
	  my $half_diff = ($diff/2);
	  my $half_diff_f = sprintf("%.0f", $half_diff);

	  #print YELLOW, "\nExpanding a region\tGene: $gene_id\tExon: $exon_region_id\tStart: $exon_region_start\tEnd: $exon_region_end\tSize: $exon_region_size", RESET;

	  $exon_region_start -= $half_diff_f;
	  $exon_region_end += $half_diff_f;
	  $exon_region_size = ($exon_region_end - $exon_region_start)+1;

	  #print YELLOW, "\n\tNew coordinates: Start: $exon_region_start\tEnd: $exon_region_end\tSize: $exon_region_size", RESET;

	}

	$bases_covered+= $exon_region_size;

	if ($target_exon_regions{$chromosome}){
	  my $regions_ref = $target_exon_regions{$chromosome}{regions};
	  $regions_ref->{$overlapping_exon_regions}->{gene_id} = $gene_id;
	  $regions_ref->{$overlapping_exon_regions}->{gene_name} = $gene_name;
	  $regions_ref->{$overlapping_exon_regions}->{chromosome} = $chromosome;
	  $regions_ref->{$overlapping_exon_regions}->{start} = $exon_region_start;
	  $regions_ref->{$overlapping_exon_regions}->{end} = $exon_region_end;
	  $regions_ref->{$overlapping_exon_regions}->{size} = $exon_region_size;

	  $exon_regions_found++;
	}else{
	  my %regions;
	  $regions{$overlapping_exon_regions}{gene_id} = $gene_id;
	  $regions{$overlapping_exon_regions}{gene_name} = $gene_name;
	  $regions{$overlapping_exon_regions}{chromosome} = $chromosome;
	  $regions{$overlapping_exon_regions}{start} = $exon_region_start;
	  $regions{$overlapping_exon_regions}{end} = $exon_region_end;
	  $regions{$overlapping_exon_regions}{size} = $exon_region_size;
	  $target_exon_regions{$chromosome}{regions} = \%regions;
	  $exon_regions_found++;
	}
      }
    }
  }


  print BLUE, "\n\nFound $exon_regions_found exon regions within the targeted genomic regions", RESET;
  print BLUE, "\n\tA total of $expanded_regions of these had to be expanded to the min size of $min_exon_region_size", RESET;
  print BLUE, "\n\tAll exon regions combined cover a total of $bases_covered bases (not taking anti-sense or overlaps into account)\n\n", RESET;

  print LOG "\n\nFound $exon_regions_found exon regions within the targeted genomic regions";
  print LOG "\n\tA total of $expanded_regions of these had to be expanded to the min size of $min_exon_region_size";
  print LOG "\n\tAll exon regions combined cover a total of $bases_covered bases (not taking anti-sense or overlaps into account)\n\n";

  return();
}


######################################################################################################################################
#identifyFinalExonRegions                                                                                                            #
######################################################################################################################################
sub identifyFinalExonRegions{
  my %args = @_;
  my $regions_processed = 0;
  my $merged_regions = 0;

  print BLUE, "\nChecking for overlaping regions and merging them\n\n", RESET;

  my $final_exon_count = 0;

  foreach my $chr (sort keys %target_exon_regions){

    print BLUE, "\n\tProcessing regions on $chr", RESET;

    #Get the regions for this chromosome
    my $regions_list_ref = $target_exon_regions{$chr}{regions};

    #Define clusters of overlapping regions for this chromosome
    my %clusters;
    my $cluster_count = 0;

    #Sort the regions by start position and go through them one by one comparing coordinates to all previously considered regions
    foreach my $region_id (sort {$regions_list_ref->{$a}->{start} <=> $regions_list_ref->{$b}->{start}} keys %{$regions_list_ref}){

      my $current_start = $regions_list_ref->{$region_id}->{start};
      my $current_end = $regions_list_ref->{$region_id}->{end};

      my $overlap_found = 0;
      my $overlaping_cluster = '';

      #Go through each cluster defined so far and check for overlaps to the current region
    CLUSTER:foreach my $cluster_id (sort {$a <=> $b} keys %clusters){

	my $ol_regions_ref = $clusters{$cluster_id}{ol_regions};

	foreach my $test_region_id (sort {$a <=> $b} keys %{$ol_regions_ref}){

	  my $test_start = $regions_list_ref->{$test_region_id}->{start};
	  my $test_end = $regions_list_ref->{$test_region_id}->{end};

	  if (($test_start > $current_start) && ($test_start < $current_end)){
	    $overlap_found = 1;
	  }
	  if (($test_end > $current_start) && ($test_end < $current_end)){
	    $overlap_found = 1;
	  }
	  if (($test_start <= $current_start) && ($test_end >= $current_end)){
	    $overlap_found = 1;
	  }
	  if ($overlap_found == 1){
	    $overlaping_cluster = $cluster_id;
	    last CLUSTER;
	  }
	}
      }


      if ($overlap_found == 1){
	#If this region does overlap one of the existing clusters, add it
	my $ol_regions_ref = $clusters{$overlaping_cluster}{ol_regions};
	$ol_regions_ref->{$region_id}->{start} = $regions_list_ref->{$region_id}->{start};
	$ol_regions_ref->{$region_id}->{end} = $regions_list_ref->{$region_id}->{end};

      }else{
	#If this region does not overlap any of the existing clusters, create a new cluster
	$cluster_count++;
	my %ol_regions;
	$ol_regions{$region_id}{start} = $regions_list_ref->{$region_id}->{start};
	$ol_regions{$region_id}{end} = $regions_list_ref->{$region_id}->{end};
	
	$clusters{$cluster_count}{ol_regions} = \%ol_regions;
	$clusters{$cluster_count}{chromosome} = $chr;
	$clusters{$cluster_count}{large_region} = "no";
	
      }

    }

    #DEBUG: Print dumper of clusters for this chromosome
    #print Dumper %clusters;

    #Go through each of the clusters defined for this chromosome and print new merged regions
    foreach my $cluster_id (sort {$a <=> $b} keys %clusters){

      #If a cluster consists of only a single region, simply keep the original record
      my $ol_regions_ref = $clusters{$cluster_id}{ol_regions};
      my $cluster_size = keys %{$ol_regions_ref};
      $regions_processed += $cluster_size;

      if ($cluster_size == 1){
	#Process cluster of only a single region
	foreach my $region_id (sort keys %{$ol_regions_ref}){
	  $final_exon_count++;

	  $final_exon_regions{$final_exon_count}{gene_id} = $regions_list_ref->{$region_id}->{gene_id};
	  $final_exon_regions{$final_exon_count}{gene_name} = $regions_list_ref->{$region_id}->{gene_name};
	  $final_exon_regions{$final_exon_count}{chromosome} = $regions_list_ref->{$region_id}->{chromosome};
	  $final_exon_regions{$final_exon_count}{start} = $regions_list_ref->{$region_id}->{start};
	  $final_exon_regions{$final_exon_count}{end} = $regions_list_ref->{$region_id}->{end};
	  $final_exon_regions{$final_exon_count}{size} = $regions_list_ref->{$region_id}->{size};

	  $merged_regions++;
	}

      }else{
	#Process actual cluster of overlapping regions
	#First get the min region ID of the cluster as well as the new start/end coordinates of the merged region

	my @ids;
	my @starts;
	my @ends;

	my %unique_gene_ids;
	my %unique_gene_names;

	foreach my $region_id (sort keys %{$ol_regions_ref}){
	  push (@ids, $region_id);
	  push (@starts, $regions_list_ref->{$region_id}->{start});
	  push (@ends, $regions_list_ref->{$region_id}->{end});
	  my $gene_id = $regions_list_ref->{$region_id}->{gene_id};
	  my $gene_name = $regions_list_ref->{$region_id}->{gene_name};
	  $unique_gene_ids{$gene_id}{tmp} = '';
	  $unique_gene_names{$gene_name}{tmp} = '';
	}
	my @unique_gene_ids = keys %unique_gene_ids;
	my @unique_gene_names = keys %unique_gene_names;;

	my $stat;

	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@ids);
	my $min_region_id = $stat->min();

	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@starts);
	my $merged_start = $stat->min();

	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@ends);
	my $merged_end = $stat->max();

	$final_exon_count++;
	$final_exon_regions{$final_exon_count}{gene_id} = "@unique_gene_ids";
	$final_exon_regions{$final_exon_count}{gene_name} = "@unique_gene_names";
	$final_exon_regions{$final_exon_count}{chromosome} = $chr;
	$final_exon_regions{$final_exon_count}{start} = $merged_start;
	$final_exon_regions{$final_exon_count}{end} = $merged_end;
	$final_exon_regions{$final_exon_count}{size} = ($merged_end - $merged_start)+1;

	$merged_regions++;
      }
    }
  }

  print BLUE, "\n\nProcessed a total of $regions_processed regions to identify overlaps", RESET;
  print BLUE, "\nIdentified a total of $final_exon_count merged exons regions\n\n", RESET;
  return();
}
