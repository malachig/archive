#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to examine coverage of exonic bases of an EnsEMBL transcriptome as a function of Illumina read depth

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use Benchmark;

#Load the ALEXA libraries
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
my $read_records_dir = '';
my $min_bit_score = '';          #Minimum bit score of BLAST hit for each read to be allowed for read-to-gene mappings
my $block_size = '';             #Number of reads before coverage is summarized
my $cutoffs = '';                #Coverage cutoffs against which genes and exonic bases will be assessed to be considered expressed/detected
my $outfile = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'read_records_dir=s'=>\$read_records_dir, 'min_bit_score=f'=>\$min_bit_score,
	    'block_size=i'=>\$block_size, 'cutoffs=s'=>\$cutoffs, 'outfile=s'=>\$outfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password", RESET;
print GREEN, "\n\tspecify a path to the read record files to be summarized using: --read_records_dir\n", RESET;
print GREEN, "\n\tSpecify the minimum bit score to a transcript for each read to be considered for the summary using: --min_bit_score", RESET;
print GREEN, "\n\tSpecify the number of reads to process each time before summarizing coverage using: --block_size", RESET;
print GREEN, "\n\tSpecify a list of minimum cutoffs using: --cutoffs", RESET;
print GREEN, "\n\tSpecify the outfile summarizing coverage after each read block using: --outfile", RESET;

print GREEN, "\n\nExample: coverageVsDepth_Gene-ExonicBase_MULTI.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=username  --password=pwd  --read_records_dir=/projects/malachig/solexa/read_records/HS04391/ENST_v53/  --min_bit_score=60.0  --block_size=100000  --cutoffs='1,5,10,50,100,500'  --outfile=/projects/malachig/solexa/read_records/HS04391/Summary/coverageVsDepth/HS04391_Lanes1-23_Gene-ExonicBases_Coverage_1-5-10-50-100-500xCoverage.txt\n\n", RESET;

unless ($database && $server && $user && $password && $read_records_dir && $min_bit_score && ($block_size =~ /\d+/) && $cutoffs && $outfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Get list of input files containing data from the specified directory
$read_records_dir = &checkDir('-dir'=>$read_records_dir, '-clear'=>"no");
my %files;
&getDataFiles('-input_dir'=>$read_records_dir);

#Get the process ID for this script
my $pid = $$;

open (OUT, ">$outfile") || die "\nCould not open output file: $outfile\n\n";

#Get the cutoffs and build a global cutoffs hash
my @cutoffs = split(",", $cutoffs);
unless(scalar(@cutoffs) > 0){
  print RED, "\nList of cutoffs not understood.  Must be supplied as a single number or a comma separated list\n\n", RESET;
  exit();
}
my %cutoffs;
foreach my $cut (@cutoffs){
  $cutoffs{$cut}{last_bases_covered} = 0;
  $cutoffs{$cut}{bases_covered} = 0;
  $cutoffs{$cut}{expressed_gene_count} = 0;
  $cutoffs{$cut}{last_expressed_gene_count} = 0;
  $cutoffs{$cut}{bases_diff} = 0;
  $cutoffs{$cut}{genes_diff} = 0;
  $cutoffs{$cut}{exonic_name} = "ExonicBases_"."$cut"."X";
  $cutoffs{$cut}{gene_name} = "Genes_"."$cut"."X";
}

#Get cutoffs string
my @tmp_e;
my @tmp_g;
foreach my $cutoff (sort {$a <=> $b} keys %cutoffs){
  push(@tmp_e, $cutoffs{$cutoff}{exonic_name});
  push(@tmp_g, $cutoffs{$cutoff}{gene_name});
}
my $exonic_string = join("\t", @tmp_e);
my $gene_string = join("\t", @tmp_g);

#0.) First get all the neccessary gene info required to perform the analysis
my $genes_ref;
my $gene_transcripts_ref;
my $gene_exon_content_ref;
my $chr_genes_ref;
my $last_gene_count = 0;
my $last_bases_count = 0;

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};
#my @gene_ids = (1 .. 100);

&getBasicGeneInfo('-gene_list'=>\@gene_ids);

#Close database connection
$alexa_dbh->disconnect();


#1.) Now initialize arrays to store coverage info
my $chr_exon_coverage_ref;
my $total_exon_bases = 0;

print OUT "#User specified parameters:\n";
print OUT "#min_bit_score = $min_bit_score\n";
print OUT "#cutoffs  = $cutoffs\n";
print OUT "#read_records_dir = $read_records_dir\n";
print OUT "#database = $database\n";

$| = 1; print BLUE, "\n\n2.) Determining the non-redundant exonic base coverage for ALL EnsEMBL genes\n", RESET; $| = 0;
&getExonCoverageObject('-gene_list'=>\@gene_ids);

print OUT "#MappedReads\t$gene_string\t$exonic_string\n";

$| = 1; print BLUE, "\n\n3-a.) Now parsing through read records files and reporting coverage after each block of $block_size reads is added\n", RESET; $| = 0;

my $total_reads_parsed = 0;
my $mapped_reads = 0;
my $grand_mapped_reads = 0;

#Report the processing time for each block of $block_size reads processed
my $t0 = new Benchmark;

#Go through each read record file and parse.
foreach my $file_count (sort {$files{$a}->{random} <=> $files{$b}->{random}} keys %files){

  print GREEN, "\nProcessing: $files{$file_count}{file_name}", RESET;

  open (READ, "zcat $files{$file_count}{file_path} |") || die "\nCould not open read records summary infile: $files{$file_count}{file_path}\n\n";
  my $header_line = 1;
  my %columns = %{$files{$file_count}{columns}};

  while(<READ>){
    chomp($_);
    my $current_line = $_;
    my @line = split("\t", $_);

    if ($header_line == 1){
      $header_line = 0;
      next();
    }

    $total_reads_parsed++;

    my $read_id = $line[$columns{Read_ID}{position}];
    my $r1_id = $line[$columns{R1_ID}{position}];
    my $r2_id = $line[$columns{R2_ID}{position}];
    my $r1_hit_type = $line[$columns{R1_HitType}{position}];
    my $r2_hit_type = $line[$columns{R2_HitType}{position}];
    my $r1_gene_id = $line[$columns{R1_GeneID}{position}];
    my $r2_gene_id = $line[$columns{R2_GeneID}{position}];
    my $r1_bit_score = $line[$columns{R1_BitScore}{position}];
    my $r2_bit_score = $line[$columns{R2_BitScore}{position}];
    my $r1_chromosome = $line[$columns{R1_Chromosome}{position}];
    my $r2_chromosome = $line[$columns{R2_Chromosome}{position}];
    my $r1_chr_start_coords = $line[$columns{R1_ChrStartCoords}{position}];
    my @r1_chr_start_coords = split(" ", $r1_chr_start_coords);
    my $r2_chr_start_coords = $line[$columns{R2_ChrStartCoords}{position}];
    my @r2_chr_start_coords = split(" ", $r2_chr_start_coords);
    my $r1_chr_end_coords = $line[$columns{R1_ChrEndCoords}{position}];
    my @r1_chr_end_coords = split(" ", $r1_chr_end_coords);
    my $r2_chr_end_coords = $line[$columns{R2_ChrEndCoords}{position}];
    my @r2_chr_end_coords = split(" ", $r2_chr_end_coords);

    #Sanity check of coordinates
    unless (scalar(@r1_chr_start_coords) == scalar(@r1_chr_end_coords)){
      print RED, "\nRead: $r1_id does not have an equal number of start and end coords!\n", RESET;
      print RED, "\nLine: $current_line\n\n", RESET;
      exit();
    }
    unless (scalar(@r2_chr_start_coords) == scalar(@r2_chr_end_coords)){
      print RED, "\nRead: $r2_id does not have an equal number of start and end coords!\n", RESET;
      print RED, "\nLine: $current_line\n\n", RESET;
      exit();
    }

    #change bit scores of 'NA' to 0
    if ($r1_bit_score eq "NA"){$r1_bit_score = 0;}
    if ($r2_bit_score eq "NA"){$r2_bit_score = 0;}

    #If both read alignments are too short, skip this record immediately
    unless ($r1_bit_score >= $min_bit_score || $r2_bit_score >= $min_bit_score){
      next();
    }

    #Test Read1 and Read2 to see if they pass the quality threshold individually
    my $read1_passes = 0;
    my $read2_passes = 0;

    if (($r1_hit_type eq "Top_Hit") && ($r1_bit_score >= $min_bit_score)){
      $read1_passes = 1;
    }
    if (($r2_hit_type eq "Top_Hit") && ($r2_bit_score >= $min_bit_score)){
      $read2_passes = 1;
    }

    #Fix chromosome formats
    if ($r1_chromosome eq "chrMT"){$r1_chromosome = "chrM";}
    if ($r2_chromosome eq "chrMT"){$r2_chromosome = "chrM";}

    my $r1_chr;
    my $r2_chr;
    if ($r1_chromosome =~ /chr(\S+)/){
      $r1_chr = $1;
    }
    if ($r2_chromosome =~ /chr(\S+)/){
      $r2_chr = $1;
    }


    #Add a count for the gene that this read hits.
    #Deal with READ1
    if ($read1_passes == 1){
      $mapped_reads++;
      $grand_mapped_reads++;

      #Make sure the gene ID found is valid
      if ($genes_ref->{$r1_gene_id}){

        #Count this read hit to the gene
        $genes_ref->{$r1_gene_id}->{quality_read_count}++;

        #Add the coverage of this read to the gene
        &addReadCoverage('-gene_id'=>$r1_gene_id, '-chr'=>$r1_chr, '-chr_starts'=>\@r1_chr_start_coords, '-chr_ends'=>\@r1_chr_end_coords);

      }else{
        print RED, "\nCould not find gene id: $r1_gene_id in Gene info object!\n\n", RESET;
        exit();
      }
    }

    #Deal with READ2
    if ($read2_passes == 1){
      $mapped_reads++;
      $grand_mapped_reads++;

      #Make sure the gene ID found is valid
      if ($genes_ref->{$r2_gene_id}){

        #Count this read hit to the gene
        $genes_ref->{$r2_gene_id}->{quality_read_count}++;

        #Add the coverage of this read to the gene
        &addReadCoverage('-gene_id'=>$r2_gene_id, '-chr'=>$r2_chr, '-chr_starts'=>\@r2_chr_start_coords, '-chr_ends'=>\@r2_chr_end_coords);

      }else{
        print RED, "\nCould not find gene id: $r2_gene_id in Gene info object!\n\n", RESET;
        exit();
      }
    }

    #After every N successfully mapped reads processed:
    # - calculate the overall coverage of the transcriptome attained
    # - calculate the number of genes with at least Y reads

    if ($mapped_reads >= $block_size){
      &summarizeCoverage();
      $mapped_reads = 0;
    }
  }
  close (READ);
}


#Calculate the exon coverage for all the exons of those genes that were actually detected as expressed at a level of 10X or greater!
$total_exon_bases = 0;
my %genes_found;
foreach my $gene_id (keys %{$genes_ref}){
  if ($genes_ref->{$gene_id}->{quality_read_count} >= 10){
    $genes_found{$gene_id}{tmp} = '';
  }
}

my $found_gene_count = keys %genes_found;
$| = 1; print BLUE, "\n\n Found a total of $found_gene_count genes expressed at >= 10X.  Determining the coverage for these genes only", RESET; $| = 0;
print OUT "#Found a total of $found_gene_count genes expressed at >= 10X.  Determining the coverage for these genes only\n";

my @found_gene_ids = keys %genes_found;
&getExonCoverageObject('-gene_list'=>\@found_gene_ids);

print OUT "#SCRIPT COMPLETE\n";
close (OUT);

exit();


###########################################################################################################
#Get data files and the columns of each                                                                   #
###########################################################################################################
sub getDataFiles{
  my %args = @_;
  my $dir = $args{'-input_dir'};

  my @required_columns = qw(Read_ID DistanceBetweenReads_Genomic DistanceBetweenReads_Transcript R1_ID R1_HitType R1_GeneID R1_BitScore R1_Chromosome R1_ChrStartCoords R1_ChrEndCoords R2_ID R2_HitType R2_GeneID R2_BitScore R2_Chromosome R2_ChrStartCoords R2_ChrEndCoords);

  my $dh = opendir(DIR, $dir) || die "\nCould not open directory: $dir\n\n";

  my @files = readdir(DIR);
  
  #Assign a random number to each file to allow random order processing of lane files
  srand();

  my $count = 0;
  foreach my $file (@files){
    my %columns;
    my $header = 1;
    
    chomp($file);
    unless ($file =~ /\.txt\.gz$/){
      next();
    }
    if (-d $file){
      next();
    }
    $count++;

    $files{$count}{file_name} = $file;
    $files{$count}{file_path} = $dir.$file;
    $files{$count}{random} = rand();

    #Get the header values for this file
    open (FILE, "zcat $dir$file |") || die "\nCould not open file: $dir$file";

    while(<FILE>){
      chomp($_);
      my @line = split("\t", $_);

      #Parse the column names and positions.  Check against a hard coded list of required columns before proceeding
      if ($header == 1){
        my $col_count = 0;
        foreach my $column (@line){
          $columns{$column}{position} = $col_count;
          $col_count++;
        }
        foreach my $req_column (@required_columns){
          unless ($columns{$req_column}){
	    print RED, "\nRequired column: $req_column was not found in the file: $file\n\n", RESET;
	    exit();
          }
        }
        last();
      }

    }
    close(FILE);
    $files{$count}{columns} = \%columns;
  }
  closedir(DIR);

  my $files_count = keys %files;
  print BLUE, "\n\nFound $files_count files to be processed (all .txt.gz files in the specified directory)", RESET;

  return();
}


############################################################################################################################################
#Get basic info for all genes from the user specified ALEXA database                                                                       #
############################################################################################################################################
sub getBasicGeneInfo{
  my %args = @_;
  my @gene_list = @{$args{'-gene_list'}};

  my %chr_genes;

  #Get gene info for all genes
  $| = 1; print BLUE, "\n1-a.) Getting gene data", RESET; $| = 0;
  my $genes_storable = "$database"."_AllGenes_GeneInfo_NoSeq.storable";
  $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_list, '-sequence'=>"no", '-storable'=>$genes_storable);

  #Get the transcript info for all transcripts of these genes
  $| = 1; print BLUE, "\n1-b.) Getting transcript data as well as exons for each transcript", RESET; $| = 0;
  my $trans_storable = "$database"."_AllGenes_TranscriptInfo_NoSeq.storable";
  $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_list, '-sequence'=>"no", '-storable'=>$trans_storable);

  #Get chromosome coordinates for all EnsEMBL transcripts
  $| = 1; print BLUE, "\n\n1-c.) Calculating chromosome coordinates for the EXONS of each gene", RESET; $| = 0;
  foreach my $gene_id (keys %{$gene_transcripts_ref}){

    #Initialize quality read_count value
    $genes_ref->{$gene_id}->{quality_read_count} = 0;


    my $chromosome = $genes_ref->{$gene_id}->{chromosome};

    if ($chromosome eq "MT"){$chromosome = "M";}

    #Add this gene to its chromosome
    if ($chr_genes{$chromosome}){
      push (@{$chr_genes{$chromosome}{gene_list}}, $gene_id);
    }else{
      my @genes;
      push (@genes, $gene_id);
      $chr_genes{$chromosome}{gene_list} = \@genes;
    }

    my $ucsc_chromosome = "chr"."$genes_ref->{$gene_id}->{chromosome}";

    #$genes_ref->{$gene_id}->{chromosome} = $ucsc_chromosome;

    $genes_ref->{$gene_id}->{ucsc_chromosome} = $ucsc_chromosome;
    my $chr_strand = $genes_ref->{$gene_id}->{chr_strand};
    my $chr_start = $genes_ref->{$gene_id}->{chr_start};
    my $chr_end = $genes_ref->{$gene_id}->{chr_end};
    my $gene_start = $genes_ref->{$gene_id}->{gene_start};
    my $gene_end = $genes_ref->{$gene_id}->{gene_end};

    my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};

    foreach my $trans_id (keys %{$transcripts_ref}){
      my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

      foreach my $exon_id (keys %{$exons_ref}){

	my $start = $exons_ref->{$exon_id}->{exon_start};
	my $end = $exons_ref->{$exon_id}->{exon_end};

        #Convert provided gene coordinates to coordinates relative to the chromosome
        my $coords_ref = &convertGeneCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$start, '-end_pos'=>$end, '-ordered'=>"yes");
        $exons_ref->{$exon_id}->{chr_start} = $coords_ref->{$gene_id}->{chr_start};
	$exons_ref->{$exon_id}->{chr_end} = $coords_ref->{$gene_id}->{chr_end};
	$exons_ref->{$exon_id}->{strand} = $coords_ref->{$gene_id}->{strand};

      }
    }
  }

  $chr_genes_ref = \%chr_genes;

  #Get exon content for all genes
  $| = 1; print BLUE, "\n\n1-d.) Getting EXON CONTENT of each gene", RESET; $| = 0;
  my $storable_name = "$database"."_AllGenes_ExonContent.storable";
  $gene_exon_content_ref = &getExonContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_list, '-storable'=>$storable_name);

  return();
}


############################################################################################################################################
#Initialize hashes to store Exon coverage for a subset of genes
############################################################################################################################################
sub getExonCoverageObject{
  my %args = @_;
  my @gene_list = @{$args{'-gene_list'}};

  my %gene_exon_coverage;
  my %chr_exon_coverage;

  #At this time, also get the chromosome coordinates for all exon-content coordinates
  $| = 1; print BLUE, "\n\na.) Calculating chromosome coordinates for the EXON CONTENT of each chromosome\n", RESET; $| = 0;

  my $counter = 0;
  foreach my $gene_id (@gene_list){

    $counter++;
    if ($counter == 100){
      $counter = 0;
      $| = 1; print BLUE, ".", RESET; $| = 0;
    }

    my $chromosome = $genes_ref->{$gene_id}->{chromosome};
    my $chr_strand = $genes_ref->{$gene_id}->{chr_strand};
    my $chr_start = $genes_ref->{$gene_id}->{chr_start};
    my $chr_end = $genes_ref->{$gene_id}->{chr_end};
    my $gene_start = $genes_ref->{$gene_id}->{gene_start};
    my $gene_end = $genes_ref->{$gene_id}->{gene_end};

    #Calculate the size of each transcript by adding up the size of its exons
    my $size = 0;
    my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};
    my %exon_content_chr;

    $gene_exon_coverage{$gene_id}{exon_content_size} = 0; #Total number of bases covered by exons of this gene

    foreach my $exon_id (keys %{$exon_content_ref}){

      my $start = $exon_content_ref->{$exon_id}->{start};
      my $end = $exon_content_ref->{$exon_id}->{end};

      my $coords_ref = &convertGeneCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$start, '-end_pos'=>$end, '-ordered'=>"yes");
      $exon_content_chr{$exon_id}{chr_start} = $coords_ref->{$gene_id}->{chr_start};
      $exon_content_chr{$exon_id}{chr_end} = $coords_ref->{$gene_id}->{chr_end};
      $exon_content_chr{$exon_id}{strand} = $coords_ref->{$gene_id}->{strand};
      $exon_content_chr{$exon_id}{size} = ($coords_ref->{$gene_id}->{chr_end} - $coords_ref->{$gene_id}->{chr_start})+1;
      $gene_exon_coverage{$gene_id}{exon_content_size} += ($coords_ref->{$gene_id}->{chr_end} - $coords_ref->{$gene_id}->{chr_start})+1;
    }
    $gene_exon_coverage{$gene_id}{exon_content} = \%exon_content_chr;
  }


  #Build a coverage hash for chromosome (keyed on chromosome positions)
  #At this time, also get the chromosome coordinates for all exon-content coordinates
  $| = 1; print BLUE, "\nb.) Initializing a COVERAGE hash for the EXON CONTENT of each chromosome\n", RESET; $| = 0;

  foreach my $current_chr (sort keys %{$chr_genes_ref}){

    $counter = 0;
    my $highest_coord = 0;
    my $lowest_coord = 100000000000000000000000000000000000000000000;

    my @chr_genes = @{$chr_genes_ref->{$current_chr}->{gene_list}};
    my $chr_gene_count = scalar(@chr_genes);

    $| = 1; print BLUE, "\n\tProcessing $current_chr (has $chr_gene_count genes)", RESET; $| = 0;

    my %coverage;

    foreach my $gene_id (@chr_genes){

      $counter++;
      if ($counter == 100){
	$counter = 0;
	$| = 1; print BLUE, ".", RESET; $| = 0;
      }

      my $exon_content_ref = $gene_exon_coverage{$gene_id}{exon_content};

      foreach my $exon_id (sort {$exon_content_ref->{$a}->{chr_start} <=> $exon_content_ref->{$b}->{chr_start}} keys %{$exon_content_ref}){

	my $start = $exon_content_ref->{$exon_id}->{chr_start};
	my $end = $exon_content_ref->{$exon_id}->{chr_end};

	if ($end > $highest_coord){$highest_coord = $end;}
	if ($start < $lowest_coord){$lowest_coord = $start;}

	#Go through each chromosome position in this exon content block and initialize that position in the hash
	for (my $i = $start; $i <= $end; $i++){
	  $coverage{$i} = 0;
	}
      }
    }

    my $chr_exon_bases = keys %coverage;
    $total_exon_bases += $chr_exon_bases;
    my $coord_span = ($highest_coord - $lowest_coord)+1;

    $chr_genes_ref->{$current_chr}->{exonic_base_count} = $chr_exon_bases;
    $| = 1; print BLUE, "\n\t\tFound $chr_exon_bases exonic bases covering a total of $coord_span chromosome bases", RESET; $| = 0;

    #initialize an empty coverage hash
    my %temp;
    $chr_genes_ref->{$current_chr}->{coverage} = \%temp;

    %coverage = ();
  }

  my $chr_count = keys %{$chr_genes_ref};
  my $gene_count = scalar(@gene_list);

  $| = 1; print BLUE, "\n\n\tFound $total_exon_bases non-redundant exonic bases in $gene_count genes on $chr_count chromosomes/contigs", RESET; $| = 0;
  print OUT "#Found $total_exon_bases non-redundant exonic bases in $gene_count genes on $chr_count chromosomes/contigs\n";

  return(\%chr_exon_coverage);
}


############################################################################################################################################
#Add the coverage of a read to a gene to the coverage hash for the exon content record for that gene
############################################################################################################################################
sub addReadCoverage{
  my %args = @_;
  my $gene_id = $args{'-gene_id'};
  my $chr = $args{'-chr'};
  my @chr_starts = @{$args{'-chr_starts'}};
  my @chr_ends = @{$args{'-chr_ends'}};

  #A.) First store coverage at the level of entire chromosomes
  my $chr_coverage_ref = $chr_genes_ref->{$chr}->{coverage};
  my @ends = @chr_ends;
  foreach my $start (@chr_starts){
    my $end = shift(@ends);

    #Go through each chromosome position in this read as it is mapped to an exon and increment that position in the hash
    for (my $i = $start; $i <= $end; $i++){
      if ($chr_coverage_ref->{$i}){
	$chr_coverage_ref->{$i}++;
      }else{
	$chr_coverage_ref->{$i} = 1;
      }
    }
  }

  #B.) Next store coverage at the level of single genes (to allow calculations for specific exons)



  return();
}


################################################################################################################################
#summarize read coverage of transcriptome and number of genes identified after given number of reads
################################################################################################################################
sub summarizeCoverage{

  #Reset coverage count for each cutoff
  foreach my $cutoff (keys %cutoffs){
    $cutoffs{$cutoff}{bases_covered} = 0;
  }

  #Calculate the exonic base coverage for the whole transcriptome - at each cutoff level
  foreach my $chr (sort keys %{$chr_genes_ref}){
    my $coverage_ref = $chr_genes_ref->{$chr}->{coverage};
    foreach my $pos (keys %{$coverage_ref}){
      my $cov = $coverage_ref->{$pos};
      foreach my $cutoff (sort {$a <=> $b} keys %cutoffs){
        if ($cov >= $cutoff){
          $cutoffs{$cutoff}{bases_covered}++;
        }
      }
    }
  }

  #Calculate the number of genes detected as expressed - at each cutoff level
  foreach my $cutoff (sort {$a <=> $b} keys %cutoffs){
    my %genes_expressed;
    foreach my $gene_id (keys %{$genes_ref}){
      if ($genes_ref->{$gene_id}->{quality_read_count} >= $cutoff){
        $genes_expressed{$gene_id}{tmp} = '';
      }
    }
    my $expressed_gene_count = keys %genes_expressed;
    $cutoffs{$cutoff}{expressed_gene_count} = $expressed_gene_count;
  }

  #Calculate the difference between this iteration and the last - at each cutoff level
  foreach my $cutoff (sort {$a <=> $b} keys %cutoffs){

    my $genes_diff = $cutoffs{$cutoff}{expressed_gene_count} - $cutoffs{$cutoff}{last_expressed_gene_count};
    $cutoffs{$cutoff}{last_expressed_gene_count} = $cutoffs{$cutoff}{expressed_gene_count};
    $cutoffs{$cutoff}{genes_diff} = $genes_diff;

    my $bases_diff = $cutoffs{$cutoff}{bases_covered} - $cutoffs{$cutoff}{last_bases_covered};
    $cutoffs{$cutoff}{last_bases_covered} = $cutoffs{$cutoff}{bases_covered};
    $cutoffs{$cutoff}{bases_diff} = $bases_diff;
  }

  #Generate display strings for each cutoff for genes and bases
  my @expressed_bases;
  my @expressed_bases_diff;
  my @expressed_genes;
  my @expressed_genes_diff;
  foreach my $cutoff (sort {$a <=> $b} keys %cutoffs){
    push(@expressed_bases, $cutoffs{$cutoff}{bases_covered});
    push(@expressed_bases_diff, $cutoffs{$cutoff}{bases_diff});
    push(@expressed_genes, $cutoffs{$cutoff}{expressed_gene_count});
    push(@expressed_genes_diff, $cutoffs{$cutoff}{genes_diff});
  }
  my $expressed_bases = join("\t", @expressed_bases);
  my $expressed_bases_diff = join ("\t", @expressed_bases_diff);
  my $expressed_genes = join("\t", @expressed_genes);
  my $expressed_genes_diff = join ("\t", @expressed_genes_diff);

  #Determine elapsed time since last report:
  my $t1 = new Benchmark;
  my $td = timediff($t1, $t0);
  $t0 = $t1;

  #Determine current %memory usage
  my $ps_query = `ps u -p $pid`;
  my @process_info = split ("\n", $ps_query);
  my $memory_usage = '';
  if ($process_info[1] =~ /\S+\s+\S+\s+\S+\s+(\S+)\s+/){
    $memory_usage = $1;
  } 

  my $string = timestr($td);
  my $seconds = '';
  if ($string =~ /(\d+)\s+wallclock/){
    $seconds = $1;
  } 

  $| = 1; print BLUE, "\n\nREADS: $grand_mapped_reads\nGENES:\n\tCumulative: $expressed_genes\n\tNew: $expressed_genes_diff\nBASES:\n\tCumulative: $expressed_bases\n\tNew: $expressed_bases_diff", RESET; $| = 0; 
  $| = 1; print BLUE, "\n[Elapsed Time = $seconds seconds]\t[% memory usage = $memory_usage]", RESET; $| = 0; 
  print OUT "$grand_mapped_reads\t$expressed_genes\t$expressed_bases\n";

  return();
}
