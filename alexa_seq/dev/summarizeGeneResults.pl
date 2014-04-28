#!/usr/bin/perl -w

#Written by Malachi Griffith
#The purpose of this script is to summarize the read hits to each gene (any transcript of that gene)
#This script takes as input a read records file created by parsing BLAST results of reads to an EnsEMBL/ALEXA transcriptome

#0.) Get neccessary gene info from a user specified ALEXA database

#Parse a read record file and do the following on the fly line-by-line
#All information required for the summary should be available on a single line of the read records summary file
#Each of these lines contains information on both members of a read pair including things like top hits, distance between reads, chromosome coordinates, etc.

#1.) Summarize Gene-level hits
#1-A.) Use this script to get the following for each gene and create a gene level summary file:
#      - The number of hits per gene
#      - % coverage of the gene
#      - The average 'X' coverage of the gene
#      - The name and description of the gene

#      - To summarize coverage, for each gene build a hash keyed on chromosome position by going through the exon content for that gene
#      - When building this hash, initialize each position to a count of 0
#      - Then go through the chromosome coordinates of each read and increment the appropriate positions in the coverage hash

#1-B.) Create custom UCSC tracks (one per chromosome) to display the coverage for each gene
#      - Create as a 'Wiggle'
#      - The wiggle track should help in these regions to see how even the coverage is...
#      - Use the same hashes generated to calculate % coverage and average 'X' coverage to generate this wig track
#1-C.) NOTE: For both (A) and (B) impose a quality filter when summarizing:
#      - Each read must be a 'Top_Hit' (not ambiguous)
#      - Each read must have an bit score of at least X (48.1 maybe...)
#1-D.) Create a gene-level summary file


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

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $library = '';
my $read_records_infile = '';
my $min_bit_score = '';    #Minimum bit score to be allowed for read-to-gene mappings
my $gene_description_file = '';
my $gene_summary_file = '';
my $ucsc_dir = '';
my $web_path = '';


GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'library=s'=>\$library, 'read_records_infile=s'=>\$read_records_infile,
	    'min_bit_score=f'=>\$min_bit_score,
	    'gene_description_file=s'=>\$gene_description_file, 'gene_summary_file=s'=>\$gene_summary_file,
	    'ucsc_dir=s'=>\$ucsc_dir, 'web_path=s'=>\$web_path);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script summarizes the data in a Read Records Summary file", RESET;
print GREEN, "\n\tThis script takes a read record summary as input, summarizes the results and creates new summary and ucsc track files", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password", RESET;
print GREEN, "\n\tMake sure this database matches the one that was used to generate the BLAST results in the first place!!\n", RESET;
print GREEN, "\n\tSpecify the library name for read data using: --library", RESET;
print GREEN, "\n\tspecify a path to the read records summary file to be summarized using: --read_records_file\n", RESET;
print GREEN, "\n\tSpecify the minimum bit score to a transcript for each read to be considered for the summary using: --min_bit_score", RESET;

print GREEN, "\n\tSupply a file containing gene descriptions using: --gene_description_file", RESET;
print GREEN, "\n\tSpecify the name of the output gene summary file using: --gene_summary_file", RESET;
print GREEN, "\n\tSpecify the target UCSC directory for custom UCSC track files using: --ucsc_dir", RESET;
print GREEN, "\n\tSpecify the html web path to this directory using: --web_path", RESET;
print GREEN, "\n\nExample: summarizeGeneResults.pl  --database=ALEXA_hs_48_36j  --server=jango.bcgsc.ca  --user=malachig  --password=pwd  --library=HS04391  --read_records_infile=/projects/malachig/solexa/read_records/HS04391/ENST_v49/HS04391_Lanes1-8_ENST_v49.txt  --min_bit_score=48.1  --gene_description_file=/projects/malachig/solexa/EnsEMBL_49_GeneDescriptions.txt  --gene_summary_file=GeneSummary.txt  --ucsc_dir=/home/malachig/www/public/htdocs/solexa/HS04391/  --web_path=http://www.bcgsc.ca/people/malachig/htdocs/solexa/HS04391/\n\n", RESET;

unless ($database && $server && $user && $password && $library && $read_records_infile && $min_bit_score && $gene_description_file && $gene_summary_file && $ucsc_dir && $web_path){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Delete the outfile if it is present to ensure a clean creation
if (-e $gene_summary_file){
  my $cmd = "rm $gene_summary_file";
  system ($cmd);
}

#Get gene descriptions for each gene
my %gene_descriptions;
&getGeneDescriptions('-infile'=>$gene_description_file);

#0.) First get all the neccessary gene info required to perform the analysis
#Establish connection with the Alternative Splicing Expression database
my $genes_ref;
my $gene_transcripts_ref;
my $gene_exon_content_ref;
my $gene_exon_coverage_ref;
&getBasicGeneInfo();


#Open the read file and begin generating the summary.
# - Gene level data can be added to the existing gene hash reference
#   - UCSC track lines for each quality read can be stored in an array for each chromosome (hash of arrays keyed on chromosome)
# - Predicted gene fusions will be added to a new hash.  Create a new key (i.e. 'gene1_gene2' where gene1_gene2 are numerically ordered)
#   - UCSC track lines for each read participating in a gene fusion can be stored in an array for each chromosome (hash of arrays keyed on chromosome)
my %ensembl_transcript_tracks;

my %columns;
my @required_columns = qw(Read_ID DistanceBetweenReads_Genomic DistanceBetweenReads_Transcript R1_ID R1_HitType R1_GeneID R1_AlignmentLength R1_Chromosome R1_ChrStartCoords R1_ChrEndCoords R2_ID R2_HitType R2_GeneID R2_AlignmentLength R2_Chromosome R2_ChrStartCoords R2_ChrEndCoords);

my $header_line = 1;
my $header;
my $total_reads_parsed = 0;
my $both_reads_fail_quality = 0;
my $both_reads_pass_quality = 0;
my $current_line = '';
my $chromosomes_processed = 0;

#Get a list of chromosomes and Gene IDs from the input file to allow chromosome-by-chromosome processing
$| = 1; print BLUE, "\n\n3.) Get a list of chromosomes and associated genes from the input file before proceeding\n", RESET; $| = 0;
my $line_counter = 0;
my %chr_list;
open (READ, "$read_records_infile") || die "\nCould not open read records summary infile: $read_records_infile\n\n";
while(<READ>){
  $line_counter++;

  if ($line_counter == 1000000){
    $| = 1; print "."; $| = 0;
    $line_counter = 0;
  }

  chomp($_);
  $current_line = $_;
  my @line = split("\t", $_);

  #Parse the column names and positions.  Check against a hard coded list of required columns before proceeding
  if ($header_line == 1){
    $header = $current_line;
    my $col_count = 0;
    foreach my $column (@line){
      $columns{$column}{position} = $col_count;
      $col_count++;
    }

    foreach my $req_column (@required_columns){
      unless ($columns{$req_column}){
	print RED, "\nRequired column: $req_column was not found in the read record file!\n\n", RESET;
	exit();
      }
    }
    $header_line = 0;
    next();
  }

  my $r1_gene_id = $line[$columns{R1_GeneID}{position}];
  my $r2_gene_id = $line[$columns{R2_GeneID}{position}];
  my $r1_chromosome = $line[$columns{R1_Chromosome}{position}];
  my $r2_chromosome = $line[$columns{R2_Chromosome}{position}];

  #Fix chromosome formats
  if ($r1_chromosome eq "chrMT"){$r1_chromosome = "chrM";}
  if ($r2_chromosome eq "chrMT"){$r2_chromosome = "chrM";}

  #If the read1 gene ID seems valid add it to the list
  if ($r1_gene_id =~ /\d+/){
    if ($chr_list{$r1_chromosome}){
      my $gene_ids_ref = $chr_list{$r1_chromosome}{gene_ids};
      $gene_ids_ref->{$r1_gene_id}->{tmp} = '';
    }else{
      my %gene_ids;
      $gene_ids{$r1_gene_id}{tmp} = '';
      $chr_list{$r1_chromosome}{gene_ids} = \%gene_ids;
    }
  }

  #If the read2 gene ID seems valid add it to the list
  if ($r2_gene_id =~ /\d+/){
    if ($chr_list{$r2_chromosome}){
      my $gene_ids_ref = $chr_list{$r2_chromosome}{gene_ids};
      $gene_ids_ref->{$r2_gene_id}->{tmp} = '';
    }else{
      my %gene_ids;
      $gene_ids{$r2_gene_id}{tmp} = '';
      $chr_list{$r2_chromosome}{gene_ids} = \%gene_ids;
    }
  }
}
close(READ);

my $chrs_found = keys %chr_list;
$| = 1; print BLUE, "\n\tFound $chrs_found chromosomes\n", RESET; $| = 0;

$line_counter = 0;

#Now begin processing data chromosome-by-chromosome
#For each chromosome, only deal with genes for which a least one read was found...
foreach my $current_chr (sort keys %chr_list){

  #Get the list of genes for this chromosome
  my @chr_genes;
  my $gene_list_ref = $chr_list{$current_chr}{gene_ids};

  foreach my $gene_id (keys %{$gene_list_ref}){
    push (@chr_genes, $gene_id);
  }
  my $chr_gene_count = scalar(@chr_genes);
  if ($chr_gene_count == 0){
    next();
  }
  $chromosomes_processed++;

  $| = 1; print YELLOW, "\n\n4.) Processing chromosome: $current_chr (contains $chr_gene_count genes with read hits in the input file)", RESET; $| = 0;

  #Get gene exon coverage for only the current set of genes
  %{$gene_exon_coverage_ref} = ();
  $gene_exon_coverage_ref = &getExonCoverageObject('-gene_ids'=>\@chr_genes);

  $| = 1; print BLUE, "\n4-c.) Parsing read record file and gathering read info for this block of reads\n", RESET; $| = 0;

  open (READ, "$read_records_infile") || die "\nCould not open read records summary infile: $read_records_infile\n\n";
  $header_line = 1;
  while(<READ>){
    $line_counter++;

    if ($line_counter == 1000000){
      $| = 1; print "."; $| = 0;
      $line_counter = 0;
    }

    chomp($_);
    $current_line = $_;
    my @line = split("\t", $_);

    #Parse the column names and positions.  Check against a hard coded list of required columns before proceeding
    if ($header_line == 1){
      $header_line = 0;
      next();
    }
    $total_reads_parsed++;

    my $r1_bit_score = $line[$columns{R1_BitScore}{position}];
    my $r2_bit_score = $line[$columns{R2_BitScore}{position}];
    my $r1_chromosome = $line[$columns{R1_Chromosome}{position}];
    my $r2_chromosome = $line[$columns{R2_Chromosome}{position}];

    #Fix chromosome formats
    if ($r1_chromosome eq "chrMT"){$r1_chromosome = "chrM";}
    if ($r2_chromosome eq "chrMT"){$r2_chromosome = "chrM";}

    #Before proceeding, make sure at least one of these reads corresponds to the current chromosome
    unless (($r1_chromosome eq "$current_chr") || ($r2_chromosome eq "$current_chr")){
      next();
    }
    #change alignment lengths of 'NA' to 0
    if ($r1_bit_score eq "NA"){$r1_bit_score = 0;}
    if ($r2_bit_score eq "NA"){$r2_bit_score = 0;}

    #Check this line to see if either read passes the basic quality filters. (Paired reads will be considered seperately for gene-level summaries)
    #If both reads fail, skip this record immediately
    unless ($r1_bit_score >= $min_bit_score || $r2_bit_score >= $min_bit_score){
      $both_reads_fail_quality++;
      next();
    }

    my $read_id = $line[$columns{Read_ID}{position}];
    my $distance_between_reads_genomic = $line[$columns{DistanceBetweenReads_Genomic}{position}];
    my $distance_between_reads_transcript = $line[$columns{DistanceBetweenReads_Transcript}{position}];
    my $r1_id = $line[$columns{R1_ID}{position}];
    my $r2_id = $line[$columns{R2_ID}{position}];
    my $r1_hit_type = $line[$columns{R1_HitType}{position}];
    my $r2_hit_type = $line[$columns{R2_HitType}{position}];
    my $r1_gene_id = $line[$columns{R1_GeneID}{position}];
    my $r2_gene_id = $line[$columns{R2_GeneID}{position}];
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

    #change distance_between read values of 'NA' to 0
    if ($distance_between_reads_transcript eq "NA"){$distance_between_reads_transcript = 0};

    #Test Read1 and Read2 to see if they pass the quality threshold individually
    my $read1_passes = 0;
    my $read2_passes = 0;

    if (($r1_hit_type eq "Top_Hit") && ($r1_bit_score >= $min_bit_score) && ($r1_chromosome eq "$current_chr")){
      $read1_passes = 1;
    }
    if (($r2_hit_type eq "Top_Hit") && ($r2_bit_score >= $min_bit_score) && ($r2_chromosome eq "$current_chr")){
      $read2_passes = 1;
    }

    #Add a count for the gene that this read hits.
    #Deal with READ1
    if ($read1_passes == 1){

      #Make sure the gene ID found is valid
      if ($genes_ref->{$r1_gene_id}){

	#Count this read hit to the gene
	$gene_exon_coverage_ref->{$r1_gene_id}->{quality_read_count}++;

	#Add the coverage of this read to the gene
	&addReadCoverage('-gene_id'=>$r1_gene_id, '-read_id'=>$r1_id, '-chr_starts'=>\@r1_chr_start_coords, '-chr_ends'=>\@r1_chr_end_coords);

      }else{
	print RED, "\nCould not find gene id: $r1_gene_id in Gene info object!\n\n", RESET;
	exit();
      }
    }

    #Deal with READ2
    if ($read2_passes == 1){

      #Make sure the gene ID found is valid
      if ($genes_ref->{$r2_gene_id}){

	#Count this read hit to the gene
	$gene_exon_coverage_ref->{$r2_gene_id}->{quality_read_count}++;

	#Add the coverage of this read to the gene
	&addReadCoverage('-gene_id'=>$r2_gene_id, '-read_id'=>$r2_id, '-chr_starts'=>\@r2_chr_start_coords, '-chr_ends'=>\@r2_chr_end_coords);

      }else{
	print RED, "\nCould not find gene id: $r2_gene_id in Gene info object!\n\n", RESET;
	exit();
      }
    }

  }
  close (READ);

  #Print out Gene Summary file
  &printGeneSummary('-outfile'=>$gene_summary_file);

  #Print out UCSC track files
  &printUCSCTracks('-ucsc_dir'=>$ucsc_dir, '-chr'=>$current_chr);
}

print BLUE, "\n\nParsed a total of $total_reads_parsed paired reads", RESET;
print BLUE, "\n\tFor $both_reads_fail_quality of these paired reads, BOTH individual reads did NOT meet the quality criteria\n\n", RESET;
print BLUE, "\n\tFor $both_reads_pass_quality of these paired reads, BOTH individual reads did PASS the quality criteria\n\n", RESET;

exit();


############################################################################################################################################
#Get gene descriptions for each gene
############################################################################################################################################
sub getGeneDescriptions{
  my %args = @_;
  my $infile = $args{'-infile'};

  $| = 1; print BLUE, "\n1.) Parsing gene descriptions from: $infile\n", RESET; $| = 0;

  open (DESC, "$infile") || die "\nCould not open gene descriptions file\n\n", RESET;

  my $first_line = 1;
  while(<DESC>){
    chomp($_);
    if ($first_line == 1){
      $first_line = 0;
      next();
    }
    my @line = split("\t", $_);

    my $ensg_id = $line[0];
    my $desc = $line[1];

    $gene_descriptions{$ensg_id}{description} = $desc;
  }

  close(DESC);

  return();
}


############################################################################################################################################
#Get basic info for all genes from the user specified ALEXA database                                                                       #
############################################################################################################################################
sub getBasicGeneInfo{
  my %args = @_;

  #Establish connection with the Alternative Splicing Expression database
  my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

  my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};

  #Get the gene info for all genes for which reads were found on the current chromosome
  $| = 1; print BLUE, "\n2-a.) Getting gene data", RESET; $| = 0;
  $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

  #Get the transcript info for all transcripts of these genes
  $| = 1; print BLUE, "\n2-b.) Getting transcript data as well as exons for each transcript", RESET; $| = 0;
  $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

  #Get chromosome coordinates for all EnsEMBL transcripts and Build UCSC tracks for all ensembl transcripts
  $| = 1; print BLUE, "\n\n2-c.) Calculating chromosome coordinates for the EXONS of each gene", RESET; $| = 0;
  foreach my $gene_id (keys %{$gene_transcripts_ref}){

    my $chromosome = $genes_ref->{$gene_id}->{chromosome};

    if ($chromosome eq "MT"){$chromosome = "M";}

    my $ucsc_chromosome = "chr"."$genes_ref->{$gene_id}->{chromosome}";
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

	#Make sure the supplied coordinates are actually within the specified gene
	unless ($start >= $gene_start-1 && $start <= $gene_end+1){
	  print RED, "\nStart coordinate ($start) provided to convertGeneCoordinates() does not appear valid for gene_id $gene_id\n\n", RESET;
	  exit();
	}
	unless ($end >= $gene_start-1 && $end <= $gene_end+1){
	  print RED, "\nEnd coordinate ($end) provided to convertGeneCoordinates() does not appear valid for gene_id $gene_id\n\n", RESET;
	  exit();
	}

	#Convert provided gene coordinates to coordinates relative to the chromosome
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

	  $exons_ref->{$exon_id}->{chr_start} = $query_chr_start;
	  $exons_ref->{$exon_id}->{chr_end} = $query_chr_end;
	  $exons_ref->{$exon_id}->{strand} = "+";

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

	  $exons_ref->{$exon_id}->{chr_start} = $query_chr_start;
	  $exons_ref->{$exon_id}->{chr_end} = $query_chr_end;
	  $exons_ref->{$exon_id}->{strand} = "-";

	}else{
	  print RED, "\nStrand format: $chr_strand not understood !\n\n", RESET;
	  exit();
	}
      }

      #Now use the exon chromosome coordinates to generate custom UCSC track records for this transcript
      my $ensembl_t_id = $transcripts_ref->{$trans_id}->{ensembl_t_id};
      foreach my $exon_id (sort {$exons_ref->{$a}->{chr_start} <=> $exons_ref->{$b}->{chr_start}} keys %{$exons_ref}){
	my $start = $exons_ref->{$exon_id}->{chr_start};
	my $end = $exons_ref->{$exon_id}->{chr_end};
	my $strand = $exons_ref->{$exon_id}->{strand};

	my $record = "\n$ucsc_chromosome\tEnsEMBL\texon\t$start\t$end\t.\t$strand\t.\t$ensembl_t_id";

	if ($ensembl_transcript_tracks{$ucsc_chromosome}){
	  my $records_ref = $ensembl_transcript_tracks{$ucsc_chromosome}{records};
	  push(@{$records_ref}, $record);
	}else{
	  my @records;
	  push(@records, $record);
	  $ensembl_transcript_tracks{$ucsc_chromosome}{records} = \@records;
	}
      }
    }
  }

  #Get exon content for all genes
  $| = 1; print BLUE, "\n\n2-d.) Getting EXON CONTENT of each gene", RESET; $| = 0;
  my $storable_name = "$database"."_AllGenes_ExonContent.storable";
  $gene_exon_content_ref = &getExonContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-storable'=>$storable_name);

  #Close database connection
  $alexa_dbh->disconnect();

  return();
}

############################################################################################################################################
#Initialize hashes to store Exon coverage for a subset of genes
############################################################################################################################################
sub getExonCoverageObject{
  my %args = @_;
  my @gene_ids = @{$args{'-gene_ids'}};

  my %gene_exon_coverage;

  #At this time, also get the chromosome coordinates for all exon-content coordinates
  $| = 1; print BLUE, "\n\n4-a.) Calculating chromosome coordinates for the EXON CONTENT of each gene (for only the current chromosome)\n", RESET; $| = 0;
  my $counter = 0;
  foreach my $gene_id (@gene_ids){

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

      #Make sure the supplied coordinates are actually within the specified gene
      unless ($start >= $gene_start-1 && $start <= $gene_end+1){
	print RED, "\nStart coordinate ($start) provided to convertGeneCoordinates() does not appear valid for gene_id $gene_id\n\n", RESET;
	exit();
      }
      unless ($end >= $gene_start-1 && $end <= $gene_end+1){
	print RED, "\nEnd coordinate ($end) provided to convertGeneCoordinates() does not appear valid for gene_id $gene_id\n\n", RESET;
	exit();
      }

      #Convert provided gene coordinates to coordinates relative to the chromosome
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
	$exon_content_chr{$exon_id}{chr_start} = $query_chr_start;
	$exon_content_chr{$exon_id}{chr_end} = $query_chr_end;
	$exon_content_chr{$exon_id}{strand} = "+";
	$exon_content_chr{$exon_id}{size} = ($query_chr_end - $query_chr_start)+1;
	$gene_exon_coverage{$gene_id}{exon_content_size} += ($query_chr_end - $query_chr_start)+1;

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

	$exon_content_chr{$exon_id}{chr_start} = $query_chr_start;
	$exon_content_chr{$exon_id}{chr_end} = $query_chr_end;
	$exon_content_chr{$exon_id}{strand} = "-";
	$exon_content_chr{$exon_id}{size} = ($query_chr_end - $query_chr_start)+1;
	$gene_exon_coverage{$gene_id}{exon_content_size} += ($query_chr_end - $query_chr_start)+1;

      }else{
	print RED, "\nStrand format: $chr_strand not understood !\n\n", RESET;
	exit();
      }
    }
    $gene_exon_coverage{$gene_id}{exon_content} = \%exon_content_chr;

  }

  #Build a coverage hash for each gene (keyed on chromosome positions)
  #At this time, also get the chromosome coordinates for all exon-content coordinates
  $| = 1; print BLUE, "\n4-b.) Initializing a COVERAGE hash for the EXON CONTENT of each gene (for only the current chromosome)\n", RESET; $| = 0;

  $counter = 0;
  foreach my $gene_id (sort keys %gene_exon_coverage){

    $counter++;
    if ($counter == 100){
      $counter = 0;
      $| = 1; print BLUE, ".", RESET; $| = 0;
    }

    my %coverage;

    my $exon_content_ref = $gene_exon_coverage{$gene_id}{exon_content};

    foreach my $exon_id (sort {$exon_content_ref->{$a}->{chr_start} <=> $exon_content_ref->{$b}->{chr_start}} keys %{$exon_content_ref}){

      my $start = $exon_content_ref->{$exon_id}->{chr_start};
      my $end = $exon_content_ref->{$exon_id}->{chr_end};

      #Go through each chromosome position in this exon content block and initialize that position in the hash
      for (my $i = $start; $i <= $end; $i++){
	$coverage{$i}{cov} = 0;
      }
    }
    $gene_exon_coverage{$gene_id}{coverage} = \%coverage;

    #Initialize quality read_count value
    $gene_exon_coverage{$gene_id}{quality_read_count} = 0;
  }

  return(\%gene_exon_coverage);
}


############################################################################################################################################
#Add the coverage of a read to a gene to the coverage hash for the exon content record for that gene
############################################################################################################################################
sub addReadCoverage{
  my %args = @_;
  my $gene_id = $args{'-gene_id'};
  my $read_id = $args{'-read_id'};
  my @chr_starts = @{$args{'-chr_starts'}};
  my @chr_ends = @{$args{'-chr_ends'}};

  my $coverage_ref = $gene_exon_coverage_ref->{$gene_id}->{coverage};

  my @ends = @chr_ends;
  foreach my $start (@chr_starts){
    my $end = shift(@ends);

    #Go through each chromosome position in this read as it is mapped to an exon and increment that position in the hash
    for (my $i = $start; $i <= $end; $i++){

      if ($coverage_ref->{$i}){
	$coverage_ref->{$i}->{cov}++;
      }else{
	print RED, "\nRead has coverage coordinates which do not correspond to exon content!", RESET;
	print RED, "\nStarts: @chr_starts\tEnds: @chr_ends", RESET;
	print Dumper $coverage_ref;
	exit();
      }
    }
  }

  return();
}


############################################################################################################################################
#Print out Gene Summary file
############################################################################################################################################
sub printGeneSummary{
  my %args = @_;
  my $outfile = $args{'-outfile'};

  $| = 1; print BLUE, "\n5.) Printing Gene Summary data to: $outfile", RESET; $| = 0;

  #Open the output file
  open (GENE_OUT, ">>$outfile") || die "\nCould not open output gene summary file: $outfile\n\n";

  if ($chromosomes_processed == 1){
    print GENE_OUT "ALEXA ID\tEnsEMBL Gene ID\tGene Name\tRead Count\tExonic Base Count\tPercent Coverage (>= 1X)\tPercent Coverage (>= 2X)\tPercent Coverage (>= 5X)\tPercent Coverage (>= 10X)\tAverage Coverage (X)\tTotal Sequence Bases\tTotal Bases Covered (>= 1X)\tGene Type\tGene Evidence\tChromosome\tStrand\tStart\tEnd\tDescription\tLink\n";
  }

  unless ($web_path =~ /.*\/$/){
    $web_path = "$web_path"."/";
  }

  #The coverage of each base of the exon content of each gene has been stored for each gene
  #Use this information to calculate the sequence coverage of each gene
  #Summarize the % of bases with at least 1X coverage and the overall average X coverage
  $| = 1; print BLUE, "\n5-a.) Summarizing sequence coverage of each gene ...\n", RESET; $| = 0;
  my $counter = 0;
  foreach my $gene_id (keys %{$gene_exon_coverage_ref}){
    $counter++;
    if ($counter == 100){
      $counter = 0;
      $| = 1; print BLUE, ".", RESET; $| = 0;
    }

    my $coverage_ref = $gene_exon_coverage_ref->{$gene_id}->{coverage};

    my $total_bases = keys %{$coverage_ref};
    my $cumulative_coverage = 0;
    my $bases_covered_1x = 0;
    my $bases_covered_2x = 0;
    my $bases_covered_5x = 0;
    my $bases_covered_10x = 0;

    foreach my $pos (sort {$a <=> $b} keys %{$coverage_ref}){
      if ($coverage_ref->{$pos}->{cov} >= 1){
	$bases_covered_1x++;
      }
      if ($coverage_ref->{$pos}->{cov} >= 2){
	$bases_covered_2x++;
      }
      if ($coverage_ref->{$pos}->{cov} >= 5){
	$bases_covered_5x++;
      }
      if ($coverage_ref->{$pos}->{cov} >= 10){
	$bases_covered_10x++;
      }
      $cumulative_coverage += $coverage_ref->{$pos}->{cov};
    }
    my $percent_coverage_1x = sprintf("%.2f", (($bases_covered_1x/$total_bases)*100));
    my $percent_coverage_2x = sprintf("%.2f", (($bases_covered_2x/$total_bases)*100));
    my $percent_coverage_5x = sprintf("%.2f", (($bases_covered_5x/$total_bases)*100));
    my $percent_coverage_10x = sprintf("%.2f", (($bases_covered_10x/$total_bases)*100));
    my $average_coverage = sprintf("%.2f", ($cumulative_coverage/$total_bases));

    $gene_exon_coverage_ref->{$gene_id}->{percent_coverage_1x} = $percent_coverage_1x;
    $gene_exon_coverage_ref->{$gene_id}->{percent_coverage_2x} = $percent_coverage_2x;
    $gene_exon_coverage_ref->{$gene_id}->{percent_coverage_5x} = $percent_coverage_5x;
    $gene_exon_coverage_ref->{$gene_id}->{percent_coverage_10x} = $percent_coverage_10x;
    $gene_exon_coverage_ref->{$gene_id}->{average_coverage} = $average_coverage;
    $gene_exon_coverage_ref->{$gene_id}->{total_sequence_bases} = $cumulative_coverage; #Total bases of sequence mapping to this gene
    $gene_exon_coverage_ref->{$gene_id}->{bases_covered_1x} = $bases_covered_1x; #Total non-redundant bases sequenced to 1x or greater depth

  }

  my $current_read_count = 1000000000000000000000000000000000000000000000000;

  #Go through each gene and print out basic data for it.  Also include a link to the custom track file for its chromosome
  $| = 1; print BLUE, "\n5-b.) Printing actual summary lines ...\n", RESET; $| = 0;
  foreach my $gene_id (sort {$gene_exon_coverage_ref->{$b}->{quality_read_count} <=> $gene_exon_coverage_ref->{$a}->{quality_read_count}} keys %{$gene_exon_coverage_ref}){

    my $description;
    my $ensg_id = $genes_ref->{$gene_id}->{ensembl_g_id};
    if ($gene_descriptions{$ensg_id}){
      $description = $gene_descriptions{$ensg_id}{description};
    }else{
      print RED, "\nCould not find a gene description for $ensg_id !! - Check gene descriptions file - ensure ensembl version match\n\n", RESET;
      exit();
    }

    #Create a link to go directly to the region of this gene (+/- 100bp) and load the correct chromosome file
    my $display_start = $genes_ref->{$gene_id}->{chr_start};
    my $display_end = $genes_ref->{$gene_id}->{chr_end};
    my $temp;
    if ($display_start > $display_end){
      $temp = $display_start;
      $display_start = $display_end;
      $display_end = $temp;
    }
    $display_start -= 100;
    $display_end += 100;

    my $chromosome = "chr"."$genes_ref->{$gene_id}->{chromosome}";

    my $link = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=$chromosome:$display_start-$display_end&hgt.customText=$web_path"."$chromosome".".txt.gz"."&ctfile_hg18=";

    print GENE_OUT "$gene_id\t$ensg_id\t$genes_ref->{$gene_id}->{gene_name}\t$gene_exon_coverage_ref->{$gene_id}->{quality_read_count}\t$gene_exon_coverage_ref->{$gene_id}->{exon_content_size}\t$gene_exon_coverage_ref->{$gene_id}->{percent_coverage_1x}\t$gene_exon_coverage_ref->{$gene_id}->{percent_coverage_2x}\t$gene_exon_coverage_ref->{$gene_id}->{percent_coverage_5x}\t$gene_exon_coverage_ref->{$gene_id}->{percent_coverage_10x}\t$gene_exon_coverage_ref->{$gene_id}->{average_coverage}\t$gene_exon_coverage_ref->{$gene_id}->{total_sequence_bases}\t$gene_exon_coverage_ref->{$gene_id}->{bases_covered_1x}\t$genes_ref->{$gene_id}->{gene_type}\t$genes_ref->{$gene_id}->{evidence}\t$genes_ref->{$gene_id}->{chromosome}\t$genes_ref->{$gene_id}->{chr_strand}\t$genes_ref->{$gene_id}->{chr_start}\t$genes_ref->{$gene_id}->{chr_end}\t$description\t$link\n";

  }

  close(GENE_OUT);

  return();
}


############################################################################################################################################
#Print out UCSC track files
############################################################################################################################################
sub printUCSCTracks{
  my %args = @_;
  my $ucsc_dir = $args{'-ucsc_dir'};
  my $chromosome = $args{'-chr'};

  unless ($ucsc_dir =~ /.*\/$/){
    $ucsc_dir = "$ucsc_dir"."/";
  }

  my $ucsc_file = "$ucsc_dir"."$chromosome".".txt";
  #Print out all the UCSC track records gathered above for this chromosome
  open (UCSC, ">$ucsc_file") || die "\nCould not open ucsc file: $ucsc_file\n\n";

  $| = 1; print BLUE, "\n6.) Printing UCSC file for $chromosome: $ucsc_file", RESET; $| = 0;

  #Browser line
  print UCSC "#Browser line";
  print UCSC "\nbrowser hide all";
  print UCSC "\nbrowser full knownGene";
  print UCSC "\nbrowser pack multiz28way";


  #1.) Track line for EnsEMBL transcripts of each gene
  #$ensembl_transcript_tracks{$ucsc_chromosome}{records} = \@records;
  print UCSC "\n\n#EnsEMBL transcripts";
  my $trans_track_name = "$database";
  my $trans_track_description = "\"EnsEMBL Transcripts used for Mapping ($database)\"";

  #Brick RED: 153,0,0
  print UCSC "\ntrack name=$trans_track_name description=$trans_track_description color=153,0,0 useScore=0 visibility=2 priority=0";
  print UCSC "\n\n#Begin DATA";

  my $trans_records_ref = $ensembl_transcript_tracks{$chromosome}{records};
  foreach my $trans_record (@{$trans_records_ref}){
    print UCSC "$trans_record";
  }

  #2.) Create a WIG track line to display the coverage level for every exonic base sequence to 1X or greater depth
  #NOTE: Currently the coverage is stored at the gene level.
  #It is possible that there is some overlap between genes that will result in overlapping coverage at chromosome level
  #Use the gene level coverage for all genes of each chromosome to create a chromosome coverage hash
  #Then use this hash to generate a WIG track for the entire chromosome
  $| = 1; print BLUE, "\n6-a.) Creating coverage hash at entire chromosome level", RESET; $| = 0;
  my %chr_coverage;
  foreach my $gene_id (keys %{$gene_exon_coverage_ref}){
    my $coverage_ref = $gene_exon_coverage_ref->{$gene_id}->{coverage};

    #For the chromosome coverage, dont worry about having 0 counts, just mark the position with some coverage
    foreach my $pos (sort {$a <=> $b} keys %{$coverage_ref}){

      if ($coverage_ref->{$pos}->{cov} > 0){
	
	#If this position is already defined in the chromosome (indicates overlapping regions between genes) then add the coverage together
	if ($chr_coverage{$pos}){
	  #If this position is already defined in the chromosome (indicates overlapping regions between genes) then add the coverage together
	  $chr_coverage{$pos}{cov} += $coverage_ref->{$pos}->{cov};
	}else{
	  #Otherwise simply initialize this position with the observed coverage level
	  $chr_coverage{$pos}{cov} = $coverage_ref->{$pos}->{cov};
	}
      }
    }
  }

  #Use the result coverage hash to generate a wig file
  print UCSC "\n\n#WIG TRACK: Exonic base coverage calculated from Illumina Paired Reads Mapped to EnsEMBL Transcripts";

  my $wig_track_name = "$database"."_Wig_"."$library";
  my $wig_track_description = "\"Exonic base coverage from Illumina Paired Reads ($library) Mapped to EnsEMBL Transcripts - UnAmbiguous Hits >= $min_bit_score bit score\"";
  #Dark Blue:  51,0,153

  print UCSC "\ntrack name=$wig_track_name description=$wig_track_description type=wiggle_0 color=51,0,153 yLineMark=0.0 yLineOnOff=on visibility=full autoScale=on graphType=bar smoothingWindow=off maxHeightPixels=80:40:2 priority=1";

  my $current_pos = -1;

  foreach my $pos (sort {$a <=> $b} keys %chr_coverage){

    #If this is a new block of covered bases, print out a new def line
    if ($pos == ($current_pos+1)){
      print UCSC "\n$chr_coverage{$pos}{cov}";
    }else{
      print UCSC "\nfixedStep chrom=$chromosome start=$pos step=1";
      print UCSC "\n$chr_coverage{$pos}{cov}";
    }
    $current_pos = $pos;
  }
  close(UCSC);

  #Finally gzip the file
  my $cmd = "gzip -f $ucsc_file";
  system($cmd);

  return();
}

