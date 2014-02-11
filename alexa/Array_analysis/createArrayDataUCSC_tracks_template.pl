#!/usr/bin/perl -w
#Written by Malachi Griffith

#The purpose of this script is to create custom UCSC tracks for blocks genes (ALEXA genes only to start)
#This script requires as input a file containing ALEXA gene IDs
#- The first three columns must be: AlexaGene_ID    Probe_ID        ProbeSet_ID
#- After this the file may have any number of columns but only those with names defined at the top of this script will be imported


#Seperate files linking these probeset IDs to PLIER or DABG expression value or splicing index values etc. will be required
#Similarly a corresponding file with gene-level expression values generated from a metaprobeset will be considered

#Sample files:
#meta-transcript-level
#/home/malachig/AlternativeSplicing/Array_analysis/Affymetrix_Exon_Arrays/MIPvs5FUR_6-20-2006/genelevel_ensembl-quant-antibgp/mip_vs_5fur/import/normalize/summarize/plier.summary.clean.txt

#probeset-level
#/home/malachig/AlternativeSplicing/Array_analysis/Affymetrix_Exon_Arrays/MIPvs5FUR_6-20-2006/probelevel-quant-antibgp/import/normalize/summarize/plier.summary.txt

use strict;
use Data::Dumper;
use Getopt::Long;
use Math::Complex;
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
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $genome_build = '';
my $map_file = '';
my $probeset_expression_file = '';
my $metaprobeset_expression_file = '';
my $ucsc_dir = '';
my $summary_file = '';
my $gzip_path = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'genome_build=s'=>\$genome_build,
	    'map_file=s'=>\$map_file, 'probeset_expression_file=s'=>\$probeset_expression_file, 'metaprobeset_expression_file=s'=>\$metaprobeset_expression_file,
	    'ucsc_dir=s'=>\$ucsc_dir, 'summary_file=s'=>\$summary_file, 'gzip_path=s'=>\$gzip_path);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the correct UCSC genome build ID using: --genome_build (e.g. 'hg18')", RESET;
print GREEN, "\n\tSpecify the input mapfile describing the maping of probesets to each Ensembl gene using --map_file", RESET;
print GREEN, "\n\tSpecify the input probeset expression data file using: --probeset_expression_file", RESET;
print GREEN, "\n\tSpecify the input meta-probeset expression data file using: --metaprobeset_expression_file", RESET;
print GREEN, "\n\tSpecify the full path to the target directory for UCSC track files using: --ucsc_dir", RESET;
print GREEN, "\n\tSpecify the output summary file name using: --summary_file", RESET;
print GREEN, "\n\tSpecify the full path to your gzip utility using: --gzip_path", RESET;
print GREEN, "\n\nExample: createArrayDataUCSC_tracks.pl  --database=ALEXA_hs_45_36g  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --genome_build=hg18  --map_file=/home/malachig/Collab/Lymphoma/mapfiles/hg18/CORE_probesets_mapped_to_ensembl45.txt  --probeset_expression_file=/home/malachig/Collab/Lymphoma/example_data_files/plier.summary.log2_exon-level.txt  --metaprobeset_expression_file=/home/malachig/Collab/Lymphoma/example_data_files/plier.summary.log2_gene-level.txt  --ucsc_dir=/home/malachig/www/public/htdocs/test   --summary_file=gene_ucsc_summary_file.txt  --gzip_path=/usr/bin/gzip\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $genome_build && $map_file && $probeset_expression_file && $metaprobeset_expression_file && $ucsc_dir && $summary_file && $gzip_path){
  print RED, "\nOptions missing!\n\n", RESET;
  exit();
}

unless ($genome_build =~ /\w{2}\d+/){
  print RED, "\nFormat of genome build you supplied ($genome_build) was not understood\n\n", RESET;
  exit();
}

my @colors_heat = ("255,242,0","255,230,0","255,217,0","255,205,0","255,192,0","255,180,0","255,167,0","255,155,0","255,142,0","255,130,0","255,117,0","255,105,0","255,92,0","255,80,0","255,67,0","255,55,0","255,42,0","255,30,0","255,12,0","255,0,0");
my @colors_green_to_blue = ("0,255,30","0,255,55","0,255,80","0,255,105","0,255,130","0,255,155","0,255,180","0,255,205","0,255,230","0,255,255","0,230,255","0,205,255","0,180,255","0,155,255","0,130,255","0,105,255","0,80,255","0,55,255","0,30,255","0,5,255");
my @colors_white_to_purple = ("200,140,255","196,133,255","193,126,255","189,119,255","186,112,255","182,105,255","179,98,255","175,91,255","172,84,255","168,77,255","165,70,255","161,63,255","158,56,255","154,49,255","151,42,255","147,35,255","144,28,255","140,21,255","137,14,255","133,7,255");
my @colors_purple_to_white = reverse(@colors_white_to_purple);

#Create a hardcoded hash of data column names mapped to human readable names
my %data_names;

#Also hardcode the default visibility (0 - hide; 1 - dense; 2 - full; 3 - pack; 4 - squish)
#Also hardcode the display order for each track (small numbers first, or at the top of the screen ...)

#hash{data_column}{required_value} = "";
#Create one of these blocks for each column of data (patient, etc.)
$data_names{'ln_ar_1_log2'}{name} = "ln_ar-1";
$data_names{'ln_ar_1_log2'}{description} = "LnCAP-AR (Androgen Treated) log2 plier data - replicate 1";
$data_names{'ln_ar_1_log2'}{colors} = "@colors_heat";
$data_names{'ln_ar_1_log2'}{visibility} = 1;
$data_names{'ln_ar_1_log2'}{display_order} = 1;

$data_names{'ln_ar_2_log2'}{name} = "ln_ar-2";
$data_names{'ln_ar_2_log2'}{description} = "LnCAP-AR (Androgen Treated) log2 plier data - replicate 2";
$data_names{'ln_ar_2_log2'}{colors} = "@colors_heat";
$data_names{'ln_ar_2_log2'}{visibility} = 1;
$data_names{'ln_ar_2_log2'}{display_order} = 1;


#Make sure the specified UCSC dir is okay before proceeding
unless ($ucsc_dir =~ /.*\/$/){
  $ucsc_dir = "$ucsc_dir"."/";
}
unless (-e $ucsc_dir && -d $ucsc_dir){
  print RED, "\nSpecified UCSC target directory: $ucsc_dir does not appear valid!\n\n", RESET;
  exit();
}
#Empty the target directory
my $cmd = "rm -f $ucsc_dir"."*.txt.gz";
print BLUE, "\n\n\tCMD: $cmd\n\tClean directory with command above (y/n)? ", RESET;
my $answer = <>;
chomp($answer);
if ($answer eq "y" || $answer eq "Y" || $answer eq "yes" || $answer eq "Yes"){
  system($cmd);
}


#Define web path (as it is accessed from the the internet!!)
#/home/malachig/www/public/htdocs/alexa_tracks/mip_v_5fur/Raw/probe_level
my $web_base = "http://www.bcgsc.ca/people/malachig/htdocs/";
my $web_loc;
if ($ucsc_dir =~ /^\/home\/malachig\/www\/public\/htdocs\/(.*)\//){
  $web_loc = $1;
}else{
  print RED, "\ntarget UCSC dir: $ucsc_dir not understood!\n\n", RESET;
  exit();
}


#1.) First get a mapping of all probes/probesets and their chromosome coordinates from the mapping file specified by the user
my %gene_info;
my %gene_probe_map;
my %master_probeset_list;
my %master_gene_list;
my $max_probeset_id = 0;  #Used to create new dummy (legend) probes

#Additional gene info for genes found in this file will also be gathered at this time
#    - Gene symbol, strand, chromosome, start pos, end pos.

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

&parseMapFile('-map_file'=>$map_file);
my $chr_stored = keys %gene_probe_map;
print BLUE, "\nFound and stored probes for: $chr_stored chromosomes\n\n", RESET;

#Close database connection
$alexa_dbh->disconnect();

#print Dumper %gene_probe_map;


#3-A.) Now get the probeset expression data to be displayed
my %probeset_data;
my @probeset_data_columns;
&parseProbesetExpressionData('-probeset_expression_file'=>$probeset_expression_file);
my $probeset_data_stored = keys %probeset_data;
print BLUE, "\nFound and stored probeset expression data for: $probeset_data_stored probesets\n\n", RESET;

#3-B.) Now get the gene-level expression data to be displayed
my %gene_data;
my @gene_data_columns;
&parseGeneExpressionData('-metaprobeset_expression_file'=>$metaprobeset_expression_file);


#4.) Now generate the custom UCSC tracks
#Establish connection with the Alternative Splicing Expression database
$alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
&generateUcscTracks('-ucsc_dir'=>$ucsc_dir);
$alexa_dbh->disconnect();

#5.) Now generate a summary file containing basic info for all genes and a link to the custom UCSC track)
&printSummaryFile('-summary_file'=>$summary_file);

#Go to the target dir and compress the track files
my $command = "$gzip_path". " $ucsc_dir"."*";
print BLUE, "\nCompressing custom track files with command: $command\n\n", RESET;
system ("$command");

exit();




#################################################################################################################################
#1.) First get a mapping of all probesets to Ensembl Genes from the mapping file specified by the user                          #
#################################################################################################################################
sub parseMapFile{
  my %args = @_;
  my $map_file = $args{'-map_file'};

  my $record_count = 0;
  my $gene_count = 0;

  print BLUE, "\nBegin importing Gene-to-Probeset mappings from: $map_file\n\n", RESET;

  open (MAPFILE, "$map_file") || die "\nCould not open mapfile: $map_file\n\n";

  while(<MAPFILE>){
    my $line = $_;
    chomp($line);

    my @line = split ("\t", $line);

    if ($line[0] =~/^\d+/){
      $record_count++;

      my $probeset_id = $line[0];
      $master_probeset_list{$probeset_id}{tmp} = '';

      if ($probeset_id > $max_probeset_id){
	$max_probeset_id = $probeset_id;
      }

      my $alexa_gene_id = $line[1];
      my $ensembl_gene_id = $line[2];
      my $chromosome = $line[4];
      my $strand = $line[5];
      my $start = $line[6];
      my $end = $line[7];
      my $probe_count = $line[8];

      if ($master_gene_list{$alexa_gene_id}){
	$master_gene_list{$alexa_gene_id}{probe_count} += $probe_count;
      }else{
	$master_gene_list{$alexa_gene_id}{probe_count} = $probe_count;
      }

      if ($gene_probe_map{$chromosome}){
	#Grab a copy of the hash reference stored in the gene_probe_map hash
	#Use this reference to directly add to the gene list of each chromosome
	my $gene_list_ref = $gene_probe_map{$chromosome}{gene_list};

	if ($gene_list_ref->{$alexa_gene_id}){

	  #Grab a copy of the hash reference stored in the gene_list hash
	  #Use this reference to directly add to the probesets hash for each gene
	  my $probesets_ref = $gene_list_ref->{$alexa_gene_id}->{probesets};
	  $probesets_ref->{$probeset_id}->{alexa_gene_id} = $alexa_gene_id;
	  $probesets_ref->{$probeset_id}->{ensembl_gene_id} = $ensembl_gene_id;
	  $probesets_ref->{$probeset_id}->{chromosome} = $chromosome;
	  $probesets_ref->{$probeset_id}->{strand} = $strand;
	  $probesets_ref->{$probeset_id}->{start} = $start;
	  $probesets_ref->{$probeset_id}->{end} = $end;
	  $probesets_ref->{$probeset_id}->{probe_count} = $probe_count;

	  #See if the start coordinate for this probeset is smaller than the one currently stored for this gene
	  if ($start < $gene_list_ref->{$alexa_gene_id}->{start_coord}){
	    $gene_list_ref->{$alexa_gene_id}->{start_coord} = $start;
	  }

	}else{
	  #Chromosome has been seen before but not this gene
	  $gene_count++;
	  my %probesets;
	  $probesets{$probeset_id}{alexa_gene_id} = $alexa_gene_id;
	  $probesets{$probeset_id}{ensembl_gene_id} = $ensembl_gene_id;
	  $probesets{$probeset_id}{chromosome} = $chromosome;
	  $probesets{$probeset_id}{strand} = $strand;
	  $probesets{$probeset_id}{start} = $start;
	  $probesets{$probeset_id}{end} = $end;
	  $probesets{$probeset_id}{probe_count} = $probe_count;
	  $gene_list_ref->{$alexa_gene_id}->{start_coord} = $start;
	  $gene_list_ref->{$alexa_gene_id}->{probesets} = \%probesets;
	}

      }else{
	#First time this chromosome (and gene) has been observed
	$gene_count++;
	my %gene_list;
	my %probesets;
	$probesets{$probeset_id}{alexa_gene_id} = $alexa_gene_id;
	$probesets{$probeset_id}{ensembl_gene_id} = $ensembl_gene_id;
	$probesets{$probeset_id}{chromosome} = $chromosome;
	$probesets{$probeset_id}{strand} = $strand;
	$probesets{$probeset_id}{start} = $start;
	$probesets{$probeset_id}{end} = $end;
	$probesets{$probeset_id}{probe_count} = $probe_count;
	$gene_list{$alexa_gene_id}{probesets} = \%probesets;
	$gene_list{$alexa_gene_id}{start_coord} = $start;  #Used to sort genes by the smallest coordinate of their probes
	$gene_probe_map{$chromosome}{gene_list} = \%gene_list;
      }
    }
  }
  close (MAPFILE);
  print BLUE, "\nProcessed: $record_count probeset records corresponding to $gene_count genes\n\n", RESET;

  #Now get neccessary additional info from the ALEXA database for these genes
  my @gene_ids = keys %master_gene_list;
  &getBasicGeneInfo('-dbh'=>$alexa_dbh, '-genes'=>\@gene_ids);

  return();
}



################################################################################################################################
#2.) Get basic gene info from ALEXA on each gene:                                                                              #
################################################################################################################################
sub getBasicGeneInfo{
  my %args = @_;
  my $dbh = $args{'-dbh'};
  my @gene_ids = @{$args{'-genes'}};

  my $gene_count = @gene_ids;

  print BLUE, "\nGetting basic gene info from the ALEXA database for $gene_count genes.", RESET;

  %gene_info = %{&getGeneInfo ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no")};
  my $gene_info_ref = \%gene_info;

  print BLUE, "\nGetting exon content info for each gene", RESET;

  #Get the exon content for the entire list of gene_ids at once
  my $gene_exon_content_ref = &getExonContent ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids);

  print BLUE, "\nConverting gene coordinates to chromosome coordinates for exon content", RESET;
  foreach my $gene_id (sort {$a <=> $b} keys %{$gene_info_ref}){

    #change format of strand info
    if ($gene_info_ref->{$gene_id}->{chr_strand} eq "-1"){
      $gene_info_ref->{$gene_id}->{chr_strand} = "-";
    }
    if ($gene_info_ref->{$gene_id}->{chr_strand} eq "1"){
      $gene_info_ref->{$gene_id}->{chr_strand} = "+";
    }

    #Get the exon content of each gene (for display of gene level expression only)
    my $exonContent_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};

    #Get genomic coordinates for these exon regions
    foreach my $exon_count (keys %{$exonContent_ref}){

      my $start = $exonContent_ref->{$exon_count}->{start};
      my $end = $exonContent_ref->{$exon_count}->{end};
      my $gene_start = $gene_info_ref->{$gene_id}->{gene_start};
      my $gene_end = $gene_info_ref->{$gene_id}->{gene_end};
      my $chr_start = $gene_info_ref->{$gene_id}->{chr_start};
      my $chr_end = $gene_info_ref->{$gene_id}->{chr_end};
      my $chr_strand = $gene_info_ref->{$gene_id}->{chr_strand};

      #Make sure the supplied coordinates are actually within the specified gene
      unless ($start >= $gene_start && $start <= $gene_end){
	print RED, "\nStart coordinate ($start) does not appear valid for gene_id $gene_id\n\n", RESET;
	exit();
      }
      unless ($end >= $gene_start && $end <= $gene_end){
	print RED, "\nEnd coordinate ($end) does not appear valid for gene_id $gene_id\n\n", RESET;
	exit();
      }

      #Convert provided gene coordinates to coordinates relative to the chromosome
      if ($chr_strand eq "+"){
	my $query_chr_start = $chr_start + $start - 1;
	my $query_chr_end = $chr_start + $end - 1;
	$exonContent_ref->{$exon_count}->{chr_start} = $query_chr_start;
	$exonContent_ref->{$exon_count}->{chr_end} = $query_chr_end;

      }elsif ($chr_strand eq "-"){
	my $query_chr_start = $chr_end - $end + 1;
	my $query_chr_end = $chr_end - $start + 1;
	$exonContent_ref->{$exon_count}->{chr_start} = $query_chr_start;
	$exonContent_ref->{$exon_count}->{chr_end} = $query_chr_end;
      }else{
	print RED, "\nStrand format: $chr_strand not understood by convertGeneCoordinates()!\n\n", RESET;
	exit();
      }
    }

    #Store the exon content data on the global gene info object
    $gene_info_ref->{$gene_id}->{exon_content} = $exonContent_ref;
  }

  #Check the number of gene entries successfully stored in the global gene-info hash
  my $gene_info_found = keys %gene_info;

  print BLUE, "\nFound basic gene info from the ALEXA database for $gene_info_found genes.\n\n", RESET;

  return();
}


#################################################################################################################################
#3-A.) Now get the probeset expression data to be displayed                                                                     #
#################################################################################################################################
sub parseProbesetExpressionData{
  my %args = @_;
  my $probeset_expression_file = $args{'-probeset_expression_file'};

  my $record_count = 0;

  print BLUE, "\nBegin importing Probeset expression data from: $probeset_expression_file\n\n", RESET;

  open (PROBESETFILE, "$probeset_expression_file") || die "\nCould not open probeset expression file: $probeset_expression_file\n\n";

  my $header_line;
  my $header = 1;
  my @headers;
  my $column_count = 0;
  my %columns;

  while(<PROBESETFILE>){
    my $line = $_;
    chomp($line);

    #Get the header line
    if ($header == 1){
      $header = 0;
      $header_line = $line;
      @headers = split ("\t", $line);

      foreach my $header (@headers){
	$columns{$header}{column_pos} = $column_count;
	$column_count++;
      }

      #Check each of the probeset_data columns to make sure they are defined in the %data_names hash
      foreach my $data_column (sort{$columns{$a}->{column_pos} <=> $columns{$b}->{column_pos}} keys %columns){
	if ($data_names{$data_column}){
	  push(@probeset_data_columns, $data_column);
	}
      }
      next();
    }

    my @line = split ("\t", $line);

    if ($line[0] =~/^\d+/){
      $record_count++;

      #get the primary id (either a single probe id or a probeset id depending on the file provided)
      my $probeset_id = $line[$columns{probeset_id}{column_pos}];

      #Only import data for those that are in the mapfile
      unless ($master_probeset_list{$probeset_id}){
	next();
      }

      foreach my $data_column (@probeset_data_columns){
	$probeset_data{$probeset_id}{$data_column} = $line[$columns{$data_column}{column_pos}];

      }
    }
  }
  close (PROBESETFILE);
  print BLUE, "\nProcessed: $record_count probeset expression records", RESET;
  print BLUE, "\nFound the following datasets: @probeset_data_columns\n\n", RESET;
  return();
}


#################################################################################################################################
#3-B.) Now get the gene-level expression data to be displayed                                                                   #
#################################################################################################################################
sub parseGeneExpressionData{
  my %args = @_;
  my $gene_expression_file = $args{'-metaprobeset_expression_file'};

  my $record_count = 0;

  print BLUE, "\nBegin importing Gene expression data from: $gene_expression_file\n\n", RESET;

  open (GENEFILE, "$gene_expression_file") || die "\nCould not open gene expression file: $gene_expression_file\n\n";

  my $header_line;
  my $header = 1;
  my @headers;
  my $column_count = 0;
  my %columns;

  while(<GENEFILE>){
    my $line = $_;
    chomp($line);

    #Get the header line
    if ($header == 1){
      $header = 0;
      $header_line = $line;
      @headers = split ("\t", $line);

      foreach my $header (@headers){
	$columns{$header}{column_pos} = $column_count;
	$column_count++;
      }

      #Check each of the gene_data columns to make sure they are defined in the %data_names hash
      foreach my $data_column (sort{$columns{$a}->{column_pos} <=> $columns{$b}->{column_pos}} keys %columns){
	if ($data_names{$data_column}){
	  push(@gene_data_columns, $data_column);
	}
      }
      next();
    }

    my @line = split ("\t", $line);

    if ($line[0] =~/^\d+/){
      $record_count++;

      my $gene_id = $line[$columns{probeset_id}{column_pos}];

      foreach my $data_column (@gene_data_columns){
	$gene_data{$gene_id}{$data_column} = $line[$columns{$data_column}{column_pos}];
      }
    }
  }
  close (GENEFILE);
  print BLUE, "\nProcessed: $record_count gene expression records", RESET;
  print BLUE, "\nFound the following datasets: @gene_data_columns\n\n", RESET;
  return();
}


#################################################################################################################################
#4.) Now generate the custom UCSC tracks                                                                                        #
#################################################################################################################################
sub generateUcscTracks{
  my %args = @_;
  my $ucsc_dir = $args{'-ucsc_dir'};

  print BLUE, "\nBegin creating custom UCSC tracks in: $ucsc_dir", RESET;
  print BLUE, "\nProcessing the following probeset data columns: @probeset_data_columns\n\n", RESET;

  #1.) First step determine the percentile ranges of the datasets ... (0th, 10th, 20th, ..., 100th percentile)
  #    Build arrays that store all the data for each data column

  my %percentiles;
  foreach my $probeset_id (keys %probeset_data){
    foreach my $probeset_data_column (@probeset_data_columns){

      #Avoid 'na' values
      if ($probeset_data{$probeset_id}{$probeset_data_column} eq "na" || $probeset_data{$probeset_id}{$probeset_data_column} eq "-na"){
	next();
      }

      if ($percentiles{$probeset_data_column}){
	push(@{$percentiles{$probeset_data_column}{data_list}}, $probeset_data{$probeset_id}{$probeset_data_column});
      }else{
	my @data_list;
	push (@data_list, $probeset_data{$probeset_id}{$probeset_data_column});
	$percentiles{$probeset_data_column}{data_list} = \@data_list;
      }
    }
  }

  #2.) Now calculate the percentiles for all of the data from each column
  #-   Define the track colors to be used - 20 colors, for 20 bins
  #-   NOTE: decided to use simple bins instead of percentiles

  my %dummy_values;

  my $number_bins = 20;

  foreach my $probeset_data_column (@probeset_data_columns){
    my @data = @{$percentiles{$probeset_data_column}{data_list}};
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@data);

    my $sample_range = $stat->sample_range();
    my $min_value = $stat->min();
    my $current_min = $min_value;

    my @colors = split (" ", $data_names{$probeset_data_column}{colors});

    #Try using linear division of sample range instead of percentiles
    print BLUE, "\n\nDividing sample range into evenly spaced bins for $probeset_data_column, Min= $min_value Range= $sample_range", RESET;
    my %percents;
    my %dummys;

    my $bin_range = $sample_range/$number_bins;
    for (my $i = 1; $i <= $number_bins; $i++){
      $current_min += $bin_range;
      $percents{$i}{percentile} = $current_min + 0.000001;
      $percents{$i}{color} = $colors[$i-1];

      my $dummy_value = ($current_min - ($bin_range / 2));
      my $dummy_value_rounded = sprintf("%.1f", $dummy_value);
      $percents{$i}{dummy_value} = $dummy_value_rounded;

      $dummys{$i}{value} = $dummy_value;

      print BLUE, "\n$i th bin of this data is less than: $current_min and the color will be $percents{$i}{color}\tLegend value = $dummy_value", RESET;
    }

    $percentiles{$probeset_data_column}{percentiles} = \%percents;
    $dummy_values{$probeset_data_column}{dummys} = \%dummys;
  }
  print BLUE, "\n", RESET;

  #3.) Go through each gene.  Foreach gene, go through each probeset for each data column and assign it to a track based on the percentiles range
  #    - Go through the genes of each chromosome and print blocks of genes to each file (instead of ~20,000 seperate files)
  #    - Once 100 genes or the next chromosome is reached, start the next block

  #Sort genes according to their start coordinate on this chromosome  (this will insure that adjacent genes are placed in the same file
  foreach my $chr (sort {$a cmp $b} keys %gene_probe_map){
    my $gene_count = 0;

    #Process the tracks for one chromosome at a time
    my %track_data;
    my $track_file_count = 1;

    #Get the gene list for this chromosome and sort it by the reference 'start_coord'
    my $gene_list_ref = $gene_probe_map{$chr}{gene_list};

    foreach my $gene_id (sort {$gene_list_ref->{$a}{start_coord} <=> $gene_list_ref->{$b}{start_coord}} keys %{$gene_list_ref}){
      $gene_count++;

      #Get all the probes for this gene and add to a track record
      my $probesets_ref = $gene_list_ref->{$gene_id}->{probesets};

      #Insert some dummy probeset values for this gene - One for each 'percentiles slot', these act as a legend for display
      #Data needed: $probesets_ref->{$probeset_id}->{start}, $probesets_ref->{$probeset_id}->{end}, $probesets_ref->{$probeset_id}->{strand},
      #             $data_value = $probeset_data{$probeset_id}{$probeset_data_column};

      #Create a dummy probeset value for each gene, for each data column, for each percentiles/bin range
      #Then add this dummy probeset to the probesets_ref object for this gene using a new unique probeset ID (continues on from the highest real one)

      #First create the new dummy probeset IDs (one for each percentiles/bin range)
      my @dummy_probeset_ids;
      for (my $i = 1; $i <= $number_bins; $i++){
	$max_probeset_id++;
	push(@dummy_probeset_ids, $max_probeset_id);
      }

      foreach my $probeset_data_column (@probeset_data_columns){
	my %dummys = %{$dummy_values{$probeset_data_column}{dummys}};

	foreach my $i (sort {$a <=> $b} keys %dummys){
	  my $dummy_probeset_id = $dummy_probeset_ids[$i-1];

	  #Create the dummy probeset record using the dummy value and coordinates for this gene
	  $probesets_ref->{$dummy_probeset_id}->{start} = ($gene_info{$gene_id}{chr_start} - 1000);
	  $probesets_ref->{$dummy_probeset_id}->{end} = ($gene_info{$gene_id}{chr_start} - 750);

	  $probesets_ref->{$dummy_probeset_id}->{strand} = $gene_info{$gene_id}{chr_strand};

	  $master_probeset_list{$dummy_probeset_id}{probe_type} = "Dummy";

	  $probeset_data{$dummy_probeset_id}{$probeset_data_column} = $dummys{$i}{value};
	}
      }

      #Only allow 50 genes per track file
      if ($gene_count == 50){
	$gene_count = 0;
	$track_file_count++;
      }

      #Set ultimate path of track file and store this info in the gene record for later reference
      my $file_name = "$chr"."_$track_file_count".".txt";
      my $file_path = "$ucsc_dir"."$file_name";

      my $ucsc_track = "$web_base"."$web_loc"."/"."$file_name".".gz";

      $gene_info{$gene_id}{track_file_path} = "$file_path";
      $gene_info{$gene_id}{ucsc_track} = "$ucsc_track";

      #Go through each probe for this gene, get intensity data for each experiment, and add to a track according to the percentiles range of intensity
      foreach my $probeset_id (sort {$a <=> $b} keys %{$probesets_ref}){


	#Go through each probeset data column for this probeset (i.e. mip101-1, 5fur-1, etc...)
	foreach my $probeset_data_column (@probeset_data_columns){
	  my $data_value = $probeset_data{$probeset_id}{$probeset_data_column};

	  #Avoid 'na' values
	  if ($data_value eq "na" || $data_value eq "-na"){
	    next();
	  }

	  #Determine which percentile range this data value falls into:
	  my %percents = %{$percentiles{$probeset_data_column}{percentiles}};

	  my $percentile_slot;
	  foreach my $i (sort {$a <=> $b} keys %percents){
	    if ($data_value <= $percents{$i}{percentile}){
	      $percentile_slot = $i;
	      last();
	    }
	  }

	  my $data_value_rounded = sprintf("%.1f", $data_value);

	  my $probeset_chr;
	  if ($chr =~ /chr/){
	    $probeset_chr = "$chr";
	  }else{
	    $probeset_chr = "chr"."$chr";
	  }


	  #Build the actual track line (or two in the case of junction probes)
	  my $feature_name = "$probeset_id"."_"."$data_value_rounded";
	  my $track_line = "$probeset_chr\talexa\tProbeSet\t$probesets_ref->{$probeset_id}->{start}\t$probesets_ref->{$probeset_id}->{end}\t.\t$probesets_ref->{$probeset_id}->{strand}\t.\t$feature_name";

	  #Build the track hash with the following levels:
	  #       track_data -> track_file_count -> experiment(probeset_data_column) -> percentile_slot -> probeset_record(actual track line) or gene_record

	  if ($track_data{$track_file_count}){

	    my $trackfile_gene_list_ref = $track_data{$track_file_count}{gene_list};
	    $trackfile_gene_list_ref->{$gene_id}->{tmp} = '';

	    my $data_column_ref = $track_data{$track_file_count}{data_columns};

	    #Track file has at least one entry, but what about this experiment (probeset_data_column)
	    if ($data_column_ref->{$probeset_data_column}){

	      my $percentile_slot_ref = $data_column_ref->{$probeset_data_column}->{percentile_slots};
	      #Experiment has at least one entry, but what about this percentile slot
	      if ($percentile_slot_ref->{$percentile_slot}){
		my $probeset_record_ref = $percentile_slot_ref->{$percentile_slot}->{probeset_records};
		$probeset_record_ref->{$probeset_id}->{track_line} = $track_line;

	      }else{
		#Initialize this probeset
		my %probeset_records;
		$probeset_records{$probeset_id}{track_line} = $track_line;
		$percentile_slot_ref->{$percentile_slot}->{probeset_records} = \%probeset_records;
	      }

	    }else{
	      #Initialize this percentile_slot, and probeset
	      my %percentile_slots;
	      my %probeset_records;

	      $probeset_records{$probeset_id}{track_line} = $track_line;
	      $percentile_slots{$percentile_slot}{probeset_records} = \%probeset_records;
	      $data_column_ref->{$probeset_data_column}->{percentile_slots} = \%percentile_slots;
	      $data_column_ref->{$probeset_data_column}->{display_order} = $data_names{$probeset_data_column}{display_order};
	    }

	  }else{
	    #Initialize this data_column, percentile_slot, and probeset - and the filepath
	    my %data_columns;
	    my %percentile_slots;
	    my %probeset_records;
	    my %trackfile_gene_list; #Keep track of the 50 genes in this trackfile

	    $probeset_records{$probeset_id}{track_line} = $track_line;

	    $trackfile_gene_list{$gene_id}{tmp} = '';

	    $percentile_slots{$percentile_slot}{probeset_records} = \%probeset_records;
	    $data_columns{$probeset_data_column}{percentile_slots} = \%percentile_slots;
	    $data_columns{$probeset_data_column}{display_order} = $data_names{$probeset_data_column}{display_order};
	    $track_data{$track_file_count}{data_columns} = \%data_columns;
	    $track_data{$track_file_count}{file_path} = $file_path;
	    $track_data{$track_file_count}{gene_list} = \%trackfile_gene_list;
	  }
	}
      }
    }#Gene loop


    #Go through the track data object and add a single track for each gene, for each experiment, place it in the next available percentile slot?
    #      track_data -> track_file_count -> experiment(data_column) -> gene_records
    print BLUE, "\n\nAdding gene-level tracks", RESET;
    foreach my $track_file_count (sort {$a <=> $b} keys %track_data){
      my $trackfile_gene_list_ref = $track_data{$track_file_count}{gene_list};
      my $gene_counter = keys %{$trackfile_gene_list_ref};
      #print GREEN, "\nDEBUG: File - $track_file_count has $gene_counter gene IDs listed", RESET;

      #Go through each gene in the list for this trackfile
      foreach my $gene_id (keys %{$trackfile_gene_list_ref}){

	#Go through each possible data column
	foreach my $probeset_data_column (@probeset_data_columns){

	  #Make sure there is gene data that corresponds to this gene and data_column
	  unless (($gene_data{$gene_id}) && ($gene_data{$gene_id}{$probeset_data_column})){
	    next();
	  }

	  #Create actual tracklines for the exons of this gene
	  my @tracks;  #Exon tracks for a single gene
	  my $data_value = $gene_data{$gene_id}{$probeset_data_column};

	  #Determine the color appropriate for this gene
	  my %percents = %{$percentiles{$probeset_data_column}{percentiles}};

	  my $percentile_slot;
	  foreach my $i (sort {$a <=> $b} keys %percents){
	    if ($data_value <= $percents{$i}{percentile}){
	      $percentile_slot = $i;
	      last();
	    }
	  }

	  my $data_value_rounded = sprintf("%.1f", $data_value);
	  my $feature_name = "$gene_id"."_"."$data_value_rounded";

	  my $gene_chr = "chr"."$gene_info{$gene_id}{chromosome}";

	  my $exon_content_ref = $gene_info{$gene_id}{exon_content};

	  foreach my $exon_count (sort {$a <=> $b} keys %{$exon_content_ref}){
	    my $chr_start = $exon_content_ref->{$exon_count}->{chr_start};
	    my $chr_end = $exon_content_ref->{$exon_count}->{chr_end};

	    my $track_line = "$gene_chr\talexa\tGene\t$chr_start\t$chr_end\t.\t$gene_info{$gene_id}{chr_strand}\t.\t$feature_name";
	    push(@tracks, $track_line);
	  }

	  #Create a track for this gene and data column
	  my $data_column_ref = $track_data{$track_file_count}{data_columns};

	  if ($data_column_ref->{$probeset_data_column}->{gene_color_group}){

	    #Continuation of genes for this data column. Check if one with the same color (expression level) had already been added
	    my $gene_color_group_ref = $data_column_ref->{$probeset_data_column}->{gene_color_group};

	    if ($gene_color_group_ref->{$percentile_slot}){
	      my $gene_records_ref = $gene_color_group_ref->{$percentile_slot}->{gene_records};
	      $gene_records_ref->{$gene_id}->{tracks} = \@tracks;

	    }else{
	      my %gene_records;
	      $gene_records{$gene_id}{tracks} = \@tracks;
	      $gene_color_group_ref->{$percentile_slot}->{gene_records} = \%gene_records;
	    }

	  }else{
	    #First gene record for this data column
	    my %gene_color_group;
	    my %gene_records; #All the genes in this block of 50 that will be colored the same
	    $gene_records{$gene_id}{tracks} = \@tracks;
	    $gene_color_group{$percentile_slot}{gene_records} = \%gene_records;
	    $data_column_ref->{$probeset_data_column}->{gene_color_group} = \%gene_color_group;
	  }
	}
      }
    }

    #4) Now that all of the genes for a chromosome have been processed, write out probeset and gene tracks organized and stored in %track_data (reset for each chromosome above)
    foreach my $track_file_count (sort {$a <=> $b} keys %track_data){

      print BLUE, "\nPrinting UCSC file: $track_data{$track_file_count}{file_path}", RESET;

      #Open output file
      open (UCSC,">$track_data{$track_file_count}{file_path}") || die "\nCould not open file for output: $track_data{$track_file_count}{file_path}\n\n";

      #Set default browser settings
      print UCSC "#Browser line";
      #print UCSC "\nbrowser position chr$genes{$gene_id}{chromosome}:$genes{$gene_id}{chr_start}-$genes{$gene_id}{chr_end}";  #Default start pos
      print UCSC "\nbrowser full knownGene refGene ensGene";
      print UCSC "\nbrowser pack multiz8way";

      my $data_column_ref = $track_data{$track_file_count}{data_columns};

      #Go through each data column - sort by the user specified display order
      foreach my $probeset_data_column (sort {$data_column_ref->{$a}->{display_order} <=> $data_column_ref->{$b}->{display_order}} keys %{$data_column_ref}){

	#A.) Print out data for each gene (if data is available for this data column)
	if ($data_column_ref->{$probeset_data_column}->{gene_color_group}){


	  #Go through each block of genes that will be colored the same
	  my $gene_color_group_ref = $data_column_ref->{$probeset_data_column}->{gene_color_group};

	  foreach my $percentile_slot (sort {$b <=> $a} keys %{$gene_color_group_ref}){

	    #Create new track
	    print UCSC "\n\n#GENE data from column: $probeset_data_column";
	    my $gene_track_name = "$data_names{$probeset_data_column}{name}"."_"."gene"."_"."$percentile_slot";
	    #my $gene_color = "137,61,190";
	    my %percents = %{$percentiles{$probeset_data_column}{percentiles}};
	    my $gene_color = $percents{$percentile_slot}{color};

	    print UCSC "\ntrack name=$gene_track_name description=\"$data_names{$probeset_data_column}{description}\" color=$gene_color useScore=0 visibility=2";

	    my $gene_records_ref = $gene_color_group_ref->{$percentile_slot}->{gene_records};

	    foreach my $gene_id (sort {$a <=> $b} keys %{$gene_records_ref}){
	      my @tracks = @{$gene_records_ref->{$gene_id}->{tracks}};

	      foreach my $gene_track (@tracks){
		print UCSC "\n$gene_track";
	      }
	    }
	  }
	}

	#B.) Print out tracks for the probeset data for each intensity bin
	my $percentile_slot_ref = $data_column_ref->{$probeset_data_column}->{percentile_slots};

	foreach my $percentile_slot (sort {$b <=> $a} keys %{$percentile_slot_ref}){

	  my %percents = %{$percentiles{$probeset_data_column}{percentiles}};
	  my $color = $percents{$percentile_slot}{color};

	  my $bin_value = $percents{$percentile_slot}{dummy_value};

	  #Create new track
	  print UCSC "\n\n#Data from column: $probeset_data_column for bin: $percentile_slot";
	  my $track_name = "$data_names{$probeset_data_column}{name}"."_"."$percentile_slot"."_"."$bin_value";
	  print UCSC "\ntrack name=$track_name description=\"$data_names{$probeset_data_column}{description}\" color=$color useScore=0 visibility=$data_names{$probeset_data_column}{visibility}";

	  my $probeset_records_ref = $percentile_slot_ref->{$percentile_slot}->{probeset_records};

	  foreach my $probeset_id (sort {$a <=> $b} keys %{$probeset_records_ref}){

	    #Populate track with actual track lines
	    print UCSC "\n$probeset_records_ref->{$probeset_id}->{track_line}";
	  }
	}
      }

      #5.) Create a final set of tracks for the genes in this track file (ORFs, protein motifs, TMDs, SPs)
      my $trackfile_gene_list_ref = $track_data{$track_file_count}{gene_list};

      #5-A.) ORFs
      print BLUE, "\nPrinting ORF and protein feature tracks\n\n", RESET;
      print UCSC "\n\n#ORF Data";
      my $track_name = "EnsEMBL_ORFs";
      print UCSC "\ntrack name=$track_name description=\"EnsEMBL ORFs\" color=51,153,0 useScore=0 visibility=3";

      my @current_gene_ids = keys %{$trackfile_gene_list_ref};
      my $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@current_gene_ids, '-sequence'=>"no");

      foreach my $gene_id (sort {$a <=> $b} keys %{$trackfile_gene_list_ref}){

	my $trans_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};

	foreach my $trans_id (sort {$a <=> $b} keys %{$trans_ref}){
	  my $cds_start = $trans_ref->{$trans_id}->{cds_start};
	  my $cds_end = $trans_ref->{$trans_id}->{cds_end};

	  if ($cds_start eq "na" || $cds_end eq "na"){next();}

	  #Get ORFs coords
	  my %cds_coords;

	  my @cds_start_coords = @{$trans_ref->{$trans_id}->{cds_start_coords}};
	  my @cds_end_coords = @{$trans_ref->{$trans_id}->{cds_end_coords}};

	  my $exon_count = 0;
	  foreach my $cds_start_coord (@cds_start_coords){
	    $exon_count++;
	    my $cds_end_coord = shift(@cds_end_coords);
	    $cds_coords{$exon_count}{start} = $cds_start_coord;
	    $cds_coords{$exon_count}{end} = $cds_end_coord;
	  }

	  #Convert ORF coords to chromosome coords
	  my $gene_info_ref = \%gene_info;
	  my $chr_start = $gene_info_ref->{$gene_id}->{chr_start};
	  my $chr_end = $gene_info_ref->{$gene_id}->{chr_end};
	  my $chr_strand = $gene_info_ref->{$gene_id}->{chr_strand};

	  foreach my $exon_count (keys %cds_coords){
	    my $start = $cds_coords{$exon_count}{start};
	    my $end = $cds_coords{$exon_count}{end};

	    #Convert provided gene coordinates to coordinates relative to the chromosome
	    if ($chr_strand eq "+"){
	      my $query_chr_start = $chr_start + $start - 1;
	      my $query_chr_end = $chr_start + $end - 1;
	      $cds_coords{$exon_count}{chr_start} = $query_chr_start;
	      $cds_coords{$exon_count}{chr_end} = $query_chr_end;

	    }elsif ($chr_strand eq "-"){
	      my $query_chr_start = $chr_end - $end + 1;
	      my $query_chr_end = $chr_end - $start + 1;
	      $cds_coords{$exon_count}{chr_start} = $query_chr_start;
	      $cds_coords{$exon_count}{chr_end} = $query_chr_end;
	    }else{
	      print RED, "\nStrand format: $chr_strand not understood by convertGeneCoordinates()!\n\n", RESET;
	      exit();
	    }
	  }

	  #Write out the track values
	  foreach my $exon_id (keys %cds_coords){
	    my $feature_name = "ORF"."_"."$trans_id";
	    my $gene_chr = "chr"."$gene_info{$gene_id}{chromosome}";
	    print UCSC "\n$gene_chr\talexa\tORF\t$cds_coords{$exon_id}{chr_start}\t$cds_coords{$exon_id}{chr_end}\t.\t$chr_strand\t.\t$feature_name";
	  }
	}
      }

      #Tracks for protein features
      my $gene_pfs_ref = &getProteinFeatures ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@current_gene_ids);

      #Get chromosome coordinates for all protein features
      my $gene_info_ref = \%gene_info;
      foreach my $gene_id (sort {$a <=> $b} keys %{$gene_pfs_ref}){
	my $chr_start = $gene_info_ref->{$gene_id}->{chr_start};
	my $chr_end = $gene_info_ref->{$gene_id}->{chr_end};
	my $chr_strand = $gene_info_ref->{$gene_id}->{chr_strand};

	my $trans_ref = $gene_pfs_ref->{$gene_id}->{transcripts};

	my @nr_pfs = @{$gene_pfs_ref->{$gene_id}->{nr_pfs}};

	foreach my $trans_id (sort {$a <=> $b} keys %{$trans_ref}){
	  my $pfs_ref = $trans_ref->{$trans_id}->{protein_features};

	  foreach my $pf_id (@nr_pfs){
	    unless ($pfs_ref->{$pf_id}){next();}  #Non-Redundant PFs were created at gene level but this causes problem as we cycle through transcripts

	    my @start_coords = @{$pfs_ref->{$pf_id}->{start_coords}};
	    my @end_coords = @{$pfs_ref->{$pf_id}->{end_coords}};

	    my @chr_starts;
	    my @chr_ends;
	    foreach my $start (@start_coords){
	      my $end = shift (@end_coords);

	      if ($chr_strand eq "+"){
		my $query_chr_start = $chr_start + $start - 1;
		my $query_chr_end = $chr_start + $end - 1;
		push (@chr_starts, $query_chr_start);
		push (@chr_ends, $query_chr_end);

	      }elsif ($chr_strand eq "-"){
		my $query_chr_start = $chr_end - $end + 1;
		my $query_chr_end = $chr_end - $start + 1;
		push (@chr_starts, $query_chr_start);
		push (@chr_ends, $query_chr_end);

	      }else{
		print RED, "\nStrand format: $chr_strand not understood by convertGeneCoordinates()!\n\n", RESET;
		exit();
	      }
	      $pfs_ref->{$pf_id}->{chr_start_coords} = \@chr_starts;
	      $pfs_ref->{$pf_id}->{chr_end_coords} = \@chr_ends;
	    }
	  }
	}
      }

      #5-B.) Signal peptides
      print UCSC "\n\n#SignalPeptide Data";
      $track_name = "EnsEMBL_SignalPeptides";
      print UCSC "\ntrack name=$track_name description=\"EnsEMBL Signal Peptide Data\" color=204,0,0 useScore=0 visibility=3";
      foreach my $gene_id (sort {$a <=> $b} keys %{$gene_pfs_ref}){
	my $gene_chr = "chr"."$gene_info{$gene_id}{chromosome}";
	my $chr_strand = $gene_info{$gene_id}{chr_strand};

	my $trans_ref = $gene_pfs_ref->{$gene_id}->{transcripts};

	my @nr_pfs = @{$gene_pfs_ref->{$gene_id}->{nr_pfs}};

	foreach my $trans_id (sort {$a <=> $b} keys %{$trans_ref}){
	  my $pfs_ref = $trans_ref->{$trans_id}->{protein_features};

	  foreach my $pf_id (@nr_pfs){
	    unless ($pfs_ref->{$pf_id}){next();}  #Non-Redundant PFs were created at gene level but this causes problem as we cycle through transcripts

	    unless ($pfs_ref->{$pf_id}->{type} eq "Signalp"){next();}

	    my $name = $pfs_ref->{$pf_id}->{name};
	    my @chr_starts = @{$pfs_ref->{$pf_id}->{chr_start_coords}};
	    my @chr_ends = @{$pfs_ref->{$pf_id}->{chr_end_coords}};

	    my $feature_name;
	    if ($name eq "na"){
	      $feature_name = "SP_"."$pf_id";
	    }else{
	      $feature_name = "SP_"."$pf_id"."_"."$name";
	    }

	    foreach my $chr_start (@chr_starts){
	      my $chr_end = shift(@chr_ends);
	      print UCSC "\n$gene_chr\talexa\tSignalP\t$chr_start\t$chr_end\t.\t$chr_strand\t.\t$feature_name";
	    }
	  }
	}
      }

      #5-C.) Transmembrane domains
      print UCSC "\n\n#TransMembrane Data";
      $track_name = "EnsEMBL_TMDs";
      print UCSC "\ntrack name=$track_name description=\"EnsEMBL Transmembrane Domain Data\" color=102,0,102 useScore=0 visibility=3";
      foreach my $gene_id (sort {$a <=> $b} keys %{$gene_pfs_ref}){
	my $gene_chr = "chr"."$gene_info{$gene_id}{chromosome}";
	my $chr_strand = $gene_info{$gene_id}{chr_strand};

	my $trans_ref = $gene_pfs_ref->{$gene_id}->{transcripts};

	my @nr_pfs = @{$gene_pfs_ref->{$gene_id}->{nr_pfs}};

	foreach my $trans_id (sort {$a <=> $b} keys %{$trans_ref}){
	  my $pfs_ref = $trans_ref->{$trans_id}->{protein_features};

	  foreach my $pf_id (@nr_pfs){
	    unless ($pfs_ref->{$pf_id}){next();}  #Non-Redundant PFs were created at gene level but this causes problem as we cycle through transcripts

	    unless ($pfs_ref->{$pf_id}->{type} eq "tmhmm"){next();}

	    my $name = $pfs_ref->{$pf_id}->{name};
	    my @chr_starts = @{$pfs_ref->{$pf_id}->{chr_start_coords}};
	    my @chr_ends = @{$pfs_ref->{$pf_id}->{chr_end_coords}};

	    my $feature_name;
	    if ($name eq "na"){
	      $feature_name = "TMD_"."$pf_id";
	    }else{
	      $feature_name = "TMD_"."$pf_id"."_"."$name";
	    }

	    foreach my $chr_start (@chr_starts){
	      my $chr_end = shift(@chr_ends);
	      print UCSC "\n$gene_chr\talexa\tTMD\t$chr_start\t$chr_end\t.\t$chr_strand\t.\t$feature_name";
	    }
	  }
	}
      }

      #5-D.) Coiled-coil domains
      print UCSC "\n\n#Coiled-coil Data";
      $track_name = "EnsEMBL_CoiledCoils";
      print UCSC "\ntrack name=$track_name description=\"EnsEMBL Coiled-coil Domain Data\" color=102,0,255 useScore=0 visibility=3";
      foreach my $gene_id (sort {$a <=> $b} keys %{$gene_pfs_ref}){
	my $gene_chr = "chr"."$gene_info{$gene_id}{chromosome}";
	my $chr_strand = $gene_info{$gene_id}{chr_strand};

	my $trans_ref = $gene_pfs_ref->{$gene_id}->{transcripts};

	my @nr_pfs = @{$gene_pfs_ref->{$gene_id}->{nr_pfs}};

	foreach my $trans_id (sort {$a <=> $b} keys %{$trans_ref}){
	  my $pfs_ref = $trans_ref->{$trans_id}->{protein_features};

	  foreach my $pf_id (@nr_pfs){
	    unless ($pfs_ref->{$pf_id}){next();}  #Non-Redundant PFs were created at gene level but this causes problem as we cycle through transcripts

	    unless ($pfs_ref->{$pf_id}->{type} eq "ncoils"){next();}

	    my $name = $pfs_ref->{$pf_id}->{name};
	    my @chr_starts = @{$pfs_ref->{$pf_id}->{chr_start_coords}};
	    my @chr_ends = @{$pfs_ref->{$pf_id}->{chr_end_coords}};

	    my $feature_name;
	    if ($name eq "na"){
	      $feature_name = "CC_"."$pf_id";
	    }else{
	      $feature_name = "CC_"."$pf_id"."_"."$name";
	    }

	    foreach my $chr_start (@chr_starts){
	      my $chr_end = shift(@chr_ends);
	      print UCSC "\n$gene_chr\talexa\tCoiledCoil\t$chr_start\t$chr_end\t.\t$chr_strand\t.\t$feature_name";
	    }
	  }
	}
      }

      #5-E.) Protein motifs
      print UCSC "\n\n#Protein Motif Data";
      $track_name = "EnsEMBL_ProtMotif";
      print UCSC "\ntrack name=$track_name description=\"EnsEMBL Protein Motif Data\" color=102,0,0 useScore=0 visibility=3";
      foreach my $gene_id (sort {$a <=> $b} keys %{$gene_pfs_ref}){
	my $gene_chr = "chr"."$gene_info{$gene_id}{chromosome}";
	my $chr_strand = $gene_info{$gene_id}{chr_strand};

	my $trans_ref = $gene_pfs_ref->{$gene_id}->{transcripts};

	my @nr_pfs = @{$gene_pfs_ref->{$gene_id}->{nr_pfs}};

	foreach my $trans_id (sort {$a <=> $b} keys %{$trans_ref}){
	  my $pfs_ref = $trans_ref->{$trans_id}->{protein_features};

	  foreach my $pf_id (@nr_pfs){
	    unless ($pfs_ref->{$pf_id}){next();}  #Non-Redundant PFs were created at gene level but this causes problem as we cycle through transcripts

	    my $type = $pfs_ref->{$pf_id}->{type};
	    unless ($type eq "Pfam" || $type eq "pfscan" || $type eq "Prints" || $type eq "scanprosite"){next();}

	    if ($type eq "scanprosite"){
	      $type = "ScanPS";
	    }

	    my $name = $pfs_ref->{$pf_id}->{name};
	    my @chr_starts = @{$pfs_ref->{$pf_id}->{chr_start_coords}};
	    my @chr_ends = @{$pfs_ref->{$pf_id}->{chr_end_coords}};

	    my $feature_name;
	    if ($name eq "na"){
	      $feature_name = "$type"."_"."$pf_id";
	    }else{
	      $feature_name = "$type"."_"."$pf_id"."_"."$name";
	    }

	    foreach my $chr_start (@chr_starts){
	      my $chr_end = shift(@chr_ends);
	      print UCSC "\n$gene_chr\talexa\tCoiledCoil\t$chr_start\t$chr_end\t.\t$chr_strand\t.\t$feature_name";
	    }
	  }
	}
      }
      close (UCSC);

    }#track file loop

  }#chromosome loop

  return();
}


#################################################################################################################################
#5.) Now generate a summary file containing basic info for all genes and a link to the custom UCSC track)                       #
#################################################################################################################################
sub printSummaryFile{

  my %args = @_;
  my $summary_file = $args{'-summary_file'};

  print BLUE, "\n\nPrinting output summary file: $summary_file\n\n", RESET;

  open (SUMMARY, ">$summary_file") || die "\nCould not open summary file: $summary_file\n\n";

  print SUMMARY "gene_id\tensembl_g_id\tgene_symbol\tprobe_count\tchromosome\tstrand\tstart\tend";

  foreach my $gene_data_column (@gene_data_columns){
    print SUMMARY "\t$gene_data_column";
  }
  print SUMMARY "\tucsc_link\n";

  foreach my $gene_id (sort {$a <=> $b} keys %gene_info){

    my $chr_corrected = "chr"."$gene_info{$gene_id}{chromosome}";

    #Shift the genome view over by 100 bp to allow the dummy legend records to be seen
    my $gene_start = ($gene_info{$gene_id}{chr_start} - 1000);

    #Note: 'ctfile_hg18=' makes sure the custom tracks are loaded fresh, preventing old tracks from displaying as well.  Make sure this matches...
    my $link = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=$genome_build&position=$chr_corrected:$gene_start-$gene_info{$gene_id}{chr_end}&hgt.customText=$gene_info{$gene_id}{ucsc_track}&ctfile_$genome_build=";

    print SUMMARY "$gene_id\t$gene_info{$gene_id}{ensembl_g_id}\t$gene_info{$gene_id}{gene_name}\t$master_gene_list{$gene_id}{probe_count}\t$gene_info{$gene_id}{chromosome}\t$gene_info{$gene_id}{chr_strand}\t$gene_info{$gene_id}{chr_start}\t$gene_info{$gene_id}{chr_end}";

    #Print out gene values in the summary file
    foreach my $gene_data_column (@gene_data_columns){
         my $data_value = $gene_data{$gene_id}{$gene_data_column};

	 #Make sure there is gene data that corresponds to this gene and data_column
	 unless (($gene_data{$gene_id}) &&($gene_data{$gene_id}{$gene_data_column})){
	   print SUMMARY "\tna";
	   next();
	 }
	 print SUMMARY "\t$gene_data{$gene_id}{$gene_data_column}";

    }

    print SUMMARY "\t$link\n";

  }
  close (SUMMARY);
  return();
}
