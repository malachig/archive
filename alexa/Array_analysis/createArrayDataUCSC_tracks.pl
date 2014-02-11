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

use lib '/home/malachig/AlternativeSplicing/perl_bin';
use Statistics::Descriptive;
use utilities::utility qw(:all);
use utilities::ALEXA_DB_35_35h qw(:all);

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $level = '';    #Level of analysis specified by user - must be 'probe' or 'probeset'
my $simple = '';   #If --simple=yes then only use Exon, Intron and Canonical probes in the UCSC tracks
my $map_file = '';
my $probeset_expression_file = '';
my $metaprobeset_expression_file = '';
my $splice_index_file = '';
my $ucsc_dir = '';
my $summary_file = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'level=s'=>\$level, 'simple=s'=>\$simple,
	    'map_file=s'=>\$map_file, 'probeset_expression_file=s'=>\$probeset_expression_file, 'metaprobeset_expression_file=s'=>\$metaprobeset_expression_file,
	    'splice_index_file=s'=>\$splice_index_file, 'ucsc_dir=s'=>\$ucsc_dir, 'summary_file=s'=>\$summary_file);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify whether the analysis should be conducted at the 'probe' or 'probeset' level using: --level", RESET;
print GREEN, "\n\tTo create tracks with only Exon, Intron, and Canonical probes or probesets use: --simple=yes", RESET;
print GREEN, "\n\tSpecify the input mapfile describing the maping of probesets to each Ensembl gene using --map_file", RESET;
print GREEN, "\n\tSpecify the input probeset expression data file using: --probeset_expression_file", RESET;
print GREEN, "\n\tSpecify the input meta-probeset expression data file using: --metaprobeset_expression_file", RESET;
print GREEN, "\n\tSpecify the input splice index data file using: --splice_index_file (OPTIONAL)", RESET;
print GREEN, "\n\tSpecify the full path to the target directory for UCSC track files using: --ucsc_dir", RESET;
print GREEN, "\n\tSpecify the output summary file name using: --summary_file", RESET;
print GREEN, "\n\nExample: createArrayDataUCSC_tracks.pl  --database=ALEXA_hs_35_35h  --server=jango.bcgsc.ca  --user=malachig  --password=pwd  --level=probeset  --simple=no  --map_file=/home/malachig/AlternativeSplicing/Array_design/ALEXA_hs_35_35h/probe_dev/NimbleGen_Design_Submission/NimbleGenArray_mipVS5FUR_coords  --probeset_expression_file=/home/malachig/AlternativeSplicing/Array_analysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/probeset_data.txt  --metaprobeset_expression_file=/home/malachig/AlternativeSplicing/Array_analysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/gene_data.txt  --splice_index_file=inclusion_index_values_1-3.txt  --ucsc_dir=/home/malachig/www/public/htdocs/alexa_tracks/mip_v_5fur/Raw/probe_level    --summary_file=gene_ucsc_summary_file_probeset-level.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $level && $simple && $map_file && $probeset_expression_file && $metaprobeset_expression_file && $ucsc_dir && $summary_file){
  print RED, "\nOptions missing!\n\n", RESET;
  exit();
}
unless ($level eq "probe" || $level eq "probeset"){
  print RED, "\nUser must specify level of analysis as 'probe' or 'probeset'\n\n";
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

#MIP - REPLICATE Log2 intensity values
$data_names{'MIP_EF_A_LOG2_I'}{name} = "mip_ef_a";
$data_names{'MIP_EF_A_LOG2_I'}{description} = "MIP101 (Sensitive) log2 intensity data - replicate EF_A";
$data_names{'MIP_EF_A_LOG2_I'}{colors} = "@colors_heat";
$data_names{'MIP_EF_A_LOG2_I'}{visibility} = 0;
$data_names{'MIP_EF_A_LOG2_I'}{display_order} = 1;

$data_names{'MIP_EF_B_LOG2_I'}{name} = "mip_ef_b";
$data_names{'MIP_EF_B_LOG2_I'}{description} = "MIP101 (Sensitive) log2 intensity data - replicate EF_B";
$data_names{'MIP_EF_B_LOG2_I'}{colors} = "@colors_heat";
$data_names{'MIP_EF_B_LOG2_I'}{visibility} = 0;
$data_names{'MIP_EF_B_LOG2_I'}{display_order} = 2;

$data_names{'MIP_GH_A_LOG2_I'}{name} = "mip_gh_a";
$data_names{'MIP_GH_A_LOG2_I'}{description} = "MIP101 (Sensitive) log2 intensity data - replicate GH_A";
$data_names{'MIP_GH_A_LOG2_I'}{colors} = "@colors_heat";
$data_names{'MIP_GH_A_LOG2_I'}{visibility} = 0;
$data_names{'MIP_GH_A_LOG2_I'}{display_order} = 3;


#FUR - REPLICATE Log2 intensity values
$data_names{'FUR_EF_A_LOG2_I'}{name} = "fur_ef_a";
$data_names{'FUR_EF_A_LOG2_I'}{description} = "5FUR (Resistant) log2 intensity data - replicate EF_A";
$data_names{'FUR_EF_A_LOG2_I'}{colors} = "@colors_heat";
$data_names{'FUR_EF_A_LOG2_I'}{visibility} = 0;
$data_names{'FUR_EF_A_LOG2_I'}{display_order} = 4;

$data_names{'FUR_EF_B_LOG2_I'}{name} = "fur_ef_b";
$data_names{'FUR_EF_B_LOG2_I'}{description} = "5FUR (Resistant) log2 intensity data - replicate EF_B";
$data_names{'FUR_EF_B_LOG2_I'}{colors} = "@colors_heat";
$data_names{'FUR_EF_B_LOG2_I'}{visibility} = 0;
$data_names{'FUR_EF_B_LOG2_I'}{display_order} = 5;

$data_names{'FUR_GH_A_LOG2_I'}{name} = "fur_gh_a";
$data_names{'FUR_GH_A_LOG2_I'}{description} = "5FUR (Resistant) log2 intensity data - replicate GH_A";
$data_names{'FUR_GH_A_LOG2_I'}{colors} = "@colors_heat";
$data_names{'FUR_GH_A_LOG2_I'}{visibility} = 0;
$data_names{'FUR_GH_A_LOG2_I'}{display_order} = 6;


#MIP vs FUR REPLICATE Log2 DE values
$data_names{'ef_a_LOG2_DE'}{name} = "DE_ef_a";
$data_names{'ef_a_LOG2_DE'}{description} = "MIP101-5FUR (Sensitive - Resistant) log2 DE data - replicate pair EF_A";
$data_names{'ef_a_LOG2_DE'}{colors} = "@colors_green_to_blue";
$data_names{'ef_a_LOG2_DE'}{visibility} = 1;
$data_names{'ef_a_LOG2_DE'}{display_order} = 7;

$data_names{'ef_b_LOG2_DE'}{name} = "DE_ef_b";
$data_names{'ef_b_LOG2_DE'}{description} = "MIP101-5FUR (Sensitive - Resistant) log2 DE data - replicate pair EF_B";
$data_names{'ef_b_LOG2_DE'}{colors} = "@colors_green_to_blue";
$data_names{'ef_b_LOG2_DE'}{visibility} = 1;
$data_names{'ef_b_LOG2_DE'}{display_order} = 8;

$data_names{'gh_a_LOG2_DE'}{name} = "DE_gh_a";
$data_names{'gh_a_LOG2_DE'}{description} = "MIP101-5FUR (Sensitive - Resistant) log2 DE data - replicate pair GH_A";
$data_names{'gh_a_LOG2_DE'}{colors} = "@colors_green_to_blue";
$data_names{'gh_a_LOG2_DE'}{visibility} = 1;
$data_names{'gh_a_LOG2_DE'}{display_order} = 9;


#MIP vs FUR GRAND MEAN (OR MEDIAN?) Log2 intensity values
#- Take mean or median of all probes in each probeset and across all replicates
$data_names{'mip_LOG2_I_U'}{name} = "mip_U";
$data_names{'mip_LOG2_I_U'}{description} = "MIP101 (Sensitive) MEAN log2 intensity data - (combines replicate 1+2+3)";
$data_names{'mip_LOG2_I_U'}{colors} = "@colors_heat";
$data_names{'mip_LOG2_I_U'}{visibility} = 1;
$data_names{'mip_LOG2_I_U'}{display_order} = 10;

$data_names{'fur_LOG2_I_U'}{name} = "fur_U";
$data_names{'fur_LOG2_I_U'}{description} = "5FUR (Resistant) MEAN log2 intensity data - (combines replicate 1+2+3)";
$data_names{'fur_LOG2_I_U'}{colors} = "@colors_heat";
$data_names{'fur_LOG2_I_U'}{visibility} = 1;
$data_names{'fur_LOG2_I_U'}{display_order} = 11;


#MIP vs FUR GRAND MEAN (OR MEDIAN?) Log2 DE values
#- Take mean or median of all DE values in each probeset and across all replicates
$data_names{'mip_v_fur_DE_U'}{name} = "DE_U";
$data_names{'mip_v_fur_DE_U'}{description} = "MIP101-5FUR (Sensitive - Resistant) MEAN log2 DE data - (combines replicate 1+2+3)";
$data_names{'mip_v_fur_DE_U'}{colors} = "@colors_green_to_blue";
$data_names{'mip_v_fur_DE_U'}{visibility} = 1;
$data_names{'mip_v_fur_DE_U'}{display_order} = 12;

#Probeset-level p-values.  All probes for each probesets are used to test for significant change between sensitive and resistant across all replicates
$data_names{'raw_wcox_pval'}{name} = "Pval_wcx";
$data_names{'raw_wcox_pval'}{description} = "Raw Wilcoxen P-values (Sens. vs Res.) - (combines all probes for each probeset across all reps)";
$data_names{'raw_wcox_pval'}{colors} = "@colors_white_to_purple";
$data_names{'raw_wcox_pval'}{visibility} = 1;
$data_names{'raw_wcox_pval'}{display_order} = 13;


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
# e.g. /home/malachig/AlternativeSplicing/Array_design/ALEXA_hs_35_35h/probe_dev/NimbleGen_Design_Submission/NimbleGenArray_mipVS5FUR_coords
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


#3-C.) Now get probeset splicing index values. (create 'na' entries for all probes that do not have a splice index value - such as non 'core' probesets)
#      - Hard code the column containing the splice index values
#if ($splice_index_file){
#  &parseProbesetSpliceIndexData('-splice_index_file'=>$splice_index_file);
#}

#4.) Now generate the custom UCSC tracks
#Establish connection with the Alternative Splicing Expression database
$alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
&generateUcscTracks('-ucsc_dir'=>$ucsc_dir);
$alexa_dbh->disconnect();

#5.) Now generate a summary file containing basic info for all genes and a link to the custom UCSC track)
&printSummaryFile('-summary_file'=>$summary_file);


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

  #First get a list of genes in the map file (all genes in the design)
  my %columns;
  open (MAPFILE, "$map_file") || die "\nCould not open mapfile: $map_file\n\n";

  my $first_line = 1;
  while(<MAPFILE>){
    chomp($_);
    my @line = split ("\t", $_);

    if ($first_line){
      #Watch for the header line which is assumed to contain neccessary columns
      my $column_count = 0;
      foreach my $column (@line){
        $columns{$column}{column_pos} = $column_count;
        $column_count++;
      }
      $first_line = 0;
      next();
    }

    #Skip negative control probes
    my $probe_type = $line[$columns{Probe_Type}{column_pos}];
    if ($probe_type eq "Control-Negative"){
      next();
    }

    my $gene_id = $line[$columns{Gene_ID}{column_pos}];
    if ($master_gene_list{$gene_id}){
      $master_gene_list{$gene_id}{probe_count}++;
    }else{
      $master_gene_list{$gene_id}{probe_count} = 1;
    }
  }
  close (MAPFILE);

  my $genes_found = keys %master_gene_list;

  print BLUE, "\nFound a total of $genes_found genes in the map file\n\n", RESET;

  #Now get neccessary additional info from the ALEXA database for these genes
  my @gene_ids = keys %master_gene_list;
  &getBasicGeneInfo('-dbh'=>$alexa_dbh, '-genes'=>\@gene_ids);

  #Now go through the mapfile again and build an object keyed on gene ids and probe or probeset_ids
  open (MAPFILE, "$map_file") || die "\nCould not open mapfile: $map_file\n\n";

  while(<MAPFILE>){
    my $line = $_;
    chomp($line);

    my @line = split ("\t", $line);

    if ($line[0] =~/^\d+/){
      $record_count++;

      #The 'probeset_id' will actually be a probe_id OR probeset_id depending on the level desired by the user

      my $probeset_id;
      if ($level eq "probeset"){
	$probeset_id = $line[$columns{ProbeSet_ID}{column_pos}];
      }elsif($level eq "probe"){
	$probeset_id = $line[$columns{Probe_Count}{column_pos}];
      }


      my $gene_id = $line[$columns{Gene_ID}{column_pos}];
      my $probe_type = $line[$columns{Probe_Type}{column_pos}];

      #Skip negative control probes
      if ($probe_type eq "Control-Negative"){
	next();
      }

      $master_probeset_list{$probeset_id}{probe_type} = $probe_type;

      if ($probeset_id > $max_probeset_id){
	$max_probeset_id = $probeset_id;
      }

      my $ensembl_gene_id = $gene_info{$gene_id}{ensembl_g_id};
      my $chromosome = $gene_info{$gene_id}{chromosome};
      my $strand = $gene_info{$gene_id}{chr_strand};

      #Determine the start/stop positions of features and probe_count (according to the level specified by the user)
      my $u1_start = $line[$columns{Unit1_start_chr}{column_pos}];
      my $u1_end = $line[$columns{Unit1_end_chr}{column_pos}];
      my $u2_start = $line[$columns{Unit2_start_chr}{column_pos}];
      my $u2_end = $line[$columns{Unit2_end_chr}{column_pos}];

      #Order the coordinates for simplicity
      my @coords;
      my ($unit1_start, $unit1_end, $unit2_start, $unit2_end);
      push (@coords, $u1_start);
      push (@coords, $u1_end);
      if ($u2_start eq "na" || $u2_end eq "na"){
	my @sorted_coords = sort {$a <=> $b} (@coords);
	$unit1_start = $sorted_coords[0];
	$unit1_end = $sorted_coords[1];
	$unit2_start = "na";
	$unit2_end = "na";

      }else{
	push (@coords, $u2_start);
	push (@coords, $u2_end);
	my @sorted_coords = sort {$a <=> $b} (@coords);

	$unit1_start = $sorted_coords[0];
	$unit1_end = $sorted_coords[1];
	$unit2_start = $sorted_coords[2];
	$unit2_end = $sorted_coords[3];
      }

      if ($gene_probe_map{$chromosome}){
	#Grab a copy of the hash reference stored in the gene_probe_map hash
	#Use this reference to directly add to the gene list of each chromosome
	my $gene_list_ref = $gene_probe_map{$chromosome}{gene_list};

	if ($gene_list_ref->{$gene_id}){

	  #Grab a copy of the hash reference stored in the gene_list hash
	  #Use this reference to directly add to the probesets hash for each gene
	  my $probesets_ref = $gene_list_ref->{$gene_id}->{probesets};

	  #Assign probe positions
	  if ($probesets_ref->{$probeset_id}){
	    if ($unit1_start < $probesets_ref->{$probeset_id}->{unit1_start}){
	      $probesets_ref->{$probeset_id}->{unit1_start} = $unit1_start;
	    }
	    if ($unit1_end > $probesets_ref->{$probeset_id}->{unit1_end}){
	      $probesets_ref->{$probeset_id}->{unit1_end} = $unit1_end;
	    }

	    if ($unit2_start eq "na" || $unit2_end eq "na"){
	      $probesets_ref->{$probeset_id}->{unit2_start} = $unit2_start;
	      $probesets_ref->{$probeset_id}->{unit2_end} = $unit2_end;
	    }else{
	      if ($unit2_start < $probesets_ref->{$probeset_id}->{unit2_start}){
		$probesets_ref->{$probeset_id}->{unit2_start} = $unit2_start;
	      }
	      if ($unit2_end > $probesets_ref->{$probeset_id}->{unit2_end}){
		$probesets_ref->{$probeset_id}->{unit2_end} = $unit2_end;
	      }
	    }
	  }else{
	    $probesets_ref->{$probeset_id}->{unit1_start} = $unit1_start;
	    $probesets_ref->{$probeset_id}->{unit1_end} = $unit1_end;
	    $probesets_ref->{$probeset_id}->{unit2_start} = $unit2_start;
	    $probesets_ref->{$probeset_id}->{unit2_end} = $unit2_end;
	  }

	  #Assign other values to this probeset
	  $probesets_ref->{$probeset_id}->{alexa_gene_id} = $gene_id;
	  $probesets_ref->{$probeset_id}->{ensembl_gene_id} = $ensembl_gene_id;
	  $probesets_ref->{$probeset_id}->{chromosome} = $chromosome;
	  $probesets_ref->{$probeset_id}->{strand} = $strand;


	  #See if the start coordinate for this probeset is smaller than the one currently stored for this gene
	  if ($unit1_start < $gene_list_ref->{$gene_id}->{start_coord}){
	    $gene_list_ref->{$gene_id}->{start_coord} = $unit1_start;
	  }

	}else{
	  #Chromosome has been seen before but not this gene
	  $gene_count++;
	  my %probesets;
	  $probesets{$probeset_id}{alexa_gene_id} = $gene_id;
	  $probesets{$probeset_id}{ensembl_gene_id} = $ensembl_gene_id;
	  $probesets{$probeset_id}{chromosome} = $chromosome;
	  $probesets{$probeset_id}{strand} = $strand;
	  $probesets{$probeset_id}{unit1_start} = $unit1_start;
	  $probesets{$probeset_id}{unit1_end} = $unit1_end;
	  $probesets{$probeset_id}{unit2_start} = $unit2_start;
	  $probesets{$probeset_id}{unit2_end} = $unit2_end;

	  $gene_list_ref->{$gene_id}->{start_coord} = $unit1_start;
	  $gene_list_ref->{$gene_id}->{probesets} = \%probesets;
	}

      }else{
	#First time this chromosome (and gene) has been observed
	$gene_count++;
	my %gene_list;
	my %probesets;
	$probesets{$probeset_id}{alexa_gene_id} = $gene_id;
	$probesets{$probeset_id}{ensembl_gene_id} = $ensembl_gene_id;
	$probesets{$probeset_id}{chromosome} = $chromosome;
	$probesets{$probeset_id}{strand} = $strand;
	$probesets{$probeset_id}{unit1_start} = $unit1_start;
	$probesets{$probeset_id}{unit1_end} = $unit1_end;
	$probesets{$probeset_id}{unit2_start} = $unit2_start;
	$probesets{$probeset_id}{unit2_end} = $unit2_end;

	$gene_list{$gene_id}{probesets} = \%probesets;
	$gene_list{$gene_id}{start_coord} = $unit1_start;  #Used to sort genes by the smallest coordinate of their probes
	$gene_probe_map{$chromosome}{gene_list} = \%gene_list;
      }
    }
  }

  close (MAPFILE);
  print BLUE, "\nProcessed: $record_count probeset records corresponding to $gene_count genes\n\n", RESET;
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
      my $probeset_id;
      if ($level eq "probeset"){
	$probeset_id = $line[$columns{ProbeSet_ID}{column_pos}];
      }elsif($level eq "probe"){
	$probeset_id = $line[$columns{Probe_ID}{column_pos}];
      }

      #Only import data for those that are in the mapfile
      unless ($master_probeset_list{$probeset_id}){
	next();
      }

      #Get the exons-skipped value if it is available
      if ($columns{Exons_Skipped}{column_pos}){
	#Add the exons-skipped value to the master_probeset_list hash for future reference
	$master_probeset_list{$probeset_id}{exons_skipped} = $line[$columns{Exons_Skipped}{column_pos}];
      }else{
	print RED, "\nExons_Skipped column not defined!\n\n", RESET;
	exit();
      }

      foreach my $data_column (@probeset_data_columns){
	$probeset_data{$probeset_id}{$data_column} = $line[$columns{$data_column}{column_pos}];

	#For p-values, change to -ve values so they get displayed in a sensible order
	if ($data_column eq "raw_wcox_pval"){
	  $probeset_data{$probeset_id}{$data_column} = -($line[$columns{$data_column}{column_pos}]);
	}
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

      my $gene_id = $line[$columns{AlexaGene_ID}{column_pos}];

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


################################################################################################################################
#3-C.) Now get probeset splicing index values. (create 'na' entries for all probes that do not have a splice index value       #
#      - Hard code the column containing the splice index values                                                               #
################################################################################################################################
sub parseProbesetSpliceIndexData{
  my %args = @_;
  my $splice_index_file = $args{'-splice_index_file'};

  my %splice_index_data;
  my $record_count = 0;

  print BLUE, "\nBegin importing splice index values from: $splice_index_file\n\n", RESET;

  open (SPLICEFILE, "$splice_index_file") || die "\nCould not open splice index file: $splice_index_file\n\n";

  my $header_line;
  my $header = 1;
  my $data_name;
  my @headers;
  my $column_count;

  while(<SPLICEFILE>){
    my $line = $_;
    chomp($line);

    #Skip comment lines if present
    if ($line =~ /^#/){
      next();
    }

    #Get the header line
    if ($header == 1){
      $header = 0;
      $header_line = $line;
      @headers = split ("\t", $header_line);
      $column_count = @headers;
    }

    my @line = split ("\t", $line);


    if ($line[0] =~/^\d+/){
      $record_count++;

      my $probeset_id = $line[0];

      #Note: the probesets mapped to Ensembl genes are NOT all probesets.
      #Only import data for those that are in the mapfile
      unless ($master_probeset_list{$probeset_id}){
	next();
      }

      for (my $i = 1; $i < $column_count; $i++){
	$splice_index_data{$probeset_id}{$headers[$i]} = $line[$i];
      }
    }
  }
  close (SPLICEFILE);
  print BLUE, "\nProcessed: $record_count splice index records", RESET;
  my $splice_index_count = keys %splice_index_data;

  #remove the probe_id column before adding column data to the probeset_data hash
  shift (@headers);

  print BLUE, "\nStored $splice_index_count splice index values", RESET;
  print BLUE, "\nAdding columns: @headers to the probeset_data_columns, probesets without a splice index value will be set to 'na'\n\n", RESET;

  #Go through the master probeset list and add records for the splice index values.  Probesets that did not have a splice index value will be set to 'na' 
  #This is neccessary because probeset expression data will be stored for all probes that map to an Ensembl gene BUT splice index values may only be calculated
  #for 'core' probesets (those that map to a known exon)
  foreach my $probeset_id (keys %master_probeset_list){
    if ($splice_index_data{$probeset_id}){

      #If this probeset ID has a splice index value AND is one of those included in the mapfile (maps to an Entrez gene), add data to the probeset_data hash
      foreach my $column_name (@headers){
	$probeset_data{$probeset_id}{$column_name} = $splice_index_data{$probeset_id}{$column_name};
      }

    }else{
      foreach my $column_name (@headers){
	$probeset_data{$probeset_id}{$column_name} = "na";
      }
    }
  }

  #Add this probe data column to the master list:
  foreach my $column_name (@headers){
    push (@probeset_data_columns, $column_name);
  }

  #Check each of the probeset_data columns to make sure they are defined in the %data_names hash
  foreach my $data_column (@probeset_data_columns){
    unless ($data_names{$data_column}){
      print RED, "\nData column: $data_column is not defined in the data_names object - check script\n\n", RESET;
      exit();
    }
  }
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
      #Data needed: $probesets_ref->{$probeset_id}->{start}, $probesets_ref->{$probeset_id}->{stop}, $probesets_ref->{$probeset_id}->{strand},
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
	  $probesets_ref->{$dummy_probeset_id}->{unit1_start} = ($gene_info{$gene_id}{chr_start} - 1000);
	  $probesets_ref->{$dummy_probeset_id}->{unit1_end} = ($gene_info{$gene_id}{chr_start} - 750);
	  $probesets_ref->{$dummy_probeset_id}->{unit2_start} = "na";
	  $probesets_ref->{$dummy_probeset_id}->{unit2_end} = "na";

	  $probesets_ref->{$dummy_probeset_id}->{strand} = $gene_info{$gene_id}{chr_strand};

	  $master_probeset_list{$dummy_probeset_id}{probe_type} = "Dummy";
	  $master_probeset_list{$dummy_probeset_id}{exons_skipped} = "na";

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

	#If the user specified the '--simple' option, only add probes that are 'Exon', 'Intron', or 'Canonical'
	#Skip all other probe types
	if ($simple eq "yes"){
	  my $probe_type = $master_probeset_list{$probeset_id}{probe_type};
	  my $exons_skipped =  $master_probeset_list{$probeset_id}{exons_skipped};

	  unless ($probe_type eq "Dummy" || $probe_type eq "Exon" || $probe_type eq "Intron" || ($probe_type eq "Exon-Exon" && $exons_skipped == 0)){
	    next();
	  }
	}

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

	  my $probeset_chr = "chr"."$chr";

	  #Build the actual track line (or two in the case of junction probes)
	  my $feature_name = "$probeset_id"."_"."$data_value_rounded";
	  my $track_line = "$probeset_chr\talexa\tProbeSet\t$probesets_ref->{$probeset_id}->{unit1_start}\t$probesets_ref->{$probeset_id}->{unit1_end}\t.\t$probesets_ref->{$probeset_id}->{strand}\t.\t$feature_name";

	  unless ($probesets_ref->{$probeset_id}->{unit2_start} eq "na" || $probesets_ref->{$probeset_id}->{unit2_end} eq "na"){
	    $track_line = "$track_line"."\n$probeset_chr\talexa\tProbeSet\t$probesets_ref->{$probeset_id}->{unit2_start}\t$probesets_ref->{$probeset_id}->{unit2_end}\t.\t$probesets_ref->{$probeset_id}->{strand}\t.\t$feature_name";
	  }

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

      my $cmd = "/usr/bin/gzip $track_data{$track_file_count}{file_path}";
      system($cmd);

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

    #Note: 'ctfile_hg17=' makes sure the custom tracks are loaded fresh, preventing old tracks from displaying as well.  Make sure the 'hg17' matches...
    my $link = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg17&position=$chr_corrected:$gene_start-$gene_info{$gene_id}{chr_end}&hgt.customText=$gene_info{$gene_id}{ucsc_track}&ctfile_hg17=";

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
