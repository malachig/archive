#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to build an exonic bases database (as a Berkley DB on disk) with one entry for every non-redundant exonic base in EnsEMBL
#This will essentially be a data structure consisting of one file per EnsEMBL chromosome
#Each file will contain an hash keyed on chromosome position (simple integer)
#The value of each element of the hash will contain useful info such as 36-mer mappability score for that base
#Other useful info specific to individual exonic base positions can be added as needed
#The value of each element of the hash will be a tab delimited string

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

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $mappability_dir = '';
my $working_dir = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'mappability_dir=s'=>\$mappability_dir, 'working_dir=s'=>\$working_dir);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script builds a Berkley DB (one file per chromosome) and stores a record for every ensembl exonic base position", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify a directory containing mapability scores for each exonic base position", RESET;
print GREEN, "\n\tSpecify the working directory to write binary tree files to using: --working_dir", RESET;
print GREEN, "\n\nExample: buildExonicBasesDatabase.pl  --database=ALEXA_hs_49_36k  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --working_dir=/projects/malachig/sequence_databases/ensembl_exonic_bases_hs_49_36k_berkDB  --mappability_dir=/projects/malachig/solexa/mapability/ENST_v49_Exons_36mers/\n\n", RESET;

unless ($database && $server && $user && $password && $mappability_dir && $working_dir){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}


#Check and clean this working dir before proceeding
$working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"yes");

#Check the mappability dir to make sure it exists
$mappability_dir = &checkDir('-dir'=>$mappability_dir);


#Get gene info for all EnsEMBL genes
my $genes_ref;
my $gene_exon_content_ref;
my %exonic_bases_chr;

&getBasicGeneInfo();

#Now add mapability values to the Berkley DB files
&addMappabilityScores('-dir'=>$mappability_dir);


#Untie all the berkley db files
foreach my $chr (keys %exonic_bases_chr){
  my $temp_ref = $exonic_bases_chr{$chr};
  untie %{$temp_ref};
}

exit();


############################################################################################################################################
#Get basic info for all genes from the user specified ALEXA database                                                                       #
############################################################################################################################################
sub getBasicGeneInfo{
  my %args = @_;

  #Establish connection with the Alternative Splicing Expression database
  my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

  my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};

  #Get the gene info for all genes
  $| = 1; print BLUE, "\n1-a.) Getting gene data", RESET; $| = 0;
  my $g_storable_name = "$database"."_AllGenes_GeneInfo_NoSeq.storable";
  $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no", '-storable'=>$g_storable_name);


  #Get exon content for all genes
  $| = 1; print BLUE, "\n\n1-b.) Getting EXON CONTENT of each gene", RESET; $| = 0;
  my $ec_storable_name = "$database"."_AllGenes_ExonContent.storable";
  $gene_exon_content_ref = &getExonContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-storable'=>$ec_storable_name);

  my $total_gene_count = scalar (@gene_ids);

  #At this time, also get the chromosome coordinates for all exon-content coordinates
  $| = 1; print BLUE, "\n\n1-c.) Calculating chromosome coordinates for the EXON CONTENT of each gene (for only the current chromosome)\n", RESET; $| = 0;
  my $counter = 0;
  my $h_ref;
  my $last_chr;
  my $gene_count = 0;

  foreach my $gene_id (sort {$genes_ref->{$a}->{chromosome} cmp $genes_ref->{$b}->{chromosome}} keys %{$genes_ref}){

    my $chromosome = $genes_ref->{$gene_id}->{chromosome};

    if ($chromosome eq "MT"){$chromosome = "M"};

    ##########
    #Debug - load in a single chromosome only for testing purposes
    #unless ($chromosome eq "Y"){
    #  next();
    #}
    ##########

    $gene_count++;

    $counter++;
    if ($counter == 100){
      $counter = 0;
      $| = 1; print BLUE, ".", RESET; $| = 0;
    }
    my $chr_strand = $genes_ref->{$gene_id}->{chr_strand};
    my $chr_start = $genes_ref->{$gene_id}->{chr_start};
    my $chr_end = $genes_ref->{$gene_id}->{chr_end};
    my $gene_start = $genes_ref->{$gene_id}->{gene_start};
    my $gene_end = $genes_ref->{$gene_id}->{gene_end};

    #Calculate the size of each transcript by adding up the size of its exons
    my $size = 0;
    my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};
    my %exon_content_chr;

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

      }else{
	print RED, "\nStrand format: $chr_strand not understood !\n\n", RESET;
	exit();
      }
    }

    #If this genes corresponds to a new chromosome, create a berkley DB file for it
    #Otherwise simply get the neccessary hash reference
    if ($exonic_bases_chr{$chromosome}){
      #If the current chromosome has been observed before get the appropriate hash reference
      $h_ref = $exonic_bases_chr{$chromosome};

    }else{
      #If the current chromosome_strand combo has NOT been observed before, create the new hash and store as a berkleyDB file
      my %h;
      my $file_name = "$working_dir"."$chromosome".".exonic_bases.btree";
      print YELLOW, "\nCreating binary tree file: $file_name", RESET;

      tie(%h, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $file_name , -Flags => DB_CREATE) or die "can't open file $file_name: $! $BerkeleyDB::Error\n";

      $exonic_bases_chr{$chromosome} = \%h;
      $h_ref = $exonic_bases_chr{$chromosome};
    }


    #Populate it with the chromosome positions of the exon content of this gene

    foreach my $exon_id (sort {$exon_content_chr{$a}->{chr_start} <=> $exon_content_chr{$b}->{chr_start}} keys %exon_content_chr){

      my $start = $exon_content_chr{$exon_id}{chr_start};
      my $end = $exon_content_chr{$exon_id}{chr_end};

      #Go through each chromosome position in this exon content block and initialize that position in the hash
      #Initialize mappability scores to 'na'
      for (my $i = $start; $i <= $end; $i++){
	if ($h_ref->{$i}){
	  my $data = $h_ref->{$i};
	  my @data = split(" ", $data);
	  $data[0]++;
	  $h_ref->{$i} = "@data"
	}else{

	  #Data format.  key{position} value{exonic_base_coverage 36mer_mappability}
	  $h_ref->{$i} = "1 na";
	}
      }
    }

    #If the current chr. is different than the last one, since we are processing genes by chr. order, we must be done with the last chr.
    #In that case untie the berkeley db for now
    #if (($last_chr) && ($last_chr ne $chromosome || $gene_count == $total_gene_count)){
    #  my $temp_ref = $exonic_bases_chr{$last_chr};
    #  untie %{$temp_ref};
    #}

    $last_chr = $chromosome;
  }

  #Close database connection
  $alexa_dbh->disconnect();

  return();
}


############################################################################################################################################
#add Mappability Scores
############################################################################################################################################
sub addMappabilityScores{
  my %args = @_;
  my $dir = $args{'-dir'};

  opendir(DIRHANDLE, "$dir") || die "\nCannot open directory: $dir\n\n";
  my @temp = readdir(DIRHANDLE);
  closedir(DIRHANDLE);

  print BLUE, "\n\nAdding mappability scores to previously created Berkley DB files\n\n", RESET;

  my %map_files;
  foreach my $temp_file (@temp){

    if ($temp_file =~ /(\w+Chr)(\w+)(\.txt)/){
      my $file_path = "$dir"."$1"."$2"."$3";
      $map_files{$2}{file_path} = $file_path;
      $map_files{$2}{file_name} = "$1"."$2"."$3";
    }
  }

  #Go through each berkeley DB file created and add mappability scores if they are available for that chromosome/contig
  my $chr_defined_bases = 0;
  my $chr_undefined_bases = 0;

  foreach my $chr (sort keys %exonic_bases_chr){

    my $h_ref;
    if ($map_files{$chr}){

      #If a berkely DB file has been created for the current chromosome, get the appropriate hash reference
      print YELLOW, "\n\tFound a mapfile AND berkeley DB file for chr: $chr", RESET;
      $h_ref = $exonic_bases_chr{$chr};

    }else{
      #If a mapfile was not found for the current berkely DB file was NOT found, warn the user
      print YELLOW, "\n\tCould not find a mapfile for berkeley DB file chr: $chr - mappability positions for this chr will remain 'na'", RESET;
      next();
    }

    print BLUE, "\nProcessing mappability file: $chr ($map_files{$chr}{file_path})", RESET;


    #Open the mappability file
    open (MAP, "$map_files{$chr}{file_path}") || die "\nCould not open mappability file: $map_files{$chr}{file_path}\n\n";

    while(<MAP>){
      chomp ($_);
      my @line = split ("\t", $_);
      my $pos = ($line[1])+1;
      my $map_val = $line[2];

      #Shorten the mappability score to 4 decimals
      my $map_val_f = sprintf("%.4f", $map_val);

      if ($h_ref->{$pos}){
	my $data = $h_ref->{$pos};
	my @data = split(" ", $data);
	$data[1] = $map_val_f;
	$h_ref->{$pos} = "@data";
	$chr_defined_bases++;

      }else{
	print YELLOW, "\n\t\tWarning: position $pos found in MAP file is not defined in the corresponding exonic bases Berkeley DB!", RESET;
	$chr_undefined_bases++;

      }
    }
    close (MAP);
    print BLUE, "\n\nA total of $chr_defined_bases bases from the mappability file matched known exonic bases ($chr_undefined_bases bases did NOT)\n\n", RESET;
  }

  return();
}
