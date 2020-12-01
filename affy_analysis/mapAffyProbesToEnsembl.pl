#!/usr/bin/perl -w
#Written by Malachi Griffith

#This script takes Affy probe and probeset files as input and attempts to map every probeset to an ALEXA/Ensembl exon and gene

#1.)  Get a hash of all ALEXA genes and their gene coordinates
#     - Organize hash according to chromosome.
#2.)  Get probe/probeset data.  For now use the coordinates for the entire probeset.  This means that all probes in the probe set must be completely
#     bounded by an Ensembl Gene.  Probeset data is provided for both NCBI build 34 and 35.  Individual probe coordinates are only provided for build 34
#     - To use the probe coordinates will therefore require converting the coordinates.

#2a.) For each probeset, get the coordinates (chr, strand, start, stop).  Make note of the number of probes in each probeset
#2b.) Make sure all (or most?) probes in the probe set lie within the boundaries of an Ensembl Gene
#     - Will try using the probeset coordinates first.
#2c.) Some probesets may map to multiple genes (in the case of overlaping genes on the same strand), make note of such ambiguous cases
#3.)  Create an output file with the following data:
#     - Probe_set_ID, EnsemblGeneID, ProbeCount

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
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $probe_file = '';
my $filter = '';
my $out_file = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'probe_file=s'=>\$probe_file, 'filter=s'=>\$filter, 'out_file=s'=>\$out_file);

#Provide instruction to the user
print BLUE, "\n\nUsage:", RESET;
print BLUE, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print BLUE, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print BLUE, "\n\tSpecify the input Affy probe and probeset coordinate data file using: --probe_file", RESET;
print BLUE, "\n\tSpecify the filter to apply using: --filter ('all', 'known_non_pseudo', or 'non_pseudo')", RESET;
print BLUE, "\n\t\tThe 'non_pseudo' choice is recommended because it allows all genes (known or predicted) except pseudo genes", RESET;
print BLUE, "\n\tSpecify the output file using: --out_file", RESET;
print BLUE, "\n\nExample: mapAffyProbesToEnsembl.pl  --database=ALEXA_hs_35_35h  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --probe_file=HuEx-1_0-st-v2.annot.hg17.csv  --filter=non_pseudo  --out_file=test.out\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $probe_file && $filter && $out_file){
  print RED, "\nOptions missing!\n\n", RESET;
  exit();
}

unless ($filter eq "all" || $filter eq "known_non_pseudo" || $filter eq "non_pseudo"){
  print RED, "\nFilter option not understood\n\n", RESET;
  exit();
}

#1.) Get a hash of all ALEXA genes and their gene coordinates
#    - Create hash to organize according to chromosome.

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);


#For each gene, go through all the probes of that gene and get genomic coordinates
my %genes;
my %genes_chr;
&getGeneCoordinates('-dbh'=>$alexa_dbh);

#Close database connection
$alexa_dbh->disconnect();

#2.) Parse probe input data and map probesets to Ensembl genes
#    For each probeset, get the coordinates (chr, strand, start, stop) for each of the probes.  Make note of the number of probes in each probeset
#    - To save memory, immediately check the coordinates and only store those that correspond to an EnsEMBL gene
#    - Also print output to file immediately
&parseProbeData('-probe_file'=>$probe_file, '-out_file'=>$out_file);


exit();


#######################################################################################################################
#Get a hash of all ALEXA genes and their gene coordinates                                                             #
#######################################################################################################################
sub getGeneCoordinates{
  my %args = @_;
  my $dbh = $args{'-dbh'};

  #Get IDs for every gene in the database BUT omit those that are pseudogenes!
  #Also only allow genes that are 'Known Genes'?

  #($filter eq "all" || $filter eq "" || $filter eq "non_pseudo"){

  my @gene_ids;
  if ($filter eq "all"){
    print BLUE, "\nGetting ALL genes\n\n", RESET;
    @gene_ids = @{&getAllGenes ('-dbh'=>$dbh, '-gene_type'=>'All')};
  }elsif($filter eq "known_non_pseudo"){
    print BLUE, "\nGetting ALL genes except UNKNOWN and PSEUDO genes\n\n", RESET;
    @gene_ids = @{&getAllGenes ('-dbh'=>$dbh, '-gene_type'=>'Non-pseudo', '-evidence'=>'Known Gene')};
  }elsif($filter eq "non_pseudo"){
    print BLUE, "\nGetting ALL genes except PSEUDO genes\n\n", RESET;
    @gene_ids = @{&getAllGenes ('-dbh'=>$dbh, '-gene_type'=>'Non-pseudo')};
  }else{
    print RED, "\nFilter option not understood!\n\n", RESET;
    exit();
  }

  my $genes_found = @gene_ids;

  print GREEN, "\n\nFound $genes_found genes, excluding Pseudo-Genes\n", RESET;

  #Get gene info for each gene
  %genes = %{&getGeneInfo ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no")};
  my $gene_info_count = keys %genes;
  print GREEN, "\nFound gene info for $gene_info_count of these genes\n", RESET;

  #Create a second hash to Organize these genes according to their chromosome
  print GREEN, "\nCreating gene hash organized by chromosome\n", RESET;
  foreach my $gene_id (keys %genes){
    my $chr = $genes{$gene_id}{chromosome};

    if ($genes_chr{$chr}){
      my %gene_list = %{$genes_chr{$chr}{gene_list}};
      $gene_list{$gene_id}{chromosome} = $chr;
      $genes_chr{$chr}{gene_list} = \%gene_list;
    }else{
      my %gene_list;
      $gene_list{$gene_id}{chromosome} = $chr;
      $genes_chr{$chr}{gene_list} = \%gene_list;
    }
  }
  return();
}


#############################################################################################################################
#2.) Parse probe input data and map probesets to Ensembl genes                                                              #
############################################################################################################################
sub parseProbeData{
  my %args = @_;
  my $probe_file = $args{'-probe_file'};
  my $out_file = $args{'-out_file'};

  print GREEN, "\n\nParsing probe data and attempting to map probesets to Ensembl Genes\n", RESET;

  open (OUTFILE, ">$out_file") || die "\n\nCould not open output file: $out_file\n\n";
  print OUTFILE "probe_set_id\talexa_gene_id\tensembl_gene_id\ttranscript_cluster_id\tchromosome\tstrand\tstart\tstop\tprobe_count\n";

  open (PROBE, "$probe_file") || die "\n\nCould not open probe file: $probe_file\n\n";

  my $probe_record_count = 0;      #Number of probeset records examined
  my $probesets_mapped = 0;        #Number of probesets that map unambiguously to a single Ensembl gene id
  my $multiple_gene_mappings = 0;  #Number of probesets that map to multiple genes
  my $probesets_unmapped = 0;      #Number of probesets that did not map to any Ensembl gene id
  my $unknown_chromosomes = 0;     #Number of probesets with weird chromosome names

  while(<PROBE>){
    chomp($_);
    my $line = $_;
    $line =~ tr/"//d;

    #Remove the quote characters from each line
    my @line = split (",", $line);

    #Unless the first element of this line array contains a number, skip the line (header lines, etc.)
    unless ($line[0] =~ /^\d+/){
      next();
    }
    $probe_record_count++;

    #Grab the following data: probeset_id,transcript_cluster_id,seqname,strand,start,stop,probe_count
    my $probeset_id = $line[0];
    my $transcript_cluster_id = $line[3];
    my $chr = $line[4];
    my $strand = $line[5];
    my $start = $line[6];
    my $stop = $line[7];
    my $probe_count = $line[8];

    #Correct the format of the chromosome name
    my $chr_corrected;
    if ($chr =~ /^chr(.*)/){
      $chr_corrected = $1;
    }

    #Attempt to map this probeset to an Ensembl gene
    my $mapped_genes = 0; #Number of genes that the current probeset maps to
    my $target_gene_id;
    my @gene_ids;

    if ($genes_chr{$chr_corrected}){
      my %gene_list = %{$genes_chr{$chr_corrected}{gene_list}};

      foreach my $gene_id (keys %gene_list){
	#First make sure the probes are on the right strand.
	unless (($strand eq "+" && $genes{$gene_id}{chr_strand} eq "1") || ($strand eq "-" && $genes{$gene_id}{chr_strand} eq "-1")){
	  next();
	}

	if ($start >= $genes{$gene_id}{chr_start} && $start <= $genes{$gene_id}{chr_end} && $stop >= $genes{$gene_id}{chr_start} && $stop <= $genes{$gene_id}{chr_end}){
	  $target_gene_id = $gene_id;
	  push(@gene_ids, $gene_id);
	  $mapped_genes++;
	}
      }

      #Check if the probeset mapped to more than one gene
      if ($mapped_genes > 1){
	$multiple_gene_mappings++;
	print RED, "\nProbeset: $probeset_id maps to multiple ALEXA/Ensembl genes: @gene_ids", RESET;
      }elsif($mapped_genes == 1){
	$probesets_mapped++;

	#Use the one target_gene_id stored previously for this probeset
	print BLUE, "\n$probe_record_count: $probeset_id ($chr:$start-$stop $strand) maps to $genes{$target_gene_id}{ensembl_g_id} ($target_gene_id)", RESET;

	print OUTFILE "$probeset_id\t$target_gene_id\t$genes{$target_gene_id}{ensembl_g_id}\t$transcript_cluster_id\t$chr\t$strand\t$start\t$stop\t$probe_count\n";

      }else{
	#Probeset did not map to any gene
	$probesets_unmapped++;
      }

    }else{
      $unknown_chromosomes++;
      print YELLOW, "\nProbeset: $probeset_id has unknown chromosome name: $chr", RESET;
    }

    #DEBUG
    #if ($probesets_mapped == 100){
    #  last();
    #}

  }

  #Print summary statistics for the mapping of probesets
  print BLUE, "\n\nTotal probeset records processed: $probe_record_count", RESET;
  print BLUE, "\nProbesets successfully mapped to a single Ensembl gene: $probesets_mapped", RESET;
  print BLUE, "\nProbesets that mapped to multiple Ensembl genes: $multiple_gene_mappings", RESET;
  print BLUE, "\nProbesets that could not be mapped to any Ensembl gene: $probesets_unmapped", RESET;
  print BLUE, "\nProbesets that have weird chromosome names and could not be processed:  $unknown_chromosomes", RESET;
  close PROBE;
  close OUTFILE;

  return();
}
