#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Essentially a wrapper for Martin K's mapability tool
#Uses this tool to calculate 36-mer read mapability for every non-redundant exonic position for the user specified version of EnsEMBL (ALEXA database)

#NOTE: this script is hardcoded to call Martin K's script and assumes you want to use hg18 and 36-mers

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
my $outfile = '';
my $chr_limit = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'chr_limit=s'=>\$chr_limit,
	    'outfile=s'=>\$outfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the outfile summarizing coverage after each read block using: --outfile", RESET;
print GREEN, "\n\tTo process only a single chromosome, include the --chr_limit option (e.g. 1 or MT or X)", RESET;
print GREEN, "\n\nExample: calculateExonicMapability.pl  --database=ALEXA_hs_49_36k  --server=jango.bcgsc.ca  --user=malachig  --password=pwd  --chr_limit='1'  --outfile=/projects/malachig/solexa/mapability/ENST_v49_36mer_mapability_Chr1.txt\n\n", RESET;

unless ($database && $server && $user && $password && $outfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}


#1.) First get all the neccessary gene info required to perform the analysis
my $genes_ref;
my $gene_transcripts_ref;
my $gene_exon_content_ref;
my $chr_genes_ref;

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
my @gene_ids;
if ($chr_limit){
  @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All', '-chromosome'=>$chr_limit)};
}else{
  @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};
}
#my @gene_ids = (1..100);

my $total_gene_count = scalar(@gene_ids);
print BLUE, "\nWorking with a total of $total_gene_count genes\n\n", RESET;

&getBasicGeneInfo('-gene_list'=>\@gene_ids);

#Close database connection
$alexa_dbh->disconnect();

$| = 1; print BLUE, "\n\n2.) Determining chromosome coordinates for the non-redundant exonic base coverage for ALL EnsEMBL genes\n"; $| = 0;
&getExonContentCoords('-gene_list'=>\@gene_ids);

#print Dumper $gene_exon_content_ref;

#Now calculate the mappability for every exonic position of every gene
#Note that there are overlapping positions from exons of adjacent genes on a single chromosome
#Resolve this after the fact by storing all positions for a chromosome and the corresponding mappability

my $nmermappability_bin = "/home/martink/cvs/nmermappability/current/nmermappability ";

#Example command:
#/home/martink/cvs/nmermappability/current/nmermappability -genome hg18 -hash 36 -nmer 36 -chr 1 -start 10000090 -end 10000100

#my $nmermappability_base_cmd = "$nmermappability_bin"."-genome hg18 -hash 36 -nmer 36 ";
my $nmermappability_base_cmd = "$nmermappability_bin"."-genome hg18 -hash 42 -nmer 42 ";

$| = 1; print BLUE, "\n\nNow calculating chromosome coordinates for each chromosome (skipping unassembled contigs!)", RESET; $| = 0;

foreach my $chr (keys %{$chr_genes_ref}){
  my %mapability;

  $| = 1; print BLUE, "\n\tProcessing chromosome: $chr\n\t", RESET; $| = 0;

  unless ($chr =~ /^\d+$|^Y$|^X$|^M$/){
    print YELLOW, "\n\tSkipping chromosome: $chr", RESET;
    next();
  }


  my @chr_genes = @{$chr_genes_ref->{$chr}->{gene_list}};

  foreach my $gene_id (@chr_genes){
    $| = 1; print BLUE, ".", RESET; $| = 0;
    my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};

    foreach my $exon_id (sort {$a <=> $b} keys %{$exon_content_ref}){

      #Exon content coordinates should already be 0 indexed?
      my $start = ($exon_content_ref->{$exon_id}->{chr_start});
      my $end = ($exon_content_ref->{$exon_id}->{chr_end});

      my $cmd = "$nmermappability_base_cmd"."-chr $chr -start $start -end $end";

      my $result = `$cmd`;

      my @result = split ("\n", $result);

      foreach my $line (@result){
	chomp($line);
	my @vals = split (" ", $line);
	my $pos = $vals[0];
	my $map = $vals[6];

	$mapability{$pos} = $map;
	
      }
    }
  }

  my $gene_count = scalar(@chr_genes);
  my $base_count = keys %mapability;
  $| = 1; print BLUE, "\n\tCalculated the mapability for $base_count bases corresponding to $gene_count genes on chr $chr\n", RESET; $| = 0;

  #Now print the mapability values to an output file
  open (OUT, ">$outfile") || die "\nCould not open output file: $outfile\n\n", RESET;
  foreach my $pos (sort {$a <=> $b} keys %mapability){
    print OUT "$chr\t$pos\t$mapability{$pos}\n";
  }
  close OUT;
}


exit();


############################################################################################################################################
#Get basic info for all genes from the user specified ALEXA database                                                                       #
############################################################################################################################################
sub getBasicGeneInfo{
  my %args = @_;
  my @gene_list = @{$args{'-gene_list'}};

  my %chr_genes;

  #Get the gene info for all genes for which reads were found on the current chromosome
  $| = 1; print BLUE, "\n1-a.) Getting gene data", RESET; $| = 0;
  $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_list, '-sequence'=>"no");

  #Get the transcript info for all transcripts of these genes
  $| = 1; print BLUE, "\n1-b.) Getting transcript data as well as exons for each transcript", RESET; $| = 0;
  $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_list, '-sequence'=>"no");

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
    }
  }

  $chr_genes_ref = \%chr_genes;

  #Get exon content for all genes
  $| = 1; print BLUE, "\n\n1-d.) Getting EXON CONTENT of each gene", RESET; $| = 0;
  $gene_exon_content_ref = &getExonContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_list);

  return();
}



############################################################################################################################################
#Initialize hashes to store Exon coverage for a subset of genes
############################################################################################################################################
sub getExonContentCoords{
  my %args = @_;
  my @gene_list = @{$args{'-gene_list'}};

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
	$exon_content_ref->{$exon_id}->{chr_start} = $query_chr_start;
	$exon_content_ref->{$exon_id}->{chr_end} = $query_chr_end;
	$gene_exon_content_ref->{$gene_id}->{exon_content_size} += ($query_chr_end - $query_chr_start)+1;

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

	$exon_content_ref->{$exon_id}->{chr_start} = $query_chr_start;
	$exon_content_ref->{$exon_id}->{chr_end} = $query_chr_end;
	$gene_exon_content_ref->{$gene_id}->{exon_content_size} += ($query_chr_end - $query_chr_start)+1;

      }else{
	print RED, "\nStrand format: $chr_strand not understood !\n\n", RESET;
	exit();
      }
    }
  }

  return();
}

