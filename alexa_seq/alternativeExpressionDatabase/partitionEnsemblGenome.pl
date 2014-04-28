#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to divide the EnsEMBL annotated genome into pieces, each with a particular size or (number of genes)
#This must be done in way that prevents the region breakpoints from overlaping any genes
#To do this, gene content regions will first be created by merging overlaping genes

#Steps
#1.) Get all genes for all chromosomes
#2.) Merge overlaping genes on each chromosome into gene blocks (note the number of genes contained within each block)
#3.) Get the size of each chromosome from the EnsEMBL API
#4.) Go through each chromosome and divide it into pieces so that each piece contains approximately the same number of genes

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use Benchmark;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);

#Initialize command line options
my $ensembl_api_version = ''; #Version of EnsEMBL to use
my $species = '';
my $ensembl_database = '';
my $ensembl_server = '';
my $ensembl_user = '';
my $ensembl_password = '';
my $alexa_database = '';
my $alexa_server = '';
my $alexa_user = '';
my $alexa_password = '';
my $target_gene_count = '';
my $outfile = '';
my $chr_filter = '';

GetOptions ('ensembl_api_version=s'=>\$ensembl_api_version, 'species=s'=>\$species,
            'ensembl_database=s'=>\$ensembl_database, 'ensembl_server=s'=>\$ensembl_server, 'ensembl_user=s'=>\$ensembl_user, 'ensembl_password=s'=>\$ensembl_password, 
            'alexa_database=s'=>\$alexa_database, 'alexa_server=s'=>\$alexa_server, 'alexa_user=s'=>\$alexa_user, 'alexa_password=s'=>\$alexa_password, 
            'target_gene_count=i'=>\$target_gene_count, 'outfile=s'=>\$outfile, 'chr_filter=s'=>\$chr_filter);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the correct EnsEMBL API version using: --ensembl_api_version (41, 42, etc.)", RESET;
print GREEN, "\n\tSpecify the SOURCE Ensembl Database, Server, User and Password using: --ensembl_database  --ensembl_server  --ensembl_user  and  --ensembl_password", RESET;
print GREEN, "\n\tThe species using: --species (e.g. --species=Human or --species='Homo sapiens')", RESET;
print GREEN, "\n\t\tMake sure the species you supply matches the EnsEMBL database you supply!!", RESET;
print GREEN, "\n\tSpecify the TARGET Database and Server to query using: --alexa_database and --alexa_server", RESET;
print GREEN, "\n\tSpecify the User and Password for access using: --alexa_user and --alexa_password", RESET;
print GREEN, "\n\tSpecify the target number of genes per region defined using: --target_gene_count  (will not come out exactly because of overlaping genes!)", RESET;
print GREEN, "\n\tSpecify the output file using: --outfile", RESET;
print GREEN, "\n\tSpecify which chromosome to process using:  --chr_filter", RESET;

print GREEN, "\n\nExample: partitionEnsemblGenome.pl  --ensembl_api_version=49  --species=Human  --ensembl_database=homo_sapiens_core_49_36k  --ensembl_server=ensembl01.bcgsc.ca  --ensembl_user=ensembl  --ensembl_password=ensembl  --alexa_database=ALEXA_hs_49_36k  --alexa_server=jango.bcgsc.ca  --alexa_user=viewer  --alexa_password=viewer  --target_gene_count=100  --outfile=/projects/malachig/solexa/batch_jobs/temp/EnsEMBL_49_Regions_100genes_Y.txt  --chr_filter=Y\n\n", RESET;

#Make sure all options were specified
unless ($ensembl_api_version && $species && $ensembl_database && $ensembl_server && $ensembl_user && defined($ensembl_password) && $alexa_database && $alexa_server && $alexa_user && $alexa_password && $target_gene_count && $outfile && $chr_filter){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

&loadEnsemblApi('-api'=>$ensembl_api_version);


#1.) Get gene models from ALEXA for all genes
#2.) Merge overlaping genes on each chromosome into gene blocks (note the number of genes contained within each block)

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$alexa_database, '-server'=>$alexa_server, '-user'=>$alexa_user, '-password'=>$alexa_password);
my $genes_ref;

$| = 1; print BLUE, "\n\n1-2.) Getting basic gene info and merging overlaping genes for chr$chr_filter", RESET; $| = 0;
my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All', '-chromosome'=>$chr_filter)};
my $chr_genic_regions_ref = &getBasicGeneInfo ('-gene_ids'=> \@gene_ids);

#Close database connection
$alexa_dbh->disconnect();


#3.) Get the size of each chromosome from the EnsEMBL API
$| = 1; print BLUE, "\n\n3.) Getting size of this chromosome/supercontig", RESET; $| = 0;

##CONNECT TO ENSEMBL SERVER
my $registry = "Bio::EnsEMBL::Registry";
$registry->load_registry_from_db(-host =>$ensembl_server, -user =>$ensembl_user, -pass=>$ensembl_password, -db_version =>$ensembl_api_version);

#Get a connection to the local Ensembl CORE database
my $ensembl_core_api = Bio::EnsEMBL::Registry->get_DBAdaptor($species, "core");

my $legacy = 0;
unless($ensembl_core_api){
  print YELLOW, "\nLoad registry failed.  Try legacy DBAdaptor\n", RESET;
  $ensembl_core_api = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $ensembl_server,
    						         -user => $ensembl_user,
  						         -dbname => $ensembl_database,
						         -pass => $ensembl_password);
  $legacy=1;
}

#Get a slice for the entire chromosome
my $slice_adaptor;
if ($legacy){
  $slice_adaptor = $ensembl_core_api->get_SliceAdaptor(); #Alternate way of getting a 'slice adaptor'
}else{
  $slice_adaptor = $registry->get_adaptor($species, 'Core', 'Slice');
}

foreach my $chr (sort {$a cmp $b} keys %{$chr_genic_regions_ref}){
#  my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr);
#  unless($slice){
#    $slice = $slice_adaptor->fetch_by_region('supercontig', $chr);
#  }
  my $slice = $slice_adaptor->fetch_by_region('toplevel', $chr);
  my $chr_length = length($slice->seq());
  print YELLOW, "\nChr$chr = $chr_length bp", RESET;
  $chr_genic_regions_ref->{$chr}->{size} = $chr_length;
}

my $i_size = 100;

#open the output file
open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
print OUT "Chromosome\tRegion\tStart_chr\tEnd_chr\tSize\tGeneCount\n";

#4.) Go through each chromosome and divide it into pieces so that each piece contains approximately the same number of genes
$| = 1; print BLUE, "\n\n4.) Definining partitioned regions for this chromosome/supercontig\n", RESET; $| = 0;
foreach my $chr (sort {$a cmp $b} keys %{$chr_genic_regions_ref}){
  $| = 1; print BLUE, "\n\nChr: $chr ($chr_genic_regions_ref->{$chr}->{size} bp in size)\n", RESET; $| = 0;
  my $chr_start = 1;
  my $chr_end = $chr_genic_regions_ref->{$chr}->{size};

  my $block_count = 1;
  my $block_start = 1;
  my $block_end;

  my %regions;

  #Define the start of the first region
  $regions{$block_count}{start} = $block_start;

  my $gr_ref = $chr_genic_regions_ref->{$chr}->{genic_regions};
  my $gr_count = keys %{$gr_ref};
  my $grand_gr_counter = 0;
  my %contained_grs;
  my $c = 0;
  my $i;

  #Go through chromosome in chuncks of $i_size.  Note the number of genic regions encompassed.  Once the target number of genes has been encompassed, start a new block and continue
  for ($i = 1; $i <= $chr_end; $i+=$i_size){
    $c++;
    $block_end = $i;

    if ($c == 100){
      $c = 0;
      $| = 1; print CYAN, ".", RESET; $| = 0;
    }
    #print CYAN, "\n\t\tTEST: $block_start - $block_end", RESET;

    #Determine the number of genic regions contained within the current growing region
    my $overlap = 0;
    foreach my $gr (sort {$gr_ref->{$a}->{chr_start} <=> $gr_ref->{$b}->{chr_start}} keys %{$gr_ref}){
      my $gr_start = $gr_ref->{$gr}->{chr_start};
      my $gr_end = $gr_ref->{$gr}->{chr_end};

      if ($gr_start >= $block_start && $gr_start <= $block_end && $gr_end >= $block_start && $gr_end <= $block_end){
        #Completely contained within
        $contained_grs{$gr}{tmp} = '';

      }elsif(($gr_start >= $block_start && $gr_start <= $block_end) || ($gr_end >= $block_start && $gr_end <= $block_end) || ($gr_start < $block_start && $gr_end > $block_end)){
        #Not contained but overlaps at one end or the other or flanks completely
        $overlap = 1;
      }
    }

    #Watch our for cases where a genic region overlaps but is not completely contained within the current coords
    #We dont want blocks to end within a genic region so in these cases just proceed to the next iteration until we are out of it
    if ($overlap == 1){
      next();
    }

    #Now check if the desired number of contained genic regions has achieved or exceed the target
    my $contained_gr_count = keys %contained_grs;
    if ($contained_gr_count >= $target_gene_count){
      $regions{$block_count}{start} = $block_start;
      $regions{$block_count}{end} = $block_end;
      $regions{$block_count}{genic_region_count} = $contained_gr_count;

      print CYAN, "\nDefined region: chr$chr:$block_start-$block_end (containing $contained_gr_count genic regions)\n", RESET;

      $grand_gr_counter += $contained_gr_count;

      #Reset counters
      $block_count++;
      %contained_grs = ();
      $block_start = $block_end+1;
    }
  }

  #Define the end of the last region - unless the last region just happened to go to the end of the chr
  $regions{$block_count}{start} = $block_start;
  $regions{$block_count}{end} = $chr_end;
  $regions{$block_count}{genic_region_count} = keys %contained_grs;

  #Store the regions object for this chromosome
  $chr_genic_regions_ref->{$chr}->{regions} = \%regions;

  #Determine the number of individual genes that are actually contained within each region defined and print out a summary
  foreach my $r (sort {$a <=> $b} keys %regions){
    my $r_start = $regions{$r}{start};
    my $r_end = $regions{$r}{end};

    $regions{$r}{size} = ($r_end-$r_start)+1;

    my $contained_gene_count = 0;
    foreach my $g (sort {$genes_ref->{$a}->{chr_start} <=> $genes_ref->{$b}->{chr_start}} keys %{$genes_ref}){
      my $g_start = $genes_ref->{$g}->{chr_start};
      my $g_end = $genes_ref->{$g}->{chr_end};
      my $gene_chr = $genes_ref->{$g}->{chromosome};
      unless ($gene_chr eq $chr){
        next();
      }

      if ($g_start >= $r_start && $g_start <= $r_end && $g_end >= $r_start && $g_end <= $r_end){
        $contained_gene_count++;
      }
    }
    $regions{$r}{gene_count} = $contained_gene_count;
    $| = 1; print MAGENTA, "\n\tRegion: $r\tStart: $r_start\tEnd: $r_end\tSize: $regions{$r}{size}\tGenicRegionCount: $regions{$r}{genic_region_count}\tGeneCount: $regions{$r}{gene_count}", RESET; $| = 0;
    print OUT "$chr\t$r\t$r_start\t$r_end\t$regions{$r}{size}\t$regions{$r}{gene_count}\n";
  }
}
print "\n\nSCRIPT COMPLETE\n\n";

exit();


################################################################################################
#Get gene models, masked gene sequence, intron content, etc.                                   #
################################################################################################
sub getBasicGeneInfo{
  my %args = @_;
  my @gene_ids = @{$args{'-gene_ids'}};

  #Get basic info about each gene
  $| = 1; print BLUE, "\n\nGetting gene models for chr$chr_filter", RESET; $| = 0;

  #DEBUG
  $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");
  #@gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All', '-chromosome'=>22)};
  #$genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

  my $gene_count = keys %{$genes_ref};
  $| = 1; print BLUE, "\n\nFound $gene_count genes in $alexa_database", RESET; $| = 0; 

  #Identify all genic regions by merging overlaping genes (using there gene coordinates)
  $| = 1; print BLUE, "\n\nUsing chromosome coordinates to determine chromosome-wide genic regions\n", RESET; $| = 0; 
  my %chr_genic_regions;
  my $gr_count = 0;

  #First make sure the chr_start is smaller than chr_end
  foreach my $g (keys %{$genes_ref}){
    my $tmp = $genes_ref->{$g}->{chr_start};
    if ($genes_ref->{$g}->{chr_start} > $genes_ref->{$g}->{chr_end}){
      $genes_ref->{$g}->{chr_start} = $genes_ref->{$g}->{chr_end};
      $genes_ref->{$g}->{chr_end} = $tmp;
    }
  }

  my $count = 0;
  foreach my $g (sort {$genes_ref->{$a}->{chr_start} <=> $genes_ref->{$b}->{chr_start}} keys %{$genes_ref}){
    my $g_start = $genes_ref->{$g}->{chr_start};
    my $g_end = $genes_ref->{$g}->{chr_end};
    my $chr = $genes_ref->{$g}->{chromosome};

    my $genic_regions_ref;
    if ($chr_genic_regions{$chr}){
      $genic_regions_ref = $chr_genic_regions{$chr}{genic_regions};
    }else{
      my %tmp;
      $chr_genic_regions{$chr}{genic_regions} = \%tmp;
      $genic_regions_ref = $chr_genic_regions{$chr}{genic_regions};
    }

    $count++;
    if ($count == 100){
      $count = 0;
      $| = 1; print BLUE, ".", RESET;  $| = 0;
    }

    #For this gene, go through all genic regions thus far and look for overlaps
    my $overlap = 0;
    foreach my $gr (sort {$genic_regions_ref->{$a}->{chr_start} <=> $genic_regions_ref->{$b}->{chr_start}} keys %{$genic_regions_ref}){
      my $gr_start = $genic_regions_ref->{$gr}->{chr_start};
      my $gr_end = $genic_regions_ref->{$gr}->{chr_end};

      #Test overlap
      if (($g_start >= $gr_start && $g_start <= $gr_end) || ($g_end >= $gr_start && $g_end <= $gr_end) || ($g_start <= $gr_start && $g_end >= $gr_end)){
        $overlap = 1;
      }

      if ($overlap == 1){
        #Merge this gene into an existing genic region and get the adjusted coordinates
        my @coords = ($g_start, $g_end, $gr_start, $gr_end);
        my @coords_sort = sort {$a <=> $b} @coords;
        $genic_regions_ref->{$gr_count}->{chr_start} = $coords_sort[0];
        $genic_regions_ref->{$gr_count}->{chr_end} = $coords_sort[3];
        last();
      }
    }

    if ($overlap == 0){
      #Enter a new genic region
      $gr_count++;
      $genic_regions_ref->{$gr_count}->{chr_start} = $g_start;
      $genic_regions_ref->{$gr_count}->{chr_end} = $g_end;
    }
  }

  my $chr_count = keys %chr_genic_regions;
  $| = 1; print BLUE, "\n\nMerged these genes (from $chr_count chromosomes) into $gr_count genic_regions\n", RESET; $| = 0; 

  return(\%chr_genic_regions);
}
