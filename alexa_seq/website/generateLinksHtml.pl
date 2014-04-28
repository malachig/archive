#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Generate index.htm files so that every single gene page is linked from somewhere.
#The main links page can be linked from the gene search page itself

#links page should described the chromosome and coordinates covered followed by a list of EnsEMBL gene IDs which link to individual gene pages
#The main links page will just be links to individual links pages for each genomic region


use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use website::WEB qw(:all);
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $annotation_dir = '';
my $junction_seq_size = '';
my $boundary_seq_size = '';
my $partition_file = '';
my $genes_dir = '';
my $search_page_url = '';
my $alexa_home_path = '';
my $alexa_seq_path = '';
my $google_analytics_id = '';

GetOptions('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'annotation_dir=s'=>\$annotation_dir, 
           'junction_seq_size=i'=>\$junction_seq_size, 'boundary_seq_size=i'=>\$boundary_seq_size,
           'partition_file=s'=>\$partition_file, 'genes_dir=s'=>\$genes_dir,
           'search_page_url=s'=>\$search_page_url, 'alexa_home_path=s'=>\$alexa_home_path, 'alexa_seq_path=s'=>\$alexa_seq_path, 'google_analytics_id=s'=>\$google_analytics_id);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script generates HTML links pages to all gene pages", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the base path to annotation files using:  --annotation_dir", RESET;
print GREEN, "\n\tSpecify the junction database sequence length using:  --junction_seq_size", RESET;
print GREEN, "\n\tSpecify the boundary database sequence length using:  --boundary_seq_size", RESET;
print GREEN, "\n\tSpecify a file containing the coordinates of genome partitions used for the analysis using: --partition_file", RESET;
print GREEN, "\n\tSpecify the top directory where gene links pages will be written using:  --genes_dir", RESET;
print GREEN, "\n\tSpecify the URL to your Xapian-Omega search page using:  --search_page_url", RESET;
print GREEN, "\n\tSpecify the URL to your ALEXA home page using: --alexa_home_path", RESET;
print GREEN, "\n\tSpecify the URL to your ALEXA-Seq results page using: --alexa_seq_path", RESET;
print GREEN, "\n\tSpecify your Google Analytics ID using: --google_analytics_id", RESET;

print GREEN, "\n\nUsage: generateLinksHtml.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --annotation_dir=/projects/malachig/sequence_databases/hs_53_36o/  --junction_seq_size=62  --boundary_seq_size=62  --partition_file=/projects/malachig/solexa/batch_jobs/EnsEMBL_53_Regions_50genes.txt  --genes_dir=/projects/malachig/solexa/temp/genes/  --search_page_url=http://www.bcgsc.ca/xapian-search/omega  --alexa_home_path=http://www.alexaplatform.org/index.htm  --alexa_seq_path=http://www.alexaplatform.org/alexa_seq/results.htm  --google_analytics_id=UA-xxxxxx-x\n\n", RESET;

#Check user supplied options
unless ($database && $server && $user && $password && $annotation_dir && $junction_seq_size && $boundary_seq_size && $partition_file && $genes_dir && $search_page_url && $alexa_home_path && $alexa_seq_path && $google_analytics_id){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}
$genes_dir = &checkDir('-dir'=>$genes_dir, '-clear'=>"no");


#1.) Get all genes
my @types_list = qw (Gene); 
my %types_list;
foreach my $type (@types_list){
  $types_list{$type}=1;
}
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
my $a_ref;
my @gene_id_list = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};
$a_ref = &getAnnotations('-annotation_dir'=>$annotation_dir, '-dbh'=>$alexa_dbh, '-database'=>$database, '-types_list'=>\%types_list, '-junction_seq_size'=>$junction_seq_size, '-boundary_seq_size'=>$boundary_seq_size);
$alexa_dbh->disconnect();

my $genes_ref = $a_ref->{'Gene'};


#2.) Import genome partitions from the specified file
my $partitions_ref = &getPartitions('-file'=>$partition_file);

#3.) Map genes to genome partitions
my $assigned_partitions = &assignPartitions('-a_ref'=>$a_ref, '-partitions_ref'=>$partitions_ref);
print BLUE, "\n\t\tAssigned $assigned_partitions records to a single genome partition", RESET;

#4.) Determine the number of partitions comprising each region
foreach my $chr (sort {$a cmp $b} keys %{$partitions_ref}){
  my $regions_ref = $partitions_ref->{$chr}->{partitions};
  $partitions_ref->{$chr}->{count} = keys %{$regions_ref};
  $partitions_ref->{$chr}->{name} = $chr;
}

#5.) Go through each partition and generate a master links page
#print Dumper $partitions_ref;
print BLUE, "\n\nGenerating master links page:\n\n", RESET;
#foreach my $chr (sort {{$a cmp $b} || {$partitions_ref->{$a}->{start} <=> $partitions_ref->{$b}->{start}}} keys %{$partitions_ref}){
my $ml_content = '';

foreach my $chr (sort {$partitions_ref->{$b}->{count} <=> $partitions_ref->{$a}->{count} || $partitions_ref->{$a}->{name} cmp $partitions_ref->{$b}->{name}} keys %{$partitions_ref}){
  my $regions_ref = $partitions_ref->{$chr}->{partitions};

  my $region_count = keys %{$regions_ref};

  $ml_content = "$ml_content"."<P CLASS=\"Indented12LR_s16_Bold\">Chromosome: chr$chr is divided into $region_count regions</P>\n";

  foreach my $part (sort {$regions_ref->{$a}->{start} <=> $regions_ref->{$b}->{start}} keys %{$regions_ref}){

    my $region = $regions_ref->{$part}->{region};
    my $start = $regions_ref->{$part}->{start};
    my $end = $regions_ref->{$part}->{end};
    my $gene_count =  $regions_ref->{$part}->{gene_count};
    my $size = $regions_ref->{$part}->{size};
    print BLUE, "\n\tchr$chr"."_"."$region:$start-$end ($size bp containing $gene_count genes)", RESET;
    my $region_name = "chr$chr"."_"."$region";
    my $region_desc = "$start-$end ($size bp containing $gene_count genes)";

    $ml_content = "$ml_content"."\t<P CLASS=\"Indented12LR_s16\"><A HREF=\"$region_name/index.html\">$region_name</A>: $region_desc</P>\n";
  }
  $ml_content = "$ml_content"."<BR>\n\n";
}
#Write out the page
my $web_path = "$genes_dir"."index.html";
my $title = "Links to regions of the genome, each containing links to 50-75 gene pages (ordered by number of genes)";
my $meta_description = "Provides a list of chromosomes broken down into sub-regions, each containing 50-75 genes";
my $meta_keywords = "Chromosomes, Gene Regions, Gene Lists, Chromosome Coordinates";
&writePage('-path'=>$web_path, '-title'=>$title, '-content'=>\$ml_content, '-css_path'=>"../../ALEXA2.css", '-alexa_home_path'=>"$alexa_home_path", '-alexa_seq_home_path'=>"$alexa_seq_path", '-summary_path'=>"../Summary.htm", '-genes_path'=>"index.html", '-search_path'=>$search_page_url, '-meta_description'=>$meta_description, '-meta_keywords'=>$meta_keywords, '-google_analytics'=>0, '-google_analytics_id'=>$google_analytics_id, '-div_count'=>12);


#6.) Make sure the chr_start is always smaller than chr_end for each gene
foreach my $gene_id (keys %{$genes_ref}){
  my $chr_start = $genes_ref->{$gene_id}->{chr_start};
  my $chr_end = $genes_ref->{$gene_id}->{chr_end};
  if ($chr_start < $chr_end){
    my $tmp = $chr_start;
    $chr_start = $chr_end;
    $chr_end = $tmp;
    $genes_ref->{$gene_id}->{chr_start} = $chr_start;
    $genes_ref->{$gene_id}->{chr_end} = $chr_end;
  }
}

#7.) Create a links page for each gene neighbourhood
print BLUE, "\n\nWriting pages for each genomic region", RESET;
foreach my $chr (sort {$partitions_ref->{$b}->{count} <=> $partitions_ref->{$a}->{count} || $partitions_ref->{$a}->{name} cmp $partitions_ref->{$b}->{name}} keys %{$partitions_ref}){
  my $regions_ref = $partitions_ref->{$chr}->{partitions};
  my $region_count = keys %{$regions_ref};

  foreach my $part (sort {$regions_ref->{$a}->{start} <=> $regions_ref->{$b}->{start}} keys %{$regions_ref}){
    my $gene_counter = 0;
    my $region = $regions_ref->{$part}->{region};
    my $start = $regions_ref->{$part}->{start};
    my $end = $regions_ref->{$part}->{end};
    my $gene_count =  $regions_ref->{$part}->{gene_count};
    my $size = $regions_ref->{$part}->{size};

    my $region_name = "chr$chr"."_"."$region";
    my $region_desc = "$start-$end ($size bp containing $gene_count genes)";
 
    my $gl_content = '';
    $gl_content = "$gl_content"."<P CLASS=\"Indented12LR_s16_Bold\">$region_name: $region_desc</P>\n";

    #Now add records for each gene that belongs to this region
    my $current_partition = "chr"."$chr"."_"."$region";
    foreach my $gene_id (sort {$genes_ref->{$a}->{chr_start} <=> $genes_ref->{$b}->{chr_start}} keys %{$genes_ref}){
      my $partition = $genes_ref->{$gene_id}->{partition};

      unless($partition eq $current_partition){
        next();
      }
      my $ensembl_g_id = $genes_ref->{$gene_id}->{ensembl_g_id};

      #Make sure the gene page is actually there - dont create broken links
      my $gene_page_path = "$genes_dir"."$current_partition"."/$ensembl_g_id".".htm";
      unless(-e $gene_page_path){
        print YELLOW, "\n\t\tCould not find gene page: $gene_page_path - skipping", RESET;
        next();
      }
      $gene_counter++;

      my $chr_start = $genes_ref->{$gene_id}->{chr_start};
      my $chr_end = $genes_ref->{$gene_id}->{chr_end};
      my $chr_strand = $genes_ref->{$gene_id}->{chr_strand};
      my $gene_name = $genes_ref->{$gene_id}->{gene_name};
      my $link = "$ensembl_g_id".".htm";
        
      $gl_content = "$gl_content"."\t<P CLASS=\"Indented12LR_s16\"><A HREF=\"$link\">$ensembl_g_id</A> a.k.a. '$gene_name' at chr$chr:$chr_start-$chr_end ($chr_strand)</P>\n";

    }
    $gl_content = "$gl_content"."<BR>\n\n";

    #Write out this page
    #print Dumper $gl_content;
    my $web_path = "$genes_dir"."$current_partition"."/index.html";
    print BLUE, "\n\tWriting page: $web_path ($current_partition)", RESET;
    my $title = "Links to individual gene pages for $gene_counter genes within the genomic region: $current_partition:$start-$end";
    my $meta_description = "Provides links to single gene pages for the genomic region: $current_partition:$start-$end which contains $gene_counter gene records";
    my $meta_keywords = "Gene list, $current_partition:$start-$end";
    &writePage('-path'=>$web_path, '-title'=>$title, '-content'=>\$gl_content, '-css_path'=>"../../../ALEXA2.css", '-alexa_home_path'=>"$alexa_home_path", '-alexa_seq_home_path'=>"$alexa_seq_path", '-summary_path'=>"../../Summary.htm", '-genes_path'=>"../index.html", '-search_path'=>$search_page_url, '-meta_description'=>$meta_description, '-meta_keywords'=>$meta_keywords, '-google_analytics'=>0, '-google_analytics_id'=>$google_analytics_id, '-div_count'=>12);
  }
}

exit();















