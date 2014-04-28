#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009,2010 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;

my $analysis_dir = '';
my $shell = '';
my $gene_script = '';
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $annotation_dir = '';
my $junction_seq_size = '';
my $boundary_seq_size = '';
my $entrez_annotations = '';
my $paths_file = '';
my $ensembl_version = '';
my $ucsc_build = '';
my $species = '';
my $species_name = '';
my $project_name = '';
my $partition_file = '';
my $partition_file_version = '';
my $track_url = '';
my $temp_dir = '';
my $search_page_url = '';
my $alexa_home_path = '';
my $alexa_seq_path = '';
my $google_analytics_id = '';
my $bash_file = '';

GetOptions ('analysis_dir=s'=>\$analysis_dir,
            'shell=s'=>\$shell, 
            'gene_script=s'=>\$gene_script,
            'database=s'=>\$database,
            'server=s'=>\$server,
            'user=s'=>\$user,
            'password=s'=>\$password,
            'annotation_dir=s'=>\$annotation_dir,
            'junction_seq_size=i'=>\$junction_seq_size,
            'boundary_seq_size=i'=>\$boundary_seq_size,
            'entrez_annotations=s'=>\$entrez_annotations,
            'paths_file=s'=>\$paths_file,
            'ensembl_version=s'=>\$ensembl_version,
            'ucsc_build=s'=>\$ucsc_build,
            'species=s'=>\$species,
            'species_name=s'=>\$species_name,
            'project_name=s'=>\$project_name,
            'partition_file=s'=>\$partition_file,
            'partition_file_version=s'=>\$partition_file_version,
            'track_url=s'=>\$track_url,
            'temp_dir=s'=>\$temp_dir,
            'search_page_url=s'=>\$search_page_url,
            'alexa_home_path=s'=>\$alexa_home_path,
            'alexa_seq_path=s'=>\$alexa_seq_path,
            'google_analytics_id=s'=>\$google_analytics_id,
            'bash_file=s'=>\$bash_file);

unless ($analysis_dir && $shell && $gene_script && $database && $server && $user && $password && $annotation_dir && $junction_seq_size && $boundary_seq_size && $entrez_annotations && $paths_file && $ensembl_version && $ucsc_build && $species && $species_name && $project_name && $partition_file && $partition_file_version && $track_url && $temp_dir && $search_page_url && $alexa_home_path && $alexa_seq_path && $google_analytics_id && $bash_file){
  print RED, "\nRequired input parameter(s) missing (createGeneHtmlJobs.pl)\n\n", RESET;
  exit();
}

my %regions;
open(REGIONS, "$partition_file") || die "\n\nCould not open region file: $partition_file\n\n";
my $header = 1;
my $r = 0;
while(<REGIONS>){
  if ($header == 1){$header = 0; next();}
  $r++;
  chomp($_);
  my @line = split("\t", $_);
  $regions{$r}{chromosome} = $line[0];
  $regions{$r}{region} = $line[1];
  $regions{$r}{start_chr} = $line[2];
  $regions{$r}{end_chr} = $line[3];
  $regions{$r}{size} = $line[4];
  $regions{$r}{gene_count} = $line[5];
}
close(REGIONS);

open(BASH, ">$bash_file") || die "\n\nCould not open bash output file: $bash_file\n\n";
foreach my $r (sort {$a <=> $b} keys %regions){
  my $chromosome = $regions{$r}{chromosome};
  my $region = $regions{$r}{region};
  my $start_chr = $regions{$r}{start_chr};
  my $end_chr = $regions{$r}{end_chr};

  my $annotation_sub_dir = "$annotation_dir"."part/$partition_file_version/chr$chromosome"."_"."$region/";

  print BASH "source $shell; $gene_script  --database=$database  --server=$server  --user=$user  --password=$password  --annotation_dir=$annotation_sub_dir  --junction_seq_size=$junction_seq_size  --boundary_seq_size=$boundary_seq_size  --entrez_annotations=$entrez_annotations  --paths_file=$paths_file  --ensembl_version=$ensembl_version  --ucsc_build=$ucsc_build  --species=\"$species\"  --species_name=$species_name  --project_name=$project_name  --partition_file=$partition_file  --track_url=$track_url  --temp_dir=$temp_dir  --search_page_url=$search_page_url  --chr_filter=\"$chromosome:$region:$start_chr-$end_chr\"  --logfile=$analysis_dir/logs/website/$project_name/genes/generateGeneHtml_"."$chromosome"."_"."$region"."_LOG.txt  --alexa_home_path=$alexa_home_path  --alexa_seq_path=$alexa_seq_path  --google_analytics_id=$google_analytics_id\n";
}
close(BASH);

exit();

