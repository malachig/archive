#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Create batch job commands for calculating expression estimates for:
#TYPES="Repeats ENST_v53 Junctions_v53 Boundaries_v53 Introns_v53 Intergenics_v53 Transcripts_v53"

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
use utilities::utility qw(:all);

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $ensembl_version = '';
my $ucsc_build = '';
my $lib_names_file = '';
my $script_dir = '';
my $batch_dir = '';
my $analysis_dir = '';
my $annotation_dir = '';
my $regions_file_version = '';
my $ucsc_dir = '';
my $web_path = '';
my $cluster_commands = '';
my $junction_size = '';
my $boundary_size = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
            'ensembl_version=s'=>\$ensembl_version, 'ucsc_build=s'=>\$ucsc_build, 'lib_names_file=s'=>\$lib_names_file, 'script_dir=s'=>\$script_dir, 
            'batch_dir=s'=>\$batch_dir, 'analysis_dir=s'=>\$analysis_dir, 'annotation_dir=s'=>\$annotation_dir, 'regions_file_version=i'=>\$regions_file_version, 'ucsc_dir=s'=>\$ucsc_dir, 'web_path=s'=>\$web_path,
            'cluster_commands=s'=>\$cluster_commands, 'junction_size=i'=>\$junction_size, 'boundary_size=i'=>\$boundary_size);

unless ($database && $server && $user && $password && $ensembl_version && $ucsc_build && $lib_names_file && $script_dir && $batch_dir && $analysis_dir && $annotation_dir && $regions_file_version && $ucsc_dir && $web_path && $junction_size && $boundary_size){
  print GREEN, "\n\nExample usage:\n\ncreateExpressionJobs.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --ensembl_version=53  --ucsc_build=hg18  --lib_names_file=/projects/malachig/solexa/batch_jobs/Neuroblastoma/Neuroblastoma_Lib_Names.txt  --script_dir=/home/malachig/svn/solexa_analysis/  --batch_dir=/projects/malachig/solexa/batch_jobs/Neuroblastoma/  --analysis_dir=/projects/malachig/solexa/  --annotation_dir=/projects/malachig/sequence_databases/hs_53_36o/  --regions_file_version=50  --ucsc_dir=/home/malachig/www/public/htdocs/solexa  --web_path=http://www.bcgsc.ca/people/malachig/htdocs/solexa  --junction_size=62  --boundary_size=62\n\n", RESET;
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

$script_dir = &checkDir('-dir'=>$script_dir, '-clear'=>"no");
$batch_dir = &checkDir('-dir'=>$batch_dir, '-clear'=>"no");
$analysis_dir = &checkDir('-dir'=>$analysis_dir, '-clear'=>"no");
$annotation_dir = &checkDir('-dir'=>$annotation_dir, '-clear'=>"no");
$ucsc_dir = &checkDir('-dir'=>$ucsc_dir, '-clear'=>"no");


#Get the repeats database
my $rep_db = "$annotation_dir"."repeats/*.txt";
my $result = `ls $rep_db`;
chomp($result);
my $repeats_db;
if ($result =~ /(.*\.txt)/){
  $repeats_db = $1;
}else{
  print RED, "\nCould not find repeats database in: $rep_db\n\n", RESET;
  exit();
}


#Get list of library IDs and corresponding names from input file
#print "\n\n#Getting library names ...", RESET;
open(LIBS, "$lib_names_file") || die "\nCould not open lib names file: $lib_names_file\n\n";
my $l = 0;
my %libs;
while(<LIBS>){
  $l++;
  chomp($_);
  my @data = split("\t", $_);
  if ($_ =~ /\w+/){
    $libs{$l}{id} = $data[0];
    $libs{$l}{name} = $data[1];
    $libs{$l}{color_set} = $data[2];
    $libs{$l}{correction_factor} = $data[3];
  }
}
close(LIBS);

#Get list of genomic regions to be processed
#print "\n\n#Getting regions ...", RESET;
my $regions_file = "$annotation_dir"."Regions_"."$regions_file_version"."_Genes.txt";

open(REGIONS, "$regions_file") || die "\nCould not open regions file: $regions_file\n\n";
my $r = 0;
my %regions;
my $header = 1;
while(<REGIONS>){
  if ($header == 1){$header = 0; next();}
  $r++;
  chomp($_);
  my @data = split("\t", $_);
  $regions{$r}{chromosome} = $data[0];
  $regions{$r}{region} = $data[1];
  $regions{$r}{chr_start} = $data[2];
  $regions{$r}{chr_end} = $data[3];
  $regions{$r}{size} = $data[4];
  $regions{$r}{gene_count} = $data[5];
}
close(REGIONS);

#Create batch files containing expression calculation commands and print out the commands needed to use them
print "\n\n#Jobs to be submitted as follows:";


#1.) REPEATS
my $outfile = "$batch_dir"."generateExpressionValues_Repeats.sh";
open (OUT, ">$outfile") || die "\nCould not open output file: $outfile\n\n";
foreach my $l (sort {$a <=> $b} keys %libs){
  my $lib_id = $libs{$l}{id};
  my $lib_name = $libs{$l}{name};
  
  my $cmd = "$script_dir"."expression/generateExpressionValues_Repeats.pl  --repeat_library_file=$repeats_db  --library=$lib_id  --read_record_dir=$analysis_dir"."read_records/$lib_id/  --mapped_reads_dir=$analysis_dir"."read_records/$lib_id/Repeats/  --working_dir=$analysis_dir"."read_records/$lib_id/Repeats/temp/  --min_bit_score=60.0  --results_dir=$analysis_dir"."read_records/$lib_id/Repeats/Summary/  --cutoffs_file=$analysis_dir"."figures_and_stats/$lib_id/Expression_v$ensembl_version/$lib_id"."_NORM1_average_coverage_cutoffs.txt  --log_file=$analysis_dir"."logs/$lib_id/generateExpressionValues/Repeats/generateExpressionValues_Repeats_LOG.txt\n";

  print OUT "source ~/.bashrc; $cmd";
}
close(OUT);
print "\n\n#---REPEATS---#";
if ($cluster_commands){
  print "\ncd $batch_dir\nrm -fr GEV_R\nmqsub  --file $outfile  --name GEV_R  --mkdir  --delay 5";
}else{
  print "\nbash $outfile";
}

#2.) GENES AND EXON REGIONS
$outfile = "$batch_dir"."generateExpressionValues_ENST_v"."$ensembl_version".".sh";
open (OUT, ">$outfile") || die "\nCould not open output file: $outfile\n\n";
foreach my $l (sort {$a <=> $b} keys %libs){
  my $lib_id = $libs{$l}{id};
  my $lib_name = $libs{$l}{name};
  my $color_set = $libs{$l}{color_set};
  foreach my $r (sort {$a <=> $b} keys %regions){
    my $chr_filter = "$regions{$r}{chromosome}".":"."$regions{$r}{region}".":"."$regions{$r}{chr_start}-$regions{$r}{chr_end}";
    my $chr_region = "$regions{$r}{chromosome}"."_"."$regions{$r}{region}";

    my $chromosome = $regions{$r}{chromosome};
    my $region = $regions{$r}{region};
    my $annotation_sub_dir = "$annotation_dir"."part/$regions_file_version/chr$chromosome"."_"."$region/";

    my $cmd = "$script_dir"."expression/generateExpressionValues_GeneExon.pl  --database=$database  --server=$server  --user=$user  --password=$password  --library=$lib_id  --library_name=$lib_name  --chr_filter=\"$chr_filter\"  --working_dir=$analysis_dir"."read_records/$lib_id/ENST_v$ensembl_version/Summary/temp/  --min_bit_score=60.0  --min_seq_coverage=75.0  --mapped_reads_dir=$analysis_dir"."read_records/$lib_id/ENST_v$ensembl_version/part/$chr_region/  --read_record_dir=$analysis_dir"."read_records/$lib_id/  --annotation_dir=$annotation_sub_dir   --results_dir=$analysis_dir"."read_records/$lib_id/ENST_v$ensembl_version/Summary/results/  --ucsc_dir=$ucsc_dir"."$lib_id/  --ucsc_build=$ucsc_build  --web_path=$web_path/$lib_id/  --color_set=$color_set  --cutoffs_file=$analysis_dir"."figures_and_stats/$lib_id/Expression_v"."$ensembl_version"."/$lib_id"."_NORM1_average_coverage_cutoffs.txt  --log_file=$analysis_dir"."logs/$lib_id/generateExpressionValues/ENST_v$ensembl_version/generateExpressionValues_GeneExon_chr"."$chr_region".".txt\n";

   print OUT "source ~/.bashrc; $cmd";
  }
}
close(OUT);
print "\n\n#---GENES AND EXON REGIONS---#";
if ($cluster_commands){
  print "\ncd $batch_dir\nrm -fr GEV_GE\nmqsub  --file $outfile  --name GEV_GE  --mkdir  --delay 5";
}else{
  print "\nbash $outfile";
}


#3.) EXON JUNCTIONS
$outfile = "$batch_dir"."generateExpressionValues_Junctions_v"."$ensembl_version".".sh";
open (OUT, ">$outfile") || die "\nCould not open output file: $outfile\n\n";
foreach my $l (sort {$a <=> $b} keys %libs){
  my $lib_id = $libs{$l}{id};
  my $lib_name = $libs{$l}{name};
  my $color_set = $libs{$l}{color_set};
  my $correction_factor = $libs{$l}{correction_factor};
  foreach my $r (sort {$a <=> $b} keys %regions){
    my $chr_filter = "$regions{$r}{chromosome}".":"."$regions{$r}{region}".":"."$regions{$r}{chr_start}-$regions{$r}{chr_end}";
    my $chr_region = "$regions{$r}{chromosome}"."_"."$regions{$r}{region}";
    my $chromosome = $regions{$r}{chromosome};
    my $region = $regions{$r}{region};
    my $annotation_sub_dir = "$annotation_dir"."part/$regions_file_version/chr$chromosome"."_"."$region/";

    my $cmd = "$script_dir"."expression/generateExpressionValues_JunctionBoundary.pl  --seq_type=junction  --database=$database  --server=$server  --user=viewer  --password=viewer  --library=$lib_id  --library_name=$lib_name  --chr_filter=\"$chr_filter\"  --mapped_reads_dir=$analysis_dir"."read_records/$lib_id/Junctions_v$ensembl_version/part/$chr_region/  --read_record_dir=$analysis_dir"."read_records/$lib_id/  --seq_database=$annotation_sub_dir/exonJunctions/exonJunctions_"."$junction_size"."mers_annotated.txt.gz  --working_dir=$analysis_dir"."read_records/$lib_id/Junctions_v$ensembl_version/Summary/temp/  --min_bit_score=69.9  --min_seq_coverage=70.0  --correction_factor=$correction_factor  --results_dir=$analysis_dir"."read_records/$lib_id/Junctions_v$ensembl_version/Summary/results/  --ucsc_dir=$ucsc_dir"."$lib_id/  --ucsc_build=$ucsc_build  --web_path=$web_path/$lib_id/  --color_set=$color_set  --cutoffs_file=$analysis_dir"."figures_and_stats/$lib_id/Expression_v$ensembl_version/$lib_id"."_NORM1_average_coverage_cutoffs.txt  --log_file=$analysis_dir"."logs/$lib_id/generateExpressionValues/Junctions_v$ensembl_version/generateExpressionValues_Junction_chr"."$chr_region".".txt\n";

    print OUT "source ~/.bashrc; $cmd";

  }
}
close(OUT);
print "\n\n#---EXON JUNCTIONS---#";
if ($cluster_commands){
  print "\ncd $batch_dir\nrm -fr GEV_J\nmqsub  --file $outfile  --name GEV_J  --mkdir  --delay 5";
}else{
  print "\nbash $outfile";
}

#4.) EXON BOUNDARIES
$outfile = "$batch_dir"."generateExpressionValues_Boundaries_v"."$ensembl_version".".sh";
open (OUT, ">$outfile") || die "\nCould not open output file: $outfile\n\n";
foreach my $l (sort {$a <=> $b} keys %libs){
  my $lib_id = $libs{$l}{id};
  my $lib_name = $libs{$l}{name};
  my $color_set = $libs{$l}{color_set};
  my $correction_factor = $libs{$l}{correction_factor};

  foreach my $r (sort {$a <=> $b} keys %regions){
    my $chr_filter = "$regions{$r}{chromosome}".":"."$regions{$r}{region}".":"."$regions{$r}{chr_start}-$regions{$r}{chr_end}";
    my $chr_region = "$regions{$r}{chromosome}"."_"."$regions{$r}{region}";
    my $chromosome = $regions{$r}{chromosome};
    my $region = $regions{$r}{region};
    my $annotation_sub_dir = "$annotation_dir"."part/$regions_file_version/chr$chromosome"."_"."$region/";

    my $cmd = "$script_dir"."expression/generateExpressionValues_JunctionBoundary.pl  --seq_type=boundary  --database=$database  --server=$server  --user=viewer  --password=viewer  --library=$lib_id  --library_name=$lib_name  --chr_filter=\"$chr_filter\"  --mapped_reads_dir=$analysis_dir"."read_records/$lib_id/Boundaries_v$ensembl_version/part/$chr_region/  --read_record_dir=$analysis_dir"."read_records/$lib_id/  --seq_database=$annotation_sub_dir/exonBoundaries/exonBoundaries_"."$boundary_size"."mers_annotated.txt.gz  --working_dir=$analysis_dir"."read_records/$lib_id/Boundaries_v$ensembl_version/Summary/temp/  --min_bit_score=69.9  --min_seq_coverage=75.0  --correction_factor=$correction_factor  --results_dir=$analysis_dir"."read_records/$lib_id/Boundaries_v$ensembl_version/Summary/results/  --ucsc_dir=$ucsc_dir"."$lib_id/  --ucsc_build=$ucsc_build  --web_path=$web_path/$lib_id/  --color_set=$color_set  --cutoffs_file=$analysis_dir"."figures_and_stats/$lib_id/Expression_v$ensembl_version/$lib_id"."_NORM1_average_coverage_cutoffs.txt  --log_file=$analysis_dir"."logs/$lib_id/generateExpressionValues/Boundaries_v$ensembl_version/generateExpressionValues_Boundary_chr"."$chr_region".".txt\n";

    print OUT "source ~/.bashrc; $cmd";
  }
}
close(OUT);
print "\n\n#---EXON BOUNDARIES---#";
if ($cluster_commands){
  print "\ncd $batch_dir\nrm -fr GEV_B\nmqsub  --file $outfile  --name GEV_B  --mkdir  --delay 5";
}else{
  print "\nbash $outfile";
}

#5.) INTRONS
$outfile = "$batch_dir"."generateExpressionValues_Introns_v"."$ensembl_version".".sh";
open (OUT, ">$outfile") || die "\nCould not open output file: $outfile\n\n";
foreach my $l (sort {$a <=> $b} keys %libs){
  my $lib_id = $libs{$l}{id};
  my $lib_name = $libs{$l}{name};
  my $color_set = $libs{$l}{color_set};
  foreach my $r (sort {$a <=> $b} keys %regions){
    my $chr_filter = "$regions{$r}{chromosome}".":"."$regions{$r}{region}".":"."$regions{$r}{chr_start}-$regions{$r}{chr_end}";
    my $chr_region = "$regions{$r}{chromosome}"."_"."$regions{$r}{region}";
    my $chromosome = $regions{$r}{chromosome};
    my $region = $regions{$r}{region};
    my $annotation_sub_dir = "$annotation_dir"."part/$regions_file_version/chr$chromosome"."_"."$region/";

    my $cmd = "$script_dir"."expression/generateExpressionValues_Intron.pl  --database=$database  --server=$server  --user=$user  --password=$password  --library=$lib_id  --library_name=$lib_name  --chr_filter=\"$chr_filter\"  --working_dir=$analysis_dir"."read_records/$lib_id/Introns_v$ensembl_version/Summary/temp/  --min_bit_score=60.0  --min_seq_coverage=75.0  --mapped_reads_dir=$analysis_dir"."read_records/$lib_id/Introns_v$ensembl_version/part/$chr_region/  --read_record_dir=$analysis_dir"."read_records/$lib_id/  --annotation_dir=$annotation_sub_dir"."introns/   --results_dir=$analysis_dir"."read_records/$lib_id/Introns_v$ensembl_version/Summary/results/  --ucsc_dir=$ucsc_dir"."$lib_id/  --ucsc_build=$ucsc_build  --web_path=$web_path/$lib_id/  --color_set=$color_set  --cutoffs_file=$analysis_dir"."figures_and_stats/$lib_id/Expression_v$ensembl_version/$lib_id"."_NORM1_average_coverage_cutoffs.txt  --log_file=$analysis_dir"."logs/$lib_id/generateExpressionValues/Introns_v$ensembl_version/generateExpressionValues_Intron_chr$chr_region".".txt\n";

    print OUT "source ~/.bashrc; $cmd";
  }
}
close(OUT);
print "\n\n#---INTRONS---#";
if ($cluster_commands){
  print "\ncd $batch_dir\nrm -fr GEV_I\nmqsub  --file $outfile  --name GEV_I  --mkdir  --delay 5";
}else{
  print "\nbash $outfile";
}

#6.) INTERGENIC REGIONS
$outfile = "$batch_dir"."generateExpressionValues_Intergenics_v"."$ensembl_version".".sh";
open (OUT, ">$outfile") || die "\nCould not open output file: $outfile\n\n";
foreach my $l (sort {$a <=> $b} keys %libs){
  my $lib_id = $libs{$l}{id};
  my $lib_name = $libs{$l}{name};
  my $color_set = $libs{$l}{color_set};
  foreach my $r (sort {$a <=> $b} keys %regions){
    my $chr_filter = "$regions{$r}{chromosome}".":"."$regions{$r}{region}".":"."$regions{$r}{chr_start}-$regions{$r}{chr_end}";
    my $chr_region = "$regions{$r}{chromosome}"."_"."$regions{$r}{region}";
    my $chromosome = $regions{$r}{chromosome};
    my $region = $regions{$r}{region};
    my $annotation_sub_dir = "$annotation_dir"."part/$regions_file_version/chr$chromosome"."_"."$region/";

    my $cmd = "$script_dir"."expression/generateExpressionValues_Intergenic.pl  --database=$database  --server=$server  --user=$user  --password=$password  --library=$lib_id  --library_name=$lib_name  --chr_filter=\"$chr_filter\"  --working_dir=$analysis_dir"."read_records/$lib_id/Intergenics_v$ensembl_version/Summary/temp/  --min_bit_score=60.0  --min_seq_coverage=75.0  --mapped_reads_dir=$analysis_dir"."read_records/$lib_id/Intergenics_v$ensembl_version/part/$chr_region/  --read_record_dir=$analysis_dir"."read_records/$lib_id/  --annotation_dir=$annotation_sub_dir"."intergenics/  --results_dir=$analysis_dir"."read_records/$lib_id/Intergenics_v$ensembl_version/Summary/results/  --ucsc_dir=$ucsc_dir"."$lib_id/  --ucsc_build=$ucsc_build  --web_path=$web_path/$lib_id/  --color_set=$color_set  --cutoffs_file=$analysis_dir"."figures_and_stats/$lib_id/Expression_v$ensembl_version/$lib_id"."_NORM1_average_coverage_cutoffs.txt  --log_file=$analysis_dir"."logs/$lib_id/generateExpressionValues/Intergenics_v$ensembl_version/generateExpressionValues_Intergenic_chr$chr_region".".txt\n";

    print OUT "source ~/.bashrc; $cmd";
  }
}
close(OUT);
print "\n\n#---INTERGENIC REGIONS---#";
if ($cluster_commands){
  print "\ncd $batch_dir\nrm -fr GEV_IG\nmqsub  --file $outfile  --name GEV_IG  --mkdir  --delay 5";
}else{
  print "\nbash $outfile";
}
print "\n\n";

exit();

