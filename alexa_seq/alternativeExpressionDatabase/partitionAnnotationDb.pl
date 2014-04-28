#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009,2010 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to divide all ALEXA-Seq annotation DB data into partitions according to genome region
#Based on an input file defining the boundaries of genome regions, a new set of annotation db files will be created region-by-region
#Use of these files in downstream analyses that are also being performed region-by-region will improve runtime and file I/O loads

#Example complete annotation dir:
#/projects/malachig/sequence_databases/hs_53_36o/

#Make a partitioned version (assumining '250' regions file was used)
#/projects/malachig/sequence_databases/hs_53_36o/part/250/
# - this dir will contain subdirs for each region
# - each region subdir will look just like the complete annotation db with the following 'types': genes/ transcripts/ exonRegions/ exonJunctions/ exonBoundaries/ introns/ intergenics/

#Steps.
#1.) Create the 'part' dir.
#2.) Create the 'part/250/' dir.  If it already exists, ask the user if they want to create it cleanly
#3.) Load the regions to be created from an input file
#4.) Create subdirs for each region
#5.) Create 'type' subdirs within each region subdir
#6.) For each 'type' load the main annotation file.  Write the header to all region specific files for that type
#7.) Go through the main annotation file.  Figure out the region corresponding to a particular entry and write it to the file
#    - If a feature has any overlap at all with the target region, write it to the file (this may result in some intergenic regions being written to two files)...
#8.) Compress all results files
#9.) Check that all every region subdir has all the expected type subdirs and files.


use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);
use website::WEB qw(:all);


my $region_file_version = '';
my $annotation_dir = '';

GetOptions ('region_file_version=i'=>\$region_file_version, 'annotation_dir=s'=>\$annotation_dir);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the genome region file version using: --region_file_version (e.g. 250)", RESET;
print GREEN, "\n\tSpecify the target annotation dir using:  --annotation_dir", RESET;

print GREEN, "\n\nExample: partitionAnnotationDir.pl  --region_file_version=250  --annotation_dir=/projects/malachig/sequence_databases/hs_53_36o/\n\n", RESET;

unless ($region_file_version && $annotation_dir){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Define the type names
my @types = qw (genes transcripts exonRegions exonJunctions exonBoundaries introns intergenics);

#Check the specific annotation dir
$annotation_dir = &checkDir('-dir'=>$annotation_dir, '-clear'=>"no");

#Check for the existence of the corresponding region file
my $regions_file = "$annotation_dir"."Regions_"."$region_file_version"."_Genes.txt";
unless(-e $regions_file){
  print RED, "\n\nCould not find a region file corresponding to the specified version ($region_file_version): $regions_file\n\n", RESET;
  exit();
}

#1.) Create the 'part' dir.
my $part_dir = "$annotation_dir"."part/";
if (-e $part_dir && -d $part_dir){
  print BLUE, "\n\nMain partition dir already exists", RESET;
  $part_dir = &checkDir('-dir'=>$part_dir, '-clear'=>"no");
}else{
  print BLUE, "\n\nMain partition dir not found - creating", RESET;
  mkdir($part_dir);
  $part_dir = &checkDir('-dir'=>$part_dir, '-clear'=>"no");
}

#2.) Create the 'part/250/' dir.  If it already exists, ask the user if they want to create it cleanly
my $part_regions_dir = "$part_dir"."$region_file_version/";
if (-e $part_regions_dir && -d $part_regions_dir){
  print BLUE, "\n\nPartition regions dir already exists", RESET;
  $part_regions_dir = &checkDir('-dir'=>$part_regions_dir, '-clear'=>"yes", '-recursive'=>"yes");
}else{
  print BLUE, "\n\nPartition regions dir not found - creating", RESET;
  mkdir($part_regions_dir);
  $part_regions_dir = &checkDir('-dir'=>$part_regions_dir, '-clear'=>"no");
}

#3.) Load the regions to be created from an input file
my $regions_ref = &getPartitions('-file'=>$regions_file);
my $region_count = keys %{$regions_ref};


#4.) Create subdirs for each region and store the path to the newly created dirs associated with each region
print BLUE, "\n\nCreating the subdirs for each region", RESET;
foreach my $chr (keys %{$regions_ref}){
  my $part_ref = $regions_ref->{$chr}->{partitions};
  foreach my $r (keys %{$part_ref}){
    my $region = $part_ref->{$r}->{region};
    my $dir = "$part_regions_dir"."chr$chr"."_"."$region/";
    mkdir($dir);
    $part_ref->{$r}->{dir} = $dir;
  }
}

#5.) Create 'type' subdirs within each region subdir
print BLUE, "\n\nCreating the subdirs for each type within each region subdir", RESET;
foreach my $chr (keys %{$regions_ref}){
  my $part_ref = $regions_ref->{$chr}->{partitions};
  foreach my $r (keys %{$part_ref}){
    my $region = $part_ref->{$r}->{region};
    my $r_dir = $part_ref->{$r}->{dir};
    foreach my $type (@types){
      my $dir = "$r_dir"."$type/";
      mkdir($dir);
    }
  }
}

#Get a list of all annotation files to be processed (if there are multiple junction/boundary files) get these now.
my %a_files;
my $j_test = "ls $annotation_dir"."exonJunctions/exonJunctions_*mers_annotated.txt.gz";
my $b_test = "ls $annotation_dir"."exonBoundaries/exonBoundaries_*mers_annotated.txt.gz";
my @j_files = `$j_test`;
my @b_files = `$b_test`;
my $c = 0;
$c++; $a_files{$c}{name} = "genes_annotated.txt"; $a_files{$c}{dir} = "genes";
$c++; $a_files{$c}{name} = "transcripts_annotated.txt"; $a_files{$c}{dir} = "transcripts";
$c++; $a_files{$c}{name} = "exonRegions_annotated.txt"; $a_files{$c}{dir} = "exonRegions";
foreach my $jfile (@j_files){
  chomp($jfile);
  if ($jfile =~ /(exonJunctions_\d+mers_annotated\.txt)\.gz/){
    $c++; $a_files{$c}{name} = "$1"; $a_files{$c}{dir} = "exonJunctions";
  }
}
foreach my $bfile (@b_files){
  chomp($bfile);
  if ($bfile =~ /(exonBoundaries_\d+mers_annotated\.txt)\.gz/){
    $c++; $a_files{$c}{name} = "$1"; $a_files{$c}{dir} = "exonBoundaries";
  }
}
$c++; $a_files{$c}{name} = "introns_annotated.txt"; $a_files{$c}{dir} = "introns";
$c++; $a_files{$c}{name} = "activeIntronRegions.txt"; $a_files{$c}{dir} = "introns";
$c++; $a_files{$c}{name} = "silentIntronRegions.txt"; $a_files{$c}{dir} = "introns";
$c++; $a_files{$c}{name} = "intergenics_annotated.txt"; $a_files{$c}{dir} = "intergenics";
$c++; $a_files{$c}{name} = "activeIntergenicRegions.txt"; $a_files{$c}{dir} = "intergenics";
$c++; $a_files{$c}{name} = "silentIntergenicRegions.txt"; $a_files{$c}{dir} = "intergenics";


#6.) For each annotation file.  Write the header to all region specific files for that type
#    - Note that for exon junctions and boundaries there may be multiple files for multiple database lengths...
#    - At the same time determine and store the identity of the columns containing the chromosome, chr_start, and chr_end values
print BLUE, "\n\nChecking the annotation files, getting column positions...", RESET;
foreach my $c (sort {$a <=> $b} keys %a_files){
  my $a_file_name = $a_files{$c}{name};
  my $name_gz = "$a_file_name".".gz";
  my $a_dir_name = $a_files{$c}{dir};
  my $a_path = "$annotation_dir"."$a_dir_name/$a_file_name";
  my $a_path_gz = "$annotation_dir"."$a_dir_name/$a_file_name".".gz";
  $a_files{$c}{name_gz} = "$name_gz";
  $a_files{$c}{path} = "$a_path";
  $a_files{$c}{path_gz} = "$a_path_gz";

  print BLUE, "\n\t$a_path_gz", RESET;
  open (ANN, "zcat $a_path_gz | head -n 1 |") || die "\n\nCould not open annotation file: $a_path_gz\n\n";
  while (<ANN>){
    chomp($_);
    my $header = $_;
    $a_files{$c}{header} = $_;
    my @cols = split("\t", $_);
    my $col_pos = 0;
    foreach my $col (@cols){
      if ($col =~ /^chromosome$/i){
        $a_files{$c}{chr_col} = $col_pos;
      }
      if ($col =~ /^Unit1\_start\_chr$/i){
        $a_files{$c}{start_col} = $col_pos;
      }
      if ($col =~ /^Unit2\_end\_chr$/i){
        $a_files{$c}{end_col} = $col_pos;
      }elsif($col =~ /Unit1\_end\_chr$/i){
        $a_files{$c}{end_col} = $col_pos;
      }
      $col_pos++;
    }
  }
  close(ANN);
  #Make sure all required columns were found
  unless (defined($a_files{$c}{chr_col}) && defined($a_files{$c}{start_col}) && defined($a_files{$c}{end_col})){
    print RED, "\nMissing required column.\tchr_col = $a_files{$c}{chr_col}\tstart_col = $a_files{$c}{start_col}\tend_col = $a_files{$c}{end_col}", RESET;
    exit();
  }
}

#Now initialize all files with the header
print BLUE, "\n\nCreating the subdirs for each type within each region subdir", RESET;
foreach my $chr (keys %{$regions_ref}){
  my $part_ref = $regions_ref->{$chr}->{partitions};
  foreach my $r (keys %{$part_ref}){
    my $region = $part_ref->{$r}->{region};
    my $r_dir = $part_ref->{$r}->{dir};
    foreach my $c (sort {$a <=> $b} keys %a_files){
      my $header = $a_files{$c}{header};
      my $typedir = $a_files{$c}{dir};
      my $name = $a_files{$c}{name};

      my $out = "$r_dir"."$typedir/$name";
      open (OUT, ">$out") || die "\n\nCould not open out file: $out\n\n";
      print OUT "$header";
      close(OUT);
    }
  }
}


#7.) Go through each main annotation file.  Figure out the region corresponding to a particular entry and write it to the appropriate file
#    - If a feature has any overlap at all with the target region, write it to the file (this may result in some intergenic regions being written to two files)...
print BLUE, "\n\nParsing annotation records by genome region...", RESET;
foreach my $c (sort {$a <=> $b} keys %a_files){
  my $counter = 0;
  my $typedir = $a_files{$c}{dir};
  my $name = $a_files{$c}{name};
  my $a_path_gz = $a_files{$c}{path_gz};
  my $chr_col = $a_files{$c}{chr_col};
  my $start_col = $a_files{$c}{start_col};
  my $end_col = $a_files{$c}{end_col};

  print BLUE, "\n\t$a_path_gz  ", RESET;

  #Open a filehandle for each output file
  foreach my $chr (keys %{$regions_ref}){
    my $part_ref = $regions_ref->{$chr}->{partitions};
    foreach my $r (keys %{$part_ref}){
      my $region = $part_ref->{$r}->{region};
      my $r_dir = $part_ref->{$r}->{dir};
      my $out = "$r_dir"."$typedir/$name";
      my $fh = IO::File->new(">>$out") || die "\nCould not open output file handle: $out\n\n";
      $part_ref->{$r}->{fh} = $fh;
    }
  }

  open (ANN, "zcat $a_path_gz |") || die "\n\nCould not open annotation file: $a_path_gz\n\n";
  my $header = 1;
  ANNO:while(<ANN>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header == 1){
      $header = 0;
      next();
    }
    $counter++;
    if ($counter == 10000){
      $| = 1; print BLUE, ".", RESET; $| = 0;
      $counter = 0;
    }
    my $chr = $line[$chr_col];
    if ($chr eq "M"){$chr = "MT";}
    my $start = $line[$start_col];
    my $end = $line[$end_col];

    unless ($regions_ref->{$chr}){
      print RED, "\n\nCould not find a matching chromosome for $chr in the Regions object\n\n", RESET;
      exit();
    }
    my $part_ref = $regions_ref->{$chr}->{partitions};
    foreach my $r (keys %{$part_ref}){
      my $region = $part_ref->{$r}->{region};
      my $r_dir = $part_ref->{$r}->{dir};
      my $r_start = $part_ref->{$r}->{start};
      my $r_end = $part_ref->{$r}->{end};

      if ($typedir eq "intergenics"){
        #Look for overlap
        if (($start >= $r_start && $start <= $r_end) || ($end >= $r_start && $end <= $r_end)){
          #my $out = "$r_dir"."$typedir/$name";
          #open (OUT, ">>$out") || die "\n\nCould not open out file: $out\n\n";
          #print OUT "\n$_";
          #close(OUT);
          my $fh = $part_ref->{$r}->{fh};
          print $fh "\n$_";
        }
      }else{
        if (($start >= $r_start && $start <= $r_end) || ($end >= $r_start && $end <= $r_end)){
          #my $out = "$r_dir"."$typedir/$name";
          #open (OUT, ">>$out") || die "\n\nCould not open out file: $out\n\n";
          #print OUT "\n$_";
          #close(OUT);
          my $fh = $part_ref->{$r}->{fh};
          print $fh "\n$_";
          next(ANNO);
        }
      }
    }
  }
  close(ANN);

  #Close the filehandles
  foreach my $chr (keys %{$regions_ref}){
    my $part_ref = $regions_ref->{$chr}->{partitions};
    foreach my $r (keys %{$part_ref}){
      my $region = $part_ref->{$r}->{region};
      my $r_dir = $part_ref->{$r}->{dir};
      my $out = "$r_dir"."$typedir/$name";
      my $fh = $part_ref->{$r}->{fh};
      $fh->close;
    }
  }


  #8.) Compress all results files for the current type
  foreach my $chr (sort keys %{$regions_ref}){
    my $part_ref = $regions_ref->{$chr}->{partitions};
    foreach my $r (sort keys %{$part_ref}){
      my $region = $part_ref->{$r}->{region};
      my $r_dir = $part_ref->{$r}->{dir};
      my $out = "$r_dir"."$typedir/$name";
      my $gzip_cmd = "gzip -f $out";
      print BLUE, "\n\t\t$gzip_cmd", RESET;
      system($gzip_cmd);
    }
  }
}

#9.) Check that every region subdir has all the expected type subdirs and files.
print BLUE, "\n\nChecking for presence of all partioned files ...", RESET;
my $missing_files = 0;
foreach my $c (sort {$a <=> $b} keys %a_files){
  my $typedir = $a_files{$c}{dir};
  my $name = $a_files{$c}{name};

  foreach my $chr (keys %{$regions_ref}){
    my $part_ref = $regions_ref->{$chr}->{partitions};
    foreach my $r (keys %{$part_ref}){
      my $region = $part_ref->{$r}->{region};
      my $r_dir = $part_ref->{$r}->{dir};
      my $out = "$r_dir"."$typedir/$name".".gz";
      unless (-e $out){
        $missing_files++;
        print RED, "\n\tCould not find file: $out (missing file count: $missing_files)", RESET;
      }
    }
  }
}

if ($missing_files == 0){
  print BLUE, "\n\nLooks good.  Should be safe to proceed ...", RESET;
}else{
  print RED, "\n\nMissing file count was greater than 0!  - Do not proceed.  (Rerun job)\n\n", RESET;
  exit();
}

print "\n\n";

exit();


