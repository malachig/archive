#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to fix read names that were incorrectly generated and are missing the lane number
#READ IDS SHOULD HAVE THE FORMAT: FLOWCELL_LANE_TILE_X_Y

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

my $library_id = '';
my $analysis_dir = '';
my $ensembl_version = '';

GetOptions ('library_id=s'=>\$library_id, 'analysis_dir=s'=>\$analysis_dir, 'ensembl_version=i'=>\$ensembl_version);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the main analysis directory using: --analysis_dir", RESET;
print GREEN, "\n\tSpecify the project name for this project using:  --library_id", RESET;
print GREEN, "\n\tSpecify the ensembl version used for this analysis using: --ensembl_version", RESET;
print GREEN, "\n\nExample:  checkParsingJobs.pl  --library_id=HS1441  --analysis_dir=/projects/malachig/alexa_seq/  --ensembl_version=53\n\n", RESET;

unless ($library_id && $analysis_dir && $ensembl_version){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}
$analysis_dir = &checkDir('-dir'=>$analysis_dir, '-clear'=>"no");
my $rr_dir = "$analysis_dir"."read_records/$library_id/";
$rr_dir = &checkDir('-dir'=>$rr_dir, '-clear'=>"no");
my $br_dir = "$analysis_dir"."blast_results/$library_id/";
$br_dir = &checkDir('-dir'=>$br_dir, '-clear'=>"no");


#Define directory types
my %mr_types;
$mr_types{'1'}{name} = "Repeat";
$mr_types{'1'}{dir} = "Repeats";
$mr_types{'1'}{file} = "Repeats";
$mr_types{'1'}{columns_expected} = 3;
$mr_types{'2'}{name} = "GeneExon";
$mr_types{'2'}{dir} = "ENST_v"."$ensembl_version";
$mr_types{'2'}{file} = "ENST_v"."$ensembl_version";
$mr_types{'2'}{columns_expected} = 3;
$mr_types{'3'}{name} = "Junction";
$mr_types{'3'}{dir} = "Junctions_v"."$ensembl_version";
$mr_types{'3'}{file} = "Junctions_v"."$ensembl_version";
$mr_types{'3'}{columns_expected} = 1;
$mr_types{'4'}{name} = "Boundary";
$mr_types{'4'}{dir} = "Boundaries_v"."$ensembl_version";
$mr_types{'4'}{file} = "Boundaries_v"."$ensembl_version";
$mr_types{'4'}{columns_expected} = 1;
$mr_types{'5'}{name} = "Intron";
$mr_types{'5'}{dir} = "Introns_v"."$ensembl_version";
$mr_types{'5'}{file} = "Introns";
$mr_types{'5'}{columns_expected} = 3;
$mr_types{'6'}{name} = "Intergenic";
$mr_types{'6'}{dir} = "Intergenics_v"."$ensembl_version";
$mr_types{'6'}{file} = "Intergenics";
$mr_types{'6'}{columns_expected} = 3;

my %br_types;
$br_types{'1'}{dir} = "repeats";
$br_types{'2'}{dir} = "transcripts";
$br_types{'3'}{dir} = "exonJunctions";
$br_types{'4'}{dir} = "exonBoundaries";
$br_types{'5'}{dir} = "introns";
$br_types{'6'}{dir} = "intergenics";


#1.) Get the read record files
print BLUE, "\n\nGet read record files:", RESET;
my %rr_files;
opendir(DIRHANDLE, "$rr_dir") || die "\nCannot open directory: $rr_dir\n\n";
my @files = readdir(DIRHANDLE);

#Skip directories within the specified directory
my $rr = 0;
foreach my $file (sort @files){
  my $file_path = "$rr_dir"."$file";
  if (-e $file_path && -d $file_path){
    next();
  }
  unless($file =~ /\.txt\.gz$/){
    next();
  }
  my $base_name = '';
  if ($file =~ /(.*)\.txt\.gz$/){
    $base_name = $1;
  }

  #Get the lane number
  my $lane = '';
  if ($file =~ /Lane(\d+)/){
    $lane = $1;
  }else{
    print RED, "\n\nCould not get lane number from file name: $file", RESET;
    exit();
  }
  $rr++;
  print BLUE, "\n\t$rr: $base_name ($file_path) (LANE: $lane)", RESET;;
  $rr_files{$rr}{dir} = $rr_dir;
  $rr_files{$rr}{name} = $file;
  $rr_files{$rr}{base_name} = $base_name;
  $rr_files{$rr}{path} = $file_path;
  $rr_files{$rr}{lane} = $lane;
  $rr_files{$rr}{invalid_read_ids} = 0;
}
closedir(DIRHANDLE);

#2.) Get the mapping results files
print BLUE, "\n\nGet the mapping results files:", RESET;
my %mr_files;
my $mr = 0;
foreach my $mrt (sort {$a <=> $b} keys %mr_types){
  my $mr_name = $mr_types{$mrt}{name};
  my $mr_dir = $mr_types{$mrt}{dir};
  my $mr_file = $mr_types{$mrt}{file};
  my $cols_expected = $mr_types{$mrt}{columns_expected};
  
  my $path = "$rr_dir"."$mr_dir/";
  print BLUE, "\n\t$mrt: $path", RESET;

  opendir(DIRHANDLE, "$path") || die "\nCannot open directory: $path\n\n";
  my @files = readdir(DIRHANDLE);

  foreach my $file (sort @files){
    my $file_path = "$path"."$file";
    if (-e $file_path && -d $file_path){
      next();
    }
    unless($file =~ /\.txt\.gz$/){
      next();
    }
    my $base_name = '';
    if ($file =~ /(.*)\.txt\.gz$/){
      $base_name = $1;
    }

    #Get the lane number
    my $lane = '';
    if ($file =~ /Lane(\d+)/){
      $lane = $1;
    }else{
      print RED, "\n\nCould not get lane number from file name: $file", RESET;
      exit();
    }
    $mr++;
    print BLUE, "\n\t\t$mr: $base_name ($file_path) (LANE: $lane)", RESET;
    $mr_files{$mr}{dir} = $path;
    $mr_files{$mr}{name} = $file;
    $mr_files{$mr}{base_name} = $base_name;
    $mr_files{$mr}{path} = $file_path;
    $mr_files{$mr}{lane} = $lane;
    $mr_files{$mr}{invalid_read_ids} = 0;
    $mr_files{$mr}{columns_expected} = $cols_expected;
  }
}

#3.) Get the Blast results files
print BLUE, "\n\nGet the blast results files:", RESET;
my %br_dirs;
my %br_files;
my $brd = 0;
my $br = 0;

my $path = "$br_dir";
print BLUE, "\n\t$path", RESET;
opendir(DIRHANDLE, "$path") || die "\nCannot open directory: $path\n\n";
my @dirs = readdir(DIRHANDLE);
foreach my $dir (sort @dirs){
  my $dir_path = "$path"."$dir/";
  unless($dir =~ /Lane\d+$/){
    next();
  }
   
  #Get the lane number
  my $lane = '';
  if ($dir =~ /Lane(\d+)/){
    $lane = $1;
  }else{
    print RED, "\n\nCould not get lane number from dir name: $dir", RESET;
    exit();
  }
  $brd++;
  print BLUE, "\n\t\t$brd: $dir ($dir_path) (LANE: $lane)", RESET;
  $br_dirs{$brd}{name} = $dir;
  $br_dirs{$brd}{path} = $dir_path;
  $br_dirs{$brd}{lane} = $lane;

  #For this dir, get the files in each sub dir
  foreach my $brt (sort {$a <=> $b} keys %br_types){
    my $brt_dir = $br_types{$brt}{dir};
    my $brt_dir_path = "$dir_path"."$brt_dir/";
    print BLUE, "\n\t\t\t$brt: $brt_dir ($brt_dir_path)", RESET;

    opendir(DIRHANDLE, "$brt_dir_path") || die "\nCannot open directory: $brt_dir_path\n\n";
    my @files = readdir(DIRHANDLE);
    foreach my $file (sort @files){
      my $file_path = "$brt_dir_path"."$file";
      if (-e $file_path && -d $file_path){
        next();
      }
      unless($file =~ /\.gz$/){
        next();
      }
      my $base_name = '';
      if ($file =~ /(.*)\.gz$/){
        $base_name = $1;
      }

      $br++;
      #print BLUE, "\n\t\t\t\t$br: $base_name ($file_path) (LANE: $lane)", RESET;
      $br_files{$br}{dir} = $brt_dir_path;
      $br_files{$br}{name} = $file;
      $br_files{$br}{base_name} = $base_name;
      $br_files{$br}{path} = $file_path;
      $br_files{$br}{lane} = $lane;
      $br_files{$br}{invalid_read_ids} = 0;
    }
    closedir(DIRHANDLE);
  }
}
closedir(DIRHANDLE);

#4.) Summarize all the files found
my $rr_count = keys %rr_files;
my $mr_count = keys %mr_files;
my $br_count = keys %br_files;

print BLUE, "\n\nFound $rr_count ReadRecord files, $mr_count MappingResults files, and $br_count BlastResults files", RESET;

#5.) Check each file to see whether the read IDs need to be repaired
print BLUE, "\n\nChecking all files to see if they need to be fixed\n", RESET;
my $files_to_fix = 0;
foreach my $rr (sort {$a <=> $b} keys %rr_files){
  my $path = $rr_files{$rr}{path};
  my $lane = $rr_files{$rr}{lane};
  my $name = $rr_files{$rr}{name};
  open (TEMP, "zcat $path |") || die "\n\nCould not open file: $path\n\n";
  my $header = 1;
  my %columns;
  while(<TEMP>){
    chomp($_);
    my @line = split("\t", $_);
    my $cc = 0;
    if ($header == 1){
      $header = 0;
      foreach my $col (@line){
        $cc++;
        $columns{$col}{col_pos} = $cc;
      }
      my $rid_col = $columns{'Read_ID'}{col_pos};
      my $rid1_col = $columns{'Read1_ID'}{col_pos};
      my $rid2_col = $columns{'Read2_ID'}{col_pos};
      unless($rid_col && $rid1_col && $rid2_col){
        print RED, "\n\nCould not find neccessary read ID columns\n\n", RESET;
        exit();
      }
      next();
    }
    #Check for invalid read id
    if ($line[0] =~ /(\w+)\_(\d+)\_(\d+)\_(\d+)\_(\d+)/){
      unless($2 == $lane){
        $rr_files{$rr}{invalid_read_ids} = 1;
        $files_to_fix++;
        $| = 0; print BLUE, ".", RESET; $| = 1;
        last();
      }
    }else{
      print RED, "\n\nFormat of read id not understood: $line[0]\n\n", RESET;
      exit();
    }
  }
  close(TEMP);
}
foreach my $mr (sort {$a <=> $b} keys %mr_files){
  my $path = $mr_files{$mr}{path};
  my $lane = $mr_files{$mr}{lane};
  my $name = $mr_files{$mr}{name};
  my $cols_expected = $mr_files{$mr}{columns_expected};

  open (TEMP, "zcat $path |") || die "\n\nCould not open file: $path\n\n";
  my $header = 1;
  my %columns;
  while(<TEMP>){
    chomp($_);
    my @line = split("\t", $_);
    my $cc = 0;
    if ($header == 1){
      $header = 0;
      foreach my $col (@line){
        $cc++;
        $columns{$col}{col_pos} = $cc;
      }
      $mr_files{$mr}{columns} = \%columns;
      my $rid_col = $columns{'Read_ID'}{col_pos};
      my $rid1_col = $columns{'R1_ID'}{col_pos};
      my $rid2_col = $columns{'R2_ID'}{col_pos};
      my $cols_found = 0;
      if ($rid_col){$cols_found++;}
      if ($rid1_col){$cols_found++;}
      if ($rid2_col){$cols_found++;}
      unless($cols_found == $cols_expected){
        print RED, "\n\nCould not find neccessary read ID columns\n\n", RESET;
        exit();
      }
      next();
    }
    #Check for invalid read id
    if ($line[0] =~ /(\w+)\_(\d+)\_(\d+)\_(\d+)\_(\d+)/){
      unless($2 == $lane){
        $mr_files{$mr}{invalid_read_ids} = 1;
        $files_to_fix++;
        $| = 0; print BLUE, ".", RESET; $| = 1;
        last();
      }
    }else{
      print RED, "\n\nFormat of read id not understood: $line[0]\n\n", RESET;
      exit();
    }
  }
  close(TEMP);
}
foreach my $br (sort {$a <=> $b} keys %br_files){
  my $path = $br_files{$br}{path};
  my $lane = $br_files{$br}{lane};
  my $name = $br_files{$br}{name};
  open (TEMP, "zcat $path |") || die "\n\nCould not open file: $path\n\n";
  while(<TEMP>){
    chomp($_);
    my @line = split("\t", $_);
    #Check for invalid read id
    if ($line[0] =~ /(\w+)\_(\d+)\_(\d+)\_(\d+)\_(\d+)/){
      unless($2 == $lane){
        $br_files{$br}{invalid_read_ids} = 1;
        $files_to_fix++;
        $| = 0; print BLUE, ".", RESET; $| = 1;
        last();
      }
    }else{
      print RED, "\n\nFormat of read id not understood: $line[0]\n\n", RESET;
      exit();
    }
  }
  close(TEMP);
}
print BLUE, "\n\nFound $files_to_fix files with invalid read IDs that need to be fixed", RESET;


#6.) Now actually fix files that need to be fixed (i.e. if invalid_read_ids==1)
#    - Create a temp file and once complete, compress it, then overwrite the garbage file


#READ RECORD FILES
foreach my $rr (sort {$a <=> $b} keys %rr_files){
  my $path = $rr_files{$rr}{path};
  my $lane = $rr_files{$rr}{lane};
  my $name = $rr_files{$rr}{name};
  my $base_name = $rr_files{$rr}{base_name};
  my $dir = $rr_files{$rr}{dir};
  my $tmp_file = "$dir"."$base_name".".tmp";
  my $invalid_read_ids = $rr_files{$rr}{invalid_read_ids};
  unless($invalid_read_ids == 1){
    next();
  }

  print MAGENTA, "\n\nRepairing: $path", RESET;
  open (IN, "zcat $path |") || die "\n\nCould not open in file: $path\n\n";
  open (OUT, ">$tmp_file") || die "\n\nCould not open out file: $tmp_file\n\n";
  my $header = 1;
  while(<IN>){
    chomp($_);
    my $line = $_;
    my @line = split("\t", $_);
    if ($header == 1){
      $header = 0;
      print OUT "$line\n";
      next();
    }
    #Check for invalid read id
    if ($line[0] =~ /(\w+)\_(\d+)\_(\d+)\_(\d+)\_(\d+)/){
      my $flowcell = $1;
      my $tile = $2;
      my $x = $3;
      my $y = $4;
      my $fixed_rid = "$flowcell"."_"."$lane"."_"."$tile"."_"."$x"."_"."$y";
      my $fixed_r1_id = "$fixed_rid"."_R1";
      my $fixed_r2_id = "$fixed_rid"."_R2";
      $line[0] = $fixed_rid;
      $line[1] = $fixed_r1_id;
      $line[2] = $fixed_r2_id;
      #print YELLOW, "$line\n", RESET;
      #my $tmp = $"; $" = "\t"; print GREEN, "@line\n", RESET; $" = $tmp;

      my $tmp = $"; $" = "\t"; print OUT "@line\n"; $" = $tmp;
    }else{
      print RED, "\n\nFirst column in READ RECORD file not understood\n\n", RESET;
      exit();
    }

  }
  close(IN);
  close(OUT);

  #Compress the output file
  my $cmd_gz = "gzip -f $tmp_file";
  print BLUE, "\n$cmd_gz";
  system($cmd_gz);

  #Overwrite the original file with the new tmp file
  my $cmd_mv = "mv -f $tmp_file".".gz"." $path";
  print BLUE, "\n$cmd_mv", RESET;
  system($cmd_mv);
}


#MAPPING RESULTS FILES
foreach my $mr (sort {$a <=> $b} keys %mr_files){
  my $path = $mr_files{$mr}{path};
  my $lane = $mr_files{$mr}{lane};
  my $name = $mr_files{$mr}{name};
  my $base_name = $mr_files{$mr}{base_name};
  my $dir = $mr_files{$mr}{dir};
  my $tmp_file = "$dir"."$base_name".".tmp";
  my $invalid_read_ids = $mr_files{$mr}{invalid_read_ids};
  unless($invalid_read_ids == 1){
    next();
  }
  
  my %columns = %{$mr_files{$mr}{columns}};
  my $cols_expected = $mr_files{$mr}{columns_expected};

  print MAGENTA, "\n\nRepairing: $path", RESET;
  open (IN, "zcat $path |") || die "\n\nCould not open in file: $path\n\n";
  open (OUT, ">$tmp_file") || die "\n\nCould not open out file: $tmp_file\n\n";
  my $header = 1;
  while(<IN>){
    chomp($_);
    my $line = $_;
    my @line = split("\t", $_);
    if ($header == 1){
      $header = 0;
      print OUT "$line\n";
      next();
    }
    #Check for invalid read id
    if ($cols_expected == 3){
      if ($line[0] =~ /(\w+)\_(\d+)\_(\d+)\_(\d+)\_(\d+)/){
        my $flowcell = $1;
        my $tile = $2;
        my $x = $3;
        my $y = $4;
        my $fixed_rid = "$flowcell"."_"."$lane"."_"."$tile"."_"."$x"."_"."$y";
        my $fixed_r1_id = "$fixed_rid"."_R1";
        my $fixed_r2_id = "$fixed_rid"."_R2";
        my $rid_col = $columns{'Read_ID'}{col_pos};
        my $rid1_col = $columns{'R1_ID'}{col_pos};
        my $rid2_col = $columns{'R2_ID'}{col_pos};
        $line[0] = $fixed_rid;
        $line[$rid1_col-1] = $fixed_r1_id;
        $line[$rid2_col-1] = $fixed_r2_id;
      }else{
        print RED, "\n\nFirst column in MAPPING results file not understood\n\n", RESET;
        exit();
      }
    }else{
      #42HRWAAXX_100_1000_1136_0_R1
      if ($line[0] =~ /(\w+)\_(\d+)\_(\d+)\_(\d+)\_(\d+)\_(\w+)/){
        my $flowcell = $1;
        my $tile = $2;
        my $x = $3;
        my $y = $4;
        my $r = $6;
        my $fixed_rid = "$flowcell"."_"."$lane"."_"."$tile"."_"."$x"."_"."$y"."_"."$r";
        $line[0] = $fixed_rid;
      }else{
        print RED, "\n\nFirst column in MAPPING results file not understood\n\n", RESET;
        exit();
      }
    }
    #print YELLOW, "$line\n", RESET;
    #my $tmp = $"; $" = "\t"; print GREEN, "@line\n", RESET; $" = $tmp;

    my $tmp = $"; $" = "\t"; print OUT "@line\n"; $" = $tmp;
   
  }
  close(IN);
  close(OUT);

  #Compress the output file
  my $cmd_gz = "gzip -f $tmp_file";
  print BLUE, "\n$cmd_gz";
  system($cmd_gz);

  #Overwrite the original file with the new tmp file
  my $cmd_mv = "mv -f $tmp_file".".gz"." $path";
  print BLUE, "\n$cmd_mv", RESET;
  system($cmd_mv);

}


#BLAST RESULTS FILES
foreach my $br (sort {$a <=> $b} keys %br_files){
  my $path = $br_files{$br}{path};
  my $lane = $br_files{$br}{lane};
  my $name = $br_files{$br}{name};
  my $base_name = $br_files{$br}{base_name};
  my $dir = $br_files{$br}{dir};
  my $tmp_file = "$dir"."$base_name".".tmp";
  my $invalid_read_ids = $br_files{$br}{invalid_read_ids};
  unless($invalid_read_ids == 1){
    next();
  }

  print MAGENTA, "\n\nRepairing: $path", RESET;
  open (IN, "zcat $path |") || die "\n\nCould not open in file: $path\n\n";
  open (OUT, ">$tmp_file") || die "\n\nCould not open out file: $tmp_file\n\n";
  my $header = 1;
  while(<IN>){
    chomp($_);
    my $line = $_;
    my @line = split("\t", $_);
    if ($header == 1){
      $header = 0;
      print OUT "$line\n";
      next();
    }
    #Check for invalid read id
    if ($line[0] =~ /(\w+)\_(\d+)\_(\d+)\_(\d+)\_(\d+)\_(\w+)/){
      my $flowcell = $1;
      my $tile = $2;
      my $x = $3;
      my $y = $4;
      my $r = $6;
      my $fixed_rid = "$flowcell"."_"."$lane"."_"."$tile"."_"."$x"."_"."$y"."_"."$r";
      $line[0] = $fixed_rid;
    }else{
      print RED, "\n\nFirst column in BLAST results file not understood\n\n", RESET;
      exit();
    }

    #print YELLOW, "$line\n", RESET;
    #my $tmp = $"; $" = "\t"; print GREEN, "@line\n", RESET; $" = $tmp;

    my $tmp = $"; $" = "\t"; print OUT "@line\n"; $" = $tmp;

  }
  close(IN);
  close(OUT);

  #Compress the output file
  my $cmd_gz = "gzip -f $tmp_file";
  print BLUE, "\n$cmd_gz";
  system($cmd_gz);

  #Overwrite the original file with the new tmp file
  my $cmd_mv = "mv -f $tmp_file".".gz"." $path";
  print BLUE, "\n$cmd_mv", RESET;
  system($cmd_mv);
}


print "\n\n";

exit();

