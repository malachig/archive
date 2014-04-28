#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Examine a directory containing summarized quality statistics already parsed with readQualityStats.pl
#These files contain quality values summarized according to their position as well as their score value
#i.e. distribution of scores according to position and distribution of the frequency of each score
#Seperate files were created for both R1 and R2 of all lanes of data.  
#Grand summary files are also provided for all R1 and R2 as well as R1/R2 combined

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Load the ALEXA libraries
my $script_dir;
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
    $script_dir = $1;
  }
}
use utilities::utility qw(:all);

#Initialize command line options
my $data_dir = '';
my $working_dir = '';

GetOptions('data_dir=s'=>\$data_dir, 'working_dir=s'=>\$working_dir);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the complete path to a data directory containing summarized quality stats using: --data_dir", RESET;
print GREEN, "\n\tSpecify the complete path to a results directory using: --working_dir", RESET;

print GREEN, "\n\nUsage: summarizeQualityScores.pl  --data_dir=~/quality_scores/  --working_dir=~/results/\n", RESET;

unless ($data_dir && $working_dir){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}  

$data_dir = &checkDir('-dir'=>$data_dir, '-clear'=>"no");
$working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"yes");

my %data_list;

#get all data files
my $max_qual = 0;
my $max_pos = 0;

my $files_ref = &getDataFiles('-input_dir'=>$data_dir);

#Summarize grand R1 file
print BLUE, "\nProcessing: $files_ref->{R1_all}", RESET;
my $data_ref = &parseDataFile('-infile'=>$files_ref->{R1_all});
$data_list{R1_all} = $data_ref;
&processData('-data_list'=>\%data_list, '-working_dir'=>$working_dir, '-data_name'=>'R1_mean');

#Summarize grand R2 file
%data_list = ();
print BLUE, "\n\nProcessing: $files_ref->{R2_all}", RESET;
$data_ref = &parseDataFile('-infile'=>$files_ref->{R2_all});
$data_list{R2_all} = $data_ref;
&processData('-data_list'=>\%data_list, '-working_dir'=>$working_dir, '-data_name'=>'R2_mean');

#Summarize grand R1/R2 file
%data_list = ();
print BLUE, "\n\nProcessing: $files_ref->{all}", RESET;
$data_ref = &parseDataFile('-infile'=>$files_ref->{all});
$data_list{all} = $data_ref;
&processData('-data_list'=>\%data_list, '-working_dir'=>$working_dir, '-data_name'=>'R12_mean');

#Summarize individual R1 files
my %data_list_r1;
my %data_list_r2;
my %data_list_all;

my $r1_files = $files_ref->{R1};
my $r2_files = $files_ref->{R2};

foreach my $file_num (sort {$r1_files->{$a}->{name} cmp $r1_files->{$b}->{name}} keys %{$r1_files}){

  my $file = $r1_files->{$file_num}->{path};
  my $name = $r1_files->{$file_num}->{name};
  
  my $data_ref = &parseDataFile('-infile'=>$file);
  $data_list_r1{$name} = $data_ref;
  $data_list_all{$name} = $data_ref;
  
}

foreach my $file_num (sort {$r2_files->{$a}->{name} cmp $r2_files->{$b}->{name}} keys %{$r2_files}){
  my $file = $r2_files->{$file_num}->{path};
  my $name = $r2_files->{$file_num}->{name};

  my $data_ref = &parseDataFile('-infile'=>$file);
  $data_list_r2{$name} = $data_ref;
  $data_list_all{$name} = $data_ref;

}


print BLUE, "\n\nProcessing: R1 files", RESET;
&processData('-data_list'=>\%data_list_r1, '-working_dir'=>$working_dir, '-data_name'=>'R1_lanes');

print BLUE, "\n\nProcessing: R2 files", RESET;
&processData('-data_list'=>\%data_list_r2, '-working_dir'=>$working_dir, '-data_name'=>'R2_lanes');

print BLUE, "\n\nProcessing: ALL files", RESET;
&processData('-data_list'=>\%data_list_all, '-working_dir'=>$working_dir, '-data_name'=>'R12_lanes');



print "\n\n";

exit();


###########################################################################################################
#get data files
###########################################################################################################
sub getDataFiles{
  my %args = @_;
  my $dir = $args{'-input_dir'};


  my $dh = opendir(DIR, $dir) || die "\nCould not open directory: $dir\n\n";

  my @files = readdir(DIR);
  my %files;
  my %files_R1;
  my %files_R2;
  
  my $count = 0;
  foreach my $file (@files){
    chomp($file);
    unless ($file =~ /\.txt$/){
      next();
    }
    if (-d $file){
      next();
    }
    $count++;
    #print "\n$file";

    if ($file =~ /(\w+)_(\w+)_(\w+)\.(.*)/){
      #print YELLOW, "\n$1\t$2\t$3\t$4", RESET;
      if ($3 eq "R1"){
        $files_R1{$count}{file_name} = $file;
        $files_R1{$count}{path} = "$dir"."$file";
        $files_R1{$count}{flowcell} = $1;
        $files_R1{$count}{lane} = $2;
        $files_R1{$count}{read_num} = $3;
        $files_R1{$count}{name} = "$1"."_"."$2"."_"."$3";
        $files_R1{$count}{extension} = $4;
      }elsif($3 eq "R2"){
        $files_R2{$count}{file_name} = $file;
        $files_R2{$count}{path} = "$dir"."$file";
        $files_R2{$count}{flowcell} = $1;
        $files_R2{$count}{lane} = $2;
        $files_R2{$count}{read_num} = $3;
        $files_R2{$count}{name} = "$1"."_"."$2"."_"."$3";
        $files_R2{$count}{extension} = $4;
      }
    }
  }

  $files{R1} = \%files_R1;
  $files{R2} = \%files_R2;
  $files{R1_all} = "$dir"."ALL_R1.qual_stats.txt";
  $files{R2_all} = "$dir"."ALL_R2.qual_stats.txt";
  
  $files{all} = "$dir"."ALL.qual_stats.txt";
  closedir(DIR);
  
  return(\%files);
}

##################################################################################################################
#Parse data file                                                                                                 #
##################################################################################################################
sub parseDataFile{
  my %args = @_;
  my $file = $args{'-infile'};
  my %data;
  my %pos;
  my %qual;
 
  open (DATA, "$file") || die "\nCould not open input file: $file";

  while(<DATA>){
    chomp($_);

    if ($_ =~ /^#/){
      next();
    }

    if ($_ =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
      #Get the quality distribution values (5 columns)
      #Quality from 1 to 40
      my $quality = $1;
      if ($quality > $max_qual){$max_qual = $quality;}
      
      $qual{$quality}{read_count} = $2; #Number of reads with this quality score
      $qual{$quality}{percent} = $3;    #Percentage of all reads with this quality score
      $qual{$quality}{cumulative_percent_forward} = $4;
      $qual{$quality}{cumulative_percent_reverse} = $5;
      
    }elsif($_ =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
      #Get the position distribution values (4 columns)
      my $position = $1;
      if ($position > $max_pos){$max_pos = $position;}
      
      $pos{$position}{average_quality} = $4;
      
    }
  }
  $data{qual} = \%qual;
  $data{pos} = \%pos;
  
  close(DATA);
  
  return(\%data);
}


#####################################################################################################
#Print data to temp file
#####################################################################################################
sub processData{
  my %args = @_;
  my $data_list = $args{'-data_list'};
  my $working_dir = $args{'-working_dir'};
  my $data_name = $args{'-data_name'};

  my @quals = (1 .. $max_qual);
  my @pos = (1 .. $max_pos);

  my $temp_file_qual_by_dataset = "$working_dir"."$data_name".".qual_by_data.txt";
  my $temp_file_pos_by_dataset = "$working_dir"."$data_name".".pos_by_data.txt";
  my $temp_file_dataset_by_qual = "$working_dir"."$data_name".".data_by_qual.txt";
  my $temp_file_dataset_by_pos = "$working_dir"."$data_name".".data_by_pos.txt";
    
  open (QUAL_BY_DATA, ">$temp_file_qual_by_dataset") || die "\nCould not open temp file: $temp_file_qual_by_dataset\n\n";
  open (POS_BY_DATA, ">$temp_file_pos_by_dataset") || die "\nCould not open temp file: $temp_file_pos_by_dataset\n\n";
  open (DATA_BY_QUAL, ">$temp_file_dataset_by_qual") || die "\nCould not open temp file: $temp_file_dataset_by_qual\n\n";
  open (DATA_BY_POS, ">$temp_file_dataset_by_pos") || die "\nCould not open temp file: $temp_file_dataset_by_pos\n\n";
  
  foreach my $data_set (sort keys %{$data_list}){
    
    print BLUE, "\n\tProcessing data: $data_set", RESET;
    
    my $data_ref = $data_list->{$data_set};

    my $qual_ref = $data_ref->{qual};
    my $pos_ref = $data_ref->{pos};
  
    foreach my $qual (@quals){
      unless($qual_ref->{$qual}){
        $qual_ref->{$qual}->{percent} = "NA";
      }
    }
    
    foreach my $pos (@pos){
      unless($pos_ref->{$pos}){
        $pos_ref->{$pos}->{average_quality} = "NA";
      }
    }
  }

  #Now print columns of data (one for each lane of data) - each row corresponds to a position or quality score 
  my $i = 0;
  foreach my $data_set (sort keys %{$data_list}){
    if ($i == 0){$i = 1; print QUAL_BY_DATA "Quality\t"; print POS_BY_DATA "Position\t";}
    print QUAL_BY_DATA "$data_set\t";
    print POS_BY_DATA "$data_set\t";
  }
  print QUAL_BY_DATA "\n";
  print POS_BY_DATA "\n";

  foreach my $qual (@quals){
    $i = 0; 
    foreach my $data_set (sort keys %{$data_list}){
      my $data_ref = $data_list->{$data_set};
      my $qual_ref = $data_ref->{qual};

      if ($i == 0){$i = 1; print QUAL_BY_DATA "$qual\t";}
      print QUAL_BY_DATA "$qual_ref->{$qual}->{percent}\t"; 
    }
    print QUAL_BY_DATA "\n"; 
  }

  foreach my $pos (@pos){
    $i = 0;

    foreach my $data_set (sort keys %{$data_list}){
      my $data_ref = $data_list->{$data_set};
      my $pos_ref = $data_ref->{pos};
     
      if ($i == 0){$i = 1; print POS_BY_DATA "$pos\t";}
      print POS_BY_DATA "$pos_ref->{$pos}->{average_quality}\t"; 
    }
    print POS_BY_DATA "\n"; 
  }


  #Now print columns of data (one for each position or quality score) - each row corresponds to a lane of data
  $i = 0;
  foreach my $qual (@quals){
    if ($i == 0){$i = 1; print DATA_BY_QUAL "Data_set\t";}
    print DATA_BY_QUAL "$qual\t";
  }
  print DATA_BY_QUAL "\n";

  $i = 0;
  foreach my $pos (@pos){
    if ($i == 0){$i = 1; print DATA_BY_POS "Data_set\t";}
    print DATA_BY_POS "$pos\t";
  }
  print DATA_BY_POS "\n";


  foreach my $data_set (sort keys %{$data_list}){
    $i = 0;
    my $data_ref = $data_list->{$data_set};
    my $qual_ref = $data_ref->{qual};

    foreach my $qual (@quals){
      if ($i == 0){$i = 1; print DATA_BY_QUAL "$data_set\t";}
      print DATA_BY_QUAL "$qual_ref->{$qual}->{percent}\t"; 
    }
    print DATA_BY_QUAL "\n"; 
  }


  foreach my $data_set (sort keys %{$data_list}){
    $i = 0;
    my $data_ref = $data_list->{$data_set};
    my $pos_ref = $data_ref->{pos};

    foreach my $pos (@pos){
      if ($i == 0){$i = 1; print DATA_BY_POS "$data_set\t";}
      print DATA_BY_POS "$pos_ref->{$pos}->{average_quality}\t"; 
    }
    print DATA_BY_POS "\n"; 
  }

  close(QUAL_BY_DATA);
  close(POS_BY_DATA);
  close(DATA_BY_QUAL);
  close(DATA_BY_POS);
  return();
}




