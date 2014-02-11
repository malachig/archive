#!/usr/bin/perl -w
#Written by Malachi Griffith

#Throw away script to rename misc. clones and add to existing assemblies based on clone names


use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

my $assembly_dir = "/home/malachig/AlternativeSplicing/UMPS_Clone_Assemblies/";
my $input_dir = "/home/malachig/AlternativeSplicing/UMPS_Clone_Assemblies/temp_reads/";

opendir(DIRHANDLE, "$input_dir") || die "\nCannot open directory: $input_dir\n\n";
my @files = readdir(DIRHANDLE);

my %trace_files;

my $count = 0;
foreach my $file (@files){
  unless ($file =~ /\.ab1/){
    print YELLOW, "\nFile: $file is not an ab1 file - skipping", RESET;
    next();
  }
  $count++;

  $trace_files{$count}{file_name} = "$file";
  $trace_files{$count}{file_path} = "$input_dir"."$file";

  my $plate;
  my $read_direction;
  my $clone_name;
  my $well;
  if ($file =~ /(MG\d+)_(\w)_(\w+_\d+_\d+)_(\w+)\.ab1/){
    $plate = $1;
    $read_direction = $2;
    $clone_name = $3;
    $well = $4;
  }else{
    print RED, "\nCould not determine info from clone name\n\n", RESET;
    exit();
  }
  $trace_files{$count}{plate} = $plate;
  $trace_files{$count}{read_direction} = $read_direction;
  $trace_files{$count}{clone_name} = $clone_name;
  $trace_files{$count}{well} = $well;
  $trace_files{$count}{assembly_path} = "$assembly_dir"."$clone_name"."/";
  $trace_files{$count}{chromat_dir_path} = "$assembly_dir"."$clone_name"."/chromat_dir/";

  #create new trace name
  my $new_name;
  if ($read_direction eq "F"){
    $new_name = "$plate".".C21_"."$well".".ab1";

  }elsif ($read_direction eq "R"){
    $new_name = "$plate".".CR_"."$well".".ab1";

  }else{
    print RED, "\nRead direction does not appear valid\n\n", RESET;
    exit();
  }

  #make sure the target chromat dir path exists!
  unless (-e $trace_files{$count}{chromat_dir_path} && -d $trace_files{$count}{chromat_dir_path}){
    print RED, "\nTarget chromat dir does not exist: $trace_files{$count}{chromat_dir_path}\n\n", RESET;
    exit();
  }

  $trace_files{$count}{new_name} = $new_name;
  $trace_files{$count}{new_path} = "$trace_files{$count}{chromat_dir_path}"."$new_name";

}

print Dumper %trace_files;

#Copy files to new location
foreach my $file (sort {$a <=> $b} keys %trace_files){

  my $cp_file_cmd = "cp $trace_files{$file}{file_path} $trace_files{$file}{new_path}";

  print BLUE, "\n$cp_file_cmd", RESET;
  system("$cp_file_cmd");

}


exit();
