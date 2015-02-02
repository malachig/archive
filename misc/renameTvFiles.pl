#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

#my $alexa_seq_config_file = '';
#my $project_config_file = '';
#my $commands_file = '';
#my $cluster_commands = '';
#my $clean = '';
#my $demo = '';

#GetOptions ('alexa_seq_config_file=s'=>\$alexa_seq_config_file, 'project_config_file=s'=>\$project_config_file, 'commands_file=s'=>\$commands_file, 'clean=s'=>\$clean, 'cluster_commands=s'=>\$cluster_commands, 'demo=s'=>\$demo);

#print "\n\nHello world\n\n";

#ls the current dir and capture the outcome
my @ls = `ls`;
chomp(@ls);
#print Dumper @ls;

#Get the current dir
my $current_dir = `pwd`;
chomp($current_dir);

#How many files found?
my $file_count = scalar (@ls);
print "#Found $file_count files in $current_dir\n\n";

#Go through each file found
foreach my $file (@ls){
  
  #Quote the file to deal with idiotic characters (space, ', ", etc)
  my $old_fn = $file;
  $old_fn =~ s/\s/\\ /g;
  $old_fn =~ s/'/\\'/g;
  $old_fn =~ s/\(/\\\(/g;
  $old_fn =~ s/\)/\\\)/g;
  $old_fn =~ s/\&/\\&/g;
 
  #Skip directories
  if (-d $file){
    print "#File: is a directory!\n";
    next();
  }

  #If the file is already in the correct format, move on
  if ($file =~ /^s\d{2}e\d{2}\.\w+$/){
    print "#File: $file looks correct\n";
    next();
  }

  #Capture the extension
  my $ext;
  if ($file =~ /(\.\w+)$/){
    $ext = $1;
  }else{
    print RED, "\n\nCould not determine file extension for file:\n$file\n\n", RESET;
    exit();
  }

  #Capture the season and episode number
  my $season;
  my $episode;

  if ($file =~ /(\d{2})\w(\d{2})/){
    #Simpsons 07x01 - Who Shot Mr Burns (Part 2) [rl].avi
    $season = $1; $episode = $2;
  }elsif($file =~ /(\d{1})(\d{2})/){
    #family.guy.909.hdtv-lol.avi
    $season = "0$1"; $episode = $2;
     
  }else{
    print RED, "\n\nCould not determine season/episode numbers for: $file\n\n", RESET;
    exit();
  }


  #Print the mv command
  print "mv $old_fn s$season"."e$episode"."$ext\n";


}


print "\n\n";
exit();




