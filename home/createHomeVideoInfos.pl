#!/usr/bin/perl -w

#create .nfo files for home videos of the form:

#<?xml version=1.0 encoding=UTF-8 standalone=yes ?>
#<musicvideo>
#	<title>2008.01.01 09:23 2008</title>
#	<album>2008 (2008)</album>
#	<year>2008</year>
#</musicvideo>

#If the video is "filename.mp4" then the info file should be "filename.nfo"

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

my $dir = '';

GetOptions ('dir=s'=>\$dir);

#Make sure parameters were provided
unless($dir1){
  print RED, "\n\nRequired parameters missing\n\n", RESET;
  print GREEN, "Usage:  createHomeVideoInfos.pl --dir=/Volumes/Volume_1/DLINK_NFS1/Home Videos/", RESET;
  exit();
}

#Make sure these are actually directories that exist
unless (-e $dir && -d $dir){
  print RED, "\n\nDir ($dir) could not be found or is not a directory\n\n", RESET;
  exit();
}

#ls the current dir and capture the outcome
my @ls = `ls`;
chomp(@ls);
#print Dumper @ls;

#Go through each file found
foreach my $file (@ls){
  
  #Capture the extension
  my $ext;
  if ($file =~ /(\.\w+)$/){
    $ext = $1;
  }else{
    print RED, "\n\nCould not determine file extension for file:\n$file\n\n", RESET;
    exit();
  }



}


print "\n\n";
exit();

