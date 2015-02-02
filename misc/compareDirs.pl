#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

my $dir1 = '';
my $dir2 = '';

GetOptions ('dir1=s'=>\$dir1, 'dir2=s'=>\$dir2);

#Make sure parameters were provided
unless($dir1 && $dir2){
  print RED, "\n\nRequired parameters missing\n\n", RESET;
  print GREEN, "Usage:  compareDirs.pl --dir1=/Volumes/Volume_1/DLINK_NFS1/Movies/  --dir2=/Volumes/Volume_1-1/DLINK_NFS2/Movies/\n\n", RESET;
  exit();
}

#Make sure these are actually directories that exist
unless (-e $dir1 && -d $dir1){
  print RED, "\n\nDir1 ($dir1) could not be found or is not a directory\n\n", RESET;
  exit();
}
unless (-e $dir2 && -d $dir2){
  print RED, "\n\nDir2 ($dir2) could not be found or is not a directory\n\n", RESET;
  exit();
}

#Make sure both dir paths have a trailing '/'
unless ($dir1 =~ /\/$/){
  $dir1 .= "/";
}
unless ($dir2 =~ /\/$/){
  $dir2 .= "/";
}

#my $test1 = "/Volumes/Volume_1/DLINK_NFS1/TV/30 Rock/Season 1";
#$test1 =~ s/ /\\ /g;
#print "\n\ntest1: $test1\n\n";
#if (-d $test1){
#  print "\n\nDir Found\n\n";
#}else{
#  print "\n\nDir Not found\n\n";
#}
#opendir(DIR, $test1);
#my @result = readdir(DIR);
#my @result = `ls $test1`;
#print YELLOW, "\n\nls result: @result\n\n", RESET;
#exit();

#Escape spaces in input paths
#$dir1 =~ s/ /\\\\ /g;
#$dir2 =~ s/ /\\\\ /g;

print BLUE, "\n\nLooking for files in Dir1 ('old dir') and then seeing if all these files exist in Dir2 ('new_dir')", RESET;


#ls the dir1 and capture the outcome
print BLUE, "\n\nFirst gathering files from dir1: $dir1", RESET;
my $dir1_temp = $dir1;
$dir1_temp =~ s/\\//g;
opendir(DIR, "$dir1_temp") || die "\n\nCould not open dir for reading: $dir1_temp";
my @ls = readdir(DIR);
chomp(@ls);


#How many files found?
my $file_count = scalar (@ls);
print BLUE, "\n\tFound $file_count files/directories in top level";

#Go through each file found/directory found and add the prefix path - but store the original base path as well
my %files;
my $fc = 0;
foreach my $file (@ls){

  #Skip files/dir starting with "."
  if ($file =~ /^\./){
    next();
  }

  $fc++;
  #$file =~ s/ /\\\\ /g;
  my $path = "$dir1"."$file";
  $files{$fc}{base} = $dir1;       #The original base dir
  $files{$fc}{full_path} = $path;  #The full path
  $files{$fc}{path} = $file;       #The path with everything after the base path
  $files{$fc}{size} = -s $path;
}

my $dirs_found = 1;
while($dirs_found = &expandPath('-files'=>\%files)){}

my $total_files_found = keys %files;
print BLUE, "\n\tFound $total_files_found files in all levels";

my $found_files = 0;
my $missing_files = 0;
foreach my $fc (sort {$files{$a}{full_path} cmp $files{$b}{full_path}} keys %files){
  my $full_path1 = $files{$fc}{full_path};
  my $path1 = $files{$fc}{path};
  my $full_path2 = "$dir2"."$path1";
  my $size1 = $files{$fc}{size};

  if ((-e $full_path2) && (-s $full_path2 == $size1)){
    $found_files++;
    #print GREEN, "\n\tFile found: $full_path1 -> $full_path2", RESET;
  }else{
    $missing_files++;
    print YELLOW, "\n\tFile missing: $full_path1 -> $full_path2", RESET;
  }
}

print BLUE, "\n\nFound $found_files in both locations, and $missing_files are missing from Dir2", RESET;

#print Dumper %files;

print "\n\n";
exit();



############################################################################################################
#A sub routine that takes a list of file paths and expands those that are actually directories             #
#The number of expanded dirs is returned.  If all paths resolve to files already, 0 is returned            #
############################################################################################################
sub expandPath{
  my %args = @_;
  my $files_ref = $args{'-files'};

  my $dir_count = 0;

  print YELLOW, "\n\t\tExpanding dirs where appropriate... ", RESET;
 
  foreach my $file (sort {$a <=> $b} keys %{$files_ref}){
    my $old_base = $files_ref->{$file}->{base};
    my $old_full_path = $files_ref->{$file}->{full_path};
    my $old_path = $files_ref->{$file}->{path};

    if (-d $old_full_path){
      $dir_count++;
      unless ($old_full_path =~ /\/$/){
        $old_full_path .= "/";
      }
      unless ($old_path =~ /\/$/){
        $old_path .= "/";
      }

      my $old_full_path_temp = $old_full_path;

      opendir(DIR, $old_full_path_temp) || die "\n\nCould not open dir for reading: $old_full_path_temp\n\n";
      my @ls = readdir(DIR);
      chomp(@ls);

      foreach my $tf (@ls){
        #Skip files/dir starting with "."
        if ($tf =~ /^\./){
          next();
        }

        $fc++;
        #$tf =~ s/ /\\ /g;
        my $new_full_path = "$old_full_path"."$tf";
        my $new_path = "$old_path"."$tf";
        $files{$fc}{base} = $old_base;            #The original base dir
        $files{$fc}{full_path} = $new_full_path;  #The full path
        $files{$fc}{path} = $new_path;                #The path with everything after the base path
        $files{$fc}{size} = -s $new_full_path;
      }

      #Remove the dir from the file hash
      delete $files_ref->{$file};
    }
  }

  print YELLOW, "found $dir_count dirs", RESET;
  return ($dir_count);
}




