#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use File::Copy qw(copy);

my $dirty_dir = '';
my $clean_dir = '';

GetOptions ('dirty_dir=s'=>\$dirty_dir, 'clean_dir=s'=>\$clean_dir);

#Make sure parameters were provided
unless($clean_dir && $dirty_dir){
  print RED, "\n\nRequired parameters missing\n\n", RESET;
  print GREEN, "Usage:  removeDuplicateMusicFiles.pl --clean_dir=/Users/mgriffit/Music/Import/Music/clean  --dirty_dir=/Users/mgriffit/Music/Import/Music/dirty  \n\n", RESET;
  exit();
}

#Make sure these are actually directories that exist
unless (-e $clean_dir && -d $clean_dir){
  print RED, "\n\nClean dir ($clean_dir) could not be found or is not a directory\n\n", RESET;
  exit();
}
unless (-e $dirty_dir && -d $dirty_dir){
  print RED, "\n\nDirty dir ($dirty_dir) could not be found or is not a directory\n\n", RESET;
  exit();
}

#Make sure both dir paths have a trailing '/'
unless ($clean_dir =~ /\/$/){
  $clean_dir .= "/";
}
unless ($dirty_dir =~ /\/$/){
  $dirty_dir .= "/";
}

#Escape spaces in input paths
#$dir1 =~ s/ /\\\\ /g;
#$dir2 =~ s/ /\\\\ /g;

#Get the artist dirs
print BLUE, "\n\nFirst gathering artist directories from clean_dir: $clean_dir", RESET;
opendir(DIR, "$dirty_dir") || die "\n\nCould not open dir for reading: $dirty_dir\n\n";
my @artists = readdir(DIR);
chomp(@artists);

#How many artist directories were found?
my $artist_count = scalar (@artists);
print BLUE, "\n\tFound $artist_count artist directories in top level", RESET;


foreach my $artist_dir (sort @artists){
  next if ($artist_dir =~ /^\./);

  print BLUE, "$artist_dir\n", RESET;
  my $full_artist_dir_path = $dirty_dir . $artist_dir;
  my $full_artist_dir_path_clean = $clean_dir . $artist_dir;
  mkdir $full_artist_dir_path_clean || die "\n\nCould not create new clean dir: $full_artist_dir_path_clean\n\n";

  opendir(DIR, "$full_artist_dir_path") || die "\n\nCould not open dir for reading: $full_artist_dir_path\n\n";
  my @albums = readdir(DIR);
  chomp(@albums);

  foreach my $album_dir (sort @albums){
    next if ($album_dir =~ /^\./);

    print YELLOW, "\t$album_dir\n", RESET;
    my $full_album_dir_path = $full_artist_dir_path . "/" . $album_dir;
    my $full_album_dir_path_clean = $full_artist_dir_path_clean . "/" . $album_dir;
    mkdir $full_album_dir_path_clean || die "\n\nCould not create new clean dir: $full_album_dir_path_clean\n\n";;

    opendir(DIR, "$full_album_dir_path") || die "\n\nCould not open dir for reading: $full_album_dir_path\n\n";
    my @songs = readdir(DIR);
    chomp(@songs);
    
    #remove files starting with "."
    my @songs_clean;
    foreach my $song_file (@songs){
      next if ($song_file =~ /^\./);
      push(@songs_clean, $song_file);
    }
    my $c_expected = scalar(@songs_clean);
    
    my %files;
    my $c = 0;

    #store files that do NOT end with a number
    foreach my $song_file (sort @songs_clean){
      my $size = -s $full_album_dir_path . "/" . $song_file || die "\n\nCould not get file size\n\n";
      print RED, "\t\t$song_file\n", RESET;
      my $base;
      if ($song_file =~ /(.*)\.\w+$/){$base = $1;}else{die "\n\nCould not get base name for song: $song_file\n\n";}
      unless($song_file =~ /\d+\.\w+$/){
        if ($files{$base}){
          my $list = $files{$base}{list};
          $list->{$song_file}->{size} = $size;
          $list->{$song_file}->{c} = $c;
          $c++;
        }else{
          my %list;
          $list{$song_file}{size} = $size;
          $list{$song_file}{c} = $c;
          $files{$base}{list} = \%list;
          $c++;
        }
      }
    }
      
    #check for files that are identical except for trailing number
    #this time skip files unless they DO end with a number
    foreach my $song_file (sort @songs_clean){
      my $size = -s $full_album_dir_path . "/" . $song_file || die "\n\nCould not get file size\n\n";
      my $base;
      next unless($song_file =~ /\d+\.\w+$/); #song must end with a number
      if ($song_file =~ /(.*)\s+\d+\.\w+$/){ #grab base without trailing number for duplicate files "something 1.mp3"
        $base = $1;
      }elsif($song_file =~ /(.*)\.\w+$/){ #grap base for other names styles that happen to end with a number
        $base = $1;
      }else{
        die "\n\nCould not get base name for song: $song_file\n\n";
      }

      if ($files{$base}){
        my $list = $files{$base}{list};
        $list->{$song_file}->{size} = $size;
        $list->{$song_file}->{c} = $c;
        $c++;
      }else{
        my %list;
        $list{$song_file}{size} = $size;
        $list{$song_file}{c} = $c;
        $files{$base}{list} = \%list;
        $c++;
      }
    }

    #print Dumper %files;

    # make sure all files are accounted for
    unless ($c == $c_expected){die "\n\nExpected $c_expected files but only stored $c\n\n"};

    print "\n";

    # now go through the list of files found and chose a single file (the largest) when multiple files are available
    my @chosen_files;
    foreach my $base (sort keys %files){
      my $list = $files{$base}{list};
      my $chosen_file;

      my $max_size = 0;
      foreach my $song_file (sort {$list->{$b}->{size} <=> $list->{$a}->{size}}keys %{$list}){
        my $size = $list->{$song_file}->{size};
        $max_size = $size if ($size > $max_size);
        #print GREEN, "\t\t$song_file\t$size\n", RESET;
      }

      #choose the one file that is largest, or observed first in the event of ties
      my $found = 0;      
      foreach my $song_file (sort {$list->{$a}->{c} <=> $list->{$b}->{c}}keys %{$list}){
        my $size = $list->{$song_file}->{size};
        next unless $size == $max_size;
        next if ($found > 0);
        #print MAGENTA, "\t\t$song_file\t$size\n", RESET;
        push (@chosen_files, $song_file);;
        $found++;
      }
     
    }

    #copy the chosen files to the new clean dir, maintaining the artist / album / song structure
    foreach my $song_file (sort @chosen_files){
      print MAGENTA, "\t\t$song_file\n", RESET;
      my $full_song_path_dirty = $full_album_dir_path . "/" . $song_file;
      my $full_song_path_clean = $full_album_dir_path_clean . "/" . $song_file;
      #print MAGENTA, "\t\t\t$full_song_path_dirty $full_song_path_clean\n", RESET;
      copy $full_song_path_dirty, $full_song_path_clean;
    }

  }
}




print "\n\n";

exit;

