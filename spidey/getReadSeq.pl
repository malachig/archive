#!/usr/bin/perl

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my @seqs;
my $seq_count;
my %read_id;
my $analysis_dir = '';
my $library_id = '';
my $test = '';

#User specified option for:
#$analysis_dir (e.g.  /projects/alexa2/alexa_seq/)
#$library_id (e.g. HS04391/)
GetOptions( 'analysis_dir=s' => \$analysis_dir,
            'library_id=s' => \$library_id,
            'test=i' => \$test );

unless ($analysis_dir && $library_id)
{
  print MAGENTA, "Please specify the analysis directory, library ID, using the following format:\n";
  print "--analysis_dir=/projects/alexa2/alexa_seq/ --library_id=HS04391 --test=1\n";
  print "getReadSeq.pl terminated\n\n", RESET;
  exit();

}

unless ($analysis_dir =~ /\/$/)
{
  #add a slash (/) at the end of the dir string
  $analysis_dir = "$analysis_dir"."\/";
}


while (my $line = <STDIN>)
{
  chomp $line;
  #print $line; print "\n";
  $read_id{$line}{defined} = 1;

}


my $seq_dir = "$analysis_dir"."read_records\/$library_id\/";

opendir (DIR, "$seq_dir") || die MAGENTA, "\n\nCould not open directory. getReadSeq.pl terminated\n\n", RESET;
my @dir_list = readdir(DIR);
closedir (DIR);


#Filtering out the non-txt.gz files
#Pushing the txt.gz files into the @reads array

foreach my $item (@dir_list)
{
  if ($item =~ /\.txt\.gz/) #using 20836AAXX_Lane5.txt.gz for testing purposes
  {
    push @seqs, $item;
  }
}
my @seqs_sorted = sort {$a cmp $b} @seqs;

$seq_count = scalar(@seqs_sorted);

#test purposes:
if ($test)
{
  $seq_count = $test;
}

my @seq_line;

for (my $m = 0; $m < $seq_count; $m++)
{
  #Opening the file
  open (SEQ, "zcat $seq_dir"."$seqs_sorted[$m] |");

  while (<SEQ>)
  {
    chomp($_);
    #Searching for read IDs that are defined in hash %read_id
    @seq_line = 0;
    (@seq_line) = split("\t");
    
    if ($read_id{$seq_line[0]})
    {
      print "\>$seq_line[1]\n";
      print "$seq_line[5]\n";
      print "\>$seq_line[2]\n";
      print "$seq_line[6]\n";
    }
  }
  close (SEQ);
}




