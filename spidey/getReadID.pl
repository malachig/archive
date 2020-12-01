#!/usr/bin/perl


#
#Command line: /home/bcho/Perl/GeneJunction/getReadID.pl --analysis_dir='/projects/alexa2/alexa_seq/' --library_id='HS04391'
#
#
#

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my @reads;
my $read_count;
my $analysis_dir;
my $library_id;
my $gene_name;
my $ensembl_version;
my $test;

#User specified option for:
#$analysis_dir (e.g.  /projects/alexa2/alexa_seq/)
#$library_id (e.g. HS04391/)
GetOptions( 'analysis_dir=s' => \$analysis_dir,
            'library_id=s' => \$library_id,
            'gene_name=s' => \$gene_name,
            'ensembl_version=s' => \$ensembl_version,
            'test=i' => \$test );

unless ($analysis_dir && $library_id && $gene_name && $ensembl_version)
{
  print MAGENTA, "Please specify the analysis directory, library ID, using the following format:\n";
  print "--analysis_dir=/projects/alexa2/alexa_seq/ --library_id=HS04391 --gene_name=CCNK --ensembl_version=53 --test=1\n";
  print "getReadID.pl terminated\n\n", RESET;
  exit();
}

unless ($analysis_dir =~ /\/$/)
{
  #add a slash (/) at the end of the dir string
  $analysis_dir = "$analysis_dir"."\/";
}

my $ENST_dir = "$analysis_dir"."read_records\/$library_id\/ENST_v$ensembl_version\/";

opendir (DIR, "$ENST_dir") || die MAGENTA, "\n\nCould not open directory. getReadID.pl terminated\n\n", RESET;
my @dir_list = readdir(DIR);
closedir (DIR);

#Filtering out the non-txt.gz files
#Pushing the txt.gz files into the @reads array

foreach my $item (@dir_list)
{
  if ($item =~ /\.txt\.gz/) #using 20836AAXX_Lane5_ENST_v53 for testing purposes
  {
    push @reads, $item;
  }

}

my @reads_sorted = sort {$a cmp $b} @reads;

$read_count = scalar(@reads_sorted);

#test purposes
if ($test)
{
  $read_count = $test;
}

my @ENST_line;

for (my $m = 0; $m < $read_count; $m++)
{
  #Opening the file
  open (MAP, "zcat $ENST_dir"."$reads_sorted[$m] |");

  while (<MAP>)
  {
    chomp($_);
    #Searching for user-specified gene name
    if ($_ =~ /\s$gene_name\s/)
    {
      @ENST_line = 0;
      (@ENST_line) = split("\t");
      print $ENST_line[0]; print "\n"
    }
  }

  close (MAP);

}


