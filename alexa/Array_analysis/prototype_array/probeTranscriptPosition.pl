#!/usr/bin/perl -w
#Written by Malachi Griffith
#This script takes any tab delimited probe file with ALEXA probe IDs in the first column, finds the position of the probe within the gene
#  - Position is with respect to the exon content and thus gives an estimate of its position in possible transcripts from that gene
#The output is simply the input file with the distances added as an extra column
#Distances are expressed as the distance from the 3' end of the transcript (although could easily be expressed relative to 5' end or as percents)

#To do this will involve the following steps
#Open the input file and for each line containing a probe ID, do the following
#1.) Get the probe ID
#2.) Get the unit1 start position (start of the probe), target gene ID, and probe type for this probe, using the probe ID
#    Note: positions will only be calculated for Exon, and Canonical junction probes (0 exons skipped)
#3.) Using the unit1 start position and geneID determine the relative position of the probe
#4.) Output the record

use strict;
use Data::Dumper;
use Getopt::Long;

#Always use the seqdev folder as the library root and then specify the packages to use.
use lib '/usr/local/ulib/beta/seqdev';
use utilities::utility qw(:all);

use lib '/home/malachig/AlternativeSplicing/Array_design/perl_bin/';
use utilities::ALEXA_DB qw(:all);

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $probe_file = '';
my $out_file = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'probe_file=s'=>\$probe_file, 'out_file=s'=>\$out_file);

#Provide instruction to the user
print "\n\nUsage:";
print "\n\tSpecify the database and server to query using: --database and --server";
print "\n\tSpecify the user and password for access using: --user and --password";
print "\n\tSpecify the input probe file using: --probe_file";
print "\n\tSpecify the output file using: --out_file";
print "\n\nExample: probeTrannscriptPosition.pl --database ALEXA_hs_31_35d --server jango.bcgsc.ca --user malachig --password pwd --probe_file test.txt --out_file test.out\n\n";

#Make sure all options were specified
unless ($database && $server && $user && $password && $probe_file && $out_file){
  print "\nOptions missing!\n\n";
  exit();
}


#1.) Open input file and probes to be processed 
my %probes;
my $header;
&getProbes('-probe_file'=>$probe_file);

#2.) Determine the distance from the 3' end of the transcript for each probe
#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = connectDB($database, $server, $user, $password);

&getDistances('-dbh'=>$alexa_dbh);

#Close database connection
$alexa_dbh->disconnect();

#3.) Write a the new output file, appending the calculated values onto the last column
&writeFile('-out_file'=>$out_file);

exit();


###############################################################################################################
#1.) getProbes                                                                                                #
###############################################################################################################
sub getProbes{
  my %args = @_;
  my $probe_file = $args{'-probe_file'};

  print "\nImporting data from input file\n\n";

  open (PROBES, "$probe_file") || die "\nCould not open input probe file: $probe_file\n\n";

  my $probe_count = 0;
  my $header_line = 1;
  while (<PROBES>){
    chomp($_);

    if ($header_line == 1){
      $header = $_;
      $header_line = 0;
      next();
    }

    my @line = split("\t", $_);

    #Skip lines that do not have a probe ID
    my $probe_id;
    if ($line[0] =~ /^(\d+)$/){
      $probe_id = $1;
    }else{
      next();
    }
    #For every line that contains a probe ID, process and output
    $probe_count++;

    #Determine probe type from file to reduce queries to database
    #Check the probe type and number of exons skipped
    #Skip all probes except exon and exon-exon with exons_skipped=0
    my $p_type = $line[3];
    my $e_skipped = $line[6];

    $probes{$probe_id}{count} = $probe_count;
    $probes{$probe_id}{line} = $_;
    $probes{$probe_id}{type} = $p_type;
    $probes{$probe_id}{exons_skipped} = $e_skipped;
  }
  close (PROBES);
  print "\n\nFound $probe_count probes in the input file: $probe_file\n\n";
  return();
}



###############################################################################################################
#getDistances                                                                                                 #
###############################################################################################################
sub getDistances{
  my %args = @_;
  my $dbh = $args{'-dbh'};

  print "\nRetrieving data from ALEXA\n\n";

  my $exon_canonical_probes = 0;
  foreach my $probe_id (sort keys %probes){

    #Check the probe type and number of exons skipped
    #Skip all probes except exon and exon-exon with exons_skipped=0
    my $probe_type = $probes{$probe_id}{type};
    my $exons_skipped = $probes{$probe_id}{exons_skipped};

    unless ($probe_type eq "Exon" || $probe_type eq "Exon-Exon"){
      $probes{$probe_id}{dist_3_prime} = "na";
      $probes{$probe_id}{dist_5_prime} = "na";
      $probes{$probe_id}{total_exon_bases} = "na";
      next();
    }
    if ($probe_type eq "Exon-Exon" && $exons_skipped > 0){
      $probes{$probe_id}{dist_3_prime} = "na";
      $probes{$probe_id}{dist_5_prime} = "na";
      $probes{$probe_id}{total_exon_bases} = "na";
      next();
    }
    $exon_canonical_probes++;

    my %probe_info = %{&getProbeInfo('-dbh'=>$alexa_dbh, '-probe_count_id'=>$probe_id)};

    my $gene_id = $probe_info{gene_id};
    my $position = $probe_info{unit1_start};
    my %gene_pos = %{&relativeGenePosition ('-dbh'=>$alexa_dbh, '-gene_id'=>$gene_id, '-position'=>$position)};

    my $distance_from_end = $gene_pos{$gene_id}{dist_3prime};
    my $distance_from_start = $gene_pos{$gene_id}{dist_5prime};
    my $total_exon_bases = $gene_pos{$gene_id}{total_exon_bases};

    $probes{$probe_id}{dist_3_prime} = $distance_from_end;
    $probes{$probe_id}{dist_5_prime} = $distance_from_start;
    $probes{$probe_id}{total_exon_bases} = $total_exon_bases;
    print "\n$probe_id\t$gene_id\t$distance_from_start\t$distance_from_end\t$total_exon_bases";
  }

  print "\nFound $exon_canonical_probes that are exon or canonical probes\n\n";

  return();
}


###############################################################################################################
#writeFile                                                                                                    #
###############################################################################################################
sub writeFile{
  my %args = @_;
  my $out_file = $args{'-out_file'};

  print "\nWriting output file: $out_file\n\n";

  open (OUT, ">$out_file") || die "\nCould not open output file: $out_file\n\n";

  print OUT "$header\tDist_5prime\tDist_3prime\tTranscriptSize\n";
  foreach my $probe_id (sort {$probes{$a}->{count} <=> $probes{$b}->{count}} keys %probes){
    print OUT "$probes{$probe_id}{line}\t$probes{$probe_id}{dist_5_prime}\t$probes{$probe_id}{dist_3_prime}\t$probes{$probe_id}{total_exon_bases}\n";
  }

  close (OUT);

  return();
}
