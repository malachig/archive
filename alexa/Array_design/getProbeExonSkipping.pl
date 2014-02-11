#!/usr/bin/perl -w
#Written by Malachi Griffith

#The purpose of this script is to analyze a file of exon junction probes and determine for each probe how many exons would be skipped
#in the theoretical exon-exon join it is attempting to interrogate
#Specifically it will look at each probe record in an input probe file and get the Exon IDs for Unit1 and Unit2.
#It will then determine the number of exons that are completely skipped by connecting Unit1 and Unit2 (a theoretical splice event)
#To do this, the script will require a hash of all genes and their exons (with coordinates).
#This information will then be added onto the end of each probe record.  The information will be used in the array design stage.
#By looking at the splicing events currently contained in Ensembl I have determined that 95% of all exon-skipping events involve a maximum
#of 5 skipped exons.  For this reason, it is probably safe to filter probes that interrogate skips of more than 6 exons.  This should allow
#more genes to be included on the array

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

#ALEXA libraries
#When a script is initiated, use the full path of the script location at execution to add the perl module libraries to @INC
#This should allow this scripts to work regardless of the current working directory or the script location (where it was unpacked).
#The /utilities directory must remain in the same directory as this script but the entire code directory can be moved around
BEGIN {
  my $script_dir = &File::Basename::dirname($0);
  push (@INC, $script_dir);
}
use utilities::ALEXA_DB qw(:all);
use utilities::utility qw(:all);

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $probe_file = '';
my $outfile = '';
my $logfile = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'probe_file=s'=>\$probe_file, 'outfile=s'=>\$outfile, 'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script determines the number of exons skipped by each exon-exon junction probe specified in an input file", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password\n", RESET;
print GREEN, "\n\tSpecify a tab-delimited input probe file using: --probe_file", RESET;
print GREEN, "\n\t\tIt is assume the the 1st column contains a unique probe ID", RESET;
print GREEN, "\n\tSpecify the name of the resulting outfile using: --outfile", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of this script: --logfile", RESET;
print GREEN, "\n\nExample: getProbeExonSkipping.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --probe_file=/home/user/alexa/ALEXA_version/unfiltered_probes/exonJunctionProbes.txt  --outfile=/home/user/alexa/ALEXA_version/unfiltered_probes/exonJunctionProbes_ExonsSkipped.txt  --logfile=/home/user/alexa/ALEXA_version/logs/getProbeExonSkipping_LOG.txt\n\n", RESET;

unless ($database && $server && $user && $password && $probe_file && $outfile && $logfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#Print out the parameters supplied by the user to the logfile for future reference
print LOG "\nUser Specified the following options:\ndatabase = $database\nprobe_file = $probe_file\noutfile = $outfile\nlogfile = $logfile\n\n";

#1.) First get all the probe IDs from the input file, as well as their Unit1 and Unit2 exon IDs
#Remember that in the case of overlapping exons, there may be multiple unit1 and Unit2 exons
my %probes;
my %master_gene_list;

print BLUE, "\nParsing probe file: $probe_file for probe data\n\n", RESET;
print LOG "\nParsing probe file: $probe_file for probe data\n\n";

open (PROBES, "$probe_file") || die "\nCould not open input probe file: $probe_file\n\n";

my $header = 1;
my %columns;
while (<PROBES>){
  chomp($_);
  my @line = split("\t", $_);

  if ($header == 1){
    my $column_count = 0;

    foreach my $column (@line){
      $columns{$column}{column_pos} = $column_count;
      $column_count++;
    }

    #Check for critical columns and their names
    unless ($columns{'Gene_ID'} && $columns{'Probe_Type'} && $columns{'Unit1_end'} && $columns{'Unit2_start'}){
      print RED, "\nCritical column missing or named incorrectly in probe file: $probe_file, check input file\n\n", RESET;
      close (LOG);
      exit();
    }

    $header = 0;
    next();
  }

  my $probe_id = $line[0];
  $master_gene_list{$line[2]}{tmp} = '';

  $probes{$probe_id}{gene_id} = $line[$columns{'Gene_ID'}{column_pos}];
  $probes{$probe_id}{probe_type} = $line[$columns{'Probe_Type'}{column_pos}];
  $probes{$probe_id}{unit1_end} = $line[$columns{'Unit1_end'}{column_pos}];
  $probes{$probe_id}{unit2_start} = $line[$columns{'Unit2_start'}{column_pos}];
}
close (PROBES);


#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

#First get all the protein coding genes in ALEXA
my @gene_ids = keys %master_gene_list;
my $gene_count = @gene_ids;
print BLUE, "\nFound probes for a total of $gene_count genes\n\n", RESET;
print LOG "\nFound probes for a total of $gene_count genes\n\n";

#Now for each of these genes, get the transcripts and exons
print BLUE, "\nExtracting Gene and Exon data from ALEXA\n\n", RESET;
print LOG "\nExtracting Gene and Exon data from ALEXA\n\n";
my $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

print BLUE, "\nGenerating reference exon objects\n\n", RESET;
print LOG "\nGenerating reference exon objects\n\n";

foreach my $gene_id (@gene_ids){

  my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};
  my $trans_count = keys %{$transcripts_ref};

  #2.) Assemble a reference set of exons (a superset of all non-redundant exons)
  my %reference_exons;
  foreach my $trans_id (keys %{$transcripts_ref}){
    my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

    foreach my $trans_exon_id (keys %{$exons_ref}){
      my $trans_exon_start = $exons_ref->{$trans_exon_id}->{exon_start};
      my $trans_exon_end = $exons_ref->{$trans_exon_id}->{exon_end};

      #Check each of the reference exons to see if one of these is the same, otherwise add it to the list
      my $redundant_exon = 0;
      foreach my $ref_exon_id (keys %reference_exons){
	if ($trans_exon_start == $reference_exons{$ref_exon_id}{exon_start} && $trans_exon_end == $reference_exons{$ref_exon_id}{exon_end}){
	  $redundant_exon = 1;
	}
      }
      #Unless the current transcript exon was found to be redundant, add it to the list
      unless ($redundant_exon == 1){
	$reference_exons{$trans_exon_id}{exon_start} = $trans_exon_start;
	$reference_exons{$trans_exon_id}{exon_end} = $trans_exon_end;
      }
    }
  }

  #3.) Get arrays to represent the reference exons
  my $ref_exon_count = keys %reference_exons;

  my @reference_exon_starts;
  my @reference_exon_ends;
  foreach my $ref_exon_id (sort {$reference_exons{$a}->{exon_start} <=> $reference_exons{$b}->{exon_start}} keys %reference_exons){
    push (@reference_exon_starts, $reference_exons{$ref_exon_id}{exon_start});
    push (@reference_exon_ends, $reference_exons{$ref_exon_id}{exon_end});
  }

  unless($ref_exon_count == 1){
    $gene_transcripts_ref->{$gene_id}->{ref_exon_starts} = \@reference_exon_starts;
    $gene_transcripts_ref->{$gene_id}->{ref_exon_ends} = \@reference_exon_ends;
    $gene_transcripts_ref->{$gene_id}->{ref_exon_count} = $ref_exon_count;
  }
}

#4.) Go through each probe and if it is an exon-exon probe, get the number of exons skipped for its theoretical target
print BLUE, "\nDetermining the number of exon skips for each exon junction probe\n\n", RESET;
print LOG "\nDetermining the number of exon skips for each exon junction probe\n\n";


my $counter = 0;
foreach my $probe_id (sort {$a <=> $b} keys %probes){

  #If this probe is not an exon-exon probe, set the skipped_exons value to 'na' and continue
  unless ($probes{$probe_id}{probe_type} eq "Exon-Exon"){
    $probes{$probe_id}{skipped_exons} = "na";
    next();
  }

  $counter++;
  if ($counter == 10000){
    print BLUE ".";
    $counter = 0;
  }

  my $gene_id = $probes{$probe_id}{gene_id};
  my $exon1_end = $probes{$probe_id}{unit1_end};
  my $exon2_start = $probes{$probe_id}{unit2_start};

  #Sanity check 
  unless ($exon2_start > $exon1_end){
    print RED, "\nExon1 end and Exon2 start positions do not make sense!\n\n", RESET;
    close (LOG);
    $alexa_dbh->disconnect();
    exit();
  }
  my $skipped_exons = 0;

  #Compare the coordinate of this probe to the reference exons and identify the number of skipped exons
  my @ref_exon_starts = @{$gene_transcripts_ref->{$gene_id}->{ref_exon_starts}};
  my @ref_exon_ends = @{$gene_transcripts_ref->{$gene_id}->{ref_exon_ends}};
  my $ref_exon_count = $gene_transcripts_ref->{$gene_id}->{ref_exon_count};

  #print "\n\nREF EXON STARTS: @ref_exon_starts";
  #print "\nREF EXON ENDS:   @ref_exon_ends";
  #print "\nPROBE: $probe_id\tExon1_end: $exon1_end\tExon2_start: $exon2_start";

  my $found_exon1 = 0;

  #First identify the last exon which corresponds to the exon1_end position
  for (my $i = 0; $i < $ref_exon_count-1; $i++){

    #Look for the exon which corresponds to the exon2_start position
    if ($exon2_start >= $ref_exon_starts[$i] && $exon2_start <= $ref_exon_ends[$i]){
      last();
    }

    #Look for the exon which corresponds to the exon1_end position
    if ($exon1_end >= $ref_exon_starts[$i] && $exon1_end <= $ref_exon_ends[$i]){
      $found_exon1 = 1;
      next();
    }

    #If exon1 has been found in the ordered reference set but exon2 has not been found yet, this is a skipped exon
    #The boundaries of the current reference exon should also be completely between the two exons targetted
    if ($found_exon1 == 1 && $exon1_end < $ref_exon_starts[$i] && $exon2_start > $ref_exon_starts[$i] && $exon1_end < $ref_exon_ends[$i] && $exon2_start > $ref_exon_ends[$i]){
      $skipped_exons++;
      next();
    }
  }
  $probes{$probe_id}{skipped_exons} = $skipped_exons;
}

#Finally print out the results
print BLUE, "\nPrinting output to file: $outfile\n\n", RESET;
print LOG "\nPrinting output to file: $outfile\n\n";

open (PROBES, "$probe_file") || die "\nCould not open input probe file: $probe_file\n\n";
open (OUTFILE, ">$outfile") || die "\nCould not open output probe file: $outfile\n\n";

$header = 1;
while (<PROBES>){
  if ($header == 1){
    $header = 0;
    chomp($_);
    print OUTFILE "$_\tExons_Skipped\n";
    next();
  }

  my @line = split("\t", $_);
  my $probe_id = $line[0];

  unless ($probes{$probe_id}){
    print RED, "\nProbe_id not found in hash!\n\n", RESET;
    close (LOG);
    $alexa_dbh->disconnect();
    exit();
  }
  chomp ($_);
  print OUTFILE "$_\t$probes{$probe_id}{skipped_exons}\n";

}
close (PROBES);
close (OUTFILE);
close (LOG);

#Close database connection
$alexa_dbh->disconnect();

exit();
