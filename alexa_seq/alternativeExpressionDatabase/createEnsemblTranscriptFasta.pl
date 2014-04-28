#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#This script takes all exons for each ensembl transcript, assembles them into the transcript sequence and creates
#a fasta file for all transcripts for every ensembl gene stored in ALEXA.  This file will be used to create a blast database
#which will then be used to test the specificity of probes for the target gene locus
use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::ALEXA_DB qw(:all);
use utilities::utility qw(:all);

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $outfile = '';
my $logfile = '';
my $chr = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'outfile=s'=>\$outfile, 'logfile=s'=>\$logfile, 'chr=s'=>\$chr);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password\n", RESET;
print GREEN, "\n\tSpecify an output fasta file using: --outfile", RESET;
print GREEN, "\n\tSpecify a logfile for output using: --logfile", RESET;
print GREEN, "\n\tOptional, use --chr to limit to transcripts from a particular chromosome.  e.g.  --chr=1", RESET;
print GREEN, "\n\nExample: createEnsemblTranscriptFasta.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --outfile=/home/user/alexa/ALEXA_version/data_sources/ensembl_transcripts/alexa_transcripts.fa  --logfile=/home/user/alexa/ALEXA_version/logs/specificity/createEnsemblTranscriptsFasta_LOG.txt\n\n", RESET;

unless ($database && $server && $user && $password && $outfile && $logfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Open logfile for output
open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";
print LOG "\nUser Specified the following options:\ndatabase = $database\noutfile = $outfile\nlogfile = $logfile\n\n";

open (OUTFILE, ">$outfile") || die "\nCould not open output fasta file $outfile\n\n";

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

#Get all ensembl genes
my @gene_ids;
if ($chr){
  @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All', '-chromosome' => $chr)}
}else{
  @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};
}

my $genes_found = @gene_ids;
my $transcripts_found = 0;

print BLUE, "\n\nFound $genes_found genes in $database\n\n", RESET;
print LOG "\n\nFound $genes_found genes in $database\n\n";

my $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"yes");
my $gene_count = 0;

foreach my $gene_id (@gene_ids){
  $gene_count++;

  #Get a transcripts object for this gene
  my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};

  #Go through each transcript, build the transcript sequence and dump it to the fasta file
  foreach my $trans_id (sort {$a <=> $b} keys %{$transcripts_ref}){
    my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

    my $transcript_seq = '';
    #Sort the exons for this transcript by their start position
    foreach my $exon_id (sort {$exons_ref->{$a}->{exon_start} <=> $exons_ref->{$b}->{exon_start}} keys %{$exons_ref}){
      $transcript_seq = "$transcript_seq"."$exons_ref->{$exon_id}->{sequence}"
    }
    $transcripts_found++;
    print OUTFILE ">$trans_id\n$transcript_seq\n";
  }
}
close (OUTFILE);

print BLUE, "\n\nPrinted fasta records for $transcripts_found transcripts to $outfile\n\n", RESET;
print LOG "\n\nPrinted fasta records for $transcripts_found transcripts to $outfile\n\n";

close (LOG);

#Close database connection
$alexa_dbh->disconnect();

exit();

