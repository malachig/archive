#!/usr/bin/perl -w
#Written by Malachi Griffith

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
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $gene_type = '';
my $verbose = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'gene_type=s'=>\$gene_type, 'verbose=s'=>\$verbose);

print GREEN, "\n\nThis script determines the size of a transcriptome (bp occupied by the exons of all transcripts of all genes)", RESET;
print GREEN, "\nOverlaping bases in overlaping exons are only counted once", RESET;
print GREEN, "\nThe user can limit the analysis to only a specific type of EnsEMBL gene\n\n", RESET;
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the gene type using: --gene_type=protein_coding (or 'All' or 'miRNA' etc)", RESET;
print GREEN, "\n\tUse --verbose=yes for detailed output", RESET;
print GREEN, "\n\nExample: determineTranscriptomeSize.pl  --database=ALEXA_hs_45_36g  --server=server_name  --user=username  --password=pwd  --gene_type=protein_coding  --verbose=no\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $gene_type && $verbose){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

#Get all the genes of the type specified by the user
my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>$gene_type)};
my $gene_count = @gene_ids;
print BLUE, "\nFound $gene_count genes of the type: $gene_type\n\n", RESET;

#Get non-redundant exon content for all of the exons of all of the transcripts of each gene
#Note that every base expressed from a gene locus is only counted once by this approach (exons are merged)
my $genes_ref = &getExonContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);

#Calculate the transcriptome size by adding up all of the merged exon regions of every gene
my $transcriptome_size = 0;
foreach my $gene_id (keys %{$genes_ref}){

  my $gene_content_size = 0;
  my $exonContent_ref = $genes_ref->{$gene_id}->{exon_content}; #Get exon content for each gene

  foreach my $exon_region (keys %{$exonContent_ref}){
    my $exon_region_size = ($exonContent_ref->{$exon_region}->{end} - $exonContent_ref->{$exon_region}->{start})+1;
    $gene_content_size += $exon_region_size;
  }
  if ($verbose eq "yes"){
    print BLUE, "\nThe exons of the transcripts of gene: $gene_id consist of $gene_content_size bp", RESET;
  }
  $transcriptome_size += $gene_content_size;
}

#Close database connection
$alexa_dbh->disconnect();

print BLUE, "\n\nFound $gene_count genes of the type: $gene_type", RESET;
print BLUE, "\nTotal transcriptome size of the $gene_type genes is $transcriptome_size\n\n", RESET;

exit();
