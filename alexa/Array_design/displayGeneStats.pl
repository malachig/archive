#!/usr/bin/perl -w
#Written by Malachi Griffith
#The purpose of this script is to provide simple summaries of genes in an ALEXA database.
#   For each gene the following data will be displayed:
#   - alexa_id, ensembl_id, gene_type, gene_name, evidence, chromosome, coordinates, size, external IDs, transcripts, and exons

#It can also be used to generate summary output files to be used for generating statistics for all genes in an ALEXA database
#   This will allow a summary of the size of exons and transcripts as well as the number of transcripts per gene and exons per transcript

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
my $target_dir = '';
my $alexa_id = '';
my $ensembl_id = '';
my $all_genes = '';
my $summarize = '';
my $allow_predicted_genes = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'target_dir=s'=>\$target_dir,
	    'alexa_id=i'=>\$alexa_id, 'ensembl_id=s'=>\$ensembl_id, 'all_genes=i'=>\$all_genes, 'summarize=i'=>\$summarize,
	    'allow_predicted_genes=s'=>\$allow_predicted_genes);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify a target directory for output files using: --target_dir", RESET;
print GREEN, "\n\tIf you want to display info for a single Alexa Gene ID, use: --alexa_id", RESET;
print GREEN, "\n\tIf you want to display info for a single Ensembl Gene ID, use: --ensembl_id", RESET;
print GREEN, "\n\tIf you want to display info for all genes in the database, use: --all_genes", RESET;
print GREEN, "\n\tIf you wish to generate summary files, use: --summarize", RESET;
print GREEN, "\n\tFinally if you wish to allow predicted genes (for species with few known genes), use: --allow_predicted_genes=yes", RESET;
print GREEN, "\n\nExample: displayGeneStats.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=mysql_user  --password=mysql_pwd  --target_dir=/home/user/alexa/ALEXA_version/stats/genes --alexa_id=1  --allow_predicted_genes=no", RESET;
print GREEN, "\nExample: displayGeneStats.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=mysql_user  --password=mysql_pwd  --target_dir=/home/user/alexa/ALEXA_version/stats/genes  --ensembl_id=ENSG00000000003  --allow_predicted_genes=no", RESET;
print GREEN, "\nExample: displayGeneStats.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=mysql_user  --password=mysql_pwd  --target_dir=/home/user/alexa/ALEXA_version/stats/genes  --all_genes=1  --allow_predicted_genes=no", RESET;
print GREEN, "\nExample: displayGeneStats.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=mysql_user  --password=mysql_pwd  --target_dir=/home/user/alexa/ALEXA_version/stats/genes  --summarize=1  --allow_predicted_genes=no\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $target_dir && $allow_predicted_genes){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

unless ($alexa_id || $ensembl_id || $all_genes || $summarize){
  print RED, "\nChose one of the following options: (--alexa_id=1, --ensembl_id=ENSG00000000003, --all_genes=1, --summarize=1)\n\n", RESET;
  exit();
}

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);


#The user may test a single gene by ALEXA or Ensembl gene ID or get info on all of them
if ($alexa_id){

  &displayGeneStats ('-dbh'=>$alexa_dbh, '-alexa_id'=>$alexa_id);

}elsif ($ensembl_id){

  &displayGeneStats ('-dbh'=>$alexa_dbh, '-ensembl_id'=>$ensembl_id);

}elsif ($all_genes){

  my @gene_ids;
  if ($allow_predicted_genes eq "yes"){
    @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'Non-pseudo')};
  }else{
    @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'Non-pseudo', '-evidence'=>"Known Gene")};
  }

  #Display statistics for each gene
  foreach my $gene_id (@gene_ids){
    &displayGeneStats ('-dbh'=>$alexa_dbh, '-alexa_id'=>$gene_id);
  }

}elsif ($summarize){
  #Summarize the number of exons and transcripts per protein coding gene.  Also dump the transcript and exon sizes as columns to allow histograms to be made
  my @gene_ids;

  if ($allow_predicted_genes eq "yes"){
    @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'Non-pseudo')};
  }else{
    @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'Non-pseudo', '-evidence'=>"Known Gene")};
  }

  print BLUE, "\n\n", RESET;

  my $gene_count = @gene_ids;

  my $trans_file = "$target_dir"."/"."transcript_counts.txt";
  my $exon_file = "$target_dir"."/"."exon_counts.txt";
  my $trans_l_file = "$target_dir"."/"."transcript_lengths.txt";
  my $exon_l_file = "$target_dir"."/"."exon_lengths.txt";

  #Open output files for each of the data types to be summarized 
  open (TRANS, ">$trans_file") || die "\nCould not open transcript counts file: $trans_file";
  open (EXONS, ">$exon_file") || die "\nCould not open exon counts file: $exon_file";
  open (TRANS_LENGTH, ">$trans_l_file") || die "\nCould not open transcript lengths file: $trans_l_file";
  open (EXONS_LENGTH, ">$exon_l_file") || die "\nCould not open exon lengths file: $exon_l_file";

  my $gene_transcript_ref = &getTranscripts ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

  foreach my $gene_id (keys %{$gene_transcript_ref}){

    my $transcripts_ref = $gene_transcript_ref->{$gene_id}->{transcripts};

    my $number_transcripts = keys %{$transcripts_ref};

    print BLUE, "\nGENE: $gene_id has $number_transcripts transcripts", RESET;
    print TRANS "$number_transcripts\n";

    foreach my $trans_id (sort keys %{$transcripts_ref}){

      my $ensembl_t_id = $transcripts_ref->{$trans_id}->{ensembl_t_id};

      my $exons_ref = $transcripts_ref->{$trans_id}->{exons};
      my $number_exons = keys %{$exons_ref};
      print EXONS "$number_exons\n";

      my $transcript_size = 0;

      foreach my $exon_id (sort keys %{$exons_ref}){
	my $exon_size = ($exons_ref->{$exon_id}->{exon_end} - $exons_ref->{$exon_id}->{exon_start})+1;
	print EXONS_LENGTH "$exon_size\t$exons_ref->{$exon_id}->{ensembl_e_id}\t$ensembl_t_id\n";

	$transcript_size += $exon_size;
      }
      print BLUE, "\n\tTRANS: $ensembl_t_id is $transcript_size bp long and has $number_exons exons", RESET;
      print TRANS_LENGTH "$transcript_size\n";
    }
  }

  close TRANS;
  close EXONS;
  close TRANS_LENGTH;
  close EXONS_LENGTH;

}

#Close database connection
$alexa_dbh->disconnect();

exit();
