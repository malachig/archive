#!/usr/bin/perl -w
#Written by Malachi Griffith

#The purpose of this script is to consider all Ensembl genes which have multiple transcripts and identify exon skipping events
#Hopefully this will provide some sense of the normal amount of exons that get skipped.  In other words, is it reasonable to skip 1,2,3 ..10?
#This information can be used to decide how many probes to include.  For example, there may be no point in including a probe that interrogates the possibility
#of exon 1 being connected to exon 50. The existing distribution of exon skipping events in EnsEMBL can be used to chose a rationale cutoff point

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
my $allow_predicted_genes = '';
my $logfile = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'target_dir=s'=>\$target_dir, 
	    'allow_predicted_genes=s'=>\$allow_predicted_genes, 'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify a target directory for output files using: --target_dir", RESET;
print GREEN, "\n\tIf you wish to allow predicted genes (for species with few known genes), use: --allow_predicted_genes=yes", RESET;
print GREEN, "\n\tSpecify a logfile for output using: --logfile", RESET;
print GREEN, "\n\nExample: summarizeEnsemblExonSkipping.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=mysql_user  --password=mysql_pwd  --target_dir=/home/user/alexa/ALEXA_version/stats/genes  --allow_predicted_genes=no  --logfile=/home/user/alexa/ALEXA_version/logs/summarizeEnsemblExonSkipping_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $target_dir && $allow_predicted_genes && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";
#Print out the parameters supplied by the user to the logfile for future reference
print LOG "\nUser Specified the following options:\ndatabase = $database\ntarget_dir = $target_dir\nallow_predicted_genes = $allow_predicted_genes\nlogfile = $logfile\n\n";

#Establish connection with the Alternative Splicing Expression database using details provided by user
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

#First get all the protein coding genes in ALEXA
my @gene_ids;

if ($allow_predicted_genes eq "yes"){
  @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'Non-pseudo')};
}else{
  @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'Non-pseudo', '-evidence'=>"Known Gene")};
}

my $multi_transcript_genes = 0;
my $multi_transcript_count = 0;
my $exon_skipping_events = 0;
my @exon_skip_counts;

#For each of these genes, get the transcripts and exons
my $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

#Print out the parameters supplied by the user to the logfile for future reference
print LOG "\nUser Specified the following options:\ndatabase = $database\ntarget_dir = $target_dir\nlogfile = $logfile\n\n";

foreach my $gene_id (@gene_ids){

  my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};

  #Unless there are multiple transcripts, skip this gene
  my $trans_count = keys %{$transcripts_ref};
  unless ($trans_count > 1){
    next();
  }
  $multi_transcript_genes++;
  $multi_transcript_count += $trans_count;

  print BLUE, "\n\nGENE: $gene_id", RESET;
  print LOG "\n\nGENE: $gene_id";

  #1.) Assemble a reference set of exons (a superset of all non-redundant exons)
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

  #2.) Get arrays to represent the reference exons
  my $ref_exon_count = keys %reference_exons;

  my @reference_exon_starts;
  my @reference_exon_ends;
  foreach my $ref_exon_id (sort {$reference_exons{$a}->{exon_start} <=> $reference_exons{$b}->{exon_start}} keys %reference_exons){
    push (@reference_exon_starts, $reference_exons{$ref_exon_id}{exon_start});
    push (@reference_exon_ends, $reference_exons{$ref_exon_id}{exon_end});
  }
  print BLUE, "\nREF EXON STARTS: @reference_exon_starts", RESET;
  print BLUE, "\nREF EXON ENDS:   @reference_exon_ends", RESET;

  print LOG "\nREF EXON STARTS: @reference_exon_starts";
  print LOG "\nREF EXON ENDS:   @reference_exon_ends";


  #3.) Now go through each exon of each transcript, compare to the reference exons and identiy exon skipping events
  foreach my $trans_id (keys %{$transcripts_ref}){
    my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

    #Sort the exons for this transcript and consider create arrays of their start and end positions
    my @trans_exon_starts;
    my @trans_exon_ends;
    foreach my $trans_exon_id (sort {$exons_ref->{$a}->{exon_start} <=> $exons_ref->{$b}->{exon_start}} keys %{$exons_ref}){
      push (@trans_exon_starts, $exons_ref->{$trans_exon_id}->{exon_start});
      push (@trans_exon_ends, $exons_ref->{$trans_exon_id}->{exon_end});
    }
    #Determine the number of exons in this transcript
    my $trans_exon_count = keys %{$exons_ref};

    print BLUE, "\n\tTRANS: $trans_id ($trans_exon_count exons)", RESET;
    print BLUE, "\n\t\tTRANS EXON STARTS: @trans_exon_starts", RESET;
    print BLUE, "\n\t\tTRANS EXON ENDS:   @trans_exon_ends", RESET;

    print LOG "\n\tTRANS: $trans_id ($trans_exon_count exons)";
    print LOG "\n\t\tTRANS EXON STARTS: @trans_exon_starts";
    print LOG "\n\t\tTRANS EXON ENDS:   @trans_exon_ends";

    #If this transcript has only 1 exon, skip it
    if ($trans_exon_count == 1){
      next();
    }

    #Now go through each exon, get the next exon in the transcript and see if these two exons are consecutive in the reference exon set
    TRANS_EXON: for (my $i = 0; $i < $trans_exon_count-1; $i++){
	my $trans_exon1_start = $trans_exon_starts[$i];
	my $trans_exon1_end = $trans_exon_ends[$i];
	my $trans_exon2_start = $trans_exon_starts[$i+1];
	my $trans_exon2_end = $trans_exon_ends[$i+1];

	#Now look for 'valid' skipped exons between exon1 and exon2 of this transcript in the reference exon set
	my $skipped_exons = 0;
	my $found_first_exon = 0;
	for (my $j = 0; $j <= $ref_exon_count-1; $j++){
	  #Find the first exon
	  if ($reference_exon_starts[$j] == $trans_exon1_start && $reference_exon_ends[$j] == $trans_exon1_end){
	    $found_first_exon = 1;
	    next(); #next exon in the reference exons
	  }
	  #If the second exon is found, stop searching for this transcript exon pair
	  if ($reference_exon_starts[$j] == $trans_exon2_start && $reference_exon_ends[$j] == $trans_exon2_end){
	    if ($skipped_exons > 0){
	      #If evidence for exon skipping was observed for this exon pair in this transcript, report it
	      $exon_skipping_events++;
	      push (@exon_skip_counts, $skipped_exons);
	      print BLUE, "\nObserved $skipped_exons skipped exons between ($trans_exon1_start - $trans_exon1_end) and ($trans_exon2_start - $trans_exon2_end)", RESET;
	      print LOG "\nObserved $skipped_exons skipped exons between ($trans_exon1_start - $trans_exon1_end) and ($trans_exon2_start - $trans_exon2_end)";
	    }

	    next TRANS_EXON;
	  }

	  #If the first trans exon is already found and the second hasn't been found yet, start counting skipped exons
	  #But only count ones that would be valid connections.  ie. the reference exon is completely within the coordinates of exon1 and exon2
	  #Otherwise it may be cases of alternate splice site usage, which is a different issue
	  if ($found_first_exon == 1 && $trans_exon1_end < $reference_exon_starts[$j] && $trans_exon2_start > $reference_exon_ends[$j]){
	    $skipped_exons++;
	    next();
	  }
	}
      }

  }
}

my $outfile = "$target_dir"."/"."exon_skip_counts.txt";
print BLUE, "\n\nSUMMARY:", RESET;
print BLUE, "\nObserved a total of $multi_transcript_genes genes with multiple transcripts", RESET;
print BLUE, "\nThese multi-transcript genes are comprised of a total of $multi_transcript_count transcripts", RESET;
print BLUE, "\nObserved a total of $exon_skipping_events exon skipping events", RESET;
print BLUE, "\nThe exon skip counts were written to $outfile\n\n", RESET;

print LOG "\n\nSUMMARY:";
print LOG "\nObserved a total of $multi_transcript_genes genes with multiple transcripts";
print LOG "\nThese multi-transcript genes are comprised of a total of $multi_transcript_count transcripts";
print LOG "\nObserved a total of $exon_skipping_events exon skipping events";
print LOG "\nThe exon skip counts were written to $outfile\n\n";

open (OUTFILE, ">$outfile") || die "\nCould not open output file $outfile\n\n";

foreach my $count (@exon_skip_counts){
  print OUTFILE "$count\n";
}
close (OUTFILE);

close (LOG);

#Close database connection
$alexa_dbh->disconnect();

exit();
