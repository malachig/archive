#!/usr/bin/perl -w

#Written by Malachi Griffith
#The purpose of this script is to test the specificity of probes in an input file of probe records.  It relies on the existence of BLAST results 
#prepared in advance using megaBLAST or blastall (likely on the cluster) and stored in a results directory.  Before running this script, first use:
#   blastProbesVersus_ensembl.pl to generate probe fasta files and then run the blast jobs that it produces on the cluster.
#This script will thus rely on a probe file from  /home/user/alexa/ALEXA_version/specificity/job_name/probes and will look for the blast results for
#these probes in /home/user/alexa/ALEXA_version/specificity/job_name/blast_results

#Script process:
#1.) Load the probes of interest which have already been blasted against all ensembl transcripts.  Once the test of specificity has been completed
#    the results will be appended to this file and written to a new output file
#    Create a 'Probes' hash (keyed on probe_ID) containing all info from the input file
#    Confirm that no probe ID is encountered more than once
#2.) Load all BLAST results from the blast-result file specified.
#    Create a 'Blast_results' hash (keyed on probe ID) which contains references to arrays or hashes of blast hits
#    Only store the minimum amount of info in the hash to save memory
#3.) Load map information from ALEXA describing all the transcripts for each gene and store in a hash
#4.) Foreach 'Probe' get the 'Blast_results' found and identify cases where a significant hit occurs to a transcript from a gene other than the 
#    target gene for that probe.  This will be noted as a failure.
#    - The result will be recorded as the number of significant hits to the transcripts from non-target genes.

#Note1: A hit on either strand to a sequence that maps outside of the targeted genomic region will fail a probe sequence.  Perhaps this is too conservative.
#-Actually, if you are hybridizing double-stranded cDNA samples to the array, this is a really good idea!

#Note2: Since the strand of the BLAST hit is not considered in this script, if a probe hits a transcript from the opposite strand of the target gene,
#       even though it is in the same genomic location, it will still be associated with a different gene ID and this hit will be marked as a non-target
#       hit.  This should effectively allow for the filtering of probes that correspond to antisense regions of the genome.

#Note3: This script assumes that the blast results were generated in tabular format!!  (using the -m 8 option)

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
my $probe_infile = '';
my $blast_results_infile = '';
my $outfile = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'probe_infile=s'=>\$probe_infile, 'blast_results_infile=s'=>\$blast_results_infile, 'outfile=s'=>\$outfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script test the specificity of probes against all Ensembl transcripts in the specified ALEXA database", RESET;
print GREEN, "\n\tThis script takes a probe file and blast results file as input, parses the results and creates and updated output file", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password\n", RESET;
print GREEN, "\n\tSpecify a tab-delimited input probe file using: --probe_infile", RESET;
print GREEN, "\n\tSpecify a tab-delimited blast results file using: --blast_results_file", RESET;
print GREEN, "\n\tSpecify the name of the resulting outfile using: --outfile", RESET;
print GREEN, "\n\nExample: testProbeSpecificity_ensembl.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --probe_infile=exonJunctionProbes.txt  --blast_results_infile=blast_results.txt  --outfile=exonJunctionProbes_SpecEnsembl.txt\n\n", RESET;

unless ($database && $server && $user && $password && $probe_infile && $blast_results_infile && $outfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);


#1.) Load the probes of interest which have already been blasted against mRNA and EST databases
my $header_line;
my $probes_ref = &loadProbeInfo('-input_file'=>$probe_infile);

#2.) Load all BLAST results from the specified blast-results file
my $blast_results_ref = &loadBlastResults('-blast_results_file'=>$blast_results_infile);

#3.) Load map information from ALEXA describing all the transcripts for each gene and store in a hash (key on transcript id)
my $transcripts_ref = &loadTranscriptInfo('-dbh'=>$alexa_dbh);

#4.) Using the infomation gathered in steps 1-3, conduct the actual specificity test
&testProbeSpecificity('-probe_object'=>$probes_ref, '-blast_object'=>$blast_results_ref, '-transcript_object'=>$transcripts_ref);

#5.) Print the probe info out to a new file and append the results of the specificity test
&printResults('-output_file'=>$outfile, '-probe_object'=>$probes_ref, '-header_line'=>$header_line);

#Close database connection
$alexa_dbh->disconnect();

exit();


#########################################################################################################
#1.) Load the probes of interest which have already been blasted against mRNA and EST databases         #
#########################################################################################################
sub loadProbeInfo{
  my %args = @_;
  my $probe_file = $args{'-input_file'};

  my %probe_object;
  my $probe_count = 0;

  #Open the probe file and read the neccessary probe data into a hash keyed on probe ID
  open (PROBES, "$probe_file") || die "\nCould not open input probe file: $probe_file\n\n";

  my $first_line = 1;
  while (<PROBES>){
    if ($first_line == 1){
      $header_line = $_;
      chomp ($header_line);
      $first_line = 0;
      next();
    }

    my $line_record = $_;
    my @line = split("\t", $_);
    my $probe_id = $line[0];
    my $gene_id = $line[2];

    if ($probe_object{$probe_id}){
      print RED, "\nNon unique Probe_id found! ($probe_id).  Check input file\n\n", RESET;
      exit();
    }
    chomp($line_record);
    $probe_object{$probe_id}{alexa_gene_id} = $gene_id;
    $probe_object{$probe_id}{line_record} = $line_record;
    $probe_count++;
  }

  close (PROBES);

  print BLUE, "\n\nFound $probe_count probes in the input file", RESET;
  return (\%probe_object);
}


#########################################################################################################
#2.) Load all BLAST results from the blast-results files                                                #
#########################################################################################################
sub loadBlastResults{
  my %args = @_;
  my $results_file = $args{'-blast_results_file'};

  my %blast_results_object;
  my $blast_results_processed = 0;
  my $blast_results_stored = 0;
  my $less_significant_hits = 0;

  #open the file provided
  chomp ($results_file);
  print BLUE, "\n\nProcessing $results_file", RESET;
  open (BLAST_RESULTS, "$results_file") || die "\nCould not open blast results file: $results_file\n\n";

  while (<BLAST_RESULTS>){
    my @line = split ("\t", $_);
    my $probe_id = $line[0];
    my $hit_accession = $line[1];
    my $percent_identity = $line[2];
    my $alignment_length = $line[3];
    my $subject_start = $line[8]; #Position on the mRNA/EST sequence at which the probe alignment begins
    my $subject_end = $line[9]; #Position on the mRNA/EST sequence at which the probe alignment ends

    #NOTE: The query sequence start and end positions are always written with the start less than the end.
    #However, the subject sequences can have the start position greater than the end, if the alignment is on the opposite DNA strand
    # from the way the subject sequence was originally put into the database.
    # - Note that although they are listed in reverse order they are still relative to the start of the subject sequence on the top strand!
    # - So you should swap them before doing coordinate tests below
    my $original_start;
    if ($subject_start > $subject_end){
      $original_start = $subject_start;
      $subject_start = $subject_end;
      $subject_end = $original_start;
    }

    $blast_results_processed++;

    #Note: using the following datastructure will not allow multiple hits by one probe to the same accession to be recorded
    #Only the more significant hit (longer) will be stored when this this occurs

    #See if there are any hits recorded for this probe already
    if ($blast_results_object{$probe_id}){
      my $hits_ref = $blast_results_object{$probe_id}{blast_hits};

      #See if there has already been a hit for this accession (ie. multiple hits of one probe to the same sequence)
      if ($hits_ref->{$hit_accession}){
	my $test_length = $hits_ref->{$hit_accession}->{al};
	#print "\nMultiple hits found for a single probe to a single database sequence";

	#Replace the old hit record for this sequence if this hit is longer than the one already recorded
	if ($alignment_length > $test_length){
	  $hits_ref->{$hit_accession}->{pi} = $percent_identity;
	  $hits_ref->{$hit_accession}->{al} = $alignment_length;
	  $hits_ref->{$hit_accession}->{ss} = $subject_start;
	  $hits_ref->{$hit_accession}->{se} = $subject_end;
	  $blast_results_stored++;
	}else{
	  $less_significant_hits++;
	  next(); #Skip less significant hits to the same sequence
	}

      }else{
	#This hit is to a new accession ID
	$hits_ref->{$hit_accession}->{pi} = $percent_identity;
	$hits_ref->{$hit_accession}->{al} = $alignment_length;
	$hits_ref->{$hit_accession}->{ss} = $subject_start;
	$hits_ref->{$hit_accession}->{se} = $subject_end;
	$blast_results_stored++;
      }

    }else{
      #create the first hit record for this probe
      my %hits;
      $hits{$hit_accession}{pi} = $percent_identity;
      $hits{$hit_accession}{al} = $alignment_length;
      $hits{$hit_accession}{ss} = $subject_start;
      $hits{$hit_accession}{se} = $subject_end;
      $blast_results_object{$probe_id}{blast_hits} = \%hits;
      $blast_results_stored++;
    }
  }
  close BLAST_RESULTS;

  print BLUE, "\n\tProcessed $blast_results_processed blast results", RESET;
  print BLUE, "\n\tFound $less_significant_hits that were of the same or shorter length to another hit on the same sequence", RESET;
  print BLUE, "\n\tActually stored $blast_results_stored blast results", RESET;

  return(\%blast_results_object);
}


#########################################################################################################
#3.)  Load map information from ALEXA describing all the transcripts for each gene and store in a hash  #
#     (key on transcript id)                                                                            #
#########################################################################################################
sub loadTranscriptInfo{
  my %args = @_;

  my $dbh = $args{'-dbh'};

  my %trans_object;

  #Get gene/transcript mappings from ALEXA for all genes
  my $sql = "SELECT id,fk_Gene__id FROM Transcript";
  my $sth = $dbh->prepare("$sql");
  $sth->execute();

  while (my ($trans_id,$gene_id) = $sth->fetchrow_array()){
    $trans_object{$trans_id}{gene_id} = $gene_id;
  }
  $sth->finish();

  my $transcripts_found = keys %trans_object;
  print BLUE, "\n\nImported $transcripts_found Ensembl Transcripts and their corresponding gene IDs from ALEXA", RESET;

  return(\%trans_object);
}


#########################################################################################################
#4.) Using the infomation gathered in steps 1-3, conduct the actual specificity test                    #
#########################################################################################################
sub testProbeSpecificity{
  my %args = @_;
  my $probe_object_ref = $args{'-probe_object'};
  my $blast_object_ref = $args{'-blast_object'};
  my $trans_object_ref = $args{'-transcript_object'};

  print BLUE, "\n\nBegin actual specificity test", RESET;

  #Go through each probe ID and test its specificity to the target locus
  foreach my $probe_id (sort {$a <=> $b} keys %{$probe_object_ref}){

    #Count the number of times a probe hits a transcript sequence from the target gene and also hits to non-target genes
    my $non_target_hits = 0;
    my $target_hits = 0;

    #Initialize variables
    $probe_object_ref->{$probe_id}->{non_target_hits} = $non_target_hits;
    $probe_object_ref->{$probe_id}->{target_hits} = $target_hits;

    #Keep track of the longest non-target and target alignment lengths for each probe
    $probe_object_ref->{$probe_id}->{largest_non_target_al} = 0;
    $probe_object_ref->{$probe_id}->{largest_target_al} = 0;

    #Get the gene ID for this probe
    my $probe_gene_id = $probe_object_ref->{$probe_id}->{alexa_gene_id};

    #Get the blast hits for this probe (as a hash keyed on ALEXA transcript ID)
    #Check for probes with no blast hits
    my $blast_hits_ref;
    if ($blast_object_ref->{$probe_id}->{blast_hits}){
      $blast_hits_ref = $blast_object_ref->{$probe_id}->{blast_hits};
    }else{
      next();
    }

    #Go through each BLAST hit and make sure it is to a transcript of the target gene and not a non-target gene
    BLAST_HIT:foreach my $trans_hit_id (sort keys %{$blast_hits_ref}){

	#Make sure this transcript_id is in the hash containing transcript to gene mappings
	unless ($trans_object_ref->{$trans_hit_id}){
	  print RED, "\nThe transcript ID for this blast hit could not be found in the map hash from $database!\n\n", RESET;
	  exit();
	}

	#Is the gene_id associated with this transcript hit the same as the gene targeted by the current probe?
	if ($trans_object_ref->{$trans_hit_id}->{gene_id} == $probe_gene_id){
	  $target_hits++;

	  #update the largest target hit found so far if neccessary
	  if ($blast_hits_ref->{$trans_hit_id}->{al} > $probe_object_ref->{$probe_id}->{largest_target_al}){
	    $probe_object_ref->{$probe_id}->{largest_target_al} = $blast_hits_ref->{$trans_hit_id}->{al};
	  }
	}else{
	  $non_target_hits++;

	  #update the largest non-target hit found so far if neccessary
	  if ($blast_hits_ref->{$trans_hit_id}->{al} > $probe_object_ref->{$probe_id}->{largest_non_target_al}){
	    $probe_object_ref->{$probe_id}->{largest_non_target_al} = $blast_hits_ref->{$trans_hit_id}->{al};
	  }
	}
      }
    $probe_object_ref->{$probe_id}->{non_target_hits} = $non_target_hits;
    $probe_object_ref->{$probe_id}->{target_hits} = $target_hits;
  }
  return();
}


#########################################################################################################
#5.) Print the probe info out to a new file and append the results of the specificity test              #
#########################################################################################################
sub printResults{
  my %args = @_;

  my $output_file = $args{'-output_file'};
  my $probe_object_ref = $args{'-probe_object'};
  my $header_line = $args{'-header_line'};

  print BLUE, "\n\nBegin printing results to $output_file\n", RESET;

  my $new_column1 = "enst_Non-TargetHits";
  my $new_column2 = "enst_TargetHits";
  my $new_column3 = "enst_largestNon-TargetHitLength";
  my $new_column4 = "enst_largestTargetHitLength";

  open (PROBE_OUT, ">$output_file") || die "\nCould not open output file: $output_file";

  print PROBE_OUT "Probe_Count\t$new_column1\t$new_column2\t$new_column3\t$new_column4\n";

  foreach my $probe_id (sort {$a <=> $b} keys %{$probe_object_ref}){

    #check the format of each variable - final sanity check
    unless ($probe_object_ref->{$probe_id}->{non_target_hits} =~ /\d+/){
      print RED, "\nUndefined non-target_hits", RESET;
      exit();
    }
    unless ($probe_object_ref->{$probe_id}->{target_hits} =~ /\d+/){
      print RED, "\nUndefined target_hits", RESET;
      exit();
    }
    unless ($probe_object_ref->{$probe_id}->{largest_non_target_al} =~ /\d+/){
      print RED, "\nUndefined largest_non_target_al", RESET;
      exit();
    }
    unless ($probe_object_ref->{$probe_id}->{largest_target_al} =~ /\d+/){
      print RED, "\nUndefined largest_target_al", RESET;
      exit();
    }

    print PROBE_OUT "$probe_id\t$probe_object_ref->{$probe_id}->{non_target_hits}\t$probe_object_ref->{$probe_id}->{target_hits}\t$probe_object_ref->{$probe_id}->{largest_non_target_al}\t$probe_object_ref->{$probe_id}->{largest_target_al}\n";
  }

  close PROBE_OUT;

  return();
}
