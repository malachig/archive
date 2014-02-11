#!/usr/bin/perl -w

#Written by Malachi Griffith
#The purpose of this script is to test the specificity of probes in an input file of probe records.
#In particular it tests for probes within the extracted probe set which are identical or nearly identical to each other.
#The other specificity tests look for probes that do not appear to be specific to mRNAs or ensembl transcripts of the target locus.
#Such probes may be unique but still not specific to the target locus (becuase probes may have not been designed for the non-target locus)
#This script is need because the other specificity test may miss non-specific probes, for example, those that come from repeated motifs within the target gene.
#Such probe are undesirable because you will not be able to tell which exon is being expressed in that gene.
#Note: It relies on the existence of BLAST results prepared in advance using megaBLAST (likely on the cluster)
# and stored in a results directory.  Before running this script, first use:
#   blastProbesVersus_Probes.pl to generate probe fasta files and then run the blast jobs that it produces on the cluster.
#This script will thus rely on a probe file from /home/user/alexa/ALEXA_version/specificity/job_name/probes and will look for the blast results for
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
#4.) Foreach 'Probe' get the 'Blast_results' found and identify cases where a significant hit occurs to probe that is not the same probe
#    Only self-self hits are acceptable (same probe, or probe of the same probeset).  Look for perfect matches or near perfect matches to other probes
#    - The result will be recorded as the number of significant hits to non self probes

#Note1: A hit to any probe other than probes of the same probeset will fail a probe sequence..
#Note2: The problem with this is that you have probesets consisting of an entire exon or intron.  It is possible that two probes from different parts of
#       an exon or intron could be very similar not becuase they are overlapping, but because there is a repeated sequence within the intron or exon...
#       - This will most likely be a rare event, but in theory could allow non-unique or nearly identical probes to be in the non-filtered probe list
#       - The introns are not expected to be expressed anyway, so this should not be a big issue
#       - For large exons with duplicated segments, such probes could be detecting multiple parts of the same target exon, resulting in a higher than expected
#         expression estimate.  The probability of such events is probably low enough to ignore this eventuality for now ...
#Note3: This script assumes that the blast results were generated in tabular format!!  (using the -m 8 option)

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my $probe_infile = '';
my $blast_results_infile = '';
my $outfile = '';

GetOptions ('probe_infile=s'=>\$probe_infile, 'blast_results_infile=s'=>\$blast_results_infile, 'outfile=s'=>\$outfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script test the specificity of probes against all experimental probes", RESET;
print GREEN, "\n\tThis script takes a probe file and blast results file as input, parses the results and creates and updated output file", RESET;
print GREEN, "\n\tSpecify a tab-delimited input probe file using: --probe_infile", RESET;
print GREEN, "\n\tSpecify a tab-delimited blast results file using: --blast_results_file", RESET;
print GREEN, "\n\tSpecify the name of the resulting outfile using: --outfile", RESET;
print GREEN, "\n\nExample: testProbeSpecificity_Probes.pl --probe_infile=exonJunctionProbes.txt --blast_results_infile=blast_results.txt --outfile=exonJunctionProbes_Spec_Probes.txt\n\n", RESET;

unless ($probe_infile && $blast_results_infile && $outfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#1.) Load the probes of interest which have already been blasted against the experimental Probes database
my $header_line;
my $probes_ref = &loadProbeInfo('-input_file'=>$probe_infile);

#2.) Load all BLAST results from the specified blast-results file
my $blast_results_ref = &loadBlastResults('-blast_results_file'=>$blast_results_infile);

#3.) Using the infomation gathered in steps 1-2, conduct the actual specificity test
&testProbeSpecificity('-probe_object'=>$probes_ref, '-blast_object'=>$blast_results_ref);

#4.) Print the probe info out to a new file and append the results of the specificity test
&printResults('-output_file'=>$outfile, '-probe_object'=>$probes_ref, '-header_line'=>$header_line);

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
    my $probeset_id = $line[1];
    my $gene_id = $line[2];

    if ($probe_object{$probe_id}){
      print RED, "\nNon unique Probe_id found! ($probe_id).  Check input file\n\n", RESET;
      exit();
    }
    chomp($line_record);
    $probe_object{$probe_id}{probeset_id} = $probeset_id;
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
  print BLUE, "\n\tProcessing $results_file", RESET;
  open (BLAST_RESULTS, "$results_file") || die "\nCould not open blast results file: $results_file\n\n";

  while (<BLAST_RESULTS>){
    my @line = split ("\t", $_);
    my $probe_id = $line[0];  #Query probe name
    my $hit_probe_name = $line[1]; #Hit probe name (format = '$probe_id_$probeset_id)'
    my $hit_probe_id;
    my $hit_probeset_id;
    if ($hit_probe_name =~ /^(\d+)_(\d+)/){
      $hit_probe_id = $1;
      $hit_probeset_id = $2;
    }else{
      print RED, "\nFormat of hit probe name: $hit_probe_name from blast results file is not understood\n\n";
      exit();
    }

    my $percent_identity = $line[2];
    my $alignment_length = $line[3];
    my $subject_start = $line[8]; #Position on the probe sequence at which the probe alignment begins
    my $subject_end = $line[9]; #Position on the probe sequence at which the probe alignment ends 

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
      my %hits = %{$blast_results_object{$probe_id}{blast_hits}};

      #See if there has already been a hit for this probe (ie. multiple hits of one probe to another probe)
      if ($hits{$hit_probe_id}){
	my $test_length = $hits{$hit_probe_id}{al};
	#print "\nMultiple hits found for a single probe to a single database sequence";

	#Replace the old hit record for this sequence if this hit is longer than the one already recorded
	if ($alignment_length > $test_length){
	  $hits{$hit_probe_id}{pi} = $percent_identity;
	  $hits{$hit_probe_id}{al} = $alignment_length;
	  $hits{$hit_probe_id}{ss} = $subject_start;
	  $hits{$hit_probe_id}{se} = $subject_end;
	  $blast_results_object{$probe_id}{blast_hits} = \%hits;
	  $blast_results_stored++;
	}else{
	  $less_significant_hits++;
	  next(); #Skip less significant hits to the same sequence
	}

      }else{
	#This hit is to a new accession ID
	$hits{$hit_probe_id}{pi} = $percent_identity;
	$hits{$hit_probe_id}{al} = $alignment_length;
	$hits{$hit_probe_id}{ss} = $subject_start;
	$hits{$hit_probe_id}{se} = $subject_end;
	$hits{$hit_probe_id}{probeset_id} = $hit_probeset_id;
	$blast_results_object{$probe_id}{blast_hits} = \%hits;
	$blast_results_stored++;
      }

    }else{
      #create the first hit record for this probe
      my %hits;
      $hits{$hit_probe_id}{pi} = $percent_identity;
      $hits{$hit_probe_id}{al} = $alignment_length;
      $hits{$hit_probe_id}{ss} = $subject_start;
      $hits{$hit_probe_id}{se} = $subject_end;
      $hits{$hit_probe_id}{probeset_id} = $hit_probeset_id;
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
#3.) Using the infomation gathered in steps 1-2, conduct the actual specificity test                    #
#########################################################################################################
sub testProbeSpecificity{
  my %args = @_;
  my $probe_object_ref = $args{'-probe_object'};
  my $blast_object_ref = $args{'-blast_object'};

  print BLUE, "\n\nBegin actual specificity test", RESET;

  #Go through each probe ID and test its specificity compared to all other probes
  foreach my $probe_id (sort {$a <=> $b} keys %{$probe_object_ref}){
    my $probeset_id = $probe_object_ref->{$probe_id}->{probeset_id};

    #Count the number of times a probe hits a transcript sequence from the target gene and also hits to non-target genes
    my $non_target_hits = 0;
    my $target_hits = 0;

    #Initialize variables
    $probe_object_ref->{$probe_id}->{non_target_hits} = $non_target_hits;
    $probe_object_ref->{$probe_id}->{target_hits} = $target_hits;

    #Keep track of the longest non-target and target alignment lengths for each probe
    $probe_object_ref->{$probe_id}->{largest_non_target_al} = 0;
    $probe_object_ref->{$probe_id}->{largest_target_al} = 0;

    my $probe_gene_id = $probe_object_ref->{$probe_id}->{alexa_gene_id};

    #Get the blast hits for this probe (as a hash keyed on the probe ID)
    #Check for probes with no blast hits
    my $blast_hits_ref;
    if ($blast_object_ref->{$probe_id}->{blast_hits}){
      $blast_hits_ref = $blast_object_ref->{$probe_id}->{blast_hits};
    }else{
      print RED, "\nDid not find any hits for this probe ($probe_id).  It should at least hit itself?!\n\n", RESET;
      next();
    }

    #Go through each BLAST hit and make sure it is to a probe of the same probeset and not a non-target probeset
    BLAST_HIT:foreach my $probe_hit_id (sort keys %{$blast_hits_ref}){
	#1.) First eliminate the cases where the probe_hit_id is not defined (it is not even one of the probes in the block input file)
	#    - This is most definitely a non-target hit
	unless ($probe_object_ref->{$probe_hit_id}){
	  $non_target_hits++;
	  #update the largest non-target hit found so far if neccessary
	  if ($blast_hits_ref->{$probe_hit_id}->{al} > $probe_object_ref->{$probe_id}->{largest_non_target_al}){
	    $probe_object_ref->{$probe_id}->{largest_non_target_al} = $blast_hits_ref->{$probe_hit_id}->{al};
	  }
	  next BLAST_HIT;
	}

	#2.) Is the probe ID for the hit the same as the probe ID of the current probe (ie. a self-self hit)?
	#    - Also check if the probeset ID of the hit is the same as the probeset ID of the current probe
	if (($probe_hit_id == $probe_id) || ($probeset_id == $blast_hits_ref->{$probe_hit_id}->{probeset_id})){
	  $target_hits++;

	  #update the largest target hit found so far if neccessary
	  if ($blast_hits_ref->{$probe_hit_id}->{al} > $probe_object_ref->{$probe_id}->{largest_target_al}){
	    $probe_object_ref->{$probe_id}->{largest_target_al} = $blast_hits_ref->{$probe_hit_id}->{al};
	  }
	  next BLAST_HIT;
	}

	#3.) If the probe ID for this hit is from the block file, but is NOT a self-self hit or the same probeset
	#    - Then this must also be a non-target hit
	$non_target_hits++;

	#update the largest non-target hit found so far if neccessary
	if ($blast_hits_ref->{$probe_hit_id}->{al} > $probe_object_ref->{$probe_id}->{largest_non_target_al}){
	  $probe_object_ref->{$probe_id}->{largest_non_target_al} = $blast_hits_ref->{$probe_hit_id}->{al};
	}
	next BLAST_HIT;

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

  print "\n\nBegin printing results to $output_file\n";

  my $new_column1 = "probe_Non-TargetHits";
  my $new_column2 = "probe_TargetHits";
  my $new_column3 = "probe_largestNon-TargetHitLength";
  my $new_column4 = "probe_largestTargetHitLength";

  open (PROBE_OUT, ">$output_file") || die "\nCould not open output file: $output_file";

  print PROBE_OUT "Probe_Count\t$new_column1\t$new_column2\t$new_column3\t$new_column4\n";

  foreach my $probe_id (sort {$a <=> $b} keys %{$probe_object_ref}){

    #check the format of each variable - final sanity check
    unless ($probe_object_ref->{$probe_id}->{non_target_hits} =~ /\d+/){
      print "\nUndefined non-target_hits";
      exit();
    }
    unless ($probe_object_ref->{$probe_id}->{target_hits} =~ /\d+/){
      print "\nUndefined target_hits";
      exit();
    }
    unless ($probe_object_ref->{$probe_id}->{largest_non_target_al} =~ /\d+/){
      print "\nUndefined largest_non_target_al";
      exit();
    }
    unless ($probe_object_ref->{$probe_id}->{largest_target_al} =~ /\d+/){
      print "\nUndefined largest_target_al";
      exit();
    }

    print PROBE_OUT "$probe_id\t$probe_object_ref->{$probe_id}->{non_target_hits}\t$probe_object_ref->{$probe_id}->{target_hits}\t$probe_object_ref->{$probe_id}->{largest_non_target_al}\t$probe_object_ref->{$probe_id}->{largest_target_al}\n";
  }

  close PROBE_OUT;

  return();
}
