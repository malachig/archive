#!/usr/bin/perl -w

#Written by Malachi Griffith
#The purpose of this script is to test the specificity of probes in an input file of probe records.  It relies on the existence of BLAST results 
#prepared in advance using megaBLAST or BLAST (likely on the cluster) and stored in a results directory.  Before running this script, first use:
#   blastProbesVersus_Database.pl (with --type=genomic option) to generate probe fasta files and then run the blast jobs that it produces on the cluster.
#This script will thus rely on a probe file from  /home/user/alexa/ALEXA_version/specificity/job_name/probes and will look for the blast results for
#these probes in /home/user/alexa/ALEXA_version/specificity/job_name/blast_results

#Script process:
#1.) Load the probes of interest which have already been blasted against the complete genome.  Once the test of specificity has been completed
#    the results will be appended to this file and written to a new output file
#    Create a 'Probes' hash (keyed on probe_ID) containing all info from the input file
#    - Store the source genomic coordinates for each probe so that this can be used to identify target vs. non-target hits
#    Confirm that no probe ID is encountered more than once
#2.) Load all BLAST results from the blast-result file specified.
#    Create a 'Blast_results' hash (keyed on probe ID) which contains references to arrays or hashes of blast hits
#    Only store the minimum amount of info in the hash to save memory
#3.) Foreach 'Probe' get the 'Blast_results' found and identify cases where a significant hit occurs to region of the genome other than where 
#    the probe was extracted from.
#    - The result will be recorded as the number of significant hits to the target region and non-target regions of the complete genome.
#    - Also note the max hit size for both target and non-target hits
#    - Also note the number of non-target hits over 50% of probe length and over 75% of probe length


#Note1: A hit on either strand outside of the targeted genomic region will fail a probe sequence.

#Note2: This script is only set up to deal with probes that have a Unit1_start and Unit1_end position (region probes, exon probes, intron probes)
#       - It will require modification to work with junction or boundary probes

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

my $probe_infile = '';
my $blast_results_infile = '';
my $outfile = '';

GetOptions ('probe_infile=s'=>\$probe_infile, 'blast_results_infile=s'=>\$blast_results_infile, 'outfile=s'=>\$outfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script test the specificity of probes against the entire genome", RESET;
print GREEN, "\n\tThis script takes a probe file and blast results file as input, parses the results and creates and updated output file", RESET;
print GREEN, "\n\tSpecify a tab-delimited input probe file using: --probe_infile", RESET;
print GREEN, "\n\tSpecify a tab-delimited blast results file using: --blast_results_file", RESET;
print GREEN, "\n\tSpecify the name of the resulting outfile using: --outfile", RESET;
print GREEN, "\n\nExample: testProbeSpecificity_Genomic.pl  --probe_infile=regionProbes.txt  --blast_results_infile=blast_results.txt  --outfile=regionProbes_SpecGenomic.txt\n\n", RESET;

unless ($probe_infile && $blast_results_infile && $outfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#1.) Load the probes of interest which have already been blasted against mRNA and EST databases
my %columns;
my $header_line;
my $probes_ref = &loadProbeInfo('-input_file'=>$probe_infile);

#2.) Load all BLAST results from the specified blast-results file
my $blast_results_ref = &loadBlastResults('-blast_results_file'=>$blast_results_infile);

#3.) Using the infomation gathered in steps 1-3, conduct the actual specificity test
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
  my $column_count = 0;
  my $probe_count = 0;

  #Open the probe file and read the neccessary probe data into a hash keyed on probe ID
  open (PROBES, "$probe_file") || die "\nCould not open input probe file: $probe_file\n\n";

  my $first_line = 1;
  while (<PROBES>){
    if ($first_line == 1){
      $header_line = $_;
      chomp ($header_line);
      $first_line = 0;

      my @column_names = split ("\t", $header_line);
      foreach my $column (@column_names){
	$columns{$column}{position} = $column_count;
	$column_count++;
      }

      next();
    }

    my $line_record = $_;
    my @line = split("\t", $_);

    #Store the following info for each probe: Probe_Count, Probe_length, Chromosome, Unit1_start, Unit1_end
    my $probe_id = $line[$columns{'Probe_Count'}{position}];
    my $probe_length = $line[$columns{'Probe_length'}{position}];

    my $chromosome;
    #Watch for NC probes which do not have a target chromosome
    if ($columns{'Probe_Type'}{position}){
      my $probe_type = $line[$columns{'Probe_Type'}{position}];
      if ($probe_type eq "Control-Negative"){
	$chromosome = "na";
      }
    }else{
      $chromosome = $line[$columns{'Chromosome'}{position}];
    }
    my $probe_start = $line[$columns{'Unit1_start'}{position}];
    my $probe_end = $line[$columns{'Unit1_end'}{position}];

    if ($probe_object{$probe_id}){
      print RED, "\nNon unique Probe_id found! ($probe_id).  Check input file\n\n", RESET;
      exit();
    }
    chomp($line_record);
    $probe_object{$probe_id}{line_record} = $line_record;
    $probe_object{$probe_id}{probe_length} = $probe_length;
    $probe_object{$probe_id}{chromosome} = $chromosome;
    $probe_object{$probe_id}{probe_start} = $probe_start;
    $probe_object{$probe_id}{probe_end} = $probe_end;

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
  my $blast_hit_count = 0;

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

    #Sanity check
    if ($subject_end < $subject_start){
      print RED, "\nSubject start is larger than subject end! - makes no sense\n\n", RESET;
      exit();
    }

    $blast_hit_count++;

    #See if there are any hits recorded for this probe already
    if ($blast_results_object{$probe_id}){
      my $hits_ref = $blast_results_object{$probe_id}{blast_hits};

      #See if there has already been a hit for this accession (ie. multiple hits of one probe to the same sequence)
      #This hit is to a new accession ID (subject ID - in this case chromosome name)
      $hits_ref->{$blast_hit_count}->{sid} = $hit_accession;
      $hits_ref->{$blast_hit_count}->{pi} = $percent_identity;
      $hits_ref->{$blast_hit_count}->{al} = $alignment_length;
      $hits_ref->{$blast_hit_count}->{ss} = $subject_start;
      $hits_ref->{$blast_hit_count}->{se} = $subject_end;

    }else{
      #create the first hit record for this probe
      my %hits;
      $hits{$blast_hit_count}{sid} = $hit_accession;
      $hits{$blast_hit_count}{pi} = $percent_identity;
      $hits{$blast_hit_count}{al} = $alignment_length;
      $hits{$blast_hit_count}{ss} = $subject_start;
      $hits{$blast_hit_count}{se} = $subject_end;
      $blast_results_object{$probe_id}{blast_hits} = \%hits;
    }
  }
  close BLAST_RESULTS;

  print BLUE, "\n\tProcessed and stored $blast_hit_count blast results", RESET;

  return(\%blast_results_object);
}


#########################################################################################################
#3.) Using the infomation gathered in steps 1-3, conduct the actual specificity test                    #
#########################################################################################################
sub testProbeSpecificity{
  my %args = @_;
  my $probe_object_ref = $args{'-probe_object'};
  my $blast_object_ref = $args{'-blast_object'};

  print BLUE, "\n\nBegin actual specificity test", RESET;

  #Go through each probe ID and test its specificity to the target locus
  foreach my $probe_id (sort {$a <=> $b} keys %{$probe_object_ref}){

    #Count the number of times a probe hits a the target genomic sequence as well as non-target regions of the genome
    my $non_target_hits = 0;         #Non-target hits any size
    my $non_target_hits_over50 = 0;  #Non-target hits over 50% of the length of the probe
    my $non_target_hits_over75 = 0;  #Non-target hits over 75% of the length of the probe
    my $target_hits = 0;             #Hits to the region of the genome the probe actually came from

    my $probe_seq_name = $probe_object_ref->{$probe_id}->{chromosome};
    my $probe_start = $probe_object_ref->{$probe_id}->{probe_start};
    my $probe_end = $probe_object_ref->{$probe_id}->{probe_end};
    my $probe_length = $probe_object_ref->{$probe_id}->{probe_length};

    #Initialize variables
    $probe_object_ref->{$probe_id}->{non_target_hits} = $non_target_hits;
    $probe_object_ref->{$probe_id}->{non_target_hits_over50} = $non_target_hits_over50;
    $probe_object_ref->{$probe_id}->{non_target_hits_over75} = $non_target_hits_over75;
    $probe_object_ref->{$probe_id}->{target_hits} = $target_hits;

    #Keep track of the longest non-target and target alignment lengths for each probe
    $probe_object_ref->{$probe_id}->{largest_non_target_al} = 0;
    $probe_object_ref->{$probe_id}->{largest_non_target_pi} = 0;
    $probe_object_ref->{$probe_id}->{largest_target_al} = 0;

    #Get the blast hits for this probe (as a hash keyed on ALEXA transcript ID)
    #Check for probes with no blast hits
    my $blast_hits_ref;
    if ($blast_object_ref->{$probe_id}->{blast_hits}){
      $blast_hits_ref = $blast_object_ref->{$probe_id}->{blast_hits};
    }else{
      next();
    }

    #Go through each BLAST hit and determine whether it corresponds to the target genomic region of the probe or a non-target region
    BLAST_HIT:foreach my $blast_hit_id (sort keys %{$blast_hits_ref}){

	#Test for overlap of this hit's genomic coordinates with the coordinates actually targeted by the probe
	my $overlap = 0;

	my $hit_seq_name = $blast_hits_ref->{$blast_hit_id}->{sid};
	my $hit_start = $blast_hits_ref->{$blast_hit_id}->{ss};
	my $hit_end = $blast_hits_ref->{$blast_hit_id}->{se};

	#Check for hits to the target genome region - except for negative control probes which have no target region!
	unless ($probe_start eq "na" && $probe_end eq "na"){

	  if (($probe_seq_name eq $hit_seq_name) && ($hit_start >= $probe_start) && ($hit_start <= $probe_end)){
	    $overlap = 1;
	  }
	  if (($probe_seq_name eq $hit_seq_name) && ($hit_end >= $probe_start) && ($hit_end <= $probe_end)){
	    $overlap = 1;
	  }
	  if (($probe_seq_name eq $hit_seq_name) && ($hit_start <= $probe_start) && ($hit_end >= $probe_end)){
	    $overlap = 1;
	  }
	}

	if ($overlap == 1){
	  $target_hits++;

	  #update the largest target hit found so far if neccessary
	  if ($blast_hits_ref->{$blast_hit_id}->{al} > $probe_object_ref->{$probe_id}->{largest_target_al}){
	    $probe_object_ref->{$probe_id}->{largest_target_al} = $blast_hits_ref->{$blast_hit_id}->{al};
	  }
	}else{
	  $non_target_hits++;

	  #update the largest non-target hit alignment length found so far if neccessary
	  if ($blast_hits_ref->{$blast_hit_id}->{al} > $probe_object_ref->{$probe_id}->{largest_non_target_al}){
	    $probe_object_ref->{$probe_id}->{largest_non_target_al} = $blast_hits_ref->{$blast_hit_id}->{al};
	  }

	  #update the largest non-target percent identity found so far if neccessary
	  if ($blast_hits_ref->{$blast_hit_id}->{pi} > $probe_object_ref->{$probe_id}->{largest_non_target_pi}){
	    $probe_object_ref->{$probe_id}->{largest_non_target_pi} = $blast_hits_ref->{$blast_hit_id}->{pi};
	  }

	  #count non-target hits over 50 and 75% of probe length
	  my $percent_probe_length = ($blast_hits_ref->{$blast_hit_id}->{al} / $probe_length)*100;
	  if ($percent_probe_length > 50){
	    $non_target_hits_over50++;
	  }
	  if ($percent_probe_length > 75){
	    $non_target_hits_over75++;
	  }
	}
      }

    #Note the number of each type of hits found for this probe
    $probe_object_ref->{$probe_id}->{non_target_hits} = $non_target_hits;
    $probe_object_ref->{$probe_id}->{non_target_hits_over50} = $non_target_hits_over50;
    $probe_object_ref->{$probe_id}->{non_target_hits_over75} = $non_target_hits_over75;
    $probe_object_ref->{$probe_id}->{target_hits} = $target_hits;

  }
  return();
}


#########################################################################################################
#4.) Print the probe info out to a new file and append the results of the specificity test              #
#########################################################################################################
sub printResults{
  my %args = @_;

  my $output_file = $args{'-output_file'};
  my $probe_object_ref = $args{'-probe_object'};
  my $header_line = $args{'-header_line'};

  print BLUE, "\n\nBegin printing results to $output_file\n", RESET;

  my $new_column1 = "genomic_TargetHits";
  my $new_column2 = "genomic_largestTargetHitLength";
  my $new_column3 = "genomic_Non-TargetHits";
  my $new_column4 = "genomic_Non-TargetHits_over50";
  my $new_column5 = "genomic_Non-TargetHits_over75";
  my $new_column6 = "genomic_largestNon-TargetHitLength";
  my $new_column7 = "genomic_largestNon-TargetPercentIdentity";

  open (PROBE_OUT, ">$output_file") || die "\nCould not open output file: $output_file";

  print PROBE_OUT "Probe_Count\t$new_column1\t$new_column2\t$new_column3\t$new_column4\t$new_column5\t$new_column6\t$new_column7\n";

  foreach my $probe_id (sort {$a <=> $b} keys %{$probe_object_ref}){

    #check the format of each variable - final sanity check
    unless ($probe_object_ref->{$probe_id}->{target_hits} =~ /\d+/){
      print RED, "\nUndefined target_hits", RESET;
      exit();
    }
    unless ($probe_object_ref->{$probe_id}->{largest_target_al} =~ /\d+/){
      print RED, "\nUndefined largest_target_al", RESET;
      exit();
    }
    unless ($probe_object_ref->{$probe_id}->{non_target_hits} =~ /\d+/){
      print RED, "\nUndefined non-target_hits", RESET;
      exit();
    }
    unless ($probe_object_ref->{$probe_id}->{non_target_hits_over50} =~ /\d+/){
      print RED, "\nUndefined non-target_hits_over50", RESET;
      exit();
    }
    unless ($probe_object_ref->{$probe_id}->{non_target_hits_over75} =~ /\d+/){
      print RED, "\nUndefined non-target_hits_over75", RESET;
      exit();
    }
    unless ($probe_object_ref->{$probe_id}->{largest_non_target_al} =~ /\d+/){
      print RED, "\nUndefined largest_non_target_al", RESET;
      exit();
    }
    unless ($probe_object_ref->{$probe_id}->{largest_non_target_pi} =~ /\d+/){
      print RED, "\nUndefined largest_non_target_pi", RESET;
      exit();
    }

    print PROBE_OUT "$probe_id\t$probe_object_ref->{$probe_id}->{target_hits}\t$probe_object_ref->{$probe_id}->{largest_target_al}\t$probe_object_ref->{$probe_id}->{non_target_hits}\t$probe_object_ref->{$probe_id}->{non_target_hits_over50}\t$probe_object_ref->{$probe_id}->{non_target_hits_over75}\t$probe_object_ref->{$probe_id}->{largest_non_target_al}\t$probe_object_ref->{$probe_id}->{largest_non_target_pi}\n";
  }

  close PROBE_OUT;

  return();
}
