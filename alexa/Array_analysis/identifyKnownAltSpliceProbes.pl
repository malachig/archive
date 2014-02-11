#!/usr/bin/perl -w
#Written by Malachi Griffith

#This script considers a particular microarray design (ex. 'Druggable Genome' and the probes present on that array
#It requires that probe, microarray, and array spot data be imported into ALEXA

#For a particular microarray design get all the of the probes present on this design, identify all the genes targeted by these probes, and then:
#identify known alternative splicing events for these genes and the probes corresponding to them

#1.) Get probe info for the specified array design
#    - Info needed probe_type,unit1_start, unit1_end, unit2_start, unit2_end
#2.) Reorganize probe data so that it is grouped by gene ID (hash of hashes keyed on gene id)
#3.) Get the transcripts for every gene in this hash
#4.) Start going through the genes and identifying alternative splicing events for every gene with more than one transcript
#  A.) Identify exon skipping events.  Determine the probe corresponding to this skipping event
#      - Find the exon-exon probe with the same 'unit1_end' and 'unit2_start' as the exon skipping event
#  B.) Identify alternate exon boundary events.  Determine the probes that overlap the alternate boundary
#      - Find exon-intron probes that span one exon's boundary but are contained within another exon
#  C.) Identify exon probes that are specific to a single isoform??
#      - Or should I identify probes that are do not correspond to all known isoforms??  Potentially variant (non-constitutive) exons ...
#  D.) Identify alternate transcription start and end sites.  Determine the probes that are specific to the end of a transcript.
#      - Identify exons that are unique to a single isoform, then identify exon and canonical probes that are unique to these exons
#5.) Print a summary of known AS events to file.  A tab-delimited file containing info about AS event, probe that corresponds, etc.

use strict;
use Data::Dumper;
use Getopt::Long;
use DBI;

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
my $array_name = '';
my $out_file = '';

GetOptions ('database=s'=>\$database, 'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'array_name=s'=>\$array_name,
	    'out_file=s'=>\$out_file);

#Provide instruction to the user
print "\n\nThis script takes a gene file and microarray design name and identifies probes in this design that correspond to known AS events";
print "\nUsage:";
print "\n\tSpecify the target server using: --server";
print "\n\tSpecify the target database using: --database";
print "\n\tSpecify the user and password using: --user and --password";
print "\n\t\tMake sure you chose the correct database and server!!";
print "\n\tSpecify the array name using: --array_name";
print "\n\tSpecify an output file for the summary of known AS events";
print "\n\nExample: identifyKnownAltSpliceProbes.pl --server jango.bcgsc.ca --database ALEXA_hs_31_35d --user malachig --password pwd --array_name Druggable_Genome --out_file test.out\n\n";

#Make sure all options were specified
unless ($database && $server && $user && $password && $array_name && $out_file){
  print "\nOptions missing!\n\n";
  exit();
}

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = connectDB($database, $server, $user, $password);

#1.) Get probe info for the specified microarray design
#Info needed
print "\nRetrieving probes and associated data for Microarray design: $array_name\n\n";
my %microarray_probes = %{&getMicroarrayProbes ('-dbh'=>$alexa_dbh, '-name'=>"Druggable_Genome", '-id_type'=>"probe_id")};

#2.) Organize probe info by gene ID
my %genes;
&sortProbesByGene();

#3.) Get the transcripts for every gene in this hash
&getGeneTranscripts('-dbh'=>$alexa_dbh);

#Define a hash to store a summary of all known splicing events encountered
#Store: AS event count, AS event type (skip, alternate boundary, etc.), gene ID, probe ID corresponding to event,
#       Sequence ID where event was observed (Ensembl transcript, EST, mRNA, etc.), number of exons skipped if applicable
my $as_event_count;
my %as_events;

#4.) Start going through the genes and identifying alternative splicing events for every gene with more than one transcript
#  A.) Identify exon skipping events.  Determine the probe corresponding to this skipping event
#      - Find the exon-exon probe with the same 'unit1_end' and 'unit2_start' as the exon skipping event
#      - Identify exon probes that are specific to a single isoform
&identifyExonSkipProbes();

#  B.) Identify alternate exon boundary events.  Determine the probes that overlap the alternate boundary
#      - Find exon-intron probes that span one exon's boundary but are contained within another exon
&identifyAltBoundaryProbes();

#  C.) Identify exon probes that are specific to a single isoform??
#      - Or should I identify probes that are do not correspond to all known isoforms??  Potentially variant (non-constitutive) exons ...
&identifyNonConstitutiveExonProbes();


#5.) Print a summary of known AS events using the information gathered in %AS_events
&printAS_Events('-out_file'=>$out_file);


#Close database connection
$alexa_dbh->disconnect();

exit();



################################################################################################################
#2.) Organize probe info by Gene                                                                               #
################################################################################################################
sub sortProbesByGene{

  print "\nOrganizing probe info into gene units\n\n";
  foreach my $probe_id (keys %microarray_probes){
    my $gene_id = $microarray_probes{$probe_id}{gene_id};

    if ($genes{$gene_id}{probes}){
      my %gene_probes = %{$genes{$gene_id}{probes}};

      $gene_probes{$probe_id}{gene_id} = $microarray_probes{$probe_id}{gene_id};
      $gene_probes{$probe_id}{type} = $microarray_probes{$probe_id}{type};
      $gene_probes{$probe_id}{unit1_start} = $microarray_probes{$probe_id}{unit1_start};
      $gene_probes{$probe_id}{unit1_end} = $microarray_probes{$probe_id}{unit1_end};
      $gene_probes{$probe_id}{unit2_start} = $microarray_probes{$probe_id}{unit2_start};
      $gene_probes{$probe_id}{unit2_end} = $microarray_probes{$probe_id}{unit2_end};

      $genes{$gene_id}{probes} = \%gene_probes;
    }else{
      my %gene_probes;

      $gene_probes{$probe_id}{gene_id} = $microarray_probes{$probe_id}{gene_id};
      $gene_probes{$probe_id}{type} = $microarray_probes{$probe_id}{type};
      $gene_probes{$probe_id}{unit1_start} = $microarray_probes{$probe_id}{unit1_start};
      $gene_probes{$probe_id}{unit1_end} = $microarray_probes{$probe_id}{unit1_end};
      $gene_probes{$probe_id}{unit2_start} = $microarray_probes{$probe_id}{unit2_start};
      $gene_probes{$probe_id}{unit2_end} = $microarray_probes{$probe_id}{unit2_end};

      $genes{$gene_id}{probes} = \%gene_probes;
    }
  }

  my $gene_count = keys %genes;

  print "\nFound $gene_count genes associated with the probes of the microarray: $array_name\n\n";

  return();
}


################################################################################################################
#3.) Get the Transcripts for each Gene                                                                         #
################################################################################################################
sub getGeneTranscripts{
  my %args = @_;
  my $dbh = $args{'-dbh'};

  foreach my $gene_id (keys %genes){
    my %transcripts = %{&getTranscripts ('-dbh'=>$dbh, '-gene_id'=>$gene_id)};

    $genes{$gene_id}{transcripts} = \%transcripts;
    $genes{$gene_id}{trans_count} = keys %transcripts;
  }
  return();
}


################################################################################################################
#4-A.) Identify exon skip probes                                                                               #
################################################################################################################
sub identifyExonSkipProbes{

  print "\nIdentifying probes representing exon skipping events\n\n";

  my $multi_transcript_genes = 0;
  my $exon_skipping_events = 0;
  my @exon_skip_counts;

  foreach my $gene_id (sort {$a <=> $b} keys %genes){

    my %transcripts = %{$genes{$gene_id}{transcripts}};

    #Unless there are multiple transcripts, skip this gene
    unless ($genes{$gene_id}{trans_count} > 1){
      next();
    }
    $multi_transcript_genes++;

    print "\n\nGENE: $gene_id";

    #1.) Assemble a reference set of exons (a superset of all non-redundant exons)
    my %reference_exons;
    foreach my $trans_id (keys %transcripts){
      my %exons = %{$transcripts{$trans_id}{exons}};

      foreach my $trans_exon_id (keys %exons){
	my $trans_exon_start = $exons{$trans_exon_id}{exon_start};
	my $trans_exon_end = $exons{$trans_exon_id}{exon_end};

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
    print "\nREF EXON STARTS: @reference_exon_starts";
    print "\nREF EXON ENDS:   @reference_exon_ends";

    #3.) Now go through each exon of each transcript, compare to the reference exons and identiy exon skipping events
    foreach my $trans_id (keys %transcripts){
      my %exons = %{$transcripts{$trans_id}{exons}};

      my $ensembl_t_id = $transcripts{$trans_id}{ensembl_t_id};

      #Sort the exons for this transcript and create arrays of their start and end positions
      my @trans_exon_starts;
      my @trans_exon_ends;
      foreach my $trans_exon_id (sort {$exons{$a}->{exon_start} <=> $exons{$b}->{exon_start}} keys %exons){
	push (@trans_exon_starts, $exons{$trans_exon_id}{exon_start});
	push (@trans_exon_ends, $exons{$trans_exon_id}{exon_end});
      }
      #Determine the number of exons in this transcript
      my $trans_exon_count = keys %exons;

      print "\n\tTRANS: $trans_id ($trans_exon_count exons)";
      print "\n\t\tTRANS EXON STARTS: @trans_exon_starts";
      print "\n\t\tTRANS EXON ENDS:   @trans_exon_ends";

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
	      print "\n\t\t\tObserved $skipped_exons skipped exons between ($trans_exon1_start - $trans_exon1_end) and ($trans_exon2_start - $trans_exon2_end)";

	      #Now attempt to identify an exon skipping probe that corresponds to this event
	      my %gene_probes = %{$genes{$gene_id}{probes}};
	      foreach my $probe_id (keys %gene_probes){
		#Skip probes unless they are 'Exon-Exon' probes
		unless ($gene_probes{$probe_id}{type} eq "Exon-Exon"){
		  next();
		}

		#See if this probe's 'unit1_end' and unit2_start' matches that of the observed exon skipping event
		if ($gene_probes{$probe_id}{unit1_end} == $trans_exon1_end && $gene_probes{$probe_id}{unit2_start} == $trans_exon2_start){
		  print "\n\t\t\tPROBE: $probe_id matches this exon skipping event";
		  $as_event_count++;
		  $as_events{$probe_id}{type} = "Exon-Exon";
		  $as_events{$probe_id}{gene_id} = $gene_id;
		  if ($as_events{$probe_id}{sequence_ids}){
		    my @tmp = @{$as_events{$probe_id}{sequence_ids}};
		    push (@tmp, $ensembl_t_id);
		    $as_events{$probe_id}{sequence_ids} = \@tmp;
		  }else{
		    my @tmp;
		    push (@tmp, $ensembl_t_id);
		    $as_events{$probe_id}{sequence_ids} = \@tmp;
		  }
		}		
	      }
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
  return();
}

###########################################################################################################################
#4-B.) Identify alternate exon boundary events.  Determine the probes that overlap the alternate boundary                 #
#      - Find exon-intron probes that span one exon's boundary but are contained within another exon                      #
###########################################################################################################################
sub identifyAltBoundaryProbes{

  print "\nIdentifying probes representing alternate exon boundaries\n\n";

  foreach my $gene_id (sort {$a <=> $b} keys %genes){

    my %transcripts = %{$genes{$gene_id}{transcripts}};

    #Unless there are multiple transcripts, skip this gene
    unless ($genes{$gene_id}{trans_count} > 1){
      next();
    }

    print "\n\nGENE: $gene_id";

    #Go through each exon of each transcript, check each exon-intron probe to see if it is completely contained within one of these exons

    foreach my $trans_id (keys %transcripts){
      my %exons = %{$transcripts{$trans_id}{exons}};

      my $ensembl_t_id = $transcripts{$trans_id}{ensembl_t_id};

      #Go through each exon for this transcript
      foreach my $trans_exon_id (sort {$exons{$a}->{exon_start} <=> $exons{$b}->{exon_start}} keys %exons){

	my $trans_exon_start = $exons{$trans_exon_id}{exon_start};
	my $trans_exon_end = $exons{$trans_exon_id}{exon_end};

	#Get the Exon-Intron and Intron-Exon probes
	my %gene_probes = %{$genes{$gene_id}{probes}};

	foreach my $probe_id (keys %gene_probes){
	  #Skip probes unless they are 'Exon-Intron' or 'Intron-Exon' probes
	  unless ($gene_probes{$probe_id}{type} eq "Exon-Intron" || $gene_probes{$probe_id}{type} eq "Intron-Exon"){
	    next();
	  }

	  #See if this probe's 'unit1_start' and unit2_end' are contained within this exon
	  if ($gene_probes{$probe_id}{unit1_start} >= $trans_exon_start && $gene_probes{$probe_id}{unit2_end} <= $trans_exon_end){
	    print "\n\t\t\tPROBE: Exon boundary probe: $probe_id is contained within an exon ($trans_exon_start - $trans_exon_end) of $ensembl_t_id";

	    $as_events{$probe_id}{type} = $gene_probes{$probe_id}{type};
	    $as_events{$probe_id}{gene_id} = $gene_id;
	    if ($as_events{$probe_id}{sequence_ids}){
	      my @tmp = @{$as_events{$probe_id}{sequence_ids}};
	      push (@tmp, $ensembl_t_id);
	      $as_events{$probe_id}{sequence_ids} = \@tmp;
	    }else{
	      my @tmp;
	      push (@tmp, $ensembl_t_id);
	      $as_events{$probe_id}{sequence_ids} = \@tmp;
	    }
	  }		
	}
      }
    }
  }

  return();
}


##############################################################################################################################################
#  C.) Identify exon probes that are specific to a single isoform??                                                                          #
#      - Or should I identify probes that are do not correspond to all known isoforms??  Potentially variant (non-constitutive) exons ...    #
##############################################################################################################################################
sub identifyNonConstitutiveExonProbes{

  print "\nIdentifying exon probes that are not present in every transcript of a gene (non-constitutive probes)\n\n";

  foreach my $gene_id (sort {$a <=> $b} keys %genes){

    my %transcripts = %{$genes{$gene_id}{transcripts}};

    #Unless there are multiple transcripts, skip this gene
    unless ($genes{$gene_id}{trans_count} > 1){
      next();
    }

    print "\n\nGENE: $gene_id";

    my $trans_count = $genes{$gene_id}{trans_count};

    #Go through each probe for this gene, attempt to find it in each transcript, if it is not found in every transcript, make note of it
    my %gene_probes = %{$genes{$gene_id}{probes}};

    foreach my $probe_id (keys %gene_probes){
      #Skip probes unless they are 'Exon' probes
      unless ($gene_probes{$probe_id}{type} eq "Exon"){
	next();
      }
      my $probe_trans_count = 0;  #Number of transcripts containing this probe sequence
      my $probe_start = $gene_probes{$probe_id}{unit1_start};
      my $probe_end = $gene_probes{$probe_id}{unit1_end};

      #Go through all of the exons of each transcript and count the number of times the probe sequence is within an exon
      my @transcripts;
    TRANS: foreach my $trans_id (keys %transcripts){
	  my %exons = %{$transcripts{$trans_id}{exons}};

	  my $ensembl_t_id = $transcripts{$trans_id}{ensembl_t_id};

	  foreach my $trans_exon_id (sort {$exons{$a}->{exon_start} <=> $exons{$b}->{exon_start}} keys %exons){

	    my $trans_exon_start = $exons{$trans_exon_id}{exon_start};
	    my $trans_exon_end = $exons{$trans_exon_id}{exon_end};

	    #Test if the start AND end position of the probe are contained within this exon 
	    if ($probe_start >= $trans_exon_start && $probe_start <= $trans_exon_end && $probe_end >= $trans_exon_start && $probe_end <= $trans_exon_end){
	      $probe_trans_count++;
	      push (@transcripts, $ensembl_t_id);
	      next TRANS;
	    }
	  }
	}

      #If this exon probe was not found in every transcript, make note of it as a non-constitutive probe
      if ($probe_trans_count < $trans_count){
	print "\n\t\t\tPROBE: Exon probe: $probe_id is non-constitutive ($probe_start - $probe_end) but occurs in @transcripts";

	$as_events{$probe_id}{type} = $gene_probes{$probe_id}{type};
	$as_events{$probe_id}{gene_id} = $gene_id;
	$as_events{$probe_id}{sequence_ids} = \@transcripts;
      }

      #Make sure this probe was found in at least one transcript
      if ($probe_trans_count == 0){
	print "\n*** Could not find probe: $probe_id in any transcript of gene: $gene_id ***\n\n";
      }
    }
  }

  return();
}


######################################################################################
#5.) Print a summary of known AS events using the information gathered in %as_events #
######################################################################################
sub printAS_Events{
  my %args = @_;
  my $out_file = $args{'-out_file'};

  open (OUT, ">$out_file") || die "\nCould not open output file: $out_file\n\n";

  print OUT "Probe ID\tGene ID\tAS Type\tSequence IDs\n";

  foreach my $probe_id (sort {$as_events{$a}->{gene_id} <=> $as_events{$b}->{gene_id}} keys %as_events){

    my @seq_ids = @{$as_events{$probe_id}{sequence_ids}};

    print OUT "$probe_id\t$as_events{$probe_id}{gene_id}\t$as_events{$probe_id}{type}\t@seq_ids\n";

  }

  close OUT;

  return();
}





