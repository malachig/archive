#!/usr/bin/perl -w

#Written by Malachi Griffith
#The purpose of this script is to test the specificity of probes in an input file of probe records.  It relies on the existence of BLAST results 
#prepared in advance using megaBLAST (likely on the cluster) and stored in a results directory.  Before running this script, first use:
#   blastProbesVersus_EST-mRNA.pl to generate probe fasta files and then run the blast jobs that it produces on the cluster.
#This script will thus rely on a probe file from /home/user/alexa/ALEXA_version/specificity/job_name/probes and will look for the blast results for
#these probes in /home/user/alexa/ALEXA_version/specificity/job_name/blast_results

#Note some of the following may be very memory intensive!
#Script process:
#1.) Load the probes of interest which have already been blasted against mRNA and EST databases.  Once the test of specificity has been completed
#    the results will be appended to this file and written to a new output file. Summarize the number and types of probes loaded.
#    Create a 'Probes' hash (keyed on probe_ID) containing all info from the input file
#    Confirm that no probe ID is encountered more than once
#2.) Load all BLAST results from the blast-result file specified.  Summarize the number of blast results contained in this file
#    Create a 'Blast_results' hash (keyed on probe ID) which contains references to arrays or hashes of blast hits
#    Only store the minimum amount of info in the hash to save memory
#3.) Load genomic alignment information from UCSC tables stored in ALEXA.  Note: because it would so cumputationally intensive to map all human
#    ESTs and mRNAs to the genome myself, I am relying on UCSC to do this for me.  Tables containing this info were downloaded from the UCSC ftp site
#    and dumped as stand-alone tables into ALEXA.
#    Create a 'Alignments' hash (keyed on sequence Accession ID) which contains references to arrays or hashed of alignment positions to the genome
#4.) Load the genomic positions for each Ensembl gene in ALEXA.
#    Create a 'Genes' hash (keyed on alexa_gene_id) containing the genomic position of each gene.
#5.) Foreach 'Probe' get the 'Blast_results' found and find the corresponding genomic positions of the mRNA/ESTs in the 'Alignments'.  If a blast hit
#    corresponds to a sequence which aligns to the genome in the region of the target 'Gene' then this is okay.  However, if a hit corresponds to a 
#    sequence that aligns outside of the genomic region targeted, even if it also aligns there, the probe can not be considered to be specific to the 
#    targeted gene region.  This will be noted as a failure.
#    - The result will be recorded as the number of significant hits to the sequences that align outside the target region.
#Note1: This script does not consider the strand that each probe is hitting.  The tabular output format of BLAST does not include this info.  This means
#       that a hit on either strand to a sequence that maps outside of the targeted genomic region will fail a probe sequence.  Perhaps this is too conservative.
#       -Actually, if you are hybridizing double-stranded cDNA samples to the array, this is a really good idea!
#Note2: After initial experimentation I decided that it is neccessary to consider the position of the probe bit to the mRNA.  If a probe hits an mRNA, and that
#       mRNA aligns to a region of the genome other than the target, this is only a problem if the region that aligns to a non-target part of the genome is the
#       same region that the probe hit.  When I first wrote this a probe could be failed if it hit an mRNA at the 3' end and that mRNA has a small sequence at
#       the 5' end that aligns to a non-target region of the genome.  There may be many cases like this where an mRNA really does come from the region targetted
#       but part of its sequence is non-specific/low complexity and aligns all over the genome.  This is ONLY a problem for probes that come from that part of the sequence.
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
my $db_type = '';
my $probe_infile = '';
my $blast_results_infile = '';
my $outfile = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'db_type=s'=>\$db_type,
	    'probe_infile=s'=>\$probe_infile, 'blast_results_infile=s'=>\$blast_results_infile, 'outfile=s'=>\$outfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script test the specificity of probes against all mRNA or ESTs in the specified ALEXA database", RESET;
print GREEN, "\n\tThis script takes a probe file and blast results file as input, parses the results and creates and updated output file", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password\n", RESET;
print GREEN, "\n\tSpecify the type of blast database being used using: --db_type (EST or mRNA)", RESET;
print GREEN, "\n\tSpecify a tab-delimited input probe file using: --probe_infile", RESET;
print GREEN, "\n\tSpecify a tab-delimited blast results file using: --blast_results_file", RESET;
print GREEN, "\n\tSpecify the name of the resulting outfile using: --outfile", RESET;
print GREEN, "\n\nExample: testProbeSpecificity_EST-mRNA.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --db_type=mRNA  --probe_infile=exonJunctionProbes.txt  --blast_results_infile=blast_results.txt  --outfile=exonJunctionProbes_Spec_mRNA.txt\n\n", RESET;

unless ($database && $server && $user && $password && $db_type && $probe_infile && $blast_results_infile && $outfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Specify the blast result type to use - mRNA or EST
unless ($db_type eq "mRNA" || $db_type eq "EST"){
  print RED, "\nMust specify the type of specificity test to conduct with --db_type=mRNA or --db_type=EST\n\n", RESET;
  exit();
}

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

#1.) Load the probes of interest which have already been blasted against mRNA and EST databases
my $header_line;
my %master_gene_list;
my $probes_ref = &loadProbeInfo('-input_file'=>$probe_infile);

#2.) Load all BLAST results from the specified blast-results file
my $blast_results_ref = &loadBlastResults('-blast_results_file'=>$blast_results_infile, '-blast_type'=>$db_type);

#3.) Load genomic alignment information from UCSC tables stored in ALEXA.
my $alignments_ref = &loadAlignmentInfo('-dbh'=>$alexa_dbh, '-type'=>$db_type);

#4.) Load the genomic positions for each Ensembl gene in ALEXA
my $genes_ref = &loadGeneInfo('-dbh'=>$alexa_dbh);

#5.) Using the infomation gathered in steps 1-4, conduct the actual specificity test
&testProbeSpecificity('-probe_object'=>$probes_ref, '-blast_object'=>$blast_results_ref, '-alignment_object'=>$alignments_ref, '-gene_object'=>$genes_ref);

#6.) Print the probe info out to a new file and append the results of the specificity test
&printResults('-output_file'=>$outfile, '-probe_object'=>$probes_ref, '-header_line'=>$header_line, '-type'=>$db_type);

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

    $master_gene_list{$gene_id}{tmp} = '';

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
  my $type = $args{'-blast_type'};

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

  print BLUE, "\nProcessed $blast_results_processed blast results", RESET;
  print BLUE, "\n\tFound $less_significant_hits that were of the same or shorter length to another hit on the same sequence", RESET;
  print BLUE, "\n\tActually stored $blast_results_stored blast results", RESET;

  return(\%blast_results_object);
}


#########################################################################################################
#3.) Load genomic alignment information from UCSC tables stored in ALEXA.                               #
#########################################################################################################
sub loadAlignmentInfo{
  my %args = @_;
  my $dbh = $args{'-dbh'};
  my $type = $args{'-type'};

  my %alignment_object;

  #Specify the correct table to retrieve alignment data from depending on whether the blast results are from ESTs or mRNAs
  my $table_name;
  if ($type eq "mRNA"){
    $table_name = "all_mrna";
  }elsif($type eq "EST"){
    $table_name = "all_est";
  }else{
    print RED, "\nBlast analysis type not understood\n\n", RESET;
    exit();
  }

  my $sql = "SELECT strand,qName,qStart,qEnd,tName,tStart,tEnd,matches FROM $table_name";

  my $sth = $dbh->prepare("$sql");
  $sth->execute();

  while (my ($strand,$acc,$acc_start,$acc_end,$chr,$chr_start,$chr_end, $matches) = $sth->fetchrow_array()){

    #Eliminate small alignments!  ie. if there are less than say 100 bp of matches ignore the alignment
    #It is probably still a good idea to eliminate alignments to short sequences because they may not align to the genome very well
    unless ($matches > 100){
      next();
    }

    #If there are already alignments for this accession - get them and add to the growing list
    if ($alignment_object{$acc}){
      my $align_count = $alignment_object{$acc}{align_count};
      $align_count++;

      my $align_ref = $alignment_object{$acc}{alignments};
      $align_ref->{$align_count}->{strand} = $strand;
      $align_ref->{$align_count}->{acc_start} = $acc_start;
      $align_ref->{$align_count}->{acc_end} = $acc_end;
      $align_ref->{$align_count}->{chr} = $chr;
      $align_ref->{$align_count}->{chr_start} = $chr_start;
      $align_ref->{$align_count}->{chr_end} = $chr_end;
      $alignment_object{$acc}{align_count} = $align_count;

    }else{
      #If there are no alignments for this accession yet, create a new hash
      my $align_count = 1;
      $alignment_object{$acc}{align_count} = $align_count;
      my %align;
      $align{$align_count}{strand} = $strand;
      $align{$align_count}{acc_start} = $acc_start;
      $align{$align_count}{acc_end} = $acc_end;
      $align{$align_count}{chr} = $chr;
      $align{$align_count}{chr_start} = $chr_start;
      $align{$align_count}{chr_end} = $chr_end;

      $alignment_object{$acc}{alignments} = \%align;
    }

  }
  $sth->finish();

  my $acc_found = keys %alignment_object;
  print BLUE, "\n\nImported UCSC Genomic alignments for $acc_found $type sequences", RESET;

  return(\%alignment_object);
}


#########################################################################################################
#4.) Load the genomic positions for each Ensembl gene in ALEXA                                          #
#########################################################################################################
sub loadGeneInfo{
  my %args = @_;

  my $dbh = $args{'-dbh'};

  my %gene_object;

  #Get gene info for only the ALEXA genes targeted by the probes in the input probe file
  my @gene_ids = keys %master_gene_list;
  my $gene_info_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids,'-sequence'=>"no");

  foreach my $gene_id (@gene_ids){
    $gene_object{$gene_id}{chr} = $gene_info_ref->{$gene_id}->{chromosome};
    $gene_object{$gene_id}{chr_start} = $gene_info_ref->{$gene_id}->{chr_start};
    $gene_object{$gene_id}{chr_end} = $gene_info_ref->{$gene_id}->{chr_end};

    #correct strand format so that it can be compared to UCSC strands
    if ($gene_info_ref->{$gene_id}->{chr_strand} eq "1"){
      $gene_object{$gene_id}{strand} = "+";
    }elsif ($gene_info_ref->{$gene_id}->{chr_strand} eq "-1"){
      $gene_object{$gene_id}{strand} = "-";
    }else{
      print RED, "\nStrand format not understood\n\n", RESET;
      exit();
    }
  }

  my $genes_found = keys %gene_object;
  print BLUE, "\n\nImported $genes_found Ensembl Genes with coordinates from ALEXA", RESET;

  return(\%gene_object);
}


#########################################################################################################
#5.) Using the infomation gathered in steps 1-4, conduct the actual specificity test                    #
#########################################################################################################
sub testProbeSpecificity{
  my %args = @_;
  my $probe_object_ref = $args{'-probe_object'};
  my $blast_object_ref = $args{'-blast_object'};
  my $align_object_ref = $args{'-alignment_object'}; #UCSC alignments of mRNAs/ESTs to the genome
  my $gene_object_ref = $args{'-gene_object'};

  #Use a hash to keep track of sequence accessions (mRNA/EST) from BLAST hits that can not be found in the UCSC alignments table
  #Compare to the number of sequence accessions from BLAST hits that ARE found in the UCSC alignments table
  my %unmapped_seqs;
  my %mapped_seqs;
  my %target_seqs;
  my %antisense_seqs;

  print BLUE, "\n\nBegin actual specificity test", RESET;

  #Go through each probe ID and test its specificity to the target locus
  foreach my $probe_id (sort {$a <=> $b} keys %{$probe_object_ref}){

    #Count the number of times a probe hits a sequence which maps outside the targeted region
    my $non_target_hits = 0;
    my $target_hits = 0;

    #Initialize variables
    $probe_object_ref->{$probe_id}->{non_target_hits} = $non_target_hits;
    $probe_object_ref->{$probe_id}->{target_hits} = $target_hits;

    #Keep track of the longest non-target and target alignment lengths for each probe
    $probe_object_ref->{$probe_id}->{largest_non_target_al} = "0";
    $probe_object_ref->{$probe_id}->{largest_target_al} = "0";

    #Get the gene ID for this probe
    my $alexa_gene_id = $probe_object_ref->{$probe_id}->{alexa_gene_id};

    #Make sure gene info exists for this gene ID
    unless ($gene_object_ref->{$alexa_gene_id}){
      print RED, "\nFound a probe associated with a gene for which gene coordinates were not found\n\n", RESET;
      exit();
    }

    #Get the chromosome, chromosome start, chromosome end and strand for this gene
    my $gene_chr = $gene_object_ref->{$alexa_gene_id}->{chr};
    my $gene_chr_start = $gene_object_ref->{$alexa_gene_id}->{chr_start};
    my $gene_chr_end = $gene_object_ref->{$alexa_gene_id}->{chr_end};
    my $gene_strand = $gene_object_ref->{$alexa_gene_id}->{strand};

    #Get the blast hits for this probe (as a hash keyed on accession ID)
    #Check for probes with no blast hits
    my $blast_hits_ref;
    if ($blast_object_ref->{$probe_id}->{blast_hits}){
      $blast_hits_ref = $blast_object_ref->{$probe_id}->{blast_hits};
    }else{
      next();
    }

    #Go through each BLAST hit and get the UCSC genomic alignments for each of these sequences that the probe has a blast hit for
    BLAST_HIT:foreach my $blast_acc (sort keys %{$blast_hits_ref}){

	#Make sure the sequence accession in the blast result can be found in the UCSC alignments table
	#If it can not, perhaps this is a sequence which could not be aligned to the genome at any position - note how often this happens
	unless($align_object_ref->{$blast_acc}){
	  $unmapped_seqs{$blast_acc}{tmp} = '';
	  next BLAST_HIT;
	}
	$mapped_seqs{$blast_acc}{tmp} = '';

	my $seq_alignments_ref = $align_object_ref->{$blast_acc}->{alignments};

	#Does this sequence to which the probe has a BLAST hit align to the targetted genomic region only?
	#Go through each alignment of this sequence to the genome ...
	GENOME_MAP: foreach my $seq_align (sort keys %{$seq_alignments_ref}){

	    #Make sure the alignment of the probe to the mRNA/EST is in the same region as the alignment of the mRNA/EST to the genome
	    #Basically, make sure we are talking about the same region of the mRNA/EST with regards to the probe and genome alignments

	    #Alignment start/end positions within the mRNA/EST sequence to genome
	    my $seq_map_start = ($seq_alignments_ref->{$seq_align}->{acc_start})+1; #Start position within the mRNA/EST for its alignment to the genome.
	    my $seq_map_end = ($seq_alignments_ref->{$seq_align}->{acc_end})+1;     #End position within the mRNA/EST for its alignment to the genome.

	    #Alignment start/end positions within the mRNA/EST sequence to the probe
	    my $probe_align_start = $blast_hits_ref->{$blast_acc}->{ss}; #Start position within the mRNA/EST of the alignment to the probe
	    my $probe_align_end = $blast_hits_ref->{$blast_acc}->{se};   #End position within the mRNA/EST of the alignment to the probe

	    #Now before processing any farther with this sequence alignment to the genome, confirm that it involves the same region
	    #of the mRNA/EST as the probe alignment to this sequence

	    #If the alignment of the probe to the mRNA/EST overlaps at either end of alignment of the mRNA/EST to the genome this okay
	    if (($probe_align_start >= $seq_map_start && $probe_align_start <= $seq_map_end) || ($probe_align_end >= $seq_map_start && $probe_align_end <= $seq_map_end)){
	      #Do nothing

	      #Note that the above condition will miss cases where the probe alignment outflanks the genomic alignment region at boths ends
	      #This should not happen since I have not allowed small alignments to the genome anyway, unless the probe is very long
	    }elsif($probe_align_start <= $seq_map_start && $probe_align_end >= $seq_map_end){
	      #Do nothing
	    }else{
	      #The alignment to the genome involves a different region of the mRNA/EST sequence than the region hit by the probe.
	      #Ignore this genomic alignment
	      next GENOME_MAP;
	    }

	    #Get the chromosome, chromosome start, and chromosome end positions for this sequence (mRNA/EST) alignment to the genome
	    my $seq_chr_unformated = $seq_alignments_ref->{$seq_align}->{chr}; #Note: Chromosome names are stored in UCSC tables as 'chr1, chrX, etc.)
	    my $seq_chr;
	    if ($seq_chr_unformated =~ /^chr(.*)$/){
	      $seq_chr = $1;
	    }else{
	      print RED, "\nChromosome name format from UCSC tables not understood\n\n", RESET;
	      exit();
	    }
	    my $seq_chr_start = $seq_alignments_ref->{$seq_align}->{chr_start};
	    my $seq_chr_end = $seq_alignments_ref->{$seq_align}->{chr_end};

	    #1.) First make sure the hit is to the same chromosome as the targeted gene
	    unless ($seq_chr eq $gene_chr){
	      $non_target_hits++;

	      #update the largest non-target hit found so far if neccessary
	      if ($blast_hits_ref->{$blast_acc}->{al} > $probe_object_ref->{$probe_id}->{largest_non_target_al}){
		$probe_object_ref->{$probe_id}->{largest_non_target_al} = $blast_hits_ref->{$blast_acc}->{al};
	      }

	      #print "\n\nWrong chromosome!";
	      #print "\nProbe= $probe_id\tNonTargetHit\tGene: $gene_chr ($gene_chr_start - $gene_chr_end)\tSeq: $blast_acc: $seq_chr ($seq_chr_start - $seq_chr_end)\tBlast_hit_length= $blast_hits{$blast_acc}{al}";
	      next BLAST_HIT;
	    }

	    #2.) Then make sure the hit is within the same genomic coordinates as the targeted gene
	    #If the alignment overlaps at either end of the gene sequence or is completely within it this will be considered okay
	    #This is somewhat sloppy because technically the actual position of the hit should be considered but this would be more complicated
	    if (($seq_chr_start >= $gene_chr_start && $seq_chr_start <= $gene_chr_end) || ($seq_chr_end >= $gene_chr_start && $seq_chr_end <= $gene_chr_end)){
	      #Do nothing 
	      $target_seqs{$blast_acc}{tmp} = '';
	      #Note that the above condition will miss cases where the sequence alignment outflanks the genomic target region at boths ends
	      #Since this also indicates that the hit is to the correct targeted region I need a second condition to allow this case

	    }elsif($seq_chr_start <= $gene_chr_start && $seq_chr_end >= $gene_chr_end){
	      #Do nothing
	      $target_seqs{$blast_acc}{tmp} = '';
	    }else{
	      $non_target_hits++;

	      #update the largest non-target hit found so far if neccessary
	      if ($blast_hits_ref->{$blast_acc}->{al} > $probe_object_ref->{$probe_id}->{largest_non_target_al}){
		$probe_object_ref->{$probe_id}->{largest_non_target_al} = $blast_hits_ref->{$blast_acc}->{al};
	      }

	      #print "\n\nNon-target hit!";
	      #print "\nProbe= $probe_id\tNonTargetHit\tGene: $gene_chr ($gene_chr_start - $gene_chr_end)\tSeq: $blast_acc: $seq_chr ($seq_chr_start - $seq_chr_end)\tBlast_hit_length= $blast_hits{$blast_acc}{al}";
	      #print "\nProbeAlignTo mRNA/EST: ($probe_align_start - $probe_align_end)\tmRNA/EST Region Involved in Align to Genome: ($seq_map_start - $seq_map_end)";
	      next BLAST_HIT;
	    }

	    #3.) To make it this far, the hit has to have been to the same chromosome and within the target coordinates
	    #    So this is most likely a valid hit to the target gene, however, what if the hit is to a transcript that maps to the opposite strand??
	    #    We want to eliminate such cases because if you hybridize ds-cDNA, you will not be able to tell where the observed expression is coming from
	    #    Make note of how often this happens
	    my $seq_strand = $seq_alignments_ref->{$seq_align}->{strand};

	    unless ($seq_strand eq $gene_strand){
	      $non_target_hits++;

	      #update the largest non-target hit found so far if neccessary
	      if ($blast_hits_ref->{$blast_acc}->{al} > $probe_object_ref->{$probe_id}->{largest_non_target_al}){
		$probe_object_ref->{$probe_id}->{largest_non_target_al} = $blast_hits_ref->{$blast_acc}->{al};
	      }
	      $antisense_seqs{$blast_acc}{tmp} = '';
	      next BLAST_HIT;
	    }

	    #4.) If the hit passed all the above filters, it will be considered a valid target hit
	    #update the largest target hit found so far if neccessary
	    if ($blast_hits_ref->{$blast_acc}->{al} > $probe_object_ref->{$probe_id}->{largest_target_al}){
	      $probe_object_ref->{$probe_id}->{largest_target_al} = $blast_hits_ref->{$blast_acc}->{al};
	    }

	    $target_hits++;
	  }
      }
    $probe_object_ref->{$probe_id}->{non_target_hits} = $non_target_hits;
    $probe_object_ref->{$probe_id}->{target_hits} = $target_hits;
  }

  my $unmapped_seqs = keys %unmapped_seqs;
  my $mapped_seqs = keys %mapped_seqs;
  my $target_seqs = keys %target_seqs;
  my $antisense_seqs = keys %antisense_seqs;

  print BLUE, "\n\tFound $unmapped_seqs mRNA or EST sequences for which a probe matches but could not be found in the UCSC genomic alignment tables", RESET;
  print BLUE, "\n\tFound $mapped_seqs mRNA or EST sequences which could be successfully mapped to the genome in the UCSC genomic alignment tables", RESET;
  print BLUE, "\n\tFound $target_seqs mRNA or EST sequences (with probe matches) which map to the genomic region targeted by the probe", RESET;
  print BLUE, "\n\tOf these, $antisense_seqs mRNA or EST sequences are from the target genomic region of a probe BUT on the opposite strand", RESET;

  return();
}


#########################################################################################################
#6.) Print the probe info out to a new file and append the results of the specificity test              #
#########################################################################################################
sub printResults{
  my %args = @_;

  my $output_file = $args{'-output_file'};
  my $probe_object_ref = $args{'-probe_object'};
  my $header_line = $args{'-header_line'};
  my $type = $args{'-type'};

  print BLUE, "\n\nBegin printing results to $output_file\n", RESET;

  my $new_column1 = "$type"."_Non-TargetHits";
  my $new_column2 = "$type"."_TargetHits";
  my $new_column3 = "$type"."_largestNon-TargetHitLength";
  my $new_column4 = "$type"."_largestTargetHitLength";

  open (PROBE_OUT, ">$output_file") || die "\nCould not open output file: $output_file\n\n";

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
