#!/usr/bin/perl -w

#Written by Malachi Griffith
#This script blasts probes against the complete gene sequence of their target gene (i.e. the gene they were extracted from)
#- It attempts to identify significant hits for regions other than where the probe was extracted from.
#- The point of this is to identify probes that are derived from regions that are repeated WITHIN the target gene

#1.) Import probes from an input probe file
#    - Group probes according to their target gene
#    - Store probe coordinate mapping info for use in determining target vs. non-target hits

#Next, go through each gene and:
#2.) Create a query fasta file of all probe sequences in a temp directory
#3.) Get the full genomic sequence for this target gene and create a fasta file for in a temp directory
#4.) Use formatDB to create a blastable database for this single genomic sequence
#5.) Blast all probes in the query sequence against this temp database and create a blast results file in the temp directory
#6.) Parse the blast results file
#7.) Test the position of hits to the target gene and compare to those expected
#    - Identify hits that do not overlap with the coordinates where the probe was actually extracted from as non-target hits
#    - Hits that do overlap are target hits
#8.) Print the probe specificity test results to the master out file

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
my $blast_bin_dir = '';
my $word_size = '';
my $multiple_cpus = '';
my $probe_file = '';
my $temp_dir = '';
my $job_name = '';
my $outfile = '';
my $logfile = '';
my $force = '';
my $big_mem = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'word_size=i'=>\$word_size, 'multiple_cpus=i'=>\$multiple_cpus, 'blast_bin_dir=s'=>\$blast_bin_dir,
	    'probe_file=s'=>\$probe_file, 'temp_dir=s'=>\$temp_dir, 'job_name=s'=>\$job_name,
	    'outfile=s'=>\$outfile, 'logfile=s'=>\$logfile, 'force=s'=>\$force, 'big_mem=s'=>\$big_mem);

#Provide instruction to the user
print GREEN, "\n\nThis script blasts probes against their target gene and parses the results to identify non-specific probes", RESET;
print GREEN, "\nThe goal is to identify probes corresponding to repeated regions WITHIN the target gene", RESET;
print GREEN, "\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password\n", RESET;
print GREEN, "\n\tSpecify the full path to a directory containing blast binaries using: --blast_bin_dir", RESET;
print GREEN, "\n\tSpecify the blast wordsize using: --word_size", RESET;
print GREEN, "\n\tTo use multiple cpus for blast commands use:  --multiple_cpus (e.g. --multiple_cpus=2)", RESET;
print GREEN, "\n\tSpecify the input probe file using: --probe_file", RESET;
print GREEN, "\n\tSpecify the temp directory to use for processing each gene target using: --temp_dir", RESET;
print GREEN, "\n\tSpecify a descriptive job name using: --job_name", RESET;
print GREEN, "\n\tSpecify the name of the resulting outfile using: --outfile", RESET;
print GREEN, "\n\tSpecify a logfile for output using: --logfile", RESET;
print GREEN, "\n\tTo prevent user prompts, use: --force=yes", RESET;
print GREEN, "\n\tIf you have access to a large memory cpu, use --big_mem=yes to get all gene sequences at once", RESET;
print GREEN, "\n\nExample: blastProbesVersus_TargetGene.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --blast_bin_dir=/home/user/BioSw/BLAST2/blast2.2.15/bin  --word_size=11  --multiple_cpus=1  --probe_file=/home/user/alexa/ALEXA_version/unfiltered_probes/exonJunctionProbes.txt  --temp_dir=/tmp/  --job_name=within_gene_hs_ej  --outfile=/home/user/alexa/ALEXA_version/unfiltered_probes/exonJunctionProbes_Spec_withinGene.txt  --logfile=/home/user/alexa/ALEXA_version/logs/specificity/blastProbesVersus_TargetGene_exonJunction_LOG.txt\n\n", RESET;

unless ($database && $server && $user && $password && $word_size && $blast_bin_dir && $probe_file && $temp_dir && $job_name && $outfile && $logfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Open logfile for output
open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";
print LOG "\nUser Specified the following options:\ndatabase = $database\nblast_bin_dir = $blast_bin_dir\nword_size = $word_size\nprobe_file = $probe_file\ntemp_dir = $temp_dir\njob_name = $job_name\noutfile = $outfile\nlogfile = $logfile\n\n";

#Make sure the blast binary path exists and is a directory
unless ($blast_bin_dir =~ /.*\/$/){
  $blast_bin_dir = "$blast_bin_dir"."/";
}
unless (-e $blast_bin_dir && -d $blast_bin_dir){
  print RED, "\nBlast binary directory: $blast_bin_dir does not appear valid!\n\n", RESET;
  close (LOG);
  exit();
}

#Create a multiple cpus command if specified by the user
my $multi_cpu_option = '';
if ($multiple_cpus){
  $multi_cpu_option = " -a $multiple_cpus";
}

#Create a working temp directory
my $working_dir;
if ($force){
  if ($force eq "yes"){
    $working_dir = &createNewDir('-path'=>$temp_dir, '-new_dir_name'=>$job_name, '-force'=>$force);
  }else{
    print RED, "\nValue supplied for option: --force ($force) not understood!!", RESET;
    close (LOG);
    exit();
  }
}else{
  $working_dir = &createNewDir('-path'=>$temp_dir, '-new_dir_name'=>$job_name);
}

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

#If the user specified a big mem job, get all gene sequences at once
my $all_genes_ref;
my %gene_list;
unless ($big_mem){
  $big_mem = "no";
}
if ($big_mem eq "yes"){
  print BLUE, "\n\nPreparing to get all gene sequences at once - getting gene_id list from input probe file", RESET;
  open (PROBE, "$probe_file") || die "\nCould not open input probe file: $probe_file\n\n";
  while(<PROBE>){
    chomp($_);
    my @line = split ("\t", $_);
    if ($line[2] =~ /^(\d+)/){
      $gene_list{$1}{tmp} = '';
    }
  }
  close (PROBE);
  my $genes_found = keys %gene_list;
  my @gene_ids = keys %gene_list;

  print BLUE, "\n\tFound $genes_found genes, getting all gene sequences at once\n", RESET;
  $all_genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"yes", '-silent'=>"no");
  $alexa_dbh->disconnect();
}

#Open master output file
my $new_column1 = "withinGene_Non-TargetHits";
my $new_column2 = "withinGene_TargetHits";
my $new_column3 = "withinGene_largestNon-TargetHitLength";
my $new_column4 = "withinGene_largestTargetHitLength";
open (MASTER_OUT, ">$outfile") || die "\nCould not open output file: $outfile";

#1.) Import probes from an input probe file
#    - NOTE this script assumes that probes are ordered according to to their target gene ID
#    - Get all the probes for a single gene, then proceed with the processing of it
#    - This is done to save the memory required to store all probe sequences at once
#    - Store probe coordinate mapping info for use in determining target vs. non-target hits
my %columns;

my $header_line;
my $first_line = 1;

my $probe_id_count = 0;
my $gene_count = 0;

my $probes_processed = 0;
my $genes_processed = 0;

#First count the total number of probes
print BLUE, "\nCounting the total number of probes in the input file\n\n", RESET;
print LOG "\nCounting the total number of probes in the input file\n\n";

open (PROBE, "$probe_file") || die "\nCould not open input probe file: $probe_file\n\n";
my $total_number_of_probes = 0;
while(<PROBE>){
  if ($first_line == 1){
    $first_line = 0;
    next();
  }
  $total_number_of_probes++;
}
close (PROBE);

print BLUE, "\nFound $total_number_of_probes probes in the input file\n\n", RESET;
print LOG "\nFound $total_number_of_probes probes in the input file\n\n";

#Start parsing the probe file
print BLUE, "\nBegin processing the input probe file: $probe_file\n\n", RESET;
print LOG "\nBegin processing the input probe file: $probe_file\n\n";

#Open the probe files and read the neccessary probe data into a hash keyed gene_id.  Store probe and probeset IDs
$first_line = 1;
open (PROBE, "$probe_file") || die "\nCould not open input probe file: $probe_file\n\n";

#Change working directory to the temp dir specified by the user
chdir($working_dir);

my $working_gene_id = 0;
my $total_probe_count = 0;
my %probes;
while(<PROBE>){

  #Deal with the header line
  if ($first_line == 1){
    $header_line = $_;
    chomp ($header_line);

    my @columns = split("\t", $header_line);
    my $col_count = 0;
    foreach my $column (@columns){
      $columns{$column}{column_pos} = $col_count;
      $col_count++;
    }

    #Check for critical columns and their names
    unless ($columns{'Probe_Count'} && $columns{'Gene_ID'} && $columns{'Sequence'} && $columns{'Unit1_start'} && $columns{'Unit1_end'} && $columns{'Unit2_start'} && $columns{'Unit2_end'}){
      print RED, "\nCritical column missing or named incorrectly in probe file: $probe_file, check input file\n\n", RESET;
      close (LOG);
      $alexa_dbh->disconnect();
      exit();
    }

    #Print the header line
    print MASTER_OUT "$header_line\t$new_column1\t$new_column2\t$new_column3\t$new_column4\n";

    $first_line = 0;
    next();
  }
  $total_probe_count++;
  chomp($_);
  my @data_line = split ("\t", $_);

  my $line_record = $_;
  my $gene_id = $data_line[$columns{'Gene_ID'}{column_pos}];
  my $probe_id = $data_line[$columns{'Probe_Count'}{column_pos}];
  my $sequence = $data_line[$columns{'Sequence'}{column_pos}];
  my $unit1_start = $data_line[$columns{'Unit1_start'}{column_pos}];
  my $unit1_end = $data_line[$columns{'Unit1_end'}{column_pos}];
  my $unit2_start = $data_line[$columns{'Unit2_start'}{column_pos}];
  my $unit2_end = $data_line[$columns{'Unit2_end'}{column_pos}];

  unless ($gene_id =~ /\d+/){
    print RED, "\nGene ID: $gene_id is not valid",RESET;
    print RED, "\nRemember that this test should not be applied to negative control probes\n\n", RESET;
    $alexa_dbh->disconnect();
    exit();
  }

  #Add probe to hash if it is part of the current working gene
  if($working_gene_id eq '0'){
    #First probe entry
    $working_gene_id = $gene_id;
    $probes{$probe_id}{seq} = $sequence;
    $probes{$probe_id}{unit1_start} = $unit1_start;
    $probes{$probe_id}{unit1_end} = $unit1_end;
    $probes{$probe_id}{unit2_start} = $unit2_start;
    $probes{$probe_id}{unit2_end} = $unit2_end;
    $probes{$probe_id}{line} = $line_record;
    $gene_count++;
    $probe_id_count++;
  }elsif($total_probe_count == $total_number_of_probes){
    #Last probe entry, process this last gene
    $probes{$probe_id}{seq} = $sequence;
    $probes{$probe_id}{unit1_start} = $unit1_start;
    $probes{$probe_id}{unit1_end} = $unit1_end;
    $probes{$probe_id}{unit2_start} = $unit2_start;
    $probes{$probe_id}{unit2_end} = $unit2_end;
    $probes{$probe_id}{line} = $line_record;

    &processGene('-probes'=>\%probes, '-gene_id'=>$working_gene_id);

    $probe_id_count++;

  }elsif ($gene_id == $working_gene_id){
    #Continuing with the same gene
    $probes{$probe_id}{seq} = $sequence;
    $probes{$probe_id}{unit1_start} = $unit1_start;
    $probes{$probe_id}{unit1_end} = $unit1_end;
    $probes{$probe_id}{unit2_start} = $unit2_start;
    $probes{$probe_id}{unit2_end} = $unit2_end;
    $probes{$probe_id}{line} = $line_record;

    $probe_id_count++;

  }else{
    #Not the last probe but this probe is for a new gene, so process the old one
    &processGene('-probes'=>\%probes, '-gene_id'=>$working_gene_id);

    #Initialize variables and start next gene
    %probes = ();
    $working_gene_id = $gene_id;
    $probes{$probe_id}{seq} = $sequence;
    $probes{$probe_id}{unit1_start} = $unit1_start;
    $probes{$probe_id}{unit1_end} = $unit1_end;
    $probes{$probe_id}{unit2_start} = $unit2_start;
    $probes{$probe_id}{unit2_end} = $unit2_end;
    $probes{$probe_id}{line} = $line_record;

    $gene_count++;
    $probe_id_count++;
  }
}
close (PROBE);

print BLUE, "\nFound $probe_id_count probes corresponding to $gene_count genes\n\n", RESET;
print BLUE, "\nA total of $probes_processed probes were written to fasta files, and these corresponded to $genes_processed genes processed\n\n", RESET;
print LOG "\nFound $probe_id_count probes corresponding to $gene_count genes\n\n";
print LOG "\nA total of $probes_processed probes were written to fasta files, and these corresponded to $genes_processed genes processed\n\n";

#Clean up the working directory
my $command = "rm -fr $working_dir";
system ($command);

#Close database connection
unless ($big_mem eq "yes"){
  $alexa_dbh->disconnect();
}

close (MASTER_OUT);
close (LOG);

exit();


#############################################################################################################################
#Process gene
#############################################################################################################################
sub processGene{
  my %args = @_;
  my $probes_ref = $args{'-probes'};
  my $working_gene_id = $args{'-gene_id'};

  my $probe_count = keys %{$probes_ref};

  print CYAN, "\n\nGene count: $gene_count", RESET;
  print LOG "\n\nGene count: $gene_count";

  #Process the probes found for the last gene, initialize variables and start the next gene
  print BLUE, "\n\t1.) Processing gene: $working_gene_id (with $probe_count probes)", RESET;
  print LOG "\n\t1.) Processing gene: $working_gene_id";

  #Process probes for last gene (working_gene_id)

  #2.) Create a query fasta file of all probe sequences in a temp directory
  print BLUE "\n\t2.) Creating query fasta file", RESET;

  my $probe_fasta_file = "$working_dir"."probes.fa";
  open (PROBE_FASTA, ">$probe_fasta_file") || die "\nCould not open fasta file: $probe_fasta_file\n\n";
  foreach my $p_id (sort {$a <=> $b} keys %probes){
    print PROBE_FASTA ">$p_id\n$probes{$p_id}{seq}\n";
    $probes_processed++;
  }
  close (PROBE_FASTA);

  #3.) Get the full genomic sequence for this target gene and create a fasta file for in a temp directory
  my $gene_seq;
  if ($big_mem eq "yes"){
    print BLUE, "\n\t3.) Getting gene sequence from previously created hash", RESET;
    $gene_seq = $all_genes_ref->{$working_gene_id}->{sequence};
  }else{
    print BLUE, "\n\t3.) Getting gene sequence from ALEXA", RESET;
    my @gene_ids;
    push (@gene_ids, $working_gene_id);
    my $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids,'-sequence'=>"yes", '-silent'=>"yes");
    $gene_seq = $genes_ref->{$working_gene_id}->{sequence};
  }

  #Make sure a valid gene_seq was returned or die - in case of bizarre mysql or DBI error
  unless ($gene_seq){
    print RED, "\nGene seq was not defined.  Mysql connection error??\n\n", RESET;
    exit();
  }

  my $gene_length = length($gene_seq);

  my $gene_fasta_file = "$working_dir"."gene.fa";
  open (GENE_FASTA, ">$gene_fasta_file") || die "\nCould not open fasta file: $gene_fasta_file\n\n";
  print GENE_FASTA ">$working_gene_id\n$gene_seq\n";
  close (GENE_FASTA);

  #4.) Use formatDB to create a blastable database for this single genomic sequence
  print BLUE, "\n\t4.) Creating blast database with formatdb (gene length is $gene_length bp)", RESET;
  my $formatdb_command = "$blast_bin_dir"."formatdb -t gene_seq -i $gene_fasta_file -p F -o T -n gene_seq";
  system($formatdb_command);

  #5.) Blast all probes in the query sequence against this temp database and create a blast results file in the temp directory
  print BLUE, "\n\t5.) Running blast command", RESET;
  my $blast_database = "$working_dir"."gene_seq";
  my $blast_results_file = "$working_dir"."blast_results.txt";
  my $blastall_command = "$blast_bin_dir"."blastall -p blastn -d $blast_database -i $probe_fasta_file -o $blast_results_file -F F -W $word_size -m 8"."$multi_cpu_option";
  system($blastall_command);

  #6.) Parse the blast results file
  print BLUE, "\n\t6.) Parsing blast results file", RESET;
  my $blast_results_ref = &loadBlastResults('-blast_results_file'=>$blast_results_file);

  #7.) Test the position of hits to the target gene and compare to those expected
  #    - Identify hits that do not overlap with the coordinates where the probe was actually extracted from as non-target hits
  #    - Hits that do overlap are target hits
  #    - Summarize these findings to a master output file
  print BLUE, "\n\t7.) Conducting specificity test", RESET;
  &testProbeSpecificity('-probe_object'=>$probes_ref, '-blast_object'=>$blast_results_ref);

  #8.) Print the probe specificity test results to the master out file
  print BLUE, "\n\t8.) Printing results", RESET;
  &printResults('-probe_object'=>$probes_ref);

  $genes_processed++;
  return();
}


#########################################################################################################
#6.) Load all BLAST results from the blast-results files                                                #
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
  print BLUE, "\n\t\tProcessing $results_file", RESET;
  print LOG "\n\t\tProcessing $results_file";
  open (BLAST_RESULTS, "$results_file") || die "\nCould not open blast results file: $results_file\n\n";

  while (<BLAST_RESULTS>){
    my @line = split ("\t", $_);
    my $probe_id = $line[0];
    my $hit_accession = $line[1];
    my $percent_identity = $line[2];
    my $alignment_length = $line[3];
    my $subject_start = $line[8]; #Position on the mRNA/EST sequence at which the probe alignment begins
    my $subject_end = $line[9]; #Position on the mRNA/EST sequence at which the probe alignment ends

    $blast_results_processed++;

    #Note: using the following datastructure will not allow multiple hits by one probe to the same accession to be recorded
    #Only the more significant hit (longer) will be stored when this this occurs

    #See if there are any hits recorded for this probe to this gene already
    if ($blast_results_object{$probe_id}){
      my $hits_ref = $blast_results_object{$probe_id}{blast_hits};
      $blast_results_object{$probe_id}{hit_count}++;

      my $hit_count = $blast_results_object{$probe_id}{hit_count};

      $hits_ref->{$hit_count}->{pi} = $percent_identity;
      $hits_ref->{$hit_count}->{al} = $alignment_length;
      $hits_ref->{$hit_count}->{ss} = $subject_start;
      $hits_ref->{$hit_count}->{se} = $subject_end;
      $blast_results_stored++;

    }else{
      #create the first hit record for this probe
      my %hits;
      my $hit_count = 1;
      $hits{$hit_count}{pi} = $percent_identity;
      $hits{$hit_count}{al} = $alignment_length;
      $hits{$hit_count}{ss} = $subject_start;
      $hits{$hit_count}{se} = $subject_end;
      $blast_results_object{$probe_id}{blast_hits} = \%hits;
      $blast_results_object{$probe_id}{hit_count} = $hit_count;
      $blast_results_stored++;
    }
  }
  close BLAST_RESULTS;

  print BLUE, "\n\t\tProcessed $blast_results_processed blast results", RESET;
  print BLUE, "\n\t\tActually stored $blast_results_stored blast results", RESET;
  print LOG "\n\t\tProcessed $blast_results_processed blast results";
  print LOG "\n\t\tActually stored $blast_results_stored blast results";

  return(\%blast_results_object);
}


#########################################################################################################
#7.) Using the infomation gathered in step 6, conduct the actual specificity test                       #
#########################################################################################################
sub testProbeSpecificity{
  my %args = @_;
  my $probe_object_ref = $args{'-probe_object'};
  my $blast_object_ref = $args{'-blast_object'};

  print BLUE, "\n\t\tBegin actual specificity test", RESET;
  print LOG "\n\t\tBegin actual specificity test";

  #Go through each probe ID and test its specificity to the target locus
  foreach my $probe_id (sort {$a <=> $b} keys %{$probe_object_ref}){

    #Count the number of times a probe hits part of the target gene where the probe was actually extracted from
    #Also count hits to parts of the target gene with no-overlap to the expected region
    my $non_target_hits = 0;
    my $target_hits = 0;

    #Initialize variables
    $probe_object_ref->{$probe_id}->{non_target_hits} = $non_target_hits;
    $probe_object_ref->{$probe_id}->{target_hits} = $target_hits;

    #Keep track of the longest non-target and target alignment lengths for each probe
    $probe_object_ref->{$probe_id}->{largest_non_target_al} = 0;
    $probe_object_ref->{$probe_id}->{largest_target_al} = 0;

    #Get the blast hits for this probe (as a hash keyed on probe ID)
    #Check for probes with no blast hits
    my $blast_hits_ref;
    if ($blast_object_ref->{$probe_id}->{blast_hits}){
      $blast_hits_ref = $blast_object_ref->{$probe_id}->{blast_hits};
    }else{
      next();
    }

    #Get the region(s) actually targetted by this probe
    my $unit1_start = $probe_object_ref->{$probe_id}->{unit1_start};
    my $unit1_end = $probe_object_ref->{$probe_id}->{unit1_end};
    my $unit2_start = $probe_object_ref->{$probe_id}->{unit2_start};
    my $unit2_end = $probe_object_ref->{$probe_id}->{unit2_end};

    #Go through each BLAST hit and make sure it overlaps the region actually targeted
    BLAST_HIT:foreach my $hit_count (sort keys %{$blast_hits_ref}){
	my $overlap = 0;

	my $subject_start = $blast_hits_ref->{$hit_count}->{ss};
	my $subject_end = $blast_hits_ref->{$hit_count}->{se};

	#Unit1: make sure the coordinates are actually defined - then look for overlap between this hit and the target coordinates
	if ($unit1_start =~ /\d+/ && $unit1_end =~ /\d+/){

	  if ($subject_start >= $unit1_start && $subject_start <= $unit1_end) {$overlap = 1;}
	  if ($subject_end >= $unit1_start && $subject_end <= $unit1_end) {$overlap = 1;}
	  if ($subject_start <= $unit1_start && $subject_end >= $unit1_end) {$overlap = 1;}

	}else{
	  print RED, "\nUnit1 coordinates do not seem valid for probe: $probe_id\n\n", RESET;
	  close (LOG);
	  $alexa_dbh->disconnect();
	  exit();
	}

	#Unit2: make sure the coordinates are actually defined (only junction probes have unit2 coordinates
	if ($unit2_start =~ /\d+/ && $unit2_end =~ /\d+/){

	  if ($subject_start >= $unit2_start && $subject_start <= $unit2_end) {$overlap = 1;}
	  if ($subject_end >= $unit2_start && $subject_end <= $unit2_end) {$overlap = 1;}
	  if ($subject_start <= $unit2_start && $subject_end >= $unit2_end) {$overlap = 1;}

	}

	if ($overlap == 1){
	#If overlap was found, this is a target hit
	  $target_hits++;

	  #update the largest target hit found so far if neccessary
	  if ($blast_hits_ref->{$hit_count}->{al} > $probe_object_ref->{$probe_id}->{largest_target_al}){
	    $probe_object_ref->{$probe_id}->{largest_target_al} = $blast_hits_ref->{$hit_count}->{al};
	  }
	}else{
	  #Otherwise it is a non-target hit
	  $non_target_hits++;

	  #update the largest non-target hit found so far if neccessary
	  if ($blast_hits_ref->{$hit_count}->{al} > $probe_object_ref->{$probe_id}->{largest_non_target_al}){
	    $probe_object_ref->{$probe_id}->{largest_non_target_al} = $blast_hits_ref->{$hit_count}->{al};
	  }
	}
      }
    $probe_object_ref->{$probe_id}->{non_target_hits} = $non_target_hits;
    $probe_object_ref->{$probe_id}->{target_hits} = $target_hits;
  }
  return();
}


#########################################################################################################
#8.) Print the probe specificity test results to the master out file                                    #
#########################################################################################################
sub printResults{
  my %args = @_;

  my $probe_object_ref = $args{'-probe_object'};

  print BLUE, "\n\t\tPrinting results to $outfile\n", RESET;
  print LOG "\n\t\tPrinting results to $outfile\n";

  foreach my $probe_id (sort {$a <=> $b} keys %{$probe_object_ref}){

    #check the format of each variable - final sanity check
    unless ($probe_object_ref->{$probe_id}->{non_target_hits} =~ /\d+/){
      print RED, "\nUndefined non-target_hits", RESET;
      close (LOG);
      $alexa_dbh->disconnect();
      exit();
    }
    unless ($probe_object_ref->{$probe_id}->{target_hits} =~ /\d+/){
      print RED, "\nUndefined target_hits", RESET;
      close (LOG);
      $alexa_dbh->disconnect();
      exit();
    }
    unless ($probe_object_ref->{$probe_id}->{largest_non_target_al} =~ /\d+/){
      print RED, "\nUndefined largest_non_target_al", RESET;
      close (LOG);
      $alexa_dbh->disconnect();
      exit();
    }
    unless ($probe_object_ref->{$probe_id}->{largest_target_al} =~ /\d+/){
      print RED, "\nUndefined largest_target_al", RESET;
      close (LOG);
      $alexa_dbh->disconnect();
      exit();
    }

    print MASTER_OUT "$probe_object_ref->{$probe_id}->{line}\t$probe_object_ref->{$probe_id}->{non_target_hits}\t$probe_object_ref->{$probe_id}->{target_hits}\t$probe_object_ref->{$probe_id}->{largest_non_target_al}\t$probe_object_ref->{$probe_id}->{largest_target_al}\n";
  }

  return();
}

