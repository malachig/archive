#!/usr/bin/perl -w
#Written by Malachi Griffith

#The purpose of this script is to consider a large set of probes and eliminate probes that do not pass certain quality checks
#These probes are taken from an input file and certain header naming formats are assumed

#Note probes that have been extracted to probe files have already passed certain basic requirements such as avoiding unknown genomic bases, avoiding
#excessive repeatMasked bases (more than 1/4 of the probe length) and avoiding extraction of a sequence from an exon that is too short to accomodate it.
#Also the probes may have been pre-filtered for Tm (in particular for exon and intron probes where a larger selection region means we can be more stringent)

#1.) Load input probe file

#2.) The following quality checks will be considered (in this order):
#2-a.) Probes must have a Tm equal to '--target_tm' +/- '--tm_range'  (for example, 67C +/- 3C) - determined by nearest neighbour Tm calculation
#2-b.) Probes must not contain low complexity regions larger than '--complexity_region_limit' (mono-, di-, or tri-nucleotide repeats) - determined by mdust
#2-c.) Probes must have self-self folding free energy scores of less than the user specified threshhold '--pairfold_score_limit' - Determined by Pairfold
#2-d.) Probes must have internal folding free energy scores of less than the use specified threshhold '--simfold_score_limit' - Determined by SimFold
#2-e.) Probes must be specific to the targeted locus - determined by blast against mRNAs
#      - The allowable hit size between a probe and non-target mRNA is specified by '--mrna_hit_percent_limit'
#      - Since the probes can vary in length, only fail probes if they have a non-specific hit of say 50% of their length or greater
#2-f.) Probes must be specific to the targeted locus - determined by blast against ESTs
#      - The allowable hit size between a probe and non-target EST is specified by '--est_hit_percent_limit'
#      - Since the probes can vary in length, only fail probes if they have a non-specific hit of say 50% of their length or greater
#2-g.) Probes must be specific to the targeted locus - determined by blast against all ensembl transcripts
#      - The allowable hit size between a probe and non-target Ensembl transcript is specified by '--enst_hit_percent_limit'
#      - Since the probes can vary in length, only fail probes if they have a non-specific hit of say 50% of their length or greater
#2-h.) Eliminate identical or near idential probes (for example 30 bp match over a 36 bp probe or more will fail)
#      - The allowable hit size between a probe and another probe from a different probeset is specified by '--probe_hit_percent_limit'
#      - Since the probes can vary in length, only fail probes if they have a non-specific hit of say 75% of their length or greater
#2-i)  Eliminate probes with significant hits to the entire human genome sequence (only used for negative-control probes)
#2-j.) Eliminate probes with significant hits to the target gene that are outside the targeted region (i.e. caused by within gene repeats)
#2-k.) Test probe uniqueness within each probeset.  Since all probes are blasted against the complete set of experimental probes the probe specificity
#      test will eliminate identical or nearly identical probes.  However, this test does not count probes within the same probeset
#      - Probes within the same probeset may be highly similar to each other, this is expected, however we want to make sure that none are actually identical

#5.) Print a probe file containing only the filtered probes
#As each of these steps is accomplished - summarize the number of probes that pass/fail each test

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Initialize input parameters needed for filtering from the user
my $infile = '';
my $target_length = '';
my $target_tm = '';
my $tm_range = '';
my $complexity_region_limit = '';
my $simfold_score_limit = '';
my $pairfold_score_limit = '';
my $mrna_hit_percent_limit = '';
my $est_hit_percent_limit = '';
my $enst_hit_percent_limit = '';
my $probe_hit_percent_limit = '';
my $genomic_hit_percent_limit = '';
my $withinGene_hit_percent_limit = '';
my $probeset_size = '';
my $outfile = '';
my $logfile = '';
my $verbose = '';

GetOptions ('infile=s'=>\$infile, 'target_length=i'=>\$target_length,
	    'target_tm=f'=>\$target_tm, 'tm_range=f'=>\$tm_range, 'complexity_region_limit=i'=>\$complexity_region_limit,
	    'simfold_score_limit=f'=>\$simfold_score_limit, 'pairfold_score_limit=f'=>\$pairfold_score_limit,
	    'mrna_hit_percent_limit=f'=>\$mrna_hit_percent_limit, 'est_hit_percent_limit=f'=>\$est_hit_percent_limit,
	    'enst_hit_percent_limit=f'=>\$enst_hit_percent_limit, 'probe_hit_percent_limit=f'=>\$probe_hit_percent_limit,
	    'genomic_hit_percent_limit=f'=>\$genomic_hit_percent_limit, 'withinGene_hit_percent_limit=f'=>\$withinGene_hit_percent_limit,
	    'probeset_size=i'=>\$probeset_size, 'outfile=s'=>\$outfile, 'logfile=s'=>\$logfile, 'verbose=s'=>\$verbose);

#Provide instruction to the user
print GREEN, "\n\nSpecify the input probe file using: --infile", RESET;
print GREEN, "\nSpecify the desired name of the filtered output probe file using: --outfile", RESET;
print GREEN, "\nSpecify the desired name of the out logfile which will act as a record of the filtering using: --logfile", RESET;
print GREEN, "\nThis script conducts a series of quality checks on probe records from this input file and produces a filtered output file", RESET;
print GREEN, "\nProbes that do not meet the required thresholds are eliminated", RESET;
print GREEN, "\nYou must specify the following parameters for filtering the probes:", RESET;
print GREEN, "\n\tTarget length using: --target_length (say 36 bp)", RESET;
print GREEN, "\n\tTarget Tm using: --target_tm (say 67.0 degrees)", RESET;
print GREEN, "\n\tAcceptable range of Tm around this target using: --tm_range (say 3.0 degrees)", RESET;
print GREEN, "\n\tThe minimum sized region of low complexity that will cause a probe to fail using: --complexity_region_limit (say 6 bp)", RESET;
print GREEN, "\n\tThe cutoff SimFold score (within-probe/internal folding) that will cause a probe to fail using: --simfold_score_limit (say -7.0)", RESET;
print GREEN, "\n\tThe cutoff PairFold score (self-self/between probe folding) that will cause a probe to fail using: --pairfold_score_limit (say -19.5)", RESET;
print GREEN, "\n\tThe mRNA specifity hit percent length cutoff using: --mrna_hit_percent_limit (say 50.0%)", RESET;
print GREEN, "\n\tThe EST specifity hit percent length cutoff using: --est_hit_percent_limit (say 50.0%)", RESET;
print GREEN, "\n\tThe Ensembl specifity hit percent length cutoff using: --enst_hit_percent_limit (say 50.0%)", RESET;
print GREEN, "\n\tThe Probe specifity hit percent length cutoff using: --probe_hit_percent_limit  (say 75.0%)", RESET;
print GREEN, "\n\tThe Genomic specificity hit percent length cutoff using: --genomic_hit_percent_limit (say 50.0%)", RESET;
print GREEN, "\n\tThe Within gene specificity hit percent length cutoff using:  --withinGene_hit_percent_limit (say 50.0%)", RESET;
print GREEN, "\n\tThe desired probeset size (number of probes per single exon, junction, etc.) using: --probeset_size (say 2)", RESET;
print GREEN, "\n\tFor verbose output, use: verbose=yes\n\n", RESET;
print GREEN, "\n\nExample filterProbes.pl --infile=exonJunctionProbes.txt  --target_length=36  --target_tm=67.0  --tm_range=3.5  --complexity_region_limit=6  --simfold_score_limit=-7.0  --pairfold_score_limit=-19.5  --mrna_hit_percent_limit=65.0  --est_hit_percent_limit=65.0  --enst_hit_percent_limit=65.0  --probe_hit_percent_limit=65.0  --genomic_hit_percent_limit=65.0  --withinGene_hit_percent_limit=65.0  --probeset_size=2  --outfile=exonJunctionProbes_filtered.txt --logfile=filter_exonJunction_LOG.txt  --verbose=no\n\n", RESET;

unless ($infile && $target_length && $target_tm && $tm_range && $complexity_region_limit && $simfold_score_limit && $pairfold_score_limit && $mrna_hit_percent_limit && $est_hit_percent_limit && $enst_hit_percent_limit && $probe_hit_percent_limit && $genomic_hit_percent_limit && $withinGene_hit_percent_limit && $probeset_size && $outfile && $logfile && $verbose){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Open logfile
open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile for output\n\n";

#Print out the parameters selected by the user so that this information can be explicitly stored in a LOG file
print BLUE, "\n\nPARAMETERS USED:\ninfile = $infile\ntarget_length = $target_length\ntarget_tm = $target_tm\ntm_range = $tm_range\ncomplexity_region_limit = $complexity_region_limit\nsimfold_score_limit = $simfold_score_limit\npairfold_score_limit = $pairfold_score_limit\nmrna_hit_percent_limit = $mrna_hit_percent_limit\nest_hit_percent_limit = $est_hit_percent_limit\nenst_hit_percent_limit = $enst_hit_percent_limit\nprobe_hit_percent_limit = $probe_hit_percent_limit\ngenomic_hit_percent_limit = $genomic_hit_percent_limit\nwithinGene_hit_percent_limit = $withinGene_hit_percent_limit\nprobeset_size = $probeset_size\noutfile = $outfile\nlogile = $logfile\nverbose = $verbose\n\n", RESET;
print LOG "\n\nPARAMETERS USED:\ninfile = $infile\ntarget_length = $target_length\ntarget_tm = $target_tm\ntm_range = $tm_range\ncomplexity_region_limit = $complexity_region_limit\nsimfold_score_limit = $simfold_score_limit\npairfold_score_limit = $pairfold_score_limit\nmrna_hit_percent_limit = $mrna_hit_percent_limit\nest_hit_percent_limit = $est_hit_percent_limit\nenst_hit_percent_limit = $enst_hit_percent_limit\nprobe_hit_percent_limit = $probe_hit_percent_limit\ngenomic_hit_percent_limit = $genomic_hit_percent_limit\nwithinGene_hit_percent_limit = $withinGene_hit_percent_limit\nprobeset_size = $probeset_size\noutfile = $outfile\nlogile = $logfile\nverbose = $verbose\n\n";

my $probe_count = 0;

#1.) Load the input probes from an input file or from ALEXA
my $header_line;
my %probes;
my %probesets;
my %columns;

&loadProbeInfo_file('-input_file'=>$infile);

#2.) Conduct each test on each probe, update the 'status' field for probes that fail and summarize the number of each type of failure
print BLUE, "\n\nProceeding with quality tests on all probes:", RESET;
print LOG "\n\nProceeding with quality tests on all probes:";

#2-a.) Tm test
if ($columns{Probe_Tm}){
  &testProbeTms('-target_tm'=>$target_tm, '-tm_range'=>$tm_range);
}else{
  print YELLOW, "\n\tProbe_Tm values not found in input file, skipping the Tm test!\n", RESET;
  print LOG "\n\tProbe_Tm values not found in input file, skipping the Tm test!\n";
}

#2-b) Low complexity test
if ($columns{MdustBases_cutoff_11}){
  &testProbeComplexity('-low_complexity_limit'=>$complexity_region_limit);
}else{
  print YELLOW, "\n\tMdustBases_cutoff_11 values not found in input file, skipping the low complexity test!\n", RESET;
  print LOG "\n\tMdustBases_cutoff_11 values not found in input file, skipping the low complexity test!\n";
}

#2-c) Pairfold - self-self folding test
if ($columns{PairFold_score}){
  &testPairfoldScore('-pairfold_score_limit'=>$pairfold_score_limit);
}else{
  print YELLOW, "\n\tPairFold_score values not found in input file, skipping the pairfold score test!\n", RESET;
  print LOG "\n\tPairFold_score values not found in input file, skipping the pairfold score test!\n";
}

#2-d) Simfold - internal folding test
if ($columns{SimFold_score}){
  &testSimfoldScore('-simfold_score_limit'=>$simfold_score_limit);
}else{
  print YELLOW, "\n\tSimFold_score values not found in input file, skipping the simfold score test!\n", RESET;
  print LOG "\n\tSimFold_score values not found in input file, skipping the simfold score test!\n";
}

#2-e) Specificity test - Do not allow any hits to mRNAs that map outside the target region
if ($columns{'mRNA_largestNon-TargetHitLength'}){
  &testProbeSpecificity_mRNA('-non_target_hit_percent_limit'=>$mrna_hit_percent_limit);
}else{
  print YELLOW, "\n\tmRNA_largestNon-TargetHitLength values not found in input file, skipping the mRNA specificity test!\n", RESET;
  print LOG "\n\tmRNA_largestNon-TargetHitLength values not found in input file, skipping the mRNA specificity test!\n";
}

#2-f) Specificity test - Do not allow any hits to ESTs that map outside the target region (not used for all probe types)
if ($columns{'EST_largestNon-TargetHitLength'}){
  &testProbeSpecificity_EST('-non_target_hit_percent_limit'=>$est_hit_percent_limit);
}else{
  print YELLOW, "\n\tEST_largestNon-TargetHitLength values not found in input file, skipping the EST specificity test!\n", RESET;
  print LOG "\n\tEST_largestNon-TargetHitLength values not found in input file, skipping the EST specificity test!\n";
}

#2-g) Specificity test - Do not allow any hits to ensembl transcripts other than those associated with the target gene
if ($columns{'enst_largestNon-TargetHitLength'}){
  &testProbeSpecificity_enst('-non_target_hit_percent_limit'=>$enst_hit_percent_limit);
}else{
  print YELLOW, "\n\tenst_largestNon-TargetHitLength values not found in input file, skipping the Ensembl transcript specificity test!\n", RESET;
  print LOG "\n\tenst_largestNon-TargetHitLength values not found in input file, skipping the Ensembl transcript specificity test!\n";
}

#2-h) Specificity test - Do not allow any hits to other probes (ie. fail identical or near identical probes)
if ($columns{'probe_largestNon-TargetHitLength'}){
  &testProbeSpecificity_probe('-non_target_hit_percent_limit'=>$probe_hit_percent_limit);
}else{
  print YELLOW, "\n\tprobe_largestNon-TargetHitLength values not found in input file, skipping the Probe specificity test!\n", RESET;
  print LOG "\n\tprobe_largestNon-TargetHitLength values not found in input file, skipping the Probe transcript specificity test!\n";
}

#2-i) Specificity test - Do not allow any hits to the human genome (only used for negative-control probes)
if ($columns{'genomic_largestNon-TargetHitLength'}){
  &testProbeSpecificity_genomic('-non_target_hit_percent_limit'=>$genomic_hit_percent_limit);
}else{
  print YELLOW, "\n\tgenomic_largestNon-TargetHitLength values not found in input file, skipping the genomic specificity test!\n", RESET;
  print LOG "\n\tgenomic_largestNon-TargetHitLength values not found in input file, skipping the genomic specificity test!\n";
}

#2-j.) Specificity test - Do not allow any hits to a region of the targeted gene other than the part actually targeted
if ($columns{'withinGene_largestNon-TargetHitLength'}){
  &testProbeSpecificity_withinGene('-non_target_hit_percent_limit'=>$withinGene_hit_percent_limit);
}else{
  print YELLOW, "\n\twithinGene_largestNon-TargetHitLength values not found in input file, skipping the withinGene specificity test!\n", RESET;
  print LOG "\n\twithinGene_largestNon-TargetHitLength values not found in input file, skipping the withinGene specificity test!\n";
}

#2-k.) Test probe uniquenes within probeset
&testProbeUniquenessWithinProbeset();

#3.) Probeset size test - Probesets with less than the desired minimum number of probes will be noted.
#  - Will not be used to pass or fail probes, simply used to give an idea how the filtering is working                   #
&testProbesetSize('-probeset_size'=>$probeset_size);


#7.) Print out files containing the filtered probes
#Go back to the original input file and print each line out to a new file - if the probe contained on that line passed the tests
&printOutputFile('-input_file'=>$infile, '-output_file'=>$outfile);

close (LOG);

exit();


#########################################################################################################
#1.) Load all probes info from an input probe file and build a probe object as a hash                   #
#########################################################################################################
sub loadProbeInfo_file{
  my %args = @_;
  my $probe_file = $args{'-input_file'};

  my $progress_count = 0;
  my $blocks_imported = 0;

  #Open the probe file and read the neccessary probe data into a hash keyed on probe ID
  open (PROBES, "$probe_file") || die "\nCould not open input probe file: $probe_file\n\n";

  #Possible column headings:
  #Probe_Count,ProbeSet_ID,Gene_ID,Sequence,Probe_length,Probe_Tm,Probe_Type
  #Exon1_IDs,Unit1_start,Unit1_end,Exon2_IDs,Unit2_start,Unit2_end,
  #masked_bases,SimFold_score,PairFold_score,MdustBases_cutoff_11
  #mRNA_Non-TargetHits,mRNA_TargetHits,mRNA_largestNon-TargetHitLength,mRNA_largestTargetHitLength
  #EST_Non-TargetHits,EST_TargetHits,EST_largestNon-TargetHitLength,EST_largestTargetHitLength
  #enst_Non-TargetHits,enst_TargetHits,enst_largestNon-TargetHitLength,enst_largestTargetHitLength
  #probe_Non-TargetHits,probe_TargetHits,probe_largestNon-TargetHitLength,probe_largestTargetHitLength
  #genomic_Non-TargetHits,genomic_TargetHits,genomic_largestNon-TargetHitLength,genomic_largestTargetHitLength
  #Exons_Skipped

  my $first_line = 1;

  print BLUE, "\nImporting probe scores from: $probe_file\n\n", RESET;
  print LOG "\nImporting probe scores from: $probe_file\n\n";

  while (<PROBES>){
    $progress_count++;
    if ($progress_count == 10000){
      $blocks_imported++;
      print BLUE, "\n\tImported $blocks_imported blocks of 10,000 probes", RESET;
      print LOG "\n\tImported $blocks_imported blocks of 10,000 probes";
      $progress_count = 0;
    }

    #Get the header line and identify column names and their positions
    if ($first_line == 1){
      $header_line = $_;
      chomp ($header_line);

      my @columns = split("\t", $header_line);
      my $col_count = 0;
      foreach my $column (@columns){
	$columns{$column}{column_pos} = $col_count;
	$col_count++;
      }
      $first_line = 0;

      #Check for critical columns and their names
      unless ($columns{Probe_Count} && $columns{ProbeSet_ID} && $columns{Probe_Type} && $columns{Gene_ID}){
	print RED, "\nCritical column missing or named incorrectly, check input file", RESET;
	exit();
      }
      next();
    }
    #Get the values of interest from each line (probe record)
    chomp($_);
    my @probe_line = split ("\t", $_);

    my $probe_id = $probe_line[$columns{Probe_Count}{column_pos}];
    my $probeset_id = $probe_line[$columns{ProbeSet_ID}{column_pos}];

    unless ($probe_id =~ /^\d+/ && $probeset_id =~ /^\d+/){
      print RED, "\nInvalid probe or probeset ID\n\nLINE: $_\n\n", RESET;
      exit();
    }
    $probe_count++;

    #Store probe in the probe object hash
    $probes{$probe_id}{gene_id} = $probe_line[$columns{Gene_ID}{column_pos}];
    $probes{$probe_id}{probe_type} = $probe_line[$columns{Probe_Type}{column_pos}];
    $probes{$probe_id}{sequence} = $probe_line[$columns{Sequence}{column_pos}];

    #Store probe-to-probeset mapping
    if ($probesets{$probeset_id}){
      my $probe_list_ref = $probesets{$probeset_id}{probes};
      $probe_list_ref->{$probe_id}->{t} = '';
    }else{
      my %probe_list;
      $probe_list{$probe_id}{t} = '';
      $probesets{$probeset_id}{probe_type} = $probe_line[$columns{Probe_Type}{column_pos}];
      $probesets{$probeset_id}{probes} = \%probe_list;
    }

    #Initialize all probes to status='Pass', if they fail any single test, this will be changed to 'Fail'
    $probes{$probe_id}{status} = "Pass";

    #Conduct check of data types to catch missing values or incorrect formats!
    if ($columns{Probe_length}){
      my $probe_length = $probe_line[$columns{Probe_length}{column_pos}];
      if ($probe_length =~ /^\d+/){
	$probes{$probe_id}{probe_length} = $probe_length;
      }else{
	print RED, "\nProbe length value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{Unit1_start}){
      my $unit1_start = $probe_line[$columns{Unit1_start}{column_pos}];
      if ($unit1_start =~ /^\d+/ || $unit1_start =~ /na/){
	$probes{$probe_id}{unit1_start} = $unit1_start;
      }else{
	print RED, "\nUnit1_start value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{Unit1_end}){
      my $unit1_end = $probe_line[$columns{Unit1_end}{column_pos}];
      if ($unit1_end =~ /^\d+/ || $unit1_end =~ /na/){
	$probes{$probe_id}{unit1_end} = $unit1_end;
      }else{
	print RED, "\nUnit1_end value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{Probe_Tm}){
      my $tm_celsius = $probe_line[$columns{Probe_Tm}{column_pos}];
      if ($tm_celsius =~ /\d+\.\d+/){
	$probes{$probe_id}{tm_celsius} = $tm_celsius;
      }else{
	print RED, "\nTm value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{'mRNA_largestNon-TargetHitLength'}){
      my $mrna_non_target_hit_length = $probe_line[$columns{'mRNA_largestNon-TargetHitLength'}{column_pos}];
      if ($mrna_non_target_hit_length =~ /\d+/){
	$probes{$probe_id}{mrna_nt_hl} = $mrna_non_target_hit_length;
      }else{
	print RED, "\nmRNA Non-Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{'EST_largestNon-TargetHitLength'}){
      my $est_non_target_hit_length = $probe_line[$columns{'EST_largestNon-TargetHitLength'}{column_pos}];
      if ($est_non_target_hit_length =~ /\d+/){
	$probes{$probe_id}{est_nt_hl} = $est_non_target_hit_length;
      }else{
	print RED, "\nEST Non-Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{'enst_largestNon-TargetHitLength'}){
      my $enst_non_target_hit_length = $probe_line[$columns{'enst_largestNon-TargetHitLength'}{column_pos}];
      if ($enst_non_target_hit_length =~ /\d+/){
	$probes{$probe_id}{enst_nt_hl} = $enst_non_target_hit_length;
      }else{
	print RED, "\nEnst Non-Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{'probe_largestNon-TargetHitLength'}){
      my $probe_non_target_hit_length = $probe_line[$columns{'probe_largestNon-TargetHitLength'}{column_pos}];
      if ($probe_non_target_hit_length =~ /\d+/){
	$probes{$probe_id}{probe_nt_hl} = $probe_non_target_hit_length;
      }else{
	print RED, "\nProbe Non-Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{'genomic_largestNon-TargetHitLength'}){
      my $genomic_non_target_hit_length = $probe_line[$columns{'genomic_largestNon-TargetHitLength'}{column_pos}];
      if ($genomic_non_target_hit_length =~ /\d+/){
	$probes{$probe_id}{genomic_nt_hl} = $genomic_non_target_hit_length;
      }else{
	print RED, "\ngenomic Non-Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{'withinGene_largestNon-TargetHitLength'}){
      my $withinGene_non_target_hit_length = $probe_line[$columns{'withinGene_largestNon-TargetHitLength'}{column_pos}];
      if ($withinGene_non_target_hit_length =~ /\d+/){
	$probes{$probe_id}{withinGene_nt_hl} = $withinGene_non_target_hit_length;
      }else{
	print RED, "\nwithinGene Non-Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{SimFold_score}){
      my $simfold_score = $probe_line[$columns{SimFold_score}{column_pos}];
      if ($simfold_score =~ /\d+\.\d+/ || $simfold_score =~ /^0$/){
	$probes{$probe_id}{simfold_score} = $simfold_score;
      }else{
	print RED, "\nSimfold Score missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{PairFold_score}){
      my $pairfold_score = $probe_line[$columns{PairFold_score}{column_pos}];
      if ($pairfold_score =~ /\d+\.\d+/ || $pairfold_score =~ /^0$/){
	$probes{$probe_id}{pairfold_score} = $pairfold_score;
      }else{
	print RED, "\nPairfold Score missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{MdustBases_cutoff_11}){
      my $mdust_bases =  $probe_line[$columns{MdustBases_cutoff_11}{column_pos}];
      if ($mdust_bases =~ /^\d+/){
	$probes{$probe_id}{mdust_bases} = $mdust_bases;
      }else{
	print RED, "\nMdust Score missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{Exons_Skipped}){
      my $exons_skipped = $probe_line[$columns{Exons_Skipped}{column_pos}];
      if ($exons_skipped =~ /\d+/ || $exons_skipped eq "na"){
	$probes{$probe_id}{Exons_Skipped} = $exons_skipped;
      }else{
	print RED, "\nExons skipped Score missing - or incorrect format: $exons_skipped\tProbe = $probe_id", RESET;
	exit();
      }
    }
  }

  close (PROBES);

  print BLUE, "\n\nFound $probe_count probes in the input file\n", RESET;
  print LOG "\n\nFound $probe_count probes in the input file\n";
  return ();
}


#########################################################################################################
#2-a.) Tm test                                                                                          #
#########################################################################################################
sub testProbeTms{
  my %args = @_;
  my $tm = $args{'-target_tm'};
  my $range = $args{'-tm_range'};

  my $upper_limit = $tm + $range;
  my $lower_limit = $tm - $range;

  my $tm_fail_count = 0;

  #Go through each probe and check if it is within the range specified from the target Tm specified
  foreach my $probe_id (keys %probes){
    my $probe_tm = $probes{$probe_id}{tm_celsius};

    if ($probe_tm > $upper_limit || $probe_tm < $lower_limit){
      $tm_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($tm_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);

  print BLUE, "\n\tA total of $tm_fail_count ($percent_fail%) Probes failed the Tm test (Tm outside the range: $lower_limit - $upper_limit)\n\n", RESET;
  print LOG "\n\tA total of $tm_fail_count ($percent_fail%) Probes failed the Tm test (Tm outside the range: $lower_limit - $upper_limit)\n\n";

  return();
}


#########################################################################################################
#2-b.) Low Complexity Test                                                                              #
#########################################################################################################
sub testProbeComplexity{
  my %args = @_;
  my $complexity_limit = $args{'-low_complexity_limit'};

  my $complexity_fail_count = 0;

  #Go through each probe and check if has more than the maximum number of low complexity bases according to mdust
  foreach my $probe_id (keys %probes){
    my $low_complexity_bases = $probes{$probe_id}{mdust_bases};

    if ($low_complexity_bases >= $complexity_limit){
      $complexity_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($complexity_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);

  print BLUE, "\n\tA total of $complexity_fail_count ($percent_fail%) Probes failed the low complexity test ($complexity_limit or more low complexity bases)\n\n", RESET;
  print LOG "\n\tA total of $complexity_fail_count ($percent_fail%) Probes failed the low complexity test ($complexity_limit or more low complexity bases)\n\n";

  return();
}


#########################################################################################################
#2-c.) Pairfold - self-self folding test                                                                #
#########################################################################################################
sub testPairfoldScore{
  my %args = @_;
  my $pairfold_score_limit = $args{'-pairfold_score_limit'};

  my $pairfold_fail_count = 0;

  #Go through each probe and check if has more than the maximum number of low complexity bases according to mdust
  foreach my $probe_id (keys %probes){
    my $pairfold_score = $probes{$probe_id}{pairfold_score};

    if ($pairfold_score <= $pairfold_score_limit){
      $pairfold_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($pairfold_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);
  print BLUE, "\n\tA total of $pairfold_fail_count ($percent_fail%) Probes failed the pairfold test (PairFold score smaller than $pairfold_score_limit)\n\n", RESET;
  print LOG "\n\tA total of $pairfold_fail_count ($percent_fail%) Probes failed the pairfold test (PairFold score smaller than $pairfold_score_limit)\n\n";

  return();
}


#########################################################################################################
#2-d.) Simfold - internal folding test                                                                  #
#########################################################################################################
sub testSimfoldScore{
  my %args = @_;

  my $simfold_score_limit = $args{'-simfold_score_limit'};

  my $simfold_fail_count = 0;

  #Go through each probe and check if has more than the maximum number of low complexity bases according to mdust
  foreach my $probe_id (keys %probes){
    my $simfold_score = $probes{$probe_id}{simfold_score};

    if ($simfold_score <= $simfold_score_limit){
      $simfold_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($simfold_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);
  print BLUE, "\n\tA total of $simfold_fail_count ($percent_fail%) Probes failed the simfold test (SimFold score smaller than $simfold_score_limit)\n\n", RESET;
  print LOG "\n\tA total of $simfold_fail_count ($percent_fail%) Probes failed the simfold test (SimFold score smaller than $simfold_score_limit)\n\n";

  return();
}


#########################################################################################################
#2-e.) Specificity test - Do not allow any hits to mRNAs that map outside the target region             #
#########################################################################################################
sub testProbeSpecificity_mRNA{
  my %args = @_;
  my $mrna_hit_percent_limit = $args{'-non_target_hit_percent_limit'};

  my $mrna_specificity_fail_count = 0;

  #Go through each probe and check if has one or more significant non-target hits
  foreach my $probe_id (keys %probes){
    my $mrna_non_target_hit_length = $probes{$probe_id}{mrna_nt_hl};
    my $probe_length = $probes{$probe_id}{probe_length};

    if ((($mrna_non_target_hit_length/$probe_length)*100) >= $mrna_hit_percent_limit){
      $mrna_specificity_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($mrna_specificity_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);
  print BLUE, "\n\tA total of $mrna_specificity_fail_count ($percent_fail%) probes failed the mRNA specificity test (hit of >= $mrna_hit_percent_limit% length to mRNAs that map to a gene other than the target)\n\n", RESET;
  print LOG "\n\tA total of $mrna_specificity_fail_count ($percent_fail%) probes failed the mRNA specificity test (hit of >= $mrna_hit_percent_limit% length to mRNAs that map to a gene other than the target)\n\n";

  return();
}


#########################################################################################################
#2-f.) Specificity test - Do not allow any hits to ESTs that map outside the target region              #
#########################################################################################################
sub testProbeSpecificity_EST{
  my %args = @_;
  my $est_hit_percent_limit = $args{'-non_target_hit_percent_limit'};

  my $est_specificity_fail_count = 0;

  #Go through each probe and check if has one or more significant non-target hits
  foreach my $probe_id (keys %probes){
    my $est_non_target_hit_length = $probes{$probe_id}{est_nt_hl};
    my $probe_length = $probes{$probe_id}{probe_length};

    if ((($est_non_target_hit_length/$probe_length)*100) >= $est_hit_percent_limit){
      $est_specificity_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($est_specificity_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);
  print BLUE, "\n\tA total of $est_specificity_fail_count ($percent_fail%) probes failed the EST specificity test (hit of >= $est_hit_percent_limit% length to ESTs that map to a gene other than the target)\n\n", RESET;
  print LOG "\n\tA total of $est_specificity_fail_count ($percent_fail%) probes failed the EST specificity test (hit of >= $est_hit_percent_limit% length to ESTs that map to a gene other than the target)\n\n";

  return();
}


#######################################################################################################################
#2-g.) Specificity test - Do not all any hits to ensembl transcripts other than those associated with the target gene #
#######################################################################################################################
sub testProbeSpecificity_enst{
  my %args = @_;
  my $enst_hit_percent_limit = $args{'-non_target_hit_percent_limit'};

  my $enst_specificity_fail_count = 0;

  #Go through each probe and check if has one or more significant non-target hits
  foreach my $probe_id (keys %probes){
    my $enst_non_target_hit_length = $probes{$probe_id}{enst_nt_hl};
    my $probe_length = $probes{$probe_id}{probe_length};

    if ((($enst_non_target_hit_length/$probe_length)*100) >= $enst_hit_percent_limit){
      $enst_specificity_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($enst_specificity_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);
  print BLUE, "\n\tA total of $enst_specificity_fail_count ($percent_fail%) probes failed the ensembl specificity test (hit of >= $enst_hit_percent_limit% length to ensembl transcripts from a non-target gene)\n\n", RESET;
  print LOG "\n\tA total of $enst_specificity_fail_count ($percent_fail%) probes failed the ensembl specificity test (hit of >= $enst_hit_percent_limit% length to ensembl transcripts from a non-target gene)\n\n";

  return();
}


#######################################################################################################################
#2-h.) Specificity test - Do not allow any hits to other probes (ie. fail idential or near identical probes)          #
#######################################################################################################################
sub testProbeSpecificity_probe{
  my %args = @_;
  my $probe_hit_percent_limit = $args{'-non_target_hit_percent_limit'};

  my $probe_specificity_fail_count = 0;

  #Go through each probe and check if has one or more significant non-target hits
  foreach my $probe_id (keys %probes){
    my $probe_non_target_hit_length = $probes{$probe_id}{probe_nt_hl};
    my $probe_length = $probes{$probe_id}{probe_length};

    if ((($probe_non_target_hit_length/$probe_length)*100) >= $probe_hit_percent_limit){
      $probe_specificity_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($probe_specificity_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);
  print BLUE, "\n\tA total of $probe_specificity_fail_count ($percent_fail%) probes failed the probes specificity test (hit of >= $probe_hit_percent_limit% length to a probe from a different probeset)\n\n", RESET;
  print LOG "\n\tA total of $probe_specificity_fail_count ($percent_fail%) probes failed the probes specificity test (hit of >= $probe_hit_percent_limit% length to a probe from a different probeset)\n\n";

  return();
}


#######################################################################################################################
#2-i.) Specificity test - Do not allow any hits to the genome (only used for negative-control probes)                 #
#######################################################################################################################
sub testProbeSpecificity_genomic{
  my %args = @_;
  my $genomic_hit_percent_limit = $args{'-non_target_hit_percent_limit'};

  my $genomic_specificity_fail_count = 0;

  #Go through each probe and check if has one or more significant non-target hits
  foreach my $probe_id (keys %probes){
    my $genomic_non_target_hit_length = $probes{$probe_id}{genomic_nt_hl};
    my $probe_length = $probes{$probe_id}{probe_length};

    if ((($genomic_non_target_hit_length/$probe_length)*100) >= $genomic_hit_percent_limit){
      $genomic_specificity_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($genomic_specificity_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);
  print BLUE, "\n\tA total of $genomic_specificity_fail_count ($percent_fail%) probes failed the genomic specificity test (hit of >= $genomic_hit_percent_limit% length to the complete genome sequence)\n\n", RESET;
  print LOG "\n\tA total of $genomic_specificity_fail_count ($percent_fail%) probes failed the genomic specificity test (hit of >= $genomic_hit_percent_limit% length to the complete genome sequence)\n\n";

  return();
}


#######################################################################################################################
#2-j.) Specificity test - Do not allow any hits to a region of the targeted gene other than the part actually targeted#
#######################################################################################################################
sub testProbeSpecificity_withinGene{
  my %args = @_;
  my $withinGene_hit_percent_limit = $args{'-non_target_hit_percent_limit'};

  my $withinGene_specificity_fail_count = 0;

  #Go through each probe and check if has one or more significant non-target hits
  foreach my $probe_id (keys %probes){
    my $withinGene_non_target_hit_length = $probes{$probe_id}{withinGene_nt_hl};
    my $probe_length = $probes{$probe_id}{probe_length};

    if ((($withinGene_non_target_hit_length/$probe_length)*100) >= $withinGene_hit_percent_limit){
      $withinGene_specificity_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($withinGene_specificity_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);
  print BLUE, "\n\tA total of $withinGene_specificity_fail_count ($percent_fail%) probes failed the withinGene specificity test (hit of >= $withinGene_hit_percent_limit% length to a non-target region of the target gene)\n\n", RESET;
  print LOG "\n\tA total of $withinGene_specificity_fail_count ($percent_fail%) probes failed the withinGene specificity test (hit of >= $withinGene_hit_percent_limit% length to a non-target region of the target gene)\n\n";

  return();
}


##########################################################################################################################
#2-k.) Test probe uniquenes within probeset                                                                              #
##########################################################################################################################
sub testProbeUniquenessWithinProbeset{
  my $duplicate_probes = 0;

  #Go through each probeset and consider all the probe sequences
  foreach my $probeset_id (keys %probesets){
    my $probes_ref = $probesets{$probeset_id}{probes};

    my %probe_seqs;
    foreach my $probe_id (keys %{$probes_ref}){
      my $seq = $probes{$probe_id}{sequence};

      #If an identical probe to this one already exists for this probeset, fail them both!
      if ($probe_seqs{$seq}){
	my $previous_probe_id = $probe_seqs{$seq}{probe_id};
	$probes{$probe_id}{status} = "Fail";
	$probes{$previous_probe_id}{status} = "Fail";
	$duplicate_probes++;

	if ($verbose eq "yes"){
	  print YELLOW "\nFound a duplicate probe in probeset $probeset_id: ($probe_id and $previous_probe_id)", RESET;
	  print LOG "\nFound a duplicate probe in probeset $probeset_id: ($probe_id and $previous_probe_id)";
	}

      }else{
	$probe_seqs{$seq}{probe_id} = $probe_id;
      }
    }
  }

  print BLUE, "\n\tFound a duplicate probe $duplicate_probes times\n\n", RESET;
  print LOG "\n\tFound a duplicate probe $duplicate_probes times\n\n";

  return();
}


##########################################################################################################################
#Probeset size test - Use as a test to see how well the probe design is working                                          #
#  - Will not be used to pass or fail probes, simply used to give an idea how the filtering is working                   #
##########################################################################################################################
sub testProbesetSize{
  my %args = @_;
  my $probeset_size = $args{'-probeset_size'};

  my $full_probeset_count = 0;
  my $partial_probeset_count = 0;
  my $probeset_size_fail_count = 0;
  my $total_probesets = keys %probesets;

  foreach my $probeset_id (keys %probesets){

    #Note that the probeset size rule does not apply to Negative Control probes
    if ($probesets{$probeset_id}{probe_type} eq "Control-Negative"){
      next();
    }

    my $probes_ref = $probesets{$probeset_id}{probes};

    #Create a counter for the number of passed probes for each probesets (will be useful later)
    $probesets{$probeset_id}{passed_probe_count} = 0;

    #Count all the probes of this probeset that passed all the previous tests
    my $passed_probeset_probes = 0;
    foreach my $probe_id (keys %{$probes_ref}){
      if ($probes{$probe_id}{status} eq "Pass"){
	$passed_probeset_probes++;
      }
    }

    #Was the desired number of passing probes found for this probeset?
    if ($passed_probeset_probes >= $probeset_size){
      $probesets{$probeset_id}{passed_probe_count} = $passed_probeset_probes;
      $full_probeset_count++;
    }elsif($passed_probeset_probes >= 1){
      $probesets{$probeset_id}{passed_probe_count} = $passed_probeset_probes;
      $partial_probeset_count++;
    }else{
      $probeset_size_fail_count++;

      #If the probeset did not have enough passing probes, make sure all members of this probeset are marked as failed
      #By definition, they should already be!
      foreach my $probe_id (keys %{$probes_ref}){
	$probes{$probe_id}{status} = "Fail";
      }
    }
  }
  my $fail_rate = ($probeset_size_fail_count / $total_probesets) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);

  my $full_pass_rate = ($full_probeset_count / $total_probesets) * 100;
  my $percent_full_pass = sprintf("%.2f", $full_pass_rate);

  my $partial_pass_rate = ($partial_probeset_count / $total_probesets) * 100;
  my $percent_partial_pass = sprintf("%.2f", $partial_pass_rate);

  print BLUE, "\n\tA total of $full_probeset_count ($percent_full_pass%) have a full probeset ($probeset_size), while $partial_probeset_count ($percent_partial_pass%) have a partial probeset\n\n", RESET;
  print LOG "\n\tA total of $full_probeset_count ($percent_full_pass%) have a full probeset ($probeset_size), while $partial_probeset_count ($percent_partial_pass%) have a partial probeset\n\n";

  print BLUE, "\n\tA total of $probeset_size_fail_count ($percent_fail%) Probesets do not have any probes\n\n", RESET;
  print LOG "\n\tA total of $probeset_size_fail_count ($percent_fail%) Probesets do not have any probes\n\n";

  return();
}



###########################################################################################################################
#6.) Print out files containing the filtered probes                                                                       #
#Go back to the original input file and print each line out to a new file - if the probe on that line passed the tests    #
###########################################################################################################################
sub printOutputFile{
  my %args = @_;
  my $input_file = $args{'-input_file'};
  my $output_file = $args{'-output_file'};

  print BLUE, "\nPrinting filtered probes to output file: $output_file\n\n", RESET;
  print LOG "\nPrinting filtered probes to output file: $output_file\n\n";

  my $passed_probe_count = 0;

  #Open the input and output files
  open(INFILE, "$input_file") || die "\nCould not open input file: $input_file\n\n";
  open(OUTFILE, ">$output_file") || die "\nCould not open output file: $output_file\n\n";

  my $header = 1;
  while(<INFILE>){
    my $probe_line = $_;
    #Print out the header line
    if ($header == 1){
      print OUTFILE "$_";
      $header = 0;
      next();
    }

    my @line = split ("\t", $_);

    if ($line[0] =~ /^(\d+)/){
      my $probe_id = $1;

      if ($probes{$probe_id}{status} eq "Pass"){
	$passed_probe_count++;

	#If this probe_id corresponds to a probe that passed, print it out
	print OUTFILE "$_";
      }

    }else{
      print RED, "\nDid not find a valid probe ID on the line: $_", RESET;
      exit();
    }
  }

  print BLUE, "\nPrinted a total of $passed_probe_count probes to the filtered probe file\n\n", RESET;
  print LOG "\nPrinted a total of $passed_probe_count probes to the filtered probe file\n\n";

  close (INFILE);
  close (OUTFILE);
  return();
}

