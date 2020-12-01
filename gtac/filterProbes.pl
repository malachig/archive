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
#2-b.) Probes must have less than x% of their length comprised of RepeatMasked bases
#      - Probes are already pre-filtered during extraction to eliminate those with more than 25% masked bases.
#      - Use this filter to make this more stringent if desired
#2-c.) Probes must have self-self folding free energy scores of less than the user specified threshhold '--pairfold_score_limit' - Determined by Pairfold
#2-d.) Probes must have internal folding free energy scores of less than the use specified threshhold '--simfold_score_limit' - Determined by SimFold
#2-e.) Probes must not contain low complexity regions larger than '--complexity_region_limit' (mono-, di-, or tri-nucleotide repeats) - determined by mdust

#2-f.) Probes must be specific to the targeted genomic regions - determined by blast against the complete genome
#      - Fail probes that have more than one hit overlapping the targeted region (indicates repeats or low complexity regions)

#2-g.) Fail probes that have more than Y hits to non-target regions that are larger than 50% of the probe length?
#2-h.) Fail probes that have more than Z hits to non-target regions that are larger than 75% of the probe length?
#2-i.) Fail probes that have a hit to a region other than that targeted if they are larger than x% of the probe length
#2-j.) Fail probes that require more than N cycles to synthesize.  Set N to max allowed by NimbleGen (N = 180?)

#3.) Summarize the number of successful probes per probeset

#4.) Print a probe file containing only the filtered probes
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
my $masked_percent_limit = '';
my $simfold_score_limit = '';
my $pairfold_score_limit = '';
my $complexity_region_percent_limit = '';
my $genomic_50percent_hits_limit = '';
my $genomic_75percent_hits_limit = '';
my $genomic_hit_percent_limit = '';
my $linker_hit_percent_limit = '';
my $nimblegen_cycles_limit = '';
my $probeset_size = '';
my $outfile = '';
my $logfile = '';
my $verbose = '';
my $region_list = '';

GetOptions ('infile=s'=>\$infile, 'target_length=i'=>\$target_length,
	    'target_tm=f'=>\$target_tm, 'tm_range=f'=>\$tm_range, 'masked_percent_limit=f'=>\$masked_percent_limit, 
	    'simfold_score_limit=f'=>\$simfold_score_limit, 'pairfold_score_limit=f'=>\$pairfold_score_limit,
	    'complexity_region_percent_limit=f'=>\$complexity_region_percent_limit,
	    'genomic_50percent_hits_limit=f'=>\$genomic_50percent_hits_limit, 'genomic_75percent_hits_limit=f'=>\$genomic_75percent_hits_limit,
	    'genomic_hit_percent_limit=f'=>\$genomic_hit_percent_limit, 'linker_hit_percent_limit=f'=>\$linker_hit_percent_limit,
	    'nimblegen_cycles_limit=i'=>\$nimblegen_cycles_limit,
	    'probeset_size=i'=>\$probeset_size, 'outfile=s'=>\$outfile, 'logfile=s'=>\$logfile, 'verbose=s'=>\$verbose,
	    'region_list=s'=>\$region_list);

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
print GREEN, "\n\tThe percent repeatmasked bases cutoff using: --masked_percent_limit (say 15.0%)", RESET;
print GREEN, "\n\tThe cutoff SimFold score (within-probe/internal folding) that will cause a probe to fail using: --simfold_score_limit (say -10.0)", RESET;
print GREEN, "\n\tThe cutoff PairFold score (self-self/between probe folding) that will cause a probe to fail using: --pairfold_score_limit (say -22.0)", RESET;
print GREEN, "\n\tThe maximum allowed percent of low complexity bases allowed using: --complexity_region_percent_limit (say 10.0%)", RESET;
print GREEN, "\n\tThe maximum number of non-specific genomic hits of 50% of probe length or greater using: --genomic_50percent_hits_limit (say 5)", RESET;
print GREEN, "\n\tThe maximum number of non-specific genomic hits of 75% of probe length or greater using: --genomic_75percent_hits_limit (say 2)", RESET;
print GREEN, "\n\tThe Genomic specificity hit percent length cutoff using: --genomic_hit_percent_limit (say 80.0%)", RESET;
print GREEN, "\n\tThe maximum hit percent length to Solexa linkers allowed using --linker_hit_percent_limit (say 0% to reject any hit)", RESET
print GREEN, "\n\tThe maximum number of cycles required for synthesis using: --nimblegen_cycles_limit (say 180)", RESET;
print GREEN, "\n\tThe desired probeset size (number of probes per target region) using: --probeset_size (say 8)", RESET;
print GREEN, "\n\tFor verbose output, use: verbose=yes\n\n", RESET;
print GREEN, "\n\tTo generate a filtered set for only a subset of target region, supply a region list with: --region_list (OPTIONAL)", RESET;
print GREEN, "\n\nExample filterProbes.pl  --infile=RegionProbes.txt  --target_length=54  --target_tm=76.0  --tm_range=5  --masked_percent_limit=15.0  --simfold_score_limit=-10.0  --pairfold_score_limit=-22.0  --complexity_region_percent_limit=10.0  --genomic_50percent_hits_limit=5  --genomic_75percent_hits_limit=2  --genomic_hit_percent_limit=80.0  --linker_hit_percent_limit=0  --nimblegen_cycles_limit=180  --probeset_size=8  --outfile=RegionProbes_filtered.txt  --logfile=filter_RegionProbes_LOG.txt  --verbose=no\n\n", RESET;

unless ($infile && $target_length && $target_tm && $tm_range && $masked_percent_limit && $simfold_score_limit && $pairfold_score_limit && $complexity_region_percent_limit && $genomic_50percent_hits_limit && $genomic_75percent_hits_limit && $genomic_hit_percent_limit && $linker_hit_percent_limit && $nimblegen_cycles_limit && $probeset_size && $outfile && $logfile && $verbose){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Open logfile
open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile for output\n\n";

#Print out the parameters selected by the user so that this information can be explicitly stored in a LOG file
print BLUE, "\n\nPARAMETERS USED:\ninfile = $infile\ntarget_length = $target_length\ntarget_tm = $target_tm\ntm_range = $tm_range\nmasked_percent_limit = $masked_percent_limit\nsimfold_score_limit = $simfold_score_limit\npairfold_score_limit = $pairfold_score_limit\ncomplexity_region_percent_limit = $complexity_region_percent_limit\ngenomic_50percent_hits_limit = $genomic_50percent_hits_limit\ngenomic_75percent_hits_limit = $genomic_75percent_hits_limit\ngenomic_hit_percent_limit = $genomic_hit_percent_limit\nlinkers_hit_percent_limit = $linker_hit_percent_limit\nnimblegen_cycles_limit = $nimblegen_cycles_limit\nprobeset_size = $probeset_size\noutfile = $outfile\nlogile = $logfile\nverbose = $verbose\n\n", RESET;

print LOG "\n\nPARAMETERS USED:\ninfile = $infile\ntarget_length = $target_length\ntarget_tm = $target_tm\ntm_range = $tm_range\nmasked_percent_limit = $masked_percent_limit\nsimfold_score_limit = $simfold_score_limit\npairfold_score_limit = $pairfold_score_limit\ncomplexity_region_percent_limit = $complexity_region_percent_limit\ngenomic_50percent_hits_limit = $genomic_50percent_hits_limit\ngenomic_75percent_hits_limit = $genomic_75percent_hits_limit\ngenomic_hit_percent_limit = $genomic_hit_percent_limit\nlinkers_hit_percent_limit = $linker_hit_percent_limit\nnimblegen_cycles_limit = $nimblegen_cycles_limit\nprobeset_size = $probeset_size\noutfile = $outfile\nlogile = $logfile\nverbose = $verbose\n\n";

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

##2-a.) Tm test
#      - Probes must have a Tm equal to '--target_tm' +/- '--tm_range'  (for example, 67C +/- 3C) - determined by nearest neighbour Tm calculation
if ($columns{Probe_Tm}){
  &testProbeTms('-target_tm'=>$target_tm, '-tm_range'=>$tm_range);
}else{
  print YELLOW, "\n\tProbe_Tm values not found in input file, skipping the Tm test!\n", RESET;
  print LOG "\n\tProbe_Tm values not found in input file, skipping the Tm test!\n";
}

##2-b.) Repeatmask test
#      - Probes must have less than x% of their length comprised of RepeatMasked bases
#      - Probes are already pre-filtered during extraction to eliminate those with more than 25% masked bases.
#      - Use this filter to make this more stringent if desired
if ($columns{masked_bases}){
  &testMaskedBases('-masked_percent_limit'=>$masked_percent_limit);
}else{
  print YELLOW, "\n\tmasked_bases values not found in input file, skipping the masked bases test!\n", RESET;
  print LOG "\n\tmasked_bases values not found in input file, skipping the masked bases test!\n";
}

##2-c) Pairfold - self-self folding test
#     - Probes must have self-self folding free energy scores of less than the user specified threshhold '--pairfold_score_limit' - Determined by Pairfold
if ($columns{PairFold_score}){
  &testPairfoldScore('-pairfold_score_limit'=>$pairfold_score_limit);
}else{
  print YELLOW, "\n\tPairFold_score values not found in input file, skipping the pairfold score test!\n", RESET;
  print LOG "\n\tPairFold_score values not found in input file, skipping the pairfold score test!\n";
}

##2-d) Simfold - internal folding test
#     - Probes must have internal folding free energy scores of less than the use specified threshhold '--simfold_score_limit' - Determined by SimFold
if ($columns{SimFold_score}){
  &testSimfoldScore('-simfold_score_limit'=>$simfold_score_limit);
}else{
  print YELLOW, "\n\tSimFold_score values not found in input file, skipping the simfold score test!\n", RESET;
  print LOG "\n\tSimFold_score values not found in input file, skipping the simfold score test!\n";
}

##2.e.) Low complexity test
#      - Probes must not contain low complexity regions larger than '--complexity_region_limit' (mono-, di-, or tri-nucleotide repeats) - determined by mdust
if ($columns{MdustBases_cutoff_11}){
  &testProbeComplexity('-low_complexity_percent_limit'=>$complexity_region_percent_limit);
}else{
  print YELLOW, "\n\tMdustBases_cutoff_11 values not found in input file, skipping the low complexity test!\n", RESET;
  print LOG "\n\tMdustBases_cutoff_11 values not found in input file, skipping the low complexity test!\n";
}

##2-f.) Target region specificity test
#      - Probes must be specific to the targeted genomic regions - determined by blast against the complete genome
#      - Fail probes that have more than one hit overlapping the targeted region (indicates repeats or low complexity regions)
if ($columns{'genomic_TargetHits'}){
  &testProbeSpecificity_TargetRegion();
}else{
  print YELLOW, "\n\tgenomic_TargetHits values not found in input file, skipping the genomic Target Region specificity test!\n", RESET;
  print LOG "\n\tgenomic_TargetHits values not found in input file, skipping the genomic Target Region specificity test!\n";
}

##2-g.) 50% non-target hits test
#      - Fail probes that have more than Y hits to non-target regions that are larger than 50% of the probe length?
if ($columns{'genomic_Non-TargetHits_over50'}){
  &testProbeSpecificity_genomic_over50('-genomic_50percent_hits_limit'=>$genomic_50percent_hits_limit);
}else{
  print YELLOW, "\n\tgenomic_Non-TargetHits_over50 values not found in input file, skipping the genomic over50 specificity test!\n", RESET;
  print LOG "\n\tgenomic_Non-TargetHits_over50 values not found in input file, skipping the genomic over50 specificity test!\n";
}


##2-h.) 75% non-target hits test
#      - Fail probes that have more than Z hits to non-target regions that are larger than 75% of the probe length?
if ($columns{'genomic_Non-TargetHits_over75'}){
  &testProbeSpecificity_genomic_over75('-genomic_75percent_hits_limit'=>$genomic_75percent_hits_limit);
}else{
  print YELLOW, "\n\tgenomic_Non-TargetHits_over75 values not found in input file, skipping the genomic over75 specificity test!\n", RESET;
  print LOG "\n\tgenomic_Non-TargetHits_over75 values not found in input file, skipping the genomic over75 specificity test!\n";
}

##2-i.) Genomic Specificity test
#      - Fail probes that have a hit to a region other than that targeted if they are larger than x% of the probe length
if ($columns{'genomic_largestNon-TargetHitLength'}){
  &testProbeSpecificity_genomic('-non_target_hit_percent_limit'=>$genomic_hit_percent_limit);
}else{
  print YELLOW, "\n\tgenomic_largestNon-TargetHitLength values not found in input file, skipping the genomic specificity test!\n", RESET;
  print LOG "\n\tgenomic_largestNon-TargetHitLength values not found in input file, skipping the genomic specificity test!\n";
}

##2-j.) Solexa linkers specificity test - (linkers_hit_percent_limit)
if ($columns{'linker_largestNon-TargetHitLength'}){
  &testProbeSpecificity_linker('-linker_hit_percent_limit'=>$linker_hit_percent_limit);
}else{
  print YELLOW, "\n\tlinker_largestNon-TargetHitLength values not found in input file, skipping the linker specificity test!\n", RESET;
  print LOG "\n\tlinker_largestNon-TargetHitLength values not found in input file, skipping the linker specificity test!\n";
}


##2-k.) NimbleGen cycles test
#      - Fail probes that require more than N cycles to synthesize.  Set N to max allowed by NimbleGen (N = 180?)
if ($columns{'NimbleGenCycles'}){
  &testProbeCycles('-nimblegen_cycles_limit'=>$nimblegen_cycles_limit);
}else{
  print YELLOW, "\n\tNimbleGenCycles values not found in input file, skipping the genomic specificity test!\n", RESET;
  print LOG "\n\tNimbleGenCycles values not found in input file, skipping the genomic specificity test!\n";
}



#3.) Summarize the number of successful probes per probeset
#    - Probeset size test - Probesets with less than the desired minimum number of probes will be noted.
#    - Will not be used to pass or fail probes, simply used to give an idea how the filtering is working
&testProbesetSize('-probeset_size'=>$probeset_size);


#4.) Print out files containing the filtered probes
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
  #Probe_Count, ProbeSet_ID, Region_ID, Probe_length, Probe_Tm, masked_bases, SimFold_score, PairFold_score, MdustBases_cutoff_11,
  #genomic_TargetHits, genomic_Non-TargetHits, genomic_Non-TargetHits_over50, genomic_Non-TargetHits_over75
  #genomic_largestNon-TargetHitLength, CycleCount

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
      unless ($columns{Probe_Count} && $columns{ProbeSet_ID}){
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

    my $region_id;
    if ($columns{Region_ID}){
      $region_id = $probe_line[$columns{Region_ID}{column_pos}];
    }else{
      $region_id = "na";
    }
    $probes{$probe_id}{region} = $region_id;

    unless ($probe_id =~ /^\d+/ && $probeset_id =~ /^\d+/){
      print RED, "\nInvalid probe or probeset ID\n\nLINE: $_\n\n", RESET;
      exit();
    }
    $probe_count++;

    #Store probe-to-probeset mapping
    if ($probesets{$probeset_id}){
      my $probe_list_ref = $probesets{$probeset_id}{probes};
      $probe_list_ref->{$probe_id}->{t} = '';
    }else{
      my %probe_list;
      $probe_list{$probe_id}{t} = '';
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
    if ($columns{Probe_Tm}){
      my $tm_celsius = $probe_line[$columns{Probe_Tm}{column_pos}];
      if ($tm_celsius =~ /\d+\.\d+/){
	$probes{$probe_id}{tm_celsius} = $tm_celsius;
      }else{
	print RED, "\nTm value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{masked_bases}){
      my $masked_bases = $probe_line[$columns{masked_bases}{column_pos}];
      if ($masked_bases =~ /^\d+/){
	$probes{$probe_id}{masked_bases} = $masked_bases;
      }elsif($masked_bases eq "na"){
	$probes{$probe_id}{masked_bases} = 0;
      }else{
	print RED, "\nMasked bases value missing - or incorrect format\tProbe = $probe_id", RESET;
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
    if ($columns{'genomic_TargetHits'}){
      my $genomic_target_hits = $probe_line[$columns{'genomic_TargetHits'}{column_pos}];
      if ($genomic_target_hits =~ /\d+/){
	$probes{$probe_id}{genomic_target_hits} = $genomic_target_hits;
      }else{
	print RED, "\ngenomic Target Hits value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{'genomic_Non-TargetHits'}){
      my $genomic_non_target_hits = $probe_line[$columns{'genomic_Non-TargetHits'}{column_pos}];
      if ($genomic_non_target_hits =~ /\d+/){
	$probes{$probe_id}{genomic_non_target_hits} = $genomic_non_target_hits;
      }else{
	print RED, "\nGenomic Non-Target Hits value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{'genomic_Non-TargetHits_over50'}){
      my $genomic_non_target_hits_over50 = $probe_line[$columns{'genomic_Non-TargetHits_over50'}{column_pos}];
      if ($genomic_non_target_hits_over50 =~ /\d+/){
	$probes{$probe_id}{genomic_non_target_hits_50} = $genomic_non_target_hits_over50;
      }else{
	print RED, "\nGenomic Non-Target Hits over50 value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{'genomic_Non-TargetHits_over75'}){
      my $genomic_non_target_hits_over75 = $probe_line[$columns{'genomic_Non-TargetHits_over75'}{column_pos}];
      if ($genomic_non_target_hits_over75 =~ /\d+/){
	$probes{$probe_id}{genomic_non_target_hits_75} = $genomic_non_target_hits_over75;
      }else{
	print RED, "\nGenomic Non-Target Hits over75 value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{'genomic_largestNon-TargetHitLength'}){
      my $genomic_non_target_hit_length = $probe_line[$columns{'genomic_largestNon-TargetHitLength'}{column_pos}];
      if ($genomic_non_target_hit_length =~ /\d+/){
	$probes{$probe_id}{genomic_non_target_hit_length} = $genomic_non_target_hit_length;
      }else{
	print RED, "\nGenomic Non-Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{'linker_largestNon-TargetHitLength'}){
      my $linker_hit_length = $probe_line[$columns{'linker_largestNon-TargetHitLength'}{column_pos}];
      if ($linker_hit_length =~ /\d+/){
	$probes{$probe_id}{linker_non_target_hit_length} = $linker_hit_length;
      }else{
	print RED, "\nLinker Non-target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	exit();
      }
    }
    if ($columns{'NimbleGenCycles'}){
      my $cycle_count = $probe_line[$columns{'NimbleGenCycles'}{column_pos}];
      if ($cycle_count =~ /\d+/){
	$probes{$probe_id}{cycle_count} = $cycle_count;
      }else{
	print RED, "\nCycle count value missing - or incorrect format\tProbe = $probe_id", RESET;
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

  print YELLOW, "\n\tA total of $tm_fail_count ($percent_fail%) Probes failed the Tm test (Tm outside the range: $lower_limit - $upper_limit)\n\n", RESET;
  print LOG "\n\tA total of $tm_fail_count ($percent_fail%) Probes failed the Tm test (Tm outside the range: $lower_limit - $upper_limit)\n\n";

  return();
}



#########################################################################################################
#2-b.) Repeatmask test                                                                                  #
#########################################################################################################
sub testMaskedBases{
  my %args = @_;
  my $masked_percent_limit = $args{'-masked_percent_limit'};

  my $masked_bases_fail_count = 0;

  #Go through each probe and check if it has too many masked bases
  foreach my $probe_id (keys %probes){
    my $masked_bases = $probes{$probe_id}{masked_bases};
    my $probe_length = $probes{$probe_id}{probe_length};

    if ((($masked_bases/$probe_length)*100) > $masked_percent_limit){
      $masked_bases_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($masked_bases_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);
  print YELLOW, "\n\tA total of $masked_bases_fail_count ($percent_fail%) probes failed the masked bases test (> $masked_percent_limit% of probe is masked bases)\n\n", RESET;
  print LOG "\n\tA total of $masked_bases_fail_count ($percent_fail%) probes failed the masked bases test (> $masked_percent_limit% of probe is masked bases)\n\n";

  return();
}


#########################################################################################################
#2-c.) Pairfold - self-self folding test                                                                #
#########################################################################################################
sub testPairfoldScore{
  my %args = @_;
  my $pairfold_score_limit = $args{'-pairfold_score_limit'};

  my $pairfold_fail_count = 0;

  #Go through each probe and check if exceeds the allowed pairfold score
  foreach my $probe_id (keys %probes){
    my $pairfold_score = $probes{$probe_id}{pairfold_score};

    if ($pairfold_score <= $pairfold_score_limit){
      $pairfold_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($pairfold_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);
  print YELLOW, "\n\tA total of $pairfold_fail_count ($percent_fail%) Probes failed the pairfold test (PairFold score smaller than $pairfold_score_limit)\n\n", RESET;
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

  #Go through each probe and check if exceeds the allowed pairfold score
  foreach my $probe_id (keys %probes){
    my $simfold_score = $probes{$probe_id}{simfold_score};

    if ($simfold_score <= $simfold_score_limit){
      $simfold_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($simfold_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);
  print YELLOW, "\n\tA total of $simfold_fail_count ($percent_fail%) Probes failed the simfold test (SimFold score smaller than $simfold_score_limit)\n\n", RESET;
  print LOG "\n\tA total of $simfold_fail_count ($percent_fail%) Probes failed the simfold test (SimFold score smaller than $simfold_score_limit)\n\n";

  return();
}


#########################################################################################################
#2-e.) Low Complexity Test                                                                              #
#########################################################################################################
sub testProbeComplexity{
  my %args = @_;
  my $complexity_percent_limit = $args{'-low_complexity_percent_limit'};

  my $complexity_fail_count = 0;

  #Go through each probe and check if has more than the maximum number of low complexity bases according to mdust
  foreach my $probe_id (keys %probes){
    my $low_complexity_bases = $probes{$probe_id}{mdust_bases};
    my $probe_length = $probes{$probe_id}{probe_length};


    if ((($low_complexity_bases/$probe_length)*100) >= $complexity_percent_limit){
      $complexity_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($complexity_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);

  print YELLOW, "\n\tA total of $complexity_fail_count ($percent_fail%) Probes failed the low complexity test ($complexity_percent_limit % or more low complexity bases)\n\n", RESET;
  print LOG "\n\tA total of $complexity_fail_count ($percent_fail%) Probes failed the low complexity test ($complexity_percent_limit % or more low complexity bases)\n\n";

  return();
}


#########################################################################################################
#2-f.) Target region specificity test                                                                   #
#########################################################################################################
sub testProbeSpecificity_TargetRegion{

  my $target_region_fail_count = 0;
  my $single_hit_count = 0;
  my $no_hit_count = 0;

  #Go through each probe and check if has more than 1 hit overlapping the target region
  foreach my $probe_id (keys %probes){

    my $genomic_target_hits = $probes{$probe_id}{genomic_target_hits};

    if ($genomic_target_hits > 1){
      $target_region_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }elsif($genomic_target_hits == 1){
      $single_hit_count++;
    }else{
      unless($probes{$probe_id}{region} eq "na"){	
        print RED, "\nProbe: $probe_id did not hit its own target genomic region!", RESET;	
        $no_hit_count++;
      }
    }

  }
  my $fail_rate = ($target_region_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);

  my $percent_single = sprintf("%.2f", (($single_hit_count/$probe_count)* 100));
  my $percent_none = sprintf("%.2f", (($no_hit_count/$probe_count)* 100));

  print YELLOW, "\n\tA total of $target_region_fail_count ($percent_fail%) Probes have more than 1 overlaping hits to the target region", RESET;
  print YELLOW, "\n\tA total of $single_hit_count ($percent_single%) had a single hit and $no_hit_count ($percent_none%) had 0 hits to target region\n\n", RESET;
  print LOG "\n\tA total of $target_region_fail_count ($percent_fail%) Probes have more than 1 overlaping hits to the target region";
  print LOG "\n\tA total of $single_hit_count ($percent_single%) had a single hit and $no_hit_count ($percent_none%) had 0 hits to target region\n\n";

  return();
}


#########################################################################################################
#2-g.) 50% non-target hits test
#########################################################################################################
sub testProbeSpecificity_genomic_over50{
  my %args = @_;
  my $genomic_50percent_hits_limit = $args{'-genomic_50percent_hits_limit'};
  my $genomic50_fail_count = 0;

  #Go through each probe and check if has more than the maximum number of low complexity bases according to mdust
  foreach my $probe_id (keys %probes){
    my $genomic_non_target_hits_over50 = $probes{$probe_id}{genomic_non_target_hits_50};

    if ($genomic_non_target_hits_over50 > $genomic_50percent_hits_limit){
      $genomic50_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($genomic50_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);
  print YELLOW, "\n\tA total of $genomic50_fail_count ($percent_fail%) Probes failed the genomic hits over 50% test (more than $genomic_50percent_hits_limit non-specific hits over 50% of probe length)\n\n", RESET;
  print LOG "\n\tA total of $genomic50_fail_count ($percent_fail%) Probes failed the genomic hits over 50% test (more than $genomic_50percent_hits_limit non-specific hits over 50% of probe length)\n\n";

  return();
}


#########################################################################################################
#2-h.) 75% non-target hits test
#########################################################################################################
sub testProbeSpecificity_genomic_over75{
  my %args = @_;
  my $genomic_75percent_hits_limit = $args{'-genomic_75percent_hits_limit'};
  my $genomic75_fail_count = 0;

  #Go through each probe and check if has more than the maximum number of low complexity bases according to mdust
  foreach my $probe_id (keys %probes){
    my $genomic_non_target_hits_over75 = $probes{$probe_id}{genomic_non_target_hits_75};

    if ($genomic_non_target_hits_over75 > $genomic_75percent_hits_limit){
      $genomic75_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($genomic75_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);
  print YELLOW, "\n\tA total of $genomic75_fail_count ($percent_fail%) Probes failed the genomic hits over 75% test (more than $genomic_75percent_hits_limit non-specific hits over 75% of probe length)\n\n", RESET;
  print LOG "\n\tA total of $genomic75_fail_count ($percent_fail%) Probes failed the genomic hits over 75% test (more than $genomic_75percent_hits_limit non-specific hits over 75% of probe length)\n\n";

  return();
}


#########################################################################################################
#2-i.) Genomic Specificity test                                                                         #
#########################################################################################################
sub testProbeSpecificity_genomic{
  my %args = @_;
  my $genomic_hit_percent_limit = $args{'-non_target_hit_percent_limit'};

  my $genomic_specificity_fail_count = 0;

  #Go through each probe and check if has one or more significant non-target hits
  foreach my $probe_id (keys %probes){
    my $genomic_non_target_hit_length = $probes{$probe_id}{genomic_non_target_hit_length};
    my $probe_length = $probes{$probe_id}{probe_length};

    if ((($genomic_non_target_hit_length/$probe_length)*100) > $genomic_hit_percent_limit){
      $genomic_specificity_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($genomic_specificity_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);
  print YELLOW, "\n\tA total of $genomic_specificity_fail_count ($percent_fail%) probes failed the genomic specificity test (hit of > $genomic_hit_percent_limit% length to the complete genome sequence)\n\n", RESET;
  print LOG "\n\tA total of $genomic_specificity_fail_count ($percent_fail%) probes failed the genomic specificity test (hit of > $genomic_hit_percent_limit% length to the complete genome sequence)\n\n";

  return();
}


#########################################################################################################
#2-j.) Solexa linkers specificity test - (linkers_hit_percent_limit)                                    #
#########################################################################################################
sub testProbeSpecificity_linker{
  my %args = @_;
  my $linker_hit_percent_limit = $args{'-linker_hit_percent_limit'};

  my $linker_specificity_fail_count = 0;

  #Go through each probe and check if has one or more significant non-target hits
  foreach my $probe_id (keys %probes){
    my $linker_non_target_hit_length = $probes{$probe_id}{linker_non_target_hit_length};
    my $probe_length = $probes{$probe_id}{probe_length};

    if ((($linker_non_target_hit_length/$probe_length)*100) > $linker_hit_percent_limit){
      $linker_specificity_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($linker_specificity_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);
  print YELLOW, "\n\tA total of $linker_specificity_fail_count ($percent_fail%) probes failed the linker specificity test (hit of > $linker_hit_percent_limit% length to one or more Solexa linkers)\n\n", RESET;
  print LOG "\n\tA total of $linker_specificity_fail_count ($percent_fail%) probes failed the linker specificity test (hit of > $linker_hit_percent_limit% length to one or more Solexa linkers)\n\n";

  return();
}



#########################################################################################################
#2-k.) NimbleGen cycles test                                                                            #
#########################################################################################################
sub testProbeCycles{
  my %args = @_;
  my $nimblegen_cycles_limit = $args{'-nimblegen_cycles_limit'};

  my $cycles_fail_count = 0;

  #Go through each probe and check if has more than the maximum number of low complexity bases according to mdust
  foreach my $probe_id (keys %probes){
    my $cycle_count = $probes{$probe_id}{cycle_count};

    if ($cycle_count >= $nimblegen_cycles_limit){
      $cycles_fail_count++;
      $probes{$probe_id}{status} = "Fail";
    }
  }
  my $fail_rate = ($cycles_fail_count / $probe_count) * 100;
  my $percent_fail = sprintf("%.2f", $fail_rate);

  print YELLOW, "\n\tA total of $cycles_fail_count ($percent_fail%) Probes failed the NimbleGen cycles test (>= $nimblegen_cycles_limit required)\n\n", RESET;
  print LOG "\n\tA total of $cycles_fail_count ($percent_fail%) Probes failed the NimbleGen cycles test (>= $nimblegen_cycles_limit required)\n\n";

  return();
}


##########################################################################################################################
#3.) Probeset size test - Use as a test to see how well the probe design is working                                      #
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

  print YELLOW, "\n\tA total of $full_probeset_count ($percent_full_pass%) have a full probeset ($probeset_size), while $partial_probeset_count ($percent_partial_pass%) have a partial probeset", RESET;
  print LOG "\n\tA total of $full_probeset_count ($percent_full_pass%) have a full probeset ($probeset_size), while $partial_probeset_count ($percent_partial_pass%) have a partial probeset";

  print YELLOW, "\n\tA total of $probeset_size_fail_count ($percent_fail%) Probesets do not have any probes\n\n", RESET;
  print LOG "\n\tA total of $probeset_size_fail_count ($percent_fail%) Probesets do not have any probes\n\n";

  return();
}



###########################################################################################################################
#4.) Print out files containing the filtered probes                                                                       #
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

      my $region_id;
      if ($columns{Region_ID}){
	$region_id = $line[$columns{Region_ID}{column_pos}];
      }else{
	$region_id = "na";
      }

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

