#!/usr/bin/perl -w

#Written by Malachi Griffith
#The purpose of this script is to select probe sequences within ensembl introns for a specified gene or list of genes
#The user will specify the desired target probe length (--target_length) and the degree of overlap/tiling (--overlap) (i.e. how closely to design probes to each other)
#If the value for --overlap is less than the --target_length, the designed probes will overlap each other
#At every design position, the length of the probe will be varied up to the --max_length_variance specified by the user to achieve the specified --target_tm
#If you want anisothermal probes of fixed length, simply set --max_variance_length=0
#The desired melting temperature is the median Tm of all exon-exon and intron-exon junction probes

#In this way, many probes will be designed for each intron of each gene.  From this pool, a selection of representative probes will be chosen
#In general I want at least a few probes per intron and they should be near to the centre of the intron and non-overlaping if possible.
#Introns will be defined relative to all of the exon content for the gene (define exon content for the gene and the intron content is everything left over)

#Note that in general we can be a bit more stringent in the selection of intron probes because unlike for junctions there is a larger region to select from

#In general the neccessary code can be summarized as follows:

#(1.) Get the single gene ID or list of IDs from the user
#     - I will design intron probes only for a control set of genes which correspond to the housekeeping control genes used on Affy's arrays

#(2.) For each gene, get all exons and their gene coordinates via ALEXA
#     Determine the exon content of the gene and thereby define the constitutively intronic regions of the gene. Then for each intron:
#     - Extract probe sequences of length '--target_length' +/- '--max_length_variance', starting every '--overlap' bases as specified by the user
#     - Make sure all probe sequences are within the boundaries of the intron (also check for introns that are too short)
#     - Store the probe sequence, start position, and end position as a hash
#     - Eliminate probes that correspond to Genomic N's
#     - Eliminate probes that are comprised of too many repeat masked bases
#     - Eliminate probes that are outside the desired target Tm +/- the acceptable range
#     - Eliminate probes with simple repeats such as monopolymers and dipolymers

#-The end result of this script will be a file containing many probes for each intron region of each gene.  A seperate script will then evaluate these:
#   - Essentially we want the 'best' probe(s) for each intron region
#   - We want the probes with the highest specificity (ie. unique in the transcriptome if possible -avoid common sequences from gene families
#   - Also the 'best' probes are those that have a Tm that is as close as possible to the target Tm and the probe is close to the center of the intron
#   - If there are multiple probes within x degrees of the Tm (where x is some reasonable range, say 5 degrees C) then chose the one closest to the centre
#   - Or alternatively come up with a scoring system to determine the probe that is close in Tm and close to the centre
#   - Store the Tm, distance from centre, coordinates, etc. and summarize for all chosen probe

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
use utilities::Probes qw(:all);

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $overlap = '';
my $target_tm = '';
my $tm_range = '';
my $target_length = '';
my $max_length_variance = '';
my $ensembl_gene_id = '';
my $gene_file = '';
my $all_genes = '';
my $allow_predicted_genes = '';
my $verbose = '';
my $probe_dir = '';
my $outfile = '';
my $logfile = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'overlap=i'=>\$overlap, 'target_tm=f'=>\$target_tm, 'tm_range=f'=>\$tm_range, 'target_length=i'=>\$target_length,
	    'max_length_variance=i'=>\$max_length_variance, 'ensembl_gene_id=s'=>\$ensembl_gene_id, 'gene_file=s'=>\$gene_file,
	    'all_genes=s'=>\$all_genes, 'allow_predicted_genes=s'=>\$allow_predicted_genes, 'verbose=s'=>\$verbose,
	    'probe_dir=s'=>\$probe_dir, 'outfile=s'=>\$outfile, 'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the tiling distance or overlap between adjacent probe selection regions using: --overlap (eg. 5)", RESET;
print GREEN, "\n\tSpecify the target Tm for probe sequences using: --target_tm (eg. 67.0)", RESET;
print GREEN, "\n\tSpecify the Tm range allowed in case varying the length does not achieve the target Tm using: --tm_range (eg. 2.0)", RESET;
print GREEN, "\n\tSpecify the target probe length (must be a multiple of 2) using: --target_length (eg. 36)", RESET;
print GREEN, "\n\tSpecify the maximum this probe length will be allowed to vary using: --max_length_variance (eg. 5)", RESET;
print GREEN, "\n\t\tIf you want to test with a single gene ID use: --ensembl_gene_id (ENSG00000000003)", RESET;
print GREEN, "\n\t\tIf you want to design probes for all genes, use: --all_genes=yes", RESET;
print GREEN, "\n\t\tIf you wish to allow predicted genes (for species with few known genes), use: --allow_predicted_genes=yes", RESET;
print GREEN, "\n\t\tIf you want to design probes for a specific list of genes, specify the file containing this list using: --gene_file", RESET;
print GREEN, "\n\t\tNote: This file must be tab-delimited, with a header line and contain valid EnsEMBL Gene IDs in the first column", RESET;
print GREEN, "\n\t\tIf you want verbose output, use: --verbose=yes", RESET;
print GREEN, "\n\tSpecify the path to the directory of current probe files (used to determine starting probe ID value) using: --probe_dir", RESET;
print GREEN, "\n\tSpecify the name of the output probe file: --outfile", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of the design process using: --logfile", RESET;
print GREEN, "\n\nExample: generate_IntronProbes.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --overlap=5  --target_tm=67.0  --tm_range=3.5  --target_length=36  --max_length_variance=8  --allow_predicted_genes=no  --all_genes=yes  --probe_dir=/home/user/alexa/ALEXA_version/unfiltered_probes/  --outfile=/home/user/alexa/ALEXA_version/unfiltered_probes/intronProbes.txt  --logfile=/home/user/alexa/ALEXA_version/logs/generate_IntronProbes_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $overlap && $target_tm && $tm_range && $target_length && ($max_length_variance || $max_length_variance eq "0") && $allow_predicted_genes && $probe_dir && $outfile && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Get the starting probe_id and probeset_id by examining the specified input probe file directory.  If it is empty start with 1
my @current_ids = &getCurrentProbeIds('-probe_dir'=>$probe_dir);
my $current_probe_id = $current_ids[0];          #Unique probe ID used for each successful probe
my $current_probeset_id = $current_ids[1];       #Unique count of exon-exon junctions with successful probes

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#Print out the parameters supplied by the user to the logfile for future reference
print LOG "\nUser Specified the following options:\ndatabase = $database\noverlap = $overlap\ntarget_tm = $target_tm\ntm_range = $tm_range\ntarget_length = $target_length\nmax_length_variance = $max_length_variance\nensembl_gene_id = $ensembl_gene_id\nall_genes = $all_genes\nallow_predicted_genes = $allow_predicted_genes\ngene_file = $gene_file\nprobe_dir = $probe_dir\noutfile = $outfile\nlogfile = $logfile\n\n";

my @gene_ids;
my $probe_number = 0 ;

#The user may test a single gene or get all of them according to certain options
if ($ensembl_gene_id){
  my @ids;
  push (@ids, $ensembl_gene_id);
  my $gene_id_ref = &getGeneIds ('-dbh'=>$alexa_dbh, '-ensembl_g_ids'=>\@ids);
  my $gene_id = $gene_id_ref->{$ensembl_gene_id}->{alexa_gene_id};
  push (@gene_ids, $gene_id);

}elsif ($gene_file){
  @gene_ids = &parse_gene_file('-gene_file'=>$gene_file, '-dbh'=>$alexa_dbh);
  my $gene_count = @gene_ids;

}elsif ($all_genes eq "yes"){

  if ($allow_predicted_genes eq "yes"){
    @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'Non-pseudo')};
    my $gene_count = @gene_ids;
    print BLUE, "\nFound $gene_count genes that meet the criteria: 'Non-pseudo'\n\n", RESET;
    print LOG "\nFound $gene_count genes that meet the criteria: 'Non-pseudo'\n\n";

  }else{
    @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'Non-pseudo', '-evidence'=>"Known Gene")};
    my $gene_count = @gene_ids;
    print BLUE, "\nFound $gene_count genes that meet the criteria: 'Non-pseudo' and 'Known Gene'\n\n", RESET;
    print LOG "\nFound $gene_count genes that meet the criteria: 'Non-pseudo' and 'Known Gene'\n\n";
  }


}else{
  print RED, "\nSelect a single ensembl gene (--ensembl_gene_id=ENSG00000000003), provide an input file (--gene_file=), or use all genes (--all_genes=yes)\n\n", RESET;
  close (LOG);
  $alexa_dbh->disconnect();
  exit();
}

#Keep track of the number of exon regions for which probes are designed.  This script will design many probes spanning each of these regions.
#A second script will try to choose the 'best' single probe from each region
my $total_intron_regions = 0;
my $total_successful_intron_regions = 0;
my $gene_count = 0;

#Open the probe output file
print BLUE, "\nAll data will be written to $outfile\n\n", RESET;
print LOG "\nAll data will be written to $outfile\n\n";

open (OUTFILE, ">$outfile") || die "\nCould not open $outfile";
print OUTFILE "Probe_Count\tProbeSet_ID\tGene_ID\tSequence\tProbe_length\tProbe_Tm\tProbe_Type\tExon1_IDs\tUnit1_start\tUnit1_end\tExon2_IDs\tUnit2_start\tUnit2_end\tmasked_bases\n";

print BLUE, "\nGetting gene sequence data", RESET;
print LOG "\nGetting gene sequence data";
my $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids,'-sequence'=>"yes");

#Get the complete masked sequence for all genes
print BLUE, "\nGetting masked gene sequence data", RESET;
print LOG "\nGetting masked gene sequence data";
my $masked_genes_ref = &getMaskedGene ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);

#create an intron content object, containing the intronic regions of each gene
print BLUE, "\nGetting intron content positions for each gene", RESET;
print LOG "\nGetting intron content positions for each gene";
my $introns_ref = &getIntronContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);

foreach my $gene_id (@gene_ids){

  $gene_count++;

  #Go through each exon cluster - for clusters of a single exon get the probe sequences, for clusters of multiple exons, define regions
  my $intron_content_ref = $introns_ref->{$gene_id}->{intron_content};

  #IF a gene has no introns (single exon genes) - skip it
  if ($introns_ref->{$gene_id}->{intron_count} < 1){next;}

  my $number_introns = keys %{$intron_content_ref};
  my $intron_regions_successfully_targeted = 0;

  print CYAN, "\nGENE: $gene_count\tALEXA: $gene_id\tENSG: $genes_ref->{$gene_id}->{ensembl_g_id} has $number_introns introns. Extracting probe sequences", RESET;
  print LOG "\nGENE: $gene_count\tALEXA: $gene_id\tENSG: $genes_ref->{$gene_id}->{ensembl_g_id} has $number_introns introns. Extracting probe sequences";

  foreach my $intron_count (sort {$a <=> $b} keys %{$intron_content_ref}){
    $total_intron_regions++;

    print BLUE, "\n\n\tIntron: $intron_count\tSTART: $intron_content_ref->{$intron_count}->{start}\tEND: $intron_content_ref->{$intron_count}->{end}", RESET;
    print LOG "\n\n\tIntron: $intron_count\tSTART: $intron_content_ref->{$intron_count}->{start}\tEND: $intron_content_ref->{$intron_count}->{end}";

    my $probes_ref;

    $probes_ref = &getRawProbes('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id,
				'-start_pos'=>$intron_content_ref->{$intron_count}->{start}, '-end_pos'=>$intron_content_ref->{$intron_count}->{end});
    my $probes_found = 0;

    if (%{$probes_ref}){
      $probes_found = keys %{$probes_ref};
    }
    print BLUE, "\n\t\tFound $probes_found probes for this region", RESET;
    print LOG "\n\t\tFound $probes_found probes for this region";

    if ($probes_found >= 1){
      $intron_regions_successfully_targeted++;
      $total_successful_intron_regions++;
      &printProbeInfo('-probe_object'=>$probes_ref);
    }
  }

  print CYAN, "\n\nIntron regions found: $number_introns\tIntron regions with at least one successful probe: $intron_regions_successfully_targeted", RESET;
  print CYAN, "\n\n****************************************************************************************************************************************\n", RESET;
  print LOG "\n\nIntron regions found: $number_introns\tIntron regions with at least one successful probe: $intron_regions_successfully_targeted";
  print LOG "\n\n****************************************************************************************************************************************\n";

}

print CYAN, "\n\nA total of $probe_number probes were printed to the output file", RESET;
print CYAN, "\nTotal Intron Regions Found: $total_intron_regions.\tTotal with at least one probe found: $total_successful_intron_regions\n\n", RESET;
print LOG "\n\nA total of $probe_number probes were printed to the output file";
print LOG "\nTotal Intron Regions Found: $total_intron_regions.\tTotal with at least one probe found: $total_successful_intron_regions\n\n";

#Close database connection
$alexa_dbh->disconnect();

#Close the output file
close (OUTFILE);
close (LOG);
exit();


################################################################################################
#Parse input list of gene IDs to be used for intron probe selection                            #
################################################################################################
sub parse_gene_file{
  my %args = @_;
  my $gene_file = $args{'-gene_file'};
  my $dbh = $args{'-dbh'};

  print BLUE, "\nGetting gene IDs to use for intron probe design from the file: $gene_file", RESET;
  print LOG "\nGetting gene IDs to use for intron probe design from the file: $gene_file";

  my $genes_found = 0;
  my @gene_ids;
  my @ensembl_ids;

  open (GENES, "$gene_file") || die "\nCould not open input file: $gene_file\n\n";

  my $first_line = 1;
  while (<GENES>){

    if ($first_line == 1){
      $first_line = 0;
      next();
    }

    my @line = split ("\t", $_);
    if ($line[0] =~ /\d+/){
      push (@ensembl_ids, $line[0]);
      $genes_found++;
    }
  }
  close (GENES);

  foreach my $ensembl_g_id (@ensembl_ids){
    my @ids;
    push (@ids, $ensembl_g_id);
    my $gene_id_ref = &getGeneIds ('-dbh'=>$alexa_dbh, '-ensembl_g_ids'=>\@ids);
    my $gene_id = $gene_id_ref->{$ensembl_gene_id}->{alexa_gene_id};
    push (@gene_ids, $gene_id);

  }

  my @ordered_gene_ids = sort {$a <=> $b} @gene_ids;

  print BLUE, "\nFound $genes_found Gene IDs in this file\n\n", RESET;
  print LOG "\nFound $genes_found Gene IDs in this file\n\n";

  return(@ordered_gene_ids);
}


################################################################################################
#Given a gene object and a sequence start and stop position, return valid probes.
#Disqualify probes that have genomic N's, too many repeat masked N's, and Tm outside
#the accepted range.  Return a hash of probes sequences with start/end coordinates, Tm, etc.
################################################################################################
sub getRawProbes{

  my %args = @_;
  my $gene_id = $args{'-gene_id'};
  my $start_pos = $args{'-start_pos'};
  my $end_pos = $args{'-end_pos'};

  my $probe_count = 0;
  my %probe_object;

  #Probe characteristics defined globally ($target_length, $target_tm, $max_length_variance, $overlap)

  #If the desired probe length is longer than the sequence provided there is no point in continuing
  my $sequence_length = ($end_pos+1) - $start_pos;
  if (($target_length + $max_length_variance) > $sequence_length){
    return(\%probe_object);
  }

  #Get the source gene sequence to work with
  my $gene_seq = $genes_ref->{$gene_id}->{sequence};
  my $masked_seq = $masked_genes_ref->{$gene_id}->{sequence};

  #Define last position in the exon region that it is safe to select probes from without running past
  my $start_cutoff = ($end_pos - ($target_length + $max_length_variance + 1));

  #Parse probe sequences of $probe_length, starting every $overlap bases
 POS: for (my $i = $start_pos-1; $i < $start_cutoff; $i += $overlap){

    my %position_probes_object;
    my $position_probe_count = 0;  #Count of probes with varying length at a single position

  LENGTH: for (my $j = $target_length-$max_length_variance; $j <= $target_length+$max_length_variance; $j++){
      # Note: $i is the current probe start position and $j is the current probe length

      my $probe_seq = substr ($gene_seq, $i, $j);

      #Get the start and end position of the probe within the gene
      my $probe_seq_start = $i + 1;
      my $probe_seq_end = $i + $j;

      #Each probe must now pass through a number of filters to make it into the final list of candidates
      #If a probe fails any of the following basic checks - skip to the next probe

      #1.) Check for presence of genomic N's and other non-valid letters such as Ambiguiety codes
      if($probe_seq =~ /[^ATCG]/){
	if ($verbose eq "yes"){
	  print YELLOW, "\n\t\tThis probe: $probe_seq contains an invalid bases such as an N or R", RESET;
	  print LOG "\n\t\tThis probe: $probe_seq contains an invalid bases such as an N or R";
	}
	next(LENGTH); #Note if these are found, making the probe longer wont help, but we dont want to abandon any probes that were already successful
      }

      #2.) Check for presence of excessive repeatMasked bases
      #Specify the number of allowable RepeatMasked bases - Experiment with thresholds 
      my $allowed_masked_bases = ($j/4);  #If more than 1/4 of probe is masked, it will be rejected
      my $masked_seq = substr ($masked_seq, $i, $j);
      my @n_count = $masked_seq =~ m/N/g;
      my $n_count = @n_count;
      if ($n_count > $allowed_masked_bases){
	if ($verbose eq "yes"){
	  print YELLOW, "\n\t\tToo many masked bases ($n_count)", RESET;
	  print LOG "\n\t\tToo many masked bases ($n_count)";
	}
	next(LENGTH); #Note making the probe longer might help in this case, but we dont want to abandon any probes that were already successful
      }

      #3.) Calculate the Tm of this probe
      my $temp_k = &tmCalc('-sequence'=>$probe_seq, '-silent'=>1);
      my $Tm_celsius = &tmConverter('-tm'=>$temp_k, '-scale'=>'Kelvin');

      #4.) Check for probes with simple repeats
      my $test_seq = $probe_seq;
      #Search for Monopolymeric repeats of 6 or longer
      if ($test_seq =~ /(A{6,}|T{6,}|C{6,}|G{6,})/){
	if ($verbose eq "yes"){
	  print YELLOW, "\n\t\tMonopolymer detected ($1): $probe_seq", RESET;
	  print LOG "\n\t\tMonopolymer detected ($1): $probe_seq";
	}
	next(LENGTH); #Note if these are found, making the probe longer wont help, but we dont want to abandon any probes that were already successful
      }
      #Search for dipolymeric repeats of 4 or longer
      my $match = "AT";
      my @matches = qw (AT AG AC TA TG TC CA CT CG GA GT GC);
      foreach my $match (@matches){
	if ($test_seq =~ /($match){4,}/g){
	  if ($verbose eq "yes"){
	    print YELLOW, "\n\t\tDipolymer detected ($1): $probe_seq", RESET;
	    print LOG "\n\t\tDipolymer detected ($1): $probe_seq";
	  }
	  next(LENGTH); #Note if these are found, making the probe longer wont help, but we dont want to abandon any probes that were already successful
	}
      }
      #Store this probe that passed basic quality check temporarily until the best one for this probe/length for this spot can be chosen
      $position_probe_count++;
      $position_probes_object{$position_probe_count}{sequence} = $probe_seq;
      $position_probes_object{$position_probe_count}{start_pos} = $probe_seq_start;
      $position_probes_object{$position_probe_count}{end_pos} = $probe_seq_end;
      $position_probes_object{$position_probe_count}{Tm} = $Tm_celsius;
      $position_probes_object{$position_probe_count}{masked_bases} = $n_count;
      $position_probes_object{$position_probe_count}{length} = $j;

    }#LENGTH variance loop - altering probe length to achieve target Tm

    #Make sure at least some probes were found for this position by varying the length
    my $position_probes_found = keys %position_probes_object;
    unless ($position_probes_found >= 1){
      if ($verbose eq "yes"){
	print YELLOW, "\n\t\tNo suitable probes found for start position: $i", RESET;
	print LOG "\n\t\tNo suitable probes found for start position: $i";
      }
      next(POS);
    }

    #Go through each probe for this position, calculate absolute Tm difference from target, and note the 'best' probe
    my $best_probe_id;
    my $closest_tm = 1000;  #Arbitrarily large tm diff
    foreach my $pp_count (keys %position_probes_object){
      my $tm_diff = abs($position_probes_object{$pp_count}{Tm} - $target_tm);
      if ($tm_diff < $closest_tm){
	$closest_tm = $tm_diff;
	$best_probe_id = $pp_count;
      }
    }

    #Note: Even after selecting the best probe for each position, we still want to eliminate probes that are outside the desired Tm range
    #This will eliminate probes where even altering the length was not sufficient to get the target Tm.
    my $Tm_lower_limit = ($target_tm - $tm_range);
    my $Tm_upper_limit = ($target_tm + $tm_range);
    unless ($position_probes_object{$best_probe_id}{Tm} >= $Tm_lower_limit &&  $position_probes_object{$best_probe_id}{Tm} <= $Tm_upper_limit){
      if ($verbose eq "yes"){
	print YELLOW, "\n\t\tBest probe at this position has Tm: $position_probes_object{$best_probe_id}{Tm} outside of range ($Tm_lower_limit - $Tm_upper_limit)", RESET;
	print LOG "\n\t\tBest probe at this position has Tm: $position_probes_object{$best_probe_id}{Tm} outside of range ($Tm_lower_limit - $Tm_upper_limit)";
      }
      next(POS);
    }

    #If the probe made it this far, add it to the list of successful probes for this exon region and proceed to the next position within the region
    #Only select one probe for each tiling position - the one with the length that results in a Tm closest to the target_Tm
    $probe_count++;
    $probe_object{$probe_count}{gene_id} = $gene_id;
    $probe_object{$probe_count}{sequence} = $position_probes_object{$best_probe_id}{sequence};
    $probe_object{$probe_count}{start_pos} = $position_probes_object{$best_probe_id}{start_pos};
    $probe_object{$probe_count}{end_pos} = $position_probes_object{$best_probe_id}{end_pos};
    $probe_object{$probe_count}{Tm} = $position_probes_object{$best_probe_id}{Tm};
    $probe_object{$probe_count}{masked_bases} = $position_probes_object{$best_probe_id}{masked_bases};
    $probe_object{$probe_count}{length} = $position_probes_object{$best_probe_id}{length};

    #print "\nProgbeSEQ: $probe_seq";

  }#POS (Position) Tiling, Overlap loop - shifting position within the target exon region

  return(\%probe_object);
}


################################################################################################
#Given a probe object containing probe sequences for an exon region, print to a probe file.
################################################################################################
sub printProbeInfo{
  my %args = @_;
  my $probe_object_ref = $args{'-probe_object'};

  #Info that needs to be printed
  #Probe_Count\tProbeSet_ID\tGene_ID\tSequence\tProbe_length\tProbe_Tm\tProbe_Type\tExon1_IDs\tUnit1_start\tUnit1_end\tExon2_IDs (na)\tUnit2_start (na)\tUnit2_end (na)\tmasked_bases\n";

  #Global file handle = OUTFILE

  $current_probeset_id++;

  my $probe_type = "Intron";

  foreach my $probe (sort {$a<=>$b} keys %{$probe_object_ref}){
    $probe_number++;
    $current_probe_id++;

    my $Tm = $probe_object_ref->{$probe}->{Tm};
    my $Tm_rounded = sprintf("%.3f", $Tm);

   print OUTFILE "$current_probe_id\t$current_probeset_id\t$probe_object_ref->{$probe}->{gene_id}\t$probe_object_ref->{$probe}->{sequence}\t$probe_object_ref->{$probe}->{length}\t$Tm_rounded\t$probe_type\tna\t$probe_object_ref->{$probe}->{start_pos}\t$probe_object_ref->{$probe}->{end_pos}\tna\tna\tna\t$probe_object_ref->{$probe}->{masked_bases}\n";

  }
  return();
}
