#!/usr/bin/perl -w

#Written by Malachi Griffith
#The purpose of this script is to generate all possible exon-exon junction probes for each ensembl gene
#The user can specify a target Tm and probe length (which must be an even number for simplicity)
#Probes will be centered on the junction and the length will be modified to achieve the closest Tm to the target
#The probe length will not be allowed to extend or shorten beyond a limit set by the user (max_length_variance)
#This allows the user to retrieve isothermal probes for each junction.
#However, if the user wishes to have probes of a set length regardless of the Tm, they can simply specify a max_length_variance of 0

#The user will also be allowed to specify the number of probes returned for each junction, generated by increasing or decreasing the probe length
#These will be returned in order according to how closely they match the target Tm
#This is only applicable if the user has allowed the probe length to vary

#To allow for testing of this script, the user can also specify a single target Ensembl gene ID or design probes for all Ensembl Genes
#NOTE: If the user specifies to get probes for all genes, only 'Known' and 'Non-Pseudo' genes will be analyzed
#      - If you wish to design probes for predicted genes or pseudogenes, you will have to modifify this script

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
my $target_tm = '';
my $target_length = '';
my $max_length_variance = '';
my $probes_per_junction = '';
my $ensembl_gene_id = '';
my $all_genes = '';
my $allow_predicted_genes = '';
my $verbose = '';
my $probe_dir = '';
my $outfile = '';
my $logfile = '';
my $ignore_masking = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'target_tm=f'=>\$target_tm, 'target_length=i'=>\$target_length, 'max_length_variance=i'=>\$max_length_variance,
	    'ensembl_gene_id=s'=>\$ensembl_gene_id, 'all_genes=s'=>\$all_genes, 'allow_predicted_genes=s'=>\$allow_predicted_genes,
	    'probes_per_junction=i'=>\$probes_per_junction, 'verbose=s'=>\$verbose, 'probe_dir=s'=>\$probe_dir, 'outfile=s'=>\$outfile,
	    'logfile=s'=>\$logfile, 'ignore_masking=s'=>\$ignore_masking);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the target Tm for probe sequences using: --target_tm (eg. 67.0)", RESET;
print GREEN, "\n\tSpecify the target probe length (must be a multiple of 2) using: --target_length (eg. 36)", RESET;
print GREEN, "\n\tSpecify the maximum this probe length will be allowed to vary using: --max_length_variance (eg. 10)", RESET;
print GREEN, "\n\tSpecify the desired number of probes per junction to return using: --probes_per_junction (eg. 2)", RESET;
print GREEN, "\n\t\tIf you want to test with a single gene ID use: --ensembl_gene_id (ENSG00000000003)", RESET;
print GREEN, "\n\t\tIf you want to design probes for all genes, use: --all_genes=yes", RESET;
print GREEN, "\n\t\tIf you wish to allow predicted genes (for species with few known genes), use: --allow_predicted_genes=yes", RESET;
print GREEN, "\n\t\tIf you want verbose output, use: --verbose=yes", RESET;
print GREEN, "\n\t\tIf you want to disregard masking for some reason (not recommended) use: --ignore_masking=yes", RESET;
print GREEN, "\n\tSpecify the path to the directory of current probe files (used to determine starting probe ID value) using: --probe_dir", RESET;
print GREEN, "\n\tSpecify the name of the output probe file: --outfile", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of the design process using: --logfile", RESET;
print GREEN, "\n\nExample: generate_ExonJunctionProbes.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --target_tm=67.0  --target_length=36  --max_length_variance=10  --probes_per_junction=3  --allow_predicted_genes=no  --all_genes=yes  --probe_dir=/home/user/alexa/ALEXA_version/unfiltered_probes/  --outfile=/home/user/alexa/ALEXA_version/unfiltered_probes/exonJunctionProbes.txt  --logfile=/home/user/alexa/ALEXA_version/logs/generate_ExonJunctionProbes_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $target_tm && $target_length && ($max_length_variance || $max_length_variance eq "0") && $probes_per_junction && $allow_predicted_genes && $probe_dir && $outfile && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}
if ($max_length_variance eq "0" && $probes_per_junction > 1){
  print RED, "\nIf max variance in length is specified as 0, only 1 probe can be found!!\n\n", RESET;
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
print LOG "\nUser Specified the following options:\ndatabase = $database\ntarget_tm = $target_tm\ntarget_length = $target_length\nmax_length_variance = $max_length_variance\nprobes_per_junction = $probes_per_junction\nallow_predicted_genes = $allow_predicted_genes\nensembl_gene_id = $ensembl_gene_id\nall_genes = $all_genes\noutfile = $outfile\nprobe_dir = $probe_dir\nlogfile = $logfile\n\n";

my @gene_ids;

#The user may test a single gene or get all of them according to certain options
if ($ensembl_gene_id){
  my @ids;
  push (@ids, $ensembl_gene_id);
  my $gene_id_ref = &getGeneIds ('-dbh'=>$alexa_dbh, '-ensembl_g_ids'=>\@ids);
  my $gene_id = $gene_id_ref->{$ensembl_gene_id}->{alexa_gene_id};
  push (@gene_ids, $gene_id);

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
  print RED, "\nMust select either a single ensembl gene (--ensembl_gene_id=ENSG00000000003), or all genes option (--all_genes=yes)\n\n", RESET;
  close (LOG);
  $alexa_dbh->disconnect();
  exit();
}

my $gene_count = 0;
my $total_successful_probes = 0;
my $total_possible_probes = 0;

#Open the probe output file
print BLUE, "\nAll data will be written to $outfile\n\n", RESET;
print LOG "\nAll data will be written to $outfile\n\n";

open (OUTFILE, ">$outfile") || die "\nCould not open $outfile";
print OUTFILE "Probe_Count\tProbeSet_ID\tGene_ID\tSequence\tProbe_length\tProbe_Tm\tProbe_Type\tExon1_IDs\tUnit1_start\tUnit1_end\tExon2_IDs\tUnit2_start\tUnit2_end\tmasked_bases\n";

#Get the gene sequence and other gene info for all genes
print BLUE, "\nGetting gene sequence data", RESET;
print LOG "\nGetting gene sequence data";
my $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids,'-sequence'=>"yes");

#Get the complete masked sequence for all genes
my $masked_genes_ref;
unless ($ignore_masking eq "yes"){
  print BLUE, "\nGetting masked gene sequence data", RESET;
  print LOG "\nGetting masked gene sequence data";
  $masked_genes_ref = &getMaskedGene ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);
}

#Get the exons for all genes
print BLUE, "\nGetting exon coordinate data", RESET;
print LOG "\nGetting exon coordinate data";
my $gene_exons_ref = &getExons ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

#Get the total theoretical probes possible for all genes
print BLUE, "\nGetting theoretical probe counts data", RESET;
print LOG "\nGetting theoretical probe counts data";
my $probe_counts_ref = &junctionProbeCombinations('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);

foreach my $gene_id (@gene_ids){
  $gene_count++;

  #Keep track of exons, successful probes and unsuccessful probes for this gene
  my $successful_probes = 0;

  print CYAN, "\n\n\n*****************************************************************************", RESET;
  print CYAN, "\nGENE: $gene_count\tALEXA: $gene_id\tENSG: $genes_ref->{$gene_id}->{ensembl_g_id}", RESET;
  print LOG "\n\n\n*****************************************************************************";
  print LOG "\nGENE: $gene_count\tALEXA: $gene_id\tENSG: $genes_ref->{$gene_id}->{ensembl_g_id}";

  my $gene_seq = $genes_ref->{$gene_id}->{sequence};

  my $maskedGeneSeq;
  if ($ignore_masking eq "yes"){
    $maskedGeneSeq = $gene_seq;
  }else{
    $maskedGeneSeq = $masked_genes_ref->{$gene_id}->{sequence};
  }

  #Get all exons for this gene
  my $exons_ref = $gene_exons_ref->{$gene_id}->{exons};
  my $exon_count = keys %{$exons_ref};

  #Go through each exon and find the total number of unique end positions
  #Keep track of the exon IDs associated with each start and end positions (store as an array to keep track of redundant start/stop positions)
  my %start_positions;
  my %end_positions;
  foreach my $exon_id (sort keys %{$exons_ref}){

    #Unique start positions and associated exon IDs
    my $exon_start = $exons_ref->{$exon_id}->{exon_start};
    if ($start_positions{$exon_start}{exon_ids}){
      my @tmp = @{$start_positions{$exon_start}{exon_ids}};
      push (@tmp, $exon_id);
      $start_positions{$exon_start}{exon_ids} = \@tmp;
    }else{
      my @tmp;
      push (@tmp, $exon_id);
      $start_positions{$exon_start}{exon_ids} = \@tmp;
    }

    #Unique end positions and associated exon IDs
    my $exon_end = $exons_ref->{$exon_id}->{exon_end};
    if ($end_positions{$exon_end}{exon_ids}){
      my @tmp = @{$end_positions{$exon_end}{exon_ids}};
      push (@tmp, $exon_id);
      $end_positions{$exon_end}{exon_ids} = \@tmp;
    }else{
      my @tmp;
      push (@tmp, $exon_id);
      $end_positions{$exon_end}{exon_ids} = \@tmp;
    }
  }

  #If the gene has only one exon - skip
  if ($exon_count == 1){
    print CYAN, "\n\tGene has only a single exon - skipping", RESET;
    next();
  }

  #NOTE: Some exons have the same start or end position!
  #Probes should be tracked primarily with respect to their gene coordinates as specified in the GeneProbe table
  #We will also attempt to track which exons are associated with each probe

  #Find every valid NON-Redundant exon-exon pair, consider every possible permutation and check for valid connections by start/end position
  #i.e. Find the number of valid connections between each end position and all start positions
  foreach my $exon1_end_pos (sort {$a <=> $b} keys %end_positions){

    #For each end position, compare to all unique start positions of exons, any that are greater represent valid connections
    foreach my $exon2_start_pos (sort {$a <=> $b} keys %start_positions){

      #Consider whether this combination of exon1 and exon2 is valid
      if ($exon2_start_pos > $exon1_end_pos){
	my %junction_probes;
	my $junction_probe_count = 0;

	for (my $i = $target_length-$max_length_variance; $i <= $target_length+$max_length_variance; $i+=2){

	  #Probe length is altered in each iteration, but the probe sequence is alway centered on the junction
	  my $probe_length = $i;

	  #Determine the start and end position of the probe halves based on the length of the probe required and the probe length variance specified
	  #First 1/2 of probe
	  my $unit1_start_pos = ($exon1_end_pos-($probe_length/2));
	  my $unit1_length = ($probe_length/2);
	  my $unit1_end_pos = $exon1_end_pos;
	  my $unit1_seq = substr($gene_seq, $unit1_start_pos, $unit1_length);
	  my $seq1_length = length($unit1_seq);

	  #Second 1/2 of probe
	  my $unit2_start_pos = $exon2_start_pos;
	  my $unit2_length = ($probe_length/2);
	  my $unit2_end_pos = ($unit2_start_pos + $unit2_length);
	  my $unit2_seq = substr($gene_seq, $unit2_start_pos-1, $unit2_length);
	  my $seq2_length = length($unit2_seq);
	
	  my $probe_seq = "$unit1_seq"."$unit2_seq";
	  my $probe_seq_length = length($probe_seq);

	  #Still need to consider the length of each exon.  There may be several exons that end/start at a particular position
	  #But not all of them are neccessarily long enough to allow a probe to be designed
	
	  #The exon_end_pos and exon_start_pos combination used for this probe are associated with at least two exons
	  #Get all the exons involved in this specific junction and check if the resulting unit_1_start position is within the range of the exon1 options
	  #Similarly, check if the resulting unit_2_end position is within the range of the exon2 possibilities
	  #If the probe exceeds the boundaries of all possible exons at each side it must be rejected, otherwise, the exons for which
	  #it is valid will be noted
	  my @exon1_ids = @{$end_positions{$exon1_end_pos}{exon_ids}};
	  my @exon2_ids = @{$start_positions{$exon2_start_pos}{exon_ids}};

	  my @valid_exon1_ids;
	  foreach my $exon1_id (@exon1_ids){
	    if ($unit1_start_pos <= $exons_ref->{$exon1_id}->{exon_start}){
	      if ($verbose eq "yes"){
		print YELLOW, "\n\tProbe exceeds 5' boundary of exon: $exon1_id\tPROBE = $probe_seq", RESET;
	      }
	      print LOG "\n\tProbe exceeds 5' boundary of exon: $exon1_id\tPROBE = $probe_seq";
	    }else{
	      push (@valid_exon1_ids, $exon1_id);
	    }
	  }
	  my @valid_exon2_ids;
	  foreach my $exon2_id (@exon2_ids){
	    if ($unit2_end_pos >= $exons_ref->{$exon2_id}->{exon_end}){
	      if ($verbose eq "yes"){
		print YELLOW, "\n\tProbe exceeds 3' boundary of exon: $exon2_id\tPROBE = $probe_seq", RESET;
	      }
	      print LOG "\n\tProbe exceeds 3' boundary of exon: $exon2_id\tPROBE = $probe_seq";
	    }else{
	      push (@valid_exon2_ids, $exon2_id);
	    }
	  }
	  my $valid_exon1_ids = @valid_exon1_ids;
	  my $valid_exon2_ids = @valid_exon2_ids;

	  #Make sure the probe is valid for at least one of the pairs of exons with redundant start/stop positions
	  unless ($valid_exon1_ids > 0 && $valid_exon2_ids > 0){
	    next();
	  }

	  #Sanity check of resulting probe length
	  my $result_length = length($probe_seq);
	  unless ($result_length == $probe_length){
	    print RED, "\n\tUnexpected probe length ($result_length)!!\n\n", RESET;
	    print RED, "\nPROBE: $probe_seq", RESET;
	    print RED, "\nUNIT1: $unit1_seq\tUNIT2: $unit2_seq", RESET;
	    print RED, "\nExon1: @exon1_ids\tExon2: @exon2_ids", RESET;
	    print RED, "\nUnit1_start=$unit1_start_pos\tUnit1_end=$unit1_end_pos", RESET;
	    print RED, "\nUnit2_start=$unit2_start_pos\tUnit2_end=$unit2_end_pos", RESET;
	    close (LOG);
	    $alexa_dbh->disconnect();
	    exit();
	  }

	  #NOTE: Some exons actually contain N's from the underlying genomic sequence in ensembl!  
	  #For simplicity, probes that incorporate these unknown bases should be skipped!

	  #Check for presence of genomic N's and other non-valid letters such as Ambiguiety codes
	  if($probe_seq =~ /[^ATCG]/){
	    if ($verbose eq "yes"){
	      print YELLOW, "\n\tThis probe: $probe_seq contains an invalid bases such as an N or R", RESET;
	    }
	    print LOG "\n\tThis probe: $probe_seq contains an invalid bases such as an N or R";
	    next();
	  }

	  #Check for excessive RepeatMasked bases in this genomic region
	  #Specify the number of allowable RepeatMasked bases - Experiment with thresholds 
	  my $allowed_masked_bases = ($probe_length/4);  #If more than 1/4 of the target probe is masked, it will be rejected

	  my $masked_unit1_seq = substr($maskedGeneSeq, $unit1_start_pos, $unit1_length);
	  my $masked_unit2_seq = substr($maskedGeneSeq, $unit2_start_pos-1, $unit2_length);
	  my $masked_probe_seq = "$masked_unit1_seq"."$masked_unit2_seq";
	  my @n_count = $masked_probe_seq =~ m/N/g;
	  my $n_count = @n_count;

	  if ($n_count > $allowed_masked_bases){
	    if ($verbose eq "yes"){
	      print YELLOW, "\nExon1: @exon1_ids\tExon2: @exon2_ids", RESET;
	      print YELLOW, "\n\tToo many masked bases ($n_count)", RESET;
	      print YELLOW, "\n\tProbe:  $probe_seq", RESET;
	      print YELLOW, "\n\tMasked: $masked_probe_seq", RESET;
	      print LOG "\nExon1: @exon1_ids\tExon2: @exon2_ids";
	      print LOG "\n\tToo many masked bases ($n_count)";
	      print LOG "\n\tProbe:  $probe_seq";
	      print LOG "\n\tMasked: $masked_probe_seq";
	    }
	    next();
	  }

	  $junction_probe_count++;

	  #If this probe passes all the simple filters here, calculate its Tm and add it to list of potential list of probes for this junction
	  my $temp_k = &tmCalc('-sequence'=>$probe_seq, '-silent'=>1);
	  my $Tm_celsius = &tmConverter('-tm'=>$temp_k, '-scale'=>'Kelvin');

	  $junction_probes{$junction_probe_count}{sequence} = $probe_seq;
	  $junction_probes{$junction_probe_count}{probe_seq_length} = $probe_seq_length;
	  $junction_probes{$junction_probe_count}{length} = $probe_length;
	  $junction_probes{$junction_probe_count}{tm} = $Tm_celsius;
	  $junction_probes{$junction_probe_count}{valid_exon1_ids} = \@valid_exon1_ids;
	  $junction_probes{$junction_probe_count}{exon1_end} = $exon1_end_pos;
	  $junction_probes{$junction_probe_count}{valid_exon2_ids} = \@valid_exon2_ids;
	  $junction_probes{$junction_probe_count}{exon2_start} = $exon2_start_pos;
	  $junction_probes{$junction_probe_count}{unit1_start} = $unit1_start_pos;
	  $junction_probes{$junction_probe_count}{unit1_end} = $unit1_end_pos;
	  $junction_probes{$junction_probe_count}{unit1_seq} = $unit1_seq;
	  $junction_probes{$junction_probe_count}{unit1_length} = $unit1_length;
	  $junction_probes{$junction_probe_count}{seq1_length} = $seq1_length;
	  $junction_probes{$junction_probe_count}{unit2_start} = $unit2_start_pos;
	  $junction_probes{$junction_probe_count}{unit2_end} = $unit2_end_pos;
	  $junction_probes{$junction_probe_count}{unit2_seq} = $unit2_seq;
	  $junction_probes{$junction_probe_count}{unit2_length} = $unit2_length;
	  $junction_probes{$junction_probe_count}{seq2_length} = $seq2_length;
	  $junction_probes{$junction_probe_count}{masked_bases} = $n_count;

	} #Probe Length variance loop

	#Now that probes have been generated for all possible lengths of probe for this junction, rank them according to Tm and select the desired number

	#First make sure the desired number of probes was found
	my $junction_probes_found = keys %junction_probes;
	unless ($junction_probes_found >= $probes_per_junction){
	  if ($verbose eq "yes"){
	    print YELLOW, "\n\tCould not find the desired number of probes for this junction", RESET;
	  }
	  print LOG "\n\tCould not find the desired number of probes for this junction";
	  next();
	}
	$current_probeset_id++;
	#Calculate absolute difference from target
	foreach my $jp_count (keys %junction_probes){
	  $junction_probes{$jp_count}{abs_target_tm_diff} = abs($junction_probes{$jp_count}{tm} - $target_tm);
	}

	my $probes_selected = 0;
	foreach my $jp_count (sort {$junction_probes{$a}->{abs_target_tm_diff} <=> $junction_probes{$b}->{abs_target_tm_diff}} keys %junction_probes){
	  #Print probe info to Output probe file
	  $probes_selected++;

	  if ($probes_selected > $probes_per_junction){
	    last();
	  }

	  my $tm_rounded = sprintf("%.3f", $junction_probes{$jp_count}{tm});

	  #Verbose Printing Section - Helpful for debugging ...
	  if ($verbose eq "yes"){
	    print YELLOW, "\n\nEXON1_end=$junction_probes{$jp_count}{exon1_end}\tEXON2_start=$junction_probes{$jp_count}{exon2_start}", RESET;
	    print YELLOW, "\nUNIT1: Start=$junction_probes{$jp_count}{unit1_start}\tEnd=$junction_probes{$jp_count}{unit1_end}\tLength=$junction_probes{$jp_count}{unit1_length}\tSeq=$junction_probes{$jp_count}{unit1_seq} ($junction_probes{$jp_count}{seq1_length})", RESET;
	    print YELLOW, "\nUNIT2: Start=$junction_probes{$jp_count}{unit2_start}\tEnd=$junction_probes{$jp_count}{unit2_end}\tLength=$junction_probes{$jp_count}{unit2_length}\tSeq=$junction_probes{$jp_count}{unit2_seq} ($junction_probes{$jp_count}{seq2_length})", RESET;
	    print YELLOW, "\nPROBE: $junction_probes{$jp_count}{sequence} ($junction_probes{$jp_count}{probe_seq_length}) (Tm = $tm_rounded)", RESET;
	  }

	  #Probe_Count,ProbeSet_ID,Gene_ID,Sequence,Probe_length,Probe_Tm,Probe_Type,Exon1_IDs,Unit1_start,Unit1_end,Exon2_IDs,Unit2_start,Unit2_end,masked_bases
	  $successful_probes++;
	  $total_successful_probes++;
	  $current_probe_id++;
	  my $probe_type = "Exon-Exon";

	  print OUTFILE "$current_probe_id\t$current_probeset_id\t$gene_id\t$junction_probes{$jp_count}{sequence}\t$junction_probes{$jp_count}{length}\t$tm_rounded\t$probe_type\t@{$junction_probes{$jp_count}{valid_exon1_ids}}\t$junction_probes{$jp_count}{unit1_start}\t$junction_probes{$jp_count}{unit1_end}\t@{$junction_probes{$jp_count}{valid_exon2_ids}}\t$junction_probes{$jp_count}{unit2_start}\t$junction_probes{$jp_count}{unit2_end}\t$junction_probes{$jp_count}{masked_bases}\n";

	}#Print probes for this junction loop
      }#Valid junction loop
    }#Junction start position loop (unit_2_start)
  }#Junction end position loop (unit_1_end)

  my $possible_probes = $probe_counts_ref->{$gene_id}->{exon_exon};

  #Possible probes must account for the number of desired probe per junction
  $possible_probes = ($possible_probes * $probes_per_junction);

  $total_possible_probes += $possible_probes;

  print CYAN, "\n\nSUMMARY for Gene ID: $gene_id", RESET;
  print CYAN, "\nNumber Exons = $exon_count\tPossible Exon-Exon Probes = $possible_probes\tSuccessful Probes = $successful_probes", RESET;
  print LOG "\n\nSUMMARY for Gene ID: $gene_id";
  print LOG "\nNumber Exons = $exon_count\tPossible Exon-Exon Probes = $possible_probes\tSuccessful Probes = $successful_probes";

}#Gene loop

print CYAN, "\n\nTotal Possible Probes = $total_possible_probes\nTotal Successful Probes = $total_successful_probes\n\n", RESET;
print LOG "\n\nTotal Possible Probes = $total_possible_probes\nTotal Successful Probes = $total_successful_probes\n\n";

close (OUTFILE);
close (LOG);

#Close database connection
$alexa_dbh->disconnect();

exit();
