#!/usr/bin/perl -w

#Written by Malachi Griffith
#The purpose of this script is to select probe sequences within ensembl exons
#The user will specify the desired target probe length (--target_length) and the degree of overlap/tiling (--overlap) (i.e. how closely to design probes to each other)
#If the value for --overlap is less than the --target_length, the designed probes will overlap each other
#At every design position, the length of the probe will be varied up to the --max_length_variance specified by the user to achieve the specified --target_tm
#If you want anisothermal probes of fixed length, simply set --max_variance_length=0
#The desired melting temperature is the median Tm of all exon-exon and intron-exon junction probes

#In this way, many probes will be designed for each exon of each gene.  From this pool, a selection of representative probes will be chosen
#In general I want at least a few probes per exon and they should be near to the centre of the exon and non-overlapin if possible.
#Exons that overlap each other complicate the process of probe selection and certain efforts must be made to get probes that represent each of the exon regions present

#Note that in general we can be a bit more stringent in the selection of exon probes because unlike for junctions there is a larger region to select from

#In general the neccessary code can be summarized as follows:

#- For each gene, get all exons and their gene coordinates via ALEXA
#(1.) For exons that are completely isolated from other exons (no overlapping coordinates)
#     - Extract probe sequences of length '--target_length' +/- '--max_length_variance', starting every '--overlap' bases as specified by the user
#     - Make sure all probe sequences are within the boundaries of the exon (also check for exons that are too short)
#     - Store the probe sequence, start position, and end position as a hash
#     - Eliminate probes that correspond to Genomic N's
#     - Eliminate probes that are comprised of too many repeat masked bases
#     - Eliminate probes that are outside the desired target Tm +/- the acceptable range
#     - Eliminate probes with simple repeats such as monopolymers and dipolymers

#(2.) For all exons that have some overlap, create clusters of these exons.
#     - (a) Identify regions of these exons that are unique.  Treat these regions as 'exons' and try to design probes as in (1.)
#     - (b) Also identify all junction points for the exons in the cluster - This will be a non-redundant list of the start/end positions
#           - Treat the regions between these as 'exons' and try to design probes as in (1.)
#     - These two approaches will sometimes fail.
#           - (a) may fail if there is not a unique portion of each of the overlapping exons or it is too small
#                 - Any probes that were successful should be useful though
#           - (b) may fail if the junctions are too close to each other
#                 - Any probes that were successful should be useful though
#     - (c) Identify exons that were not successfully targeted in (a) or (b) and simply design probes for each of these exons as in (1.)?
#     Note: When designing probes for overlapping exons I want to avoid probes that span the junctions!  These probes should have been designed already
#           as intron-exon probes

#Basically the main goal of this script is to identify exon regions that will interogate the expression of individual exons in the most informative way
#possible.  Many probes will be designed across each exon region (either a whole exon or parts of exons from multi-exon clusters)
#-The end result of this script will be a file containing many probes for each exon region of each gene.  A seperate script will then evaluate these:
#   - Essentially we want the 'best' probe(s) for each exon region
#   - We want the probes with the highest specificity (ie. unique in the transcriptome if possible -avoid common sequences from gene families
#   - Also the 'best' probes are those that have a Tm that is as close as possible to the target Tm and the probe is close to the center of the exon
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
my $all_genes = '';
my $allow_predicted_genes = '';
my $verbose = '';
my $probe_dir = '';
my $outfile = '';
my $logfile = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'overlap=i'=>\$overlap, 'target_tm=f'=>\$target_tm, 'tm_range=f'=>\$tm_range, 'target_length=i'=>\$target_length,
	    'max_length_variance=i'=>\$max_length_variance, 'ensembl_gene_id=s'=>\$ensembl_gene_id, 'all_genes=s'=>\$all_genes,
	    'allow_predicted_genes=s'=>\$allow_predicted_genes, 'verbose=s'=>\$verbose, 'probe_dir=s'=>\$probe_dir, 'outfile=s'=>\$outfile,
	    'logfile=s'=>\$logfile);

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
print GREEN, "\n\t\tIf you want verbose output, use: --verbose=yes", RESET;
print GREEN, "\n\tSpecify the path to the directory of current probe files (used to determine starting probe ID value) using: --probe_dir", RESET;
print GREEN, "\n\tSpecify the name of the output probe file: --outfile", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of the design process using: --logfile", RESET;
print GREEN, "\n\nExample: generate_ExonProbes.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --overlap=5  --target_tm=67.0  --tm_range=3.5  --target_length=36  --max_length_variance=8  --all_genes=yes  --allow_predicted_genes=no  --probe_dir=/home/user/alexa/ALEXA_version/unfiltered_probes/  --outfile=/home/user/alexa/ALEXA_version/unfiltered_probes/exonProbes.txt  --logfile=/home/user/alexa/ALEXA_version/logs/generate_ExonProbes_LOG.txt\n\n", RESET;

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
print LOG "\nUser Specified the following options:\nmy database = $database\noverlap = $overlap\ntarget_tm = $target_tm\ntm_range = $tm_range\ntarget_length = $target_length\nmax_length_variance = $max_length_variance\nensembl_gene_id = $ensembl_gene_id\nall_genes = $all_genes\nallow_predicted_genes = $allow_predicted_genes\nverbose = $verbose\nprobe_dir = $probe_dir\noutfile = $outfile\nlogfile = $logfile\n\n";

my @gene_ids;
my $probe_number = 0 ;

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

#Keep track of the number of exon regions for which probes are designed.  This script will design many probes spanning each of these regions.
#A second script will try to choose the 'best' single probe from each region
my $total_exon_regions = 0;
my $total_successful_exon_regions = 0;
my $gene_count = 0;

#Open the probe output file
print BLUE, "\nAll data will be written to $outfile\n\n", RESET;
print LOG "\nAll data will be written to $outfile\n\n";

open (OUTFILE, ">$outfile") || die "\nCould not open $outfile";
print OUTFILE "Probe_Count\tProbeSet_ID\tGene_ID\tSequence\tProbe_length\tProbe_Tm\tProbe_Type\tExon1_IDs\tUnit1_start\tUnit1_end\tExon2_IDs\tUnit2_start\tUnit2_end\tmasked_bases\n";

#Create a gene object containing the gene sequence, masked sequence, exon positions, exon clusters, etc. for all genes
my $genes_ref = &createGeneObject ('-gene_ids'=> \@gene_ids);

foreach my $gene_id (@gene_ids){

  $gene_count++;
  my $exon_regions_found = 0;
  my $exon_regions_successfully_targeted = 0;

  #IF a gene has only one exon - skip it
  #if ($genes_ref->{$gene_id}->{number_exons} <= 1){next;} #Not neccessary

  #Go through each exon cluster - for clusters of a single exon get the probe sequences, for clusters of multiple exons, define regions
  my $exons_ref = $genes_ref->{$gene_id}->{exons};

  my $number_exons = keys %{$exons_ref};

  print CYAN, "\nGENE: $gene_count\tALEXA: $gene_id\tENSG: $genes_ref->{$gene_id}->{ensembl_g_id} has $number_exons exons. Extracting probe sequences", RESET;
  print LOG "\nGENE: $gene_count\tALEXA: $gene_id\tENSG: $genes_ref->{$gene_id}->{ensembl_g_id} has $number_exons exons. Extracting probe sequences";

  my $exon_clusters_ref = $genes_ref->{$gene_id}->{clusters};

  foreach my $cluster (sort {$a <=> $b} keys %{$exon_clusters_ref}){
    my $cluster_exons_aref = $exon_clusters_ref->{$cluster}->{exons};
    my $cluster_size = @{$cluster_exons_aref};

    print BLUE, "\n\n\tProcessing cluster: $cluster consisting of $cluster_size exons", RESET;
    print LOG "\n\n\tProcessing cluster: $cluster consisting of $cluster_size exons";

    my $probes_ref;

    #1.) Single-exon clusters
    if ($cluster_size == 1){
      my $exon_id =@{$cluster_exons_aref}[0];

      print BLUE, "\n\t\tSingle Exon Region: ($exons_ref->{$exon_id}->{exon_start} - $exons_ref->{$exon_id}->{exon_end}) covers exon ids: $exon_id", RESET;
      print LOG "\n\t\tSingle Exon Region: ($exons_ref->{$exon_id}->{exon_start} - $exons_ref->{$exon_id}->{exon_end}) covers exon ids: $exon_id";
      $probes_ref = &getRawProbes('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$exons_ref->{$exon_id}->{exon_start},
				  '-end_pos'=>$exons_ref->{$exon_id}->{exon_end}, '-exon_ids'=>$cluster_exons_aref);

      $exon_regions_found++;
      my $probes_found = 0;

      if (%{$probes_ref}){
	$probes_found = keys %{$probes_ref};
      }
      print BLUE, "\n\t\tFound $probes_found probes for this region", RESET;
      print LOG "\n\t\tFound $probes_found probes for this region";

      if ($probes_found >= 1){
	$exon_regions_successfully_targeted++;
	&printProbeInfo('-probe_object'=>$probes_ref);
      }

    }else{
      #2.) Multi-exon clusters
      #Get the most informative regions from the cluster of overlapping exons (non-overlaping where possible)
      my $exon_regions_ref = &defineExonRegions('-exon_object'=>$exons_ref, '-exon_ids'=>$cluster_exons_aref);

      #For each region identified attempt to extract probe sequences as done for single-exon clusters
      foreach my $region (sort {$a <=> $b} keys %{$exon_regions_ref}){
	my $region_start = $exon_regions_ref->{$region}->{region_start};
	my $region_end = $exon_regions_ref->{$region}->{region_end};
	my $region_exons_aref = $exon_regions_ref->{$region}->{exons};

	print BLUE, "\n\t\tExon Region: ($exon_regions_ref->{$region}->{region_start} - $exon_regions_ref->{$region}->{region_end}) covers exon ids: @{$region_exons_aref}", RESET;
	print LOG "\n\t\tExon Region: ($exon_regions_ref->{$region}->{region_start} - $exon_regions_ref->{$region}->{region_end}) covers exon ids: @{$region_exons_aref}";
	$probes_ref = &getRawProbes('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$region_start,
				    '-end_pos'=>$region_end, '-exon_ids'=>$region_exons_aref);

	$exon_regions_found++;
	my $probes_found = 0;

	if (%{$probes_ref}){
	  $probes_found = keys %{$probes_ref};
	}
	print BLUE, "\n\t\tFound $probes_found probes for this region", RESET;
	print LOG "\n\t\tFound $probes_found probes for this region";

	if ($probes_found >= 1){
	  $exon_regions_successfully_targeted++;
	  &printProbeInfo('-probe_object'=>$probes_ref);
	}
      }
    }
  }
  print CYAN, "\n\nExon regions found: $exon_regions_found\tExon regions with at least one successful probe: $exon_regions_successfully_targeted", RESET;
  print CYAN, "\n****************************************************************************************************************************************\n", RESET;
  print LOG "\n\nExon regions found: $exon_regions_found\tExon regions with at least one successful probe: $exon_regions_successfully_targeted";
  print LOG "\n****************************************************************************************************************************************\n";
  $total_exon_regions += $exon_regions_found;
  $total_successful_exon_regions += $exon_regions_successfully_targeted;
}

print CYAN, "\n\nA total of $probe_number probes were printed to the output file", RESET;
print CYAN, "\nTotal Exon Regions Found: $total_exon_regions\tTotal with at least one probe found: $total_successful_exon_regions\n\n", RESET;
print LOG "\n\nA total of $probe_number probes were printed to the output file";
print LOG "\nTotal Exon Regions Found: $total_exon_regions\tTotal with at least one probe found: $total_successful_exon_regions\n\n";

#Close database connection
$alexa_dbh->disconnect();

#Close the output file
close (OUTFILE);
close (LOG);
exit();


################################################################################################
#Create gene object - return hash with gene sequence, masked sequence, exon info, etc.         #
#Identify overlapping exons and define exon clusters                                           #
################################################################################################
sub createGeneObject{
  my %args = @_;
  my @gene_ids = @{$args{'-gene_ids'}};

  my %gene_object;

  #Get the raw gene sequence
  print BLUE, "\nGetting gene sequence data", RESET;
  print LOG "\nGetting gene sequence data";

  my $gene_info_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids,'-sequence'=>"yes");

  #Get the complete masked sequence for this gene
  print BLUE, "\nGetting masked gene sequence data", RESET;
  print LOG "\nGetting masked gene sequence data";
  my $masked_genes_ref = &getMaskedGene ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);

  #Get all exons for this gene
  print BLUE, "\nGetting exon coordinate data", RESET;
  print LOG "\nGetting exon coordinate data";
  my $gene_exons_ref = &getExons ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

  #Process each gene
  foreach my $gene_id (@gene_ids){

    my $exons_ref = $gene_exons_ref->{$gene_id}->{exons};

    my $exon_count = keys %{$exons_ref};

    #Go through each exon and create clusters of overlapping exons.  Any exon that shares at least one bp with another exon will be placed in the same cluster
    #An exon that is completely isolated from others will comprise a cluster of one
    my $cluster_count = 0;
    my %clusters;

  EXON: foreach my $exon_id (sort {$exons_ref->{$a}->{exon_start} <=> $exons_ref->{$b}->{exon_start}} keys %{$exons_ref}){

      #Get the start and end positions of the current exon
      my $exon_start = $exons_ref->{$exon_id}->{exon_start};
      my $exon_end = $exons_ref->{$exon_id}->{exon_end};

      #Go through each cluster and see if the current exon overlaps with one of the exons already present in one of these clusters
      foreach my $cluster_id (sort keys %clusters){
	my @tmp_exons = @{$clusters{$cluster_id}{exons}};
	foreach my $exon_id_test (@tmp_exons){

	  my $exon_start_test = $exons_ref->{$exon_id_test}->{exon_start};
	  my $exon_end_test = $exons_ref->{$exon_id_test}->{exon_end};

	  #See if the start of the current exon is within the range of the test exon
	  if ($exon_start >= $exon_start_test && $exon_start <= $exon_end_test){
	    push (@tmp_exons, $exon_id);
	    $clusters{$cluster_id}{exons} = \@tmp_exons;
	    next(EXON);
	  }
	  #See if the end of the current exon is within the range of the test exon
	  if ($exon_end >= $exon_start_test && $exon_end <= $exon_end_test){
	    push (@tmp_exons, $exon_id);
	    $clusters{$cluster_id}{exons} = \@tmp_exons;
	    next(EXON);
	  }
	  #See if the current exon completely flanks the test exon - if so it should be added to the cluster
	  #New condition to fix bug!!
	  if ($exon_start <= $exon_start_test && $exon_end >= $exon_end_test){
	    push (@tmp_exons, $exon_id);
	    $clusters{$cluster_id}{exons} = \@tmp_exons;
	    next(EXON);
	  }
	}
      }
      #If the current exon was not added to any of the current clusters - create a new cluster
      $cluster_count++;
      my @tmp_exons;
      push (@tmp_exons, $exon_id);
      $clusters{$cluster_count}{exons} = \@tmp_exons;
    }

    #Build the gene object from the info gathered above
    $gene_object{$gene_id}{sequence} = $gene_info_ref->{$gene_id}->{sequence};
    $gene_object{$gene_id}{ensembl_g_id} = $gene_info_ref->{$gene_id}->{ensembl_g_id};
    $gene_object{$gene_id}{masked_sequence} = $masked_genes_ref->{$gene_id}->{sequence};
    $gene_object{$gene_id}{gene_start} = $gene_info_ref->{$gene_id}->{gene_start};
    $gene_object{$gene_id}{gene_end} = $gene_info_ref->{$gene_id}->{gene_end};
    $gene_object{$gene_id}{number_exons} = $exon_count;
    $gene_object{$gene_id}{exons} = $exons_ref;
    $gene_object{$gene_id}{clusters} = \%clusters;
  }

  return(\%gene_object);
}


################################################################################################
#Given a gene object and a sequence start and stop position, return valid probes.
#Disqualify probes that have genomic N's, too many repeat masked N's, and Tm outside
#the accepted range.  Return a hash of probes sequences with start/end coordinates, Tm, etc.
################################################################################################
sub getRawProbes{

  my %args = @_;
  my $gene_object_ref = $args{'-gene_object'};
  my $gene_id = $args{'-gene_id'};
  my $exon_ids_ref = $args{'-exon_ids'};
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
  my $gene_seq = $gene_object_ref->{$gene_id}->{sequence};
  my $masked_seq = $gene_object_ref->{$gene_id}->{masked_sequence};

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
    $probe_object{$probe_count}{exon_ids} = $exon_ids_ref;
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
#For exon clusters that have been derived from overlapping exons, identify informative regions #
#of these exons to use for probe design.                                                       #
#For each of these regions, make note of which of the exons from the cluster of overlapping    #
#exons are involved in each of the defined regions                                             #
################################################################################################
sub defineExonRegions{
  my %args = @_;

  my $exon_object_ref = $args{'-exon_object'};
  my $exon_ids_aref = $args{'-exon_ids'};

  #The exons in @exon_ids comprise a single cluster of overlapping exons
  #First try to identify regions of these exons that are unique to each exon

  #Compile a non-redundant list of start/end positions
  #Note that if two exons are in a cluster and have exactly the same start/end positions (which does happen in Ensembl!) then when
  #converted to a non-redundant list they will resolve to a single region which is good.
  my %junctions;
  my @junctions;
  foreach my $exon_id (@{$exon_ids_aref}){
    #Use a hash to compile a non-redundant list of start/end positions
    my $exon_start = $exon_object_ref->{$exon_id}->{exon_start};
    my $exon_end = $exon_object_ref->{$exon_id}->{exon_end};
    $junctions{$exon_start}{tmp}='na';
    $junctions{$exon_end}{tmp}='na';
    #print "\nEXON: $exon_id\tStart: $exon_start\tEnd: $exon_end";
 }
  #Create a sorted array of these junctions
  foreach my $junct (sort {$a <=> $b} keys %junctions){
    push (@junctions, $junct);
  }

  #Now consider the regions between the junctions and identify the exons associated with each
  my $number_junctions = @junctions;
  my %regions;
  my $region_count;

  for (my $i = 0; $i < $number_junctions-1; $i++){
    my $region_start = ($junctions[$i])+1;
    my $region_end = ($junctions[$i+1])-1;

    #print "\nCompare: $region_start - $region_end";

    #Confirm that the region selected is actually valid and larger than the required probe size
    if ($region_end <= $region_start){
      #print "\nExon junctions are too close to each other to allow probe design";
      next();
    }
    my $region_size = ($region_end - $region_start);
    if ($region_size <= ($target_length + $max_length_variance)){
      #print "\nExon junctions are too close to each other to allow probe design";
      next();
    }

    #Check each exon to see which overlap within this region at either end - or flank it completely
    my @exons_within_region;
    foreach my $exon_id (@{$exon_ids_aref}){
      my $exon_start = $exon_object_ref->{$exon_id}->{exon_start};
      my $exon_end = $exon_object_ref->{$exon_id}->{exon_end};

      #Is the start position of the selected region within this exon?
      if ($region_start >= $exon_start && $region_start <= $exon_end){
	push (@exons_within_region, $exon_id);
	next();
      }
      #Is the end position of the selected region within this exon?
      if ($region_end >= $exon_start && $region_end <= $exon_end){
	push (@exons_within_region, $exon_id);
	next();
      }
      #Does the selected region completely flank this exon?
      #New condition to fix bug!!
      if ($region_start <= $exon_start && $region_end >= $exon_end){
	push (@exons_within_region, $exon_id);
	next();
      }

    }
    $region_count++;
    $regions{$region_count}{region_start} = $region_start;
    $regions{$region_count}{region_end} = $region_end;
    $regions{$region_count}{exons} = \@exons_within_region;
    if ($verbose eq "yes"){
      print YELLOW, "\n\t\tRegion: $region_count ($region_start - $region_end) covers exons: @exons_within_region", RESET;
    }
    print LOG "\n\t\tRegion: $region_count ($region_start - $region_end) covers exons: @exons_within_region";
  }

  #print Dumper %regions;

  #Keep track of the exons that are successfully covered by the regions defined
  #Check each region defined to see if it is completely within an exon.
  #For every exon, I want a region that is completely within the boundaries of the exon.  For those exons where this is not true,
  #define a region within the exon and note which other exons it covers - for this evaluation, if a probe covers any amount of sequence
  #of another exon it will be noted.  This is a conservative approach, if there is any ambiguity at all regarding which exons a probe covers,
  #it must be noted.  Nevertheless I want at least one probe that is completely within each exon.
  my %exons_covered;
  foreach my $region (sort {$a <=> $b} keys %regions){

    my $region_start = $regions{$region}{region_start};
    my $region_end = $regions{$region}{region_end};

    #Check this region against the initial list of exons for this exon cluster
    foreach my $exon_id (@{$exon_ids_aref}){
      my $exon_start = $exon_object_ref->{$exon_id}->{exon_start};
      my $exon_end = $exon_object_ref->{$exon_id}->{exon_end};

      #Which exons completely cover this region?
      if ($region_start >= $exon_start && $region_end <= $exon_end){
	$exons_covered{$exon_id}{tmp} = 'na';
      }
    }
  }
  #Now go through the original list of exons for this exon cluster and see which have not been successfully covered
  foreach my $exon_id (@{$exon_ids_aref}){
    unless ($exons_covered{$exon_id}{tmp}){
      if ($verbose eq "yes"){
	print YELLOW, "\n\t\tInitial attempts to define a probe region for Exon: $exon_id failed.  Last attempt ...", RESET;
	print LOG "\n\t\tInitial attempts to define a probe region for Exon: $exon_id failed.  Last attempt ...";
      }

      #Still need a probe region for this exon!
      my $region_start = ($exon_object_ref->{$exon_id}->{exon_start})+1;
      my $region_end = ($exon_object_ref->{$exon_id}->{exon_end})-1;

      #Make sure this region fits the required probe length
      my $region_size = $region_end - $region_start;
      if ($region_size <= ($target_length + $max_length_variance)){
	if ($verbose eq "yes"){
	  print YELLOW, "\n\t\tExon is too small to allow specific probe design.  Exon: $exon_id", RESET;
	  print LOG "\n\t\tExon is too small to allow specific probe design.  Exon: $exon_id";
	}
	next();
      }

      #Define the region and note the exons that it covers.
      my @exons_within_region;
      foreach my $exon_id (@{$exon_ids_aref}){
	my $exon_start = $exon_object_ref->{$exon_id}->{exon_start};
	my $exon_end = $exon_object_ref->{$exon_id}->{exon_end};

	#Is the start position of the selected region within this exon?
	if ($region_start >= $exon_start && $region_start <= $exon_end){
	  push (@exons_within_region, $exon_id);
	  next();
	}
	#Is the end position of the selected region within this exon?
	if ($region_end >= $exon_start && $region_end <= $exon_end){
	  push (@exons_within_region, $exon_id);
	  next();
	}
      }
      $region_count++;
      $regions{$region_count}{region_start} = $region_start;
      $regions{$region_count}{region_end} = $region_end;
      $regions{$region_count}{exons} = \@exons_within_region;
      if ($verbose eq "yes"){
	print YELLOW, "\n\t\tRegion: $region_count ($region_start - $region_end) for exon: $exon_id covers exons: @exons_within_region", RESET;
	print LOG "\n\t\tRegion: $region_count ($region_start - $region_end) for exon: $exon_id covers exons: @exons_within_region";
      }
    }
  }

  return(\%regions);
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

  my $probe_type = "Exon";

  foreach my $probe (sort {$a<=>$b} keys %{$probe_object_ref}){
    $probe_number++;
    $current_probe_id++;

    my $Tm = $probe_object_ref->{$probe}->{Tm};
    my $Tm_rounded = sprintf("%.3f", $Tm);

    my @exon1_IDs = @{$probe_object_ref->{$probe}->{exon_ids}};

   print OUTFILE "$current_probe_id\t$current_probeset_id\t$probe_object_ref->{$probe}->{gene_id}\t$probe_object_ref->{$probe}->{sequence}\t$probe_object_ref->{$probe}->{length}\t$Tm_rounded\t$probe_type\t@exon1_IDs\t$probe_object_ref->{$probe}->{start_pos}\t$probe_object_ref->{$probe}->{end_pos}\tna\tna\tna\t$probe_object_ref->{$probe}->{masked_bases}\n";

  }
  return();
}
