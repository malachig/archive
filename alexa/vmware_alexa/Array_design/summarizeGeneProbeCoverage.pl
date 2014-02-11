#!/usr/bin/perl -w

#Written by Malachi Griffith
#This script will determine the number of probes that are theoretically possible for each targeted by the probes found in input probe files
#It will also determine the number of probes that have successfully been designed based on input files and output these values
#as a file contain two columns, and one row for every gene.  The two columns will contain the theoretical and actual number of 
#probes for each gene.
#This script is also capable of considering the exon-skipping parameter when calculating the theoretical number of probes expected as well
#as when it counts the probes actually designed

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
my $probe_dir = '';
my $skipped_exon_threshold = '';
my $junction_probeset_size = '';
my $exon_probeset_size = '';
my $target_length = '';
my $max_length_variance = '';
my $outfile = '';
my $logfile = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'probe_dir=s'=>\$probe_dir, 'skipped_exon_threshold=i'=>\$skipped_exon_threshold,
	    'junction_probeset_size=i'=>\$junction_probeset_size, 'exon_probeset_size=i'=>\$exon_probeset_size,
	    'target_length=i'=>\$target_length, 'max_length_variance=s'=>\$max_length_variance,
	    'outfile=s'=>\$outfile, 'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script determines the number of exons skipped by each exon-exon junction probe specified in an input file", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password\n", RESET;
print GREEN, "\n\tSpecify the full path to a directory containing probe files using: --probe_dir", RESET;
print GREEN, "\n\tSpecify the threshold for exons skipped to allow using: --skipped_exon_threshold (use 0 for no limit)", RESET;
print GREEN, "\n\tSpecify the desired number of probes per junction probeset (probes per junction) using: --junction_probeset_size", RESET;
print GREEN, "\n\tSpecify the desired number of probes per exon probeset (probes per exon) using: --exon_probeset_size", RESET;
print GREEN, "\n\tSpecify the target length used to originally generate exon probes using: --target_length (say 36)", RESET;
print GREEN, "\n\tSpecify the max length variance used to originally generate exon probes using: --max_length_variance (say 8)", RESET;
print GREEN, "\n\n\tSpecify the name of the resulting outfile using: --outfile", RESET;
print GREEN, "\n\tSpecify the name of an output logfile using: --logfile", RESET;
print GREEN, "\n\nExample: summarizeGeneProbeCoverage.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --probe_dir=/home/user/alexa/ALEXA_version/unfiltered_probes  --skipped_exon_threshold=3  --junction_probeset_size=2  --exon_probeset_size=3  --target_length=36  --max_length_variance=8  --outfile=/home/user/alexa/ALEXA_version/stats/gene_probe_coverage_beforeFilter.txt  --logfile=/home/user/alexa/ALEXA_version/logs/summarizeGeneProbeCoverage_beforeFilter_LOG.txt\n\n", RESET;

unless ($database && $server && $user && $password && $probe_dir && $junction_probeset_size && $exon_probeset_size && $target_length && $max_length_variance && $outfile && $logfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}
unless ($skipped_exon_threshold || $skipped_exon_threshold eq '0'){
  print RED, "\nOption --skipped_exon_threshold must be a positive integer or 0\n\n", RESET;
  exit();
}

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";
#Print out parameters selected
print LOG "\n\nUser Specified the following options:\ndatabase = $database\nprobe_dir = $probe_dir\nskipped_exon_threshold = $skipped_exon_threshold\njunction_probeset_size = $junction_probeset_size\nexon_probeset_size = $exon_probeset_size\ntarget_length = $target_length\nmax_length_variance = $max_length_variance\noutfile = $outfile\nlogfile = $logfile\n\n";

#1.) Go through all three input files and build an object containing all probes and probesets for each gene ID found
#    - Skip negative control and intron probes
my $gene_probes_ref = &parseProbeFiles('-probe_dir'=>$probe_dir, '-skipped_exon_threshold'=>$skipped_exon_threshold);


#2.) For each of these genes, get the number of exons and exon-clusters and add to the gene_probes object (keyed on gene id)
#    Also calculate the theoretical number of probes needed to target all exon regions
print BLUE, "\nGet theoretical number of probes for all of the exon regions of each gene\n\n", RESET;
print LOG "\nGet theoretical number of probes for all of the exon regions of each gene\n\n";

&createGeneObject('-dbh'=>$alexa_dbh, '-gene_probes_object'=>$gene_probes_ref);

#3.) Now, get the theoretical number of exon-exon and exon-intron junction probes that would ideally be designed for each gene
print BLUE, "\nGet theoretical number of exon junction and intron junction probes for each gene", RESET;
print LOG "\nGet theoretical number of exon junction and intron junction probes for each gene";

my @multi_exon_genes;
foreach my $gene_id (keys %{$gene_probes_ref}){
  unless ($gene_probes_ref->{$gene_id}->{exon_count} == 1){
    push (@multi_exon_genes, $gene_id);
  }
}
my $probe_counts_ref = &junctionProbeCombinations ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@multi_exon_genes, '-exon_skip_limit'=>$skipped_exon_threshold);

foreach my $gene_id (@multi_exon_genes){

  my $exon_exon = (($probe_counts_ref->{$gene_id}->{exon_exon}) * $junction_probeset_size);        #Multiply by desired probeset size
  my $intron_exon = (($probe_counts_ref->{$gene_id}->{intron_exon}) * $junction_probeset_size);    #Multiply by desired probeset size

  $gene_probes_ref->{$gene_id}->{theoretical_exon_junction_probes} = $exon_exon;
  $gene_probes_ref->{$gene_id}->{theoretical_intron_junction_probes} = $intron_exon;
}

#4.) Now go through each gene and count the probes of each type
#Only count probes if the desired number of probes per probeset is achieved
#Do not count exon-junction probes that exceed the skipped_exon threshold
#Intron probes and negative control probes are not considered in these calculations
foreach my $gene_id (sort {$a <=> $b} keys %{$gene_probes_ref}){

  my $probesets_ref = $gene_probes_ref->{$gene_id}->{probesets};

  foreach my $probeset_id (sort {$a <=> $b} keys %{$probesets_ref}){

    my $probe_type = $probesets_ref->{$probeset_id}->{probe_type};

    my $probes_ref = $probesets_ref->{$probeset_id}->{probes};

    #make sure the desired number of probes was found for this probeset
    my $probe_count = keys %{$probes_ref};


    if ($probe_type eq "Exon-Exon"){

      if ($probe_count >= $junction_probeset_size){
	$gene_probes_ref->{$gene_id}->{actual_exon_junction_probes} += $junction_probeset_size;
      }else{
	$gene_probes_ref->{$gene_id}->{actual_exon_junction_probes} += $probe_count;
      }

    }elsif (($probe_type eq "Exon-Intron") || ($probe_type eq "Intron-Exon")){
      if ($probe_count >= $junction_probeset_size){
	$gene_probes_ref->{$gene_id}->{actual_intron_junction_probes} += $junction_probeset_size;
      }else{
	$gene_probes_ref->{$gene_id}->{actual_intron_junction_probes} += $probe_count;
      }

    }elsif ($probe_type eq "Exon"){
      if ($probe_count >= $exon_probeset_size){
	$gene_probes_ref->{$gene_id}->{actual_exon_probes} += $exon_probeset_size;
      }else{
	$gene_probes_ref->{$gene_id}->{actual_exon_probes} += $probe_count;
      }
    }else{
      print RED, "\nProbe type: $probe_type not understood!!\n\n";
      close (LOG);
      $alexa_dbh->disconnect();
      exit();
    }
  }
}

#5.) Print an output file containing: gene_ID, ideal probe count, actual probe count
&printOutputFile('-gene_probes_object'=>$gene_probes_ref, '-output_file'=>$outfile);

close (LOG);

#Close database connection
$alexa_dbh->disconnect();

exit();


#######################################################################################
#Parse the input probe files
#######################################################################################
sub parseProbeFiles{
  my %args = @_;
  my $probe_dir = $args{'-probe_dir'};
  my $skipped_exon_threshold = $args{'-skipped_exon_threshold'};

  my %gene_probes;

  #1.) Get files from this directory
  my %files;

  unless ($probe_dir =~ /.*\/$/){
    $probe_dir = "$probe_dir"."/";
  }
  #First make sure the specified base path exists and is a directory
  unless (-e $probe_dir && -d $probe_dir){
    print RED, "\nSpecified directory: $probe_dir does not appear valid!\n\n", RESET;
    close (LOG);
    $alexa_dbh->disconnect();
    exit();
  }

  #Get all the input files in the repeat masked result directory
  print BLUE, "\nSearching $probe_dir for files\n", RESET;
  print LOG "\nSearching $probe_dir for files\n";

  my $file_count = 0;
  opendir(DIRHANDLE, "$probe_dir") || die "\nCannot open directory: $probe_dir\n\n";
  my @test_files = readdir(DIRHANDLE);

 FILE:foreach my $test_file (@test_files){
    my $file_path = "$probe_dir"."$test_file";

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      print BLUE, "\n\t$file_path  is a directory - skipping", RESET;
      print LOG "\n\t$file_path  is a directory - skipping";
      next();
    }

    #Only process the first line of each file for now
    my %columns;
    open (PROBE_FILE, "$file_path") || die "\nCould not open probe file: $file_path\n\n";
    while (<PROBE_FILE>){
      my $line_record = $_;
      chomp ($line_record);

      my @line = split ("\t", $line_record);

      #Watch for the header line which is assumed to contain neccessary columns
      my $column_count = 0;
      foreach my $column (@line){
	$columns{$column}{column_pos} = $column_count;
	$column_count++;
      }

      #If the file has the neccessary columns, add it to the list of files to be processed
      if ($columns{'Probe_Count'} && $columns{'ProbeSet_ID'} && $columns{'Gene_ID'} && $columns{'Probe_Type'}){
	$file_count++;
	$files{$file_count}{file_path} = $file_path;
	$files{$file_count}{file_name} = $test_file;
	$files{$file_count}{source_dir} = $probe_dir;
	$files{$file_count}{columns} = \%columns;

      }else{
	print YELLOW, "\nFile: $file_path does not appear to be a valid probe file (expected column missing) - skipping\n\n", RESET;
	print LOG "\nFile: $file_path does not appear to be a valid probe file (expected column missing) - skipping\n\n";
      }
      close (PROBE_FILE);
      next FILE;  #skip to next file after processing only the first line
    }
  }

  #2.) Summarize the files found - list old and new new names
  print BLUE, "\n\nFile summary:", RESET;
  print LOG "\n\nFile summary:";
  foreach my $file_count (sort {$a <=> $b} keys %files){
    print BLUE, "\n\tFile: $file_count\n\t\tINFILE: $files{$file_count}{file_path}", RESET;
    print LOG "\n\tFile: $file_count\n\t\tINFILE: $files{$file_count}{file_path}";
  }
  print BLUE, "\n\n", RESET;
  print LOG "\n\n";


  #3.) Go through each file and get probe data
  foreach my $file_count (sort {$a <=> $b} keys %files){

    print BLUE, "\nBegin processing the input probe file: $files{$file_count}{file_path}\n\n", RESET;
    print LOG "\nBegin processing the input probe file: $files{$file_count}{file_path}\n\n";

    open (PROBE_FILE, "$files{$file_count}{file_path}") || die "\nCould not open probe file: $files{$file_count}{file_path}\n\n";
    my $header = 1;
    while (<PROBE_FILE>){

      #skip the first line
      if ($header == 1){
	$header = 0;
	next();
      }

      my $line_record = $_;
      chomp ($line_record);
      my @line = split ("\t", $line_record);

      my %columns = %{$files{$file_count}{columns}};
      my $probe_id = $line[$columns{'Probe_Count'}{column_pos}];
      my $probeset_id = $line[$columns{'ProbeSet_ID'}{column_pos}];
      my $gene_id = $line[$columns{'Gene_ID'}{column_pos}];
      my $probe_type = $line[$columns{'Probe_Type'}{column_pos}];

      #If this probe is an Exon-Exon probe, deal with the exon-skipping values
      if ($probe_type eq "Exon-Exon"){
	#Deal with exon junction probes

	unless ($columns{'Exons_Skipped'}){
	  print RED, "\nExon junction probes can not be processed unless Exons_Skipped column is defined!!\n\n", RESET;
	  close (LOG);
	  $alexa_dbh->disconnect();
	  exit();
	}

	#Get the number of exons-skipped by the exon-exon junction of this probe
	my $exons_skipped = $line[$columns{'Exons_Skipped'}{column_pos}];

	unless ($exons_skipped =~ /^\d+/){
	  print RED, "\nFormat of exons_skipped value: $exons_skipped, not understood\n\n", RESET;
	  close (LOG);
	  $alexa_dbh->disconnect();
	  exit();
	}

	#If the number exons skipped by this exon exceeds the max, do not count this probe
	if ($skipped_exon_threshold == 0 || $exons_skipped <= $skipped_exon_threshold){

	  #Add probe to hash
	  if ($gene_probes{$gene_id}){
	    my $probesets_ref = $gene_probes{$gene_id}{probesets};

	    if ($probesets_ref->{$probeset_id}){
	      my $probes_ref = $probesets_ref->{$probeset_id}->{probes};
	      $probes_ref->{$probe_id}->{tmp} = '';

	    }else{
	      my %probes;
	      $probes{$probe_id}{tmp} = '';
	      $probesets_ref->{$probeset_id}->{probes} = \%probes;
	      $probesets_ref->{$probeset_id}->{probe_type} = $probe_type;
	    }
	  }else{
	    my %probesets;
	    my %probes;
	    $probes{$probe_id}{tmp} = '';
	    $probesets{$probeset_id}{probes} = \%probes;
	    $probesets{$probeset_id}{probe_type} = $probe_type;
	    $gene_probes{$gene_id}{probesets} = \%probesets;
	  }
	}

      }elsif ($probe_type eq "Control-Negative" || $probe_type eq "Intron"){
	#Skip negative control and intron probes
	next();

      }else{
	#Add all other probe types
	#Add probe to hash
	if ($gene_probes{$gene_id}){
	  my $probesets_ref = $gene_probes{$gene_id}{probesets};

	  if ($probesets_ref->{$probeset_id}){
	    my $probes_ref = $probesets_ref->{$probeset_id}->{probes};
	    $probes_ref->{$probe_id}->{tmp} = '';

	  }else{
	    my %probes;
	    $probes{$probe_id}{tmp} = '';
	    $probesets_ref->{$probeset_id}->{probes} = \%probes;
	    $probesets_ref->{$probeset_id}->{probe_type} = $probe_type;
	  }
	}else{
	  my %probesets;
	  my %probes;
	  $probes{$probe_id}{tmp} = '';
	  $probesets{$probeset_id}{probes} = \%probes;
	  $probesets{$probeset_id}{probe_type} = $probe_type;
	  $gene_probes{$gene_id}{probesets} = \%probesets;
	}

      }
    }
    close (PROBE_FILE);
  }

  return(\%gene_probes);
}


################################################################################################
#Create gene object - return hash with gene sequence, masked sequence, exon info, etc.         #
################################################################################################
sub createGeneObject{
  my %args = @_;
  my $alexa_dbh = $args{'-dbh'};
  my $gene_object_ref = $args{'-gene_probes_object'};

  my @gene_ids = keys %{$gene_object_ref};

  #Get all exons for these genes
  my $gene_exons_ref = &getExons ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

  foreach my $gene_id (keys %{$gene_object_ref}){

    my $exons_ref = $gene_exons_ref->{$gene_id}->{exons};

    my $unique_exons = keys %{$exons_ref};

    #Go through each exon and create clusters of overlapping exons.  Any exon that shares at least one bp with another exon will be placed in the same cluster
    #An exon that is completely isolated from others will comprise a cluster of one
    my $cluster_count = 0;
    my %clusters;
  EXON: foreach my $exon_id (sort {$exons_ref->{$a}->{exon_start} <=> $exons_ref->{$b}->{exon_start}} keys %{$exons_ref}){

      #Get the start and end positions of the current exon
      my $exon_start = $exons_ref->{$exon_id}->{exon_start};
      my $exon_end = $exons_ref->{$exon_id}->{exon_end};

      #Go through each cluster and see if the current exon overlaps with one of the exons already present
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


    #Now using the exon clusters identified here, determine the number of exon regions exactly as was done when generating the exon probes
    my $total_exon_region_count = 0;
    foreach my $cluster (sort {$a <=> $b} keys %clusters){
      my $cluster_exons_aref = $clusters{$cluster}{exons};
      my $cluster_size = @{$cluster_exons_aref};

      if ($cluster_size == 1){
	#If cluster consists of only a single exon - there is only one exon region
	$total_exon_region_count++;

      }else{
	#If the cluster consists of multiple exons, determine the number of regions to be targeted
	my $exon_regions_ref = &defineExonRegions('-exon_object'=>$exons_ref, '-exon_ids'=>$cluster_exons_aref);
	my $exon_region_count = keys %{$exon_regions_ref};
	$total_exon_region_count += $exon_region_count;
      }
    }

    #Add the info gathered above to the gene object
    my $final_cluster_count = keys %clusters;
    $gene_object_ref->{$gene_id}->{clusters} = \%clusters;
    $gene_object_ref->{$gene_id}->{cluster_count} = $final_cluster_count;
    $gene_object_ref->{$gene_id}->{exon_count} = $unique_exons;
    $gene_object_ref->{$gene_id}->{exon_region_count} = $total_exon_region_count;
    $gene_object_ref->{$gene_id}->{theoretical_exon_probes} = ($total_exon_region_count * $exon_probeset_size);

    #Initialize exon-exon and exon-intron theoretical probe counts
    $gene_object_ref->{$gene_id}->{theoretical_exon_junction_probes} = 0;
    $gene_object_ref->{$gene_id}->{theoretical_intron_junction_probes} = 0;

    #Initialize exon-exon, exon-intron, and exon actual probe counts
    $gene_object_ref->{$gene_id}->{actual_exon_junction_probes} = 0;
    $gene_object_ref->{$gene_id}->{actual_intron_junction_probes} = 0;
    $gene_object_ref->{$gene_id}->{actual_exon_probes} = 0;

  }
  return();
}


##################################################################################################
#Print gene IDs, probe counts and theoretical probe counts                                       #
##################################################################################################
sub printOutputFile{
  my %args = @_;
  my $gene_object_ref = $args{'-gene_probes_object'};
  my $output_file = $args{'-output_file'};

  open (OUTFILE, ">$output_file") || die "\nCould not open output file: $output_file\n\n";

  print OUTFILE "Gene_ID\tTheoreticalProbeCount\tActualProbeCount\n";

  my $grand_a_probes = 0;
  my $grand_a_ej_probes = 0;
  my $grand_a_ij_probes = 0;
  my $grand_a_e_probes = 0;
  my $grand_t_probes = 0;
  my $grand_t_ej_probes = 0;
  my $grand_t_ij_probes = 0;
  my $grand_t_e_probes = 0;


  foreach my $gene_id (sort {$a <=> $b} keys %{$gene_object_ref}){
    my $combined_theoretical_count = $gene_object_ref->{$gene_id}->{theoretical_exon_junction_probes} + $gene_object_ref->{$gene_id}->{theoretical_intron_junction_probes} + $gene_object_ref->{$gene_id}->{theoretical_exon_probes};

    my $combined_actual_count = $gene_object_ref->{$gene_id}->{actual_exon_junction_probes} + $gene_object_ref->{$gene_id}->{actual_intron_junction_probes} + $gene_object_ref->{$gene_id}->{actual_exon_probes};

    $grand_a_probes += $combined_actual_count;
    $grand_a_ej_probes += $gene_object_ref->{$gene_id}->{actual_exon_junction_probes};
    $grand_a_ij_probes += $gene_object_ref->{$gene_id}->{actual_intron_junction_probes};
    $grand_a_e_probes += $gene_object_ref->{$gene_id}->{actual_exon_probes};
    $grand_t_probes += $combined_theoretical_count;
    $grand_t_ej_probes += $gene_object_ref->{$gene_id}->{theoretical_exon_junction_probes};
    $grand_t_ij_probes += $gene_object_ref->{$gene_id}->{theoretical_intron_junction_probes};
    $grand_t_e_probes += $gene_object_ref->{$gene_id}->{theoretical_exon_probes};

    print BLUE, "\nGene: $gene_id\tExons: $gene_object_ref->{$gene_id}->{exon_count}\tCombinedProbes: ($combined_actual_count / $combined_theoretical_count)\tExonJunction: ($gene_object_ref->{$gene_id}->{actual_exon_junction_probes} / $gene_object_ref->{$gene_id}->{theoretical_exon_junction_probes})\tIntronJunction: ($gene_object_ref->{$gene_id}->{actual_intron_junction_probes} /  $gene_object_ref->{$gene_id}->{theoretical_intron_junction_probes})\tExon: ($gene_object_ref->{$gene_id}->{actual_exon_probes} / $gene_object_ref->{$gene_id}->{theoretical_exon_probes})", RESET;
    print LOG "\nGene: $gene_id\tExons: $gene_object_ref->{$gene_id}->{exon_count}\tCombinedProbes: ($combined_actual_count / $combined_theoretical_count)\tExonJunction: ($gene_object_ref->{$gene_id}->{actual_exon_junction_probes} / $gene_object_ref->{$gene_id}->{theoretical_exon_junction_probes})\tIntronJunction: ($gene_object_ref->{$gene_id}->{actual_intron_junction_probes} /  $gene_object_ref->{$gene_id}->{theoretical_intron_junction_probes})\tExon: ($gene_object_ref->{$gene_id}->{actual_exon_probes} / $gene_object_ref->{$gene_id}->{theoretical_exon_probes})";

    print OUTFILE "$gene_id\t$combined_theoretical_count\t$combined_actual_count\n";
  }
  close OUTFILE;

  print BLUE, "\n\nGRAND TOTALS:", RESET;
  print BLUE, "\nExonJunction Probes: $grand_a_ej_probes / $grand_t_ej_probes", RESET;
  print BLUE, "\nIntronJunction Probes: $grand_a_ij_probes / $grand_t_ij_probes", RESET;
  print BLUE, "\nExon Probes: $grand_a_e_probes / $grand_t_e_probes", RESET;
  print BLUE, "\nCombined Probes: $grand_a_probes / $grand_t_probes\n\n", RESET;

  print LOG "\n\nGRAND TOTALS:", RESET;
  print LOG "\nExonJunction Probes: $grand_a_ej_probes / $grand_t_ej_probes";
  print LOG "\nIntronJunction Probes: $grand_a_ij_probes / $grand_t_ij_probes";
  print LOG "\nExon Probes: $grand_a_e_probes / $grand_t_e_probes";
  print LOG "\nCombined Probes: $grand_a_probes / $grand_t_probes\n\n";

  return();
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
      if ($region_start <= $exon_start && $region_end >= $exon_end){
	push (@exons_within_region, $exon_id);
	next();
      }

    }
    $region_count++;
    $regions{$region_count}{region_start} = $region_start;
    $regions{$region_count}{region_end} = $region_end;
    $regions{$region_count}{exons} = \@exons_within_region;
    #print BLUE, "\n\t\tRegion: $region_count ($region_start - $region_end) covers exons: @exons_within_region", RESET;
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

      #Still need a probe region for this exon!
      my $region_start = ($exon_object_ref->{$exon_id}->{exon_start})+1;
      my $region_end = ($exon_object_ref->{$exon_id}->{exon_end})-1;

      #Make sure this region fits the required probe length
      my $region_size = $region_end - $region_start;
      if ($region_size <= ($target_length + $max_length_variance)){
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
    }
  }

  return(\%regions);
}
