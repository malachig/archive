#!/usr/bin/perl -w
#Written by Malachi Griffith

#The purpose of this script is to summarize DE or SI events detected at the probeset level (exons, boundaries, junctions, etc.)
#Stats will be generated for each gene but information for each event will also be gathered and printed to a new file

#EVENT LEVEL SUMMARY
#1.) Import list of all probesets
#    - For example:
#    - ~/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/Probeset_List.txt

#2.) For each of these events get info for the corresponding probesets.  Get this info from probe files stored here:
#    - /home/malachig/AlternativeSplicing/Array_design/ALEXA_hs_35_35h/probe_dev/filtered_probes/
#    - Go through these files and grab info for probes corresponding to the probesets imported in (1)
#    - Store:
#            - start and end coordinate for entire probeset (max and min coord observed for any probe)
#            - max probe_length of probes in this probeset
#            - max 'enst_largestTargetHitLength'.  This tells us whether the event is 'Known' or 'Novel'

#3.) Identify features for each probeset.
#    A.) Does the probeset overlap the ORF (or is contained within it), or is it entirely limited to the 5' UTR or 3' UTR
#        - To determine this, first get the CDS start and end for ALL transcripts for the target gene.
#        - Then identify the MIN CDS start and MAX CDS end coords for the gene.  Use these to classify the probeset according to its coords
#    B.) Does the probeset overlap any transmembrane domains? How many?
#    C.) Does the probeset overlap any signal peptides? How many?
#    D.) Does the probeset overlap any Coiled-coil domains? How many?
#    E.) Does the probeset overlap any protein domains?  How many?  Store the complete list of their names!
#    F.) Does the probeset correspond to a known transcript sequence?  'Known' vs 'Novel'
#        - Consider it to be known if the max EnsEMBL transcript hit from the target gene is 95% or more of the length of the probe.

#4.) Print out an appended file. Same an input file except it has the info gathered in (2) and (3) added


use strict;
use Data::Dumper;
use Getopt::Long;
use Math::Complex;
use Term::ANSIColor qw(:constants);

use lib '/home/malachig/AlternativeSplicing/perl_bin';
use Statistics::Descriptive;
use utilities::utility qw(:all);
use utilities::ALEXA_DB_35_35h qw(:all);

#1.) Import lists of significant DE and SI events
my %all_probesets_list;
my %all_genes_list;

my $all_EventsFile = "/home/malachig/AlternativeSplicing/Array_analysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/Probeset_List.txt";

#Create a hash of these events, keyed on probeset ID
my $all_events_ref = &parse_EventsFile('-all_file'=>$all_EventsFile);

#2.) For each of these events get info for the corresponding probesets.
my $probe_dir = "/home/malachig/AlternativeSplicing/Array_design/ALEXA_hs_35_35h/probe_dev/filtered_probes/";
my @required_columns = qw(Probe_Count ProbeSet_ID Gene_ID Probe_length Probe_Tm Probe_Type Unit1_start Unit1_end Unit2_start Unit2_end enst_largestTargetHitLength mRNA_largestTargetHitLength EST_largestTargetHitLength);
my %files = %{&getProbeFiles('-probe_dir'=>$probe_dir)};

my %probesets;
foreach my $file_count (sort {$files{$a}->{min_probe_id} <=> $files{$b}->{min_probe_id}} keys %files){

  print BLUE, "\nProcessing file: $files{$file_count}{file_path}\n\n", RESET;
  my %columns = %{$files{$file_count}{columns}};
  &importProbeData('-input_file'=>$files{$file_count}{file_path}, '-columns'=>\%columns);
}

#3.) Identify features for each probeset.
my $alexa_dbh = &connectDB('-database'=>'ALEXA_hs_35_35h', '-server'=>'jango.bcgsc.ca', '-user'=>'malachig', '-password'=>'gEEnom$');
my @gene_ids = keys %all_genes_list;

print BLUE, "\nGetting gene transcript and orf info\n\n", RESET;
my $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

#3-A.) Does the probeset overlap the ORF (or is contained within it), or is it entirely limited to the 5' UTR or 3' UTR
print BLUE, "\nIdentifying the position of each event with respect to the CDS\n\n", RESET;
&identifyEventPositions();

#3-B-E.) Does the probeset overlap any protein features? How many of each type (TMDs, SPs, CCs, and protein domains)?
print BLUE, "\nGetting protein features for all genes with significant events\n\n", RESET;
my $gene_protein_features_ref = &getProteinFeatures ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);

my $gene_count1 = @gene_ids;
my $gene_count2 = keys %{$gene_transcripts_ref};
my $gene_count3 = keys %{$gene_protein_features_ref};

&identifyOverlapingProteinFeatures();

#3-F.) Does the probeset correspond to a known transcript sequence?  'Known' vs 'Novel'
&identifyNovelKnownProbesets();

#4.) Print out an appended file. Same an input file except it has the info gathered in (2) and (3) added
my $out_file = "/home/malachig/AlternativeSplicing/Array_analysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/Probeset_List_ANNOTATED.txt";

printOutputFiles('-out_file'=>$out_file, '-out_file'=>$out_file);

$alexa_dbh->disconnect();

exit();


#########################################################################################################
#1.) Import lists of significant DE and SI events
#########################################################################################################
sub parse_EventsFile{
  my %args = @_;
  my $all_file = $args{'-all_file'};

  my %all_events;
  print BLUE, "\n\n1-A.) Processing input file: $all_file\n", RESET;

  #ProbeSet_ID, AlexaGene_ID, Probe_Type, Exons_Skipped
  open (ALL, "$all_file") || die "\nCould not open file: $all_file\n\n";

  my $line_count = 0;
  while(<ALL>){
    chomp($_);
    my @line = split("\t", $_);

    unless ($line[0] =~ /^\d+/){
      next();
    }
    $line_count++;
    my $probeset_id = $line[0];

    $all_probesets_list{$probeset_id}{gene_id} = $line[1];
    $all_probesets_list{$probeset_id}{tmd_count} = 0;
    $all_probesets_list{$probeset_id}{sp_count} = 0;
    $all_probesets_list{$probeset_id}{cc_count} = 0;
    $all_probesets_list{$probeset_id}{domain_count} = 0;
    my @domain_names;
    $all_probesets_list{$probeset_id}{domain_names} = \@domain_names;
    $all_probesets_list{$probeset_id}{novel_known} = 'na';
    $all_probesets_list{$probeset_id}{est_support} = 'na';

    unless ($line[1] eq "0"){
      $all_genes_list{$line[1]}{tmp} = '';
    }
    $all_events{$probeset_id}{gene_id} = $line[1];
    $all_events{$probeset_id}{probe_type} = $line[3];
    $all_events{$probeset_id}{exons_skipped} = $line[4];
    $all_events{$probeset_id}{line_count} = $line_count;
  }
  close (ALL);

  my $all_events_count = keys %all_events;

  print BLUE, "\nImported $all_events_count events\n\n", RESET;

  return(\%all_events);
}



###########################################################################################################
#Get probe files and order according to probeset ranges                                                   #
###########################################################################################################
sub getProbeFiles{
  my %args = @_;
  my $probe_dir = $args{'-probe_dir'};

  unless ($probe_dir =~ /.*\/$/){
    $probe_dir = "$probe_dir"."/";
  }

  #Make sure the specified directory is valid
  unless (-e $probe_dir && -d $probe_dir){
    print RED, "\nSpecified directory does not appear to be valid: $probe_dir\n\n", RESET;
    exit();
  }

  #Get files from this directory
  print BLUE, "\nSearching $probe_dir for probe files", RESET;

  my %possible_files;
  my $possible_file_count = 0;
  opendir(DIRHANDLE, "$probe_dir") || die "\nCannot open directory: $probe_dir\n\n";
  my @test_files = readdir(DIRHANDLE);

  foreach my $test_file (@test_files){
    my $file_path = "$probe_dir"."$test_file";

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      print BLUE, "\n\t$file_path  is a directory - skipping", RESET;
      next();
    }
    $possible_file_count++;

    $possible_files{$possible_file_count}{file_name} = $test_file;
    $possible_files{$possible_file_count}{source_dir} = $probe_dir;
    $possible_files{$possible_file_count}{file_path} = $file_path;
  }

  my $file_num = keys %possible_files;
  print BLUE, "\n\nFound $file_num probe files in the specified directory\n", RESET;

  #Check each file for expected columns, get probe id ranges to allow ordering of files
  my %files;
 FILE: foreach my $file_count (sort {$a <=> $b} keys %possible_files){
    my $probe_file = $possible_files{$file_count}{file_path};

    print BLUE, "\n\nExamining probe file: $probe_file\n", RESET;

    open (PROBE_FILE, "$probe_file") || die "\nCould not open probe file: $probe_file\n\n";

    my %columns;
    my $max_probe_id = 0;
    my $min_probe_id = 0;
    my $max_probeset_id = 0;
    my $min_probeset_id = 0;
    my $header = 1;
    my $first_data = 1;

    while (<PROBE_FILE>){
      my $line_record = $_;
      chomp ($line_record);

      my @line = split ("\t", $line_record);

      #Watch for the header line
      if ($header == 1){
	$header = 0;
	my $column_count = 0;
	foreach my $column (@line){
	  $columns{$column}{column_pos} = $column_count;
	  $column_count++;
	}

	#If the file has the neccessary columns, add it to the list of files to be processed
	my $required_columns_found = 1;
	my @missing_columns;

	foreach my $req_col (@required_columns){
	  unless ($columns{$req_col}){
	    $required_columns_found = 0;
	    push(@missing_columns, $req_col);
	  }
	}

	if ($required_columns_found == 1){

	  $files{$file_count}{file_path} = $probe_file;
	  $files{$file_count}{file_name} = $possible_files{$file_count}{file_name};
	  $files{$file_count}{source_dir} = $possible_files{$file_count}{source_dir};
	  $files{$file_count}{columns} = \%columns;

	}else{
	  print RED, "\nFile: $possible_files{$file_count}{file_name} does not appear to be a complete probe file.", RESET;
	  print RED, "\n\tColumns missing: @missing_columns", RESET;
	  print RED, "\n\tComplete appropriate tests before attempting to import probe records\n\n", RESET;
	  exit();
	}
	next();
      }
      #Process data lines
      if ($first_data == 1){

	#Skip files with negative control probes
	#if ($line[$columns{'Probe_Type'}{column_pos}] eq "Control-Negative"){
	#  print YELLOW, "\nFile: $possible_files{$file_count}{file_name} appears to contain negative control probes - skipping\n\n", RESET;
	#  my $test = delete ($files{$file_count});
	#  next FILE;
	#}

	$max_probe_id = $line[0];
	$min_probe_id = $line[0];
	$first_data = 0;
	next();
      }
      if ($line[0] > $max_probe_id){$max_probe_id = $line[0];}
      if ($line[0] < $min_probe_id){$min_probe_id = $line[0];}
    }
    close (PROBE_FILE);

    $files{$file_count}{min_probe_id} = $min_probe_id;
    $files{$file_count}{max_probe_id} = $max_probe_id;
  }

  #Summarize the files found
  print BLUE, "\n\nFile summary:", RESET;
  foreach my $file_count (sort {$files{$a}->{min_probe_id} <=> $files{$b}->{min_probe_id}} keys %files){
    print BLUE, "\n\tFile: $file_count\n\t\tINFILE: $files{$file_count}{file_path} (probe_ids: $files{$file_count}{min_probe_id} - $files{$file_count}{max_probe_id})", RESET;
  }
  print BLUE, "\n\n", RESET;

  return(\%files);
}


###########################################################################################################
#Import Probe Data from an input probe file
###########################################################################################################
sub importProbeData{
  my %args = @_;
  my $probe_file = $args{'-input_file'};
  my %columns = %{$args{'-columns'}};

  #Open the probe file and read the neccessary probe data into a hash keyed on probe ID
  open (PROBES, "$probe_file") || die "\nCould not open input probe file: $probe_file";

  my $first_line = 1;

  while (<PROBES>){

    #Skip the header line
    if ($first_line == 1){
      $first_line = 0;
      next();
    }

    #Get the values of interest from each line (probe record)
    chomp($_);
    my @probe_line = split ("\t", $_);

    #Probe_Count
    my $probe_id = $probe_line[$columns{Probe_Count}{column_pos}];
    unless ($probe_id =~ /^\d+/){
      print RED, "\nInvalid probe ID: $probe_id\n\n", RESET;
      exit();
    }

    #ProbeSet_ID
    my $probeset_id = $probe_line[$columns{ProbeSet_ID}{column_pos}];

    #Unless this probeset is in the list of significant events, skip it
    unless ($all_probesets_list{$probeset_id}){
      next();
    }

    my $gene_id = $probe_line[$columns{Gene_ID}{column_pos}];
    my $probe_length = $probe_line[$columns{Probe_length}{column_pos}];
    my $probe_type = $probe_line[$columns{Probe_Type}{column_pos}];
    my $unit1_start = $probe_line[$columns{Unit1_start}{column_pos}];
    my $unit1_end = $probe_line[$columns{Unit1_end}{column_pos}];
    my $unit2_start = $probe_line[$columns{Unit2_start}{column_pos}];
    my $unit2_end = $probe_line[$columns{Unit2_end}{column_pos}];
    my $enst_target_hit_length = $probe_line[$columns{'enst_largestTargetHitLength'}{column_pos}];
    my $mrna_target_hit_length = $probe_line[$columns{'mRNA_largestTargetHitLength'}{column_pos}];
    my $est_target_hit_length = $probe_line[$columns{'EST_largestTargetHitLength'}{column_pos}];

    if ($probesets{$probeset_id}){

      if ($probe_length > $probesets{$probeset_id}{max_probe_length}){
	$probesets{$probeset_id}{max_probe_length} = $probe_length;
      }

      if ($probe_type eq "Control-Negative"){
	next();
      }

      if ($unit1_start < $probesets{$probeset_id}{unit1_start}){
	$probesets{$probeset_id}{unit1_start} = $unit1_start;
      }
      if ($unit1_end > $probesets{$probeset_id}{unit1_end}){
	$probesets{$probeset_id}{unit1_end} = $unit1_end;
      }
      unless ($unit2_start eq "na" || $unit2_end eq "na"){
	if ($unit2_start < $probesets{$probeset_id}{unit2_start}){
	  $probesets{$probeset_id}{unit2_start} = $unit2_start;
	}
	if ($unit2_end > $probesets{$probeset_id}{unit2_end}){
	  $probesets{$probeset_id}{unit2_end} = $unit2_end;
	}
      }
      if ($enst_target_hit_length > $probesets{$probeset_id}{largest_enst_hit}){
	$probesets{$probeset_id}{largest_enst_hit} = $enst_target_hit_length;
      }
      if ($mrna_target_hit_length > $probesets{$probeset_id}{largest_mrna_hit}){
	$probesets{$probeset_id}{largest_mrna_hit} = $mrna_target_hit_length;
      }
      if ($est_target_hit_length > $probesets{$probeset_id}{largest_est_hit}){
	$probesets{$probeset_id}{largest_est_hit} = $est_target_hit_length;
      }

    }else{
      $probesets{$probeset_id}{gene_id} = $gene_id;
      $probesets{$probeset_id}{max_probe_length} = $probe_length;
      $probesets{$probeset_id}{probe_type} = $probe_type;
      $probesets{$probeset_id}{unit1_start} = $unit1_start;
      $probesets{$probeset_id}{unit1_end} = $unit1_end;
      $probesets{$probeset_id}{unit2_start} = $unit2_start;
      $probesets{$probeset_id}{unit2_end} = $unit2_end;
      $probesets{$probeset_id}{largest_enst_hit} = $enst_target_hit_length;
      $probesets{$probeset_id}{largest_mrna_hit} = $mrna_target_hit_length;
      $probesets{$probeset_id}{largest_est_hit} = $est_target_hit_length;
    }

  }
  close (PROBES);

  return ();
}

############################################################################################################################
#3.) Identify features for each probeset.
############################################################################################################################

############################################################################################################################
#3-A.) Does the probeset overlap the ORF (or is contained within it), or is it entirely limited to the 5' UTR or 3' UTR
############################################################################################################################
sub identifyEventPositions{

  foreach my $probeset_id (keys %all_probesets_list){

    my $gene_id = $all_probesets_list{$probeset_id}{gene_id};
    my $probe_type = $probesets{$probeset_id}{probe_type};

    $all_probesets_list{$probeset_id}{position} = 'na';

    if ($probe_type eq "Control-Negative"){
      next();
    }

    my $grand_cds_start = $gene_transcripts_ref->{$gene_id}->{grand_cds_start};
    my $grand_cds_end = $gene_transcripts_ref->{$gene_id}->{grand_cds_end};

    my $test_start = $probesets{$probeset_id}{unit1_start};
    my $test_end = $probesets{$probeset_id}{unit1_end};

    if ($probesets{$probeset_id}{unit2_end} =~ /\d+/){
      $test_end = $probesets{$probeset_id}{unit2_end};
    }


    if ($test_end < $grand_cds_start){
      $all_probesets_list{$probeset_id}{position} = "5prime_UTR";
    }
    if ($test_start > $grand_cds_end){
      $all_probesets_list{$probeset_id}{position} = "3prime_UTR";
    }
    if (($test_end >= $grand_cds_start && $test_end <= $grand_cds_end) || ($test_start >= $grand_cds_start && $test_start <= $grand_cds_end)){
      $all_probesets_list{$probeset_id}{position} = "ORF";
    }

  }
  return();
}


#############################################################################################################################
#3-B-E.) Does the probeset overlap any protein features? How many of each type (TMDs, SPs, CCs, and protein domains)?
#############################################################################################################################
sub identifyOverlapingProteinFeatures{

  foreach my $probeset_id (keys %all_probesets_list){
    my $gene_id = $all_probesets_list{$probeset_id}{gene_id};

    #Watch out for cases where no protein features were found for a gene
    unless($gene_protein_features_ref->{$gene_id}){
      print YELLOW, "\nNo protein features found for this gene: $gene_id", RESET;
      next();
    }

    #Get a list of non-redundant protein features for this gene
    my @nr_pfs = @{$gene_protein_features_ref->{$gene_id}->{nr_pfs}};

    #Get the transcripts and protein features for the gene this probeset belongs to
    my $trans_ref = $gene_protein_features_ref->{$gene_id}->{transcripts};

    #Go through each transcript
    foreach my $trans_id (keys %{$trans_ref}){

      my $pfs_ref = $trans_ref->{$trans_id}->{protein_features};

      #Go through each NR protein feature and see if the current probset overlaps with it
      foreach my $pf_id (@nr_pfs){

	unless ($pfs_ref->{$pf_id}){next();}  #Non-Redundant PFs were created at gene level but this causes problem as we cycle through transcripts

	my $type = $pfs_ref->{$pf_id}->{type};
	my $name = $pfs_ref->{$pf_id}->{name};

	my $combo_name = "$pf_id"."_"."$type"."_"."$name";
	my @start_coords = @{$pfs_ref->{$pf_id}->{start_coords}};
	my @end_coords = @{$pfs_ref->{$pf_id}->{end_coords}};

	my $unit1_start = $probesets{$probeset_id}{unit1_start};
	my $unit1_end = $probesets{$probeset_id}{unit1_end};
	my $unit2_start = $probesets{$probeset_id}{unit2_start};
	my $unit2_end = $probesets{$probeset_id}{unit2_end};

	my $overlap_found = 0;

	#Go through the coordinate blocks of the protein feature and see if the probeset coords overlap
	foreach my $pf_start (@start_coords){
	  my $pf_end = shift (@end_coords);

	  if (($unit1_start > $pf_start && $unit1_start < $pf_end) || ($unit1_end > $pf_start && $unit1_end < $pf_end)){
	    $overlap_found = 1;
	  }
	  #probe completely flanks domain
	  if (($unit1_start < $pf_start) && ($unit1_end > $pf_end)){
	    $overlap_found = 1;
	  }

	  if ($unit2_start =~ /\d+/ && $unit2_end =~ /\d+/){
	    if (($unit2_start > $pf_start && $unit2_start < $pf_end) || ($unit2_end > $pf_start && $unit2_end < $pf_end)){
	      $overlap_found = 1;
	    }
	    #probe completely flanks domain
	    if (($unit2_start < $pf_start) && ($unit2_end > $pf_end)){
	      $overlap_found = 1;
	    }
	  }
	}

	#If overlap between the probeset and this protein feature was found, make note of it in the appropriate category
	#Type Categories: TMDs='tmhmm', SPs='Signalp', CCs='ncoils', protein_domains = 'Pfam','pfscan','Prints','scanprosite'
	if ($overlap_found == 1){
	  if ($type eq "tmhmm"){
	    $all_probesets_list{$probeset_id}{tmd_count}++;
	  }
	  if ($type eq "Signalp"){
	    $all_probesets_list{$probeset_id}{sp_count}++;
	  }
	  if ($type eq "ncoils"){
	    $all_probesets_list{$probeset_id}{cc_count}++;
	  }
	  if ($type eq "Pfam" || $type eq "pfscan" || $type eq "Prints" || $type eq "scanprosite"){
	    $all_probesets_list{$probeset_id}{domain_count}++;
	    push(@{$all_probesets_list{$probeset_id}{domain_names}}, $combo_name);
	  }
	}
      }

    }

  }
  return();
}


#################################################################################################################################
#3-F.) Does the probeset correspond to a known transcript sequence?  'Known' vs 'Novel'.  Does it have EST support?
#################################################################################################################################
sub identifyNovelKnownProbesets{

  foreach my $probeset_id (keys %all_probesets_list){

    my $max_probe_length = $probesets{$probeset_id}{max_probe_length};
    my $enst_target_hit_length = $probesets{$probeset_id}{largest_enst_hit};
    my $mrna_target_hit_length = $probesets{$probeset_id}{largest_mrna_hit};
    my $est_target_hit_length = $probesets{$probeset_id}{largest_est_hit};

    my $enst_match_percent = ($enst_target_hit_length / $max_probe_length)*100;
    my $mrna_match_percent = ($mrna_target_hit_length / $max_probe_length)*100;
    my $est_match_percent = ($est_target_hit_length / $max_probe_length)*100;

    $all_probesets_list{$probeset_id}{novel_known} = "novel";
    $all_probesets_list{$probeset_id}{est_support} = "NO";

    if ($enst_match_percent >= 95){
      $all_probesets_list{$probeset_id}{novel_known} = "known";
    }
    if ($mrna_match_percent >= 95){
      $all_probesets_list{$probeset_id}{novel_known} = "known";
    }
    if ($est_match_percent >= 95){
      $all_probesets_list{$probeset_id}{est_support} = "YES";
    }
  }

  return();
}


#################################################################################################################################
#4.) Print out an appended file. Same an input file except it has the info gathered in (2) and (3) added
#################################################################################################################################
sub printOutputFiles{
  my %args = @_;
  my $out_file = $args{'-out_file'};

  #DE OUT FILE
  open (OUT, ">$out_file") || die "\nCould not open output file: $out_file\n\n";

  print OUT "ProbeSet_ID\tAlexaGene_ID\tProbe_Type\tExons_Skipped\tProbePosition\tTMD_count\tSignalPeptide_count\tCoiledCoil_count\tProteinDomainCount\tNovel_or_Known\tEST_Support\tProteinDomainNames\n";

  foreach my $probeset_id (sort {$all_events_ref->{$a}->{line_count} <=> $all_events_ref->{$b}->{line_count}} keys %{$all_events_ref}){
    print OUT "$probeset_id\t$all_events_ref->{$probeset_id}->{gene_id}\t$all_events_ref->{$probeset_id}->{probe_type}\t$all_events_ref->{$probeset_id}->{exons_skipped}\t$all_probesets_list{$probeset_id}{position}\t$all_probesets_list{$probeset_id}{tmd_count}\t$all_probesets_list{$probeset_id}{sp_count}\t$all_probesets_list{$probeset_id}{cc_count}\t$all_probesets_list{$probeset_id}{domain_count}\t$all_probesets_list{$probeset_id}{novel_known}\t$all_probesets_list{$probeset_id}{est_support}\t@{$all_probesets_list{$probeset_id}{domain_names}}\n";
  }
  close (OUT);

  return();
}

