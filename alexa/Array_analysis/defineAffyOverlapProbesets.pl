#!/usr/bin/perl -w

#Written by Malachi Griffith
#The purpose of this script is to grab all of the raw probe values from the Affy platform for genes targeted on the ALEXA design
#These are gathered regardless of whether the particular exon or intron was covered in the ALEXA platform
#These are therefore different from AffyAlexa (AA) probesets which correspond to only exons that are covered in BOTH platforms

#NOTE: This entire analysis is based on build hg17

#1.) First get the list of ALEXA genes of interest from a design file containing info on all ALEXA probes
#    - At this time also grab the strand and chromosome for the genes targeted by these probes from the appropriate ALEXA database
#    /home/malachig/AlternativeSplicing/Array_design/ALEXA_hs_35_35h/probe_dev/NimbleGen_Design_Submission/NimbleGenArray_mipVS5FUR_coords
#
#2.) The Affy exon array probesets which map to within an EnsEMBL GENE are already known and defined in a mapfile
#    The Affy exon array probesets which map to within an EnsEMBL EXON are also already known and defined in a mapfile
#    - For those genes identified above, get all the relevant probesets
#    - Create a hash of AFFY_PROBES(genes -> probesets -> coordinates)
#    - Probesets in the CORE file are 'Exon' probesets and the remainder that were in the COMPLETE file are 'Intron' probesets
#    /home/malachig/AlternativeSplicing/perl_bin/Affy_analysis/mapfiles/hg17/COMPLETE_probesets_mapped_to_ensembl.txt
#    /home/malachig/AlternativeSplicing/perl_bin/Affy_analysis/mapfiles/hg17/CORE_probesets_mapped_to_ensembl.txt
#
#3.) Print out a file summarizing each of the AO probesets identified
#    - Print: Affy_probeset_ID, chr, start, end, strand, Affy_probe_count, probe_type
#
#4.) Create datafiles
#    - For each platform, import a raw datafile containing triplicate MIP vs 5FUR data
#    - Print out each line to a new data file only if it belongs to an AO probeset
#    - Keep the original probe ID but print out the new AO probeset ID for each probe data line
#Affy_data:
#/home/malachig/AlternativeSplicing/Array_analysis/Affymetrix_Exon_Arrays/MIP_vs_5FUR/MIP_vs_5FUR_8-26-2006/probelevel/HuEX-1_001-006_RawProbe.txt

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

use lib '/home/malachig/AlternativeSplicing/perl_bin/';
use utilities::utility qw(:all);
use utilities::ALEXA_DB_35_35h qw(:all);

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $design_file = '';
my $core_map_file = '';
my $complete_map_file = '';
my $affy_data_file = '';
my $summary_file = '';
my $AO_affy_data_file = '';

GetOptions ('database=s'=>\$database, 'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'design_file=s'=>\$design_file, 'core_map_file=s'=>\$core_map_file, 'complete_map_file=s'=>\$complete_map_file,
	    'affy_data_file=s'=>\$affy_data_file, 'summary_file=s'=>\$summary_file,
	    'AO_affy_data_file=s'=>\$AO_affy_data_file);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\nExample: defineAffyOverlapProbesets.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --design_file=/home/malachig/AlternativeSplicing/Array_design/ALEXA_hs_35_35h/probe_dev/NimbleGen_Design_Submission/NimbleGenArray_mipVS5FUR_coords  --core_map_file=/home/malachig/AlternativeSplicing/perl_bin/Affy_analysis/mapfiles/hg17/CORE_probesets_mapped_to_ensembl.txt  --complete_map_file=/home/malachig/AlternativeSplicing/perl_bin/Affy_analysis/mapfiles/hg17/COMPLETE_probesets_mapped_to_ensembl.txt  --affy_data_file=/home/malachig/AlternativeSplicing/Array_analysis/Affymetrix_Exon_Arrays/MIP_vs_5FUR/MIP_vs_5FUR_8-26-2006/probelevel/HuEX-1_001-006_RawProbe.txt  --summary_file=AO_Probesets_Summary.txt  --AO_affy_data_file=AO_affy_data.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $design_file && $core_map_file && $complete_map_file && $affy_data_file && $summary_file && $AO_affy_data_file){
  print RED, "\nOptions missing!\n\n", RESET;
  exit();
}

#1.) First get the list of ALEXA genes of interest from a design file containing info on all ALEXA probes
#    - At this time also grab the strand and chromosome for the genes targeted by these probes from the appropriate ALEXA database
#    /home/malachig/AlternativeSplicing/Array_design/ALEXA_hs_35_35h/probe_dev/NimbleGen_Design_Submission/NimbleGenArray_mipVS5FUR_coords
my %alexa_genes;

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

&getAlexaGenes('-design_file'=>$design_file, '-dbh'=>$alexa_dbh);

#Close database connection
$alexa_dbh->disconnect();

#print Dumper %alexa_genes;


#2.) The Affy exon array probesets which map to within an EnsEMBL GENE are already known and defined in a mapfile
#    The Affy exon array probesets which map to within an EnsEMBL EXON are also already known and defined in a mapfile
#    - For those genes identified above, get all the relevant probesets
#    - Create a hash of AFFY_PROBES(genes -> probesets -> coordinates)
#    - Probesets in the CORE file are 'Exon' probesets and the remainder that were in the COMPLETE file are 'Intron' probesets
#    /home/malachig/AlternativeSplicing/perl_bin/Affy_analysis/mapfiles/hg17/COMPLETE_probesets_mapped_to_ensembl.txt
#    /home/malachig/AlternativeSplicing/perl_bin/Affy_analysis/mapfiles/hg17/CORE_probesets_mapped_to_ensembl.txt

my %affy_probesets;
my %affy_probeset_list;
&getAffyProbesets('-core_map_file'=>$core_map_file, '-complete_map_file'=>$complete_map_file);

#3.) Print out a file summarizing each of the AO probesets identified
#    - Print: Affy_probeset_ID, chr, start, end, strand, Affy_probe_count, probe_type
#
&summarize_AO_probesets('-summary_file'=>$summary_file);

#4.) Create datafiles
#    - For each platform, import a raw datafile containing triplicate MIP vs 5FUR data
#    - Print out each line to a new data file only if it belongs to an AO probeset
#    - Keep the original probe ID but print out the new AO probeset ID for each probe data line
#Affy_data:
#/home/malachig/AlternativeSplicing/Array_analysis/Affymetrix_Exon_Arrays/MIP_vs_5FUR/MIP_vs_5FUR_8-26-2006/probelevel/HuEX-1_001-006_RawProbe.txt
&createAO_datafile('-affy_data_file'=>$affy_data_file, '-AO_affy_data_file'=>$AO_affy_data_file);


exit();


###################################################################################################################
#1.) First get the list of ALEXA genes of interest from a design file containing info on all ALEXA probes         #
###################################################################################################################
sub getAlexaGenes{
  my %args = @_;
  my $design_file = $args{'-design_file'};
  my $dbh = $args{'-dbh'};

  open (DESIGN, "$design_file") || die "\nCould not open design file: $design_file\n\n";

  my $first_line = 1;
  my %columns;

  while(<DESIGN>){
    chomp($_);
    my @line = split("\t", $_);

    #Watch for the header line
    if ($first_line == 1){
      $first_line = 0;
      my $column_count = 0;
      foreach my $column (@line){
	$columns{$column}{column_pos} = $column_count;
	$column_count++;
      }

      #Check for required columns
      unless($columns{'Probe_Count'} && $columns{'Gene_ID'}){
	print RED, "\nExpected column missing\n\n", RESET;
	exit();
      }
      next();
    }

    my $gene_id = $line[$columns{'Gene_ID'}{column_pos}];
    my $probe_id = $line[$columns{'Probe_Count'}{column_pos}];

    $alexa_genes{$gene_id}{chromosome} = '';
  }
  close (DESIGN);

  #Now, for all genes gathered, get additional neccessary gene info from the alexa database
  my @gene_ids = keys %alexa_genes;

  my $gene_info_ref = &getGeneInfo ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids,'-sequence'=>"no");

  foreach my $gene_id (keys %alexa_genes){

    my $temp_chr = $gene_info_ref->{$gene_id}->{chromosome};
    $alexa_genes{$gene_id}{chromosome} = "chr"."$temp_chr";

    my $strand = $gene_info_ref->{$gene_id}->{chr_strand};
    if ($strand eq "-1"){
      $strand = "-";
    }else{
      $strand = "+";
    }

    $alexa_genes{$gene_id}{strand} = $strand;
  }

  my $gene_count = keys %alexa_genes;
  print BLUE "\nFound $gene_count genes in the ALEXA design\n\n", RESET;

  return();
}


###################################################################################################################
#2.) The Affy exon array probesets which map to within an EnsEMBL gene are already known and defined in a mapfile #
###################################################################################################################
sub getAffyProbesets{
  my %args = @_;
  my $core_map_file = $args{'-core_map_file'};
  my $complete_map_file = $args{'-complete_map_file'};

  my $exon_probesets = 0;
  my $intron_probesets = 0;

  #1.) First identify the 'Exon' (CORE) probesets which match those genes targeted by the ALEXA design 
  open (CORE_MAP, "$core_map_file") || die "\nCould not open core map file: $core_map_file\n\n";

  my $first_line = 1;
  my %columns;

  while(<CORE_MAP>){
    chomp($_);
    my @line = split("\t", $_);

    #Watch for the header line
    if ($first_line == 1){
      $first_line = 0;
      my $column_count = 0;
      foreach my $column (@line){
	$columns{$column}{column_pos} = $column_count;
	$column_count++;
      }

      #Check for required columns
      unless($columns{'probe_set_id'} && $columns{'alexa_gene_id'} && $columns{'chromosome'} && $columns{'strand'} && $columns{'start'} && $columns{'stop'} && $columns{'probe_count'} && $columns{'ensembl_gene_id'}){
	print RED, "\nExpected column missing\n\n", RESET;
	exit();
      }
      next();
    }

    my $gene_id = $line[$columns{'alexa_gene_id'}{column_pos}];

    #Unless this probeset belongs to one of the genes targeted on the ALEXA design, skip it.
    unless ($alexa_genes{$gene_id}){
      next();
    }
    $exon_probesets++;

    my $probeset_id = $line[$columns{'probe_set_id'}{column_pos}];
    my $ensembl_g_id = $line[$columns{'ensembl_gene_id'}{column_pos}];
    my $chromosome = $line[$columns{'chromosome'}{column_pos}];
    my $strand = $line[$columns{'strand'}{column_pos}];
    my $start = $line[$columns{'start'}{column_pos}];
    my $end = $line[$columns{'stop'}{column_pos}];
    my $probe_count = $line[$columns{'probe_count'}{column_pos}];
    my $type = "Exon";

    $affy_probeset_list{$probeset_id}{gene_id} = $gene_id;
    $affy_probeset_list{$probeset_id}{type} = $type;
    $affy_probeset_list{$probeset_id}{ensembl_g_id} = $ensembl_g_id;

    if ($affy_probesets{$gene_id}){
      #If this gene was previously observed 
      my $probesets_ref = $affy_probesets{$gene_id}{probesets};

      $affy_probesets{$gene_id}{exon_probeset_count}++;  #Number of exon probes for this gene

      if ($probesets_ref->{$probeset_id}){
	#If this probeset was previously observed
	print RED, "\n\tDuplicate probeset!!", RESET;
      }else{
	#First time this probeset has been observed
	$probesets_ref->{$probeset_id}->{ensembl_g_id} = $ensembl_g_id;
	$probesets_ref->{$probeset_id}->{start} = $start;
	$probesets_ref->{$probeset_id}->{end} = $end;
	$probesets_ref->{$probeset_id}->{type} = $type;
	$probesets_ref->{$probeset_id}->{probe_count} = $probe_count;
	$probesets_ref->{$probeset_id}->{chromosome} = $chromosome;
	$probesets_ref->{$probeset_id}->{strand} = $strand;
      }
    }else{
      #First time gene was observed
      my %probesets;
      $probesets{$probeset_id}{ensembl_g_id} = $ensembl_g_id;
      $probesets{$probeset_id}{start} = $start;
      $probesets{$probeset_id}{end} = $end;
      $probesets{$probeset_id}{type} = $type;
      $probesets{$probeset_id}{probe_count} = $probe_count;
      $probesets{$probeset_id}{chromosome} = $chromosome;
      $probesets{$probeset_id}{strand} = $strand;

      $affy_probesets{$gene_id}{probesets} = \%probesets;
      $affy_probesets{$gene_id}{exon_probeset_count} = 1;  #Number of probes that will be used to estimate gene expression
    }
  }
  my $exon_probeset_count = keys %affy_probeset_list;
  my $gene_count = keys %affy_probesets;

  print BLUE "\nFound $gene_count overlaping genes and $exon_probeset_count Exon probesets in the AFFY design\n\n", RESET;

  close (CORE_MAP);


  #2.) Now identify the 'Intron' probesets which match those genes targeted by the ALEXA design 
  open (COMPLETE_MAP, "$complete_map_file") || die "\nCould not open complete map file: $complete_map_file\n\n";

  $first_line = 1;
  %columns = ();

  while(<COMPLETE_MAP>){
    chomp($_);
    my @line = split("\t", $_);

    #Watch for the header line
    if ($first_line == 1){
      $first_line = 0;
      my $column_count = 0;
      foreach my $column (@line){
	$columns{$column}{column_pos} = $column_count;
	$column_count++;
      }

      #Check for required columns
      unless($columns{'probe_set_id'} && $columns{'alexa_gene_id'} && $columns{'chromosome'} && $columns{'strand'} && $columns{'start'} && $columns{'stop'} && $columns{'probe_count'} && $columns{'ensembl_gene_id'}){
	print RED, "\nExpected column missing\n\n", RESET;
	exit();
      }
      next();
    }

    my $gene_id = $line[$columns{'alexa_gene_id'}{column_pos}];

    #Unless this probeset belongs to one of the genes targeted on the ALEXA design, skip it.
    unless ($alexa_genes{$gene_id}){
      next();
    }

    #If this probeset has already been identified as a CORE probeset, skip it
    my $probeset_id = $line[$columns{'probe_set_id'}{column_pos}];
    if ($affy_probeset_list{$probeset_id}){
      next();
    }

    $intron_probesets++;
    my $ensembl_g_id = $line[$columns{'ensembl_gene_id'}{column_pos}];
    my $chromosome = $line[$columns{'chromosome'}{column_pos}];
    my $strand = $line[$columns{'strand'}{column_pos}];
    my $start = $line[$columns{'start'}{column_pos}];
    my $end = $line[$columns{'stop'}{column_pos}];
    my $probe_count = $line[$columns{'probe_count'}{column_pos}];
    my $type = "Intron";

    $affy_probeset_list{$probeset_id}{gene_id} = $gene_id;
    $affy_probeset_list{$probeset_id}{type} = $type;
    $affy_probeset_list{$probeset_id}{ensembl_g_id} = $ensembl_g_id;

    if ($affy_probesets{$gene_id}){
      #If this gene was previously observed 
      my $probesets_ref = $affy_probesets{$gene_id}{probesets};

      $affy_probesets{$gene_id}{intron_probeset_count}++;  #Number of exon probes for this gene

      if ($probesets_ref->{$probeset_id}){
	#If this probeset was previously observed
	print RED, "\n\tDuplicate probeset!!", RESET;
      }else{
	#First time this probeset has been observed
	$probesets_ref->{$probeset_id}->{ensembl_g_id} = $ensembl_g_id;
	$probesets_ref->{$probeset_id}->{start} = $start;
	$probesets_ref->{$probeset_id}->{end} = $end;
	$probesets_ref->{$probeset_id}->{type} = $type;
	$probesets_ref->{$probeset_id}->{probe_count} = $probe_count;
	$probesets_ref->{$probeset_id}->{chromosome} = $chromosome;
	$probesets_ref->{$probeset_id}->{strand} = $strand;
      }
    }else{
      #First time gene was observed
      my %probesets;
      $probesets{$probeset_id}{ensembl_g_id} = $ensembl_g_id;
      $probesets{$probeset_id}{start} = $start;
      $probesets{$probeset_id}{end} = $end;
      $probesets{$probeset_id}{type} = $type;
      $probesets{$probeset_id}{probe_count} = $probe_count;
      $probesets{$probeset_id}{chromosome} = $chromosome;
      $probesets{$probeset_id}{strand} = $strand;

      $affy_probesets{$gene_id}{probesets} = \%probesets;
      $affy_probesets{$gene_id}{exon_probeset_count} = 1;  #Number of probes that will be used to estimate gene expression
    }
  }
  my $total_probeset_count = keys %affy_probeset_list;
  my $intron_probeset_count = $total_probeset_count - $exon_probeset_count;
  $gene_count = keys %affy_probesets;

  print BLUE "\nFound $gene_count overlaping genes and $intron_probeset_count Intron probesets in the AFFY design\n\n", RESET;

  close (COMPLETE_MAP);

  return();
}


#############################################################################################################
#3.) Print out a file summarizing each of the AO probesets identified
#############################################################################################################
sub summarize_AO_probesets{
  my %args = @_;
  my $summary_file = $args{'-summary_file'};

  open (SUMMARY, ">$summary_file") || die "\nCould not open summary output file: $summary_file\n\n";

  print SUMMARY "probe_set_ID\talexa_gene_id\tchromosome\tstart\tend\tstrand\tprobe_count\tProbe_Type\n";

  foreach my $gene_id (sort {$a <=> $b} keys %affy_probesets){

    my $probesets_ref = $affy_probesets{$gene_id}{probesets};

    foreach my $probeset_id (sort {$a <=> $b} keys %{$probesets_ref}){

      print SUMMARY "$probeset_id\t$gene_id\t$probesets_ref->{$probeset_id}->{chromosome}\t$probesets_ref->{$probeset_id}->{start}\t$probesets_ref->{$probeset_id}->{end}\t$probesets_ref->{$probeset_id}->{strand}\t$probesets_ref->{$probeset_id}->{probe_count}\t$probesets_ref->{$probeset_id}->{type}\n";

    }
  }

  close (SUMMARY);

  return();
}


####################################################################################################################################################
#4.) Create datafiles
#    - For each platform, import a raw datafile containing triplicate MIP vs 5FUR data
#    - Print out each line to a new data file only if it belongs to an AO probeset
#    - Keep the original probe ID but print out the new AO probeset ID for each probe data line
#Affy_data:
#/home/malachig/AlternativeSplicing/Array_analysis/Affymetrix_Exon_Arrays/MIP_vs_5FUR/MIP_vs_5FUR_8-26-2006/probelevel/HuEX-1_001-006_RawProbe.txt
###################################################################################################################################################
sub createAO_datafile{
  my %args = @_;
  my $affy_data_file = $args{'-affy_data_file'};
  my $AO_affy_data_file = $args{'-AO_affy_data_file'};

  #Affy data file
  print BLUE, "\nProcessing Affy Data IN/OUT files\n\n", RESET;
  open (AFFY_IN, "$affy_data_file") || die "\nCould not open affy_data_file: $affy_data_file\n\n";
  open (AFFY_OUT, ">$AO_affy_data_file") || die "\nCould not open output affy_data_file: $AO_affy_data_file\n\n";

  my $affy_header = 1;
  while (<AFFY_IN>){
    chomp($_);

    #Skip comment lines
    if ($_ =~ /^\#/){
      next();
    }

    if ($affy_header == 1){
      print AFFY_OUT "AlexaGene_ID\tEnsEMBL_Gene_ID\tProbe_Type\t$_\n";
      $affy_header = 0;
      next();
    }

    my @line = split ("\t", $_);
    my $probe_id = $line[0];
    my $probeset_id = $line[1];

    unless ($probe_id =~ /\d+/){
      next();
    }

    #Unless this probe belongs to one of the probeset IDs identified above, skip it
    unless ($affy_probeset_list{$probeset_id}){
      next();
    }

    my $gene_id = $affy_probeset_list{$probeset_id}{gene_id};
    my $probe_type = $affy_probeset_list{$probeset_id}{type};
    my $ensembl_g_id = $affy_probeset_list{$probeset_id}{ensembl_g_id};

    print AFFY_OUT "$gene_id\t$ensembl_g_id\t$probe_type\t$_\n";

  }

  close (AFFY_IN);
  close (AFFY_OUT);

  return();
}
