#!/usr/bin/perl -w
#Written by Malachi Griffith
#This script takes a tab-delimited input file Array design file and creates UCSC tracks to represent the position of selected probes
#The purpose of this is to assist a manual validation of the design

#Process:
#1.) Parse the filtered probes files in the specified directory and get coordinates and gene info for these selected probes
#2.) Use ALEXA to get basic gene info
#3.) Create UCSC tracks to display the probe positions
#4.) Create a summary file with gene_id and a link to the custom track file

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

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $filtered_probe_dir = '';
my $ucsc_dir = '';
my $web_path = '';
my $ucsc_build = '';
my $outfile = '';
my $gzip_path = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'filtered_probe_dir=s'=>\$filtered_probe_dir, 'ucsc_dir=s'=>\$ucsc_dir, 'web_path=s'=>\$web_path, 'ucsc_build=s'=>\$ucsc_build,
	    'outfile=s'=>\$outfile, 'gzip_path=s'=>\$gzip_path);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the directory containing filtered probe files using: --filtered_probe_dir", RESET;
print GREEN, "\n\tSpecify the name of the target ucsc track file directory using: --ucsc_dir", RESET;
print GREEN, "\n\t\tThis should be a directory that is accesible to the internet (i.e. served by Apache)", RESET;
print GREEN, "\n\tSpecify the web path for this directory using: --web_path", RESET;
print GREEN, "\n\tSpecify the name of the UCSC genome build corresponding to this array design using: --ucsc_build", RESET;
print GREEN, "\n\tSpecify the name of a summary gene outfile using: --outfile", RESET;
print GREEN, "\n\tSpecify the full path to your gzip utility using: --gzip_path", RESET;
print GREEN, "\n\nExample: createProbeUcscTracks.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=username  --password=pwd  --filtered_probe_dir=/home/user/alexa/ALEXA_version/filtered_probes/  --ucsc_dir=/home/user/www/public/htdocs/ALEXA_version/  --web_path='http://www.bcgsc.ca/people/user/htdocs/alexa_tracks/ALEXA_version/'  --ucsc_build=hg17  --outfile=/home/user/alexa/ALEXA_version/gene_probe_summary.txt  --gzip_path=/usr/bin/gzip\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $filtered_probe_dir && $ucsc_dir && $web_path && $ucsc_build && $outfile && $gzip_path){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

unless ($ucsc_dir =~ /.*\/$/){
  $ucsc_dir = "$ucsc_dir"."/";
}
unless ($web_path =~ /.*\/$/){
  $web_path = "$web_path"."/";
}

#Make sure the ucsc_dir provided exists and is a directory
unless (-e $ucsc_dir && -d $ucsc_dir){
  print RED, "\nSpecified directory: $ucsc_dir does not appear valid!\n\n", RESET;
  exit();
}

#Specify desired standard display tracks depending on database specified
my $desired_tracks = "mrna";
if ($database eq "ALEXA_dm_42_43"){
  $desired_tracks = "mrna flyBaseGene flyBaseNoncoding";
}
if ($database eq "ALEXA_dr_42_6c" || $database eq "ALEXA_hs_41_36c" || $database eq "ALEXA_hs_35_35h" || $database eq "ALEXA_mm_41_36b" || $database eq "ALEXA_rn_43_34m"){
  $desired_tracks = "ensGene knownGene";
}
if ($database eq "ALEXA_sc_41_1d"){
  $desired_tracks = "mrna sgdGene";
}
if($database eq "ALEXA_ce_45_170b"){
  $desired_tracks = "sangerGene sangerRnaGene";
}

my %track_names;
$track_names{'Exon-Exon'}{name} = "ExonJunct";
$track_names{'Exon-Exon'}{description} = "Exon-Exon Junction Probes";
$track_names{'Exon-Exon'}{color} = "0,0,255"; #blue

$track_names{'Intron-Junction'}{name} = "IntronJunct";
$track_names{'Intron-Junction'}{description} = "Intron-Exon Junction Probes";
$track_names{'Intron-Junction'}{color} = "255,0,0"; #red

$track_names{'Exon'}{name} = "Exon";
$track_names{'Exon'}{description} = "Exon Probes";
$track_names{'Exon'}{color} = "0,100,0"; #yellow

$track_names{'Intron'}{name} = "Intron";
$track_names{'Intron'}{description} = "Intron Probes";
$track_names{'Intron'}{color} = "47,79,79"; #dark slate

#1.) Parse the filtered probes files in the specified directory and get coordinates and gene info for each probe
my $gene_probes_ref = &importFilteredProbes('-filtered_probe_dir'=>$filtered_probe_dir);

#2.) Use ALEXA to get basic gene info

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

my $gene_info_ref;
&getBasicGeneInfo('-dbh'=>$alexa_dbh, '-gene_probes_object'=>$gene_probes_ref);

#Close database connection
$alexa_dbh->disconnect();

#3.) Create UCSC tracks to display the probe positions
&generateUCSC_tracks('-ucsc_dir'=>$ucsc_dir, '-gene_probes_object'=>$gene_probes_ref);

#4.) Create a summary file with gene_id and a link to the custom track file
&printSummaryFile('-outfile'=>$outfile, '-gene_probes_object'=>$gene_probes_ref);

#Go to the target dir and compress the track files
my $command = "$gzip_path". " $ucsc_dir"."*";
print BLUE, "\nCompressing custom track files with command: $command\n\n", RESET;
system ("$command");

exit();


###################################################################################################################
#1.)Import the probes from the filtered exon and junction probes files.                                           #
###################################################################################################################
sub importFilteredProbes{
  my %args = @_;
  my $filtered_probe_dir = $args{'-filtered_probe_dir'};
  my $design_probes_ref = $args{'-design_probes_object'};

  unless ($filtered_probe_dir =~ /.*\/$/){
    $filtered_probe_dir = "$filtered_probe_dir"."/";
  }

  #First make sure the specified base path exists and is a directory
  unless (-e $filtered_probe_dir && -d $filtered_probe_dir){
    print RED, "\nSpecified directory: $filtered_probe_dir does not appear valid!\n\n", RESET;
    exit();
  }

  print BLUE, "\nProcessing input probe files from $filtered_probe_dir to get probe info for each probe from the design file\n\n", RESET;

  opendir(DIRHANDLE, "$filtered_probe_dir") || die "\nCannot open directory: $filtered_probe_dir\n\n";
  my @test_files = readdir(DIRHANDLE);
  closedir(DIRHANDLE);

  my @files;
  foreach my $file (@test_files){

    my $file_path = "$filtered_probe_dir"."$file";

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      print YELLOW, "\n\t$file_path  is a directory - skipping", RESET;;
      next();
    }
    push (@files, $file);
  }

  print BLUE, "\n\nFound the following filtered probe files:\n\t@files\n", RESET;

  my $grand_probe_count = 0;
  my %gene_probes;
  my %columns;

  foreach my $probe_file (@files){

    chomp($probe_file);
    my $probe_file_path = "$filtered_probe_dir"."$probe_file";

    open (PROBE, "$probe_file_path") || die "\nCould not open probe file: $probe_file_path\n\n";

    print BLUE, "\n\tCounting probes from: $probe_file", RESET;
    my $header_line;
    my $first_line = 1;

    while (<PROBE>){
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
	unless ($columns{'Probe_Count'} && $columns{'Gene_ID'} && $columns{'Probe_Type'} && $columns{'Unit1_start_chr'} && $columns{'Unit1_end_chr'} && $columns{'Unit2_start_chr'} && $columns{'Unit2_end_chr'}){
	  print RED, "\nCritical column missing or named incorrectly in probe input file: $probe_file, check input file", RESET;
	  exit();
	}
	next();
      }

      #Get the values of interest from each line (gene record)
      chomp($_);
      my @data_line = split ("\t", $_);

      my $probe_id = $data_line[$columns{'Probe_Count'}{column_pos}];
      my $gene_id = $data_line[$columns{'Gene_ID'}{column_pos}];
      my $probe_type = $data_line[$columns{'Probe_Type'}{column_pos}];
      my $unit1_start = $data_line[$columns{'Unit1_start_chr'}{column_pos}];
      my $unit1_end = $data_line[$columns{'Unit1_end_chr'}{column_pos}];
      my $unit2_start = $data_line[$columns{'Unit2_start_chr'}{column_pos}];
      my $unit2_end = $data_line[$columns{'Unit2_end_chr'}{column_pos}];

      #Skip negative control probes which do not even map to the genome
      if ($probe_type eq "Control-Negative"){
	next();
      }

      if ($probe_type eq "Exon-Intron" || $probe_type eq "Intron-Exon"){
	$probe_type = "Intron-Junction";
      }

      $grand_probe_count++;

      if ($gene_probes{$gene_id}){
	my $probes_ref = $gene_probes{$gene_id}{probes};
	$probes_ref->{$probe_id}->{probe_type} = $probe_type;
	$probes_ref->{$probe_id}->{unit1_start} = $unit1_start;
	$probes_ref->{$probe_id}->{unit1_end} = $unit1_end;
	$probes_ref->{$probe_id}->{unit2_start} = $unit2_start;
	$probes_ref->{$probe_id}->{unit2_end} = $unit2_end;
	$gene_probes{$gene_id}{probe_count}++;

      }else{
	my %probes;

	$probes{$probe_id}{probe_type} = $probe_type;
	$probes{$probe_id}{unit1_start} = $unit1_start;
	$probes{$probe_id}{unit1_end} = $unit1_end;
	$probes{$probe_id}{unit2_start} = $unit2_start;
	$probes{$probe_id}{unit2_end} = $unit2_end;
	$gene_probes{$gene_id}{probe_count} = 1;
	$gene_probes{$gene_id}{probes} = \%probes;
      }
    }
    close (PROBE);
  }

  my $genes_with_probes = keys %gene_probes;

  print BLUE, "\n\nFound a total of $grand_probe_count probes corresponding to $genes_with_probes gene targets\n\n", RESET;

  return(\%gene_probes);
}


#######################################################################################################################
#2.)Get basic info for each gene                                                                                      #
#######################################################################################################################
sub getBasicGeneInfo{
  my %args = @_;

  my $dbh = $args{'-dbh'};
  my $gene_probes_ref = $args{'-gene_probes_object'};

  print BLUE, "\nGetting basic gene info for each gene\n\n", RESET;

  my @gene_ids = keys %{$gene_probes_ref};
  $gene_info_ref = &getGeneInfo ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");
  my $gene_count = 0;

  my $gene_transcripts_ref = &getTranscripts('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

  my $type = "EntrezGene";
  my $gene_ids_ref = &getGeneTerms ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids, '-id_type'=>$type);

  foreach my $gene_id (sort {$a <=> $b} keys %{$gene_probes_ref}){
    $gene_count++;

    my $chr_corrected = "chr"."$gene_info_ref->{$gene_id}->{chromosome}";
    $gene_probes_ref->{$gene_id}->{chromosome} = $chr_corrected;

    my $strand = $gene_info_ref->{$gene_id}->{chr_strand};
    if ($strand eq "1"){
      $gene_probes_ref->{$gene_id}->{chr_strand} = "+";
    }else{
      $gene_probes_ref->{$gene_id}->{chr_strand} = "-";
    }

    #Add the gene info
    $gene_probes_ref->{$gene_id}->{chr_start} = $gene_info_ref->{$gene_id}->{chr_start};
    $gene_probes_ref->{$gene_id}->{chr_end} = $gene_info_ref->{$gene_id}->{chr_end};

    my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};
    my $trans_count = keys %{$transcripts_ref};
    $gene_info_ref->{$gene_id}->{trans_count} = $trans_count;

    my %nr_exons;
    foreach my $trans_id (keys %{$transcripts_ref}){
      my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

      foreach my $exon_id (keys %{$exons_ref}){
	my $start_stop = "$exons_ref->{$exon_id}->{exon_start}"."."."$exons_ref->{$exon_id}->{exon_end}";
	$nr_exons{$start_stop}{tmp} = 'na';
      }
    }
    my $exon_count = keys %nr_exons;
    $gene_info_ref->{$gene_id}->{exon_count} = $exon_count;

    my $external_ids_ref = $gene_ids_ref->{$gene_id}->{external_ids};
    if ($external_ids_ref){
      my @entrez_ids = keys %{$external_ids_ref};
      $gene_info_ref->{$gene_id}->{entrez_ids} = \@entrez_ids;
    }else{
      my @entrez_ids = ('na');
      #push (@entrez_ids, 'na');
      $gene_info_ref->{$gene_id}->{entrez_ids} = \@entrez_ids;
    }
  }
  return();
}


#################################################################################################################################
#3.) Generate custom track files for the UCSC browser                                                                           #
#################################################################################################################################
sub generateUCSC_tracks{
  my %args = @_;
  my $ucsc_dir = $args{'-ucsc_dir'};
  my $gene_probes_ref = $args{'-gene_probes_object'};

  print BLUE "\nOrganizing genes by chromosome and creating custom tracks object\n\n", RESET;

  #First organize genes by chromosome
  my %gene_map;
  foreach my $gene_id (sort {$gene_probes_ref->{$a}->{chromosome} cmp $gene_probes_ref->{$b}->{chromosome}} keys %{$gene_probes_ref}){
    my $chromosome = $gene_probes_ref->{$gene_id}->{chromosome};

    if ($gene_map{$chromosome}){
      my $genes_ref = $gene_map{$chromosome}{genes};
      $genes_ref->{$gene_id}->{chr_start} = $gene_probes_ref->{$gene_id}->{chr_start};

    }else{
      my %genes;
      $genes{$gene_id}{chr_start} = $gene_probes_ref->{$gene_id}->{chr_start};
      $gene_map{$chromosome}{genes} = \%genes;
    }
  }

  #Now go through each chromosome and build a tracks object organize by chromosome, gene_block, gene ID, probe type, then each probe within that
  my %tracks;
  foreach my $chromosome (sort keys %gene_map){
    my $genes_ref = $gene_map{$chromosome}{genes};

    my $block = 1;
    my $gene_count = 0;

    foreach my $gene_id (sort {$genes_ref->{$a}->{chr_start} <=> $genes_ref->{$b}->{chr_start}} keys %{$genes_ref}){
      $gene_count++;
      if ($gene_count == 100){
	$block++;
	$gene_count = 0;
      }

      #Make note of the block this gene was put into and the ultimate file name:
      my $filename = "$chromosome"."_"."$block".".txt";
      my $filepath = "$ucsc_dir"."$filename";
      my $weblink = "$web_path"."$filename".".gz";

      $gene_probes_ref->{$gene_id}->{filepath} = $filepath;
      $gene_probes_ref->{$gene_id}->{weblink} = $weblink;

      my $strand = $gene_probes_ref->{$gene_id}->{chr_strand};

      #Get the probes for this gene
      my $probes_ref = $gene_probes_ref->{$gene_id}->{probes};

      foreach my $probe_id (sort {$a <=> $b} keys %{$probes_ref}){
	my $probe_type = $probes_ref->{$probe_id}->{probe_type};
	my $unit1_start_genomic = $probes_ref->{$probe_id}->{unit1_start};
	my $unit1_end_genomic = $probes_ref->{$probe_id}->{unit1_end};
	my $unit2_start_genomic = $probes_ref->{$probe_id}->{unit2_start};
	my $unit2_end_genomic = $probes_ref->{$probe_id}->{unit2_end};

	#create a track line for this probe
	my $feature_name = "$probe_id";
	my $track_line_unit1 = "$chromosome\tALEXA\tProbe\t$unit1_start_genomic\t$unit1_end_genomic\t.\t$strand\t.\t$feature_name";

	#Create a second track line if the unit2 coords are defined for this probe
	my $track_line_unit2;
	if ($unit2_start_genomic =~ /\d+/ && $unit2_end_genomic =~ /\d+/){
	  $track_line_unit2 = "$chromosome\tALEXA\tProbe\t$unit2_start_genomic\t$unit2_end_genomic\t.\t$strand\t.\t$feature_name";
	}

	if ($tracks{$chromosome}){
	  my $gene_blocks_ref = $tracks{$chromosome}{gene_blocks};

	  if ($gene_blocks_ref->{$block}){
	    my $probe_types_ref = $gene_blocks_ref->{$block}->{probe_types};

	    if ($probe_types_ref->{$probe_type}){
	      my $track_lines_ref = $probe_types_ref->{$probe_type}->{track_lines};
	      push (@{$track_lines_ref}, $track_line_unit1);
	      if ($track_line_unit2){
		push (@{$track_lines_ref}, $track_line_unit2);
	      }

	    }else{
	      #First entry for this probe type but there are some other probe types with entries already
	      my @track_lines;
	      push (@track_lines, $track_line_unit1);
	      if ($track_line_unit2){
		push (@track_lines, $track_line_unit2);
	      }
	      $probe_types_ref->{$probe_type}->{track_lines} = \@track_lines;
	    }

	  }else{
	    #First entry for this block but there are some other blocks for this chromosome already
	    my %probe_types;
	    my @track_lines;
	    push (@track_lines, $track_line_unit1);
	    if ($track_line_unit2){
	      push (@track_lines, $track_line_unit2);
	    }
	    $probe_types{$probe_type}{track_lines} = \@track_lines;
	    $gene_blocks_ref->{$block}->{probe_types} = \%probe_types;
	    $gene_blocks_ref->{$block}->{filepath} = $filepath;
	  }

	}else{
	  #First entry for this chromosome
	  my %gene_blocks;
	  my %probe_types;
	  my @track_lines;

	  push (@track_lines, $track_line_unit1);
	  if ($track_line_unit2){
	    push (@track_lines, $track_line_unit2);
	  }

	  $probe_types{$probe_type}{track_lines} = \@track_lines;
	  $gene_blocks{$block}{probe_types} = \%probe_types;
	  $gene_blocks{$block}{filepath} = $filepath;
	  $tracks{$chromosome}{gene_blocks} = \%gene_blocks;
	}
      }
    }
  }

  #Go through each track print it out
  foreach my $chromosome (sort keys %tracks){

    my $gene_blocks_ref = $tracks{$chromosome}{gene_blocks};

    foreach my $block (sort {$a <=> $b} keys %{$gene_blocks_ref}){

      my $filepath = $gene_blocks_ref->{$block}->{filepath};

      open (UCSC, ">$filepath") || die "\nCould not open ucsc file: $filepath\n\n";
      #Set default browser settings
      print UCSC "#Browser line";
      #print UCSC "\nbrowser position chr$genes{$gene_id}{chromosome}:$genes{$gene_id}{chr_start}-$genes{$gene_id}{chr_end}";  #Default start pos
      print UCSC "\nbrowser hide all";
      print UCSC "\nbrowser full $desired_tracks";

      my $probe_types_ref = $gene_blocks_ref->{$block}->{probe_types};

      foreach my $probe_type (sort keys %{$probe_types_ref}){

	#Print each track out
	print UCSC "\n\n#Probe types: $probe_type";
	print UCSC "\ntrack name=$track_names{$probe_type}{name} description=\"$track_names{$probe_type}{description}\" color=$track_names{$probe_type}{color} useScore=0 visibility=3";

	#Data lines: seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup
	print UCSC "\n#BEGIN DATA - $probe_type";
	my $track_lines_ref = $probe_types_ref->{$probe_type}->{track_lines};

	foreach my $track_line (@{$track_lines_ref}){
	  print UCSC "\n$track_line";
	}
      }
      close (UCSC);
    }

  }

  return();
}


########################################################################################################################
#4.) Create a summary file with gene_id and a link to the custom track file                                            #
########################################################################################################################
sub printSummaryFile{
  my %args = @_;
  my $outfile = $args{'-outfile'};
  my $gene_probes_ref = $args{'-gene_probes_object'};

  print BLUE, "\nPrint summary output file: $outfile\n\n", RESET;

  open (SUMMARY, ">$outfile") || die "\nCould not open summary file: $outfile\n\n";

  #alexa_gene_id, ensembl_gene_id, gene_name, entrez_id, transcript count, exon count, chromosome, strand, start, end, probe_count, ucsc link
  print SUMMARY "ALEXA_ID\tEnsEMBL_Gene_ID\tGeneName\tEntrez_ID\tTranscript_Count\tExon_Count\tChromosome\tStrand\tStart\tEnd\tProbeCount\tucsc_link\n";

  foreach my $gene_id (sort {$a <=> $b} keys %{$gene_probes_ref}){

    my $chr_start = $gene_probes_ref->{$gene_id}->{chr_start};
    my $chr_end = $gene_probes_ref->{$gene_id}->{chr_end};
    my $chromosome = $gene_probes_ref->{$gene_id}->{chromosome};
    my $chr_strand = $gene_probes_ref->{$gene_id}->{chr_strand};
    my $weblink = $gene_probes_ref->{$gene_id}->{weblink};

    my $link = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=$ucsc_build&position=$chromosome:$chr_start-$chr_end&hgt.customText=$weblink&ctfile_$ucsc_build=";

    print SUMMARY "$gene_id\t$gene_info_ref->{$gene_id}->{ensembl_g_id}\t$gene_info_ref->{$gene_id}->{gene_name}\t@{$gene_info_ref->{$gene_id}->{entrez_ids}}\t$gene_info_ref->{$gene_id}->{trans_count}\t$gene_info_ref->{$gene_id}->{exon_count}\t$gene_probes_ref->{$gene_id}->{chromosome}\t$gene_probes_ref->{$gene_id}->{chr_strand}\t$gene_probes_ref->{$gene_id}->{chr_start}\t$gene_probes_ref->{$gene_id}->{chr_end}\t$gene_probes_ref->{$gene_id}->{probe_count}\t$link\n";
    #print BLUE, "$gene_id\t$link\n", RESET;

  }

  close (SUMMARY);

  return();
}
