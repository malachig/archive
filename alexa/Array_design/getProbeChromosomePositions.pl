#!/usr/bin/perl -w

#Written by Malachi Griffith
#This script parse probe files, gets probe 'gene' coordinates and prints new probe files with 'chromosome' coordinates appended
#1.) Get input files, test for neccessary columns and make note of their position, generate new probe names
#2.) For each probe file:
#    - a.) Get all probes and gene coordinate info
#    - b.) Build a list of genes found in this probe file
#    - c.) Get gene info for these genes from ALEXA
#    - d.) Go through all probes and calculate chromosome coordinates
#    - e.) Open the original probe file again
#    - f.) Write out a new probe file with chromosome coordinates appended immediately after the gene coordinates
#    - g.) Proceed to the next probe file

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

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $probe_dir = '';
my $logfile = '';
my $batch_size = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'probe_dir=s'=>\$probe_dir, 'logfile=s'=>\$logfile, 'batch_size=i'=>\$batch_size);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the full path to a directory containing input probe files to be processed using: --probe_dir", RESET;
print GREEN, "\n\t\tOutput files will be generated automatically and will have the same names with _chr_coords appended", RESET;
print GREEN, "\n\tSpecify a logfile for output using: --logfile", RESET;
print GREEN, "\n\tSpecify the number of probes to process at a time using: --batch_size (reduce if system memory is limited)", RESET;
print GREEN, "\n\nExample: getProbeChromosomePositions.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --batch_size=100000  --probe_dir=/home/user/alexa/ALEXA_version/unfiltered_probes/  --logfile=/home/user/alexa/ALEXA_version/logs/getProbeChromosomePositions_LOG.txt\n\n", RESET;

unless ($database && $server && $user && $password && $probe_dir && $logfile && $batch_size){
  print RED, "\nRequired input parameter(s)\n\n", RESET;
  exit();
}

#Open logfile for output
open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";
print LOG "\nUser Specified the following options:\ndatabase = $database\nprobe_dir = $probe_dir\nlogfile = $logfile\nbatch_size = $batch_size\n\n";

#1.) Get input files, test for neccessary columns and make note of their position, generate new probe names
my %files = &getProbeFiles('-dir'=>$probe_dir);

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

#2.) For each probe file:
foreach my $file_count (sort {$a <=> $b} keys %files){

  print BLUE, "\n\nPROCESSING FILE $file_count: $files{$file_count}{file_path}\n", RESET;
  print LOG "\n\nPROCESSING FILE $file_count: $files{$file_count}{file_path}\n";

  #Process each probe files in batches

  #a.) Get all probes and gene coordinate info
  my %probes;
  my %genes;

  #open input and output files
  my $probe_file = $files{$file_count}{file_path};
  my $new_probe_file = $files{$file_count}{new_file_path};

  open (PROBE_FILE, "$probe_file") || die "\nCould not open probe file $probe_file\n\n";
  open (OUT_FILE, ">$new_probe_file") || die "\nCould not open new probe file $new_probe_file\n\n";


  my $columns_ref = $files{$file_count}{columns};

  my $first_line = 1;
  my $line_count = 0;
  while (<PROBE_FILE>){
    chomp($_);
    my @data_line = split("\t", $_);

    #Deal with header
    if ($first_line == 1){
      $first_line = 0;

      #print the first block of columns
      for (my $i = 0; $i <= $columns_ref->{'Unit2_end'}->{column_pos}; $i++){
	print OUT_FILE "$data_line[$i]\t";
      }
      #Insert the new block
      print OUT_FILE "Unit1_start_chr\tUnit1_end_chr\tUnit2_start_chr\tUnit2_end_chr\t";

      my $total_columns = @data_line;

      #print the second block of columns
      my $continue_pos = ($columns_ref->{'Unit2_end'}->{column_pos})+1;
      for (my $j = $continue_pos; $j <= $total_columns-2; $j++){
	print OUT_FILE "$data_line[$j]\t";
      }

      #print the last column
      print OUT_FILE "$data_line[$total_columns-1]\n";

      next();
    }
    $line_count++;

    my $probe_id = $data_line[$columns_ref->{'Probe_Count'}->{column_pos}];

    #check for non-unique probe ids
    if ($probes{$probe_id}){
      print RED, "\nFound non-unique probe id!! - exiting\n\n", RESET;
      close (LOG);
      exit();
    }

    $probes{$probe_id}{gene_id} = $data_line[$columns_ref->{'Gene_ID'}->{column_pos}];
    $probes{$probe_id}{unit1_start} = $data_line[$columns_ref->{'Unit1_start'}->{column_pos}];
    $probes{$probe_id}{unit1_end} = $data_line[$columns_ref->{'Unit1_end'}->{column_pos}];
    $probes{$probe_id}{unit2_start} = $data_line[$columns_ref->{'Unit2_start'}->{column_pos}];
    $probes{$probe_id}{unit2_end} = $data_line[$columns_ref->{'Unit2_end'}->{column_pos}];
    $probes{$probe_id}{line_count} = $line_count;
    $probes{$probe_id}{data_line} = \@data_line;


    #Build a list of genes found in this probe file
    if ($data_line[$columns_ref->{'Gene_ID'}->{column_pos}] =~ /\d+/){
      $genes{$data_line[$columns_ref->{'Gene_ID'}->{column_pos}]}{tmp} = '';
    }

    if ($line_count == $batch_size){
      $line_count = 0;

      &processProbes('-probe_object'=>\%probes, '-gene_object'=>\%genes, '-columns_ref'=>$columns_ref);
      %probes = ();
      %genes = ();

      $| = 1;
      print BLUE, ".", RESET;
      $| = 0;
      next();
    }
  }
  #Process any remaining probes that did not make up a full batch
  &processProbes('-probe_object'=>\%probes, '-gene_object'=>\%genes, '-columns_ref'=>$columns_ref);

  close (PROBE_FILE);
  close (OUT_FILE);

}

#Close database connection
$alexa_dbh->disconnect();

close (LOG);

print "\n\n";

exit();


#################################################################################################################################
#getProbeFiles                                                                                                                  #
#################################################################################################################################
sub getProbeFiles{
  my %args = @_;
  my $probe_dir = $args{'-dir'};

  my %files;

  unless ($probe_dir =~ /.*\/$/){
    $probe_dir = "$probe_dir"."/";
  }

  #First make sure the specified base path exists and is a directory
  unless (-e $probe_dir && -d $probe_dir){
    print RED, "\nSpecified directory: $probe_dir does not appear valid!\n\n", RESET;
    close (LOG);
    exit();
  }

  #Get all the input files in the repeat masked result directory
  print BLUE, "\nSearching $probe_dir for files\n", RESET;
  print LOG "\nSearching $probe_dir for files\n";

  my %possible_files;
  my $possible_file_count = 0;
  opendir(DIRHANDLE, "$probe_dir") || die "\nCannot open directory: $probe_dir\n\n";
  my @test_files = readdir(DIRHANDLE);

  foreach my $test_file (@test_files){
    my $file_path = "$probe_dir"."$test_file";

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      print BLUE, "\n\t$file_path  is a directory - skipping", RESET;
      print LOG "\n\t$file_path  is a directory - skipping";
      next();
    }
    $possible_file_count++;

    $possible_files{$possible_file_count}{file_name} = $test_file;
    $possible_files{$possible_file_count}{source_dir} = $probe_dir;
    $possible_files{$possible_file_count}{file_path} = $file_path;
  }

  my $file_num = keys %possible_files;
  print BLUE, "\n\nFound $file_num probe files in the specified directory\n", RESET;
  print LOG "\n\nFound $file_num probe files in the specified directory\n";

  #Check each file for expected columns
  FILE: foreach my $file_count (sort {$a <=> $b} keys %possible_files){
      my $probe_file = $possible_files{$file_count}{file_path};

      my %columns;

      open (PROBE_FILE, "$probe_file") || die "\nCould not open probe file: $probe_file\n\n";

      #Only process the first line
      while (<PROBE_FILE>){
	my $line_record = $_;
	chomp ($line_record);

	my @line = split ("\t", $line_record);

	#Watch for the header line which is assumed to start with "Probe_Count" and contain "Sequence"
	my $column_count = 0;
	foreach my $column (@line){
	  $columns{$column}{column_pos} = $column_count;
	  $column_count++;
	}

	#Check if this probefile already appears to have chromosome coordinates
	if ($columns{'Unit1_end_chr'}){
	  print YELLOW, "\nFile: $possible_files{$file_count}{file_name} appears to already have chromosome coordinate columns - skipping\n\n", RESET;
	  print LOG "\nFile: $possible_files{$file_count}{file_name} appears to already have chromosome coordinate columns - skipping\n\n";
	  close (PROBE_FILE);
	  next FILE;
	}

	#If the file has the neccessary columns, add it to the list of files to be processed
	if ($columns{'Probe_Count'} && $columns{'Gene_ID'} && $columns{'Unit1_start'} && $columns{'Unit1_end'} && $columns{'Unit2_start'} && $columns{'Unit2_end'}){
	  $files{$file_count}{file_path} = $probe_file;
	  $files{$file_count}{file_name} = $possible_files{$file_count}{file_name};
	  $files{$file_count}{source_dir} = $possible_files{$file_count}{source_dir};
	  $files{$file_count}{columns} = \%columns;

	  #Create a new file name and path
	  my $new_file_name = "$files{$file_count}{file_name}"."_chr_coords";
	  $files{$file_count}{new_file_name} = $new_file_name;
	  $files{$file_count}{new_file_path} = "$probe_dir"."$new_file_name";

	}else{
	  print YELLOW, "\nFile: $possible_files{$file_count}{file_name} does not appear to be a valid probe file (expected column missing) - skipping\n\n", RESET;
	  print LOG "\nFile: $possible_files{$file_count}{file_name} does not appear to be a valid probe file (expected column missing) - skipping\n\n";
	}
	close (PROBE_FILE);
	next FILE;
      }
    }

  #Summarize the files found - list old and new new names
  print BLUE, "\n\nFile summary:", RESET;
  print LOG "\n\nFile summary:";
  foreach my $file_count (sort {$a <=> $b} keys %files){
    print BLUE, "\n\tFile: $file_count\n\t\tINFILE: $files{$file_count}{file_path}\n\t\tOUTFILE: $files{$file_count}{new_file_path}", RESET;
    print LOG "\n\tFile: $file_count\n\t\tINFILE: $files{$file_count}{file_path}\n\t\tOUTFILE: $files{$file_count}{new_file_path}";
  }
  print BLUE, "\n\n", RESET;
  print LOG "\n\n";

  return(%files);
}


################################################################################################################################
#processProbes                                                                                                                 #
################################################################################################################################
sub processProbes{
  my %args = @_;
  my $probes_ref = $args{'-probe_object'};
  my $genes_ref = $args{'-gene_object'};
  my $columns_ref = $args{'-columns_ref'};

  #a.) Get gene info for these genes from ALEXA
  my @gene_ids;
  my $gene_info_ref;
  my $gene_count = keys %{$genes_ref};

  if ($gene_count > 0){
    @gene_ids = keys %{$genes_ref};
    $gene_info_ref = &getGeneInfo('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no", '-silent'=>"yes");
  }

  #b.) Go through all probes and calculate chromosome coordinates
  foreach my $probe_id (keys %{$probes_ref}){

    my $gene_id = $probes_ref->{$probe_id}->{gene_id};

    #If the probe does not have a valid gene id (negative control probe for example) set all values to 'na'
    unless ($gene_id =~ /\d+/){
      $probes_ref->{$probe_id}->{unit1_start_chr} = "na";
      $probes_ref->{$probe_id}->{unit1_end_chr} = "na";
      $probes_ref->{$probe_id}->{unit2_start_chr} = "na";
      $probes_ref->{$probe_id}->{unit2_end_chr} = "na";
      next();
    }

    #get the coords for unit1 (should be defined for all probes)
    my $unit1_start = $probes_ref->{$probe_id}->{unit1_start};
    my $unit1_end = $probes_ref->{$probe_id}->{unit1_end};
    my $coords1_ref = &convertProbeCoordinates('-gene_object'=>$gene_info_ref, '-gene_id'=>$gene_id, '-start_pos'=>$unit1_start, '-end_pos'=>$unit1_end);

    $probes_ref->{$probe_id}->{unit1_start_chr} = $coords1_ref->{$gene_id}->{chr_start};
    $probes_ref->{$probe_id}->{unit1_end_chr} = $coords1_ref->{$gene_id}->{chr_end};

    #get the coords for unit2 (should be defined for junction probes)
    my $unit2_start = $probes_ref->{$probe_id}->{unit2_start};
    my $unit2_end = $probes_ref->{$probe_id}->{unit2_end};

    if ($unit2_start eq "na" || $unit2_end eq "na"){
      $probes_ref->{$probe_id}->{unit2_start_chr} = "na";
      $probes_ref->{$probe_id}->{unit2_end_chr} = "na";
    }else{
      my $coords2_ref = &convertProbeCoordinates('-gene_object'=>$gene_info_ref, '-gene_id'=>$gene_id, '-start_pos'=>$unit2_start, '-end_pos'=>$unit2_end);
      $probes_ref->{$probe_id}->{unit2_start_chr} = $coords2_ref->{$gene_id}->{chr_start};
      $probes_ref->{$probe_id}->{unit2_end_chr} = $coords2_ref->{$gene_id}->{chr_end};
    }
  }


  #c.) Now print out the result to the out file
  foreach my $probe_id (sort {$probes_ref->{$a}->{line_count} <=> $probes_ref->{$b}->{line_count}} keys %{$probes_ref}){

    my @data_line = @{$probes_ref->{$probe_id}->{data_line}};

    #print the first block of columns
    for (my $i = 0; $i <= $columns_ref->{'Unit2_end'}->{column_pos}; $i++){
      print OUT_FILE "$data_line[$i]\t";
    }
    #Insert the new block
    print OUT_FILE "$probes_ref->{$probe_id}->{unit1_start_chr}\t$probes_ref->{$probe_id}->{unit1_end_chr}\t$probes_ref->{$probe_id}->{unit2_start_chr}\t$probes_ref->{$probe_id}->{unit2_end_chr}\t";

    my $total_columns = @data_line;

    #print the second block of columns
    my $continue_pos = ($columns_ref->{'Unit2_end'}->{column_pos})+1;
    for (my $j = $continue_pos; $j <= $total_columns-2; $j++){
      print OUT_FILE "$data_line[$j]\t";
    }

    #print the last column
    print OUT_FILE "$data_line[$total_columns-1]\n";

  }

  return();
}


############################
#convertProbeCoordinates   #
############################
sub convertProbeCoordinates{
  my %args = @_;
  my $gene_object_ref = $args{'-gene_object'};
  my $gene_id = $args{'-gene_id'};
  my $start = $args{'-start_pos'};
  my $end = $args{'-end_pos'};

  my %coords;

  #Note: All gene coordinates stored in ALEXA are relative to the coding strand (i.e. start codon near beginning, end codon near the end)
  #This means that when coverting exon, or other coordinates back to the chromosome context, the original 'strand' of the gene on the chromosome
  #needs to be considered.
  my $chromosome = $gene_object_ref->{$gene_id}->{chromosome};
  my $chr_strand = $gene_object_ref->{$gene_id}->{chr_strand};
  my $chr_start = $gene_object_ref->{$gene_id}->{chr_start};
  my $chr_end = $gene_object_ref->{$gene_id}->{chr_end};
  my $gene_start = $gene_object_ref->{$gene_id}->{gene_start};
  my $gene_end = $gene_object_ref->{$gene_id}->{gene_end};

  #Make sure the supplied coordinates are actually within the specified gene
  unless ($start >= $gene_start-1 && $start <= $gene_end+1){
    print "\nStart coordinate ($start) provided to convertGeneCoordinates() does not appear valid for gene_id $gene_id\n\n";
    close (LOG);
    exit();
  }
  unless ($end >= $gene_start-1 && $end <= $gene_end+1){
    print "\nEnd coordinate ($end) provided to convertGeneCoordinates() does not appear valid for gene_id $gene_id\n\n";
    close (LOG);
    exit();
  }

  #Convert provided gene coordinates to coordinates relative to the chromosome
  if ($chr_strand == 1){
    my $query_chr_start = $chr_start + $start - 1;
    my $query_chr_end = $chr_start + $end - 1;

    $coords{$gene_id}{chr_start} = $query_chr_start;
    $coords{$gene_id}{chr_end} = $query_chr_end;
    $coords{$gene_id}{chr_name} = $chromosome;
    $coords{$gene_id}{strand} = "+";
  }elsif ($chr_strand == -1){

    my $query_chr_start = $chr_end - $end + 1;
    my $query_chr_end = $chr_end - $start + 1;

    $coords{$gene_id}{chr_start} = $query_chr_start;
    $coords{$gene_id}{chr_end} = $query_chr_end;
    $coords{$gene_id}{chr_name} = $chromosome;
    $coords{$gene_id}{strand} = "-";

  }else{
    print "\nStrand format: $chr_strand not understood by convertGeneCoordinates()!\n\n";
    close (LOG);
    exit();
  }

  return(\%coords);
}
