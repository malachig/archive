#!/usr/bin/perl -w
#Written by Malachi Griffith

#The purpose of this script is to import probe records from probe files into the database and link them to the appropriate gene records

#Note: This script can be used to import ALL probe records in these files into the database (i.e. unfiltered)
#      It will also allow views of all possible probes within the ALEXA platform, not just those that wind up on an array
#      The probes that make up each array will consist of subsets of these probe records, indicated by a seperate relational table
#      These microarray probe sets will then be associated with the experiments conducted using particular probe sets (microarray designs)
#      This script can also be used to import only filtered probes

#For purposes of clarity, I will keep the probe_IDs used in the probe files the same as the Probe.id used in the database
#Probe IDs will only start back from '1' if they involve a different species (or possibly different EnsEMBL build).
#In any case, such probes would occur in a seperate instance of ALEXA (one for each species at the least)

#Probes should come from several seperate files, exonJunction probes, exonBoundary probes, exon probes and intron probes
#Negative control probes are not associated with any gene (by definition) and therefore will not be associated with gene records
#They will be imported for completeness though

use strict;
use Data::Dumper;
use Getopt::Long;
use DBI;
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

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $probe_dir = '';
my $block_size = '';
my $populate_database = '';
my $logfile = '';

GetOptions ('database=s'=>\$database, 'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'probe_dir=s'=>\$probe_dir, 'block_size=i'=>\$block_size, 'populate_database=s'=>\$populate_database,
	    'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nThis script parses probe files in a specfified directory and imports probe records to the specified database", RESET;
print GREEN, "\nNOTE: Import all probes for all probe types", RESET;
print GREEN, "\nUsage:", RESET;
print GREEN, "\n\tSpecify the target server using: --server", RESET;
print GREEN, "\n\tSpecify the target database using: --database", RESET;
print GREEN, "\n\tSpecify the user and password using: --user and --password", RESET;
print GREEN, "\n\t\tMake sure you chose the correct database and server!!", RESET;
print GREEN, "\n\tSpecify the directory containing probe records using: --probe_dir", RESET;
print GREEN, "\n\tSpecify the block size for import statements using: --block_size", RESET;
print GREEN, "\n\tAfter testing this script use:  --populate_database=yes to insert records to database", RESET;
print GREEN, "\n\nExample: importProbeRecords.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --probe_dir=/home/user/alexa/ALEXA_version/unfiltered_probes/  --block_size=1000  --populate_database=no  --logfile=/home/user/alexa/ALEXA_version/logs/database_population/importProbeRecords_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $probe_dir && $block_size && $populate_database && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#Print out the parameters supplied by the user to the logfile for future reference
print LOG "\nUser Specified the following options:\ndatabase = $database\nprobe_dir = $probe_dir\nblock_size = $block_size\npopulate_database = $populate_database\nlogfile = $logfile\n\n";

#Establish connection with the Alternative Splicing Expression database using details provided by user
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

#1.) Get probe files to be processed - order them by probe IDs contained within them

my @required_columns = qw(Probe_Count ProbeSet_ID Gene_ID Sequence Probe_length Probe_Tm Probe_Type Exon1_IDs Unit1_start Unit1_end Exon2_IDs Unit2_start Unit2_end Unit1_start_chr Unit1_end_chr Unit2_start_chr Unit2_end_chr masked_bases SimFold_score PairFold_score MdustBases_cutoff_11);
my @optional_columns = qw(mRNA_largestNon-TargetHitLength mRNA_largestTargetHitLength EST_largestNon-TargetHitLength EST_largestTargetHitLength enst_largestNon-TargetHitLength enst_largestTargetHitLength probe_largestNon-TargetHitLength probe_largestTargetHitLength genomic_largestNon-TargetHitLength withinGene_largestNon-TargetHitLength withinGene_largestTargetHitLength Exons_Skipped);

my %files = %{&getProbeFiles('-probe_dir'=>$probe_dir)};

my $total_probe_insert_count = 0;
my $total_geneprobe_insert_count = 0;

foreach my $file_count (sort {$files{$a}->{min_probe_id} <=> $files{$b}->{min_probe_id}} keys %files){

  print BLUE, "\nProcessing file: $files{$file_count}{file_path}\n\n", RESET;
  print LOG "\nProcessing file: $files{$file_count}{file_path}\n\n";

  #Foreach probe file do the following
  #2.) Import data from the user specified probe file
  #3.) Everytime a block of probes are acquired from the probe file, insert them into the database 
  #    (saves memory -rather than loading all probes first)

  my %columns = %{$files{$file_count}{columns}};
  &importProbeData('-input_file'=>$files{$file_count}{file_path}, '-columns'=>\%columns, '-block_size'=>$block_size);

}

unless ($populate_database eq "yes"){
  print YELLOW, "\nRecords not inserted to database!  Once testing is complete use: --populate_database=yes\n\n", RESET;
}

print BLUE, "\n\nInserted a total of $total_probe_insert_count Probe records to the 'Probe' table", RESET;
print BLUE, "\nInserted a total of $total_geneprobe_insert_count GeneProbe records to the 'GeneProbe' table\n\n", RESET;
print LOG "\n\nInserted a total of $total_probe_insert_count Probe records to the 'Probe' table";
print LOG "\nInserted a total of $total_geneprobe_insert_count GeneProbe records to the 'GeneProbe' table\n\n";

#Close database connection
$alexa_dbh->disconnect();

close (LOG);

exit();


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
    $alexa_dbh->disconnect(); close (LOG);
    exit();
  }

  #Get files from this directory
  print BLUE, "\nSearching $probe_dir for probe files", RESET;
  print LOG "\nSearching $probe_dir for probe files";

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

  #Check each file for expected columns, get probe id ranges to allow ordering of files
  my %files;
 FILE: foreach my $file_count (sort {$a <=> $b} keys %possible_files){
    my $probe_file = $possible_files{$file_count}{file_path};

    print BLUE, "\n\nExamining probe file: $probe_file\n", RESET;
    print LOG "\n\nExamining probe file: $probe_file\n";

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

	  #Summarize the optional columns that are missing from this file
	  print BLUE, "\nTesting for presence of optional columns - check if they seem okay:", RESET;
	  print LOG "\nTesting for presence of optional columns - check if they seem okay:";
	  foreach my $opt_col (@optional_columns){
	    unless ($columns{$opt_col}){
	      print YELLOW, "\n\t$opt_col column missing from this probe file", RESET;
	      print LOG "\n\t$opt_col column missing from this probe file";
	    }
	  }
	  print BLUE, "\n\n", RESET;
	  print LOG "\n\n";

	  $files{$file_count}{file_path} = $probe_file;
	  $files{$file_count}{file_name} = $possible_files{$file_count}{file_name};
	  $files{$file_count}{source_dir} = $possible_files{$file_count}{source_dir};
	  $files{$file_count}{columns} = \%columns;

	}else{
	  print RED, "\nFile: $possible_files{$file_count}{file_name} does not appear to be a complete probe file.", RESET;
	  print RED, "\n\tColumns missing: @missing_columns", RESET;
	  print RED, "\n\tComplete appropriate tests before attempting to import probe records\n\n", RESET;
	  $alexa_dbh->disconnect(); close (LOG);
	  exit();
	}
	next();
      }
      #Process data lines
      if ($first_data == 1){

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
  print LOG "\n\nFile summary:";
  foreach my $file_count (sort {$files{$a}->{min_probe_id} <=> $files{$b}->{min_probe_id}} keys %files){
    print BLUE, "\n\tFile: $file_count\n\t\tINFILE: $files{$file_count}{file_path} (probe_ids: $files{$file_count}{min_probe_id} - $files{$file_count}{max_probe_id})", RESET;
    print LOG "\n\tFile: $file_count\n\t\tINFILE: $files{$file_count}{file_path} (probe_ids: $files{$file_count}{min_probe_id} - $files{$file_count}{max_probe_id})";
  }
  print BLUE, "\n\n", RESET;
  print LOG "\n\n";

  return(\%files);
}


###########################################################################################################
#Import Probe Data from an input file to the database
###########################################################################################################

sub importProbeData{
  my %args = @_;
  my $probe_file = $args{'-input_file'};
  my %columns = %{$args{'-columns'}};
  my $block_size = $args{'-block_size'};

  #########################################################################################################
  #1.) Load all probes info from an input probe file and build a probe object as a hash                   #
  #########################################################################################################

  my $blocks_imported = 0;
  my $probe_count = 0;
  my $total_probe_count = 0;
  my %probes;

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
      $alexa_dbh->disconnect(); close (LOG);
      exit();
    }
    $probe_count++;
    $total_probe_count++;

    #Conduct check of data types to catch missing values or incorrect formats!

    #ProbeSet_ID
    if ($columns{ProbeSet_ID}){
      my $probeset_id = $probe_line[$columns{ProbeSet_ID}{column_pos}];
      if ($probeset_id =~ /\d+/){
	$probes{$probe_id}{probeset_id} = $probeset_id;
      }else{
	print RED, "\nProbeset ID value missing - or incorrect format\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }

    #Gene_ID
    $probes{$probe_id}{gene_id} = $probe_line[$columns{Gene_ID}{column_pos}];

    #Sequence
    $probes{$probe_id}{sequence} = $probe_line[$columns{Sequence}{column_pos}];

    #Probe_length
    if ($columns{Probe_length}){
      my $probe_length = $probe_line[$columns{Probe_length}{column_pos}];
      if ($probe_length =~ /\d+/){
	$probes{$probe_id}{probe_length} = $probe_length;
      }else{
	print RED, "\nProbe length value missing - or incorrect format: $probe_length\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }

    #Probe_Tm
    if ($columns{Probe_Tm}){
      my $probe_tm = $probe_line[$columns{Probe_Tm}{column_pos}];
      if ($probe_tm =~ /\d+\.\d+/){
	$probes{$probe_id}{probe_tm} = $probe_tm;
      }else{
	print RED, "\nTm value missing - or incorrect format: $probe_tm\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }

    #Probe_Type
    $probes{$probe_id}{probe_type} = $probe_line[$columns{Probe_Type}{column_pos}];

    #Exon1_IDs
    if ($columns{Exon1_IDs}){
      my $exon1_ids = $probe_line[$columns{Exon1_IDs}{column_pos}];
      if ($exon1_ids =~ /\d+/ || $exon1_ids eq "na"){
	$probes{$probe_id}{exon1_ids} = $exon1_ids;
      }else{
	print RED, "\nExon1 IDs missing - or incorrect format: $exon1_ids\tProbe = $probe_id", RESET;
	exit();
      }
    }

    #Unit1_start
    if ($columns{Unit1_start}){
      my $unit1_start = $probe_line[$columns{Unit1_start}{column_pos}];
      if ($unit1_start =~ /^\d+$/ || $unit1_start eq "na"){
	$probes{$probe_id}{unit1_start} = $unit1_start;
      }else{
	print RED, "\nUnit1 Start missing - or incorrect format: $unit1_start\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }

    #Unit1_end,
    if ($columns{Unit1_end}){
      my $unit1_end = $probe_line[$columns{Unit1_end}{column_pos}];
      if ($unit1_end =~ /^\d+$/ || $unit1_end eq "na"){
	$probes{$probe_id}{unit1_end} = $unit1_end;
      }else{
	print RED, "\nUnit1 End missing - or incorrect format: $unit1_end\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }

    #Exon2_IDs
    if ($columns{Exon2_IDs}){
      my $exon2_ids = $probe_line[$columns{Exon2_IDs}{column_pos}];
      if ($exon2_ids =~ /\d+/ || $exon2_ids eq "na"){
	$probes{$probe_id}{exon2_ids} = $exon2_ids;
      }else{
	print RED, "\nExon2 IDs missing - or incorrect format: $exon2_ids\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }

    #Unit2_start
    if ($columns{Unit2_start}){
      my $unit2_start = $probe_line[$columns{Unit2_start}{column_pos}];
      if ($unit2_start =~ /^\d+$/ || $unit2_start eq "na"){
	$probes{$probe_id}{unit2_start} = $unit2_start;
      }else{
	print RED, "\nUnit2 Start missing - or incorrect format: $unit2_start\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }

    #Unit2_end
    if ($columns{Unit2_end}){
      my $unit2_end = $probe_line[$columns{Unit2_end}{column_pos}];
      if ($unit2_end =~ /^\d+$/ || $unit2_end eq "na"){
	$probes{$probe_id}{unit2_end} = $unit2_end;
      }else{
	print RED, "\nUnit2 End missing - or incorrect format: $unit2_end\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }

    #Unit1_start_chr
    if ($columns{Unit1_start_chr}){
      my $unit1_start_chr = $probe_line[$columns{Unit1_start_chr}{column_pos}];
      if ($unit1_start_chr =~ /^\d+$/ || $unit1_start_chr eq "na"){
	$probes{$probe_id}{unit1_start_chr} = $unit1_start_chr;
      }else{
	print RED, "\nUnit1_Start_Chr missing - or incorrect format: $unit1_start_chr\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }

    #Unit1_end_chr
    if ($columns{Unit1_end_chr}){
      my $unit1_end_chr = $probe_line[$columns{Unit1_end_chr}{column_pos}];
      if ($unit1_end_chr =~ /^\d+$/ || $unit1_end_chr eq "na"){
	$probes{$probe_id}{unit1_end_chr} = $unit1_end_chr;
      }else{
	print RED, "\nUnit1_End_Chr missing - or incorrect format: $unit1_end_chr\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }

    #Unit2_start_chr
    if ($columns{Unit2_start_chr}){
      my $unit2_start_chr = $probe_line[$columns{Unit2_start_chr}{column_pos}];
      if ($unit2_start_chr =~ /^\d+$/ || $unit2_start_chr eq "na"){
	$probes{$probe_id}{unit2_start_chr} = $unit2_start_chr;
      }else{
	print RED, "\nUnit2_Start_Chr missing - or incorrect format: $unit2_start_chr\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }

    #Unit2_end_chr
    if ($columns{Unit2_end_chr}){
      my $unit2_end_chr = $probe_line[$columns{Unit2_end_chr}{column_pos}];
      if ($unit2_end_chr =~ /^\d+$/ || $unit2_end_chr eq "na"){
	$probes{$probe_id}{unit2_end_chr} = $unit2_end_chr;
      }else{
	print RED, "\nUnit2_End_Chr missing - or incorrect format: $unit2_end_chr\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }

    #masked_bases
    if ($columns{masked_bases}){
      my $masked_bases = $probe_line[$columns{masked_bases}{column_pos}];
      if ($masked_bases =~ /^\d+$/){
	$probes{$probe_id}{masked_bases} = $masked_bases;
      }elsif($masked_bases eq "na"){
	$probes{$probe_id}{masked_bases} = 0;
      }else{
	print RED, "\nmasked_bases missing - or incorrect format: $masked_bases\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }

    #SimFold_Score
    if ($columns{SimFold_score}){
      my $simfold_score = $probe_line[$columns{SimFold_score}{column_pos}];
      if ($simfold_score =~ /\d+\.\d+/ || $simfold_score =~ /^0$/){
	$probes{$probe_id}{simfold_score} = $simfold_score;
      }else{
	print RED, "\nSimfold Score missing - or incorrect format\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }

    #PairFold_score
    if ($columns{PairFold_score}){
      my $pairfold_score = $probe_line[$columns{PairFold_score}{column_pos}];
      if ($pairfold_score =~ /\d+\.\d+/ || $pairfold_score =~ /^0$/){
	$probes{$probe_id}{pairfold_score} = $pairfold_score;
      }else{
	print RED, "\nPairfold Score missing - or incorrect format\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }

    #MdustBases_cutoff_11
    if ($columns{MdustBases_cutoff_11}){
      my $mdust_bases =  $probe_line[$columns{MdustBases_cutoff_11}{column_pos}];
      if ($mdust_bases =~ /^\d+/){
	$probes{$probe_id}{mdust_bases} = $mdust_bases;
      }else{
	print RED, "\nMdust Score missing - or incorrect format\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }

    #Exons_Skipped
    if ($columns{Exons_Skipped}){
      my $exons_skipped = $probe_line[$columns{Exons_Skipped}{column_pos}];
      if ($exons_skipped =~ /\d+/ || $exons_skipped eq "na"){
	$probes{$probe_id}{exons_skipped} = $exons_skipped;
      }else{
	print RED, "\nExons skipped Score missing - or incorrect format: $exons_skipped\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }else{
      $probes{$probe_id}{exons_skipped} = 'na';
    }

    #mRNA_largestTargetHitLength
    if ($columns{'mRNA_largestTargetHitLength'}){
      my $mrna_target_hit_length = $probe_line[$columns{'mRNA_largestTargetHitLength'}{column_pos}];
      if ($mrna_target_hit_length =~ /\d+/){
	$probes{$probe_id}{mrna_t_hl} = $mrna_target_hit_length;
      }else{
	print RED, "\nmRNA Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }else{
      $probes{$probe_id}{mrna_t_hl} = 'NULL';
    }

    #mRNA_largestNon-TargetHitLength
    if ($columns{'mRNA_largestNon-TargetHitLength'}){
      my $mrna_non_target_hit_length = $probe_line[$columns{'mRNA_largestNon-TargetHitLength'}{column_pos}];
      if ($mrna_non_target_hit_length =~ /\d+/){
	$probes{$probe_id}{mrna_nt_hl} = $mrna_non_target_hit_length;
      }else{
	print RED, "\nmRNA Non-Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }else{
      $probes{$probe_id}{mrna_nt_hl} = 'NULL';
    }

    #EST_largestTargetHitLength
    if ($columns{'EST_largestTargetHitLength'}){
      my $est_target_hit_length = $probe_line[$columns{'EST_largestTargetHitLength'}{column_pos}];
      if ($est_target_hit_length =~ /\d+/){
	$probes{$probe_id}{est_t_hl} = $est_target_hit_length;
      }else{
	print RED, "\nEST Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }else{
      $probes{$probe_id}{est_t_hl} = 'NULL';
    }

    #EST_largestNon-TargetHitLength
    if ($columns{'EST_largestNon-TargetHitLength'}){
      my $est_non_target_hit_length = $probe_line[$columns{'EST_largestNon-TargetHitLength'}{column_pos}];
      if ($est_non_target_hit_length =~ /\d+/){
	$probes{$probe_id}{est_nt_hl} = $est_non_target_hit_length;
      }else{
	print RED, "\nEST Non-Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }else{
      $probes{$probe_id}{est_nt_hl} = 'NULL';
    }

    #enst_largestTargetHitLength
    if ($columns{'enst_largestTargetHitLength'}){
      my $enst_target_hit_length = $probe_line[$columns{'enst_largestTargetHitLength'}{column_pos}];
      if ($enst_target_hit_length =~ /\d+/){
	$probes{$probe_id}{enst_t_hl} = $enst_target_hit_length;
      }else{
	print RED, "\nEnst Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }else{
      $probes{$probe_id}{enst_t_hl} = 'NULL';
    }

    #enst_largestNon-TargetHitLength
    if ($columns{'enst_largestNon-TargetHitLength'}){
      my $enst_non_target_hit_length = $probe_line[$columns{'enst_largestNon-TargetHitLength'}{column_pos}];
      if ($enst_non_target_hit_length =~ /\d+/){
	$probes{$probe_id}{enst_nt_hl} = $enst_non_target_hit_length;
      }else{
	print RED, "\nEnst Non-Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }else{
      $probes{$probe_id}{enst_nt_hl} = 'NULL';
    }

    #probe_largestTargetHitLength
    if ($columns{'probe_largestTargetHitLength'}){
      my $probe_target_hit_length = $probe_line[$columns{'probe_largestTargetHitLength'}{column_pos}];
      if ($probe_target_hit_length =~ /\d+/){
	$probes{$probe_id}{probe_t_hl} = $probe_target_hit_length;
      }else{
	print RED, "\nProbe Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }else{
      $probes{$probe_id}{probe_t_hl} = 'NULL';
    }

    #probe_largestNon-TargetHitLength
    if ($columns{'probe_largestNon-TargetHitLength'}){
      my $probe_non_target_hit_length = $probe_line[$columns{'probe_largestNon-TargetHitLength'}{column_pos}];
      if ($probe_non_target_hit_length =~ /\d+/){
	$probes{$probe_id}{probe_nt_hl} = $probe_non_target_hit_length;
      }else{
	print RED, "\nProbe Non-Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }else{
      $probes{$probe_id}{probe_nt_hl} = 'NULL';
    }

    #genomic_largestNon-TargetHitLength
    if ($columns{'genomic_largestNon-TargetHitLength'}){
      my $genomic_non_target_hit_length = $probe_line[$columns{'genomic_largestNon-TargetHitLength'}{column_pos}];
      if ($genomic_non_target_hit_length =~ /\d+/){
	$probes{$probe_id}{genomic_nt_hl} = $genomic_non_target_hit_length;
      }else{
	print RED, "\nGenomic Non-Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }else{
      $probes{$probe_id}{genomic_nt_hl} = 'NULL';
    }

    #withinGene_largestTargetHitLength
    if ($columns{'withinGene_largestTargetHitLength'}){
      my $within_gene_target_hit_length = $probe_line[$columns{'withinGene_largestTargetHitLength'}{column_pos}];
      if ($within_gene_target_hit_length =~ /\d+/){
	$probes{$probe_id}{within_gene_t_hl} = $within_gene_target_hit_length;
      }else{
	print RED, "\nwithinGene Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }else{
      $probes{$probe_id}{within_gene_t_hl} = 'NULL';
    }

    #withinGene_largestNon-TargetHitLength
    if ($columns{'withinGene_largestNon-TargetHitLength'}){
      my $within_gene_non_target_hit_length = $probe_line[$columns{'withinGene_largestNon-TargetHitLength'}{column_pos}];
      if ($within_gene_non_target_hit_length =~ /\d+/){
	$probes{$probe_id}{within_gene_nt_hl} = $within_gene_non_target_hit_length;
      }else{
	print RED, "\nwithinGene Non-Target Hit Length value missing - or incorrect format\tProbe = $probe_id", RESET;
	$alexa_dbh->disconnect(); close (LOG);
	exit();
      }
    }else{
      $probes{$probe_id}{within_gene_nt_hl} = 'NULL';
    }

    #Importing probe data to the database and reseting variables
    if ($probe_count == $block_size){

      #Once a set of probe data is acquired, dump it to the database and start over again
      if ($populate_database eq "yes"){
	&insertProbeRecords('-dbh'=>$alexa_dbh, '-probe_object'=>\%probes);
      }
      $blocks_imported++;
      $probe_count = 0;
      %probes = ();

      print BLUE, "\n\t$blocks_imported: blocks of $block_size processed", RESET;
      print LOG "\n\t$blocks_imported: blocks of $block_size processed";
    }
  }
  close (PROBES);

  #Import the remaining probes that did not fill a complete block
  my $final_block_size = keys %probes;
  if ($populate_database eq "yes" && $final_block_size > 0){
    &insertProbeRecords('-dbh'=>$alexa_dbh, '-probe_object'=>\%probes);
    print BLUE, "\n\t$blocks_imported: block of $final_block_size processed", RESET;
    print LOG "\n\t$blocks_imported: block of $final_block_size processed";
  }

  #Clean-up and summary
  $blocks_imported++;
  $probe_count = 0;
  %probes = ();

  print BLUE, "\n\nFound a total of $total_probe_count probes in the input file\n\n", RESET;
  print LOG "\n\nFound a total of $total_probe_count probes in the input file\n\n";

  return ();
}

########################################################################################################################
#2.) Import probe records to the 'Probe' table                                                                         #
########################################################################################################################
sub insertProbeRecords{
  my %args = @_;
  my $dbh = $args{'-dbh'};
  my %p = %{$args{'-probe_object'}};

  #1.) Create probe records

  #Database Fields
  #id, probe_count, probe_set_id, sequence, probe_length, tm, type, masked_bases, sim_fold_score, pair_fold_score
  #mdust_bases, mrna_target_hit_length, mrna_non_target_hit_length, est_target_hit_length, est_non_target_hit_length,
  #enst_target_hit_length, enst_non_target_hit_length, probe_target_hit_length, probe_non_target_hit_length,
  #genomic_non_target_hit_length, within_gene_target_hit_length, within_gene_non_target_hit_length, exons_skipped

  my $probe_count = 0;
  my $gene_probe_count = 0;
  my $total_probes = keys %p;
  my $nc_probes_found = 0;

  my $sql = "INSERT INTO Probe (probe_count, probe_set_id, sequence, probe_length, tm, type, masked_bases, sim_fold_score, pair_fold_score, mdust_bases, mrna_target_hit_length, mrna_non_target_hit_length, est_target_hit_length, est_non_target_hit_length, enst_target_hit_length, enst_non_target_hit_length, probe_target_hit_length, probe_non_target_hit_length, genomic_non_target_hit_length, within_gene_target_hit_length, within_gene_non_target_hit_length, exons_skipped) VALUES ";

  foreach my $probe_id (sort {$a <=> $b} keys %p){
    $probe_count++;

    #Watch for a block of NC probes - Since probe files are processed one at a time this will be all or nothing (all NC or not)
    if ($p{$probe_id}{probe_type} eq "Control-Negative"){
      $nc_probes_found = 1;
    }

    #If a block is complete, or the last probe record has been reached, perform the insert
    if ($probe_count == $total_probes){

      #Build the final part of the insert statement, remember that non-numeric values should be quoted
      $sql = "$sql"."($probe_id,$p{$probe_id}{probeset_id},\'$p{$probe_id}{sequence}\',$p{$probe_id}{probe_length},$p{$probe_id}{probe_tm},\'$p{$probe_id}{probe_type}\',$p{$probe_id}{masked_bases},$p{$probe_id}{simfold_score},$p{$probe_id}{pairfold_score},$p{$probe_id}{mdust_bases},$p{$probe_id}{mrna_t_hl},$p{$probe_id}{mrna_nt_hl},$p{$probe_id}{est_t_hl},$p{$probe_id}{est_nt_hl},$p{$probe_id}{enst_t_hl},$p{$probe_id}{enst_nt_hl},$p{$probe_id}{probe_t_hl},$p{$probe_id}{probe_nt_hl},$p{$probe_id}{genomic_nt_hl},$p{$probe_id}{within_gene_t_hl},$p{$probe_id}{within_gene_nt_hl},\'$p{$probe_id}{exons_skipped}\');";

      #Actually insert the multi-insert statement
      my $sth = $dbh->prepare("$sql");
      $sth->execute();
      $sth->finish();
      $total_probe_insert_count++;

    }else{
      #otherwise, continue building the multi-insert statement
      $sql = "$sql"."($probe_id,$p{$probe_id}{probeset_id},\'$p{$probe_id}{sequence}\',$p{$probe_id}{probe_length},$p{$probe_id}{probe_tm},\'$p{$probe_id}{probe_type}\',$p{$probe_id}{masked_bases},$p{$probe_id}{simfold_score},$p{$probe_id}{pairfold_score},$p{$probe_id}{mdust_bases},$p{$probe_id}{mrna_t_hl},$p{$probe_id}{mrna_nt_hl},$p{$probe_id}{est_t_hl},$p{$probe_id}{est_nt_hl},$p{$probe_id}{enst_t_hl},$p{$probe_id}{enst_nt_hl},$p{$probe_id}{probe_t_hl},$p{$probe_id}{probe_nt_hl},$p{$probe_id}{genomic_nt_hl},$p{$probe_id}{within_gene_t_hl},$p{$probe_id}{within_gene_nt_hl},\'$p{$probe_id}{exons_skipped}\'),";
      $total_probe_insert_count++;
    }
  }

  #Unless we are dealing with negative control probes, associate each probe record with a Gene
  unless($nc_probes_found == 1){
    #2.) Create GeneProbe records
    #Get the last insert ID from the insert statement executed above.  Use this to make sure the probe IDs for gene/probe relationships
    #correspond even if they have become out of sync with the probe count values 

    #Database Fields
    #id, fk_Gene__id, fk_Probe__id, unit1_start, unit1_end, unit1_start_chr, unit1_end_chr, exon1_ids,
    #unit2_start, unit2_end, unit2_start_chr, unit2_end_chr, exon2_ids

    #Get last insert ID
    my $sql_insert_id = "SELECT LAST_INSERT_ID();";
    my $sth_insert_id = $dbh->prepare("$sql_insert_id");
    $sth_insert_id->execute();
    my $last_insert_id = $sth_insert_id->fetchrow_array(); #Note this actually gives the FIRST id of the inserted block
    $sth_insert_id->finish();

    my @ids = ($last_insert_id .. ($last_insert_id + ($total_probes-1)));

    #print "\nLast Insert ID = $last_insert_id\tTotal Probes = $total_probes\t IDS: @ids\n\n";

    #Insert the GeneProbe Records
    my $sql2 = "INSERT INTO GeneProbe (fk_Gene__id, fk_Probe__id, unit1_start, unit1_end, unit1_start_chr, unit1_end_chr, exon1_ids, unit2_start, unit2_end, unit2_start_chr, unit2_end_chr, exon2_ids) VALUES ";

    foreach my $probe_id (sort {$a <=> $b} keys %p){


      $gene_probe_count++;

      my $fk_probe_id = shift @ids;
      my $fk_gene_id = $p{$probe_id}{gene_id};

      #If a block is complete, or the last probe record has been reached, perform the insert
      if ($gene_probe_count == $total_probes){

	#Build the final part of the insert statement, remember that non-numeric values should be quoted
	$sql2 = "$sql2"."($fk_gene_id,$fk_probe_id,\'$p{$probe_id}{unit1_start}\',\'$p{$probe_id}{unit1_end}\',\'$p{$probe_id}{unit1_start_chr}\',\'$p{$probe_id}{unit1_end_chr}\',\'$p{$probe_id}{exon1_ids}\',\'$p{$probe_id}{unit2_start}\',\'$p{$probe_id}{unit2_end}\',\'$p{$probe_id}{unit2_start_chr}\',\'$p{$probe_id}{unit2_end_chr}\',\'$p{$probe_id}{exon2_ids}\');";

	#Actually insert the multi-insert statement
	my $sth2 = $dbh->prepare("$sql2");
	$sth2->execute();
	$sth2->finish();
	$total_geneprobe_insert_count++;

      }else{
	#otherwise, continue building the multi-insert statement
	$sql2 = "$sql2"."($fk_gene_id,$fk_probe_id,\'$p{$probe_id}{unit1_start}\',\'$p{$probe_id}{unit1_end}\',\'$p{$probe_id}{unit1_start_chr}\',\'$p{$probe_id}{unit1_end_chr}\',\'$p{$probe_id}{exon1_ids}\',\'$p{$probe_id}{unit2_start}\',\'$p{$probe_id}{unit2_end}\',\'$p{$probe_id}{unit2_start_chr}\',\'$p{$probe_id}{unit2_end_chr}\',\'$p{$probe_id}{exon2_ids}\'),";

	$total_geneprobe_insert_count++;
      }
    }
  }
  return();
}

