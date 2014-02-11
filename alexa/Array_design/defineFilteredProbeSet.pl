#!/usr/bin/perl -w
#Written by Malachi Griffith

#The purpose of this script is to define a set of filtered probes within the ALEXA database.
#All experimental probes are imported to the database.  This simply identifies which successfully pass the filtering steps.

#Probes should come from several seperate files, exonJunction probes, exonBoundary probes, exon probes and intron probes
#Negative control probes are not associated with any gene (by definition) and therefore will not be defined in the database

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
my $probeset_description = '';
my $logfile = '';

GetOptions ('database=s'=>\$database, 'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'probe_dir=s'=>\$probe_dir, 'block_size=i'=>\$block_size, 'populate_database=s'=>\$populate_database, 
	    'probeset_description=s'=>\$probeset_description, 'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nThis script parses filtered probe files in a specfified directory and defines them as a set in the specified database", RESET;
print GREEN, "\nNOTE: Import all probes for all probe types all quality tests are complete", RESET;
print GREEN, "\nUsage:", RESET;
print GREEN, "\n\tSpecify the target server using: --server", RESET;
print GREEN, "\n\tSpecify the target database using: --database", RESET;
print GREEN, "\n\tSpecify the user and password using: --user and --password", RESET;
print GREEN, "\n\t\tMake sure you chose the correct database and server!!", RESET;
print GREEN, "\n\tSpecify the directory containing probe records using: --probe_dir", RESET;
print GREEN, "\n\tSpecify the block size for import statements using: --block_size", RESET;
print GREEN, "\n\tAfter testing this script use:  --populate_database=yes to insert records to database", RESET;
print GREEN, "\n\tProvide a description for this set of filtered probes using: --probeset_description", RESET;
print GREEN, "\n\nExample: defineFilteredProbeSet.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --probe_dir=/home/user/alexa/ALEXA_version/filtered_probes/  --block_size=1000  --populate_database=no  --probeset_description='Standard filtering options'  --logfile=/home/user/alexa/ALEXA_version/logs/database_population/defineFilteredProbeSet_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $probe_dir && $block_size && $populate_database && $probeset_description && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#Print out the parameters supplied by the user to the logfile for future reference
print LOG "\nUser Specified the following options:\ndatabase = $database\nprobe_dir = $probe_dir\nblock_size = $block_size\npopulate_database = $populate_database\nprobeset_description = $probeset_description\nlogfile = $logfile\n\n";

#Establish connection with the Alternative Splicing Expression database using details provided by user
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

#1.) Get probe files to be processed - order them by probe IDs contained within them
my @required_columns = qw(Probe_Count);
my %files = %{&getProbeFiles('-probe_dir'=>$probe_dir)};

#2.) Define the probe set entry which all filtered probes will be added to
my $fk_Probe_set__id;
my $type = "Filtered_Probes";

if ($populate_database eq "yes"){

  my $sql = "INSERT INTO Probe_set (type, description) VALUES (\'$type\', \'$probeset_description\');";
  my $sth = $alexa_dbh->prepare("$sql");
  $sth->execute();
  $sth->finish();

  my $sql_insert_id = "SELECT LAST_INSERT_ID();";
  my $sth_insert_id = $alexa_dbh->prepare("$sql_insert_id");
  $sth_insert_id->execute();
  my $last_insert_id = $sth_insert_id->fetchrow_array();
  $sth_insert_id->finish();

  $fk_Probe_set__id = $last_insert_id;
  print BLUE "\nAdded Probe_set record and retrieved id: $fk_Probe_set__id\n\n", RESET;
  print LOG "\nAdded Probe_set record and retrieved id: $fk_Probe_set__id\n\n";
}


#3.) Define relationships between all filtered probes found and the new probe set entry

my $total_probe_insert_count = 0;

foreach my $file_count (sort {$files{$a}->{min_probe_id} <=> $files{$b}->{min_probe_id}} keys %files){

  print BLUE, "\nProcessing file: $files{$file_count}{file_path}\n\n", RESET;
  print LOG "\nProcessing file: $files{$file_count}{file_path}\n\n";

  #Foreach probe file do the following
  #  - Import probe to probe-set relationships
  #  - Everytime a block of probes are acquired from the probe file, insert them into the database 
  #  - (saves memory -rather than loading all probes first)

  my %columns = %{$files{$file_count}{columns}};
  &defineProbeSetRelationships('-input_file'=>$files{$file_count}{file_path}, '-columns'=>\%columns, '-block_size'=>$block_size);

}

unless ($populate_database eq "yes"){
  print YELLOW, "\nRecords not inserted to database!  Once testing is complete use: --populate_database=yes\n\n", RESET;
}

print BLUE, "\n\nInserted a total of $total_probe_insert_count Probe records to the 'ProbeProbe_set' table\n\n", RESET;
print LOG "\n\nInserted a total of $total_probe_insert_count Probe records to the 'ProbeProbe_set' table\n\n";

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

sub defineProbeSetRelationships{
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

  #Possible file column headings:
  #Probe_Count

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
    $probes{$probe_id}{tmp} = '';

    $probe_count++;
    $total_probe_count++;

    #Define probe-probe_set relationship in the database and reset variables
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
  #id, fk_Probe_set__id, fk_Probe__id

  my $probe_count = 0;
  my $gene_probe_count = 0;
  my $total_probes = keys %p;

  my $sql = "INSERT INTO ProbeProbe_set (fk_Probe_set__id, fk_Probe__id) VALUES ";

  foreach my $probe_id (sort {$a <=> $b} keys %p){
    $probe_count++;

    #If a block is complete, or the last probe record has been reached, perform the insert
    if ($probe_count == $total_probes){

      #Build the final part of the insert statement, remember that non-numeric values should be quoted
      $sql = "$sql"."($fk_Probe_set__id, $probe_id);";

      #Actually insert the multi-insert statement
      my $sth = $dbh->prepare("$sql");
      $sth->execute();
      $sth->finish();
      $total_probe_insert_count++;

    }else{
      #otherwise, continue building the multi-insert statement
      $sql = "$sql"."($fk_Probe_set__id, $probe_id),";
      $total_probe_insert_count++;
    }
  }

  return();
}

