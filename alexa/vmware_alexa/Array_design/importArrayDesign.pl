#!/usr/bin/perl -w
#Written by Malachi Griffith

#The purpose of this script is to import an array design using a NimbleGen array design file
#These files are tab-delimited and have the format:
#   PROBE_DESIGN_ID\tCONTAINER\tDESIGN_NOTE\tSELECTION_CRITERIA\tSEQ_ID\tPROBE_SEQUENCE\tMISMATCH\tMATCH_INDEX\tFEATURE_ID\tROW_NUM\tCOL_NUM\tPROBE_CLASS PROBE_ID
#   \tPOSITION\tDESIGN_ID\tX\tY
#
#Example file: /home/malachig/AlternativeSplicing/Array_analysis/NimbleGen_Data/Nimblegen_pilot_010606/DesignFiles/2005-10-28_Marra_Human_Splice.ndf
#Only the x and y coordinates, and row and col numbers will be imported
#NimbleGen's control probes will not be imported
#The probe_id from each line of this file will be used to associate each record to a probe record

#Tables to be populated are: 'Microarray','ArraySpot', and 'ArraySpot_negativeControl'

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use DBI;
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
my $design_file = '';
my $array_name = '';
my $array_type = '';
my $manufacturer = '';
my $probe_length = '';
my $order_date = '';
my $block_size = '';
my $update_database = '';

GetOptions ('database=s'=>\$database, 'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,'design_file=s'=>\$design_file,
	    'array_name=s'=>\$array_name, 'array_type=s'=>\$array_type, 'manufacturer=s'=>\$manufacturer, 'probe_length=i'=>\$probe_length,
	    'order_date=s'=>\$order_date, 'block_size=i'=>\$block_size, 'update_database=s'=>\$update_database);

#Provide instruction to the user
print GREEN, "\n\nThis script takes a NimbleGen design file as input and imports the array spot records contained within it to the specified database", RESET;
print GREEN, "\nNOTE: NimbleGen's internal control probes will not be imported", RESET;
print GREEN, "\nUsage:", RESET;
print GREEN, "\n\tSpecify the target server using: --server", RESET;
print GREEN, "\n\tSpecify the target database using: --database", RESET;
print GREEN, "\n\tSpecify the user and password using: --user and --password", RESET;
print GREEN, "\n\t\tMake sure you chose the correct database and server!!", RESET;
print GREEN, "\n\tSpecify the design file containing array spot records using: --design_file", RESET;
print GREEN, "\n\tSpecify the array name using: --array_name", RESET;
print GREEN, "\n\tSpecify the array type using: --array_type", RESET;
print GREEN, "\n\tSpecify the manufacturer using: --manufacturer", RESET;
print GREEN, "\n\tSpecify the probe length using: --probe_length (in case of isothermal design, specify target length)", RESET;
print GREEN, "\n\tSpecify the date the array design was ordered using: --order_date", RESET;
print GREEN, "\n\tSpecify the block size for import statements using: --block_size", RESET;
print GREEN, "\n\tOnce this script has been tested, use the --update_database option to allow a database update (ex. --update_database=yes)", RESET;
print GREEN, "\n\nExample: importArrayDesign.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --design_file=2006-09-29_MG_Human.ndf  --array_name='MIP101_vs_MIP5FUR'  --array_type=Validation_Isothermal  --manufacturer=NimbleGen  --probe_length=36  --order_date=2006-09-29  --block_size=1000  --update_database=no\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $design_file && $array_name && $array_type && $manufacturer && $probe_length && $order_date && $block_size && $update_database){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

#1.) First create a microarray record - Description of the microarray design etc.
my $microarray_id;
if ($update_database eq "yes"){
  print BLUE, "\nCreating microarray record for: Name=$array_name Type=$array_type Manufacturer=$manufacturer ProbeLength=$probe_length OrderDate=$order_date\n\n", RESET;
  &createMicroarrayRecord('-array_name'=>$array_name, '-array_type'=>$array_type, '-manufacturer'=>$manufacturer, '-probe_length'=>$probe_length,
			  '-order_date'=>$order_date, '-dbh'=>$alexa_dbh);
}

#2.) Import data from the user specified probe file
my $total_insert_count = 0;
&importArraySpotData('-input_file'=>$design_file, '-block_size'=>$block_size, '-dbh'=>$alexa_dbh);

print BLUE, "\n\nInserted a total of $total_insert_count ArraySpot records to the 'ArraySpot' table\n\n", RESET;

unless ($update_database eq "yes"){
  print YELLOW, "\n\nDesign info not imported to database!! - Run script again using -update_database=yes option\n\n", RESET;
}

#Close database connection
$alexa_dbh->disconnect();

exit();


###########################################################################################################
#Create microarray record - a description of the microarray design                                        #
###########################################################################################################
sub createMicroarrayRecord{

  my %args = @_;
  my $name = $args{'-array_name'};
  my $type = $args{'-array_type'}; 
  my $manufacturer = $args{'-manufacturer'};
  my $probe_length = $args{'-probe_length'};
  my $order_date = $args{'-order_date'};
  my $dbh = $args{'-dbh'};

  #1.) See if a microarray of this name has already been entered
  my $sql_test = "SELECT id FROM Microarray WHERE name='$name';";

  my $sth_test = $dbh->prepare("$sql_test");
  $sth_test->execute();
  my ($m_id) = $sth_test->fetchrow_array();
  $sth_test->finish();

  if ($m_id){
    print RED, "\n\nA Microarray of the name: $name is already present in the database!!\n\n", RESET;
    exit();
  }

  #2.) Enter a new microarray record to the database
  #FIELDS for 'Microarray'
  #name, type, manufacturer, probe_length, order_date
  my $sql_insert = "INSERT INTO Microarray (name,type,manufacturer,probe_length,order_date) VALUES (\'$name\',\'$type\',\'$manufacturer\',$probe_length,\'$order_date\');";
  my $sth_insert = $dbh->prepare("$sql_insert");
  $sth_insert->execute();
  $sth_insert->finish();

  #3.) Get the id of the microarray entered
  my $sql_id = "SELECT id FROM Microarray WHERE name=\'$name\';";
  my $sth_id = $dbh->prepare("$sql_id");
  $sth_id->execute();
  $m_id = $sth_id->fetchrow_array();
  $sth_id->finish();

  $microarray_id = $m_id;

  return();
}


###########################################################################################################
#2.) Import ArraySpot Data from an input file to the database                                             #
###########################################################################################################
sub importArraySpotData{
  my %args = @_;
  my $input_file = $args{'-input_file'};
  my $block_size = $args{'-block_size'};
  my $dbh = $args{'-dbh'};

  #Load all probes info from an input probe file and build a probe object as a hash
  my $line_count = 0;
  my $blocks_imported = 0;
  my $probe_count = 0;
  my $total_probe_count = 0;
  my %spots;
  my %spots_nc;

  #Open the probe file and read the neccessary probe data into a hash keyed on probe ID
  open (SPOTS, "$input_file") || die "\nCould not open input design file: $input_file";

  #Possible column headings:
  #PROBE_DESIGN_ID,CONTAINER,DESIGN_NOTE,SELECTION_CRITERIA,SEQ_ID,PROBE_SEQUENCE,MISMATCH,MATCH_INDEX,FEATURE_ID,ROW_NUM,COL_NUM,PROBE_CLASS,PROBE_ID,
  #POSITION,DESIGN_ID,X,Y

  my $first_line = 1;
  my %columns;

  print BLUE, "\n\nParsing input design file.  NimbleGen internal controls will be skipped\n\n", RESET;

  while (<SPOTS>){

    #Get the header line and identify column names and their positions
    my $header_line;
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
      next();
    }
    $line_count++;

    #Get the values of interest from each line (probe record)
    chomp($_);
    my $line_record = $_;
    my @spot_line = split ("\t", $_);

    #Check for critical columns and their names
    unless ($columns{PROBE_ID} && $columns{ROW_NUM} && $columns{COL_NUM} && $columns{X} && $columns{Y}){
      print RED, "\nCritical column missing or named incorrectly, check input file", RESET;
      exit();
    }

    my $probe_val = $spot_line[$columns{PROBE_ID}{column_pos}];

    my $probe_count_id;
    if ($probe_val =~ /^ALEXA\_(\d+)/){
      $probe_count_id = $1;
      print CYAN, "\n\t$line_count. Found probe ID: $probe_count_id - getting probe info", RESET;
    }else{
      print YELLOW, "\n\t$line_count. Invalid probe ID: $probe_val - skipping", RESET;
      next();
    }

    #Using this probe ID (a file probe ID, get the database probe ID, which should normally be the same but not neccessarily)
    my $probe_ref;

    $probe_ref = &getProbeInfo ('-dbh'=>$dbh, '-probe_count_id'=>$probe_count_id);
    my $probe_id = $probe_ref->{db_id};

    #DEBUG
    unless ($probe_id){
      print RED, "\nA probe ID could not be found for Probe_Count_ID: $probe_count_id", RESET;
      exit();
    }

    $probe_count++;
    $total_probe_count++;

    #Store probe in the probe object hash
    $spots{$probe_id}{row_num} = $spot_line[$columns{ROW_NUM}{column_pos}];
    $spots{$probe_id}{col_num} = $spot_line[$columns{COL_NUM}{column_pos}];
    $spots{$probe_id}{x} = $spot_line[$columns{X}{column_pos}];
    $spots{$probe_id}{y} = $spot_line[$columns{Y}{column_pos}];
    $spots{$probe_id}{type} = "PROBE";


    #Importing probe data to the database and reseting variables
    if ($probe_count == $block_size){

      #Once a set of probe data is acquired, dump it to the database and start over again
      if ($update_database eq "yes"){
	$blocks_imported++;
	&insertSpotRecords('-dbh'=>$dbh, '-spot_object'=>\%spots);
	print BLUE, "\n\t\t$blocks_imported: blocks of $block_size imported to database", RESET;
      }
      $probe_count = 0;
      %spots = ();

    }
  }
  close (SPOTS);

  my $final_block_size = keys %spots;

  #Import the remaining spots that did not fill a complete block
  if ($update_database eq "yes" && $final_block_size > 0){
    &insertSpotRecords('-dbh'=>$alexa_dbh, '-spot_object'=>\%spots);
    print BLUE, "\n\t\t$blocks_imported: block of $final_block_size imported to database", RESET;
  }

  #Clean-up and summary
  $blocks_imported++;
  $probe_count = 0;
  %spots = ();

  print BLUE, "\n\nFound a total of $total_probe_count probes in the input file\n\n", RESET;

  return ();
}

########################################################################################################################
#2.) Import probe records to the 'Probe' table                                                                         #
########################################################################################################################
sub insertSpotRecords{
  my %args = @_;
  my $dbh = $args{'-dbh'};
  my %spots = %{$args{'-spot_object'}};

  #2.) Create array spot records

  #FIELDS for 'ArraySpot'
  #fk_Microarray__id,fk_Probe__id,type,x_coord,y_coord,array_row,array_column

  my $spot_count = 0;
  my $total_spots = keys %spots;

  my $sql = "INSERT INTO ArraySpot (fk_Microarray__id,fk_Probe__id,type,x_coord,y_coord,array_row,array_column) VALUES ";

  foreach my $probe_id (sort {$a <=> $b} keys %spots){
    $spot_count++;

    #If a block is complete, or the last probe record has been reached, perform the insert
    if ($spot_count == $total_spots){

      #Build the final part of the insert statement, remember that non-numeric values should be quoted
      $sql = "$sql"."($microarray_id,$probe_id,\'$spots{$probe_id}{type}\',$spots{$probe_id}{x},$spots{$probe_id}{y},$spots{$probe_id}{row_num},$spots{$probe_id}{col_num});";

      #Actually insert the multi-insert statement
      my $sth = $dbh->prepare("$sql");
      $sth->execute();
      $sth->finish();
      $total_insert_count++;

    }else{
      #otherwise, continue building the multi-insert statement
      $sql = "$sql"."($microarray_id,$probe_id,\'$spots{$probe_id}{type}\',$spots{$probe_id}{x},$spots{$probe_id}{y},$spots{$probe_id}{row_num},$spots{$probe_id}{col_num}),";
      $total_insert_count++;
    }
  }
  return();
}

