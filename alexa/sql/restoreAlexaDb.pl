#!/usr/bin/perl -w
#Written by Malachi Griffith

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use DBI;

#The purpose of this script is to act as a wrapper for mysql commandline utilities
#It facilitates the restoration of a specified mysql database on a table by table basis
#Specify the location of a folder containing .sql files for an ALEXA database

#Because the ALEXA database schema is relational and has foreign key constraints, the database tables
#must be imported in a particular order.

#mysql utilities used: mysql

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $mysql_bin = '';
my $working_dir = '';
my $logfile = '';

#Hardcode a list of table names in the order they must be populated
my @tables = qw (Gene Transcript Protein_feature External_id Exon TranscriptExon Probe Probe_set ProbeProbe_set GeneProbe Microarray ArraySpot MaskedGene all_mrna all_est);

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'mysql_bin=s'=>\$mysql_bin, 'working_dir=s'=>\$working_dir, 'logfile=s'=>\$logfile);


#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tNote: Before running this script, log into your mysql server and create the target database!!", RESET;
print GREEN, "\n\t\ti.e.  CREATE DATABASE ALEXA_test;", RESET;
print GREEN, "\n\tSpecify the full path to the desired mysql client binary using: --mysql_bin", RESET;
print GREEN, "\n\tSpecify the full path to a directory containing ALEXA .sql table files using:  --working_dir", RESET;
print GREEN, "\n\t\tNote: You may need to unzip and unpack the archive first!", RESET;
print GREEN, "\n\t\ti.e.  gunzip ALEXA_version.table.tar.gz ... THEN ... tar -xvf ALEXA_version.table.tar", RESET;
print GREEN, "\n\tSpecify a logfile for output using: --logfile", RESET;
print GREEN, "\nExample: restoreAlexaDb.pl  --database=ALEXA_test  --server=server_name  --user=mysql_user  --password=mysql_pwd  --mysql_bin=/usr/bin/mysql  --working_dir=/home/user/alexa/ALEXA_version/database_backup/  --logfile=/home/user/alexa/ALEXA_version/logs/database_population/restoreAlexaDb_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $mysql_bin && $working_dir && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";
print LOG "\nUser Specified the following options:\ndatabase = $database\nserver = $server\nuser = $user\npassword = pwd\nmysql_bin = $mysql_bin\nworking_dir = $working_dir\nlogfile = $logfile\n\n";

#Make sure the provided working directory is valid
unless ($working_dir =~ /.*\/$/){
  $working_dir = "$working_dir"."/";
}

unless (-e $working_dir && -d $working_dir){
  print RED, "\nSpecified directory does not appear to be valid: $working_dir\n\n", RESET;
  exit();
}

unless (-e $mysql_bin){
  print RED, "\nSpecified mysql binary: $mysql_bin does not appear to exist!\n\n", RESET;
  exit();
}


#Establish connection with the Alternative Splicing Expression database using details provided by user
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

#Make sure the max allowed packet size is large enough to allow large genes
&checkMaxPacketSize('-dbh'=>$alexa_dbh);

#Get the version of mysqldump used to create this backup and store this info in the log file
my $version_command = "$mysql_bin -V";
my $version_info = `$version_command`;
chomp($version_info);

print BLUE, "\nThe following database restoration will be done with mysql client version: $version_info\n\n", RESET;
print LOG "\nThe following database restoration will be done with mysql client version: $version_info\n\n";

#Get a complete list of .sql files in the working directory
print BLUE, "\nGathering .sql file from the specified working directory: $working_dir\n\n", RESET;
print LOG "\nGathering .sql file from the specified working directory: $working_dir\n\n";

my %files;
opendir(DIRHANDLE, "$working_dir") || die "\nCannot open directory: $working_dir\n\n";
my @test_files = readdir(DIRHANDLE);

foreach my $test_file (@test_files){
  my $file_path = "$working_dir"."$test_file";

  #Skip directories within the specified directory
  if (-e $file_path && -d $file_path){
    print YELLOW, "\n\t$file_path is a directory - skipping", RESET;
    next();
  }
  #Skip files that do not have a '.sql' extension
  my $table_name;
  if ($test_file =~ /(.*)\.sql$/){
    $table_name = $1;
  }else{
    print YELLOW, "\n\t$test_file does not have a .sql extension - skipping", RESET;
    next();
  }
  $files{$test_file}{file_path} = $file_path;
  $files{$test_file}{table_name} = $table_name;
}

#Make sure the number of files matches the expected number of tables
my $table_count = @tables;
my $file_count = keys %files;
unless ($table_count == $file_count){
  print RED, "\nNumber of files found does not match expected number of tables - exiting\n\n", RESET;
  exit();
}

#Make sure all the table names are as expected
print BLUE, "\n\nSearching for expected table name\n", RESET;
print LOG "\n\nSearching for expected table name\n";
foreach my $table (@tables){
  my $table_found = 0;
  foreach my $file (keys %files){
    my $file_table = $files{$file}{table_name};

    if ($file_table eq $table){
      $table_found = 1;
      print BLUE, "\n\tFound table: $table", RESET;
      print LOG "\n\tFound table: $table";
    }
  }
  if ($table_found == 0){
    print RED, "\nCould not find table: $table!!\n\n", RESET;
    exit();
  }
}

#Do the actual import statements
print BLUE, "\n\nImporting table with: $mysql_bin\n", RESET;
print LOG "\n\nImporting table with: $mysql_bin\n";
TABLE:foreach my $table (@tables){

  foreach my $file (keys %files){
    my $file_table = $files{$file}{table_name};
    if ($file_table eq $table){

      my $file_path = $files{$file}{file_path};

      my $import_command = "$mysql_bin --host=$server --user=$user --password=$password --database=$database < $file_path";

      print BLUE, "\n\t$import_command", RESET;
      print LOG "\n\t$import_command";
      system($import_command);
      next TABLE;
    }
  }
}

print "\n\n";
print LOG "\n\n";

close (LOG);
exit();

###############################################################################################################
#Create mysql database connection                                                                             #
###############################################################################################################
sub connectDB {
  my %args = @_;
  my $database_name = $args{'-database'};
  my $database_host = $args{'-server'};
  my $user_name = $args{'-user'};
  my $user_pw = $args{'-password'};

  my $dbh = DBI->connect( "dbi:mysql:database=$database_name;host=$database_host", $user_name, $user_pw, { PrintError => 1 } );
  return $dbh;
}


#############################################################################################
#checkMaxPacketSize                                                                         #
#############################################################################################
sub checkMaxPacketSize{
  my %args = @_;
  my $dbh = $args{'-dbh'};

  my $sql = "SHOW GLOBAL variables LIKE 'max_allowed_packet'";
  my $sth = $dbh->prepare("$sql");
  $sth->execute();
  my $packet_size = $sth->fetchrow_array();
  unless ($packet_size > 5000000){
    print RED, "\nmax_allowed_packet = $packet_size bytes", RESET;
    print RED, "\nWARNING: max_allowed_packet variable is too small to allow insertion of large genes", RESET;
    print RED, "\nGet a DB Admin to execute the following: SET GLOBAL max_allowed_packet=10000000\n\n", RESET;
    exit();
  }

  print BLUE, "\nFound a suitable packet size ($packet_size) for the target ALEXA database\n\n", RESET;
  $sth->finish();
  return();
}


