#!/usr/bin/perl -w
#Written by Malachi Griffith

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#The purpose of this script is to act as a wrapper for mysql commandline utilities
#It facilitates the backup of a specified mysql database on a table by table basis

#Because the ALEXA database schema is relational and has foreign key constraints, the database tables
#must be created in a particular order.

#mysql utilities used: mysqldump
#unix utilities used: tar and gzip

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $mysql_dump_bin = '';
my $working_dir = '';
my $logfile = '';

#Hardcode a list of table names in the order they must be populated
my @tables = qw (Gene Transcript Protein_feature External_id Exon TranscriptExon Probe Probe_set ProbeProbe_set GeneProbe Microarray ArraySpot MaskedGene Probe_negativeControl);
my @tables_no_data = qw (all_mrna all_est);

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'mysql_dump_bin=s'=>\$mysql_dump_bin, 'working_dir=s'=>\$working_dir, 'logfile=s'=>\$logfile);


#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the full path to the desired mysqldump binary using: --mysql_dump_bin", RESET;
print GREEN, "\n\tSpecify the full path to target directory for backup file using:  --working_dir", RESET;
print GREEN, "\n\tSpecify a logfile for output using: --logfile", RESET;
print GREEN, "\n\nExample: backupAlexaDb.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=mysql_user  --password=mysql_pwd  --mysql_dump_bin=/usr/bin/mysqldump  --working_dir=/home/user/alexa/ALEXA_version/database_backup/  --logfile=/home/user/alexa/ALEXA_version/logs/database_population/backAlexaDb_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $mysql_dump_bin && $working_dir && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";
print LOG "\nUser Specified the following options:\ndatabase = $database\nserver = $server\nuser = $user\npassword = pwd\nmysql_dump_bin = $mysql_dump_bin\nworking_dir = $working_dir\nlogfile = $logfile\n\n";

#Make sure the provided working directory is valid
unless ($working_dir =~ /.*\/$/){
  $working_dir = "$working_dir"."/";
}

unless (-e $working_dir && -d $working_dir){
  print RED, "\nSpecified directory does not appear to be valid: $working_dir\n\n", RESET;
  exit();
}

unless (-e $mysql_dump_bin){
  print RED, "\nSpecified mysqldump binary: $mysql_dump_bin does not appear to exist!\n\n", RESET;
  exit();
}

#Get the version of mysqldump used to create this backup and store this info in the log file
my $version_command = "$mysql_dump_bin -V";
my $version_info = `$version_command`;
chomp($version_info);

print BLUE, "\nThe following database backup will be created with mysqldump version: $version_info\n\n", RESET;
print LOG "\nThe following database backup will be created with mysqldump version: $version_info\n\n";

#Go through each mysql table name (in order) and use mysqldump to create a backup in the working directory
print BLUE, "\nCreating table dumps with: $mysql_dump_bin\n\n", RESET;
print LOG "\nCreating table dumps with: $mysql_dump_bin\n\n";

#Dump table structure + data
foreach my $table (@tables){

  print BLUE, "\n\nProcessing table: $table", RESET;
  print LOG "\n\nProcessing table: $table";

  my $file = "$working_dir"."$table".".sql";

  my $dump_command = "$mysql_dump_bin --user=$user --password=$password --result-file=$file $database $table";
  print YELLOW, "\n\t$mysql_dump_bin --user=$user --password=$password --result-file=$file $database $table", RESET;
  print LOG "\n\t$mysql_dump_bin --user=$user --password=pwd --result-file=$file $database $table";

  system($dump_command);
}

#Dump table structure only
foreach my $table (@tables_no_data){

  print BLUE, "\n\nProcessing table: $table", RESET;
  print LOG "\n\nProcessing table: $table";

  my $file = "$working_dir"."$table".".sql";

  my $dump_command = "$mysql_dump_bin --user=$user --password=$password --no-data --result-file=$file $database $table";
  print YELLOW, "\n\t$mysql_dump_bin --user=$user --password=$password --no-data --result-file=$file $database $table", RESET;
  print LOG "\n\t$mysql_dump_bin --user=$user --password=pwd --no-data --result-file=$file $database $table";

  system($dump_command);
}


#Package all of these files into a tarball
my $tarball = "$working_dir"."$database".".tables.tar";
my $table_path = "*.sql";
my $tar_command = "tar -cf $tarball $table_path";
chdir($working_dir);

print BLUE, "\n\nCreating tarball of .sql table files", RESET;
print YELLOW, "\n\n$tar_command\n\n", RESET;
print LOG "\n\nCreating tarball of .sql table files";
print LOG "\n\n$tar_command\n\n";
system($tar_command);

#Compress the tarball
my $compress_command = "gzip -9 $tarball";
print BLUE, "\nCompressing ...", RESET;
print YELLOW, "\n$compress_command\n\n", RESET;
print LOG "\nCompressing ...";
print LOG "\n$compress_command\n\n";
system($compress_command);


#Cleanup the individual files
my $rm_path = "$working_dir"."*.sql";
my $rm_command = "rm -f $rm_path";
print BLUE, "\nRemoving individual table files ...", RESET;
print YELLOW, "\n$rm_command\n\n", RESET;
print LOG "\nRemoving individual table files ...";
print LOG "\n$rm_command\n\n";
system($rm_command);

close (LOG);

exit();




