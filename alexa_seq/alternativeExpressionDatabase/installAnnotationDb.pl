#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to download and install a complete ALEXA-Seq annotation database
#This includes the ALEXA mysql database, ALEXA-Seq annotation files and BLAST and BWA indexed databases
#All of these components are available as a single tarball for each species/build on the GSC FTP server 

#Steps:
#1.) Get the database build, annotation dir, and mysql credentials as parameters from the command line
#2.) Check if the annotation dir already contains this build - if so, exit
#3.) Download the annotation archive from the FTP site (use wget)
#4.) Unpack the annotation archive
#5.) Create the mysql database and import the schema
#6.) Unpack the mysql database files and populate the mysql database
#7.) Delete the mysql database archive

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

#Load the ALEXA modules
my $script_dir;
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
    $script_dir = $1;
  }
}
use utilities::utility qw(:all);

#Initialize command line options
my $annotation_dir = '';
my $db_build = '';
my $server = '';
my $user = '';
my $password = '';

GetOptions ('annotation_dir=s'=>\$annotation_dir, 'db_build=s'=>\$db_build, 'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a the directory to store annotation databases using: --annotation_dir", RESET;
print GREEN, "\n\tSpecify the species/build using:  --db_build", RESET;
print GREEN, "\n\tSpecify mysql server name, user, and passord using:  --server  --user  and  --password", RESET;
print GREEN, "\n\nExample: installAnnotationDb.pl  --annotation_dir=/projects/malachig/sequence_databases/  --db_build=sc_54_1i  --server=localhost  --user=alexa-seq  --password=alexa-seq\n\n", RESET;

#Make sure all options were specified
unless ($annotation_dir && $db_build && $server && $user && $password){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Get the binaries needed
my $mysql_bin = `which mysql 2>/dev/null`;
chomp($mysql_bin);
my $wget_bin = `which wget 2>/dev/null`;
chomp($wget_bin);
my $tar_bin = `which tar 2>/dev/null`;
chomp($tar_bin);

unless ($mysql_bin && $wget_bin && $tar_bin){
  print RED, "\n\nCould not find a needed system binary (mysql, wget, tar)\n\n", RESET;
  exit();
}


#Check the base annotation dir
$annotation_dir = &checkDir('-dir'=>$annotation_dir, '-clear'=>"no");

#global vars
my $build_dir = "$annotation_dir"."$db_build/";
my $alexa_db_name = "ALEXA_"."$db_build";
my $ftp_path = "ftp://ftp03.bcgsc.ca/public/ALEXA/alexa_seq/";
my $archive_name = "$db_build".".tar.gz";
my $alexa_db_dir = "$build_dir"."alexa_db/";
my $alexa_db_archive_name = "ALEXA_"."$db_build".".tables.tar.gz";

#Check current status of the installation...
my $test = &checkInstall();

if ($test){
  print YELLOW, "\n\nALEXA-Seq annotation database is already installed, should be safe to proceed with analysis\n\n", RESET;
  exit();
}

#Enter the working directory
chdir($annotation_dir);


#Download the annotation database 
if (-e "$annotation_dir$archive_name"){
  print YELLOW, "\n\nArchive file is already present, delete this file if you want to re-download it...\n\n", RESET;
}else{
  print BLUE, "\n\nAttempting to download specified ALEXA-Seq annotation database from BCGSC servers ...", RESET;
  my $wget_cmd = "$wget_bin $ftp_path$archive_name";
  print BLUE, "\n1. $wget_cmd\n\n", RESET;
  my $exit_code = system($wget_cmd);
  if ($exit_code){
    #Make second attempt
    print YELLOW, "\n\nwget failed... attempting again", RESET;
    sleep(5);
    print BLUE, "\n2. $wget_cmd\n\n", RESET;
    my $exit_code = system($wget_cmd);
    if ($exit_code){
      print RED, "\n\nDownload failed.  Make sure the 'ALEXA_SEQ_DB' parameter in your config file corresponds to a build available at the FTP server:\n$ftp_path\n\n", RESET;
      exit();
    }
  }
}


#Unpack the archive
my $unpack_cmd = "$tar_bin -zxvf $archive_name 1>/dev/null";
print BLUE, "\n\nUnpacking the downloaded archive file:\n$unpack_cmd\n", RESET;
system($unpack_cmd);


#Make sure the target dir was actually created
$build_dir = &checkDir('-dir'=>$build_dir, '-clear'=>"no");
print BLUE, "\n\nCreated the target database directory", RESET;


#Check if a mysql database with the target name aleady exists
my $db_check_cmd = "echo \"SHOW DATABASES like \'$alexa_db_name\';\" | $mysql_bin -h $server -u $user -p$password";
print BLUE, "\n\nChecking if the mysql database to be created already exists:\n$db_check_cmd", RESET;
my $db_check = `$db_check_cmd`;
if ($db_check){
  print YELLOW, "\n\nIt looks like the annotation mysql database has already been installed: $alexa_db_name", RESET;
  print YELLOW, "\n\n'Drop' the database if you wish to create it cleanly - otherwise proceed to the analysis\n\n", RESET;
}else{
  #Create the mysql database
  my $create_db_cmd = "echo \"CREATE DATABASE $alexa_db_name;\" | $mysql_bin -h $server -u $user -p$password";
  print BLUE, "\n\nCreating the mysql database:\n$create_db_cmd", RESET;
  system($create_db_cmd);

  #Unpack the mysql database files and populate the mysql database
  chdir($alexa_db_dir);
  $unpack_cmd = "$tar_bin -zxvf $alexa_db_archive_name 1>/dev/null";
  print BLUE, "\n\nUnpacking the mysql tables archive file:\n$unpack_cmd\n", RESET;
  system($unpack_cmd);

  my $populate_cmd = "$script_dir/sql/restoreAlexaDb.pl  --database=$alexa_db_name  --server=$server  --user=$user  --password=$password  --mysql_bin=$mysql_bin  --working_dir=$alexa_db_dir  --logfile=restoreAlexaDb_LOG.txt";
  print BLUE, "\n\nPopulating the mysql database with each table:\n$populate_cmd\n", RESET;
  system($populate_cmd);

  #Delete the mysql database tables
  my $rm_cmd2 = "rm -f *.sql";
  print BLUE, "\n\nDeleting the mysql table files:\n$rm_cmd2\n", RESET;
  system($rm_cmd2);
}

#Check for success of installation...
$test = &checkInstall();

if ($test){
  print BLUE, "\n\nALEXA-Seq annotation database seems to have been installed successfully, should be safe to proceed with analysis\n\n", RESET;
}else{
  print RED, "\n\nALEXA-Seq annotation database installation seems to have failed!!\n\nTry droping the mysql database AND/OR deleting the annotation dir and rebuilding:", RESET;
  print RED, "\n\ne.g. use the commands:\n\nDROP DATABASE $alexa_db_name\nrm -fr $build_dir\nrm -f $annotation_dir$archive_name\n\n", RESET;
  exit();
}


#Delete the archive
my $rm_cmd = "rm -f $annotation_dir$archive_name";
print BLUE, "\n\nDeleting the archive file:\n$rm_cmd", RESET;
system($rm_cmd);


print "\n\n";
exit();




################################################################################################################################
#Check for presence of a complete ALEXA-Seq annotation DB                                                                      #
################################################################################################################################
sub checkInstall{

  #Check for all standard directories and also that the mysql database has been created and contains the expected tables
  my $test = 1;

  #Now check for the presence of the requested annotation db within the annotation dir
  unless (-e $build_dir && -d $build_dir){$test = 0;}

  #Check for expected subdirs within this dir
  my $dir1 = "$build_dir"."alexa_db/";
  my $dir2 = "$build_dir"."exonBoundaries/";
  my $dir3 = "$build_dir"."exonJunctions/";
  my $dir4 = "$build_dir"."exonRegions/";
  my $dir5 = "$build_dir"."genes/";
  my $dir6 = "$build_dir"."intergenics/";
  my $dir7 = "$build_dir"."introns/";
  my $dir8 = "$build_dir"."repeats/";
  my $dir9 = "$build_dir"."transcripts/";
  unless (-e $dir1 && -d $dir1){$test = 0;}
  unless (-e $dir2 && -d $dir2){$test = 0;}
  unless (-e $dir3 && -d $dir3){$test = 0;}
  unless (-e $dir4 && -d $dir4){$test = 0;}
  unless (-e $dir5 && -d $dir5){$test = 0;}
  unless (-e $dir6 && -d $dir6){$test = 0;}
  unless (-e $dir7 && -d $dir7){$test = 0;}
  unless (-e $dir8 && -d $dir8){$test = 0;}
  unless (-e $dir9 && -d $dir9){$test = 0;}

  #Check for presence of mysql database and all expected tables
  my $db_check_cmd;
  my $db_check;
  $db_check_cmd = "echo \"SHOW DATABASES like \'$alexa_db_name\';\" | $mysql_bin -h $server -u $user -p$password";
  $db_check = `$db_check_cmd`;
  unless($db_check){$test = 0;}

  if ($db_check){
    $db_check_cmd = "echo \"USE \'$alexa_db_name\'; SHOW TABLES like \'Gene\';\" | $mysql_bin -h $server -u $user -p$password";
    $db_check = `$db_check_cmd`;
    unless($db_check){$test = 0;}

    $db_check_cmd = "echo \"USE \'$alexa_db_name\'; SHOW TABLES like \'Exon\';\" | $mysql_bin -h $server -u $user -p$password";
    $db_check = `$db_check_cmd`;
    unless($db_check){$test = 0;}

    $db_check_cmd = "echo \"USE \'$alexa_db_name\'; SHOW TABLES like \'External_id\';\" | $mysql_bin -h $server -u $user -p$password";
    $db_check = `$db_check_cmd`;
    unless($db_check){$test = 0;}

    $db_check_cmd = "echo \"USE \'$alexa_db_name\'; SHOW TABLES like \'MaskedGene\';\" | $mysql_bin -h $server -u $user -p$password";
    $db_check = `$db_check_cmd`;
    unless($db_check){$test = 0;}

    $db_check_cmd = "echo \"USE \'$alexa_db_name\'; SHOW TABLES like \'Protein_feature\';\" | $mysql_bin -h $server -u $user -p$password";
    $db_check = `$db_check_cmd`;
    unless($db_check){$test = 0;}

    $db_check_cmd = "echo \"USE \'$alexa_db_name\'; SHOW TABLES like \'Transcript\';\" | $mysql_bin -h $server -u $user -p$password";
    $db_check = `$db_check_cmd`;
    unless($db_check){$test = 0;}

    $db_check_cmd = "echo \"USE \'$alexa_db_name\'; SHOW TABLES like \'TranscriptExon\';\" | $mysql_bin -h $server -u $user -p$password";
    $db_check = `$db_check_cmd`;
    unless($db_check){$test = 0;}
  }

  return($test);
}

