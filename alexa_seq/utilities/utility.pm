=head1 NAME

utility.pm - library modules that contains generic utilities

=head1 SYNOPSIS

use utility qw(:all);

=head2 NOTE

currently located in '~/utilities'

=head2 RECENT CHANGES

None.  Last modified 21 December 2006

=head1 DESCRIPTION

Generic utility

=head1 EXAMPLES

use lib './';

use utilities::utility qw(:all);

=head1 SEE ALSO

None

=head1 BUGS

Contact author via email

=head1 AUTHOR

Written by Malachi Griffith (malachig@bcgsc.ca)

=head1 ACKNOWLEDGEMENTS

University of British Columbia Graduate Studies

Michael Smith Foundation for Health Research

Natural Sciences and Engineering Research Council

Genome British Columbia

=head1 AFFLIATIONS

Malachi Griffith is supervised by Marco A. Marra

Genome Sciences Centre, BC Cancer Research Centre, BC Cancer Agency, UBC Faculty of Medicine - Medical Genetics

=head1 SUBROUTINES

=cut

#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

package utilities::utility;
require Exporter;

@ISA = qw( Exporter );
@EXPORT = qw();

@EXPORT_OK = qw(&checkMaxPacketSize &connectDB &loadEnsemblApi &createNewDir &checkDir &commify &checkConfig &memoryUsage);

%EXPORT_TAGS = (
     all => [qw(&checkMaxPacketSize &connectDB &loadEnsemblApi &createNewDir &checkDir &commify &checkConfig &memoryUsage)]
);

use strict;
use Data::Dumper;
use DBI;
use Term::ANSIColor qw(:constants);


#Get path to EnsEMBL API installations from an environment variable
my $ensembl_api_path;
BEGIN{
  $ensembl_api_path=$ENV{"ENSEMBL_API_PATH"};
  unless($ensembl_api_path){
    print RED, "\n\nEnsembl_api_path undefined - set the environment variable: ENSEMBL_API_PATH\n\n", RESET;
    exit();
  }
  unless(-e $ensembl_api_path){
    print RED, "\n\nCould not find ensembl_api_path: ($ensembl_api_path) - correct the environment variable: ENSEMBL_API_PATH\n\n", RESET;
    exit();
  }
}
=head2 checkMaxPacketSize()

=over 3

=item Function:

Check the database max_allowed_packet variable to ensure that large gene entries will work

=item Return:

N/A - Will exit if packet size is not sufficient

=item Args:

'-dbh' => database handle

=item Example(s):

&checkMaxPacketSize('-dbh'=>$alexa_dbh);

=back

=cut

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


=head2 connectDB

=over 3

=item Function:

Get database connection handler

=item Return:

database handle

=item Args:

'-database' - mysql database name

'-server' - mysql host

'-user' - mysql user name

'-password' - mysql password

=item Example(s):

$dbh = &connectDB('-database'=>'database_name', '-server'=>'server_name', '-user'=>'user_name', '-password'=>'passwd')

=back

=cut

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
  my $connection_attempts = 1;

  #Check for failures to reconnect... if the connection failed, try again a few times
  my $connection_attempt_limit = 3;
  my $sleep_time = 15;
  while(!(defined($dbh)) && ($connection_attempts <= $connection_attempt_limit)){
    $connection_attempts++;
    print YELLOW, "\n\nDBH connect failed ... sleeping for $sleep_time and try attempt: $connection_attempts", RESET;
    sleep $sleep_time;
    $dbh = DBI->connect( "dbi:mysql:database=$database_name;host=$database_host", $user_name, $user_pw, { PrintError => 1 } );
  }

  return $dbh;
}


=head2 loadEnsemblApi

=over 3

=item Function:

Load API code for specified EnsEBML veriosn

=item Return:

NULL

=item Args:

'-api' - EnsEMBL API version.  e.g. '49'

=item Example(s):

loadEnsemblApi('-api'=>$ensembl_api_version)

=back

=cut


##############################################################################################################
#Load EnsEMBL API library code                                                                               #
##############################################################################################################
sub loadEnsemblApi{
  my %args = @_;
  my $ensembl_api_version = $args{'-api'};


  #**********************************************************************************************************
  #IMPORTANT NOTE: You must have the correct Ensembl API installed locally AND bioperl 1.2 or greater!!
  #Both the EnsEMBL core API as well as Compara are required
  #Refer to the ALEXA manual for additional details on how to install these
  #Then update the following paths:
  unless ($ensembl_api_path =~ /.*\/$/){
    $ensembl_api_path = "$ensembl_api_path"."/";
  }

  my $path1 = "$ensembl_api_path"."ensembl_"."$ensembl_api_version"."_perl_API/ensembl/modules"; 
  my $path2 = "$ensembl_api_path"."ensembl_"."$ensembl_api_version"."_perl_API/ensembl-variation/modules";
  my $path3 = "$ensembl_api_path"."bioperl-1.4";

  unless (-e $path1 && -e $path2 && -e $path3){
    print RED, "\n\nutility.pm could not find one of the components needed for the EnsEMBL API:", RESET;
    print RED, "\n$path1\n$path2\n$path3", RESET;
    print RED, "\n\nRun alternativeExpressionDatabase/installEnsemblAPI.pl\n\n", RESET;
    exit();
  }


  unshift(@INC, "$path1");
  unshift(@INC, "$path2");
  unshift(@INC, "$path3");
  require Bio::EnsEMBL::DBSQL::DBAdaptor; #Used for local connections to EnsEMBL core databases
  require Bio::EnsEMBL::Variation::DBSQL::DBAdaptor; #Used for local connections to EnsEMBL variation databases

  return();
}



=head2 createNewDir

=over 3

=item Function:

Create a new directory cleanly in the specified location - Prompt user for confirmation

=item Return:

Full path to new directory

=item Args:

'-path' - Full path to new directoy

'-new_dir_name' - Name of new directory

=item Example(s):

my $fasta_dir = &createNewDir('-path'=>$temp_dir, '-new_dir_name'=>"ensembl_genes_fasta");

=back

=cut

###############################################################################################################
#Create a new directory in a specified location                                                               #
###############################################################################################################
sub createNewDir{
  my %args = @_;
  my $base_path = $args{'-path'};
  my $name = $args{'-new_dir_name'};
  my $force = $args{'-force'};

  #Now make sure the desired new dir does not already exist
  unless ($base_path =~ /.*\/$/){
    $base_path = "$base_path"."/";
  }

  #First make sure the specified base path exists and is a directory
  unless (-e $base_path && -d $base_path){
    print RED, "\nSpecified working directory: $base_path does not appear valid! Create a working directory before proceeding\n\n", RESET;
    exit();
  }

  unless ($name =~ /.*\/$/){
    $name = "$name"."/";
  }

  my $new_path = "$base_path"."$name";

  if (-e $new_path && -d $new_path){

    if ($force){
      #If this directory already exists, and the -force option was provide, delete this directory and start it cleanly
      if ($force eq "yes"){
	print YELLOW, "\nForcing clean creation of $new_path\n\n", RESET;
	my $command = "rm -r $new_path";
	system ($command);
	mkdir($new_path);
      }else{
	print RED, "\nThe '-force' option provided to utility.pm was not understood!!", RESET;
	exit();
      }

    }else{

      #If this directory already exists, ask the user if they wish to erase it and start clean
      print YELLOW, "\nNew dir: $new_path already exists.\n\tDo you wish to delete it and create it cleanly (y/n)? ", RESET;
      my $answer = <>;

      chomp($answer);

      if ($answer =~ /^y$/i | $answer =~ /^yes$/i){
	my $command = "rm -r $new_path";
	system ($command);
	mkdir($new_path);
      }else{
	print YELLOW, "\nUsing existing directory, some files may be over-written and others that are unrelated to the current analysis may remain!\n", RESET;
      }
    }

  }else{
    mkdir($new_path)
  }
  return($new_path);
}


=head2 checkDir

=over 3

=item Function:

Check validity of a directory and empty if the user desires - Prompt user for confirmation

=item Return:

Path to clean,valid directory

=item Args:

'-dir' - Full path to directory to be checked

'-clear' - 'yes/no' option to clear the specified directory of files

'-force' - 'yes/no' force clear without user prompt

=item Example(s):

my $working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"yes");

=back

=cut


#############################################################################################################################
#Check dir
#############################################################################################################################
sub checkDir{
  my %args = @_;
  my $dir = $args{'-dir'};
  my $clear = $args{'-clear'};
  my $force = $args{'-force'};
  my $recursive = $args{'-recursive'};

  unless ($dir =~ /\/$/){
    $dir = "$dir"."/";
  }
  unless (-e $dir && -d $dir){
    print RED, "\nDirectory: $dir does not appear to be valid!\n\n", RESET;
    exit();
  }

  unless ($force){
    $force = "no";
  }
  unless ($clear){
    $clear = "no";
  }
  unless ($recursive){
    $recursive = "no";
  }

  #Clean up the working directory
  opendir(DIRHANDLE, "$dir") || die "\nCannot open directory: $dir\n\n";
  my @temp = readdir(DIRHANDLE);
  closedir(DIRHANDLE);

  if ($clear =~ /y|yes/i){

    if ($force =~ /y|yes/i){
      if ($recursive =~ /y|yes/i){
        my $files_present = scalar(@temp) - 2;
        my $clean_dir_cmd = "rm -fr $dir"."*";
        print YELLOW, "\n\n$clean_dir_cmd\n\n", RESET;
        system($clean_dir_cmd);
      }else{
        my $files_present = scalar(@temp) - 2;
        my $clean_dir_cmd = "rm -f $dir"."*";
        print YELLOW, "\n\n$clean_dir_cmd\n\n", RESET;
        system($clean_dir_cmd);
      }
    }else{

      my $files_present = scalar(@temp) - 2;
      my $clean_dir_cmd = "rm $dir"."*";
      if ($recursive =~ /y|yes/i){
        $clean_dir_cmd = "rm -fr $dir"."*";
      }

      unless ($files_present == 0){
	print YELLOW, "\nFound $files_present files in the specified directory ($dir)\nThis directory will be cleaned with the command:\n\t$clean_dir_cmd\n\nProceed (y/n)? ", RESET;
	my $answer = <>;
	chomp($answer);
	if ($answer =~ /y|yes/i){
          if ($recursive =~ /y|yes/i){
            system($clean_dir_cmd);
          }else{
	    system($clean_dir_cmd);
          }
	}else{
	  print YELLOW, "\nContinuing and leaving files in place then ...\n\n", RESET;
	}
      }
    }
  }
  return($dir);
}


=head2 checkConfig

=over 3

=item Function:

Check validity of configuration values iimported from system and project specific config files

=item Return:

Warnings for potential problems with config values

=item Args:

'-alexa_seq_conf' = Object containing system-wide configuration values

'-project_conf' - Object containing project specific configuration values

=item Example(s):

my $warnings = &checkConfig('-alexa_seq_conf'=>\%alexa_seq_conf, '-project_conf'=>\%project_conf);

=back

=cut


#############################################################################################################################
#Check config files - and objects created from them                                                                         #
#############################################################################################################################
sub checkConfig{
  my %args = @_;
  my $alexa_conf_ref = $args{'-alexa_seq_conf'};
  my $project_conf_ref = $args{'-project_conf'};
  my $annotation_conf_ref = $args{'-annotation_conf'};
  my $warnings = 0;
  my $wstring = '';


  #####################################################################################################################################
  #1.) First check all the alexa_seq_conf values

  if ($alexa_conf_ref){
    #SHELL_FILE
    if ($alexa_conf_ref->{SHELL_FILE}){
      unless (-e $alexa_conf_ref->{SHELL_FILE}){
        $warnings++;
        $wstring .= "\nCould not find shell file.  Check ALEXA_Seq config file parameter: SHELL_FILE";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: SHELL_FILE";
    }
    #ANALYSIS_DIR
    if ($alexa_conf_ref->{ANALYSIS_DIR}){
      unless (-e $alexa_conf_ref->{ANALYSIS_DIR} && -d $alexa_conf_ref->{ANALYSIS_DIR}){
        $warnings++;
        $wstring .= "\nAnalysis directory does not seem valid.  Check ALEXA_Seq config file parameter: ANALYSIS_DIR";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: ANALYSIS_DIR";
    }

    #SEQUENCE_DB_DIR
    if ($alexa_conf_ref->{SEQUENCE_DB_DIR}){
      unless (-e $alexa_conf_ref->{SEQUENCE_DB_DIR} && -d $alexa_conf_ref->{SEQUENCE_DB_DIR}){
        $warnings++;
        $wstring .= "\nSequence database directory does not seem valid.  Check ALEXA_Seq config file parameter: SEQUENCE_DB_DIR";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: SEQUENCE_DB_DIR";
    }

    #WEB_DIR
    if ($alexa_conf_ref->{WEB_DIR}){
      unless (-e $alexa_conf_ref->{WEB_DIR} && -d $alexa_conf_ref->{WEB_DIR}){
        $warnings++;
        $wstring .= "\nWeb directory does not seem valid.  Check ALEXA_Seq config file parameter: WEB_DIR";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: WEB_DIR";
    }

    #WEB_URL - difficult to check validity (is possible with checkLinks.pl script that is publicly available...)
    if ($alexa_conf_ref->{WEB_URL}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: WEB_URL";
    }

    #ALEXA_HOME_URL
    if ($alexa_conf_ref->{ALEXA_HOME_URL}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: ALEXA_HOME_URL";
    }

    #ALEXA_SEQ_URL - difficult to check validity (is possible with checkLinks.pl script that is publicly available...)
    if ($alexa_conf_ref->{ALEXA_SEQ_URL}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: ALEXA_SEQ_URL";
    }

    #XAPIAN_BIN
    if ($alexa_conf_ref->{XAPIAN_BIN}){
      unless (-e $alexa_conf_ref->{XAPIAN_BIN}){
        $warnings++;
        $wstring .= "\nXapian-Omega binary not found.  Check ALEXA_Seq config file parameter: XAPIAN_BIN";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: XAPIAN_BIN";
    }

    #XAPIAN_URL - difficult to check validity (is possible with checkLinks.pl script that is publicly available...)
    if ($alexa_conf_ref->{XAPIAN_URL}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: XAPIAN_URL";
    }

    #XAPIAN_DIR
    if ($alexa_conf_ref->{XAPIAN_DIR}){
      unless (-e $alexa_conf_ref->{XAPIAN_DIR} && -d $alexa_conf_ref->{XAPIAN_DIR}){
        $warnings++;
        $wstring .= "\nXapian-Omega directory does not seem valid.  Check ALEXA_Seq config file parameter: XAPIAN_DIR";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: XAPIAN_DIR";
    }

    #GOOGLE_ANALYTICS_ID
    if ($alexa_conf_ref->{GOOGLE_ANALYTICS_ID}){
      unless($alexa_conf_ref->{GOOGLE_ANALYTICS_ID} =~ /UA\-\d+\-\d+/ || $alexa_conf_ref->{GOOGLE_ANALYTICS_ID} =~ /UA\-[x]+\-[x]+/){
        $warnings++;
        $wstring .= "\nFormat of Google Analytics ID not understood.  Expecting: UA-xxxxxxx-x (where x's are integers)";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: GOOGLE_ANALYTICS_ID";
    }

    #ALEXA mySQL credentials - difficult to authenticate.
    #ALEXA_SERVER
    if ($alexa_conf_ref->{ALEXA_SERVER}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: ALEXA_SERVER";
    }

    #ALEXA_USER1
    if ($alexa_conf_ref->{ALEXA_USER1}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: ALEXA_USER1";
    }

    #ALEXA_PASSWORD1
    if ($alexa_conf_ref->{ALEXA_PASSWORD1}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: ALEXA_PASSWORD1";
    }

    #ALEXA_USER2
    if ($alexa_conf_ref->{ALEXA_USER2}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: ALEXA_USER2";
    }

    #ALEXA_PASSWORD2
    if ($alexa_conf_ref->{ALEXA_PASSWORD2}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: ALEXA_PASSWORD2";
    }

    #ENSEMBL mySQL credentials - difficult to authenticate.
    #ENSEMBL_SERVER
    if ($alexa_conf_ref->{ENSEMBL_SERVER}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: ENSEMBL_SERVER";
    }

    #ENSEMBL_USER
    if ($alexa_conf_ref->{ENSEMBL_USER}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: ENSEMBL_USER";
    }

    #ENSEMBL_PASSWORD
    if ($alexa_conf_ref->{ENSEMBL_PASSWORD}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: ENSEMBL_PASSWORD";
    }

    #MDUST_BIN
    if ($alexa_conf_ref->{MDUST_BIN}){
      unless (-e $alexa_conf_ref->{MDUST_BIN}){
        $warnings++;
        $wstring .= "\nmdust binary not found.  Check ALEXA_Seq config file parameter: MDUST_BIN";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: MDUST_BIN";
    }

    #BLAST_BIN_DIR
    if ($alexa_conf_ref->{BLAST_BIN_DIR}){
      unless (-e $alexa_conf_ref->{BLAST_BIN_DIR} && -d $alexa_conf_ref->{BLAST_BIN_DIR}){
        $warnings++;
        $wstring .= "\nBLAST binary directory does not seem valid.  Check ALEXA_Seq config file parameter: BLAST_BIN_DIR";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: BLAST_BIN_DIR";
    }

    #BWA_BIN_DIR
    if ($alexa_conf_ref->{BWA_BIN_DIR}){
      unless (-e $alexa_conf_ref->{BWA_BIN_DIR} && -d $alexa_conf_ref->{BWA_BIN_DIR}){
        $warnings++;
        $wstring .= "\nBWA binary directory does not seem valid.  Check ALEXA_Seq config file parameter: BWA_BIN_DIR";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: BWA_BIN_DIR";
    }

    #CLUSTER_HEAD_NODE - difficult to authenticate (ping possibly)
    if ($alexa_conf_ref->{CLUSTER_HEAD_NODE}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: CLUSTER_HEAD_NODE";
    }

    #JOB_SIZE
    if ($alexa_conf_ref->{JOB_SIZE}){
      unless ($alexa_conf_ref->{JOB_SIZE} =~/^\d+$/){
        $warnings++;
        $wstring .= "\nInteger format expected but not found.  Check ALEXA_Seq config file parameter: JOB_SIZE";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: JOB_SIZE";
    }

    #SHORT_DELAY
    if ($alexa_conf_ref->{SHORT_DELAY}){
      unless ($alexa_conf_ref->{SHORT_DELAY} =~/^\d+$/){
        $warnings++;
        $wstring .= "\nInteger format expected but not found.  Check ALEXA_Seq config file parameter: SHORT_DELAY";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: SHORT_DELAY";
    }

    #MEDIUM_DELAY
    if ($alexa_conf_ref->{MEDIUM_DELAY}){
      unless ($alexa_conf_ref->{MEDIUM_DELAY} =~/^\d+$/){
        $warnings++;
        $wstring .= "\nInteger format expected but not found.  Check ALEXA_Seq config file parameter: MEDIUM_DELAY";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: MEDIUM_DELAY";
    }

    #LONG_DELAY
    if ($alexa_conf_ref->{LONG_DELAY}){
      unless ($alexa_conf_ref->{LONG_DELAY} =~/^\d+$/){
        $warnings++;
        $wstring .= "\nInteger format expected but not found.  Check ALEXA_Seq config file parameter: LONG_DELAY";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: LONG_DELAY";
    }

    #VLONG_DELAY
    if ($alexa_conf_ref->{VLONG_DELAY}){
      unless ($alexa_conf_ref->{VLONG_DELAY} =~/^\d+$/){
        $warnings++;
        $wstring .= "\nInteger format expected but not found.  Check ALEXA_Seq config file parameter: VLONG_DELAY";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add ALEXA_Seq config file parameter: VLONG_DELAY";
    }
  }


  #####################################################################################################################################
  #2.) Next check all the project_conf values
  if ($project_conf_ref){
    #PROJECT_NAME
    if ($project_conf_ref->{PROJECT_NAME}){
      if ($project_conf_ref->{PROJECT_NAME} =~ /\s+/){
        $warnings++;
        $wstring .= "\nParameter should not have spaces: Check Project config file parameter: PROJECT_NAME";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Project config file parameter: PROJECT_NAME";
    }

    #PROJECT_DESCRIPTION
    if ($project_conf_ref->{PROJECT_DESCRIPTION}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Project config file parameter: PROJECT_DESCRIPTION";
    }

    #ENSEMBL_VERSION
    if ($project_conf_ref->{ENSEMBL_VERSION}){
      unless ($project_conf_ref->{ENSEMBL_VERSION} =~/^\d+$/){
        $warnings++;
        $wstring .= "\nInteger format expected but not found.  Check Project config file parameter: ENSEMBL_VERSION";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Project config file parameter: ENSEMBL_VERSION";
    }

    #SPECIES_NAME
    if ($project_conf_ref->{SPECIES_NAME}){
      if ($project_conf_ref->{SPECIES_NAME} =~ /\s+/){
        $warnings++;
        $wstring .= "\nParameter should not have spaces: Check Project config file parameter: SPECIES_NAME";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Project config file parameter: SPECIES_NAME";
    }

    #SPECIES_NAME_COMMON
    if ($project_conf_ref->{SPECIES_NAME_COMMON}){
      if ($project_conf_ref->{SPECIES_NAME_COMMON} =~ /\s+/){
        $warnings++;
        $wstring .= "\nParameter should not have spaces: Check Project config file parameter: SPECIES_NAME_COMMON";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Project config file parameter: SPECIES_NAME_COMMON";
    }

    #UCSC_BUILD
    if ($project_conf_ref->{UCSC_BUILD}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Project config file parameter: UCSC_BUILD";
    }

    #ALEXA_SEQ_DB
    if ($project_conf_ref->{ALEXA_SEQ_DB}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Project config file parameter: ALEXA_SEQ_DB";
    }

    #ALIGNMENT OPTION
    if ($project_conf_ref->{ALIGNMENT_OPTION}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Project config file parameter: ALIGNMENT_OPTION";
    }

    #REGIONS_FILE_VERSION
    if ($project_conf_ref->{REGIONS_FILE_VERSION}){
      unless ($project_conf_ref->{REGIONS_FILE_VERSION} =~/^\d+$/){
        $warnings++;
        $wstring .= "\nInteger format expected but not found.  Check Project config file parameter: REGIONS_FILE_VERSION";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Project config file parameter: REGIONS_FILE_VERSION";
    }

    #JUNCTION_SEQ_LENGTH
    if ($project_conf_ref->{JUNCTION_SEQ_LENGTH}){
      unless ($project_conf_ref->{JUNCTION_SEQ_LENGTH} =~/^\d+$/){
        $warnings++;
        $wstring .= "\nInteger format expected but not found.  Check Project config file parameter: JUNCTION_SEQ_LENGTH";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Project config file parameter: JUNCTION_SEQ_LENGTH";
    }

    #BOUNDARY_SEQ_LENGTH
    if ($project_conf_ref->{BOUNDARY_SEQ_LENGTH}){
      unless ($project_conf_ref->{BOUNDARY_SEQ_LENGTH} =~/^\d+$/){
        $warnings++;
        $wstring .= "\nInteger format expected but not found.  Check Project config file parameter: BOUNDARY_SEQ_LENGTH";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Project config file parameter: BOUNDARY_SEQ_LENGTH";
    }

    #MIN_OVERLAP
    if ($project_conf_ref->{MIN_OVERLAP}){
      unless ($project_conf_ref->{MIN_OVERLAP} =~/^\d+$/){
        $warnings++;
        $wstring .= "\nInteger format expected but not found.  Check Project config file parameter: MIN_OVERLAP";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Project config file parameter: MIN_OVERLAP";
    }

    #WEBSITE VERSION
    if ($project_conf_ref->{WEBSITE_VERSION}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Project config file parameter: WEBSITE_VERSION";
    }



    #Check format of LANE, LIBRARY, COMPARISON, GROUP, and GRP_COMPARISON entries
    my $lanes_ref = $project_conf_ref->{LANE};
    my $lib_ref = $project_conf_ref->{LIBRARY};
    my $comp_ref = $project_conf_ref->{COMPARISON};
    my $group_ref = $project_conf_ref->{GROUP};
    my $group_comp_ref = $project_conf_ref->{GRP_COMPARISON};

    foreach my $lane_count (sort {$a <=> $b} keys %{$lanes_ref}){
      if (defined($lanes_ref->{$lane_count}->{library_id}) && defined($lanes_ref->{$lane_count}->{flowcell_name}) && defined($lanes_ref->{$lane_count}->{lane_number}) && defined($lanes_ref->{$lane_count}->{data_path}) && defined($lanes_ref->{$lane_count}->{seq_file_type}) && defined($lanes_ref->{$lane_count}->{read_length}) && defined($lanes_ref->{$lane_count}->{read_trim}) && defined($lanes_ref->{$lane_count}->{max_n}) && defined($lanes_ref->{$lane_count}->{min_phred}) && defined($lanes_ref->{$lane_count}->{qual_type})){

      }else{
        $warnings++;
        $wstring .= "\nFormat of LANE entry is not correct for LANE entry $lane_count - Check Project config file LANE entry $lane_count";
      }
    }
    #Check the path for each lane of data
    foreach my $lane_count (sort {$a <=> $b} keys %{$lanes_ref}){
      if (defined($lanes_ref->{$lane_count}->{data_path})){
        my $source_dir = $lanes_ref->{$lane_count}->{data_path};
        unless (-e $source_dir && -d $source_dir){
          $warnings++;
          $wstring .= "\nCould not verify lane data path: $source_dir - Check Project config file LANE entry";
        }
      }
    }

    foreach my $lib_count (sort {$lib_ref->{$a}->{library_count} <=> $lib_ref->{$b}->{library_count}} keys %{$lib_ref}){
       if (defined($lib_ref->{$lib_count}->{library_count}) && defined($lib_ref->{$lib_count}->{library_name}) && defined($lib_ref->{$lib_count}->{group})){

      }else{
        $warnings++;
        $wstring .= "\nFormat of LIBRARY entry is not correct for LIBRARY entry $lib_count - Check Project config file LIBRARY entry $lib_count";
      }
    }

    foreach my $comp_count (sort {$a <=> $b} keys %{$comp_ref}){
       if (defined($comp_ref->{$comp_count}->{libraryA_id}) && defined($comp_ref->{$comp_count}->{libraryA_name}) && defined($comp_ref->{$comp_count}->{libraryB_id}) && defined($comp_ref->{$comp_count}->{libraryB_name}) && defined($comp_ref->{$comp_count}->{comparison_id})){

      }else{
        $warnings++;
        $wstring .= "\nFormat of COMPARISON entry is not correct for COMPARISON entry $comp_count - Check Project config file COMPARISON entry $comp_count";
      }
    }

    foreach my $group_count (sort {$a <=> $b} keys %{$group_ref}){
       if (defined($group_ref->{$group_count}->{group_name}) && defined($group_ref->{$group_count}->{library_list}) ){

      }else{
        $warnings++;
        $wstring .= "\nFormat of GROUP entry is not correct for GROUP entry $group_count - Check Project config file GROUP entry $group_count";
      }
    }

    foreach my $group_comp_count (sort {$a <=> $b} keys %{$group_comp_ref}){
       if (defined($group_comp_ref->{$group_comp_count}->{group_list}) && defined($group_comp_ref->{$group_comp_count}->{group_comparison_id}) ){

      }else{
        $warnings++;
        $wstring .= "\nFormat of GRP_COMPARISON entry is not correct for GRP_COMPARISON entry $group_comp_count - Check Project config file GRP_COMPARISON entry $group_comp_count";
      }
    }

  }

  #####################################################################################################################################
  #3.) Next check all the annotation_conf values
  if ($annotation_conf_ref){

    #ENSEMBL_VERSION
    if ($annotation_conf_ref->{ENSEMBL_VERSION}){
      unless ($annotation_conf_ref->{ENSEMBL_VERSION} =~/^\d+$/){
        $warnings++;
        $wstring .= "\nInteger format expected but not found.  Check Annotation config file parameter: ENSEMBL_VERSION";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Annotation config file parameter: ENSEMBL_VERSION";
    }

    #ALEXA_SEQ_DB
    if ($annotation_conf_ref->{ALEXA_SEQ_DB}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Annotation config file parameter: ALEXA_SEQ_DB";
    }

    #SPECIES_NAME
    if ($annotation_conf_ref->{SPECIES_NAME}){
      if ($annotation_conf_ref->{SPECIES_NAME} =~ /\s+/){
        $warnings++;
        $wstring .= "\nParameter should not have spaces: Check Annotation config file parameter: SPECIES_NAME";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Annotation config file parameter: SPECIES_NAME";
    }

    #SPECIES_NAME_COMMON
    if ($annotation_conf_ref->{SPECIES_NAME_COMMON}){
      if ($annotation_conf_ref->{SPECIES_NAME_COMMON} =~ /\s+/){
        $warnings++;
        $wstring .= "\nParameter should not have spaces: Check Annotation config file parameter: SPECIES_NAME_COMMON";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Annotation config file parameter: SPECIES_NAME_COMMON";
    }

    #ENTREZ_SPECIES_GROUP
    if ($annotation_conf_ref->{ENTREZ_SPECIES_GROUP}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Annotation config file parameter: UCSC_BUILD";
    }

    #UCSC_BUILD
    if ($annotation_conf_ref->{UCSC_BUILD}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Annotation config file parameter: UCSC_BUILD";
    }

    #JUNCTION_SEQ_LENGTH
    if ($annotation_conf_ref->{JUNCTION_SEQ_LENGTH}){
      unless ($annotation_conf_ref->{JUNCTION_SEQ_LENGTH} =~/^\d+$/){
        $warnings++;
        $wstring .= "\nInteger format expected but not found.  Check Annotation config file parameter: JUNCTION_SEQ_LENGTH";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Annotation config file parameter: JUNCTION_SEQ_LENGTH";
    }
    #BOUNDARY_SEQ_LENGTH
    if ($annotation_conf_ref->{BOUNDARY_SEQ_LENGTH}){
      unless ($annotation_conf_ref->{BOUNDARY_SEQ_LENGTH} =~/^\d+$/){
        $warnings++;
        $wstring .= "\nInteger format expected but not found.  Check Annotation config file parameter: BOUNDARY_SEQ_LENGTH";
      }
    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Annotation config file parameter: BOUNDARY_SEQ_LENGTH";
    }

    #ALIGNMENT OPTION
    if ($annotation_conf_ref->{ALIGNMENT_OPTION}){

    }else{
      $warnings++;
      $wstring .= "\nParameter not found: Add Annotation config file parameter: ALIGNMENT_OPTION";
    }

  }

  if ($warnings){
    print YELLOW, "\n\nFound $warnings warnings as follows:", RESET;
    print YELLOW, "\n$wstring\n\n", RESET;

    #Ask the user if they wish to proceed
    print YELLOW, "If you want the pipeline to run smoothly it is advisable to resolve these before proceeding ...", RESET;
    print YELLOW, "\nDo you still wish to proceed (y/n)? ", RESET;
    my $answer = <>;
    chomp($answer);
    unless($answer =~ /y|yes/i){
      print YELLOW, "\nAborting then...\n\n", RESET;
      exit();
    }
  }else{
    print YELLOW, "\n\tNo problems identified", RESET;
  }

  return($warnings);
}


#############################################################################################################################
#Add commas to number.  e.g. 1000000 to 1,000,000                                                                           #
#############################################################################################################################
sub commify {
   local $_  = shift;
   1 while s/^(-?\d+)(\d{3})/$1,$2/;
   return $_;
}


#############################################################################################################################
#Return message describing memory usage of the current process                                                              #
#############################################################################################################################
sub memoryUsage{
  my $pid = $$;
  my $ps_query = `ps -p $pid -o pmem,rss`;
  my @process_info = split ("\n", $ps_query);
  my $memory_usage = '';
  my $memory_usage_p = '';
  if ($process_info[1] =~ /(\S+)\s+(\S+)/){
    $memory_usage_p = $1;
    $memory_usage = $2;
  }
  my $memory_usage_m = sprintf("%.1f", ($memory_usage/1024));
  my $message = "Memory usage: $memory_usage_m Mb ($memory_usage_p%)";
  return($message);
}

1;




