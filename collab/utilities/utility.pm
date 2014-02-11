=head1 NAME

utility.pm - library modules that contains generic utilities

=head1 SYNOPSIS

use utility qw(:all);

=head2 NOTE

currently located in '~/Array_design/utilities'

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

package utilities::utility;
require Exporter;

@ISA = qw( Exporter );
@EXPORT = qw();

@EXPORT_OK = qw(&connectDB &createNewDir);

%EXPORT_TAGS = (
     all => [qw(&connectDB &createNewDir)]
);

use strict;
use Data::Dumper;
use DBI;
use Term::ANSIColor qw(:constants);

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
  return $dbh;
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


1;




