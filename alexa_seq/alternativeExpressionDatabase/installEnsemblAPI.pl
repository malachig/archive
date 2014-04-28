#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is the help install the EnsEMBL API (user specifies the version to install)

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

#Initialize command line options
my $install_dir = '';
my $ensembl_version = '';

GetOptions ('install_dir=s'=>\$install_dir, 'ensembl_version=i'=>\$ensembl_version);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the directory that the Ensembl API will be installed into (a sub-directory for this version will be created automatically) using: --install_dir", RESET;
print GREEN, "\n\tSpecify the version of Ensembl that you wish to install using: --version", RESET;
print GREEN, "\n\nExample: installEnsemblAPI.pl  --install_dir=/home/malachig/svn/alexa_seq/ensembl_api/  --ensembl_version=53\n\n", RESET;

#Make sure all options were specified
unless ($install_dir && $ensembl_version){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Make sure CVS is available
unless (-e "/usr/bin/cvs"){
  print RED, "\nCould not find cvs.  Either it is not installed or not in your path.  You need cvs to install the Ensembl API\n\n", RESET;
  exit();
}

#See if the user has an entry in their .cvspass file for Ensembl
my $help_string = "\n\ncvs -d :pserver:cvsuser\@cvsro.sanger.ac.uk:/cvsroot/CVSmaster login\npassword: CVSUSER\n\n";
my $home_dir = $ENV{"HOME"};
my $cvspass = "$home_dir/".".cvspass";
if (-e $cvspass){
  my $ensembl_found = 0;
  open(CVSPASS, "$cvspass") || die "\n\nCould not open .cvspass file: $cvspass\n\n";
  while(<CVSPASS>){
    if ($_ =~ /cvsro\.sanger\.ac\.uk/){
      $ensembl_found = 1;
    }
  }
  close(CVSPASS);
  unless($ensembl_found){
    print RED, "\n\nCould not find a .cvspass file.  Execute the following at the command line and then rerun this script: $help_string", RESET;
    exit();
  }
}else{
  print RED, "\n\nCould not find a .cvspass file.  Execute the following at the command line and then rerun this script: $help_string", RESET;
  exit();
}

#Check install dir
$install_dir = &checkDir('-dir'=>$install_dir, '-clear'=>"no");

#Create installation sub-dir
my $sub_dir_name = "ensembl_"."$ensembl_version"."_perl_API/";

#If the sub_dir already contains all neccessary directories, perhaps this version is already installed
my $sub_dir = "$install_dir"."$sub_dir_name";
my $found_api = 0;
if (-e $sub_dir){
  #Sub dir already exists. Check for installed API
  my $dir1 = "$sub_dir"."ensembl/";
  my $dir2 = "$sub_dir"."ensembl-compara/";
  my $dir3 = "$sub_dir"."ensembl-external/";
  my $dir4 = "$sub_dir"."ensembl-variation/";

  if ((-e $dir1) && (-e $dir2) && (-e $dir3) && (-e $dir4)){
    print YELLOW, "\n\nIt appears that this EnsEMBL API is already installed\n\n", RESET;
    $found_api = 1;
  }else{
    $sub_dir = &createNewDir('-path'=>$install_dir, '-new_dir_name'=>$sub_dir_name);
  }
}else{
  $sub_dir = &createNewDir('-path'=>$install_dir, '-new_dir_name'=>$sub_dir_name);
}


#Now move to this directory
chdir($sub_dir);

#Now execute the CVS commands to install the API
unless ($found_api == 1){
  my $cmd1 = "cvs -d :pserver:cvsuser\@cvsro.sanger.ac.uk:/cvsroot/CVSmaster checkout -r branch-ensembl-$ensembl_version  ensembl";
  my $cmd2 = "cvs -d :pserver:cvsuser\@cvsro.sanger.ac.uk:/cvsroot/CVSmaster checkout -r branch-ensembl-$ensembl_version  ensembl-external";
  my $cmd3 = "cvs -d :pserver:cvsuser\@cvsro.sanger.ac.uk:/cvsroot/CVSmaster checkout -r branch-ensembl-$ensembl_version  ensembl-compara";
  my $cmd4 = "cvs -d :pserver:cvsuser\@cvsro.sanger.ac.uk:/cvsroot/CVSmaster checkout -r branch-ensembl-$ensembl_version  ensembl-variation";

  print BLUE, "\n\n$cmd1\n", RESET;
  system($cmd1);
  print BLUE, "\n\n$cmd2\n", RESET;
  system($cmd2);
  print BLUE, "\n\n$cmd3\n", RESET;
  system($cmd3);
  print BLUE, "\n\n$cmd4\n", RESET;
  system($cmd4);
}


#Now check for BioPerl and install it if neccessary
chdir($install_dir);

my $bio_perl_name = "bioperl-1.4";
my $bio_perl_dir = "$install_dir"."$bio_perl_name/";

if (-e $bio_perl_dir && -d $bio_perl_dir){
  print YELLOW, "\n\nIt appears that the BioPerl ($bio_perl_name) is already installed\n\n", RESET;
}else{

  #Download BioPerl version 1.4
  my $source_name = "bioperl-1.4";
  my $source_name_tar = "$source_name".".tar";
  my $source_name_tar_gz = "$source_name_tar".".gz";

  my $wget_cmd = "wget http://bioperl.org/DIST/bioperl-1.4.tar.gz";
  system($wget_cmd);

  #Decompress
  my $gz_cmd = "gunzip $source_name_tar_gz";
  print BLUE "\n\n$gz_cmd";
  system($gz_cmd);

  #Unpack 
  my $tar_cmd = "tar -xvf $source_name_tar";
  print BLUE "\n\n$tar_cmd";
  system($tar_cmd);

  #Configure
  chdir($bio_perl_dir);
  my $cfg_cmd = "perl Makefile.PL";
  print BLUE "\n\n$cfg_cmd  (answer the questions that come up - or simply press <enter> for each one to accept the default)";
  system($cfg_cmd);

  #Make
  my $make_cmd = "make";
  print BLUE "\n\n$make_cmd";
  system($make_cmd);

  #Delete the tar file
  my $rm_cmd = "rm -f $install_dir"."$source_name_tar";
  print BLUE "\n\n$rm_cmd";
  system($rm_cmd);


}

print "\n\n";

exit();

