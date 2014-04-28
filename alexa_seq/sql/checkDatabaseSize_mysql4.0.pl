#!/usr/bin/perl
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.


use strict;
use Data::Dumper;
use Getopt::Std;

#Always use the seqdev folder as the library root and then specify the packages to use.
use lib '/usr/local/ulib/beta/seqdev';
use utilities::utility qw(:all);

getopts("h:d:u:p:");
use vars qw($opt_h $opt_d $opt_u $opt_p);

#Usage instructions 
unless ($opt_h && $opt_d && $opt_u && $opt_p){
  print "\nThis calculates the disk usage (size) of a specified database on a particular host";
  print "\nUse -h to specify the database host and -d to specify the database name";
  print "\nUse -u to specify the username and -p to specify your password";
  print "\nUsage: checkDatabaseSize.pl -h seqdev02 -d ALEXA -u viewer -p viewer\n\n";
  exit();
}

#Establish connection with the Alternative Splicing Expression database
my $dbh = connectDB($opt_d, $opt_h, $opt_u, $opt_p);

my $total_size = 0;

my $sql = "show table status";

my $sth = $dbh->prepare("$sql");
$sth->execute();


while (my (@row) = $sth->fetchrow_array()){

  my $table_name = $row[0];
  my $data_length = $row[5];
  my $index_length = $row[7];

  my $data_size = (($data_length/1024)/1024);
  my $index_size = (($index_length/1024)/1024);

  my $data_size_rounded = sprintf("%.3f", $data_size);
  my $index_size_rounded = sprintf("%.3f", $index_size);

  print "\nTable_name: $table_name\tData_Size: $data_size_rounded\tIndex_Size: $index_size_rounded";

  $total_size += $data_size;
  $total_size += $index_size;
}
$sth->finish();

my $total_size_rounded = sprintf("%.3f", $total_size);
print "\n\nTotal Table Size = $total_size_rounded Mb\n\n";

#Close database connection
$dbh->disconnect();

exit();
