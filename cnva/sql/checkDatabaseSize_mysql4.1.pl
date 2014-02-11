#!/usr/bin/perl

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use DBI;

my $database = '';
my $server = '';
my $user = '';
my $password = '';

GetOptions ('database=s'=>\$database, 'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password);

#Usage instructions 
print GREEN, "\nThis script calculates the disk usage (size) of a specified database on a particular host", RESET;
print GREEN, "\nUse --server to specify the database host and --database to specify the database name", RESET;
print GREEN, "\nUse --user to specify the username and --password to specify your password", RESET;
print GREEN, "\nUsage: checkDatabaseSize_mysql4.1.pl  --database=ALEXA_hs_41_36c  --server=server_name  --user=user_name  --password=pwd\n\n", RESET;

unless ($database && $server && $user && $password){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Establish connection with the Alternative Splicing Expression database
my $dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

my $total_size = 0;

my $sql = "show table status";

my $sth = $dbh->prepare("$sql");
$sth->execute();

print BLUE, "\nSummarize size of tables for $database", RESET;

while (my (@row) = $sth->fetchrow_array()){

  my $table_name = $row[0];
  my $data_length = $row[6];
  my $index_length = $row[8];

  my $data_size = (($data_length/1024)/1024);
  my $index_size = (($index_length/1024)/1024);

  my $data_size_rounded = sprintf("%.3f", $data_size);
  my $index_size_rounded = sprintf("%.3f", $index_size);

  print BLUE, "\nTable_name: $table_name\tData_Size: $data_size_rounded\tIndex_Size: $index_size_rounded", RESET;

  $total_size += $data_size;
  $total_size += $index_size;
}
$sth->finish();

my $total_size_rounded = sprintf("%.3f", $total_size);
print YELLOW, "\n\nTotal Table Size for $database = $total_size_rounded Mb\n\n", RESET;

#Close database connection
$dbh->disconnect();

exit();

##########################################################################################################################
#mysql db connection                                                                                                     #
##########################################################################################################################

sub connectDB {
  my %args = @_;
  my $database_name = $args{'-database'};
  my $database_host = $args{'-server'};
  my $user_name = $args{'-user'};
  my $user_pw = $args{'-password'};

  my $dbh = DBI->connect( "dbi:mysql:database=$database_name;host=$database_host", $user_name, $user_pw, { PrintError => 1 } );
  return $dbh;
}
