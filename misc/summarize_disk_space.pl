#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

#Griffith lab DiskGroup summary:
#https://imp-lims.gsc.wustl.edu/entity/disk-group/39999

my @dirs = qw (/gscmnt/gc2142/griffithlab/ 
               /gscmnt/sata206/griffithlab/ 
               /gscmnt/gc2547/griffithlab/ 
               /gscmnt/gc2502/griffithlab/ 
               /gscmnt/gc2602/griffithlab/ 
               /gscmnt/gc2607/griffithlab/ 
               /gscmnt/gc2543/griffithlab/
               /gscmnt/gc2736/griffithlab_gms/
               );

my $grand_total_allocated_tb = 0;
my $grand_total_used_tb = 0;
my $grand_total_available_tb = 0;

my $c = 0;
foreach my $dir (@dirs){
  $c++;
  my $result1 = `ls $dir`;
  my $result2 = `df -h $dir | grep gscmnt`;
  chomp $result2;
  if ($result2 =~ /\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+/){
    my $allocated = $1;
    my $used = $2;
    my $available = $3;

    #Format for KB, MB, GB, TB -> convert to TB
    my $allocated_tb = convert_size($allocated);
    my $used_tb = convert_size($used);
    my $available_tb = convert_size($available);

    print "\n$c.) Allocated = $allocated_tb Tb; Used = $used_tb Tb; Available = $available_tb Tb [$dir]";
    $grand_total_allocated_tb += $allocated_tb;
    $grand_total_used_tb += $used_tb;
    $grand_total_available_tb += $available_tb;

  }
}

$grand_total_allocated_tb = sprintf ("%.1f", $grand_total_allocated_tb);
$grand_total_available_tb = sprintf ("%.1f", $grand_total_available_tb);

print "\n\nGrand Total Allocated = $grand_total_allocated_tb Tb; Grand Total Used = $grand_total_used_tb Tb; Grand Total Available = $grand_total_available_tb Tb\n\n";

exit;


sub convert_size{
  my $size_string = shift;
  my $size_tb = $size_string; 


  if ($size_string =~ /^(\d+)(\w)$/){
    $size_string = $1.".0".$2;
  }
  if ($size_string =~ /^(\d+)$/){
    $size_string = $1.".0";
  }

  if ($size_string =~ /(\d+\.\d+)T/){
    $size_tb = $1
  }elsif ($size_string =~ /(\d+\.\d+)G/){
    $size_tb = $1/1000;
  }elsif ($size_string =~ /(\d+\.\d+)M/){
    $size_tb = $1/1000000;
  }elsif ($size_string =~ /(\d+\.\d+)K/){
    $size_tb = $1/1000000000;
  }elsif ($size_string =~ /(\d+\.\d+)/){
    $size_tb = $1/1000000000000;
  }else{
    die "\n\nFormat of size string not recognized: $size_string\n\n";
  }

  return $size_tb;
}



