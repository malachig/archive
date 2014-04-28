#!/usr/bin/perl

#Written by Ryan Morin

use strict;
use BerkeleyDB;

my %h;

my $chr = $ARGV[0] || die "$0 chrX\n";
my $filename = "bins/$chr.btree";
unlink $filename;



tie(%h, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $filename , -Flags => DB_CREATE) or die "can't open file $filename: $! $BerkeleyDB::Error\n";

my $n;
#open IN, "exonic_bases_gened_tabbed.txt.chr_nsort" or die "$!\n";
open IN, "exonic_bases_gened_exoned_tabbed.txt" or die "$!\n";
my @test;
while(<IN>){
    chomp;
    my ($a,$b,$c,$d) = split /\t/;
    next unless $a =~ /$chr$/;
  
    $h{$b} = "$c\t$d";
    $n++;
    unless($n%5000){
	push @test, $b;
    }
}
print "added $n lines\n";

for(@test){
    my $gene = $h{$_};
    print "$_\t$gene\n";
}

untie %h ;
