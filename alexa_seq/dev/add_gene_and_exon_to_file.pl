#!/usr/bin/perl

#Written by Ryan Morin


use strict;
use lib '/home/rmorin/lib/lib64/perl5/site_perl/5.8.5/x86_64-linux-thread-multi/';
use BerkeleyDB;
my $read_len; #only needed if format is mapview
my $file_format; #guess format using first liine
my %all;

my $add_strand = $ARGV[0];

my $base = "/projects/rmorin/common/exon_content/ensembl_human_49/bins/";
chomp(my @bins = `ls $base`);

for my $bin(@bins){
    $bin =~ /(\S+)\.btree/;
    my $chr = $1;
    my $filename = $base . $bin;
    unless(-e $filename){
	warn("hash file $filename does not exist\n");
    }
    print STDERR "tieing $chr hash\n";
    tie(my %h, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $filename, -Flags => 'DB_RDONLY' ) or die "can't open file $filename: $! $BerkeleyDB::Error\n";
    $all{$chr} = \%h;
}
my $chr_index;
my $pos_index;

my $gene_strand = "/projects/rmorin/common/exon_content/ensembl_human_49/gene_strand.txt";
open G, $gene_strand or die "$!\n";
my %gene_strand;
while(<G>){
    chomp;
    my ($gene,$strand) = split;
    $gene_strand{$gene} = $strand;
}


my $first_line = 1;
while(<STDIN>){
    chomp;
    next if /^\#/;
    if($first_line){
	$first_line = undef;
	#check format by counting columns
	my @a = split /\s/;
	my $an = @a;
	if(@a == 1){
	    if($a[0] =~ /(\S+):(\d+)/){
		$file_format = "join";
	    }
	}
	elsif($a[0] =~ /(\S+):(\d+)/){
	    $file_format = "join";
	}
	elsif($a[0] =~ /chr/ && $a[1] =~ /^\d+$/ && $a[2] =~ /^\d+$/){
	    $file_format = "bed";
	}
	elsif(@a == 16 || @a ==14){
	    $file_format = "mapview";
	}
	elsif(($an == 7 || $an == 8 || $an == 4 || $an == 5)){
	    $file_format = "pileup";
	}
	elsif($a[0] =~ /:/){
	    $file_format = "join";
	}
	elsif(@a == 9 ){
	    $file_format = "snp";
	}
	else{
	    die "unrecognized format\n$_\t($an)\n";
	}

	print STDERR "recognized file format is $file_format\t($an fields)\n";
    }
    my $line = $_;
    
    my ($chr,$s,$gene,$exon);
    my @vals = split /\s/, $line;
    if($file_format eq 'pileup'){
	$chr = $vals[0];
	$s = $vals[1];
	unless($chr =~ /chr/){
	    $chr = "chr" . $chr;
	}
	$gene = $all{$chr}{$s};
#	print STDERR "$chr $s $gene\n";
    }
    elsif($file_format eq "bed"){
	$chr = $vals[0];
	$s = $vals[1];
	my $e = $vals[2];
	my $gene_ex;
	for($s..$e){
	    $gene_ex = $all{$chr}{$_};
	    last if $gene_ex;
	}
	my ($g,$e) = split /\t/,$gene_ex;

	my $other = join "\t", @vals[4..$#vals];
	if($add_strand){
	    print "$vals[0]\t$vals[1]\t$vals[2]\t$vals[3]_$g\t$other\t$gene_strand{$g}\n";
	}
	else{
	    print "$vals[0]\t$vals[1]\t$vals[2]\t$vals[3]_$g\t$other\n";
	    next;
	}
    }
    elsif($file_format eq 'mapview'){
	$chr = $vals[1];
	unless($chr =~ /chr/){
            $chr = "chr" . $chr;
	}
	$s = $vals[2];
	my $e = $s + $read_len;
	for($s..$e){
	    $gene = $all{$chr}{$s};
	    last if $gene;  #just in case the first few nt of a read aren't mapped toa  gene
	}
	
    }
    elsif($file_format eq 'snp'){
	$chr = $vals[0];
        $s = $vals[1];
        $gene = $all{$chr}{$s};
    }
    elsif($file_format eq 'join'){
	($chr,$s) = split /:/, $vals[0];
	unless($chr =~ /chr/){
	    $chr = "chr" . $chr;
	}
	$gene = $all{$chr}{$s};
    }
    if($add_strand){
	if($gene =~ /^(\S+)/){
	    $gene = $1;
	    $exon = $2;
	}
	print "$line\t$gene\t$exon\t$gene_strand{$gene}\n";
    }else{
	print "$line\t$gene\n";
    }
}
