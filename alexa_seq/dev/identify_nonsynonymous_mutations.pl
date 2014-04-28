#!/usr/bin/perl

use lib "/home/malachig/perl/bioperl-1.4";

use strict;
use Data::Dumper;
use Bio::Seq;


my %lookup;

my %error_count; #by gene, for debugging purposes
my %grantham_distances;

load_ambig();
load_grantham();

while(<STDIN>){
    
    chomp;
    print;
    my @cols = split /\s/;
    my $gene = $cols[1];
    my $strand = $cols[2];
    my $codon = $cols[3];
    unless ($codon =~ /[ACTG]{3}/){
	print "\n";
	next;
    }
    my $codon_pos = $cols[4];
    my $ref = $cols[5];
    my $amb = $cols[6];
    my $mut = amb_lookup($amb,$ref,$cols[0]);
    next if $mut == -1;
    my $mut_cor = $mut;
    
    my $ref_cor = $ref;
    if($strand < 0){
	$mut_cor =~ tr/ACTG/TGAC/;
	$ref_cor =~ tr/ACTG/TGAC/;
    }
    my @three = split //, $codon;
    #sanity check
  
    my $codon_ind = $codon_pos -1;
    my $orig = $three[$codon_ind];
    if($orig ne $ref_cor){
	print STDERR "$cols[0]\tERROR, base $ref is not the same as codon position $codon_pos ($orig)\n";
	$error_count{$gene}++;
    }
    $three[$codon_ind] = $mut_cor;
    my $new_codon = join "", @three;
    print " $codon $new_codon ";
    my $old = Bio::Seq->new(-seq=>$codon);
    my $old_aa = $old->translate;
    my $new = Bio::Seq->new(-seq=>$new_codon);
    my $new_aa = $new->translate;
    my $old_aa_base = $old_aa->seq;
    my $new_aa_base = $new_aa->seq;
    my $type = "CODING";
    if($old_aa_base eq $new_aa_base){
	$type = "SYNONYMOUS";
    }
    my ($charge_effect,$polarity_effect) = rad_cons($old_aa_base,$new_aa_base);
    my $gdist = $grantham_distances{$old_aa_base}{$new_aa_base};
    print "$old_aa_base $new_aa_base $type $charge_effect $polarity_effect $gdist\n";
}
#print Dumper %error_count;

sub amb_lookup{
    my $amb = shift;
    my $ref = shift;
    my $pos = shift;
    if($amb =~ /[ACTG]/){
	#print "$amb\n";
	return($amb);
    }
    #print "LOOKUP: $amb\n";
    my @all = @{$lookup{$amb}};
    #print Dumper @all;
    my @some;
    for(@all){
	next if $_ eq $ref;
	push @some, $_;

    }
    if(@some > 1){
	#print Dumper @some;
	print STDERR "$pos\t$ref\t$amb\ttwo bases at this position, neither is reference\n";
	return(-1);
    }
    return($some[0]);
}
sub load_grantham{
    my $matrix_file = "/projects/rmorin/common/ucsc_files/grantham_matrix.txt";
    my @col_names = split //, "ARNDCQEGHILKMFPSTWYV";
    open M, $matrix_file or die "$!\n";
    my $i = 0;

    while(<M>){
	next if /^\#/;
	my @fields = split /\s+/;
	shift @fields;
	for (my $j = 0;$j < @fields;$j++){
	    $fields[$j] =~ s/\.//;
	    #print "d($col_names[$i],$col_names[$j]) = $fields[$j]\n";
	    $grantham_distances{$col_names[$i]}{$col_names[$j]} = $fields[$j];
	    $grantham_distances{$col_names[$j]}{$col_names[$i]} = $fields[$j];
	}
	$i++;
    }
}

sub rad_cons{
    my $before = shift;
    my $after = shift;
    my @pos = qw(R H K);
    my @neg = qw(D E);
    my @neutral = qw(A N C Q G I L M F P S T W Y V);
    my %charges;
    for(@pos){
	$charges{$_} = "POS";
    }
    for(@neg){
	$charges{$_} = "NEG";
    }
    for(@neutral){
	$charges{$_} = "NEUTRAL";
    }
    my %polarity;
    for(qw(C)){
	$polarity{$_} = "SPECIAL";
    }
    for(qw(A G P S T)){
	$polarity{$_} = "SMALL_NEUTRAL";
    }
    for(qw(N D Q E)){
	$polarity{$_} = "SMALL_POLAR";
    }
    for(qw(R H K)){
	$polarity{$_} = "LARGE_POLAR";
    }
    for(qw(I L M V)){
	$polarity{$_} = "SMALL_NONPOLAR";
    }
    for(qw(F W Y)){
	$polarity{$_} = "LARGE_NONPOLAR";
    }
    my $before_charge = $charges{$before};
    my $after_charge = $charges{$after};
    my $charge_change = "CONSERVATIVE";
    if($before_charge ne $after_charge){
	$charge_change = "RADICAL";
    }
    my $polar_change = "CONSERVATIVE";
    if($polarity{$before} ne $polarity{$after}){
	$polar_change = "RADICAL";
    }
    return($charge_change,$polar_change);
}
sub load_ambig{
    add('M',qw(A C));
    add('R',qw(A G));
    add('W',qw(A T));
    add('S',qw(C G));
    add('Y',qw(T C));
    add('K',qw(T G));
    add('V',qw(A C G));
    add('H',qw(A C T));
    add('D',qw(A G T));
    add('B',qw(C G T));
}

sub add{

    my $amb = shift;
    my @copy = @_;
    $lookup{$amb} = \@copy;
}
