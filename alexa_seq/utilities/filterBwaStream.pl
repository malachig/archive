#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Filter 'BWA samse' output on the fly via a pipe of standard out from a BWA->SAMSE job
#Filter according to bit_score:

#Example BWA Alignment command (slower more extensive result):
#Allow up to three mismatches in aligning (-n 3), search for all alignments with up to 3 mismatches (-N), show up to 100 total alignments (-n 100) 
#/home/malachig/tools/bwa/bwa-0.5.8a/bwa aln -n 3 -N db/transcripts.fa.gz reads/test2.fa 2>/dev/null | /home/malachig/tools/bwa/bwa-0.5.8a/bwa samse -n 100 db/transcripts.fa.gz - reads/test2.fa  2>/dev/null

#Example BWA Alignment command (faster result):
#Let BWA decide the max mismatches based on read length, display up 100 total alignments (-n 100)
#/home/malachig/tools/bwa/bwa-0.5.8a/bwa aln db/transcripts.fa.gz reads/fasta_0000 2>/dev/null | /home/malachig/tools/bwa/bwa-0.5.8a/bwa samse -n 100 db/transcripts.fa.gz - reads/fasta_0000  2>/dev/null | /home/malachig/svn/alexa_seq/utilities/filterBwaStream.pl --min_bit_score=40.0

#USAGE:
#BWA COMMAND or cat of BWA output file | filterBwaStream.pl  --min_bit_score=40.0 > filtered_bwa_file.txt

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;

my $min_bit_score = '';

GetOptions ('min_bit_score=f'=>\$min_bit_score);

unless ($min_bit_score =~ /^\d+/){
  print RED, "\nRequired input parameter(s) missing or incorrect format!\n\n", RESET;
  print RED, "\nExample usage: BWA (aln+samse) COMMAND or cat of BWA file | filterBwaStream.pl  --min_bit_score=40.0 > filtered_bwa_file.txt\n", RESET;
  exit();
}

#Sanity check of cutoff format
chomp($min_bit_score);
unless ($min_bit_score =~ m/^\d+.\d+$/ || $min_bit_score =~ m/^\d+$/){
  print "\nfilterSoapStream.pl was supplied a --min_bit_score value that does not appear to be a valid number!\n\n";
  exit();
}

#Define flags to check 
my $flag_read_unmapped = 0x0004;
my $flag_reverse_strand = 0x0010;


#Convert, the standard SOAP output format to BLAST tabular output format
#Using the alignment length and number of mismatches, deletions and insertions, calculate a bit score
while(<STDIN>){
  chomp($_);
  my @line = split("\t", $_);

  #Check for basic line format
  unless (scalar(@line) >= 11){
    #print RED, "\nIncomplete line: $_", RESET;
    next();
  }

  #Check the flag immediately and skip reads with no alignment
  my $flag = $line[1];
  my $aligned = ( ($flag & $flag_read_unmapped) ? 1 : 0 );
  if($aligned){
    next();
  }
  my $query_id = $line[0];
  my $subject_id = $line[2];
  my $subject_start = $line[3];
  my $cigar = $line[5];
  my $seq = $line[9];
  my $read_length = length($seq);

  #Strand
  my $strand = "+";
  my $reverse_strand = ( ($flag & $flag_reverse_strand) ? 1 : 0 );
  if($reverse_strand){
    $strand = "-";
  }

  #Grab all the optional fields as one big string
  my $opt_fields;
  if ($_ =~ /(\w+\:\w+\:\w+.*)$/){
    $opt_fields = $1;
  }else{
    print RED, "\n\nCould not identify optional fields in line:\n$_", RESET;
    exit();
  }

  #Mismatches (This is the total mismatching bases - includes substitutions + deletions + insertions!)
  my $mismatches = 0;
  if ($opt_fields =~ /NM\:\w\:(\d+)/){
    $mismatches = $1;
  }else{
    print RED, "\n\nCould not find number of mismatches in line:\n$_", RESET;
    exit();
  }

  #Store all the alignments basic info in a single array
  my @hits;
  my $first_record = "$subject_id,$strand"."$subject_start,$cigar,$mismatches";
  push(@hits, $first_record);

  #Grab the alternate alignments
  if ($opt_fields =~ /XA\:\w\:(.*)/){
    my $alt_align_string = $1;
    my @alt_hits = split(";", $alt_align_string);
    push(@hits, @alt_hits);
  }

  #Now go through each set of basic alignment info and build a BLAST record
  foreach my $align_string (@hits){

    #Basic info for each alignment
    #$subject_id, $strand, $subject_start, $cigar, $mismatches
    #e.g. "24875,+1695,51M,0"  e.g. "17471,-532,51M,0"
    my @vals = split(",", $align_string);
    my $subject_id = $vals[0];
    my $strand_pos = $vals[1];
    my $cigar = $vals[2];
    my $mismatches = $vals[3];
    my $strand;
    my $subject_start;
    if ($strand_pos =~ /([\-|\+])(\d+)/){
      $strand = $1;
      $subject_start = $2;
    }else{
      print RED, "\n\nAlignment format not understood:\n$align_string\n\n", RESET;
      exit();
    }

    #print YELLOW, "\n$subject_id $strand $subject_start $cigar $mismatches", RESET;

    #Get the alignment length, query end position, and number of matching bases from the CIGAR string
    #Get the following from the CIGAR string: query_start, query_end, subject_start, subject_end, matching_bases, deletion_count, deletion_bases, insertion_count, insertion_bases
    my $cig_info = &parseCigar('-cigar'=>$cigar, '-subject_start'=>$subject_start, '-read_length'=>$read_length);
    my $query_start = $cig_info->{query_start};
    my $query_end = $cig_info->{query_end};
    $subject_start = $cig_info->{subject_start};
    my $subject_end = $cig_info->{subject_end};
    my $matching_bases = $cig_info->{matching_bases};
    my $deletion_count = $cig_info->{deletion_count};
    my $deletion_bases = $cig_info->{deletion_bases};
    my @deletion_sizes = @{$cig_info->{deletion_sizes}};
    my $insertion_count = $cig_info->{insertion_count};
    my $insertion_bases = $cig_info->{insertion_bases};
    my @insertion_sizes = @{$cig_info->{insertion_sizes}};
    my $gap_openings = $deletion_count + $insertion_count;

    #Alignment length
    my $align_length_q = ($query_end-$query_start)+1;
    my $align_length_s = ($subject_end-$subject_start)+1;
    my $align_length = $align_length_q;
    if ($align_length_s > $align_length_q){
      $align_length = $align_length_s;
    }

    #Percent ID
    my $percent_id = sprintf("%.2f", ((($align_length-$mismatches)/$align_length)*100));

    #Switch subject start and end for reverse strand alignments
    if ($strand eq "-"){
      my $temp = $subject_start;
      $subject_start = $subject_end;
      $subject_end = $temp;
    }

    #Adjusted number of mismatches (total mismatches - number of insertions and deletions)
    #my $adjusted_mismatches = $mismatches - ($deletion_bases+$insertion_bases);
    my $adjusted_mismatches = $mismatches - ($deletion_count+$insertion_count);
    my $indel_count = $deletion_count+$insertion_count;
    if ($adjusted_mismatches < 0){
      $adjusted_mismatches = 0;
    }

    #Calculate a BLAST-like bit score using the alignment length and number of mismatches, deletions and insertions
    #Bit.Score calculated as follows :
    #Base score = (Matching Bases*2)
    #Mismatch penalty = Substitution_Count * 8
    #Deletion penalty = +13.9 for the first base of each deletion and +4.0 for each additional base of each deletion (i.e. a deletion of 2 is (13.9+4.0) and deletion of 3 is (13.9+4.0+4.0)
    #Insertion penalty = +15.5 for the first base of each insertion and +6.0 for each additional base of each insertion
    #- (Adjusted_Mismatch_Count * 8) - (Deletion_Size*16) - (InsertionSize*4)
    my $del_penalty = 0;
    my $ins_penalty = 0;
    my $adjust = 0.3;
    foreach my $del_size (@deletion_sizes){$del_penalty += (13.9+(($del_size-1)*4.0));}
    foreach my $ins_size (@insertion_sizes){$ins_penalty += (15.9+(($ins_size-1)*6.0));}
    my $starting_score = ($align_length_q*2);
    my $mismatch_penalty = ($adjusted_mismatches*8);
    my $bit_score = $starting_score - $mismatch_penalty - $del_penalty - $ins_penalty - $adjust;
    my $calculation = "$bit_score = $starting_score - $mismatch_penalty - $del_penalty - $ins_penalty - $adjust";

    #Create an output line the is like the BLAST tabular output format
    #Tabular blast output format.  Fields:
    #        (0)Query id,(1)Subject id,(2)% identity,(3) alignment length,(4) mismatches,(5) gap openings,
    #        (6)q. start, (7)q. end, (8)s. start, (9)s. end, (10)e-value, (11)bit score

    my $blast_record = "$query_id\t$subject_id\t$percent_id\t$align_length\t$adjusted_mismatches\t$gap_openings\t$query_start\t$query_end\t$subject_start\t$subject_end\t1\t$bit_score";

    #If the bit_score exceeds the specified min, print out the result line
    if ($bit_score >= $min_bit_score){
      print "\n$blast_record";
    }
  }
}
exit();


##############################################################################################################################################################
#Based on a CIGAR, Subject Start and Read Length - Get various information from a CIGAR string                                                               #
##############################################################################################################################################################
sub parseCigar{
  my %args = @_;
  my $cigar = $args{'-cigar'};
  my $pos = $args{'-subject_start'};
  my $rlen = $args{'-read_length'};

  $cigar =~ s/(\d+)([MIDSN])/$1 $2\t/g;
  my $t_start = $pos;
  my $t_end = $pos-1;
  my $q_start = 1;
  my $q_end = 0;

  my $M_bases = 0;
  my $I_count = 0;
  my $I_bases = 0;
  my @I_bases;
  my $D_count = 0;
  my $D_bases = 0;
  my @D_bases;

  my $r_used = 0;

  my $DEBUG = 0;
  for my $c (split(/\t/,   $cigar)) {
    my ($len,  $op) = split(/ /,   $c);
    if($DEBUG && $r_used >= $rlen) { print STDERR "ERROR: more bases accounted for than read length! (rlen = $rlen, r_used = $r_used)\n"; }

    if($op eq 'M') { $t_end += $len; $q_end += $len; $r_used += $len; $M_bases += $len; }
    elsif($op eq 'D') { $t_end += $len; $D_bases += $len; $D_count++; push(@D_bases, $len);}
    elsif($op eq 'N') { $t_end += $len; }
    elsif($op eq 'I') { $q_end += $len; $r_used += $len; $I_bases += $len; $I_count++; push(@I_bases, $len);}
    elsif($op eq 'S') {
      if($q_start == 1 && $q_end == 0) { $q_start = $len+1; $r_used += $len;}
      else { $r_used += $len;}
    }
    else { print STDERR "ERROR: Unknown cigar op: $op\n"; }
  }
  if($r_used != $rlen) { print STDERR "ERROR: more bases accounted for than read length! (rlen = $rlen, r_used = $r_used)\n"; }

  #Return the query_start, query_end, subject_start, subject_end, matching_bases, deletion_count, deletion_bases, insertion_count, insertion_bases
  my %cig_info;
  $cig_info{query_start}=$q_start;
  $cig_info{query_end}=$q_end;
  $cig_info{subject_start}=$t_start;
  $cig_info{subject_end}=$t_end;
  $cig_info{matching_bases}=$M_bases;
  $cig_info{deletion_count}=$D_count;
  $cig_info{deletion_bases}=$D_bases;
  $cig_info{deletion_sizes}=\@D_bases;
  $cig_info{insertion_count}=$I_count;
  $cig_info{insertion_bases}=$I_bases;
  $cig_info{insertion_sizes}=\@I_bases;

  return(\%cig_info);
}


