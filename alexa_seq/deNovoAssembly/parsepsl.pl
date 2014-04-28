#!/usr/bin/perl -w
#Written by Kim Wong (kwong@bcgsc.ca) MARCH 2008
#Adopted by Malachi Griffith JUNE 2008
#Copyright 2009 Malachi Griffith and Kim Wong
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.


# Parse PSL files from BLAT:
# Input contig (query) list, blat psl format file.
# Count exact, partial hits, average coverage of targets

use strict;

my ($list,$pslfile,$histo) = @ARGV;

my @contigs = `cat $list`;
chomp @contigs;

# blat headers for psl format


#match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
#     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count

# note: output is sorted by strand then score..

#my $nohit=0;
#my ($exact_match,$internal,$partial1,$partial2,$partial3,
#$partial_gaps1,$partial_gaps2,$partial_gaps3) = (0,0,0,0,0,0,0,0);
my %hits = ();

# hash of hashes that hold the query name
# $matches-$mismatches
# $gaps[y/n]
# best hit output line

my %targets;
my %seen;
#my $count = 0;
open (PSL, "<$pslfile") or die;
while (<PSL>) {
   next if /^#/;
   next if /^psL|^\s*match|^\s*$|----/;
   my $hit = $_;
   my @info = split "\t", $hit;
   # determine matches-mismatches and #gaps
   my $score = $info[0];
   my $qgaps = $info[5];
   my $tgaps = $info[7];
   my $query = $info[9];
   my $tname = $info[13];
   my $tsize = $info[14];
   my $tstart = $info[15];
   my $tend = $info[16];

   # query not seen yet, or query exist but current score is better
   if (!$seen{$query}{score} || $seen{$query}{score} && $score > $seen{$query}{score}) {
      $seen{$query}{score} = $score;
      $seen{$query}{qgaps} = $qgaps;
      $seen{$query}{tgaps} = $tgaps;
      $seen{$query}{hit} = $hit;
      $seen{$query}{target} = $tname;
      $seen{$query}{tcoords} = [$tstart,$tend];
      #print STDERR "target $seen{$query}{target}, $seen{$query}{tcoords}[0], $seen{$query}{tcoords}[1]\n";
     
      # track targets that are hit
      $targets{$tname} = $tsize if !$targets{$tname};
   }
   # if score is the same
   elsif ($score ==  $seen{$query}{score}) {
      if ($seen{$query}{qgaps}+$seen{$query}{tgaps} > $qgaps+$tgaps) {
         $seen{$query}{score} = $score;
         $seen{$query}{qgaps} = $qgaps;
         $seen{$query}{tgaps} = $tgaps;
         $seen{$query}{hit} = $hit;
         $seen{$query}{target} = $tname;
         $seen{$query}{tcoords} = [$tstart,$tend];
         # track targets that are hit
         $targets{$tname} = $tsize if !$targets{$tname};
      }

   }
   else {
    #print STDERR "Passing: $hit ($score)\n$seen{$query}{score}\n";
   }
   #next if $seen{$info[9]};
   #push @blat, $_;
   #$seen{$info[9]} = $_;
   #$count++;
}
close PSL;

my $count = keys %seen;

my $contig_num = @contigs;


# tally number of mismatches in partial alignments with no internal mismatches/gaps
my %diffs;
my @diffs;
my %targethits;
my $totbases;
my $tottgap;
foreach my $contig (@contigs) {
      my $hit = $seen{$contig}{hit};
      if (!$hit) {
        #print "$contig\n";  # print out no hits
	#$nohit++;
	$hits{nohit}[0]++;
	$hits{nohit}[1] .= "$contig\n";
	print "$contig\tNOHIT\n";
	next;
      }
      push @{ $targethits{$seen{$contig}{target}} },  $seen{$contig}{tcoords} ;
      #print $hit;
      # look at best hit only (???)
      #for (@hits) {
           chomp $hit;
           my @info = split "\t", $hit;
           my ($match,$mismatch,$qgap,$tgap,$contig,$qsize,$qstart,$qend) = ($info[0],$info[1],$info[5],$info[7],$info[9],$info[10],$info[11],$info[12]);
	   #$qstart++; # add 1 to start  - it's a blat thing
	   my $chunk = $qend - $qstart;
	   $totbases += $match;
	   $tottgap += $tgap;
	   #print  "$hit\n" if $match == $qsize;
	   # exact match
	   if ($match == $qsize && $qgap==0 && $tgap==0) {
             #$exact_match++;
             #$exact_match_size += $qsize;
	     $hits{exact}[0]++;
	     $hits{exact}[1]+=$qsize;
	     print "$contig\tEXACT\t$mismatch\t$qgap\t$tgap\t$qsize\n";
            # print  "contig $contig, $info[13]\n";

	   }
	   # Not exact match but aligned from qstart to qend
	   # ie: internal mismatches or gaps 
	   elsif ($qstart==0 && $qend==$qsize) {
	      #$internal++;
	      #$internal_size;
	      $hits{internal}[0]++;
	      $hits{internal}[1]+=$qsize;
               print  "$contig\tINT\t$mismatch\t$qgap\t$tgap\t$qsize\n";
	       if ($mismatch > 0) {
	          $hits{fmismatch}[0]++;
	          $hits{fmismatch}[1] += $mismatch;
	       }
	       if ($qgap > 0) {
	          $hits{fqgap}[0]++;
	          $hits{fqgap}[1] += $qgap;
	       }
	       if ($tgap > 0) {
	          $hits{ftgap}[0]++;
	          $hits{ftgap}[1] += $tgap;

	       }

	   }
	   # Alignment does not include the whole query
	   # and no mismatches or qgapbases or tgapbases
	   elsif ($match < $qsize && $info[1]==0 && $info[5]==0 && $info[7]==0) {
               # classify partial alignments
	       #if ($qsize - $match < 50) {
	       if ($qsize - $chunk < 50) {
              #    $partial1++; 
	           $hits{partial1}[0]++;
	           $hits{partial1}[1]+=$qsize;
                   print  "$contig\tINT1\t$mismatch\t$qgap\t$tgap\t$qsize\t$hit\n"
	       }
	       elsif ($qsize - $chunk < 100) {
                   #$partial2++;
	           $hits{partial2}[0]++;
	           $hits{partial2}[1]+=$qsize;
                   print  "$contig\tINT2\t$mismatch\t$qgap\t$tgap\t$qsize\n"
	       }
	       else {  
	           #$partial3++;
	           $hits{partial3}[0]++;
	           $hits{partial3}[1]+=$qsize;
		   # print STDERR "$hit\n";
                   print  "$contig\tINT3\t$mismatch\t$qgap\t$tgap\t$qsize\n"
	       
	       }
	       my $diff = $qsize - $match;
	       $diffs{$diff}++;
	       push @diffs, $diff;
	   }
	   # Alignment does not include the whole query
	   # and there are mismatches/gaps
	   
           else {
               if ($qsize - $chunk < 50) {
                   #$partial_gaps1++; 
	           $hits{partial_gaps1}[0]++;
	           $hits{partial_gaps1}[1]+=$qsize;
                   print  "$contig\tP1\t$mismatch\t$qgap\t$tgap\t$qsize\n"
	       }
	       elsif ($qsize - $chunk < 100) {
                   #$partial_gaps2++;
	           $hits{partial_gaps2}[0]++;
	           $hits{partial_gaps2}[1]+=$qsize;
                   print  "$contig\tP2\t$mismatch\t$qgap\t$tgap\t$qsize\n"
	       }
	       else { 
	          #$partial_gaps3++ ;
	           $hits{partial_gaps3}[0]++;
	           $hits{partial_gaps3}[1]+=$qsize;
                   print  "$contig\tP3\t$mismatch\t$qgap\t$tgap\t$qsize\n"
	       
	       
	       }
	       if ($mismatch > 0) {
	          $hits{mismatch}[0]++;
	          $hits{mismatch}[1] += $mismatch;
	       }
	       if ($qgap > 0) {
	          $hits{qgap}[0]++;
	          $hits{qgap}[1] += $qgap;
	       }
	       if ($tgap > 0) {
	          $hits{tgap}[0]++;
	          $hits{tgap}[1] += $tgap;

	       }
	   }

      #}

}
for ('nohit','exact','internal','partial1','partial2','partial3','partial_gaps1','partial_gaps2','partial_gaps3','mismatch','qgap','tgap','fmismatch','fqgap','ftgap' ) {
    $hits{$_}[0] = '0' if !$hits{$_}[0];   # tally of contigs 
    $hits{$_}[1] = '0' if !$hits{$_}[1];   # total sum of contig lengths
    if ($hits{$_}[1] ne '0') { # sum of contigs is not zero
        $hits{$_}[2] = int($hits{$_}[1]/$hits{$_}[0]) if $_ ne 'nohit';  # average contig length
	$hits{$_}[3] = 100*($hits{$_}[0]/$contig_num); # percent of total contigs in each category
        $hits{$_}[3] =~ s/^(\S+\.\S\S)\S*/$1/;
    }
    $hits{$_}[2] = '0' if !$hits{$_}[2];
    $hits{$_}[3] = '0' if !$hits{$_}[3];
}

print  "Total queries with hits = $count; total contigs = $contig_num\n";
if (!$histo) {
print <<END;
No hits: $hits{nohit}[0] $hits{nohit}[3]%
Exact match: $hits{exact}[0] ($hits{exact}[2]) $hits{exact}[3]%
Full alignment with internal gaps/mismatches: $hits{internal}[0] ($hits{internal}[2]) $hits{internal}[3]%
 - Contigs with mismatches: $hits{fmismatch}[0] ($hits{fmismatch}[1];$hits{fmismatch}[2]) $hits{fmismatch}[3]%;
 - Contigs with insertions: $hits{fqgap}[0] ($hits{fqgap}[1];$hits{fqgap}[2]) $hits{fqgap}[3]%;
 - Contigs with deletions: $hits{ftgap}[0] ($hits{ftgap}[1];$hits{ftgap}[2]) $hits{ftgap}[3]%;

Partial alignment (<50 bp not aligned at ends) no internal gaps/mismatches: $hits{partial1}[0] ($hits{partial1}[2]) $hits{partial1}[3]%
Partial alignment (50 to 100 bp not aligned at ends) no internal gaps/mismatches: $hits{partial2}[0] ($hits{partial2}[2]) $hits{partial2}[3]%
Partial alignment (>=100 bp not aligned at end) no internal gaps/mismatches: $hits{partial3}[0] ($hits{partial3}[2]) $hits{partial3}[3]%

Partial alignment (<50 bp not aligned at ends) with internal gaps/mismatches: $hits{partial_gaps1}[0] ($hits{partial_gaps1}[2]) $hits{partial_gaps1}[3]%
Partial alignment (50 to 100 bp not aligned at ends) with internal gaps/mismatches: $hits{partial_gaps2}[0] ($hits{partial_gaps2}[2]) $hits{partial_gaps2}[3]%
Partial alignment (>=100 bp not aligned at ends) with internal gaps/mismatches: $hits{partial_gaps3}[0] ($hits{partial_gaps3}[2]) $hits{partial_gaps3}[3]%

- Contigs with mismatches: $hits{mismatch}[0] ($hits{mismatch}[1];$hits{mismatch}[2]) $hits{mismatch}[3]%;
- Contigs with insertions: $hits{qgap}[0] ($hits{qgap}[1];$hits{qgap}[2]) $hits{qgap}[3]%;
- Contigs with deletions: $hits{tgap}[0] ($hits{tgap}[1];$hits{tgap}[2]) $hits{tgap}[3]%;

END
}

else {
   print  join "\n", @diffs;

}

# for target coverage
my $tot=0;
my $tot_targetsize=0;
my $tot_hitsize=0;
#my %targethit=();
for (keys %targethits) {
    my $target = $_;
    my @coords=();
    # sort coords 
    for (0..$#{ $targethits{$target}} ) {
        push @coords,  $targethits{$target}[$_] ;
	#print STDERR "$target, $_, $targethits{$target}[$_][0] $targethits{$target}[$_][1]\n";
    }
    next if !@coords;
    my @sorted = sort { $a->[0] <=> $b->[0] } @coords;
    my $start = 0;
    my $end = 0;
    my $totalbases = 0;
    #print "Sorting\n";
    for (0..$#sorted) {
    #   print "$sorted[$_][0],$sorted[$_][1]\n";
       if (!$start && !$end) {
         $start = $sorted[$_][0];
         $end = $sorted[$_][1];
       }
       else {
           # there is overlap
           if ($sorted[$_][0] >= $start && $sorted[$_][0] <= $end  && $sorted[$_][1] > $end) {
               $end = $sorted[$_][1];
	   }
	   # there is no overlap
	   elsif ($sorted[$_][0] > $end) {
              $totalbases += ($end - $start +1);
              $start = $sorted[$_][0];
              $end = $sorted[$_][1];
	   }

       }
    } # end of for loop
    # add final from last loop
    $totalbases += ($end - $start +1);
    my $perc = 100*$totalbases/$targets{$target};
    #print  "$target $perc $totalbases $targets{$target}\n";
    $tot += $perc;
    $tot_targetsize += $targets{$target};
   # print $targets{$target}."\n"; # print out target hit sizes
    $tot_hitsize += $totalbases;
}
my $num = scalar (keys %targethits);
my $ave = $tot/$num;
my $ave_targetsize = $tot_targetsize/$num;
my $ave_hitsize = $tot_hitsize/$num;
print "Average coverage per target is $ave\n";
print "$num targets hit - ave. target size $ave_targetsize - ave. hit size $ave_hitsize\n";
print "Total target size $tot_targetsize\n";
print "Total $totbases matching; $tottgap target bases in alignment gaps\n";

exit;

my @sorted = sort {$b <=> $a} keys %diffs;
for (1..$sorted[0]) {
  if ($diffs{$_}) {
    print "$_\t$diffs{$_}\n";
    
  } 
  else {
    print "$_\t0\n";
  }

}

