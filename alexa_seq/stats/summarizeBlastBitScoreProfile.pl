#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Based on a tab-delimited blast results file, help the user determine a suitable bit score cutoff by summarizing scores
#Display the score that results for each hit alignment length with increasing numbers of mismatches or gaps
#Obviously the outcome will be influenced by the BLAST actually conducted by the user to create the results file in the first place
#Specifically bit score is affected by the size of the search space (WORD SIZE for example)

#Note that the Bit Score for a hit with a gap depends on the size of the gap
#When summarizing, favor the smallest gaps

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my $blast_file = '';
GetOptions ('blast_file=s'=>\$blast_file);

print GREEN, "\n\nExample usage:\n\nsummarizeBlastBitScoreProfile.pl  --blast_file=blast_results.txt OR --blast_file=blast_results.txt.gz\n\n", RESET;

unless ($blast_file){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Check validity of input files
unless (-e $blast_file){
  print RED, "\nInput blast file does not appear to be invalid!\n\n", RESET;
  exit();
}

if ($blast_file =~ /\.txt$/){
  open (BLAST, "$blast_file") || die "\nCould not open BLAST input file\n\n";
}elsif($blast_file =~ /\.gz$/){
  open (BLAST, "zcat $blast_file |") || die "\nCould not open BLAST input file\n\n";
}else{
  print RED, "\n\nFile extension not recognized\n\n", RESET;
  exit();
}

my %blast;

while (<BLAST>){
  chomp($_);
  my @line = split("\t", $_);
  my $fields = scalar(@line);

  unless($fields == 12){
    print YELLOW, "\nFound an invalid blast line\n\n", RESET;
    next();
  }

  my $percent_identity = $line[2];
  my $length = $line[3];
  my $mismatches = $line[4];
  my $gaps = $line[5];
  my $bit_score = $line[11];

  #Skip hits with both mismatches and gaps so that they can be considered seperately
  if ($mismatches > 0 && $gaps > 0){
    next();
  }

  unless ($gaps > 0){
    if ($blast{$length}{mismatches}){

      my $mm_ref = $blast{$length}{mismatches};

      if ($mm_ref->{$mismatches}){

	my $bit_score_ref = $mm_ref->{$mismatches}->{bit_score};
	if ($bit_score_ref->{$bit_score}){
	  $bit_score_ref->{$bit_score}->{count}++;
	}else{
	  $bit_score_ref->{$bit_score}->{count} = 1;
	}
      }else{
	my %bit_score;
	$bit_score{$bit_score}{count} = 1;
	$mm_ref->{$mismatches}->{bit_score} = \%bit_score;
      }


    }else{
      my %mm;
      my %bit_score;

      $bit_score{$bit_score}{count} = 1;

      $mm{$mismatches}{bit_score} = \%bit_score;
      $blast{$length}{mismatches} = \%mm;
    }
  }

  unless ($mismatches > 0){
    if ($blast{$length}{gaps}){
      my $gap_ref = $blast{$length}{gaps};

      if ($gap_ref->{$gaps}){

	#Exchange the previous bit score if this hit has the same number of gaps but a higher percent id (i.e. a smaller gap)
	if ($percent_identity > $gap_ref->{$gaps}->{best_percent_identity}){
	  my %bit_score;
	  $bit_score{$bit_score}{count} = 1;

	  $gap_ref->{$gaps}->{bit_score} = \%bit_score;
	  $gap_ref->{$gaps}->{best_percent_identity} = $percent_identity;

	}else{
	  #Percent id is the same.

	  my $bit_score_ref = $gap_ref->{$gaps}->{bit_score};

	  if ($bit_score_ref->{$bit_score}){
	    $bit_score_ref->{$bit_score}->{count}++;
	  }else{
	    $bit_score_ref->{$bit_score}->{count} = 1;
	  }

	  $gap_ref->{$gaps}->{best_percent_identity} = $percent_identity;

	}

      }else{
	#First time this number of gaps has been observed for this length
	my %bit_score;
	$bit_score{$bit_score}{count} = 1;

	$gap_ref->{$gaps}->{bit_score} = \%bit_score;
	$gap_ref->{$gaps}->{best_percent_identity} = $percent_identity;

      }

    }else{
      #First time this hit length has been observed
      my %gaps;
      my %bit_score;
      $bit_score{$bit_score}{count} = 1;

      $gaps{$gaps}{bit_score} = \%bit_score;
      $gaps{$gaps}{best_percent_identity} = $percent_identity;

      $blast{$length}{gaps} = \%gaps;
    }
  }

}

#Now print out the outcome in a friendly format
print BLUE, "\n\nSummary of bit scores observed by hit length and according to number of mismatches or gaps found", RESET;
print BLUE, "\n\tNote that hits with both gaps and mismatches are ignored here", RESET; 

print BLUE, "\n\nMisMatches (# mismatches followed by score in brackets)", RESET;
foreach my $length (sort {$a <=> $b} keys %blast){
  print YELLOW, "\nLENGTH: $length", RESET;

  my $mm_ref = $blast{$length}{mismatches};
  foreach my $mm (sort {$a <=> $b} keys %{$mm_ref}){

    #Get the bit scores associated with hits of this length with this many mismatches and print only the most common one
    my $bs_ref = $mm_ref->{$mm}->{bit_score};

    foreach my $bs (sort {$bs_ref->{$b}->{count} <=> $bs_ref->{$a}->{count}} keys %{$bs_ref}){

      print YELLOW, "  $mm ($bs)", RESET;
      last();
    }
  }
}

print BLUE, "\n\nGaps (# gaps followed by score in brackets)", RESET;
foreach my $length (sort {$a <=> $b} keys %blast){
  print YELLOW, "\nLENGTH: $length", RESET;

  my $gap_ref = $blast{$length}{gaps};
  foreach my $gap (sort {$a <=> $b} keys %{$gap_ref}){

    #Get the bit scores associated with hits of this length with this many gaps and print only the most common one
    my $bs_ref = $gap_ref->{$gap}->{bit_score};

    foreach my $bs (sort {$bs_ref->{$b}->{count} <=> $bs_ref->{$a}->{count}} keys %{$bs_ref}){
      print YELLOW, "  $gap ($bs)", RESET;
      last();
    }
  }
}

print "\n\n";

close (BLAST);

exit();

