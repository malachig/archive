#!/usr/bin/perl -w
# This script will merge a series of bed files from hmmSplicer into a single BED file
# Use this script to merge the data from a series of lanes for a single library into one file
# Malachi Griffith and Rodrigo Goya

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;
use Benchmark;

my $results_dir = '';
my $dir_depth = '';
my $summary_name = '';
my $track_name = '';
my $ref_genome = '';

GetOptions ('results_dir=s'=>\$results_dir, 'dir_depth=s'=>\$dir_depth, 'summary_name=s'=>\$summary_name, 'track_name=s'=>\$track_name, 'ref_genome=s'=>\$ref_genome);

if ($results_dir && $dir_depth && $summary_name && $track_name && $ref_genome){
  #Parameters ok
}else{
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  print GREEN, "\nSpecify the results directory to be summarized using: --results_dir", RESET;
  print GREEN, "\nSpecify the depth to search in this results dir for .bed files using: --dir_depth [1 | 2]", RESET;
  print GREEN, "\nSpecify the name of the grouping of data you are summarizing using: --summary_name (e.g. library id)", RESET;
  print GREEN, "\nSpecify the name that will be used for the BED track (should be unique) using: --track_name (e.g. library name)", RESET;
  print GREEN, "\n\nExample: mergeJunctionBed.pl  --results_dir=/projects/alexa2/hmmSplicer/SA_TN_Breast/results/SA017  --dir_depth=1  --summary_name=SA017  --track_name=SA017  --ref_genome=/projects/malachig/sequence_databases/hg18_genome/all_human.fa\n\n", RESET;
  exit();
}
#Check parameters
unless (-e $results_dir && -d $results_dir){
  print RED, "\n\nResults dir specified does not appear valid: $results_dir", RESET;
  exit();
}
unless ($dir_depth =~ /^1$|^2$|^3$/){
  print RED, "\n\nSupported values for --dir_depth are 1, 2 or 3\n\n", RESET;
  exit();
}
unless ($results_dir =~ /\/$/){
  $results_dir .= "/";
}
unless (-e $ref_genome){
  print RED, "\n\nCould not find reference genome file: $ref_genome\n\n", RESET;
  exit();
}

#Get a list of bed files by searching at the appropriate depth
my $search;
if ($dir_depth == 1){
  $search = "ls $results_dir"."*/junction.final.bed*";
}elsif($dir_depth == 2){
  $search = "ls $results_dir"."*/*/junction.final.bed*";
}elsif($dir_depth == 3){
  $search = "ls $results_dir"."*/results/*junction.final.bed*";
}

my @files = `$search`;
chomp(@files);

#Go through each file and store a master list of junctions
print BLUE, "\n\nProcessing files from: $search", RESET;
my %junctions;
foreach my $bed_file (@files){
  my $lc = 0; #line count
  my $vc = 0; #valid junction line count
  my $pipe;
  if ($bed_file =~ /\.gz$/){
    $pipe = "zcat $bed_file |";
  }elsif($bed_file =~ /\.bz2$/){
    $pipe = "bzcat $bed_file |";
  }else{
    $pipe = "$bed_file";
  }

  print BLUE, "\n\tProcessing: $bed_file", RESET;
  open(BED, "$pipe") || die "\n\nCould not open junction bed file: '$bed_file'\n\n";
  while(<BED>){
    $lc++;
    chomp($_); 
    s///;

    #chrY    2647717 2648420 SOLEXA7_25:4:94:6:1221/2|junc=2 1002.27266187   -       2647717 2648420 0       2       29,24,  0,679,
    my @f = split(/\t/, $_);
    if (!($f[0] =~ m/^chr/)) {next;}
    my $blocks = $f[9];
    unless ($blocks == 2) {print RED, "\n\nNumber of blocks is not 2!\n\n", RESET; next;}
    my $start = $f[1];
    my $end = $f[2];
    my $info = $f[3];
    my $score = $f[4];
    my $strand = $f[5];
    my @b_sizes = split(/,/, $f[10]);
    my @b_offsets = split(/,/, $f[11]);

    #Get the read count from the info field
    my $read_count;
    if ($info =~ /junc\=(\d+)/){
      $read_count = $1;
    }elsif($info =~ /rc\=(\d+)\|lbp\=\d+\|rbp\=\d+/){
      $read_count = $1;      
    }else{
      $read_count = 1;
    }

    for(my $b = 1; $b < $blocks; $b++){
      my $left = $start + $b_offsets[$b-1] + $b_sizes[$b-1];
      my $right = $start + $b_offsets[$b] + 1;

      #Keep track of info needed to determine the maximum sequence covered by the reads supporting each junction
      #Multiple reads corresponding to the same junction overlap the junction differently - we want to know the total overlap observed
      #We also want to know the average score assigned to the junction (rather than arbitrarily grabbing the first one observed)
      my $b_size1 = $b_sizes[$b-1];
      my $b_size2 = $b_sizes[$b];
      my $b_offset1 = $b_offsets[$b-1];
      my $b_offset2 = $b_offsets[$b];

      #Keep track of the strand counts for each junction

      #Primary key for junctions will look like: chrX:1000-1040 (Note that we can't use strandedness since the Illumina libraries are not stranded)
      #my $j = "$f[0]($strand):$left-$right"; #Stranded
      my $j = "$f[0]:$left-$right"; #NOT Stranded


      #Determine the intron size for the junction
      my $intron_size = ($right-$left)+1;

      #Store the junction.  If it has been observed before, increment the count of the existing record
      if ($junctions{$j}){
        $junctions{$j}{read_count} += $read_count;
        push(@{$junctions{$j}{scores}}, $score);

        $junctions{$j}{strand} = ".";
        if ($start < $junctions{$j}{min_start}){
          $junctions{$j}{min_start} = $start;
          $junctions{$j}{max_b_size1} = $b_size1;
          $junctions{$j}{final_offset2} = $b_offset2;
        }
        if ($end > $junctions{$j}{max_end}){
          $junctions{$j}{max_end} = $end;
          $junctions{$j}{max_b_size2} = $b_size2;
        }
      }else{
        my @scores;
        push(@scores, $score);
        $junctions{$j}{line} = \@f;
        $junctions{$j}{read_count} = $read_count;
        $junctions{$j}{scores} = \@scores;
        $junctions{$j}{min_start} = $start;
        $junctions{$j}{max_end} = $end;
        $junctions{$j}{final_offset1} = $b_offset1;
        $junctions{$j}{final_offset2} = $b_offset2;
        $junctions{$j}{max_b_size1} = $b_size1;
        $junctions{$j}{max_b_size2} = $b_size2;
        $junctions{$j}{intron_size} = $intron_size;
      }
    }
    $vc++;
  }
  close(BED);

  my $unique_junctions = keys %junctions;
  print BLUE, "\n\t\tFound $lc records ($vc valid junction records) -> unique junctions stored thus far = $unique_junctions", RESET;

}


#Determine strand by comparison back to the reference genome...
print BLUE, "\n\nAttempting to determine correct strand for each junction", RESET;
open (REF, "$ref_genome") || die "\n\nCould not open fasta file: $ref_genome";
my $tmp = $/;
$/ = "\n>";  # read by FASTA record

while (<REF>){
  chomp $_;
  my $chr_seq = $_;
  my ($chr) = $chr_seq =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
  $chr_seq =~ s/^>*.+\n//;  # remove FASTA header
  $chr_seq =~ s/\n//g;  # remove endlines

  my $chr_length = length($_);
  print BLUE, "\n\tFound $chr sequence (length = $chr_length)", RESET;

  #Now go through the junctions found for this chromosome and look for donor/acceptor splice sites at the coordinates reported
  #SPLICE_SITES = ["GT-AG", "CT-AC", "GC-AG", "CT-GC", "AT-AC", "GT-AT"]

  foreach my $j (sort keys %junctions){
    if ($j =~ /(.*)\:(\d+)\-(\d+)/){
      my $j_chr = $1;
      my $left = $2;
      my $right = $3;
      unless($chr eq $j_chr){
        next();
      }
      my $left_dn = uc(substr($chr_seq, $left, 2));
      my $right_dn = uc(substr($chr_seq, $right-3, 2));

      #print "\n\t\tDEBUG: $left_dn ... $right_dn";

      #Assign strand...
      if ($left_dn eq "GT" && $right_dn eq "AG"){
        $junctions{$j}{strand} = "+";
        $junctions{$j}{splice_site} = "GT-AG";
      }elsif($left_dn eq "CT" && $right_dn eq "AC"){
        $junctions{$j}{strand} = "-";
        $junctions{$j}{splice_site} = "GT-AG";
      }elsif($left_dn eq "GC" && $right_dn eq "AG"){
        $junctions{$j}{strand} = "+";
        $junctions{$j}{splice_site} = "GC-AG";
      }elsif($left_dn eq "CT" && $right_dn eq "GC"){
        $junctions{$j}{strand} = "-";
        $junctions{$j}{splice_site} = "GC-AG";
      }elsif($left_dn eq "AT" && $right_dn eq "AC"){
        $junctions{$j}{strand} = "+";
        $junctions{$j}{splice_site} = "AT-AC";     
      }elsif($left_dn eq "GT" && $right_dn eq "AT"){
        $junctions{$j}{strand} = "-";
        $junctions{$j}{splice_site} = "AT-AC";
      }else{
        print RED, "\n\nUnable to determine strand!!\n\n", RESET;
        exit();
      }

      #Get the actual junction sequence while we are at it...
      my $left_size = $junctions{$j}{max_b_size1};
      my $right_size = $junctions{$j}{max_b_size2};
      my $left_seq = substr($chr_seq, $left-($left_size), $left_size);
      my $right_seq = substr($chr_seq, $right-1, $right_size);
      #my $ls = length($left_seq);
      #my $rs = length($right_seq);
      #print YELLOW, "\n$left_size ($ls)\t$right_size ($rs)", RESET;
      my $seq = "$left_seq"."$right_seq";
      $junctions{$j}{seq} = $seq;
    }else{
      print RED, "\n\nObserved junction not understood\n\n", RESET;
      exit();
    }

  }
}
close(REF);
$/ = $tmp;



#Write out a new BED file and add an updated track line using the summary name in the track name
#e.g.  track name='$track_name' description='HMMSplicer Canonical Junctions' useScore=1

#At the same time, write a junction ID file in simple tabular format with the read counts and scores only
my $outfile1 = "$results_dir"."$summary_name".".junction.final.bed";
my $outfile2 = "$results_dir"."$summary_name".".junction.list.txt";
print BLUE, "\n\nPrinting merged junction bed data to new bed file: $outfile1", RESET;
print BLUE, "\nAlso printing simple junction list to: $outfile2", RESET;
open (OUT1, ">$outfile1") || die "\n\nCould not open output file: $outfile1\n\n";
open (OUT2, ">$outfile2") || die "\n\nCould not open output file: $outfile2\n\n";
print OUT1 "browser position chr3:125931903-125946730\n";
print OUT1 "browser pack knownGene ensGene\n";
print OUT1 "track name='$track_name' description='HMMSplicer Canonical Junctions for $track_name ($summary_name)' useScore=1 visibility=3\n";
print OUT2 "JID\tRead_Count\tScore\tLeft_Size\tRight_Size\tIntron_Size\tSplice_Site\tSequence\n";
foreach my $j (sort keys %junctions){
  my @original_line = @{$junctions{$j}{line}};
  my @adjusted_line = @original_line;

  #Get average score
  my @scores = @{$junctions{$j}{scores}};
  my $n = scalar(@scores);
  my $cum_score = 0;
  foreach my $s (@scores){
    $cum_score += $s;
  }
  my $avg_score = sprintf("%.8f", ($cum_score/$n));

  #chr1    100088542       100089151       read_count=3    1024.14811807   +       100088542       100089151       0       2       17,33,  0,576,

  #Create an adjusted BED line that merges data from redundant observations of the same junction in different lanes of data

  #Replace name field with read count
  $adjusted_line[3] = "rc=$junctions{$j}{read_count}|lbp=$junctions{$j}{max_b_size1}|rbp=$junctions{$j}{max_b_size2}";
  local $" = "\t";

  #Replace score with average score
  $adjusted_line[4] = $avg_score;

  #Replace start, end, block sizes, and offsets
  $adjusted_line[1] = $junctions{$j}{min_start};
  $adjusted_line[2] = $junctions{$j}{max_end};
  $adjusted_line[5] = $junctions{$j}{strand};
  $adjusted_line[6] = $junctions{$j}{min_start};
  $adjusted_line[7] = $junctions{$j}{max_end};
  $adjusted_line[10] = "$junctions{$j}{max_b_size1},$junctions{$j}{max_b_size2}";
  $adjusted_line[11] = "$junctions{$j}{final_offset1},$junctions{$j}{final_offset2}";


  #Print out the adjusted BED line
  print OUT1 "@adjusted_line\n";

  #Print out the junction list file
  print OUT2 "$j($junctions{$j}{strand})\t$junctions{$j}{read_count}\t$avg_score\t$junctions{$j}{max_b_size1}\t$junctions{$j}{max_b_size2}\t$junctions{$j}{intron_size}\t$junctions{$j}{splice_site}\t$junctions{$j}{seq}\n";
}
close(OUT1);
close(OUT2);

#Compress the bed file
my $gzip_cmd = "gzip -f $outfile1";
print BLUE, "\n\nCompressing output bed file:\n$gzip_cmd", RESET;
system($gzip_cmd);

#Sort the junction list list
print BLUE, "\n\nSorting and overwritting the junction list file", RESET;
my $outfile3 = "$results_dir"."$summary_name".".junction.list.tmp";
my $sort_cmd = "sort $outfile2 > $outfile3";
my $mv_cmd = "mv -f $outfile3 $outfile2";
system($sort_cmd);
system($mv_cmd);

print "\n\n";
exit();
