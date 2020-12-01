#!/usr/bin/perl -w
# Malachi Griffith

# Determines the following for all observed junctions in a project (or all projects)
#1.) Junctions that can be anchored to a known junction (matching donor OR acceptor) - consider strand here 
#    - If the Donor matches to a gene on the positive strand, the non matching Acceptor should have a higher chromosome coordinate
#    - If the Donor matches to a gene on the negative strand, the non matching Acceptor should have a lower chromosome coordinate
#    - If the Acceptor matches to a gene on the positive strand, the non-matching Donor should have a lower chromosome coordinate
#    - If the Acceptor matches to a gene on the negative strand, the non-matching Donor should have a higher chromosome coordinate
#2.) Junctions that could not be anchored at all.  
#3.) Exon skipping 
#    - Number of known splice sites skipped
#    - Number of known exon clusters skipped
#4.) Transcript IDs of matching transcripts (where both Acceptor and Donor match, or if no match, then the corresponding anchored transcripts)...
#    - CCDS, Ensembl, MGC, Refseq, UCSC, Vega
#5.) Gene IDs of matching transcripts (non-redundant list)
#6.) Determine the number of exons skipped and exon donors/acceptors skipped by each junction
#7.) Join Splign results onto the junction annotations (splign validated, number of splign junctions)

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;
use Benchmark;

my $bedtools_bin_dir = '';
my $analysis_dir = '';
my $project = '';
my $ref_junction_file = '';
my $ref_ec_file = '';
my $splign_file = '';

GetOptions ('bedtools_bin_dir=s'=>\$bedtools_bin_dir, 'analysis_dir=s'=>\$analysis_dir, 'project=s'=>\$project, 
            'ref_junction_file=s'=>\$ref_junction_file, 'ref_ec_file=s'=>\$ref_ec_file, 'splign_file=s'=>\$splign_file);

if ($bedtools_bin_dir && $analysis_dir && $project && $ref_junction_file && $ref_ec_file && $splign_file){
  #Parameters ok
}else{
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  print GREEN, "\nSpecify the path to the BEDTools binary dir using: --bedtools_bin_dir", RESET;
  print GREEN, "\nSpecify the base analysis dir using: --analysis_dir", RESET;
  print GREEN, "\nSpecify the project name using: --project", RESET;
  print GREEN, "\nSpecify the reference junction file using: --ref_junction_file", RESET;
  print GREEN, "\nSpecify the reference exon content file (merge of overlapping exons on each strand) exons: --ref_ec_file", RESET;
  print GREEN, "\nSpecify the master splign results file using: --splign_file", RESET;
  print GREEN, "\n\nExample: annotateObservedJunctions.pl  --bedtools_bin_dir=/home/malachig/tools/BEDTools-Version-2.10.1/bin/  --analysis_dir=/projects/alexa2/hmmSplicer/  --project=SA_TN_Breast  --ref_junction_file=/projects/alexa2/hmmSplicer/ReferenceAnnotations/hg18/ALL.junc  --ref_ec_file=/projects/alexa2/hmmSplicer/ReferenceAnnotations/hg18/ALL.ExonContent  --splign_file=/projects/alexa2/hmmSplicer/SplignValidation/hg18/splignMasterResults.txt\n\n", RESET;
  exit();
}
unless ($analysis_dir =~ /\/$/){
  $analysis_dir .= "/";
}
unless (-e $analysis_dir && -d $analysis_dir){
  print RED, "\n\nAnalysis dir: $analysis_dir not found\n\n", RESET;
  exit();
}
unless ($bedtools_bin_dir =~ /\/$/){
  $bedtools_bin_dir .= "/";
}
unless (-e $bedtools_bin_dir && -d $bedtools_bin_dir){
  print RED, "\n\nBEDTools binary dir: $bedtools_bin_dir not found\n\n", RESET;
  exit();
}
unless (-e $ref_junction_file){
  print RED, "\n\nCould not find ref_junction_file: $ref_junction_file\n\n", RESET;
  exit();
}
unless (-e $ref_ec_file){
  print RED, "\n\nCould not find ref_ec_file: $ref_ec_file\n\n", RESET;
  exit();
}
unless (-e $splign_file){
  print RED, "\n\nCould not find splign_file: $splign_file\n\n", RESET;
  exit();
}


#Define the observed junctions file and outfile based on the analysis dir and project name
my $obs_junction_file = "$analysis_dir$project/results/$project".".junction.list.txt";
my $outfile = "$analysis_dir$project/results/$project".".junction.list.anno.txt";

#Make sure the merge project junction file can be found
unless (-e $obs_junction_file){
  print RED, "\n\nObserved junction file ($obs_junction_file) could not be found\n\n", RESET;
  exit();
}

#Import the reference junctions (from Ensembl + UCSC + MGC + Refseq + Vega + CCDS)
my %known_junctions;
my %known_donors;
my %known_acceptors;
&importRefJunctions('-ref_junction_file'=>$ref_junction_file);

#Import the reference exons (actually merged exon content blocks (from Ensembl + UCSC + MGC + Refseq + Vega + CCDS)
my %known_ec_blocks;
&importEcBlocks('-ref_ec_file'=>$ref_ec_file);

#Import all observed junctions for the project using the merged project junction file
my %observed_junctions;
&importObsJunctions('-obs_junction_file'=>$obs_junction_file);

#Determine the exon and splice site skipping of each junction observed (if it known, or anchored to a known splice site at one or both ends)
#&annotateSkipping();
&annotateSkippingBT();

# Join Splign results onto the junction annotations (splign validated, number of splign junctions)
&joinSplignResults('-infile'=>$splign_file);

open(OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
print OUT "JID\tRead_Count\tScore\tLeft_Size\tRight_Size\tIntron_Size\tSplice_Site\tAnchored\tExons_Skipped\tDonors_Skipped\tAcceptors_Skipped\tSplign_Validated\tSplign_Align_Count\tGene_Name\tSequence\n";
foreach my $jid (sort {$observed_junctions{$a}{order} <=> $observed_junctions{$b}{order}} keys %observed_junctions){
  my $string = "$jid\t$observed_junctions{$jid}{read_count}\t$observed_junctions{$jid}{score}\t$observed_junctions{$jid}{left_size}\t$observed_junctions{$jid}{right_size}\t$observed_junctions{$jid}{intron_size}\t$observed_junctions{$jid}{splice_site}\t$observed_junctions{$jid}{anchored}\t$observed_junctions{$jid}{exons_skipped}\t$observed_junctions{$jid}{donors_skipped}\t$observed_junctions{$jid}{acceptors_skipped}\t$observed_junctions{$jid}{splign_validated}\t$observed_junctions{$jid}{splign_alignment_count}\t$observed_junctions{$jid}{gid_list}\t$observed_junctions{$jid}{seq}\n";

  #print "$string";
  print OUT "$string";
}
close(OUT);

print BLUE, "\n\nPrinted resulting annotated junctions to: $outfile", RESET;

print "\n\n";

exit();


######################################################################################################################################
#Import the reference junctions (e.g. Ensembl + UCSC + MGC + Refseq + Vega + CCDS)                                                   #
######################################################################################################################################
sub importRefJunctions{
  my %args = @_;
  my $infile = $args{'-ref_junction_file'};

  print BLUE, "\n\nImporting reference junctions", RESET;
  #Inport the ref junction file.  In this file, the coordinates are always ordered regardless of strand
  #i.e. the 'left' coordinate is always a smaller number than the 'right' coordinate
  #If the strand is '+' then Donor=left and Acceptor=right
  #If the strand is '-' then Donor=right and Acceptor=left
  open(JUNC, "$infile") || die "\n\nCould not open input file: $infile\n\n";
  my $c = 0;
  while(<JUNC>){
    chomp($_);
    my @line = split("\t", $_);
    my $jid = $line[0];
    my $chr = $line[1];
    my $left = $line[2];
    my $right = $line[3];
    my $strand = $line[4];
    my $count = $line[5];
    my $gid_list = $line[6];

    $known_junctions{$jid}=$gid_list;
    $c++;

    if ($strand eq "+"){
      $known_donors{$chr}{$strand}{$left}=$gid_list;
      $known_acceptors{$chr}{$strand}{$right}=$gid_list;

    }elsif($strand eq "-"){
      $known_donors{$chr}{$strand}{$right}=$gid_list;
      $known_acceptors{$chr}{$strand}{$left}=$gid_list;
    }else{
      print RED, "\n\nUnknown strand in reference junctions file", RESET;
      exit();
    }
  }
  close(JUNC);

  print BLUE, "\n\tImported $c known junctions and associated acceptor/donor pairs\n", RESET;

  #foreach my $pos (sort {$a <=> $b} keys %{$test{'chr1'}{'+'}}){
  #  print "\n$pos";
  #}
  #exit();

  return();
}


######################################################################################################################################
#Import the reference exons (actually merged exon content blocks (from Ensembl + UCSC + MGC + Refseq + Vega + CCDS)                  #
######################################################################################################################################
sub importEcBlocks{
  my %args = @_;
  my $infile = $args{'-ref_ec_file'};

  print BLUE, "\n\nImporting reference exon content blocks", RESET;

  #Import the ref exon content block file.  In this file, the coordinates are always ordered regardless of strand
  #i.e. the 'left' coordinate is always a smaller number than the 'right' coordinate
  open(EXON, "$infile") || die "\n\nCould not open input file: $infile\n\n";
  my $c = 0;
  while(<EXON>){
    chomp($_);
    my @line = split("\t", $_);


    my $chr = $line[0];
    my $left = $line[1];
    my $right = $line[2];
    my $count = $line[3];
    my $strand = $line[4];

    $c++;

    if ($known_ec_blocks{$chr}{$strand}){
      my $ec_ref = $known_ec_blocks{$chr}{$strand};
      $ec_ref->{$c}->{left} = $left;
      $ec_ref->{$c}->{right} = $right;
    }else{
      my %tmp;
      $tmp{$c}{left} = $left;
      $tmp{$c}{right} = $right;
      $known_ec_blocks{$chr}{$strand} = \%tmp;
    }
  }
  close(EXON);
  print BLUE, "\n\tImported $c known exon content blocks\n", RESET;

  return();
}


######################################################################################################################################
#Import the observed junctions from the specified file                                                                               #
######################################################################################################################################
sub importObsJunctions{
  my %args = @_;
  my $infile = $args{'-obs_junction_file'};

  print BLUE, "\n\nImporting observed junctions", RESET;
  #Inport the obs junction file.  In this file, the coordinates are always ordered regardless of strand
  #i.e. the 'left' coordinate is always a smaller number than the 'right' coordinate
  #The strand is unknown...
  open(JUNC, "$infile") || die "\n\nCould not open input file: $infile\n\n";
  my $c = 0;
  my $header = 1;
  my %columns;

  while(<JUNC>){
    chomp($_);
    my @line = split("\t", $_);

    if ($header){
      my $pos = 0;
      foreach my $head (@line){
        $columns{$head}{pos} = $pos;
        $pos++;
      }
      $header = 0;
      next();
    }

    my $jid = $line[$columns{'JID'}{pos}];
    my $read_count = $line[$columns{'Read_Count'}{pos}];
    my $score = $line[$columns{'Score'}{pos}];
    my $left_size = $line[$columns{'Left_Size'}{pos}];
    my $right_size = $line[$columns{'Right_Size'}{pos}];
    my $intron_size = $line[$columns{'Intron_Size'}{pos}];
    my $splice_site = $line[$columns{'Splice_Site'}{pos}];
    my $seq = $line[$columns{'Sequence'}{pos}];

    my $chr;
    my $left;
    my $right;
    my $strand;
    if ($jid =~ /(.*)\:(\d+)\-(\d+)\((.*)\)/){
      $chr = $1;
      $left = $2;
      $right = $3;
      $strand = $4;
    }else{
      print RED, "\n\nObserved junction not understood\n\n", RESET;
      exit();
    }

    $observed_junctions{$jid}{read_count} = $read_count;
    $observed_junctions{$jid}{score} = $score;
    $observed_junctions{$jid}{left_size} = $left_size;
    $observed_junctions{$jid}{right_size} = $right_size;
    $observed_junctions{$jid}{intron_size} = $intron_size;
    $observed_junctions{$jid}{splice_site} = $splice_site;
    $observed_junctions{$jid}{seq} = $seq;
    $observed_junctions{$jid}{anchored} = "N";
    $observed_junctions{$jid}{gid_list} = "na";
    $c++;

    #First check for an exact match to a known junction (i.e. anchored by both Donor and Acceptor)
    $observed_junctions{$jid}{order} = $c;
    if ($known_junctions{$jid}){
      $observed_junctions{$jid}{anchored} = "DA";
      $observed_junctions{$jid}{gid_list} = $known_junctions{$jid};
    }else{
      #Now check for anchoring on one side only
      #    - If the Donor matches to a gene on the positive strand, the non matching Acceptor should have a higher chromosome coordinate
      #    - If the Donor matches to a gene on the negative strand, the non matching Acceptor should have a lower chromosome coordinate
      #    - If the Acceptor matches to a gene on the positive strand, the non-matching Donor should have a lower chromosome coordinate
      #    - If the Acceptor matches to a gene on the negative strand, the non-matching Donor should have a higher chromosome coordinate
      if ($strand eq "+"){
        if($known_donors{$chr}{'+'}{$left} && $known_acceptors{$chr}{'+'}{$right}){
          $observed_junctions{$jid}{anchored} = "NDA";
          $observed_junctions{$jid}{gid_list} = $known_donors{$chr}{'+'}{$left};
        }elsif($known_donors{$chr}{'+'}{$left}){
          $observed_junctions{$jid}{anchored} = "D";
          $observed_junctions{$jid}{gid_list} = $known_donors{$chr}{'+'}{$left};
        }elsif($known_acceptors{$chr}{'+'}{$right}){
          $observed_junctions{$jid}{anchored} = "A";
          $observed_junctions{$jid}{gid_list} = $known_acceptors{$chr}{'+'}{$right};
        }
      }elsif ($strand eq "-"){
        if($known_donors{$chr}{'-'}{$right} && $known_acceptors{$chr}{'-'}{$left}){
          $observed_junctions{$jid}{anchored} = "NDA";
          $observed_junctions{$jid}{gid_list} = $known_donors{$chr}{'-'}{$right};          
        }elsif ($known_donors{$chr}{'-'}{$right}){
          $observed_junctions{$jid}{anchored} = "D";
          $observed_junctions{$jid}{gid_list} = $known_donors{$chr}{'-'}{$right};
        }elsif($known_acceptors{$chr}{'-'}{$left}){
          $observed_junctions{$jid}{anchored} = "A";
          $observed_junctions{$jid}{gid_list} = $known_acceptors{$chr}{'-'}{$left};
        }
      }
    }

    #Initialize the skipping values
    my $anchored = $observed_junctions{$jid}{anchored};

    if ($anchored eq "N"){
      $observed_junctions{$jid}{exons_skipped} = "na";
      $observed_junctions{$jid}{donors_skipped} = "na";
      $observed_junctions{$jid}{acceptors_skipped} = "na";
    }else{
      $observed_junctions{$jid}{exons_skipped} = 0;
      $observed_junctions{$jid}{donors_skipped} = 0;
      $observed_junctions{$jid}{acceptors_skipped} = 0;
    }

  }
  close(JUNC);

  print BLUE, "\n\tImported $c observed junctions\n", RESET;
  return();
}


##############################################################################################################################################
#Determine the exon and splice site skipping of each junction observed (if it known, or anchored to a known splice site at one or both ends) #
##############################################################################################################################################
sub annotateSkippingBT{
  print BLUE, "\n\nAnnotating skipping of each junction - exons, then acceptors, then donors - Using BEDTools\n", RESET;

  #Use BEDTools to determine the exons, acceptors, or donors contained within each exon-exon observed junction (i.e. intron)
  #Do this by writing a temp BED file for each pair of coordinates (of the form: chr, start, end, strand) and then parsing the results file

  #'intersectBed -a exons.bed -b junctions.bed -f 1.0 -s -wa -wb'  
  #The '-f 1.0' option should give the exons that are entirely overlapped by junctions
  #The '-s' option should make this overlap between things on the same strand
  #The '-wa -wb' options, write the original coordinates for exons and junctions (as opposed to the merged coordinates).  Each overlaping pair will be reported as a seperate line
  #Then process the output file and count the exon entries associated with each junction entry
  
  #Example BEDTools command
  #Run BEDTools as follows and 'cut' the column containing the junction ID value.  
  #The occurence count for each of these IDs, should be the number of exons contained within the corresponding junction (i.e. skipped)
  #Unix sort and uniq can even do this counting for us...
  #Then just parse the counts out

  my $temp_obs_junctions = "$analysis_dir$project/results/ObsJunc.tmp.bed";
  my $temp_known_exons = "$analysis_dir$project/results/KnownExonContent.tmp.bed";
  my $temp_known_donors = "$analysis_dir$project/results/KnownDonors.tmp.bed";
  my $temp_known_acceptors = "$analysis_dir$project/results/KnownAcceptors.tmp.bed";
  my $temp_result = "$analysis_dir$project/results/Result.tmp.txt";

  #OBSERVED JUNCTIONS
  print BLUE, "\n\tLooking for entire exons skipped", RESET;
  open (TMP_OJ, ">$temp_obs_junctions") || die "\n\nCould not open output file: $temp_obs_junctions";
  #Print print out the observed hmmSplicer junctions (i.e. introns as a temp bed file
  foreach my $jid (sort {$observed_junctions{$a}{order} <=> $observed_junctions{$b}{order}} keys %observed_junctions){
    my $chr;
    my $left;
    my $right;
    my $strand;
    if ($jid =~ /(.*)\:(\d+)\-(\d+)\((.*)\)/){
      $chr = $1;
      $left = $2;
      $right = $3;
      $strand = $4;
    }else{
      print RED, "\n\nObserved junction not understood\n\n", RESET;
      exit();
    }
    my $score = $observed_junctions{$jid}{score};
    print TMP_OJ "$chr\t$left\t$right\t$jid\t$score\t$strand\n";
  }
  close (TMP_OJ);

  #ENTIRE EXONS
  #Print out all the known exon content blocks
  open (TMP_KE, ">$temp_known_exons") || die "\n\nCould not open output file: $temp_known_exons";
  foreach my $chr (sort keys %known_ec_blocks){
    foreach my $strand (sort keys %{$known_ec_blocks{$chr}}){
      my $ec_ref = $known_ec_blocks{$chr}{$strand};

      foreach my $ec (sort {$ec_ref->{$a}->{left} <=> $ec_ref->{$b}->{left}} keys %{$ec_ref}){
        my $ec_left = $ec_ref->{$ec}->{left};
        my $ec_right = $ec_ref->{$ec}->{right};
        print TMP_KE "$chr\t$ec_left\t$ec_right\tEC\t.\t$strand\n";

      }
    }
  }
  close (TMP_KE);

  my $bed_cmd1 = "$bedtools_bin_dir"."intersectBed -a $temp_known_exons -b $temp_obs_junctions -f 1.0 -s -wa -wb | cut -f 10 | sort | uniq -c > $temp_result";
  print BLUE, "\n\t$bed_cmd1", RESET;
  system($bed_cmd1);
  open (COUNTS, "$temp_result") || die "\n\nCould not open temp results file: $temp_result\n\n";
  while(<COUNTS>){
    chomp($_);
    if ($_ =~ /(\d+)\s+(.*)/){
      my $count = $1;
      my $jid = $2;
      if ($observed_junctions{$jid}){
        my $anchored = $observed_junctions{$jid}{anchored};
        unless ($anchored eq "N"){
          $observed_junctions{$jid}{exons_skipped} = $count;
        }
      }else{
        print RED, "\n\nUnrecognized junction id: $jid\n\n", RESET;
        exit();
      }
    }else{
      print RED, "\n\nEntry in results file not understood: $_\n\n", RESET;
      exit();
    }
  }
  close(COUNTS);


  #DONORS
  print BLUE, "\n\tLooking for donors skipped", RESET;
  open (TMP_KD, ">$temp_known_donors") || die "\n\nCould not open output file: $temp_known_donors";
  foreach my $chr (sort keys %known_donors){
    foreach my $strand (sort keys %{$known_donors{$chr}}){
      foreach my $donor (sort keys %{$known_donors{$chr}{$strand}}){
        my $donor_p = $donor+1;
        print TMP_KD "$chr\t$donor\t$donor_p\tD\t.\t$strand\n";
      }
    }
  }
  close (TMP_KD);

  my $bed_cmd2 = "$bedtools_bin_dir"."intersectBed -a $temp_known_donors -b $temp_obs_junctions -f 1.0 -s -wa -wb | cut -f 10 | sort | uniq -c > $temp_result";
  print BLUE, "\n\t$bed_cmd2", RESET;
  system($bed_cmd2);
  open (COUNTS, "$temp_result") || die "\n\nCould not open temp results file: $temp_result\n\n";
  while(<COUNTS>){
    chomp($_);
    if ($_ =~ /(\d+)\s+(.*)/){
      my $count = $1;
      my $jid = $2;
      if ($observed_junctions{$jid}){
        my $anchored = $observed_junctions{$jid}{anchored};
        unless ($anchored eq "N"){
          $observed_junctions{$jid}{donors_skipped} = $count;
        }
      }else{
        print RED, "\n\nUnrecognized junction id: $jid\n\n", RESET;
        exit();
      }
    }else{
      print RED, "\n\nEntry in results file not understood: $_\n\n", RESET;
      exit();
    }
  }
  close(COUNTS);


  #ACCEPTORS
  print BLUE, "\n\tLooking for acceptors skipped", RESET;
  open (TMP_KA, ">$temp_known_acceptors") || die "\n\nCould not open output file: $temp_known_acceptors";
  foreach my $chr (sort keys %known_acceptors){
    foreach my $strand (sort keys %{$known_acceptors{$chr}}){
      foreach my $acceptor (sort keys %{$known_acceptors{$chr}{$strand}}){
        my $acceptor_p = $acceptor+1;
        print TMP_KA "$chr\t$acceptor\t$acceptor_p\tA\t.\t$strand\n";
      }
    }
  }
  close (TMP_KA);

  my $bed_cmd3 = "$bedtools_bin_dir"."intersectBed -a $temp_known_acceptors -b $temp_obs_junctions -f 1.0 -s -wa -wb | cut -f 10 | sort | uniq -c > $temp_result";
  print BLUE, "\n\t$bed_cmd3", RESET;
  system($bed_cmd3);
  open (COUNTS, "$temp_result") || die "\n\nCould not open temp results file: $temp_result\n\n";
  while(<COUNTS>){
    chomp($_);
    if ($_ =~ /(\d+)\s+(.*)/){
      my $count = $1;
      my $jid = $2;
      if ($observed_junctions{$jid}){
        my $anchored = $observed_junctions{$jid}{anchored};
        unless ($anchored eq "N"){
          $observed_junctions{$jid}{acceptors_skipped} = $count;
        }
      }else{
        print RED, "\n\nUnrecognized junction id: $jid\n\n", RESET;
        exit();
      }
    }else{
      print RED, "\n\nEntry in results file not understood: $_\n\n", RESET;
      exit();
    }
  }
  close(COUNTS);



  
  #Clean up the temp files
  my $rm_cmd = "rm -f $temp_obs_junctions $temp_known_exons $temp_known_donors $temp_known_acceptors $temp_result";
  print BLUE, "\n\nCleaning up ...\n$rm_cmd", RESET;
  system ($rm_cmd);


  return();
}


##############################################################################################################################################
# Join Splign results onto the junction annotations (splign validated, number of splign junctions)                                           #
##############################################################################################################################################
sub joinSplignResults{
  my %args = @_;
  my $infile = $args{'-infile'};
  
  print BLUE, "\n\nJoining pre-computed splign results to the annotation of each hmmSplicer junction\n", RESET;

  #Go through the splign results and add splign values to the observed junctions object
  open (SPLIGN, "$infile") || die "\n\nCould not open splign results file\n\n";
  while(<SPLIGN>){
    chomp($_);
    my @line = split("\t", $_);
    my $jid = $line[0];
    my $splign_validated = $line[1];
    my $splign_alignment_count = $line[2];
    
    if ($observed_junctions{$jid}){
      $observed_junctions{$jid}{splign_validated} = $splign_validated;
      $observed_junctions{$jid}{splign_alignment_count} = $splign_alignment_count;
    }
  }
  close(SPLIGN);

  #Make sure splign results were found for every junction to be annotated!
  foreach my $jid (keys %observed_junctions){
    unless (defined($observed_junctions{$jid}{splign_validated})){
      print RED, "\n\nFailed to find splign results for an observed junction, update splign results for this project!!\n\n", RESET;
      exit();
    }
  }
  return()
}

