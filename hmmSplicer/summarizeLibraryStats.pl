#!/usr/bin/perl -w
# Malachi Griffith

#Summarize the following for each library of a project:
#1.) Total unique junctions found
#2.) Total junction mapping reads
#3.) Total lanes contributing to the analysis
#4.) Read lengths of each lane (average read length)
#5.) Observed junctions corresponding to a known junction (matching donor AND acceptor) (same thing for junctions observed more than once)
#    - Junction class 'DA'
#6.) Junctions that can be anchored to a known junction (matching donor OR acceptor) - (same thing for junctions observed more than once)
#    - Junction class 'NDA', 'D', 'A'
#7.) Junctions that could not be anchored at all.
#    - Junction class 'N'
#8.) Percent of all known junctions detected by one or more reads
#9.) Percent of observed junctions with 2 or more supporting reads that correspond to a known junction
#10.) Exons skipped
#     - Percent of observed DA junctions corresponding to 0 exons skipped, 1 exons skipped, 1 or more exons skipped (known exon skipping)
#     - Percent of observed NDA junctions corresponding to 0 exons skipped, 1 exons skipped, 1 or more exons skipped (novel exon skipping)
#     - Percent of observed DA+NDA junctions corresponding to 0 exons skipped, 1 exons skipped, 1 or more exons skipped (known/novel exon skipping)

#Start with a library list for the project and collect stats for each library

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;
use Benchmark;

my $analysis_dir = '';
my $project = '';
my $splign_filter = '';

GetOptions ('analysis_dir=s'=>\$analysis_dir, 'project=s'=>\$project, 'splign_filter=s'=>\$splign_filter);

if ($analysis_dir && $project && $splign_filter){
  #Parameters ok
}else{
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  print GREEN, "\nSpecify the base analysis dir using: --analysis_dir", RESET;
  print GREEN, "\nSpecify the project name using: --project", RESET;
  print GREEN, "\nTo limit the analysis to Splign validated junctions only use: --splign_filter=yes", RESET;
  print GREEN, "\n\nExample: summarizeLibraryStats.pl  --analysis_dir=/projects/alexa2/hmmSplicer/  --project=SA_TN_Breast  --splign_filter=yes\n\n", RESET;
  exit();
}
chomp($splign_filter);
#Check the analysis dir
unless ($analysis_dir =~ /\/$/){
  $analysis_dir .= "/";
}
unless (-e $analysis_dir && -d $analysis_dir){
  print RED, "\n\nAnalysis dir: $analysis_dir not found\n\n", RESET;
  exit();
}
#Find and import the library list file for the project
my $library_list = "$analysis_dir$project/jobs/LibraryList.txt";
my %libraries = %{&importLibraries('-infile'=>$library_list)};

#Find and import the list of annotated junctions for the project
#This is the superset of junction observed for the entire project in at least one library
#Each junction has been annotated with respect to known transcript annotations from various transcript annotation resources
my $anno_junc_file = "$analysis_dir$project/results/$project".".junction.list.anno.txt";
my $anno_juncs_ref = &importObservedJunctions('-infile'=>$anno_junc_file);

#Now go through each individual library and summarize the junctions observed
my $outfile = "$analysis_dir$project/results/$project".".junction.library.summary.txt";
&summarizeLibraries('-outfile'=>$outfile);


print "\n\n";

exit();




##################################################################################################################################################
#Find and import list of library for the project                                                                                                 #
##################################################################################################################################################
sub importLibraries{
  my %args = @_;
  my $infile = $args{'-infile'};

  unless (-e $infile){
    print RED, "\n\nCould not find input library list file for this project:\n\t$infile\n\n", RESET;
    exit();
  }
  my %libs;

  print BLUE, "\n\nImporting list of libraries to summarize from: $infile", RESET;
  open(LIBS, "$infile") || die "\n\nCould not open input file: $infile\n\n";
  my $c = 0;
  while(<LIBS>){
    $c++;
    chomp($_);
    my @line = split("\t", $_);
    $libs{$line[0]}{count} = $c;
  }
  close(LIBS);

  return(\%libs);
}


##################################################################################################################################################
#Find and import the list of annotated junctions for the project                                                                                 #
##################################################################################################################################################
sub importObservedJunctions{
  my %args = @_;
  my $infile = $args{'-infile'};

  unless (-e $infile){
    print RED, "\n\nCould not find input annotated junctions file for this project:\n\t$infile\n\n", RESET;
    exit();
  }
  my %anno;
  
  print BLUE, "\n\nImporting observed junction annotations for the project from: $infile", RESET;
  open(ANNO, "$infile") || die "\n\nCould not open input file: $infile\n\n";
  my $header = 1;
  my %columns;
  while(<ANNO>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header == 1){
      my $pos = 0;
      foreach my $head (@line){
        $columns{$head}{pos} = $pos;
        $pos++;
      }
      $header = 0;
      next();
    }
    my $jid = $line[$columns{'JID'}{pos}];

    $anno{$jid}{read_count} = $line[$columns{'Read_Count'}{pos}];
    $anno{$jid}{score} = $line[$columns{'Score'}{pos}];
    $anno{$jid}{anchored} = $line[$columns{'Anchored'}{pos}];
    $anno{$jid}{exons_skipped} = $line[$columns{'Exons_Skipped'}{pos}];
    $anno{$jid}{donors_skipped} = $line[$columns{'Donors_Skipped'}{pos}];
    $anno{$jid}{acceptors_skipped} = $line[$columns{'Acceptors_Skipped'}{pos}];
    $anno{$jid}{gene_name} = $line[$columns{'Gene_Name'}{pos}];
    $anno{$jid}{splign_validated} = $line[$columns{'Splign_Validated'}{pos}];
    $anno{$jid}{splign_align_count} = $line[$columns{'Splign_Align_Count'}{pos}];
  }
  close(ANNO);

  return(\%anno);
}


##################################################################################################################################################
#Now go through each individual library and summarize the junctions observed
sub summarizeLibraries{
  my %args = @_;
  my $outfile = $args{'-outfile'};

  print BLUE, "\n\nSummarizing the junctions observed for each library...\n", RESET;
  my %lib_stats;
  my @cutoffs = qw (0 1 10);
  my $class_header_string = '';
  my $first_lib = 1;

  foreach my $lib (sort {$libraries{$a}{count} <=> $libraries{$b}{count}} keys %libraries){
    my $lib_count = $libraries{$lib}{count};
    $lib_stats{$lib_count}{lib_id} = $lib;
    my %lib_class;
    $lib_class{'DA'}{order} = 1;
    $lib_class{'NDA'}{order} = 2;
    $lib_class{'D'}{order} = 3;
    $lib_class{'A'}{order} = 4;
    $lib_class{'N'}{order} = 5;
    foreach my $cutoff (@cutoffs){
      foreach my $c (sort {$lib_class{$a}{order} <=> $lib_class{$b}{order}} keys %lib_class){
        $lib_class{$c}{$cutoff} = "na";
      }
    }

    my $lib_junction_file = "$analysis_dir$project/results/$lib/$lib".".junction.list.txt";
    print BLUE, "\nProcessing $lib ... ", RESET;
    my %lib_juncs;
    if (-e $lib_junction_file){
      print BLUE, "Found", RESET;
      foreach my $cutoff (@cutoffs){
        foreach my $c (sort {$lib_class{$a}{order} <=> $lib_class{$b}{order}} keys %lib_class){
          $lib_class{$c}{$cutoff} = 0;
          if ($first_lib){
            $class_header_string .= "$c"."_"."$cutoff\t"
          }
        }
      }
      if ($first_lib == 1){
        $first_lib = 0;
      }
      #Store the junctions for this library
      open (JUNC, "$lib_junction_file") || die "\n\nCould not open junction list file for $lib: $lib_junction_file\n\n";
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
        $lib_juncs{$jid}{read_count} = $line[$columns{'Read_Count'}{pos}];
        $lib_juncs{$jid}{score} = $line[$columns{'Score'}{pos}];
      }
      close(JUNC);

      #Now determine library level stats
      $lib_stats{$lib_count}{distinct_junctions_0r} = 0;
      $lib_stats{$lib_count}{distinct_junctions_1r} = 0;
      $lib_stats{$lib_count}{distinct_junctions_10r} = 0;
      $lib_stats{$lib_count}{junction_reads} = 0;
      
      $lib_stats{$lib_count}{es_known_0s} = 0;
      $lib_stats{$lib_count}{es_known_1s} = 0;
      $lib_stats{$lib_count}{es_known_1plus} = 0;
      $lib_stats{$lib_count}{es_novel_0s} = 0;
      $lib_stats{$lib_count}{es_novel_1s} = 0;
      $lib_stats{$lib_count}{es_novel_1plus} = 0;
  
      #Get counts for each junction class
      foreach my $jid (keys %lib_juncs){
        #Make sure an annotation record exists for every observed junction
        unless ($anno_juncs_ref->{$jid}){
          print RED, "\n\nA junction was found for this library that is not defined in the junction annotations file: re-run merge and annotation jobs before proceeding...\n\n", RESET;
          exit();
        }
        my $read_count = $lib_juncs{$jid}{read_count};
        my $score = $lib_juncs{$jid}{score};
        my $class = $anno_juncs_ref->{$jid}->{anchored};
        my $exons_skipped = $anno_juncs_ref->{$jid}->{exons_skipped};
        my $splign_validated = $anno_juncs_ref->{$jid}->{splign_validated};
        my $splign_align_count = $anno_juncs_ref->{$jid}->{splign_align_count};

        #Apply the splign filter if specified by the user
        if ($splign_filter =~ /^y$|^yes$/i){
          unless ($splign_validated == 1 && $splign_align_count <= 2){
            next();
          }
        }


        #Summarize exons skipped
        #Percent of observed DA junctions corresponding to 0 exons skipped, 1 exons skipped, 1 or more exons skipped (known exon skipping)
        #Percent of observed NDA junctions corresponding to 0 exons skipped, 1 exons skipped, 1 or more exons skipped (novel exon skipping)
        #Percent of observed DA+NDA junctions corresponding to 0 exons skipped, 1 exons skipped, 1 or more exons skipped (known/novel exon skipping)
        if ($class eq "DA"){
          if ($exons_skipped =~ /\d+/){
            if ($exons_skipped == 0){$lib_stats{$lib_count}{es_known_0s}++;}
            if ($exons_skipped == 1){$lib_stats{$lib_count}{es_known_1s}++;}
            if ($exons_skipped >= 1){$lib_stats{$lib_count}{es_known_1plus}++;}
          }
        }
        if ($class eq "NDA"){
          if ($exons_skipped =~ /\d+/){
            if ($exons_skipped == 0){$lib_stats{$lib_count}{es_novel_0s}++;}
            if ($exons_skipped == 1){$lib_stats{$lib_count}{es_novel_1s}++;}
            if ($exons_skipped >= 1){$lib_stats{$lib_count}{es_novel_1plus}++;}
          }
        }
        $lib_stats{$lib_count}{junction_reads} += $read_count;

        foreach my $cutoff (@cutoffs){
          if ($read_count > $cutoff){
            $lib_class{$class}{$cutoff}++;
          }
        }
        if ($read_count > 0){$lib_stats{$lib_count}{distinct_junctions_0r}++;}
        if ($read_count > 1){$lib_stats{$lib_count}{distinct_junctions_1r}++;}
        if ($read_count > 10){$lib_stats{$lib_count}{distinct_junctions_10r}++;}
      }
      $lib_stats{$lib_count}{lib_class} = \%lib_class;

    }else{
      print YELLOW, "Not found", RESET;
      #Set values to 'na' for this library
      $lib_stats{$lib_count}{junction_reads} = "na";
      $lib_stats{$lib_count}{lib_class} = \%lib_class;
      $lib_stats{$lib_count}{distinct_junctions_0r} = "na";
      $lib_stats{$lib_count}{distinct_junctions_1r} = "na";
      $lib_stats{$lib_count}{distinct_junctions_10r} = "na";
      $lib_stats{$lib_count}{es_known_0s} = "na";
      $lib_stats{$lib_count}{es_known_1s} = "na";
      $lib_stats{$lib_count}{es_known_1plus} = "na";
      $lib_stats{$lib_count}{es_novel_0s} = "na";
      $lib_stats{$lib_count}{es_novel_1s} = "na";
      $lib_stats{$lib_count}{es_novel_1plus} = "na";
    }
  }

  open(OUT, ">$outfile") || die "\n\nCould not open outfile for writing: $outfile\n\n";
  print OUT "Library_ID\tJunction_Reads\tDistinct_Junctions_0r\tDistinct_Junctions_1r\tDistinct_Junctions_10r\tES_Known_0s\tES_Known_1s\tES_Known_1plus\tES_Novel_0s\tES_Novel_1s\tES_Novel_1plus\t$class_header_string\n";
  foreach my $lc (sort {$a <=> $b} keys %lib_stats){
    my $lib_id = $lib_stats{$lc}{lib_id};
    my %lib_class = %{$lib_stats{$lc}{lib_class}};
    
    my $class_string = '';
    foreach my $cutoff (@cutoffs){
      foreach my $c (sort {$lib_class{$a}{order} <=> $lib_class{$b}{order}} keys %lib_class){
        $class_string .= "$lib_class{$c}{$cutoff}\t"
      }
    }
    print OUT "$lib_id\t$lib_stats{$lc}{junction_reads}\t$lib_stats{$lc}{distinct_junctions_0r}\t$lib_stats{$lc}{distinct_junctions_1r}\t$lib_stats{$lc}{distinct_junctions_10r}\t$lib_stats{$lc}{es_known_0s}\t$lib_stats{$lc}{es_known_1s}\t$lib_stats{$lc}{es_known_1plus}\t$lib_stats{$lc}{es_novel_0s}\t$lib_stats{$lc}{es_novel_1s}\t$lib_stats{$lc}{es_novel_1plus}\t$class_string\n";

  }
  close(OUT);

  print BLUE, "\n\nPrinted results to: $outfile", RESET;

  return();
}

















