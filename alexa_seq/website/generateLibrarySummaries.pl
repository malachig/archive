#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this website is to create a simple webpage for each sequencing library to summarize basic statistics about that library

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use website::WEB qw(:all);
use utilities::utility qw(:all);

my $paths_file = '';  #File containing root paths to basic types of data (expression, DE and SI)
my $search_page_url = '';
my $target_web_dir = '';
my $alexa_home_path = '';
my $alexa_seq_path = '';
my $google_analytics_id = '';
my $project_name = '';

GetOptions('paths_file=s'=>\$paths_file, 'search_page_url=s'=>\$search_page_url, 'target_web_dir=s'=>\$target_web_dir, 'alexa_home_path=s'=>\$alexa_home_path, 'alexa_seq_path=s'=>\$alexa_seq_path, 
           'google_analytics_id=s'=>\$google_analytics_id, 'project_name=s'=>\$project_name);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the project name using: --project_name", RESET;
print GREEN, "\n\tSpecify a file containing root paths to data using: --paths_file", RESET;
print GREEN, "\n\tSpecify the URL to your Xapian-Omega search page using:  --search_page_url", RESET;
print GREEN, "\n\tSpecify the target web-dir to place the files in using: --target_web_dir", RESET;
print GREEN, "\n\tSpecify the URL to your ALEXA home page using: --alexa_home_path", RESET;
print GREEN, "\n\tSpecify the URL to your ALEXA-Seq results page using: --alexa_seq_path", RESET;
print GREEN, "\n\tSpecify your Google Analytics ID using: --google_analytics_id", RESET;

print GREEN, "\n\nUsage: generateLibrarySummaries.pl  --project_name=5FU  --paths_file=/projects/malachig/solexa/batch_jobs/5FU/library_paths.txt  --search_page_url=http://www.bcgsc.ca/xapian-search/omega  --target_web_dir=/gsc/www/alexaplatform.org/alexa_seq/5FU/  --alexa_home_path=http://www.alexaplatform.org/index.htm  --alexa_seq_path=http://www.alexaplatform.org/alexa_seq/results.htm  --google_analytics_id=UA-xxxxxx-x\n\n", RESET;

#Check user supplied options
unless ($project_name && $paths_file && $search_page_url && $target_web_dir && $alexa_home_path && $alexa_seq_path && $google_analytics_id){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}

$target_web_dir = &checkDir('-dir'=>$target_web_dir, '-clear'=>"no");

#Define read assignment classes and assign names
my %read_classes;
$read_classes{'Reads'}{name} = "Total read count";
$read_classes{'Reads'}{order} = 1;
$read_classes{'Duplicate'}{name} = "Read1-Read2 identical";
$read_classes{'Duplicate'}{order} = 2;
$read_classes{'Low_Complexity'}{name} = "Low complexity";
$read_classes{'Low_Complexity'}{order} = 3;
$read_classes{'Low_Quality'}{name} = "Low quality (>1 N)";
$read_classes{'Low_Quality'}{order} = 4;
$read_classes{'ENST_U'}{name} = "Ensembl transcript";
$read_classes{'ENST_U'}{order} = 5;
$read_classes{'ENST_A'}{name} = "Ensembl transcript (ambiguous)";
$read_classes{'ENST_A'}{order} = 6;
$read_classes{'NOVEL_JUNCTION_U'}{name} = "Novel exon junction";
$read_classes{'NOVEL_JUNCTION_U'}{order} = 7;
$read_classes{'NOVEL_JUNCTION_A'}{name} = "Novel exon junction (ambiguous)";
$read_classes{'NOVEL_JUNCTION_A'}{order} = 8;
$read_classes{'NOVEL_BOUNDARY_U'}{name} = "Novel exon boundary extension";
$read_classes{'NOVEL_BOUNDARY_U'}{order} = 9;
$read_classes{'NOVEL_BOUNDARY_A'}{name} = "Novel exon boundary extension (ambiguous)";
$read_classes{'NOVEL_BOUNDARY_A'}{order} = 10;
$read_classes{'INTRON_U'}{name} = "Intron";
$read_classes{'INTRON_U'}{order} = 11;
$read_classes{'INTRON_A'}{name} = "Intron (ambiguous)";
$read_classes{'INTRON_A'}{order} = 12;
$read_classes{'INTERGENIC_U'}{name} = "Intergenic";
$read_classes{'INTERGENIC_U'}{order} = 13;
$read_classes{'INTERGENIC_A'}{name} = "Intergenic (ambiguous)";
$read_classes{'INTERGENIC_A'}{order} = 14;
$read_classes{'Repeat_U'}{name} = "Repeat element";
$read_classes{'Repeat_U'}{order} = 15;
$read_classes{'Repeat_A'}{name} = "Repeat element (ambiguous)";
$read_classes{'Repeat_A'}{order} = 16;
$read_classes{'Unassigned'}{name} = "Unassigned";
$read_classes{'Unassigned'}{order} = 17;

#Get paths variables from Lib_Paths file
my $data_paths_ref = &getDataPaths('-file'=>$paths_file);
my $exp_ref = $data_paths_ref->{'Expression'};
my $exp_paths_ref = $exp_ref->{'paths'};
my $qual_ref = $data_paths_ref->{'LibraryQuality'};
my $qual_paths_ref = $qual_ref->{'paths'};
my $mapping_ref = $data_paths_ref->{'TranscriptMapping'};
my $mapping_paths_ref = $mapping_ref->{'paths'};
my $commands_ref = $data_paths_ref->{'AnalysisCommands'};
my $commands_paths_ref = $commands_ref->{'paths'};

#For each library, generate summary page
foreach my $lib (sort {$exp_paths_ref->{$a}->{line_order} <=> $exp_paths_ref->{$b}->{line_order}} keys %{$exp_paths_ref}){
  my $content_count = 0;
  my $box_id = '';
  my $images_dir = "$target_web_dir"."images/$lib/";
  unless (-e $images_dir){
    mkdir($images_dir);
  }

  my $lib_name = $exp_paths_ref->{$lib}->{name};
  my $exp_data_path = $exp_paths_ref->{$lib}->{data_path};
  my $exp_stats_path = $exp_paths_ref->{$lib}->{stats_path};
  my $qual_data_path = $qual_paths_ref->{$lib}->{data_path};
  my $qual_stats_path = $qual_paths_ref->{$lib}->{stats_path};
  my $mapping_data_path = $mapping_paths_ref->{$lib}->{data_path};
  my $mapping_stats_path = $mapping_paths_ref->{$lib}->{stats_path};
  my $commands_data_path = $commands_paths_ref->{$lib}->{data_path};
  my $commands_stats_path = $commands_paths_ref->{$lib}->{stats_path};
  
  print BLUE, "\n\n$lib ($lib_name) ($exp_stats_path) ($qual_stats_path) ($mapping_stats_path)", RESET;

  #N.) Get and summarize basic library info lane-by-lane (passing read count, read length, source dir)
  my $lane_info_file = "$commands_data_path"."$project_name"."_Lib_Data.txt";
  $content_count++;
  $box_id = "box"."$content_count"."s";
  my $li_content = "\n<!-- Table summarizing input lane data -->";
  $li_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Summary of input lanes for this library</P>";
  $li_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $li_content .= "\n<div id=\"$box_id\">";
  $li_content .= "\n<div id=\"li_legend\" style=\"width: 700px;\">";
  $li_content .= "\n<P CLASS=\"Indented24LR_s16\">Basic info is provided for each lane of data comprising this library</P><BR>";
  $li_content .= "\n</div>";

  if (-e $lane_info_file){
    my %li;
    open (LI, "$lane_info_file") || die "\nCould not open read assignment stats file: $lane_info_file\n\n";
    my $lc = 0;
    while(<LI>){
      chomp($_);
      my @line = split("\t", $_);
      my $library_id = $line[0];
      unless ($library_id eq $lib){
        next();
      }
      $lc++;
      $li{$lc}{flowcell_lane} = $line[2];
      $li{$lc}{read_length} = $line[5];
      $li{$lc}{quality_reads} = &commify($line[7]);
      $li{$lc}{source_dir} = $line[8];

      #Shorten the source dir to obscure complete path for security reasons
      if ($li{$lc}{source_dir} =~ /.*\/alexa\_seq\/(.*)/){
        $li{$lc}{source_dir} = $1;
      }

    }
    close(LI);
    $li_content .= "<TABLE CLASS=\"Data2\">\n\t<TR><TD CLASS=\"Head1\">Flowcell-Lane</TD><TD CLASS = \"Head3\">Read Length</TD><TD CLASS=\"Head3\">Quality Reads</TD><TD CLASS=\"Head3\">Source Dir</TD></TR>";
    foreach my $lc (sort {$a <=> $b} keys %li){
      $li_content .= "\n\t\t<TR><TD CLASS=\"Data3\">$li{$lc}{flowcell_lane}</TD><TD CLASS=\"Data2\">2 x $li{$lc}{read_length}</TD><TD CLASS=\"Data2\">$li{$lc}{quality_reads}</TD><TD CLASS=\"Data2\">$li{$lc}{source_dir}</TD></TR>";
    }
    $li_content .= "</TABLE>";
  }else{
    print YELLOW, "\n\n($content_count) Could not find lane info file: $lane_info_file\n\n", RESET;
    $li_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $li_content .= "\n<BR></div><BR>";


  #N.) Get and summarize basic library stats lane-by-lane
  my $lane_stats_file = "$exp_data_path"."Summary/SummaryStats.csv";
  $content_count++;
  $box_id = "box"."$content_count"."s";
  my $lane_stats_content = "\n<!-- Table summarizing sequence statistics for each lane -->";
  $lane_stats_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Summary of sequence statistics for each lane</P>";
  $lane_stats_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $lane_stats_content .= "\n<div id=\"$box_id\">";
  $lane_stats_content .= "\n<div id=\"lane_stats_legend\" style=\"width: 1000px;\">";
  $lane_stats_content .= "\n<P CLASS=\"Indented24LR_s16\">Basic statistics for each lane of sequence data are summarized below. 'Low quality' refers to reads with too many  ambiguous bases (N's). 'Low complexity' refers to homopolymeric
                                                          (e.g. polyA tails) or other sequences of low complexity identified by 'mdust'.  'R1' and 'R2' refers to the first and second read of a paired-end read.</P><BR>";
  $lane_stats_content .= "\n</div>";

  if (-e $lane_stats_file){
    my %ls;
    open(LS, "$lane_stats_file") || die "\n\nCould not open lane stats file: $lane_stats_file\n\n";
    my $header = 1;
    my %columns;
    my $col_count = 0;
    my $lc = 0;
    while(<LS>){
      chomp($_);
      my @line = split(",", $_);
      if ($header == 1){
        $header = 0;
        foreach my $value (@line){
          $columns{$value}{position} = $col_count;
          $col_count++;
        }
        next();
      }
      $lc++;
      my $total_readpairs = $line[$columns{'total_readpairs'}{position}];
      my $total_reads = $line[$columns{'total_reads'}{position}];
      my $duplicate_reads = ($line[$columns{'duplicate_readpair_reads'}{position}] + $line[$columns{'duplicate_readpair_reads_afterRC'}{position}]);
      my $r1_low_quality_reads = $line[$columns{'R1_low_quality_reads'}{position}];
      my $r1_low_complexity_reads = $line[$columns{'R1_low_complexity_reads'}{position}];
      my $r2_low_quality_reads = $line[$columns{'R2_low_quality_reads'}{position}];
      my $r2_low_complexity_reads = $line[$columns{'R2_low_complexity_reads'}{position}];
      my $r1_passing_reads = $total_readpairs-($duplicate_reads+$r1_low_quality_reads+$r1_low_complexity_reads);
      my $r2_passing_reads = $total_readpairs-($duplicate_reads+$r2_low_quality_reads+$r2_low_complexity_reads);
      my $r12_passing_reads = $r1_passing_reads+$r2_passing_reads;

      my $r12_total_bases = $line[$columns{'total_bases'}{position}];
      my $r1_total_bases = $r12_total_bases/2;
      my $r2_total_bases = $r12_total_bases/2;
      my $r1_ambiguous_bases = $line[$columns{'R1_ambiguous_bases'}{position}];
      my $r2_ambiguous_bases = $line[$columns{'R2_ambiguous_bases'}{position}];
      my $r12_ambiguous_bases = $line[$columns{'R12_ambiguous_bases'}{position}];
      my $r1_mdust_bases = $line[$columns{'R1_mdust_bases'}{position}];
      my $r2_mdust_bases = $line[$columns{'R2_mdust_bases'}{position}];
      my $r12_passing_bases = $r12_total_bases-($r1_ambiguous_bases+$r2_ambiguous_bases+$r1_mdust_bases+$r2_mdust_bases);

      $ls{$lc}{flowcell_name} = $line[$columns{'flowcell_name'}{position}];
      $ls{$lc}{lane_number} = $line[$columns{'lane_number'}{position}];
      $ls{$lc}{total_readpairs} = &commify($total_readpairs);
      $ls{$lc}{R1_low_quality_reads_p} = sprintf("%.2f", (($r1_low_quality_reads/$total_readpairs)*100));
      $ls{$lc}{R2_low_quality_reads_p} = sprintf("%.2f", (($r2_low_quality_reads/$total_readpairs)*100));
      $ls{$lc}{R1_low_complexity_reads_p} = sprintf("%.2f", (($r1_low_complexity_reads/$total_readpairs)*100));
      $ls{$lc}{R2_low_complexity_reads_p} = sprintf("%.2f", (($r2_low_complexity_reads/$total_readpairs)*100));
      $ls{$lc}{R1_passing_reads_p} = sprintf("%.2f", (($r1_passing_reads/$total_readpairs)*100));
      $ls{$lc}{R2_passing_reads_p} = sprintf("%.2f", (($r2_passing_reads/$total_readpairs)*100));
      $ls{$lc}{R12_total_bases} = &commify($r12_total_bases);
      $ls{$lc}{R1_ambiguous_bases_p} = sprintf("%.2f", (($r1_ambiguous_bases/$r1_total_bases)*100));
      $ls{$lc}{R2_ambiguous_bases_p} = sprintf("%.2f", (($r2_ambiguous_bases/$r2_total_bases)*100));
      $ls{$lc}{R12_ambiguous_bases_p} = sprintf("%.2f", (($r12_ambiguous_bases/$r12_total_bases)*100));
      $ls{$lc}{R12_passing_bases_p} =  sprintf("%.2f", (($r12_passing_bases/$r12_total_bases)*100));
    }
    close(LS);

    #Build data table
    $lane_stats_content .= "<TABLE CLASS=\"Data2\">\n\t<TR><TD CLASS=\"Head1\">Flowcell-Lane</TD><TD CLASS=\"Head3\">Paired Reads</TD><TD CLASS=\"Head3\">Low Quality<br>(R1 | R2)</TD><TD CLASS=\"Head3\">Low Complexity<br>(R1 | R2)</TD><TD CLASS=\"Head3\">Passing<br>(R1 | R2)</TD><TD CLASS=\"Head3\">Total Bases</TD><TD CLASS=\"Head3\">Passing Bases</TD></TR>";
    foreach my $lc (sort {$a <=> $b} keys %ls){
      my $flowcell_lane = "$ls{$lc}{flowcell_name}"."_Lane"."$ls{$lc}{lane_number}";
      $lane_stats_content .= "\n\t\t<TR><TD CLASS=\"Data3\">$flowcell_lane</TD><TD CLASS=\"Data2\">$ls{$lc}{total_readpairs}</TD><TD CLASS=\"Data2\">$ls{$lc}{R1_low_quality_reads_p}% | $ls{$lc}{R2_low_quality_reads_p}%</TD><TD CLASS=\"Data2\">$ls{$lc}{R1_low_complexity_reads_p}% | $ls{$lc}{R2_low_complexity_reads_p}%</TD><TD CLASS=\"Data2\">$ls{$lc}{R1_passing_reads_p}% | $ls{$lc}{R2_passing_reads_p}%</TD><TD CLASS=\"Data2\">$ls{$lc}{R12_total_bases}</TD><TD CLASS=\"Data2\">$ls{$lc}{R12_passing_bases_p}%</TD></TR>";
    }
    $lane_stats_content .= "</TABLE>";

  }else{
    print YELLOW, "\n\n($content_count) Could not find lane stats file: $lane_stats_file\n\n", RESET;
    $lane_stats_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $lane_stats_content .= "\n<BR></div><BR>";


  #N.) Get and summarize mapping summary stats lane-by-lane
  my $mapped_reads_file = "$exp_data_path"."Summary/MappedReads.txt";
  $content_count++;
  $box_id = "box"."$content_count"."s";
  
  my $mr_content = "\n<!-- Table summarizing mapping summary stats lane-by-lane -->";
  $mr_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Summary of mapping results for each lane</P>";
  $mr_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $mr_content .= "\n<div id=\"$box_id\">";
  $mr_content .= "\n<div id=\"mr_legend\" style=\"width: 1200px;\">";
  $mr_content .= "\n<P CLASS=\"Indented24LR_s16\">The mapping of reads to genome and transcriptome sequences is summarized below on a lane-by-lane basis.  Only reads mapped with high confidence are summarized here.
                                                  Reads mapped to each sequence type were assigned unambiguously.  Reads mapping ambiguously (i.e. map equally well to multiple places) are summarized in the last column.</P><BR>";
  $mr_content .= "\n</div>";


  if (-e $mapped_reads_file){
    my %mr;
    open(MR, "$mapped_reads_file") || die "\n\nCould not open mapped reads file: $mapped_reads_file\n\n";
    my $header = 1;
    my %columns;
    my $col_count = 0;
    my $mc = 0;
    while(<MR>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        $header = 0;
        foreach my $value (@line){
          $columns{$value}{position} = $col_count;
          $col_count++;
        }
        next();
      }
      $mc++;
      my $lane = $line[$columns{'LANE'}{position}];
      my $enst_reads = $line[$columns{'ENST_Reads'}{position}];
      my $junction_reads = $line[$columns{'NOVEL_JUNCTION_Reads'}{position}];
      my $boundary_reads = $line[$columns{'NOVEL_BOUNDARY_Reads'}{position}];
      my $intron_reads = $line[$columns{'INTRON_Reads'}{position}];
      my $intergenic_reads = $line[$columns{'INTERGENIC_Reads'}{position}];
      my $repeat_reads = $line[$columns{'REPEAT_Reads'}{position}];
      my $ambiguous_reads = $line[$columns{'AMBIGUOUS_Reads'}{position}];
      my $total_reads = $line[$columns{'TOTAL_Reads'}{position}];
      $mr{$mc}{flowcell_lane} = $lane;
      $mr{$mc}{total_reads} = &commify($total_reads);
      $mr{$mc}{enst_reads_p} = sprintf("%.2f", (($enst_reads/$total_reads)*100));
      $mr{$mc}{junction_reads_p} = sprintf("%.2f", (($junction_reads/$total_reads)*100));
      $mr{$mc}{boundary_reads_p} = sprintf("%.2f", (($boundary_reads/$total_reads)*100));
      $mr{$mc}{intron_reads_p} = sprintf("%.2f", (($intron_reads/$total_reads)*100));
      $mr{$mc}{intergenic_reads_p} = sprintf("%.2f", (($intergenic_reads/$total_reads)*100));
      $mr{$mc}{repeat_reads_p} = sprintf("%.2f", (($repeat_reads/$total_reads)*100));
      $mr{$mc}{ambiguous_reads_p} = sprintf("%.2f", (($ambiguous_reads/$total_reads)*100));
    }
    close(MR);

    $mr_content .= "<TABLE CLASS=\"Data2\">\n\t<TR><TD CLASS=\"Head1\">Flowcell-Lane</TD><TD CLASS = \"Head3\">Total</TD><TD CLASS=\"Head3\">Transcript</TD><TD CLASS=\"Head3\">Novel junction</TD><TD CLASS=\"Head3\">Novel boundary</TD><TD CLASS=\"Head3\">Intron</TD><TD CLASS=\"Head3\">Intergenic</TD><TD CLASS=\"Head3\">Repeat element</TD><TD CLASS=\"Head3\">Ambiguous</TD></TR>";
    foreach my $mc (sort {$a <=> $b} keys %mr){
      $mr_content .= "\n\t\t<TR><TD CLASS=\"Data3\">$mr{$mc}{flowcell_lane}</TD><TD CLASS=\"Data2\">$mr{$mc}{total_reads}</TD><TD CLASS=\"Data2\">$mr{$mc}{enst_reads_p}%</TD><TD CLASS=\"Data2\">$mr{$mc}{junction_reads_p}%</TD><TD CLASS=\"Data2\">$mr{$mc}{boundary_reads_p}%</TD><TD CLASS=\"Data2\">$mr{$mc}{intron_reads_p}%</TD><TD CLASS=\"Data2\">$mr{$mc}{intergenic_reads_p}%</TD><TD CLASS=\"Data2\">$mr{$mc}{repeat_reads_p}%</TD><TD CLASS=\"Data2\">$mr{$mc}{ambiguous_reads_p}%</TD></TR>";
    }
    $mr_content .= "</TABLE>";

  }else{
    print YELLOW, "\n\n($content_count) Could not find read assign file: $mapped_reads_file\n\n", RESET;
    $mr_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $mr_content .= "\n<BR></div><BR>";


  #1.) Get the read assignment stats - calculate percents - create display content
  my $read_assign_file = "$exp_stats_path"."ReadAssignmentStats.txt";
  $content_count++;
  $box_id = "box"."$content_count"."s";
  my $ras_content = "\n<!-- Table summarizing number of reads assigned to each class -->";
  $ras_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Summary of read assignments by assignment class</P>";
  $ras_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $ras_content .= "\n<div id=\"$box_id\">";
  $ras_content .= "\n<div id=\"mr_legend\" style=\"width: 700px;\">";
  $ras_content .= "\n<P CLASS=\"Indented24LR_s16\">The assignment of reads to each read class is summarized below for the entire library.  'Unassigned' reads could not be assigned with high confidence to any known or predicted 
                                                  transcriptome or genome sequence.  This does not mean that they had no significant similarity, only that they could not be assigned with high confidence.  In this table, reads
                                                  are summarized as individual reads (not paired reads)</P><BR>";
  $ras_content .= "\n</div>";

  if (-e $read_assign_file){
    my %ras;
    open (RAS, "$read_assign_file") || die "\nCould not open read assignment stats file: $read_assign_file\n\n";
    while(<RAS>){
      chomp($_);
      if ($_ =~ /(\S+)\:\s+(\d+)/){
        if ($read_classes{$1}){
          $ras{$1}{count} = $2;
        }
      }
    }
    close(RAS);
    foreach my $class (keys %ras){
      my $percent = ($ras{$class}{count}/$ras{'Reads'}{count})*100;
      $ras{$class}{percent} = sprintf("%.2f", $percent);
      $ras{$class}{count_pretty} = &commify($ras{$class}{count});
      $ras{$class}{order} = $read_classes{$class}{order};
      $ras{$class}{name} = $read_classes{$class}{name};
    }
    $ras{'Reads'}{percent} = 100;
    $ras_content .= "<TABLE CLASS=\"Data2\">\n\t<TR><TD CLASS=\"Head1\">Read Class</TD><TD CLASS = \"Head3\">Read Count</TD><TD CLASS=\"Head3\">Percent of Total</TD></TR>";
    foreach my $class (sort {$ras{$a}->{order} <=> $ras{$b}->{order}} keys %ras){
      $ras_content .= "\n\t\t<TR><TD CLASS=\"Data3\">$ras{$class}{name}</TD><TD CLASS=\"Data2\">$ras{$class}{count_pretty}</TD><TD CLASS=\"Data2\">$ras{$class}{percent}%</TD></TD></TR>";
    }
    $ras_content .= "</TABLE>";
  }else{
    print YELLOW, "\n\n($content_count) Could not find read assign file: $read_assign_file\n\n", RESET;
    $ras_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $ras_content .= "\n<BR></div><BR>";


  #2.) Get the grand average coverage values for each feature type
  my $coverage_file = "$exp_stats_path"."GrandAverageCoverageValues.txt";
  $content_count++;
  $box_id = "box"."$content_count"."s";
  my $gac_content = "\n<!-- Table summarizing grand average coverage statistics -->";
  $gac_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Summary of average coverage values by feature type</P>";
  $gac_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $gac_content .= "\n<div id=\"$box_id\">";
  $gac_content .= "\n<div id=\"gac_legend\" style=\"width: 800px;\">";
  $gac_content .= "\n<P CLASS=\"Indented24LR_s16\">The grand average coverage observed for each type of sequence feature is summarized below. Average coverage is calculated as the cumulative number of mapped bases divided by the total
                                                   number of base positions in the genome (for that type of sequence feature).</P><BR>";
  $gac_content .= "\n</div>";

  my %gac;
  if (-e $coverage_file){
    open (GAC, "$coverage_file") || die "\nCould not open grand average coverage file: $coverage_file\n\n";
    my $current_type = '';
    my $order = 0;
    while(<GAC>){
      chomp($_);
      if ($_ =~ /FOR\:\s+(\w+)/){
        $current_type = $1;
        $order++;
      }
      if ($_ =~ /.*\=\s+(\d+).*\=\s+(\d+).*\=\s+(\S+)/){
        $gac{$current_type}{order} = $order;
        $gac{$current_type}{cumulative_coverage}=&commify($1);
        $gac{$current_type}{base_count}=&commify($2);
        $gac{$current_type}{average_coverage}=$3;
      }
    }
    close(GAC);
    $gac_content .= "<TABLE CLASS=\"Data2\">\n\t<TR><TD CLASS=\"Head1\">Feature Type</TD><TD CLASS = \"Head3\">Average Coverage</TD><TD CLASS=\"Head3\">Total Base Count</TD><TD CLASS = \"Head3\">Cumulative Coverage</TD></TR>";
    foreach my $type (sort {$gac{$a}->{order} <=> $gac{$b}->{order}} keys %gac){
      $gac_content .= "\n\t\t<TR><TD CLASS=\"Data3\">$type</TD><TD CLASS=\"Data4\">$gac{$type}{average_coverage}</TD><TD CLASS=\"Data2\">$gac{$type}{base_count}</TD><TD CLASS=\"Data2\">$gac{$type}{cumulative_coverage}</TD></TR>";
    }
    $gac_content .= "</TABLE>";
  }else{
    print YELLOW, "\n\n($content_count) Could not find average coverage file: $coverage_file\n\n", RESET;
    $gac_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $gac_content .= "\n<BR></div><BR>";



  #3.) Summarize number of expressed events of each feature type in each library
  my $expressed_file = "$exp_stats_path"."$lib"."_ExpressedFeatures.txt";
  $content_count++;
  $box_id = "box"."$content_count"."s";
  my $ef_content = "\n<!-- Table summarizing expressed events by feature type -->";
  $ef_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Summary of expressed events by feature type</P>";
  $ef_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $ef_content .= "\n<div id=\"$box_id\">";
  $ef_content .= "\n<div id=\"ef_legend\" style=\"width: 700px;\">";
  $ef_content .= "\n<P CLASS=\"Indented24LR_s16\">The number of sequence features of each type detected above the level of background noise are summarized below.  For each feature type, the total number of features is listed, 
                                                  followed by the subset of these features that are expressed versus not expressed above background.</P><BR>";
  $ef_content .= "\n</div>";

  if (-e $expressed_file){
    my %ef;
    open(EF, "$expressed_file") || die "\nCould not open expressed features file: $expressed_file";
    my $header = 1;
    my $order = 0;
    $ef{100}{feature_type} = "GRAND TOTAL<BR>(non-redundant)";
    $ef{100}{total_events} = 0;
    $ef{100}{expressed} = 0;
    $ef{100}{not_expressed} = 0;

    while(<EF>){
      chomp($_);
      if($header){
        $header = 0;
        next();
      }
      my @data = split("\t", $_);
      $order++;
      $ef{$order}{feature_type} = $data[0];
      $ef{$order}{total_events} = &commify($data[1]); unless ($data[0] =~ /Known|Novel/){$ef{100}{total_events} += $data[1]};
      $ef{$order}{expressed} = &commify($data[2]); unless ($data[0] =~ /Known|Novel/){$ef{100}{expressed} += $data[2]};
      $ef{$order}{expressed_p} = $data[3];
      $ef{$order}{not_expressed} = &commify($data[4]); unless ($data[0] =~ /Known|Novel/){$ef{100}{not_expressed} += $data[4]};
      $ef{$order}{not_expressed_p} = $data[5];
    }
    $ef{100}{expressed_p} = sprintf("%.2f", (($ef{100}{expressed}/$ef{100}{total_events})*100)); 
    $ef{100}{not_expressed_p} = sprintf("%.2f", (($ef{100}{not_expressed}/$ef{100}{total_events})*100)); 
    $ef{100}{total_events} = &commify($ef{100}{total_events}); 
    $ef{100}{expressed} = &commify($ef{100}{expressed}); 
    $ef{100}{not_expressed} = &commify($ef{100}{not_expressed}); 
    close(EF);
    $ef_content .= "<TABLE CLASS=\"Data2\">\n\t<TR><TD CLASS=\"Head1\">Feature Type</TD><TD CLASS = \"Head3\">Feature Count</TD><TD CLASS=\"Head3\">Expressed (%)</TD><TD CLASS = \"Head3\">Not Expressed (%)</TD></TR>";

    foreach my $o (sort {$a <=> $b} keys %ef){
      $ef_content .= "\n\t\t<TR><TD CLASS=\"Data3\">$ef{$o}{feature_type}</TD><TD CLASS=\"Data2\">$ef{$o}{total_events}</TD><TD CLASS=\"Data4\">$ef{$o}{expressed} ($ef{$o}{expressed_p}%)</TD><TD CLASS=\"Data2\">$ef{$o}{not_expressed} ($ef{$o}{not_expressed_p}%)</TD></TR>";
    }
    $ef_content .= "</TABLE>";
  }else{
    print YELLOW, "\n\n($content_count) Could not find expressed features file: $expressed_file\n\n", RESET;
    $ef_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $ef_content .= "\n<BR></div><BR>";







  #3.) Summarize library expression skew
  my $skew_file = "$exp_stats_path"."LibraryExpressionSkew.txt";
  $content_count++;
  $box_id = "box"."$content_count"."s";
  my $sk_content = "\n<!-- Table summarizing library expression skew -->";
  $sk_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Summary of library expression skew (proportion of reads consumed by the top N% of expressed genes)</P>";
  $sk_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $sk_content .= "\n<div id=\"$box_id\">";
  $sk_content .= "\n<div id=\"ef_legend\" style=\"width: 700px;\">";
  $sk_content .= "\n<P CLASS=\"Indented24LR_s16\">The expression level of the most highly and lowly expressed can vary by several orders of magnitude.  In some RNA-seq libraries the most highly expressed genes can heavily dominate. 
                                                  This may genuinely reflect the nature of the RNA sample or it may be indicative of a problem with the library (excessive rRNA contamination, low complexity, etc.).
                                                  The following table summarizes the proportion of all reads mapped to genes consumed by the top N% sample genes (by raw read count).</P><BR>";
  $sk_content .= "\n</div>";

  if (-e $skew_file){
    my %sk;
    open(SK, "$skew_file") || die "\nCould not open expression skew file: $skew_file";
    my $header = 1;
    my $order = 0;

    while(<SK>){
      chomp($_);
      if($header){
        $header = 0;
        next();
      }
      my @data = split("\t", $_);
      $order++;
      $sk{$order}{percent_genes} = $data[0];
      $sk{$order}{number_genes} = &commify($data[1]); 
      $sk{$order}{read_count} = &commify($data[2]); 
      $sk{$order}{percent_reads} = $data[3];
      $sk{$order}{non_zero_genes} = &commify($data[4]);
    }
    close(SK);
    $sk_content .= "<TABLE CLASS=\"Data2\">\n\t<TR><TD CLASS=\"Head1\">Percent Genes</TD><TD CLASS = \"Head3\">Gene Count</TD><TD CLASS=\"Head3\">Read Count</TD><TD CLASS = \"Head3\">Percent Reads</TD><TD CLASS = \"Head3\">Non Zero Genes</TD></TR>";

    foreach my $o (sort {$a <=> $b} keys %sk){
      $sk_content .= "\n\t\t<TR><TD CLASS=\"Data3\">$sk{$o}{percent_genes}%</TD><TD CLASS=\"Data2\">$sk{$o}{number_genes}</TD><TD CLASS=\"Data2\">$sk{$o}{read_count}</TD> <TD CLASS=\"Data2\">$sk{$o}{percent_reads}%</TD><TD CLASS=\"Data2\">$sk{$o}{non_zero_genes}</TD></TR>";
    }
    $sk_content .= "</TABLE>";
  }else{
    print YELLOW, "\n\n($content_count) Could not find library expression skew file: $skew_file\n\n", RESET;
    $sk_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $sk_content .= "\n<BR></div><BR>";






  #4.) Calculate the (exon/silent_intron_region) and (exon/silent_intergenic_region) ratios using average coverage values for these feature types
  $content_count++;
  $box_id = "box"."$content_count"."s";
  my $signal_to_noise_content = "\n<!-- Section summarizing estimates of signal-to-noise ratios -->";
  $signal_to_noise_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Estimates of signal-to-noise ratio</P>";
  $signal_to_noise_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $signal_to_noise_content .= "\n<div id=\"$box_id\">";

  if (-e $coverage_file){
    my $er_to_si = sprintf("%.2f", ($gac{'ExonRegion'}{average_coverage}/$gac{'SilentIntronRegion'}{average_coverage}));
    my $er_to_sig = sprintf("%.2f", ($gac{'ExonRegion'}{average_coverage}/$gac{'SilentIntergenicRegion'}{average_coverage}));
    $signal_to_noise_content .= "\n<P CLASS=\"Indented24LR_s16\">(Average coverage of <i>exon regions</i> / Average coverage of <i>silent intron regions</i>) = $er_to_si</P>";
    $signal_to_noise_content .= "\n<P CLASS=\"Indented24LR_s16\">(Average coverage of <i>exon regions</i> / Average coverage of <i>silent intergenic regions</i>) = $er_to_sig</P>";
  }else{
    print YELLOW, "\n\n($content_count) Could not find signal-to-noise file: $coverage_file\n\n", RESET;
    $signal_to_noise_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $signal_to_noise_content .= "\n<BR></div><BR>";


  #5.) Get the 95th percentile value statements for both silent_intronic_region and silent_intergenic_regions
  my $exp_stats_file = "$exp_stats_path"."$lib"."_stats.txt";
  $content_count++;
  $box_id = "box"."$content_count"."s";
  my $percentiles_content = "\n<!-- Section summarizing 95th percentiles of silent intron and silent intergenic regions -->";
  $percentiles_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Estimates of intronic and intergenic noise levels (95th percentiles of silent intron and intergenic regions)</P>";
  $percentiles_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $percentiles_content .= "\n<div id=\"$box_id\">";
  
  if (-e $exp_stats_file){
    open (STATS, "$exp_stats_file") || die "\nCould not open stats file: $exp_stats_file\n\n";
    my $si_95th = '';
    my $sig_95th = '';

    while(<STATS>){
      chomp($_);
      if ($_ =~ /95th\s+percentile\s+of\s+silent\s+intron/){
        $si_95th = $_;
        $percentiles_content .= "<P CLASS=\"Indented24LR_s16\">$si_95th</P>";
      }
      if ($_ =~ /95th\s+percentile\s+of\s+silent\s+intergenic/){
        $sig_95th = $_;
        $percentiles_content .= "<P CLASS=\"Indented24LR_s16\">$sig_95th</P><BR>";
      }
    }
    close(STATS);
  }else{
    print YELLOW, "\n\n($content_count) Could not find percentiles file: $exp_stats_file\n\n", RESET;
    $percentiles_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $percentiles_content .= "\n<BR></div><BR>"; 


  #Get SVG/JPEG files for each figure to be displayed and copy them to the target images dir
  my $image_desc;
  my $image_width;
  my $image_height;

  #Figure content - Position bias
  $content_count++;
  $box_id = "box"."$content_count"."h";
  my $position_bias_content = "\n<!-- Figure - Position Bias -->";
  $image_desc = "Distribution of relative positions of reads mapping within known transcripts (i.e. position bias test)";
  $position_bias_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">$image_desc</P>";
  $position_bias_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $position_bias_content .= "\n<div id=\"$box_id\">";

  my $position_bias_file = "$mapping_stats_path"."PositionBias_LinesScaledRunMed.svg"; #750x750
  my $position_bias_file2 = "images/$lib/"."PositionBias_LinesScaledRunMed.svg";
  $image_width = 750;
  $image_height = 750;
  if (-e $position_bias_file){
    system("cp -f $position_bias_file $images_dir");
    $position_bias_content .= "\n<div id=\"position_bias_legend\" style=\"width: 750px;\">";
    $position_bias_content .= "\n<P CLASS=\"Indented24LR_s16\">The frequency of read positions are plotted against the relative position within each transcript, where 0% is the 5' end of the transcript and 100% is the 3' end.
                                                               Read positions were binned according to the size of transcripts they map within and the plot is produced for each bin (i.e. reads mapping to transcripts of 0-500 bp, 
                                                               500-1000 bp, ..., 15000-20000 bp, and > 20000 bp).</P><BR>";
    $position_bias_content .= "\n</div>";
    $position_bias_content .= "\n<P CLASS=\"Indented24LR\"><OBJECT DATA=\"$position_bias_file2\" TYPE=\"image/svg+xml\" WIDTH=\"$image_width\" HEIGHT=\"$image_height\" TITLE=\"$image_desc\" CLASS=\"Pic_bordered\" STANDBY=\"$image_desc\"><EMBED SRC=\"$position_bias_file2\" TYPE=\"image/svg+xml\" WIDTH=\"$image_width\" HEIGHT=\"$image_height\" /></OBJECT></P><BR>";

  }else{
    print YELLOW, "\n\n($content_count) Could not find figure file: $position_bias_file\n\n", RESET;
    $position_bias_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $position_bias_content .= "\n<BR></div><BR>";


  #Figure content - Fragment size
  $content_count++;
  $box_id = "box"."$content_count"."h";
  my $frag_size_content = "\n<!-- Figure - Fragment Size -->";
  $frag_size_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Distribution of fragment size for paired-end mappings to transcripts</P>";
  $frag_size_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $frag_size_content .= "\n<div id=\"$box_id\">";

  my $frag_size_file = "$mapping_stats_path"."FragmentSize.jpeg"; #750x750
  my $frag_size_file2 = "images/$lib/"."FragmentSize.jpeg";
  if (-e $frag_size_file){
    system("cp -f $frag_size_file $images_dir");
    $frag_size_content .= "\n<div id=\"fragment_size_legend\" style=\"width: 750px;\">";
    $frag_size_content .= "\n<P CLASS=\"Indented24LR_s16\">The distribution of fragment sizes is plotted for all paired-end reads in which both pairs could be aligned to the transcriptome. Pairs with reads aligning to different gene loci or different chromosomes are excluded. Fragment sizes with less than 0.01% of total reads are not plotted.</P><BR>";
    $frag_size_content .= "\n</div>";
    $frag_size_content .= "\n<P CLASS=\"Indented24LR\"><IMG SRC=\"$frag_size_file2\" CLASS=\"Pic_unbordered\" ALT=\"Fragment size distribution\" TITLE=\"Fragment size distribution\" WIDTH=\"750\" HEIGHT=\"750\"></P><BR>";
  }else{
    print YELLOW, "\n\n($content_count) Could not find figure file: $frag_size_file\n\n", RESET;
    $frag_size_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $frag_size_content .= "\n<BR></div><BR>";


  #Figure content - Library complexity
  $content_count++;
  $box_id = "box"."$content_count"."h";
  my $library_complexity_content = "\n<!-- Figure - Library Complexity -->";
  $image_desc = "Summary of library complexity - estimated by tag redundancy per million reads and compared to other libraries";
  $library_complexity_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">$image_desc</P>";
  $library_complexity_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $library_complexity_content .= "\n<div id=\"$box_id\">";

  my $library_complexity_SEQ_file = "$qual_stats_path"."LibraryComplexity_SEQ.svg"; #750x750
  my $library_complexity_SEQ_file2 = "images/$lib/"."LibraryComplexity_SEQ.svg";

  my $library_complexity_MAP_file = "$qual_stats_path"."LibraryComplexity_MAP.svg"; #750x750
  my $library_complexity_MAP_file2 = "images/$lib/"."LibraryComplexity_MAP.svg";

  $image_width = 750;
  $image_height = 750;
  if (-e $library_complexity_SEQ_file && -e $library_complexity_MAP_file){
    system("cp -f $library_complexity_SEQ_file $images_dir");
    system("cp -f $library_complexity_MAP_file $images_dir");
    $library_complexity_content .= "\n<div id=\"library_complexity_legend\" style=\"width: 850px;\">";
    $library_complexity_content .= "\n<P CLASS=\"Indented24LR_s16\">Library complexity is calculated for the sequence library by randomly sampling 1 million reads and determining the number of unique and redundant reads
                                                                    within the pool.  This sampling is repeated (with replacement) at least 3 times and average values across these samples is used.  A 'redundant' read is 
                                                                    one where both reads of a pair have either identical sequence (first figure) or mapping location (second figure) to at least one other read in the sampled pool.  
                                                                    In each panel below, the value for the current library is indicated by a red dot and compared to other libraries, summarized as box and whisker plots.  
                                                                    In the first panel (green), the percent of unique reads per million reads sampled is summarized.  
                                                                    In the second panel (orange), the percent of redundant reads per million is summarized.  In the third panel (yellow), the redundant reads are 
                                                                    further examined to determine the number of distinct redundant reads.  For example, if all redudant reads corresponded to a few sequences occuring many times
                                                                    this would result in a low number of distinct redundant reads.  A  'good' library with high complexity should have high, low, and high values for panels 
                                                                    one, two and three respectively.</P><BR>";
    $library_complexity_content .= "\n</div>";
    $library_complexity_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Library complexity assessed by read sequence identity</P>";
    $library_complexity_content .= "\n<P CLASS=\"Indented24LR\"><OBJECT DATA=\"$library_complexity_SEQ_file2\" TYPE=\"image/svg+xml\" WIDTH=\"$image_width\" HEIGHT=\"$image_height\" TITLE=\"$image_desc\" CLASS=\"Pic_bordered\" STANDBY=\"$image_desc\"><EMBED SRC=\"$library_complexity_SEQ_file2\" TYPE=\"image/svg+xml\" WIDTH=\"$image_width\" HEIGHT=\"$image_height\" /></OBJECT></P><BR>";
    $library_complexity_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Library complexity assessed by read sequence mapping positions</P>";
    $library_complexity_content .= "\n<P CLASS=\"Indented24LR\"><OBJECT DATA=\"$library_complexity_MAP_file2\" TYPE=\"image/svg+xml\" WIDTH=\"$image_width\" HEIGHT=\"$image_height\" TITLE=\"$image_desc\" CLASS=\"Pic_bordered\" STANDBY=\"$image_desc\"><EMBED SRC=\"$library_complexity_MAP_file2\" TYPE=\"image/svg+xml\" WIDTH=\"$image_width\" HEIGHT=\"$image_height\" /></OBJECT></P><BR>";

  }else{
    print YELLOW, "\n\n($content_count) Could not find figure files: $library_complexity_SEQ_file or $library_complexity_MAP_file\n\n", RESET;
    $library_complexity_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $library_complexity_content .= "\n<BR></div><BR>";


  #Figure content - Gene coverage at increasing library depth levels
  $content_count++;
  $box_id = "box"."$content_count"."h";
  my $gene_coverage_content = "\n<!-- Figure - Gene coverage -->";
  $image_desc = "Distribution % of gene bases covered for each expressed gene (at various minimum depth levels)";
  $gene_coverage_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">$image_desc</P>";
  $gene_coverage_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $gene_coverage_content .= "\n<div id=\"$box_id\">";
  my $gene_coverage_file = "$mapping_stats_path"."PercentGeneCoverage_1x-500x_Expressed.jpeg"; #700x700
  my $gene_coverage_file2 = "images/$lib/"."PercentGeneCoverage_1x-500x_Expressed.jpeg";
  $image_width = 700;
  $image_height = 700;
  if (-e $gene_coverage_file){
    system("cp -f $gene_coverage_file $images_dir");
    $gene_coverage_content .= "\n<div id=\"position_bias_legend\" style=\"width: 700px;\">";
    $gene_coverage_content .= "\n<P CLASS=\"Indented24LR_s16\">The distribution of percent coverage levels for each gene is summarized as a box and whisker plot.  These plots are produced for six minimum depth requirements.
                                                               For example, in the first plot (>=1X), the percentage of bases covered to a depth of 1X or greater is determined for each gene.  If a gene is completely covered by
                                                               reads from the first base to the last, at a depth of 1X or greater, then this gene is given a value of 100%.  The distribution of these percent values for all genes
                                                               detected above background is summarized by the box and whisker plot.</P><BR>";
    $gene_coverage_content .= "\n</div>";
    $gene_coverage_content .= "\n<P CLASS=\"Indented24LR\"><IMG SRC=\"$gene_coverage_file2\" CLASS=\"Pic_unbordered\" ALT=\"$image_desc\" TITLE=\"$image_desc\" WIDTH=\"$image_width\" HEIGHT=\"$image_width\"></P><BR>";

  }else{
    print YELLOW, "\n\n($content_count) Could not find figure file: $gene_coverage_file\n\n", RESET;
    $gene_coverage_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $gene_coverage_content .= "\n<BR></div><BR>";



  #Figure1 content
  $content_count++;
  $box_id = "box"."$content_count"."h";
  my $fig1_content = "\n<!-- Figure 1 -->";
  $fig1_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Distribution of log2 raw expression values for each feature type</P>";
  $fig1_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $fig1_content .= "\n<div id=\"$box_id\">";

  my $exp_raw_file = "$exp_stats_path"."ExpressionValues_RAW_Log2_AllTypes_bplot.jpeg"; #1000x1000
  my $exp_raw_file2 = "images/$lib/"."ExpressionValues_RAW_Log2_AllTypes_bplot.jpeg";
  if (-e $exp_raw_file){
    system("cp -f $exp_raw_file $images_dir");
    $fig1_content .= "\n<div id=\"fig1_legend\" style=\"width: 1000px;\">";
    $fig1_content .= "\n<P CLASS=\"Indented24LR_s16\">Box-and-whisker plots for log2 expression values for each feature type.</P><BR>";
    $fig1_content .= "\n</div>";
    $fig1_content .= "\n<P CLASS=\"Indented24LR\"><IMG SRC=\"$exp_raw_file2\" CLASS=\"Pic_unbordered\" ALT=\"Raw expression values by feature type\" TITLE=\"Raw expression values by feature type\" WIDTH=\"1000\" HEIGHT=\"1000\"></P><BR>";
  }else{
    print YELLOW, "\n\n($content_count) Could not find figure file: $exp_raw_file\n\n", RESET;
    $fig1_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $fig1_content .= "\n<BR></div><BR>";


  #Figure2 content
  $content_count++;
  $box_id = "box"."$content_count"."h";
  my $fig2_content = "\n<!-- Figure 2 -->";
  $fig2_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Distribution of log2 normalized expression values for each feature type</P>";
  $fig2_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $fig2_content .= "\n<div id=\"$box_id\">";

  my $exp_norm_file = "$exp_stats_path"."ExpressionValues_NORM1_Log2_AllTypes_bplot.jpeg"; #1000x1000
  my $exp_norm_file2 = "images/$lib/"."ExpressionValues_NORM1_Log2_AllTypes_bplot.jpeg";
  if (-e $exp_norm_file){
    system("cp -f $exp_norm_file $images_dir");
    $fig2_content .= "\n<div id=\"fig2_legend\" style=\"width: 1000px;\">";
    $fig2_content .= "\n<P CLASS=\"Indented24LR_s16\">Box-and-whisker plots for normalized log2 expression values for each feature type.</P><BR>";
    $fig2_content .= "\n</div>";
    $fig2_content .= "\n<P CLASS=\"Indented24LR\"><IMG SRC=\"$exp_norm_file2\" CLASS=\"Pic_unbordered\" ALT=\"Normalized expression values by feature type\" TITLE=\"Normalized expression values by feature type\" WIDTH=\"1000\" HEIGHT=\"1000\"></P><BR>";
  }else{
    print YELLOW, "\n\n($content_count) Could not find figure file: $exp_norm_file\n\n", RESET;
    $fig2_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $fig2_content .= "\n<BR></div><BR>";


  #Figure3 content
  $content_count++;
  $box_id = "box"."$content_count"."h";
  my $fig3_content = "\n<!-- Figure 3 -->";
  $fig3_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Density scatter plot of exon region versus gene expression values</P>";
  $fig3_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $fig3_content .= "\n<div id=\"$box_id\">";

  my $er_vs_gene_file = "$exp_stats_path"."ExonRegionVsGeneExpression_Log2.jpeg"; #700x700
  my $er_vs_gene_file2 = "images/$lib/"."ExonRegionVsGeneExpression_Log2.jpeg";
  if (-e $er_vs_gene_file){
    system("cp -f $er_vs_gene_file $images_dir");
    $fig3_content .= "\n<div id=\"fig3_legend\" style=\"width: 700px;\">";
    $fig3_content .= "\n<P CLASS=\"Indented24LR_s16\">Density scatter plot of log2 expression values for exon regions versus corresponding gene expression values.</P><BR>";
    $fig3_content .= "\n</div>";
    $fig3_content .= "\n<P CLASS=\"Indented24LR\"><IMG SRC=\"$er_vs_gene_file2\" CLASS=\"Pic_unbordered\" ALT=\"Exon region vs. gene expression\" TITLE=\"Exon region vs. gene expression\" WIDTH=\"700\" HEIGHT=\"700\"></P><BR>";
  }else{
    print YELLOW, "\n\n($content_count) Could not find figure file: $er_vs_gene_file\n\n", RESET;
    $fig3_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $fig3_content .= "\n<BR></div><BR>";


  #Figure4 content
  $content_count++;
  $box_id = "box"."$content_count"."h";
  my $fig4_content = "\n<!-- Figure 4 -->";
  $fig4_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Density scatter plot of silent intron region versus gene expression values</P>";
  $fig4_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $fig4_content .= "\n<div id=\"$box_id\">";

  my $si_vs_gene_file = "$exp_stats_path"."SilentIntronicVsGeneExpression_Log2.jpeg"; #700x700
  my $si_vs_gene_file2 = "images/$lib/"."SilentIntronicVsGeneExpression_Log2.jpeg";
  if (-e $si_vs_gene_file){
    system("cp -f $si_vs_gene_file $images_dir");
    $fig4_content .= "\n<div id=\"fig4_legend\" style=\"width: 700px;\">";
    $fig4_content .= "\n<P CLASS=\"Indented24LR_s16\">Density scatter plot of log2 expression values for silent intron regions versus corresponding gene expression values.</P><BR>";
    $fig4_content .= "\n</div>";
    $fig4_content .= "\n<P CLASS=\"Indented24LR\"><IMG SRC=\"$si_vs_gene_file2\" CLASS=\"Pic_unbordered\" ALT=\"Silent intron region vs. gene expression\" TITLE=\"Silent intron region vs. gene expression\" WIDTH=\"700\" HEIGHT=\"700\"></P><BR>";
  }else{
    print YELLOW, "\n\n($content_count) Could not find figure file: $si_vs_gene_file\n\n", RESET;
    $fig4_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $fig4_content .= "\n<BR></div><BR>";


  #Figure5 content
  $content_count++;
  $box_id = "box"."$content_count"."h";
  my $fig5_content = "\n<!-- Figure 5 -->";
  $fig5_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Distribution of gene-by-gene expression cutoff values</P>";
  $fig5_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $fig5_content .= "\n<div id=\"$box_id\">";

  my $cutoffs_file = "$exp_stats_path"."DistributionOfCutoffs.jpeg"; #750x750
  my $cutoffs_file2 = "images/$lib/"."DistributionOfCutoffs.jpeg";
  if(-e $cutoffs_file){
    system("cp -f $cutoffs_file $images_dir");
    $fig5_content .= "\n<div id=\"fig5_legend\" style=\"width: 750px;\">";
    $fig5_content .= "\n<P CLASS=\"Indented24LR_s16\">Histogram depicting the distribution of expression cutoff values used for each gene.  In order to be considered 'expressed above background', all features 
                      (genes, transcripts, exons, junctions, etc.) must be expressed above the level of INTERGENIC noise.
                      This INTERGENIC cutoff is the 95th percentile of expression values for all Silent Intergenic Regions in this library and is depicted below as a dotted red line.
                      The number of genes for which only INTERGENIC noise is considered is provided in the legend.
                      For features within the boundaries of highly expressed genes, additional noise is expected due to the presence of un-processed RNA contamination.  
                      For this reason, a higher INTRAGENIC cutoff is determined. These are calculated by fitting a linear model to the 95th percentile of expression values for silent intronic regions.
                      The INTRAGENIC cutoff for a gene is then determined by using the gene expression level and the coefficients of the model fit.  The distribution of the resulting gene-by-gene cutoffs is 
                      depicted below as a histogram.  The number of genes that required an INTRAGENIC cutoff is also indicated in the legend.</P><BR>";
    $fig5_content .= "\n</div>";
    $fig5_content .= "\n<P CLASS=\"Indented24LR\"><IMG SRC=\"$cutoffs_file2\" CLASS=\"Pic_unbordered\" ALT=\"Distribution of cutoffs values\" TITLE=\"Distribution of cutoffs values\" WIDTH=\"750\" HEIGHT=\"750\"></P><BR>";
  }else{
    print YELLOW, "\n\n($content_count) Could not find figure file: $cutoffs_file\n\n", RESET;
    $fig5_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $fig5_content .= "\n<BR></div><BR>";


  #Figure6 content
  $content_count++;
  $box_id = "box"."$content_count"."h";
  my $fig6_content = "\n<!-- Figure 6 -->";
  $fig6_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Histograms of expression values for each feature type</P>";
  $fig6_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $fig6_content .= "\n<div id=\"$box_id\">";

  my $exp_all_features_file = "$exp_stats_path"."ExpressionValues_Log2_AllTypes_hist.jpeg"; #1920x1200
  my $exp_all_features_file2 = "images/$lib/"."ExpressionValues_Log2_AllTypes_hist.jpeg";
  if (-e $exp_all_features_file){
    system("cp -f $exp_all_features_file $images_dir");
    $fig6_content .= "\n<div id=\"fig6_legend\" style=\"width: 1920px;\">";
    $fig6_content .= "\n<P CLASS=\"Indented24LR_s16\">Histograms depicting distribution of log2 expression values for individual feature types. The 95th percentile of expression values for silent intergenic region is depicted as a dotted line on all plots.</P><BR>";
    $fig6_content .= "\n</div>";
    $fig6_content .= "\n<P CLASS=\"Indented24LR\"><IMG SRC=\"$exp_all_features_file2\" CLASS=\"Pic_unbordered\" ALT=\"Expression distribution by feature type\" TITLE=\"Expression distribution by feature type\" WIDTH=\"1920\" HEIGHT=\"1200\"></P><BR>";
  }else{
    print YELLOW, "\n\n($content_count) Could not find figure file: $exp_all_features_file\n\n", RESET;
    $fig6_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $fig6_content .= "\n<BR></div><BR>";


  #Figure7 content
  $content_count++;
  $box_id = "box"."$content_count"."h";
  my $fig7_content = "\n<!-- Figure 7 -->";
  $fig7_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Percentiles plot for expression of exon regions, silent intron regions and silent intergenic regions</P>";
  $fig7_content .= "<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $fig7_content .= "<div id=\"$box_id\">";

  my $percentiles_file = "$exp_stats_path"."PercentilesPlot_Log2.jpeg"; #800x800
  my $percentiles_file2 = "images/$lib/"."PercentilesPlot_Log2.jpeg";
  if (-e $percentiles_file){
    system("cp -f $percentiles_file $images_dir");
    $fig7_content .= "\n<div id=\"fig7_legend\" style=\"width: 800px;\">";
    $fig7_content .= "\n<P CLASS=\"Indented24LR_s16\">Percentiles plot for exon region, silent intron region and silent intergenic region expression values.  The 95th percentiles of intronic and intergenic distributions are depicted as colored, dotted lines</P><BR>";
    $fig7_content .= "\n</div>";
    $fig7_content .= "\n<P CLASS=\"Indented24LR\"><IMG SRC=\"$percentiles_file2\" CLASS=\"Pic_unbordered\" ALT=\"Percentiles plot\" TITLE=\"Percentiles plot\" WIDTH=\"800\" HEIGHT=\"800\"></P><BR>";
  }else{
    print YELLOW, "\n\n($content_count) Could not find figure file: $percentiles_file\n\n", RESET;
    $fig7_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $fig7_content .= "\n<BR></div><BR>";


  #Figure8 content
  $content_count++;
  $box_id = "box"."$content_count"."h";
  my $fig8_content = "\n<!-- Figure 8 -->";
  $fig8_content .= "\n<P CLASS=\"Indented24LR_s16_bold\">Cumulative distribution of mean Phred score for all reads</P>";
  $fig8_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[$box_id]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $fig8_content .= "\n<div id=\"$box_id\">";

  my $phred_file = "$qual_stats_path"."$lib"."_phred.jpeg"; #750x750
  my $phred_file2 = "images/$lib/"."$lib"."_phred.jpeg";
  if (-e $phred_file){
    system("cp -f $phred_file $images_dir");
    $fig8_content .= "\n<div id=\"fig8_legend\" style=\"width: 800px;\">";
    $fig8_content .= "\n<P CLASS=\"Indented24LR_s16\">CDF plot depicting cumulative distribution of mean Phred score for all reads broken down by R1/R2 and sequencing lane</P><BR>";
    $fig8_content .= "\n</div>";
    $fig8_content .= "\n<P CLASS=\"Indented24LR\"><IMG SRC=\"$phred_file2\" CLASS=\"Pic_unbordered\" ALT=\"Phred CDF plot\" TITLE=\"Phred CDF plot\" WIDTH=\"800\" HEIGHT=\"800\"></P><BR>";
  }else{
    print YELLOW, "\n\n($content_count) Could not find figure file: $phred_file\n\n", RESET;
    $fig8_content .= "\n<P CLASS=\"Indented24LR_s16\">Not available ...</P><BR>";
  }
  $fig8_content .= "\n<BR></div><BR>";


###########################################################################################################################################################################################################
#Generate final content - Each content block below has an html comment, title line, and expandable button section
  my $summary_content = '';
$summary_content = <<"EOF";

<!--- TABLES SECTION --->

$li_content

$lane_stats_content 

$ras_content

$mr_content

$gac_content

$ef_content

$sk_content

$signal_to_noise_content

$percentiles_content

<!-- FIGURES SECTION -->

$library_complexity_content

$position_bias_content

$frag_size_content

$gene_coverage_content

$fig1_content

$fig2_content

$fig3_content

$fig4_content

$fig5_content

$fig6_content

$fig7_content

$fig8_content

<BR><BR>

EOF
#############################################################################################################################################################################################################


  #8.) Write page
  my $meta_description = "Provides a description of statistics for a single sequencing library";
  my $meta_keywords = "Paired-end sequencing, Illumina, Solexa";
  my $web_path = "$target_web_dir"."$lib".".htm";
  my $title = "Statistics and figures for sequencing library: $lib ($lib_name)";
  &writePage('-path'=>$web_path, '-title'=>$title, '-content'=>\$summary_content, '-css_path'=>"../ALEXA2.css", '-alexa_home_path'=>"$alexa_home_path", '-alexa_seq_home_path'=>"$alexa_seq_path", '-summary_path'=>"Summary.htm", '-genes_path'=>"genes/index.html", '-search_path'=>$search_page_url, '-meta_description'=>$meta_description, '-meta_keywords'=>$meta_keywords, '-google_analytics'=>1, '-google_analytics_id'=>$google_analytics_id, '-collapse_div_script'=>'../animatedcollapse.js', '-jquery_min_script'=>'../jquery.min.js', '-div_count'=>$content_count);

}

exit();

