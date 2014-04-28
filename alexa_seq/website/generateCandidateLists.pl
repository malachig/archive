#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Provide 'Master' candidate lists that summarize DE and AE events at the gene level (i.e. at most one entry per gene)
#For any comparison, provide such a list for 'ALL' events, 'Library A specific' events, and 'Library B specific' events
#This list should be ranked by fold-change of either the gene or the feature with the highest SI value

#1.) Get all DE genes and all AE features 
#    - Grand RANK according to gene DE or highest feature seq DE (or SI?)
#    - Use sign of DE values to divide into candidated into library A/B specific groups

#2.) Info to include in each output file
# Gene_id
# Ensembl_gene_id
# Ensembl_gene_name
# Entrez_gene_name
# Gene Type (i.e. protein coding)
# Gene: DE, Library A Exp, Library B Exp
# Top Feature: SI, DE, Percent Seq DE, Library A Exp, Library B Exp
# ORF affected?
# Class (e.g. Gene DE, AE-exon_skipping, AE-alternative_boundary, AE-intron_retention, AE-TSS_PolyA) - each gene can be assigned multiple classes?
# Number AE features (#for each feature type: transcript, exon, knownjunction, noveljunction, etc.)
# Total Number AE features
# Reciprocity


use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use Benchmark;

#Load the ALEXA libraries
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);
use website::WEB qw(:all);


#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $project = '';
my $comparison_id = '';
my $comparison_name = '';
my $ensembl_version = '';
my $analysis_dir = '';
my $annotation_dir = '';
my $junction_seq_size = '';
my $boundary_seq_size = '';
my $partition_file = '';
my $alexa_home_path = '';
my $alexa_seq_path = '';
my $search_page_url = '';
my $google_analytics_id = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
            'junction_seq_size=i'=>\$junction_seq_size, 'boundary_seq_size=i'=>\$boundary_seq_size,
            'project=s'=>\$project, 'comparison_id=s'=>\$comparison_id, 'comparison_name=s'=>\$comparison_name, 'ensembl_version=s'=>\$ensembl_version, 
            'analysis_dir=s'=>\$analysis_dir, 'annotation_dir=s'=>\$annotation_dir, 'partition_file=s'=>\$partition_file, 
            'alexa_home_path=s'=>\$alexa_home_path, 'alexa_seq_path=s'=>\$alexa_seq_path, 'search_page_url=s'=>\$search_page_url, 'google_analytics_id=s'=>\$google_analytics_id);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the project name using:  --project", RESET;
print GREEN, "\n\tSpecify the comparison id used for stats files using:  --comparison_id", RESET;
print GREEN, "\n\tSpecify the comparison name used for stats files using:  --comparison_name", RESET;
print GREEN, "\n\tSpecify the ensembl version using:  --ensembl_version", RESET;
print GREEN, "\n\tSpecify the root path to the analysis dir using:  --analysis_dir", RESET;
print GREEN, "\n\tSpecify the root path to the annotation dir used for the analysis using:  --annotation_dir", RESET;
print GREEN, "\n\tSpecify the junction database sequence length using:  --junction_seq_size", RESET;
print GREEN, "\n\tSpecify the boundary database sequence length using:  --boundary_seq_size", RESET;
print GREEN, "\n\tSpecify the genome partition file using:  --partition_file", RESET;
print GREEN, "\n\tSpecify the URL to your ALEXA home page using: --alexa_home_path", RESET;
print GREEN, "\n\tSpecify the URL to your ALEXA-Seq results page using: --alexa_seq_path", RESET;
print GREEN, "\n\tSpecify the URL to your Xapian-Omega search page using:  --search_page_url", RESET;
print GREEN, "\n\tSpecify your Google Analytics ID using: --google_analytics_id", RESET;

print GREEN, "\n\nExample: generateCandidateLists.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --project=Breast  --comparison_id=HS1188_vs_HS1187  --comparison_name=Myo_Epi_vs_Lum_Epi  --ensembl_version=53  --analysis_dir=/projects/malachig/alexa_seq/  --annotation_dir=/projects/alexa/sequence_databases/hs_53_36o/  --junction_seq_size=62  --boundary_seq_size=62  --partition_file=/projects/alexa/sequence_databases/hs_53_36o/Regions_250_Genes.txt  --alexa_home_path=http://www.alexaplatform.org/index.htm  --alexa_seq_path=http://www.alexaplatform.org/alexa_seq/results.htm  --search_page_url=http://www.bcgsc.ca/xapian-search/omega  --google_analytics_id=UA-xxxxxx-x\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $project && $comparison_id && $comparison_name && $ensembl_version && $analysis_dir && $annotation_dir && $junction_seq_size && $boundary_seq_size && $partition_file && $alexa_home_path && $alexa_seq_path && $search_page_url && $google_analytics_id){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}


#Form file and directory paths needed and check them
$analysis_dir = &checkDir('-dir'=>$analysis_dir, '-clear'=>"no");
my $stats_dir = "$analysis_dir"."figures_and_stats/";
$stats_dir = &checkDir('-dir'=>$stats_dir, '-clear'=>"no");
my $de_dir = "$stats_dir"."DE/"."$project";
my $si_dir = "$stats_dir"."SI/"."$project";
$de_dir = &checkDir('-dir'=>$de_dir, '-clear'=>"no");
$si_dir = &checkDir('-dir'=>$si_dir, '-clear'=>"no");
$annotation_dir = &checkDir('-dir'=>$annotation_dir, '-clear'=>"no");
my $web_dir = "$analysis_dir"."temp/website/$project/";
$web_dir = &checkDir('-dir'=>$web_dir, '-clear'=>"no");
my $summary_file = "$web_dir"."data/"."$comparison_name".".txt";
my $web_path ="$web_dir"."$comparison_id".".htm";

#Import the genome partitions
my $partitions_ref = &getPartitions('-file'=>$partition_file);


#The 'Score' used to rank candidate will be based on various criteria and subject to bonuses.  Bonuses will increment the score by the following ratio\
my $bonus_ratio = 0.25;


#Get basic gene info
print BLUE, "Getting basic gene info:", RESET;
my $genes_ref;
&getBasicGeneInfo();

#Keep a list of all gene IDs that have at least one DE or AE feature
my %master_gene_list;

#Get the DE values to be considered
print BLUE, "\n\nGetting Gene DE values:", RESET;
my $de_ref = &getDeValues('-de_dir'=>$de_dir, '-ensembl_version'=>$ensembl_version, '-comparison_name'=>$comparison_name);
my $de_count = keys %{$de_ref};
print BLUE, "\n\tFound $de_count Gene DE records", RESET;


#Get the AE values to be considered
print BLUE, "\n\nGetting all feature SI values:", RESET;
my $ae_ref = &getAeValues('-ae_dir'=>$si_dir, '-ensembl_version'=>$ensembl_version, '-comparison_name'=>$comparison_name);
my $ae_count = keys %{$ae_ref};
print BLUE, "\n\tFound $ae_count Feature AE records", RESET;


#Get annotations so that every feature can be tied to gene_id, gene name, etc...
print BLUE, "\n\nGetting feature annotations:", RESET;
my $feature_ref = &getFeatures('-annotation_dir'=>$annotation_dir, '-junction_seq_size'=>$junction_seq_size, '-boundary_seq_size'=>$boundary_seq_size);
my $feature_count = keys %{$feature_ref};
print BLUE, "\n\tFound $feature_count feature annotation records for the DE and SI records", RESET;


#Examine the adjacency of AE events to sub-classify events and assign bonus scores for events with multiple supporting features
print BLUE, "\n\nExamining adjacency to identify AE events with multiple support:", RESET;
&examineAdjacency();


#Build the ranked candidate list combining both DE and AE events into one grand ranked list
print BLUE, "\n\nBuilding ranked candidate list combining both DE and AE events into one grand ranked list:", RESET;
my $candidates_ref = &buildCandidateList('-feature_ref'=>$feature_ref, '-de_ref'=>$de_ref, '-ae_ref'=>$ae_ref);


#Generate a string of all AE feature names (ordered by chr coords)
print BLUE, "\n\nGetting name string for AE features of each gene:", RESET;
&generateFeatureNameStrings();


#Apply bonus values to scores for each candidate
print BLUE, "\n\nAdding bonus value to scores for each candidate:", RESET;
while (my ($gene_id) = each %{$candidates_ref}){
  $candidates_ref->{$gene_id}->{adjusted_score} = $candidates_ref->{$gene_id}->{score} + $candidates_ref->{$gene_id}->{bonus};
}


#Now apply the grand rank to each candidate
print BLUE, "\n\nDetermining rank of each candidate:", RESET;
my $n = 0;
foreach my $gene_id (sort {$candidates_ref->{$b}->{adjusted_score} <=> $candidates_ref->{$a}->{adjusted_score}} keys %{$candidates_ref}){
  $n++;
  $candidates_ref->{$gene_id}->{rank} = $n;
}


#Print a grand summary file that can be loaded into excel, filtered, used for downstream analysis, etc.
#Provide a link to this file in the candidate list page
print BLUE, "\n\nPrinting complete candidate list file: $summary_file", RESET;
open(SUMMARY, ">$summary_file") || die "\n\nCould not open summary output file: $summary_file\n\n";
print SUMMARY "rank\tadjusted_score\tbonus\talexa_gene_id\tensembl_g_id\tgene_name\tgene_type\tknown_transcript_count\texon_region_count\tmain_type\tdirection\tfold_change\tae_event_count\tae_exons\tae_exons_included_j\tae_exons_skipped_j\tae_alternative_boundaries\tae_retained_introns\tae_cryptic_exons\tae_codes\ttop_fid\ttop_seq_name\treciprocity\tadjacency_transcript\tadjacency_percent_transcript\n";
foreach my $gene_id (sort {$candidates_ref->{$b}->{adjusted_score} <=> $candidates_ref->{$a}->{adjusted_score}} keys %{$candidates_ref}){
  my $adjusted_score = sprintf("%.2f", $candidates_ref->{$gene_id}->{adjusted_score});
  my $bonus = sprintf("%.2f", $candidates_ref->{$gene_id}->{bonus});
  my $fold_change = sprintf("%.2f", $candidates_ref->{$gene_id}->{fold_change});
  my $reciprocity = $candidates_ref->{$gene_id}->{reciprocity};
  unless ($reciprocity eq "N/A"){
    $reciprocity = sprintf("%.2f", $candidates_ref->{$gene_id}->{reciprocity});
  }
  print SUMMARY "$candidates_ref->{$gene_id}->{rank}\t$adjusted_score\t$bonus\t$gene_id\t$candidates_ref->{$gene_id}->{ensembl_g_id}\t$candidates_ref->{$gene_id}->{gene_name}\t$genes_ref->{$gene_id}->{gene_type}\t$genes_ref->{$gene_id}->{transcript_count}\t$genes_ref->{$gene_id}->{exon_count}\t$candidates_ref->{$gene_id}->{main_type}\t$candidates_ref->{$gene_id}->{direction}\t$fold_change\t$candidates_ref->{$gene_id}->{ae_event_count}\t$candidates_ref->{$gene_id}->{ae_exons}\t$candidates_ref->{$gene_id}->{ae_exons_included_j}\t$candidates_ref->{$gene_id}->{ae_exons_skipped_j}\t$candidates_ref->{$gene_id}->{ae_alternative_boundaries}\t$candidates_ref->{$gene_id}->{ae_introns_retained}\t$candidates_ref->{$gene_id}->{ae_cryptic_exons}\t$candidates_ref->{$gene_id}->{ae_codes}\t$candidates_ref->{$gene_id}->{fid}\t$candidates_ref->{$gene_id}->{seq_name}\t$reciprocity\t$genes_ref->{$gene_id}->{adjacency_transcript}\t$genes_ref->{$gene_id}->{adjacency_percent_transcript}\n";

}
close(SUMMARY);


#Summarize the number of genes with various criteria:
#Total significant genes, DE genes, AE genes, AE exon/junction genes, AE exon skipping genes, AE alternative boundary genes, AE intron retention genes, AE cryptic exon genes
print BLUE, "\n\nSummarizing number of genes in various categories:", RESET;
my $candidate_genes_count = 0;
my $de_genes_count = keys %{$de_ref};
my $ae_genes_count = 0;
my $ae_eu_genes_count = 0;
my $ae_es_genes_count = 0;
my $ae_ab_genes_count = 0;
my $ae_ir_genes_count = 0;
my $ae_ce_genes_count = 0;
foreach my $gene_id (sort {$candidates_ref->{$b}->{adjusted_score} <=> $candidates_ref->{$a}->{adjusted_score}} keys %{$candidates_ref}){
  $candidate_genes_count++;
  if ($candidates_ref->{$gene_id}->{ae_event_count} > 0){

    $ae_genes_count++;
    if ($candidates_ref->{$gene_id}->{ae_exons} > 0 || $candidates_ref->{$gene_id}->{ae_exons_included_j} > 0){$ae_eu_genes_count++;}
    if ($candidates_ref->{$gene_id}->{ae_exons_skipped_j} > 0){$ae_es_genes_count++;}
    if ($candidates_ref->{$gene_id}->{ae_alternative_boundaries} > 0){$ae_ab_genes_count++;}
    if ($candidates_ref->{$gene_id}->{ae_introns_retained} > 0){$ae_ir_genes_count++;}
    if ($candidates_ref->{$gene_id}->{ae_cryptic_exons} > 0){$ae_ce_genes_count++;}
  }
}
$candidate_genes_count = &commify($candidate_genes_count);
$de_genes_count = &commify($de_genes_count);
$ae_genes_count = &commify($ae_genes_count);
$ae_eu_genes_count = &commify($ae_eu_genes_count);
$ae_es_genes_count = &commify($ae_es_genes_count);
$ae_ab_genes_count = &commify($ae_ab_genes_count);
$ae_ir_genes_count = &commify($ae_ir_genes_count);
$ae_ce_genes_count = &commify($ae_ce_genes_count);

#Generate link to master candidate gene list file
my $data_name = "$comparison_name".".txt";
my $data_path = "data/$data_name";
my $data_link = "<A HREF=\"$data_path\">$data_name</A>";
my $summary_content = '';
$summary_content = "<P CLASS=\"Indented12LR_s16_bold\">Download complete candidate list as tab delimited text file: $data_link</P><BR>\n";

#Summarize counts
my $total_gene_count = keys %{$genes_ref};
$total_gene_count = &commify($total_gene_count);
$summary_content .= "<P CLASS=\"Indented12LR_s16_Bold\">Summary of Differential (DE) and Alternative Expression (AE) for all gene loci:</P>\n";
$summary_content .= "<P CLASS=\"Indented12LR_s16\">Total Candidate Genes: $candidate_genes_count (of $total_gene_count possible genes)</P>\n";
$summary_content .= "<P CLASS=\"Indented12LR_s16\">DE Genes: $de_genes_count</P>\n";
$summary_content .= "<P CLASS=\"Indented12LR_s16\">AE Genes: $ae_genes_count</P>\n";
$summary_content .= "<P CLASS=\"Indented24LR_s16\">Alternative Exon Usage (EU) Genes: $ae_eu_genes_count</P>\n";
$summary_content .= "<P CLASS=\"Indented24LR_s16\">Alernative Exon Skipping (ES) Genes: $ae_es_genes_count</P>\n";
$summary_content .= "<P CLASS=\"Indented24LR_s16\">Alternative Exon Boundary (AB) Genes: $ae_ab_genes_count</P>\n";
$summary_content .= "<P CLASS=\"Indented24LR_s16\">Intron Retention (IR) Genes: $ae_ir_genes_count</P>\n";
$summary_content .= "<P CLASS=\"Indented24LR_s16\">Cryptic Exon (CE) Genes: $ae_ce_genes_count</P><BR>\n";


#Now generate html tables for display in the data viewer (for these, only the topN will be displayed):
#1.) Main candidate list
#2,) AE events only
#3.) Gains only
#4.) Losses only
print BLUE, "\n\nGenerating HTML table content:", RESET;
my $top_n = 250;
my $main_list_content = '';
my $main_list_count = 0;
my $ae_list_content = '';
my $ae_list_count = 0;
my $gain_list_content = '';
my $gain_list_count = 0;
my $loss_list_content = '';
my $loss_list_count = 0;

#Define the header row
my $header_row ="<TR><TD CLASS=\"Head3\">Rank</TD><TD CLASS=\"Head3\">Overall<BR>Rank</TD><TD CLASS=\"Head3\">Score</TD><TD CLASS=\"Head3\">Name</TD><TD CLASS=\"Head3\">Gene<BR>Type</TD><TD CLASS=\"Head3\">Trans.<BR>Count</TD><TD CLASS=\"Head3\">Exon<BR>Count</TD><TD CLASS=\"Head3\">Event<BR>Type</TD><TD CLASS=\"Head3\">Direction</TD><TD CLASS=\"Head3\">FC</TD><TD CLASS=\"Head3\"># AE<BR>Events</TD><TD CLASS=\"Head3\">AE Codes</TD><TD CLASS=\"Head3\">Top Feature</TD><TD CLASS=\"Head3\">Adjacency<BR>%</TD></TR>\n";
$main_list_content .= "$header_row";
$ae_list_content .= "$header_row";
$gain_list_content .= "$header_row";
$loss_list_content .= "$header_row";

#Define the legend
my $legend = "<div id=\"legend\" style=\"width: 1000px;\">\n";
$legend .= "<P CLASS=\"Indented24LR_s16\">The following table lists candidate genes ranked according to 'Score'.  Only the top $top_n genes appear in each list.  For a complete list download the results file provided above.  To appear in this list, a 'candidate gene' must be differentially expressed at the gene level (DE) or contain at least one alternatively expressed (AE) feature.  The 'Score' is based on the fold change as well as other factors for AE events (multiple features supporting the same event, reciprocity, etc.).  The 'Event Type' indicates whether the gene was differentially expressed ('DE') or alternatively expressed 'AE'.  'Direction' refers to the direction (+ve/-ve) of the fold change.  'Gain' (+ve) indicates higher expression in library A than library B and vice versa.  'FC' stands for fold change.  '# AE Events' is the number of significant alternative expression events observed for each gene.  'AE Codes' refers to various alternative expression subtypes (listed above).  The 'Top Feature' is the most significant AE event observed for an AE gene.  'Adjacency' refers to the number of consecutive AE features (potentially supporting the same isoform).  'Adjacency %' refers to the percent of all AE events for a gene that were adjacent to each other.</P><BR>\n";
$legend .= "</div>\n";

my $rank1 = 0;
my $rank2 = 0;
my $rank3 = 0;
my $rank4 = 0;
foreach my $gene_id (sort {$candidates_ref->{$b}->{adjusted_score} <=> $candidates_ref->{$a}->{adjusted_score}} keys %{$candidates_ref}){
  my $adjusted_score = sprintf("%.2f", $candidates_ref->{$gene_id}->{adjusted_score});
  my $fold_change = sprintf("%.2f", $candidates_ref->{$gene_id}->{fold_change});
  my $ensembl_g_id = $genes_ref->{$gene_id}->{ensembl_g_id};
  my $gene_name = $genes_ref->{$gene_id}->{gene_name};
  my $gene_record_path = "genes/$genes_ref->{$gene_id}->{partition}/$ensembl_g_id".".htm";
  my $gene_record_link = "<A HREF=\"$gene_record_path\">$gene_name</A>";

  my $adjacency_percent_transcript = $genes_ref->{$gene_id}->{adjacency_percent_transcript};
  if ($candidates_ref->{$gene_id}->{ae_event_count} == 0){
    $adjacency_percent_transcript = "N/A";
  }

  #Generate the html row $genes_ref->{$gene_id}->{transcript_count}
  my $row = "<TD CLASS=\"Data4\">$candidates_ref->{$gene_id}->{rank}</TD><TD CLASS=\"Data4\">$adjusted_score</TD><TD CLASS=\"Data1\">$gene_record_link</TD><TD CLASS=\"Data1\">'$genes_ref->{$gene_id}->{gene_type}'</TD><TD CLASS=\"Data2\">$genes_ref->{$gene_id}->{transcript_count}</TD><TD CLASS=\"Data2\">$genes_ref->{$gene_id}->{exon_count}</TD><TD CLASS=\"Data2\">$candidates_ref->{$gene_id}->{main_type}</TD><TD CLASS=\"Data2\">$candidates_ref->{$gene_id}->{direction}</TD><TD CLASS=\"Data2\">$fold_change</TD><TD CLASS=\"Data2\">$candidates_ref->{$gene_id}->{ae_event_count}</TD><TD CLASS=\"Data2\">$candidates_ref->{$gene_id}->{ae_codes}</TD><TD CLASS=\"Data2\">$candidates_ref->{$gene_id}->{seq_name}</TD><TD CLASS=\"Data2\">$adjacency_percent_transcript</TD>\n";

  #Add it to each section if appropriate
  if ($main_list_count <= $top_n){
    $rank1++;
    my $row1 = "<TR>"."<TD CLASS=\"Data4\">$rank1</TD>"."$row"."</TR>";
    $main_list_content .= $row1;
    $main_list_count++;
  }

  if ($ae_list_count <= $top_n && $candidates_ref->{$gene_id}->{ae_event_count} > 0){
    $rank2++;
    my $row2 = "<TR>"."<TD CLASS=\"Data4\">$rank2</TD>"."$row"."</TR>";
    $ae_list_content .= $row2;
    $ae_list_count++;
  }

  if ($gain_list_count <= $top_n && $candidates_ref->{$gene_id}->{direction} eq "Gain"){
    $rank3++;
    my $row3 = "<TR>"."<TD CLASS=\"Data4\">$rank3</TD>"."$row"."</TR>";
    $gain_list_content .= $row3; 
    $gain_list_count++;
  }

  if ($loss_list_count <= $top_n && $candidates_ref->{$gene_id}->{direction} eq "Loss"){
    $rank4++;
    my $row4 = "<TR>"."<TD CLASS=\"Data4\">$rank4</TD>"."$row"."</TR>";
    $loss_list_content .= $row4;
    $loss_list_count++;
  }
}


#Join the content sections together
my $div_count = 0;
my $div_name = '';
my $content = '';
$content .= $summary_content;

#Main list table
$div_count++;
$div_name="box"."$div_count"."h";
$content .= "<P CLASS=\"Indented12LR_s16_Bold\">Candidate Differential AND Alternative Expression Genes</P>\n";
$content .= "<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[$div_name]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>\n";
$content .= "<div id=\"$div_name\">\n";
$content .= "$legend";
$content .= "<TABLE CLASS=\"Data2\">\n";
$content .= $main_list_content;
$content .= "</TABLE><BR>\n</div><BR>\n";

#AE list table
$div_count++;
$div_name="box"."$div_count"."h";
$content .= "<P CLASS=\"Indented12LR_s16_Bold\">Candidate Alternative Expression Genes</P>\n";
$content .= "<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[$div_name]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>\n";
$content .= "<div id=\"$div_name\">\n";
$content .= "$legend";
$content .= "<TABLE CLASS=\"Data2\">\n";
$content .= $ae_list_content;
$content .= "</TABLE><BR>\n</div><BR>\n";

#Gain list table
$div_count++;
$div_name="box"."$div_count"."h";
$content .= "<P CLASS=\"Indented12LR_s16_Bold\">Candidate Differential and Alternative Expression Genes (Gain ONLY)</P>\n";
$content .= "<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[$div_name]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>\n";
$content .= "<div id=\"$div_name\">\n";
$content .= "$legend";
$content .= "<TABLE CLASS=\"Data2\">\n";
$content .= $gain_list_content;
$content .= "</TABLE><BR>\n</div><BR>\n";

#Loss list table
$div_count++;
$div_name="box"."$div_count"."h";
$content .= "<P CLASS=\"Indented12LR_s16_Bold\">Candidate Differential and Alternative Expression Genes (Loss ONLY)</P>\n";
$content .= "<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[$div_name]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>\n";
$content .= "<div id=\"$div_name\">\n";
$content .= "$legend";
$content .= "<TABLE CLASS=\"Data2\">\n";
$content .= $loss_list_content;
$content .= "</TABLE><BR>\n</div><BR>\n";


#Write out the page
my $title = "Summary page for comparison: '$comparison_name' ($comparison_id) - Project: $project";
my $meta_description = "Provides lists of candidate DE and AE genes for comparison '$comparison_name' of project '$project'";
my $meta_keywords = "Differential Expression, Alternative Expression, $comparison_name, $project";

&writePage('-path'=>$web_path, '-title'=>$title, '-content'=>\$content, '-css_path'=>"../ALEXA2.css", '-alexa_home_path'=>"$alexa_home_path", '-alexa_seq_home_path'=>"$alexa_seq_path", '-summary_path'=>"Summary.htm", '-genes_path'=>"genes/index.html", '-search_path'=>"$search_page_url", '-meta_description'=>$meta_description, '-meta_keywords'=>$meta_keywords, '-google_analytics'=>1, '-google_analytics_id'=>$google_analytics_id, '-collapse_div_script'=>'../animatedcollapse.js', '-jquery_min_script'=>'../jquery.min.js', '-div_count'=>$div_count);

my $mem_message = &memoryUsage();
print YELLOW, "\n\n$mem_message\n\n", RESET;
print "\n\n";



exit();

#######################################################################################################################################################################
#Get basic gene info from ALEXA DB                                                                                                                                    #
#######################################################################################################################################################################
sub getBasicGeneInfo{
  my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
  my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};
  my $gene_storable = "$database"."_AllGenes_GeneInfo_NoSeq.storable";
  $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no", '-storable'=>$gene_storable, '-silent'=>"yes");
  $alexa_dbh->disconnect();

  #Initialize a hash to store the feature list for each gene.
  while (my ($gene_id) = each %{$genes_ref}){
    my %fid_list;
    $genes_ref->{$gene_id}->{fid_list} = \%fid_list;
  }

  #Determine which region each gene corresponds to...
  while (my ($gene_id) = each %{$genes_ref}){
    my $chromosome = $genes_ref->{$gene_id}->{chromosome};
    my $chr_start = $genes_ref->{$gene_id}->{chr_start};
    my $chr_end = $genes_ref->{$gene_id}->{chr_end};
    if ($partitions_ref->{$chromosome}){
      my $part_ref = $partitions_ref->{$chromosome}->{partitions};
      foreach my $part (keys %{$part_ref}){
        my $part_start = $part_ref->{$part}->{start};
        my $part_end = $part_ref->{$part}->{end};
        if (($chr_start >= $part_start && $chr_start <= $part_end) && ($chr_end >= $part_start && $chr_end <= $part_end)){
          $genes_ref->{$gene_id}->{partition} = "chr"."$chromosome"."_"."$part_ref->{$part}->{region}";
        }
      }
    }else{
      print RED, "\n\nCould not find matching chromosome ($chromosome) in regions object for gene: $gene_id\n\n", RESET;
      exit();
    }
  }

  return();
}


#######################################################################################################################################################################
#Get the DE values to be considered                                                                                                                                   #
#######################################################################################################################################################################
sub getDeValues{
  my %args = @_;
  my $de_dir = $args{'-de_dir'};
  my $ensembl_version = $args{'-ensembl_version'};
  my $comparison_name = $args{'-comparison_name'};

  #Only get the gene DE values (all other values will come from SI files)
  my $gene_de_dir = "$de_dir"."ENST_v"."$ensembl_version";
  $gene_de_dir = &checkDir('-dir'=>$gene_de_dir, '-clear'=>"no");

  my $gene_de_file = "$gene_de_dir"."$comparison_name"."_Gene_DE_Values_Significant_MTC.txt";

  #Parse through the ranked file (note the gene DE rank) and store critical info in a hash
  my %gene_de;
  open(GENE_DE, "$gene_de_file") || die "\nCould not open file: $gene_de_file\n\n";
  my $r = 0;
  my $header = 1;
  my %columns;
  while(<GENE_DE>){
    chomp($_);
    my @line = split("\t", $_);

    if ($header == 1){
      my $column_count = 0;
      foreach my $column (@line){
        $columns{$column}{column_pos} = $column_count;
        $column_count++;
      }
      $header = 0;
      next();
    }
    $r++;
    my $fid = $line[$columns{'FID'}{column_pos}];
    my $gene_id = $line[$columns{'Gene_ID'}{column_pos}];
    $gene_de{$fid}{gene_id} = $gene_id;
    $gene_de{$fid}{seq_name} = $line[$columns{'Seq_Name'}{column_pos}];
    $gene_de{$fid}{rank} = $r;
    $gene_de{$fid}{libA_norm_coverage} = $line[$columns{'A_Norm'}{column_pos}];
    $gene_de{$fid}{libB_norm_coverage} = $line[$columns{'B_Norm'}{column_pos}];
    $gene_de{$fid}{fold_change} = $line[$columns{'Fold_Change'}{column_pos}];
    $gene_de{$fid}{log2_diff} = $line[$columns{'Log2_Diff'}{column_pos}];

    my @gene_ids = split(" ", $gene_id);
    foreach my $gid (@gene_ids){
      $master_gene_list{$gid}=1;
    }

  }
  close(GENE_DE);

  return(\%gene_de);
}


#######################################################################################################################################################################
#Get the SI values to be considered                                                                                                                                   #
#######################################################################################################################################################################
sub getAeValues{
  my %args = @_;
  my $ae_dir = $args{'-ae_dir'};
  my $ensembl_version = $args{'-ensembl_version'};
  my $comparison_name = $args{'-comparison_name'};

  my %types;
  $types{'1'}{type} = "Transcripts";
  $types{'1'}{subtype} = "Transcript";
  $types{'2'}{type} = "ENST";
  $types{'2'}{subtype} = "ExonRegion";
  $types{'3'}{type} = "Junctions";
  $types{'3'}{subtype} = "Junction";
  $types{'4'}{type} = "Boundaries";
  $types{'4'}{subtype} = "Boundary";
  $types{'5'}{type} = "Introns";
  $types{'5'}{subtype} = "Intron";
  $types{'6'}{type} = "Introns";
  $types{'6'}{subtype} = "ActiveIntronRegion";
  $types{'7'}{type} = "Introns";
  $types{'7'}{subtype} = "SilentIntronRegion";

  my %ae;

  foreach my $c (sort {$a <=> $b} keys %types){
    my $type = $types{$c}{type};
    my $subtype = $types{$c}{subtype};

    #Only get the gene DE values (all other values will come from SI files)
    my $feature_ae_dir = "$ae_dir"."$type"."_v"."$ensembl_version";
    $feature_ae_dir = &checkDir('-dir'=>$feature_ae_dir, '-clear'=>"no");
    my $feature_ae_file = "$feature_ae_dir"."$comparison_name"."_"."$subtype"."_SI_Values_Sorted_Cutoff.txt";

    open(AE, "$feature_ae_file") || die "\nCould not open file: $feature_ae_file\n\n";
    my $r = 0;
    my $header = 1;
    my %columns;
    while(<AE>){
      chomp($_);
      my @line = split("\t", $_);

      if ($header == 1){
        my $column_count = 0;
        foreach my $column (@line){
          $columns{$column}{column_pos} = $column_count;
          $column_count++;
        }
        $header = 0;
        next();
      }
      $r++;
      my $fid = $line[$columns{'FID'}{column_pos}];
      my $gene_id = $line[$columns{'Gene_ID'}{column_pos}];
      $ae{$fid}{rank} = $r;
      $ae{$fid}{libA_norm_coverage} = $line[$columns{'A_SEQ_Norm'}{column_pos}];
      $ae{$fid}{libB_norm_coverage} = $line[$columns{'B_SEQ_Norm'}{column_pos}];
      $ae{$fid}{fold_change} = $line[$columns{'SEQ_Fold_Change'}{column_pos}];
      $ae{$fid}{log2_diff} = $line[$columns{'SEQ_Log2_Diff'}{column_pos}];
      $ae{$fid}{si} = $line[$columns{'SI'}{column_pos}];
      $ae{$fid}{reciprocal} = $line[$columns{'Reciprocal'}{column_pos}];
      $ae{$fid}{reciprocity} = $line[$columns{'Reciprocity'}{column_pos}];
      $ae{$fid}{percent_seq_de} = $line[$columns{'percent_SEQ_Log2_DE'}{column_pos}];

      my @gene_ids = split(" ", $gene_id);
      foreach my $gid (@gene_ids){
        $master_gene_list{$gid}=1;
      }

    }
    close(AE);
  }

  return(\%ae);
}

#######################################################################################################################################################################
#Get the feature annotations                                                                                                                                          #
#######################################################################################################################################################################
sub getFeatures{
  my %args = @_;
  my $annotation_dir = $args{'-annotation_dir'};
  my $junction_seq_size = $args{'-junction_seq_size'};
  my $boundary_seq_size = $args{'-boundary_seq_size'};

  my %a;

  my %types;
  $types{'1'}{type} = "genes";
  $types{'1'}{subtype} = "Gene";
  $types{'1'}{file} = "genes_annotated.txt.gz";
  $types{'2'}{type} = "transcripts";
  $types{'2'}{subtype} = "Transcript";
  $types{'2'}{file} = "transcripts_annotated.txt.gz";
  $types{'3'}{type} = "exonRegions";
  $types{'3'}{subtype} = "ExonRegion";
  $types{'3'}{file} = "exonRegions_annotated.txt.gz";
  $types{'4'}{type} = "exonJunctions";
  $types{'4'}{subtype} = "Junction";
  $types{'4'}{file} = "exonJunctions_"."$junction_seq_size"."mers_annotated.txt.gz";
  $types{'5'}{type} = "exonBoundaries";
  $types{'5'}{subtype} = "Boundary";
  $types{'5'}{file} = "exonBoundaries_"."$boundary_seq_size"."mers_annotated.txt.gz";
  $types{'6'}{type} = "introns";
  $types{'6'}{subtype} = "Intron";
  $types{'6'}{file} = "introns_annotated.txt.gz";
  $types{'7'}{type} = "introns";
  $types{'7'}{subtype} = "ActiveIntronRegion";
  $types{'7'}{file} = "activeIntronRegions.txt.gz";
  $types{'8'}{type} = "introns";
  $types{'8'}{subtype} = "SilentIntronRegion";
  $types{'8'}{file} = "silentIntronRegions.txt.gz";

  foreach my $c (sort {$a <=> $b} keys %types){
    my $type = $types{$c}{type};
    my $subtype = $types{$c}{subtype};
    my $file = $types{$c}{file};

    my $dir = "$annotation_dir"."$type";
    $dir = &checkDir('-dir'=>$dir, '-clear'=>"no");

    my $file_path = "$dir"."$file";

    open(ANN, "zcat $file_path |") || die "\nCould not open annotation file: $file_path\n\n";
    my $header = 1;
    my %columns;
    while(<ANN>){
      chomp($_);
      my @line = split("\t", $_);

      if ($header == 1){
        my $column_count = 0;
        foreach my $column (@line){
          $columns{$column}{column_pos} = $column_count;
          $column_count++;
        }
        $header = 0;
        next();
      }
      my $fid = $line[$columns{'FID'}{column_pos}];

      #Get gene ID list, remember that intron features are often associated with multiple genes
      my $gene_id_list;
      if ($columns{'Gene_ID'}){
        $gene_id_list = $line[$columns{'Gene_ID'}{column_pos}];
      }elsif($columns{'Gene_ID_List'}){
        $gene_id_list = $line[$columns{'Gene_ID_List'}{column_pos}];
      }else{
        print RED, "\n\nCould not find gene id column\n\n", RESET;
      }
      my @gene_ids = split(" ", $gene_id_list);

      #Only import annotations for genes that have at least one DE or AE feature
      my $target_gene = 0;
      foreach my $gid (@gene_ids){
        if ($master_gene_list{$gid}){
          $target_gene = 1;
        }
      }
      unless($target_gene){
        next();
      }

      #Get coordinates to allow ordering of features
      my $chr_start;
      my $chr_end;
      if ($columns{'Unit1_start_chr'}){
        $chr_start = $line[$columns{'Unit1_start_chr'}{column_pos}];
      }else{
        $chr_start = 0;
      }
      if ($columns{'Unit2_end_chr'}){
        $chr_end = $line[$columns{'Unit2_end_chr'}{column_pos}];
      }elsif($columns{'Unit1_end_chr'}){
        $chr_end = $line[$columns{'Unit1_end_chr'}{column_pos}];
      }else{
        $chr_end = 0;
      }


      #Add this feature of to a list of ALL features for the corresponding gene(s) (fid_list)
      foreach my $gene_id (@gene_ids){
        unless ($genes_ref->{$gene_id}){
          print RED, "\n\nFound gene: $gene_id that is not defined in the genes_ref\n\n$_", RESET;
          exit();
        }
        my $fid_list_ref = $genes_ref->{$gene_id}->{fid_list};
        $fid_list_ref->{$fid}->{chr_start} = $chr_start;
        $fid_list_ref->{$fid}->{chr_end} = $chr_end;
        $fid_list_ref->{$fid}->{subtype} = $subtype;
      }

      #Gene specific columns
      if ($type eq "genes"){
        my $gene_id = $gene_ids[0];

        unless ($genes_ref->{$gene_id}){
          print RED, "\n\nGene id: $gene_id should already be defined\n\n", RESET;
          exit();
        }
        unless ($columns{'Description'} && $columns{'Gene_Type'} && $columns{'Gene_Evidence'} && $columns{'Transcript_Count'} && $columns{'Exon_Count'} && $columns{'Exon_Content_Count'}){
          print RED, "\n\nMissing column in gene annotation file\n\n", RESET;
          exit();
        }
        $genes_ref->{$gene_id}->{description} = $line[$columns{'Description'}{column_pos}];
        $genes_ref->{$gene_id}->{gene_type} = $line[$columns{'Gene_Type'}{column_pos}];
        $genes_ref->{$gene_id}->{evidence} = $line[$columns{'Gene_Evidence'}{column_pos}];
        $genes_ref->{$gene_id}->{transcript_count} = $line[$columns{'Transcript_Count'}{column_pos}];
        $genes_ref->{$gene_id}->{exon_count} = $line[$columns{'Exon_Count'}{column_pos}];
        $genes_ref->{$gene_id}->{exon_content_count} = $line[$columns{'Exon_Content_Count'}{column_pos}];
      }

      #Unless this FID is defined in the SI or DE objects, dont bother storing it:
      unless ($de_ref->{$fid} || $ae_ref->{$fid}){
        next();
      }

      $a{$fid}{subtype} = $subtype;
      $a{$fid}{gene_ids} = \@gene_ids;
      $a{$fid}{ensembl_g_id} = $line[$columns{'EnsEMBL_Gene_ID'}{column_pos}];
      $a{$fid}{seq_name} = $line[$columns{'Seq_Name'}{column_pos}];

      #Junction specific columns
      if ($type eq "exonJunctions"){
         $a{$fid}{exons_skipped} = $line[$columns{'Exons_Skipped'}{column_pos}];
      }
    }
    close(ANN);
  }

  #Add an order to all features within each gene based on coordinates position and according to strand
  #Generate several orderings:
  #1st ordering is 'all' (all feature types) - sort on chr_start  - sort of chr_start only (features with the same start coord will be ties)
  #
  #3rd ordering is 'exonic' (only exon and exon junction features considered) - sort of chr_start only (features with the same start coord will be ties)
  #4th ordering is 'intronic' (only exon boundaries and intron elements) - sort of chr_start only (features with the same start coord will be ties)

  while (my ($gene_id) = each %{$genes_ref}){

    my $fid_list_ref = $genes_ref->{$gene_id}->{fid_list};
    my $strand = $genes_ref->{$gene_id}->{chr_strand};
    if ($strand eq "1"){
      $strand = "+";
    }else{
      $strand = "-";
    }

    #1st ordering is 'all' (all feature types) - sort on chr_start AND chr_end
    my $order = 0;
    my $current_coord = 0;
    if ($strand eq "+"){
      foreach my $fid (sort {$fid_list_ref->{$a}->{chr_start} <=> $fid_list_ref->{$b}->{chr_start}} keys %{$fid_list_ref}){
        unless ($fid_list_ref->{$fid}->{chr_start} == $current_coord){
          $order++;
        }
        $fid_list_ref->{$fid}->{order_all} = $order;
        $current_coord = $fid_list_ref->{$fid}->{chr_start};
      }
    }else{
      foreach my $fid (sort {$fid_list_ref->{$b}->{chr_end} <=> $fid_list_ref->{$a}->{chr_end}} keys %{$fid_list_ref}){
        unless ($fid_list_ref->{$fid}->{chr_end} == $current_coord){
          $order++;
        }
        $fid_list_ref->{$fid}->{order_all} = $order;
        $current_coord = $fid_list_ref->{$fid}->{chr_end};
      }
    }

    #3rd ordering is 'transcript' (exon AND junction features considered) - sort of chr_start only (features with the same start coord will be ties)
    $order = 0;
    $current_coord = 0;
    if ($strand eq "+"){
      foreach my $fid (sort {$fid_list_ref->{$a}->{chr_start} <=> $fid_list_ref->{$b}->{chr_start}} keys %{$fid_list_ref}){
        my $subtype = $fid_list_ref->{$fid}->{subtype};
        unless ($subtype eq "ExonRegion" || $subtype eq "Junction"){
          $fid_list_ref->{$fid}->{order_transcript} = 0;
          next();
        }
        unless ($fid_list_ref->{$fid}->{chr_start} == $current_coord){
          $order++;
        }
        $fid_list_ref->{$fid}->{order_transcript} = $order;
        $current_coord = $fid_list_ref->{$fid}->{chr_start};
      }
    }else{
      foreach my $fid (sort {$fid_list_ref->{$b}->{chr_end} <=> $fid_list_ref->{$a}->{chr_end}} keys %{$fid_list_ref}){
        my $subtype = $fid_list_ref->{$fid}->{subtype};
        unless ($subtype eq "ExonRegion" || $subtype eq "Junction"){
          $fid_list_ref->{$fid}->{order_transcript} = 0;
          next();
        }
        unless ($fid_list_ref->{$fid}->{chr_end} == $current_coord){
          $order++;
        }
        $fid_list_ref->{$fid}->{order_transcript} = $order;
        $current_coord = $fid_list_ref->{$fid}->{chr_end};
      }
    }


    #3rd ordering is 'exonic' (only exon features considered) - sort of chr_start only (features with the same start coord will be ties)
    $order = 0;
    $current_coord = 0;
    if ($strand eq "+"){
      foreach my $fid (sort {$fid_list_ref->{$a}->{chr_start} <=> $fid_list_ref->{$b}->{chr_start}} keys %{$fid_list_ref}){
        my $subtype = $fid_list_ref->{$fid}->{subtype};
        unless ($subtype eq "ExonRegion"){
          $fid_list_ref->{$fid}->{order_exonic} = 0;
          next();
        }
        unless ($fid_list_ref->{$fid}->{chr_start} == $current_coord){
          $order++;
        }
        $fid_list_ref->{$fid}->{order_exonic} = $order;
        $current_coord = $fid_list_ref->{$fid}->{chr_start};
      }
    }else{
      foreach my $fid (sort {$fid_list_ref->{$b}->{chr_end} <=> $fid_list_ref->{$a}->{chr_end}} keys %{$fid_list_ref}){
        my $subtype = $fid_list_ref->{$fid}->{subtype};
        unless ($subtype eq "ExonRegion" || $subtype eq "Junction"){
          $fid_list_ref->{$fid}->{order_exonic} = 0;
          next();
        }
        unless ($fid_list_ref->{$fid}->{chr_end} == $current_coord){
          $order++;
        }
        $fid_list_ref->{$fid}->{order_exonic} = $order;
        $current_coord = $fid_list_ref->{$fid}->{chr_end};
      }
    }

    #4th ordering is 'intronic' (only exon boundaries and intron elements) - sort of chr_start only (features with the same start coord will be ties)
    $order = 0;
    $current_coord = 0;
    if ($strand eq "+"){
      foreach my $fid (sort {$fid_list_ref->{$a}->{chr_start} <=> $fid_list_ref->{$b}->{chr_start}} keys %{$fid_list_ref}){
        my $subtype = $fid_list_ref->{$fid}->{subtype};
        unless ($subtype eq "Boundary" || $subtype eq "Intron" || $subtype eq "ActiveIntronRegion" || $subtype eq "SilentIntronRegion"){
          $fid_list_ref->{$fid}->{order_intronic} = 0;
          next();
        }
        unless ($fid_list_ref->{$fid}->{chr_start} == $current_coord){
          $order++;
        }

        $fid_list_ref->{$fid}->{order_intronic} = $order;
        $current_coord = $fid_list_ref->{$fid}->{chr_start};
      }
    }else{
      foreach my $fid (sort {$fid_list_ref->{$b}->{chr_end} <=> $fid_list_ref->{$a}->{chr_end}} keys %{$fid_list_ref}){
        my $subtype = $fid_list_ref->{$fid}->{subtype};
        unless ($subtype eq "Boundary" || $subtype eq "Intron" || $subtype eq "ActiveIntronRegion" || $subtype eq "SilentIntronRegion"){
          $fid_list_ref->{$fid}->{order_intronic} = 0;
          next();
        }
        unless ($fid_list_ref->{$fid}->{chr_end} == $current_coord){
          $order++;
        }

        $fid_list_ref->{$fid}->{order_intronic} = $order;
        $current_coord = $fid_list_ref->{$fid}->{chr_end};
      }
    }

  }

  return(\%a);
}


#######################################################################################################################################################################
#Examine the adjacency of AE events to sub-classify events and assign bonus scores for events with multiple supporting features                                       #
#######################################################################################################################################################################
sub examineAdjacency{

  #For each gene and for each set of ordering values ('all', 'transcript', 'exonic' and 'intronic') retrieve the order values for AE events of that gene
  #The 'transcript' adjacency considers both exons and junctions, wheras 'exonic' considers only exons
  #Examinine these ordering value all look for evidence of adjacent AE events
  #Adjacent AE events should have either the same order value or order+1
  #Count the number of adjacencies for each gene's list of AE events

  foreach my $gene_id (keys %{$genes_ref}){
    my $fid_list_ref = $genes_ref->{$gene_id}->{fid_list};
    my $strand = $genes_ref->{$gene_id}->{chr_strand};
    my $gene_name = $genes_ref->{$gene_id}->{gene_name};

    #1.) 'All' ordering 
    my @orders_ae_all;
    my $adjacency_all = 0;
    my $adjacency_percent_all = 0;
    foreach my $fid (sort {$fid_list_ref->{$a}->{order_all} <=> $fid_list_ref->{$b}->{order_all}} keys %{$fid_list_ref}){
      if (($ae_ref->{$fid}) && ($fid_list_ref->{$fid}->{order_all})){
        push(@orders_ae_all, $fid_list_ref->{$fid}->{order_all});
      }
    }
    if (scalar(@orders_ae_all) > 1){
      $adjacency_all = &calculateAdjacency('-orders'=>\@orders_ae_all);
      $adjacency_percent_all = sprintf("%.2f", ($adjacency_all/(scalar(@orders_ae_all)-1))*100);
      #print YELLOW, "\n\nGene: $gene_name\tStrand: $strand", RESET;
      #print YELLOW, "\nAdjacency_ALL = $adjacency_all ($adjacency_percent_all%) for array: @orders_ae_all", RESET;
    }


    #2.) 'transcript' ordering 
    my @orders_ae_transcript;
    my $adjacency_transcript = 0;
    my $adjacency_percent_transcript = 0;
    foreach my $fid (sort {$fid_list_ref->{$a}->{order_transcript} <=> $fid_list_ref->{$b}->{order_transcript}} keys %{$fid_list_ref}){
      if (($ae_ref->{$fid}) && ($fid_list_ref->{$fid}->{order_transcript})){
        push(@orders_ae_transcript, $fid_list_ref->{$fid}->{order_transcript});
      }
    }
    if (scalar(@orders_ae_transcript) > 1){
      $adjacency_transcript = &calculateAdjacency('-orders'=>\@orders_ae_transcript);
      $adjacency_percent_transcript = sprintf("%.2f", ($adjacency_transcript/(scalar(@orders_ae_transcript)-1))*100);
      #print YELLOW, "\nAdjacency_TRANSCRIPT = $adjacency_transcript ($adjacency_percent_transcript%) for array: @orders_ae_transcript", RESET;
    }


    #3.) 'exonic' ordering 
    my @orders_ae_exonic;
    my $adjacency_exonic = 0;
    my $adjacency_percent_exonic = 0;
    foreach my $fid (sort {$fid_list_ref->{$a}->{order_exonic} <=> $fid_list_ref->{$b}->{order_exonic}} keys %{$fid_list_ref}){
      if (($ae_ref->{$fid}) && ($fid_list_ref->{$fid}->{order_exonic})){
        push(@orders_ae_exonic, $fid_list_ref->{$fid}->{order_exonic});
      }
    }
    if (scalar(@orders_ae_exonic) > 1){
      $adjacency_exonic = &calculateAdjacency('-orders'=>\@orders_ae_exonic);
      $adjacency_percent_exonic = sprintf("%.2f", ($adjacency_exonic/(scalar(@orders_ae_exonic)-1))*100);
      #print YELLOW, "\nAdjacency_EXONIC = $adjacency_exonic ($adjacency_percent_exonic%) for array: @orders_ae_exonic", RESET;
    }

    #4.) 'intronic' ordering 
    my @orders_ae_intronic;
    my $adjacency_intronic = 0;
    my $adjacency_percent_intronic = 0;
    foreach my $fid (sort {$fid_list_ref->{$a}->{order_intronic} <=> $fid_list_ref->{$b}->{order_intronic}} keys %{$fid_list_ref}){
      if (($ae_ref->{$fid}) && ($fid_list_ref->{$fid}->{order_intronic})){
        push(@orders_ae_intronic, $fid_list_ref->{$fid}->{order_intronic});
      }
    }
    if (scalar(@orders_ae_intronic) > 1){
      $adjacency_intronic = &calculateAdjacency('-orders'=>\@orders_ae_intronic);
      $adjacency_percent_intronic = sprintf("%.2f", ($adjacency_intronic/(scalar(@orders_ae_intronic)-1))*100);
      #print YELLOW, "\nAdjacency_INTRONIC = $adjacency_intronic ($adjacency_percent_intronic%) for array: @orders_ae_intronic", RESET;
    }

    #Store the adjacency values for each gene
    $genes_ref->{$gene_id}->{adjacency_all} = $adjacency_all;
    $genes_ref->{$gene_id}->{adjacency_transcript} = $adjacency_transcript;
    $genes_ref->{$gene_id}->{adjacency_exonic} = $adjacency_exonic;
    $genes_ref->{$gene_id}->{adjacency_intronic} = $adjacency_intronic;
    $genes_ref->{$gene_id}->{adjacency_percent_all} = $adjacency_percent_all;
    $genes_ref->{$gene_id}->{adjacency_percent_transcript} = $adjacency_percent_transcript;
    $genes_ref->{$gene_id}->{adjacency_percent_exonic} = $adjacency_percent_exonic;
    $genes_ref->{$gene_id}->{adjacency_percent_intronic} = $adjacency_percent_intronic;
  
  
  }


  return();
}


#######################################################################################################################################################################
#Calculate the numerical adjacency given an input array or sorted postion/order values                                                                                #
#######################################################################################################################################################################
sub calculateAdjacency{
  my %args = @_;
  my @orders = @{$args{'-orders'}};
  my $adjacency = 0;

  my $n = scalar(@orders);
  for (my $i = 1; $i < $n; $i++){
    my $diff = $orders[$i] - $orders[$i-1];
    if ($diff == 0 || $diff == 1){
      $adjacency++;
    }
  } 

  return($adjacency);
}


#######################################################################################################################################################################
#Build the ranked candidate list combining both DE and AE events into one grand ranked list                                                                           #
#######################################################################################################################################################################
sub buildCandidateList{
  my %args = @_;
  my $feature_ref = $args{'-feature_ref'};
  my $de_ref = $args{'-de_ref'};
  my $ae_ref = $args{'-ae_ref'};

  my %candidates;

  #Candidate list will be keyed on gene ID (only one entry allowed per gene)

  #Grab gene DE events first
  while (my ($fid) = each %{$de_ref}){
    my $gene_id = $de_ref->{$fid}->{gene_id};
    $candidates{$gene_id}{fid} = $fid;
    $candidates{$gene_id}{ae_event_count} = 0;
    $candidates{$gene_id}{diff_abs} = abs($de_ref->{$fid}->{log2_diff});
    $candidates{$gene_id}{score} = abs($de_ref->{$fid}->{log2_diff});
    $candidates{$gene_id}{bonus} = 0;
    $candidates{$gene_id}{diff} = $de_ref->{$fid}->{log2_diff};
    $candidates{$gene_id}{main_type} = "DE";
    $candidates{$gene_id}{ensembl_g_id} = $feature_ref->{$fid}->{ensembl_g_id};
    $candidates{$gene_id}{gene_name} = $genes_ref->{$gene_id}->{gene_name};
    $candidates{$gene_id}{seq_name} = $feature_ref->{$fid}->{seq_name};
    $candidates{$gene_id}{reciprocity} = "N/A";
    $candidates{$gene_id}{percent_seq_de} = "N/A";
    $candidates{$gene_id}{direction} = "Loss";
    $candidates{$gene_id}{exons_skipped} = "N/A";
    if ($de_ref->{$fid}->{log2_diff} > 0){
       $candidates{$gene_id}{direction} = "Gain";
    }
    #Calculate fold change
    my $fc = 2**(abs($de_ref->{$fid}->{log2_diff}));
    if ($de_ref->{$fid}->{log2_diff} < 0){
      $fc = $fc*-1;
    }
    $candidates{$gene_id}{fold_change} = $fc;

    my %ae_fid_list;
    $candidates{$gene_id}{ae_fid_list} = \%ae_fid_list;
  }


  #Now AE events.  If the gene is already added but AE diff for at least one feature is greater, then replace
  while (my ($fid) = each %{$ae_ref}){
    my @gene_ids = @{$feature_ref->{$fid}->{gene_ids}};
    my $subtype = $feature_ref->{$fid}->{subtype};
    my $si = $ae_ref->{$fid}->{si};
    my $si_abs = abs($ae_ref->{$fid}->{si});
    my $seq_de = $ae_ref->{$fid}->{log2_diff};
    my $seq_de_abs = abs($ae_ref->{$fid}->{log2_diff});

    my $diff = $si;
    my $diff_abs = $si_abs;
    my $score = $si_abs;

    #ALEXA-Seq 'Score' will be based on some combination of values: Gene DE, Feature DE, SI, Reciprocal?, multiple events per gene?, event type? 
    

    #If the event is reciprocal, make a note to give a bonus to this gene
    #Since the gene may have several SI events, also give a bonus to the specific feature to increase the chance that it will be the 'top' feature for this gene
    my $reciprocity_abs = abs($ae_ref->{$fid}->{reciprocity});
    my $bonus = 0;
    if ($ae_ref->{$fid}->{reciprocal} == 1){
      $bonus = $si_abs*$bonus_ratio;
      $score += $bonus;
    }

    #If the AE event has a high percent_SEQ_DE, give a bonus to this gene
    if ($ae_ref->{$fid}->{percent_seq_de} >= 80){
      $bonus += $si_abs*$bonus_ratio;
      $score += $si_abs*$bonus_ratio;
    }

    #Calculate fold change
    my $fc = 2**($seq_de_abs);
    if ($seq_de < 0){
      $fc = $fc*-1;
    }

    #If this feature is a junction make note of the number of exons skipped
    my $exons_skipped = "N/A";
    if ($subtype eq "Junction"){
      $exons_skipped = $feature_ref->{$fid}->{exons_skipped};
    }

    #Remember that each feature can be associated with multiple gene (especially introns)
    foreach my $gene_id (@gene_ids){
      if(defined($candidates{$gene_id})){
        #Already found at least one AE feature for this gene
        unless ($subtype eq "Transcript" || $subtype eq "Gene"){
          $candidates{$gene_id}{ae_event_count}++;
        }
        my $ae_fid_list_ref = $candidates{$gene_id}{ae_fid_list};
        $ae_fid_list_ref->{$fid}=1;
        if($score > $candidates{$gene_id}{score}){
          $candidates{$gene_id}{fid} = $fid;
          $candidates{$gene_id}{diff_abs} = $diff_abs;
          $candidates{$gene_id}{score} = $si_abs;
          $candidates{$gene_id}{bonus} = $bonus;
          $candidates{$gene_id}{diff} = $diff;
          $candidates{$gene_id}{main_type} = "AE";
          $candidates{$gene_id}{seq_name} = $feature_ref->{$fid}->{seq_name};
          $candidates{$gene_id}{reciprocity} = $ae_ref->{$fid}->{reciprocity};
          $candidates{$gene_id}{percent_seq_de} = $ae_ref->{$fid}->{percent_seq_de};
          $candidates{$gene_id}{direction} = "Loss";
          if ($seq_de > 0){
            $candidates{$gene_id}{direction} = "Gain";
          }
          $candidates{$gene_id}{exons_skipped} = $exons_skipped;
          $candidates{$gene_id}{fold_change} = $fc;
        }
      }else{
        #First AE feature for this gene
        $candidates{$gene_id}{ae_event_count} = 0;
        unless ($subtype eq "Transcript" || $subtype eq "Gene"){
          $candidates{$gene_id}{ae_event_count} = 1;
        }
        my %ae_fid_list;
        $ae_fid_list{$fid}=1;
        $candidates{$gene_id}{ae_fid_list} = \%ae_fid_list;

        $candidates{$gene_id}{fid} = $fid;
        $candidates{$gene_id}{diff_abs} = $diff_abs;
        $candidates{$gene_id}{score} = $si_abs;
        $candidates{$gene_id}{bonus} = $bonus;
        $candidates{$gene_id}{diff} = $diff;
        $candidates{$gene_id}{main_type} = "AE";
        $candidates{$gene_id}{ensembl_g_id} = $feature_ref->{$fid}->{ensembl_g_id};
        $candidates{$gene_id}{gene_name} = $genes_ref->{$gene_id}->{gene_name};
        $candidates{$gene_id}{seq_name} = $feature_ref->{$fid}->{seq_name};
        $candidates{$gene_id}{reciprocity} = $ae_ref->{$fid}->{reciprocity};
        $candidates{$gene_id}{percent_seq_de} = $ae_ref->{$fid}->{percent_seq_de};
        $candidates{$gene_id}{direction} = "Loss";
        if ($seq_de > 0){
          $candidates{$gene_id}{direction} = "Gain";
        }
        $candidates{$gene_id}{exons_skipped} = $exons_skipped;
        $candidates{$gene_id}{fold_change} = $fc;
      }
    }
  }

  #If there are multiple AE events supporting the same event ... give a bonus
  #These multiple AE events supporting the 'same event' should be adjacent to each other for the bonus to be applied
  while (my ($gene_id) = each %candidates){
    my $adjacency_all = $genes_ref->{$gene_id}->{adjacency_all};
    my $adjacency_transcript = $genes_ref->{$gene_id}->{adjacency_transcript};
    my $adjacency_exonic = $genes_ref->{$gene_id}->{adjacency_exonic};
    my $adjacency_intronic = $genes_ref->{$gene_id}->{adjacency_intronic};
    my $adjacency_percent_all = $genes_ref->{$gene_id}->{adjacency_percent_all};
    my $adjacency_percent_transcript = $genes_ref->{$gene_id}->{adjacency_percent_transcript};
    my $adjacency_percent_exonic = $genes_ref->{$gene_id}->{adjacency_percent_exonic};
    my $adjacency_percent_intronic = $genes_ref->{$gene_id}->{adjacency_percent_intronic};

    if ($adjacency_percent_exonic >= 75){
      my $bonus = ($candidates{$gene_id}{score})*$bonus_ratio;
      $candidates{$gene_id}{bonus} += $bonus;
    }
  }

  #Count the number of each type of ae events for each candidate gene
  while (my ($gene_id) = each %candidates){
    $candidates{$gene_id}{ae_exons} = 0;
    $candidates{$gene_id}{ae_exons_included_j} = 0;
    $candidates{$gene_id}{ae_exons_skipped_j} = 0;
    $candidates{$gene_id}{ae_alternative_boundaries} = 0;
    $candidates{$gene_id}{ae_introns_retained} = 0;
    $candidates{$gene_id}{ae_cryptic_exons} = 0;

    my $ae_fid_list_ref = $candidates{$gene_id}{ae_fid_list};
    foreach my $fid (keys %{$ae_fid_list_ref}){
      my $subtype = $feature_ref->{$fid}->{subtype};

      #ExonRegion
      if ($subtype eq "ExonRegion"){
        $candidates{$gene_id}{ae_exons}++;
      }

      #Junction
      if ($subtype eq "Junction"){
        my $exons_skipped = $feature_ref->{$fid}->{exons_skipped};
        if ($exons_skipped == 0){
          $candidates{$gene_id}{ae_exons_included_j}++;
        }else{
          $candidates{$gene_id}{ae_exons_skipped_j}++;
        }
      }
      #Boundary
      if ($subtype eq "Boundary"){
        $candidates{$gene_id}{ae_alternative_boundaries}++;
      }

      #Intron
      if ($subtype eq "Intron"){
        $candidates{$gene_id}{ae_introns_retained}++;
      }

      #ActiveIntronRegion
      if ($subtype eq "ActiveIntronRegion"){
        #Refine this step to use more information to differentiate cryptic exons from retained introns
        #If the whole intron is retained, then don't count this as a cryptic exon
        my $air_name = $feature_ref->{$fid}->{seq_name};
        if ($air_name =~ /I(\d+)\_AR\d+/){
          my $intron1_count = $1;
          my $intron_retained = 0;
          foreach my $fid (keys %{$ae_fid_list_ref}){
            my $subtype = $feature_ref->{$fid}->{subtype};
            unless($subtype eq "Intron"){
              next();
            }
            my $i_name = $feature_ref->{$fid}->{seq_name};
            if ($i_name =~ /I(\d+)/){
              my $intron2_count = $1;
              if ($intron2_count == $intron1_count){
                $intron_retained = 1;
              }
            }
          }
          unless ($intron_retained){
            $candidates{$gene_id}{ae_cryptic_exons}++;
          }
        }
      }
    }
  }

  #Assign alternative expression codes for each gene based on the types of AE events observed
  while (my ($gene_id) = each %candidates){

    $candidates{$gene_id}{ae_codes} = "N/A";
    if ($candidates{$gene_id}{ae_event_count} > 0 ){
      my $codes= '';
      if ($candidates{$gene_id}{ae_exons} > 0 || $candidates{$gene_id}{ae_exons_included_j} > 0){
        $codes .= "EU ";
      }
      if ($candidates{$gene_id}{ae_exons_skipped_j} > 0){
        $codes .= "ES ";
      }
      if ($candidates{$gene_id}{ae_alternative_boundaries} > 0){
        $codes .= "AB ";
      }
      if ($candidates{$gene_id}{ae_introns_retained} > 0){
        $codes .= "IR ";
      }
      if ($candidates{$gene_id}{ae_cryptic_exons} > 0){
        $codes .= "CE ";
      }
      $candidates{$gene_id}{ae_codes} = $codes;
    }
  }


  return(\%candidates);
}

#######################################################################################################################################################################
#Generate a string of all AE feature names (ordered by chr coords)                                                                                                    #
#######################################################################################################################################################################
sub generateFeatureNameStrings{
  
  #Go through the features of each gene, in order, and if the feature is AE, add it to the list of feature names - ALL FEATURE TYPES
  while (my ($gene_id) = each %{$genes_ref}){
    my $fid_list_ref = $genes_ref->{$gene_id}->{fid_list};

    my @feature_string_all;
    foreach my $fid (sort {$fid_list_ref->{$a}->{order_all} <=> $fid_list_ref->{$b}->{order_all}} keys %{$fid_list_ref}){
      if ($ae_ref->{$fid}){
        my $subtype = $feature_ref->{$fid}->{subtype};
        unless ($subtype eq "Transcript" || $subtype eq "Gene"){
          push(@feature_string_all, $feature_ref->{$fid}->{seq_name});
        }
      }
    }
    $genes_ref->{$gene_id}->{feature_string_all} = \@feature_string_all;
  }

  #Go through the features of each gene, in order, and if the feature is AE, add it to the list of feature names - TRANSCRIPT FEATURE TYPES
  while (my ($gene_id) = each %{$genes_ref}){
    my $fid_list_ref = $genes_ref->{$gene_id}->{fid_list};

    my @feature_string_transcript;
    foreach my $fid (sort {$fid_list_ref->{$a}->{order_transcript} <=> $fid_list_ref->{$b}->{order_transcript}} keys %{$fid_list_ref}){
      if ($ae_ref->{$fid}){
        my $subtype = $feature_ref->{$fid}->{subtype};
        unless ($subtype eq "ExonRegion" || $subtype eq "Junction"){
          next();
        }
        push(@feature_string_transcript, $feature_ref->{$fid}->{seq_name});
      }
    }
    $genes_ref->{$gene_id}->{feature_string_transcript} = \@feature_string_transcript;
  }


  #Go through the features of each gene, in order, and if the feature is AE, add it to the list of feature names - EXONIC FEATURE TYPES
  while (my ($gene_id) = each %{$genes_ref}){
    my $fid_list_ref = $genes_ref->{$gene_id}->{fid_list};

    my @feature_string_exonic;
    foreach my $fid (sort {$fid_list_ref->{$a}->{order_exonic} <=> $fid_list_ref->{$b}->{order_exonic}} keys %{$fid_list_ref}){
      if ($ae_ref->{$fid}){
        my $subtype = $feature_ref->{$fid}->{subtype};
        unless ($subtype eq "ExonRegion"){
          next();
        }
        push(@feature_string_exonic, $feature_ref->{$fid}->{seq_name});
      }
    }
    $genes_ref->{$gene_id}->{feature_string_exonic} = \@feature_string_exonic;
  }


  #Go through the features of each gene, in order, and if the feature is AE, add it to the list of feature names - INTRONIC FEATURE TYPES
  while (my ($gene_id) = each %{$genes_ref}){
    my $fid_list_ref = $genes_ref->{$gene_id}->{fid_list};

    my @feature_string_intronic;
    foreach my $fid (sort {$fid_list_ref->{$a}->{order_intronic} <=> $fid_list_ref->{$b}->{order_intronic}} keys %{$fid_list_ref}){
      if ($ae_ref->{$fid}){
        my $subtype = $feature_ref->{$fid}->{subtype};
        unless ($subtype eq "Boundary" || $subtype eq "Intron" || $subtype eq "ActiveIntronRegion" || $subtype eq "SilentIntronRegion"){
          next();
        }
        push(@feature_string_intronic, $feature_ref->{$fid}->{seq_name});
      }
    }
    $genes_ref->{$gene_id}->{feature_string_intronic} = \@feature_string_intronic;
  }

  return();
}



