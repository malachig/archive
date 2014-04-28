#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

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
my $peptide_size = '';
my $project = '';
my $comparison_id = '';
my $comparison_name = '';
my $ensembl_version = '';
my $analysis_dir = '';
my $annotation_dir = '';
my $junction_seq_size = '';
my $partition_file = '';
my $alexa_home_path = '';
my $alexa_seq_path = '';
my $search_page_url = '';
my $google_analytics_id = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'peptide_size=i'=>\$peptide_size,
            'project=s'=>\$project, 'comparison_id=s'=>\$comparison_id, 'comparison_name=s'=>\$comparison_name, 'ensembl_version=s'=>\$ensembl_version, 
            'analysis_dir=s'=>\$analysis_dir, 'annotation_dir=s'=>\$annotation_dir, 'junction_seq_size=i'=>\$junction_seq_size, 'partition_file=s'=>\$partition_file, 
            'alexa_home_path=s'=>\$alexa_home_path, 'alexa_seq_path=s'=>\$alexa_seq_path, 'search_page_url=s'=>\$search_page_url, 'google_analytics_id=s'=>\$google_analytics_id);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the desired peptide size using:  --peptide_size", RESET;
print GREEN, "\n\tSpecify the project name using:  --project", RESET;
print GREEN, "\n\tSpecify the comparison id used for stats files using:  --comparison_id", RESET;
print GREEN, "\n\tSpecify the comparison name used for stats files using:  --comparison_name", RESET;
print GREEN, "\n\tSpecify the ensembl version using:  --ensembl_version", RESET;
print GREEN, "\n\tSpecify the root path to the analysis dir using:  --analysis_dir", RESET;
print GREEN, "\n\tSpecify the root path to the annotation dir used for the analysis using:  --annotation_dir", RESET;
print GREEN, "\n\tSpecify the junction database sequence length using:  --junction_seq_size", RESET;
print GREEN, "\n\tSpecify the genome partition file using:  --partition_file", RESET;
print GREEN, "\n\tSpecify the URL to your ALEXA home page using: --alexa_home_path", RESET;
print GREEN, "\n\tSpecify the URL to your ALEXA-Seq results page using: --alexa_seq_path", RESET;
print GREEN, "\n\tSpecify the URL to your Xapian-Omega search page using:  --search_page_url", RESET;
print GREEN, "\n\tSpecify your Google Analytics ID using: --google_analytics_id", RESET;

print GREEN, "\n\nExample: generatePeptideLists.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --peptide_size=15  --project=Breast  --comparison_id=HS1188_vs_HS1187  --comparison_name=Myo_Epi_vs_Lum_Epi  --ensembl_version=53  --analysis_dir=/projects/malachig/alexa_seq/  --annotation_dir=/projects/alexa/sequence_databases/hs_53_36o/  --junction_seq_size=62  --partition_file=/projects/alexa/sequence_databases/hs_53_36o/Regions_250_Genes.txt  --alexa_home_path=http://www.alexaplatform.org/index.htm  --alexa_seq_path=http://www.alexaplatform.org/alexa_seq/results.htm  --search_page_url=http://www.bcgsc.ca/xapian-search/omega  --google_analytics_id=UA-xxxxxx-x\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $peptide_size && $project && $comparison_id && $comparison_name && $ensembl_version && $analysis_dir && $annotation_dir && $junction_seq_size && $partition_file && $alexa_home_path && $alexa_seq_path && $search_page_url && $google_analytics_id){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Load the specified ensembl API and BioPerl
&loadEnsemblApi('-api'=>$ensembl_version);
require Bio::Perl;

#Limit the number of features stored for each category (exons, known junctions, novel junctions) for each gene
my $feature_limit = 1;

#Get the library names
my $libA_name = '';
my $libB_name = '';
if ($comparison_name =~ /(.*)\_vs\_(.*)/){
  $libA_name = $2;
  $libB_name = $1;
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
my $summary_file_gain = "$web_dir"."data/"."$comparison_name"."_peptides_GAIN.txt";
my $summary_file_loss = "$web_dir"."data/"."$comparison_name"."_peptides_LOSS.txt";
my $web_path ="$web_dir"."$comparison_id"."_peptides.htm";

#Import the genome partitions
my $partitions_ref = &getPartitions('-file'=>$partition_file);

#Get DE features:
print BLUE, "\n\nGetting Gene DE values:", RESET;
my $de_ref = &getDeValues('-de_dir'=>$de_dir, '-ensembl_version'=>$ensembl_version, '-comparison_name'=>$comparison_name);
my $de_count = keys %{$de_ref};
print BLUE, "\nStored $de_count Feature DE records", RESET;

#Which of these DE features is also AE??
my $ae_count = &getAeFeatures('-si_dir'=>$si_dir, '-ensembl_version'=>$ensembl_version, '-comparison_name'=>$comparison_name);
print BLUE, "\nOf these, $ae_count are also AE", RESET;

#Load the feature annotations
print BLUE, "\n\nGetting feature annotations:", RESET;
my %master_gene_list;
my $feature_ref = &getFeatures('-annotation_dir'=>$annotation_dir);
my $feature_count = keys %{$feature_ref};
print BLUE, "\n\tFound $feature_count feature annotation records for the DE and SI records", RESET;

#Get basic gene info
print BLUE, "\n\nGetting basic gene info:", RESET;
my $genes_ref;
my $gene_transcripts_ref;
my @grand_orfs;
&getBasicGeneInfo();


#Determine the top N DE events for each gene.  
#Store junctions and exon regions separately
#Store gains and losses separately
#Disqualify exons below a certain size: ($peptide_size*3)+10
print BLUE, "\n\nBuilding candidate feature list:", RESET;
my %master_candidate_list;
my $candidates_ref = &buildCandidateList();

#Grab the actual sequences for the candidate features
print BLUE, "\n\nRetrieving candidate feature sequences:", RESET;
&getCandidateSeqs();

#Print a summary file (one for gains, one for losses) containing all candidate peptide features
#Also generate the row data for an html summary
my $nj_gain_content = '';
my $nj_loss_content = '';
my $all_gain_content = '';
my $all_loss_content = '';
my $top_n = 100;
print BLUE, "\n\nPrinting summary files:", RESET;
&printSummaryFiles('-gain_file'=>$summary_file_gain, '-loss_file'=>$summary_file_loss);

#Generate link to master candidate gene list file
my $summary_file_gain_name = "$comparison_name"."_peptides_GAIN.txt";
my $summary_file_loss_name = "$comparison_name"."_peptides_LOSS.txt";
my $data_path_gain = "data/$summary_file_gain_name";
my $data_path_loss = "data/$summary_file_loss_name";
my $data_link_gain = "<A HREF=\"$data_path_gain\">$summary_file_gain_name</A>";
my $data_link_loss = "<A HREF=\"$data_path_loss\">$summary_file_loss_name</A>";
my $download_content = '';
$download_content = "<P CLASS=\"Indented12LR_s16_bold\">Download complete peptide lists as tab delimited text file: $data_link_gain | $data_link_loss</P><BR>\n";

#List the names of the libraries
my $lib_names_content = '';
$lib_names_content .= "<P CLASS=\"Indented12LR_s16_bold\">Library A = $libA_name</P>\n";
$lib_names_content .= "<P CLASS=\"Indented12LR_s16_bold\">Library B = $libB_name</P><BR>\n";


#Define the legend
my $legend_content = "<div id=\"legend\" style=\"width: 1250px;\">\n";
$legend_content .= "<P CLASS=\"Indented24LR_s16\">The following table lists candidate peptides according to their difference in observed expression between the conditions.  Each peptide corresponds to a specific exon region or exon junction.  These may correspond to genes that are differentially expressed overall, or specific isoforms that are differentially expressed.  All candidates are divided into those 'gained' or 'lost' in library A compared to B.  The lists are further divided into peptides corresponding to novel junctions or any features type. 'Grand rank' refers to the overall ranking within each list.  'Rank' refers to the rank of peptides within a single gene if multiple candidates were identified. 'Gene Name' is the name assigned to this gene by EnsEMBL. This name links to the corresponding ALEXA-Seq gene expression record. 'Feature ID' is the ALEXA-Seq id assigned to the sequence feature used to predict the peptide sequence.  The 'Feature Name' is the name assigned to the exon region or junction. 'Feature Type' indicates whether the feature corresponds to an exon or junction.  'Fold Change' is the fold difference observed between the libraries for the feature.  'Lib A Level' and 'Lib B Level' are the expression levels observed for this feature in Library A and B respectively. 'Exons Skipped' indicates the number of exons skipped by a junction.  'Is AE?' indicates whether the feature was considered alternatively expressed.  'Expressed Count' indicates whether the feature was expressed above background noise levels in one library (1) or both (2). 'Gene ORF Count' indicates the number of known transcript ORFs for the gene.  'Gene ORF Matches' indicates the subset of these that the peptide matches.  'ORFeome Matches' indicates the number of ORFs in the entire known ORFeome that the peptide matches perfectly.  Finally the actual peptide sequence is provided.  Only the top $top_n events are shown in the following table. For a complete list and further details on each peptide, including the corresponding feature cDNA sequences, download the full data file above.</P><BR>\n";
$legend_content .= "</div>\n";

#Join the content sections together
my $div_count = 0;
my $div_name = '';
my $content = '';
$content .= $download_content;
$content .= $lib_names_content;

#Novel junctions - GAIN
$div_count++;
$div_name="box"."$div_count"."h";
$content .= "<P CLASS=\"Indented12LR_s16_Bold\">Candidate Peptides - Novel Junction Features Only (Gain Only)</P>\n";
$content .= "<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[$div_name]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>\n";
$content .= "<div id=\"$div_name\">\n";
$content .= "$legend_content";
$content .= "<TABLE CLASS=\"Data2\">\n";
$content .= $nj_gain_content;
$content .= "</TABLE><BR>\n</div><BR>\n";

#Novel junctions - LOSS
$div_count++;
$div_name="box"."$div_count"."h";
$content .= "<P CLASS=\"Indented12LR_s16_Bold\">Candidate Peptides - Novel Junction Features Only (Loss Only)</P>\n";
$content .= "<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[$div_name]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>\n";
$content .= "<div id=\"$div_name\">\n";
$content .= "$legend_content";
$content .= "<TABLE CLASS=\"Data2\">\n";
$content .= $nj_loss_content;
$content .= "</TABLE><BR>\n</div><BR>\n";

#All features - GAIN
$div_count++;
$div_name="box"."$div_count"."h";
$content .= "<P CLASS=\"Indented12LR_s16_Bold\">Candidate Peptides - All Exons and Junctions (Gain Only)</P>\n";
$content .= "<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[$div_name]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>\n";
$content .= "<div id=\"$div_name\">\n";
$content .= "$legend_content";
$content .= "<TABLE CLASS=\"Data2\">\n";
$content .= $all_gain_content;
$content .= "</TABLE><BR>\n</div><BR>\n";

#All features - LOSS
$div_count++;
$div_name="box"."$div_count"."h";
$content .= "<P CLASS=\"Indented12LR_s16_Bold\">Candidate Peptides - All Exons and Junctions (Loss Only)</P>\n";
$content .= "<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[$div_name]\" data-openimage=\"images/Minus_icon.gif\" data-closedimage=\"images/Plus_icon.gif\"><img src=\"images/Plus_icon.gif\" border=\"0\" /></a></P>\n";
$content .= "<div id=\"$div_name\">\n";
$content .= "$legend_content";
$content .= "<TABLE CLASS=\"Data2\">\n";
$content .= $all_loss_content;
$content .= "</TABLE><BR>\n</div><BR>\n";


#Write out the page
my $title = "Candidate peptides for comparison: '$comparison_name' ($comparison_id) - Project: $project";
my $meta_description = "Provides lists of candidate DE and AE peptides for comparison '$comparison_name' of project '$project'";
my $meta_keywords = "Differential Expression, Alternative Expression, $comparison_name, $project, antibody design, therapeutic antibodies, exon junction peptides";

&writePage('-path'=>$web_path, '-title'=>$title, '-content'=>\$content, '-css_path'=>"../ALEXA2.css", '-alexa_home_path'=>"$alexa_home_path", '-alexa_seq_home_path'=>"$alexa_seq_path", '-summary_path'=>"Summary.htm", '-genes_path'=>"genes/index.html", '-search_path'=>"$search_page_url", '-meta_description'=>$meta_description, '-meta_keywords'=>$meta_keywords, '-google_analytics'=>1, '-google_analytics_id'=>$google_analytics_id, '-collapse_div_script'=>'../animatedcollapse.js', '-jquery_min_script'=>'../jquery.min.js', '-div_count'=>$div_count);

my $mem_message = &memoryUsage();
print YELLOW, "\n\n$mem_message\n\n", RESET;
print "\n\n";

exit();



#################################################################################################################################################################
#Get DE values                                                                                                                                                  #
#################################################################################################################################################################
sub getDeValues{
  my %args = @_;
  my $de_dir = $args{'-de_dir'};
  my $ensembl_version = $args{'-ensembl_version'};
  my $comparison_name = $args{'-comparison_name'};

  #Get DE values for the following feature types
  my %types;
  $types{'ExonRegion'}{dir} = "ENST"."_v$ensembl_version";
  $types{'Junction'}{dir} = "Junctions_v$ensembl_version";

  my %de;

  print BLUE, "\n\nGetting DE values for DE features of each gene", RESET;
  foreach my $type (sort keys %types){
    my $de_dir = "$de_dir"."$types{$type}{dir}";
    $de_dir = &checkDir('-dir'=>$de_dir, '-clear'=>"no");
    my $de_file = "$de_dir"."$comparison_name"."_"."$type"."_DE_Values_Significant_MTC.txt";
    print YELLOW, "\n\t$de_file", RESET;

    open (DE, "$de_file") || die "\n\nCould not open DE file: $de_file\n\n";
    my $header = 1;
    my %columns;
    while(<DE>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        $header = 0;
        my $c = 0;
        foreach my $col (@line){
          $columns{$col}{pos} = $c;
          $c++;
        }
        next();
      }
      my $fid = $line[$columns{'FID'}{pos}];
      $de{$fid}{seq_name} = $line[$columns{'Seq_Name'}{pos}];
      $de{$fid}{a_expressed} = $line[$columns{'A_Expressed'}{pos}];
      $de{$fid}{b_expressed} = $line[$columns{'B_Expressed'}{pos}];
      $de{$fid}{a_norm} = $line[$columns{'A_Norm'}{pos}];
      $de{$fid}{b_norm} = $line[$columns{'B_Norm'}{pos}];
      $de{$fid}{fold_change} = $line[$columns{'Fold_Change'}{pos}];
      $de{$fid}{type} = $type;
      $de{$fid}{is_ae} = 0;
    }
    close(DE);
  }

  return(\%de);
}


#################################################################################################################################################################
#Get DE values                                                                                                                                                  #
#################################################################################################################################################################
sub getAeFeatures{
  my %args = @_;
  my $si_dir = $args{'-si_dir'};
  my $ensembl_version = $args{'-ensembl_version'};
  my $comparison_name = $args{'-comparison_name'};

  #Get DE values for the following feature types
  my %types;
  $types{'ExonRegion'}{dir} = "ENST"."_v$ensembl_version";
  $types{'Junction'}{dir} = "Junctions_v$ensembl_version";

  print BLUE, "\n\nChecking DE features to see if they are also AE", RESET;
  my $ae_count = 0;

  foreach my $type (sort keys %types){
    my $si_dir = "$si_dir"."$types{$type}{dir}";
    $si_dir = &checkDir('-dir'=>$si_dir, '-clear'=>"no");
    my $si_file = "$si_dir"."$comparison_name"."_"."$type"."_SI_Values_Sorted_Cutoff.txt";
    print YELLOW, "\n\t$si_file", RESET;

    open (SI, "$si_file") || die "\n\nCould not open DE file: $si_file\n\n";
    my $header = 1;
    my %columns;
    while(<SI>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        $header = 0;
        my $c = 0;
        foreach my $col (@line){
          $columns{$col}{pos} = $c;
          $c++;
        }
        next();
      }
      my $fid = $line[$columns{'FID'}{pos}];

      if ($de_ref->{$fid}){
        $de_ref->{$fid}->{is_ae} = 1;
        $ae_count++;
      }
    }
    close(SI);
  }
  return($ae_count);
}


#######################################################################################################################################################################
#Get basic gene info from ALEXA DB                                                                                                                                    #
#######################################################################################################################################################################
sub getBasicGeneInfo{

  my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
  my @gene_ids = keys %master_gene_list;
  $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"yes", '-silent'=>"yes");

  #Get a transcript object so that cds start/end coords will be available
  $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no", '-silent'=>"yes");
  $alexa_dbh->disconnect();

  #Assemble the ORF sequence
  #Note that the %master_gene_list should only contain protein coding genes as only protein coding features (and their corresponding genes) were allowed
  foreach my $gene_id (keys %{$gene_transcripts_ref}){
    my $trans_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};
    my $sequence = $genes_ref->{$gene_id}->{sequence};
    my $gene_name = $genes_ref->{$gene_id}->{gene_name};
    #print YELLOW, "\n\n\n$gene_name", RESET;
    foreach my $trans_id (sort keys %{$trans_ref}){

      my @cds_start_coords = @{$trans_ref->{$trans_id}->{cds_start_coords}};
      my @cds_end_coords = @{$trans_ref->{$trans_id}->{cds_end_coords}};
      my @temp = @cds_end_coords;
      my $cds_seq = "";
      
      foreach my $start (@cds_start_coords){
        my $end = shift (@temp);
        unless($end =~ /\d+/ && $start =~ /\d+/){
          next();
        }
        my $length = ($end-$start)+1;
        my $exon_cds_seq = substr($sequence, $start-1, $length);
        $cds_seq .= $exon_cds_seq;
      }
      $trans_ref->{$trans_id}->{cds_seq} = $cds_seq; 
      
      if ($cds_seq){
        #Translate the CDS sequence and store it
        my $seq_obj = Bio::Seq->new('-seq'=>$cds_seq, '-desc'=>'Sample Bio::Seq object', '-display_id' => 'something', '-accession_number' => 'accnum', '-alphabet' => 'dna' );
        my @seqs = Bio::SeqUtils->translate_3frames($seq_obj);
        my $cds_protein = $seqs[0]->seq();
        push(@grand_orfs, $cds_protein);
        $trans_ref->{$trans_id}->{cds_protein} = $cds_protein;
        #print YELLOW, "\nStarts: @cds_start_coords\nEnds: @cds_end_coords\n$cds_seq\n$cds_protein\n", RESET;
      }else{
        $trans_ref->{$trans_id}->{cds_protein} = "";
      }
    }
  }

  #Initialize a hash to store the feature list for each gene.
  foreach my $gene_id (keys %{$genes_ref}){
    my %fid_list;
    $genes_ref->{$gene_id}->{fid_list} = \%fid_list;
  }

  #Determine which region each gene corresponds to...
  foreach my $gene_id (keys %{$genes_ref}){
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
#Get the feature annotations                                                                                                                                          #
#######################################################################################################################################################################
sub getFeatures{
  my %args = @_;
  my $annotation_dir = $args{'-annotation_dir'};

  my %a;

  my %types;
  $types{'1'}{type} = "exonRegions";
  $types{'1'}{subtype} = "ExonRegion";
  $types{'1'}{file} = "exonRegions_annotated.txt.gz";
  $types{'2'}{type} = "exonJunctions";
  $types{'2'}{subtype} = "Junction";
  $types{'2'}{file} = "exonJunctions_"."$junction_seq_size"."mers_annotated.txt.gz";

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

      #Unless this FID is defined in the DE objects, dont bother storing it:
      unless ($de_ref->{$fid}){
        next();
      }

      my $base_count = $line[$columns{'Base_Count'}{column_pos}];
      my $coding_base_count = $line[$columns{'Coding_Base_Count'}{column_pos}];

      #Unless this FID corresponds to a feature that is large enough, dont bother storing it
      unless ($base_count >= (($peptide_size*3)+10)){
        delete($de_ref->{$fid});
        next();
      }

      #Unless this FID corresponds to a feature that is entirely within the CDS, don't bother storing it.
      unless ($base_count == $coding_base_count){
        delete($de_ref->{$fid});
        next();
      }

      #Skip features defined for overlapping genes for now?
      my $gene_count = scalar(@gene_ids);
      if ($gene_count > 1){
        delete($de_ref->{$fid});
        next();
      }

      $a{$fid}{id} = $line[0];
      $a{$fid}{subtype} = $subtype;
      $a{$fid}{gene_ids} = \@gene_ids;
      $a{$fid}{supporting_ensembl_count} = $line[$columns{'Supporting_EnsEMBL_Count'}{column_pos}];

      #Junction specific columns
      if ($type eq "exonJunctions"){
         $a{$fid}{exons_skipped} = $line[$columns{'Exons_Skipped'}{column_pos}];
      }else{
         $a{$fid}{exons_skipped} = "N/A";
      }

      #Store gene IDs associated with DE junctions and exon regions
      if ($type eq "exonRegions" || $type eq "exonJunctions"){
        foreach my $gene_id (@gene_ids){
          $master_gene_list{$gene_id}=1;
        }
      }
    }
    close(ANN);
  }

  return(\%a);
}


#######################################################################################################################################################################
#Build candidate list                                                                                                                                                 #
#######################################################################################################################################################################
sub buildCandidateList{

  my %candidates;

  #Go through all DE features (sorted by fold change) - first for gains, then for losses
  #Get Gene ID(s) associated with feature
  #Store top fold change observed (in each direction) and a list of top N features associated with each gene
  #Store exon, known junction and novel junction features seperately

  #Initialize arrays
  foreach my $gene_id (keys %master_gene_list){
    my @tmp1; my @tmp2; my @tmp3; my @tmp4; my @tmp5; my @tmp6;
    $candidates{$gene_id}{fc_gain_overall} = 0;
    $candidates{$gene_id}{fc_gain_known_exon} = 0;
    $candidates{$gene_id}{fc_gain_known_junction} = 0;
    $candidates{$gene_id}{fc_gain_novel_junction} = 0;
    $candidates{$gene_id}{fc_loss_overall} = 0;
    $candidates{$gene_id}{fc_loss_known_exon} = 0;
    $candidates{$gene_id}{fc_loss_known_junction} = 0;
    $candidates{$gene_id}{fc_loss_novel_junction} = 0;
    $candidates{$gene_id}{gain_known_exons} = \@tmp1;
    $candidates{$gene_id}{gain_known_junctions} = \@tmp2;
    $candidates{$gene_id}{gain_novel_junctions} = \@tmp3;
    $candidates{$gene_id}{loss_known_exons} = \@tmp4;
    $candidates{$gene_id}{loss_known_junctions} = \@tmp5;
    $candidates{$gene_id}{loss_novel_junctions} = \@tmp6;
  }

  #Gains
  foreach my $fid (sort {$de_ref->{$b}->{fold_change} <=> $de_ref->{$a}->{fold_change}} keys %{$de_ref}){
    my $fc = $de_ref->{$fid}->{fold_change};
    my $id = $feature_ref->{$fid}->{id};
    if ($fc < 0){
      last();
    }
    my $type = $de_ref->{$fid}->{type};
    my $supporting_ensembl_count = $feature_ref->{$fid}->{supporting_ensembl_count};
    my @gene_ids = @{$feature_ref->{$fid}->{gene_ids}};

    foreach my $gene_id (@gene_ids){
      if ($type eq "ExonRegion"){
        $de_ref->{$fid}->{subtype} = "Exon";
        if ($fc > $candidates{$gene_id}{fc_gain_known_exon}){
          $candidates{$gene_id}{fc_gain_known_exon} = $fc;
        }
        if ($fc > $candidates{$gene_id}{fc_gain_overall}){
          $candidates{$gene_id}{fc_gain_overall} = $fc;
        }
        unless(scalar(@{$candidates{$gene_id}{gain_known_exons}}) >= $feature_limit){
          push(@{$candidates{$gene_id}{gain_known_exons}}, $fid);
          $master_candidate_list{$id}{seq}="";
        }

      }elsif($type eq "Junction" && $supporting_ensembl_count > 0){
        $de_ref->{$fid}->{subtype} = "KnownJunction";
        if ($fc > $candidates{$gene_id}{fc_gain_known_junction}){
          $candidates{$gene_id}{fc_gain_known_junction} = $fc;
        }
        if ($fc > $candidates{$gene_id}{fc_gain_overall}){
          $candidates{$gene_id}{fc_gain_overall} = $fc;
        }
        unless(scalar(@{$candidates{$gene_id}{gain_known_junctions}}) >= $feature_limit){
          push(@{$candidates{$gene_id}{gain_known_junctions}}, $fid);
          $master_candidate_list{$id}{seq}="";
        }

      }elsif($type eq "Junction" && $supporting_ensembl_count == 0){
        $de_ref->{$fid}->{subtype} = "NovelJunction";
        if ($fc > $candidates{$gene_id}{fc_gain_novel_junction}){
          $candidates{$gene_id}{fc_gain_novel_junction} = $fc;
        }
        if ($fc > $candidates{$gene_id}{fc_gain_overall}){
          $candidates{$gene_id}{fc_gain_overall} = $fc;
        }
        unless(scalar(@{$candidates{$gene_id}{gain_novel_junctions}}) >= $feature_limit){
          push(@{$candidates{$gene_id}{gain_novel_junctions}}, $fid);
          $master_candidate_list{$id}{seq}="";
        }
      }
    }
  }

  #Losses
  foreach my $fid (sort {$de_ref->{$a}->{fold_change} <=> $de_ref->{$b}->{fold_change}} keys %{$de_ref}){
    my $fc = $de_ref->{$fid}->{fold_change};
    my $id = $feature_ref->{$fid}->{id};
    if ($fc > 0){
      last();
    }
    my $type = $de_ref->{$fid}->{type};
    my $supporting_ensembl_count = $feature_ref->{$fid}->{supporting_ensembl_count};
    my @gene_ids = @{$feature_ref->{$fid}->{gene_ids}};

    foreach my $gene_id (@gene_ids){
      if ($type eq "ExonRegion"){
        $de_ref->{$fid}->{subtype} = "Exon";
        if ($fc < $candidates{$gene_id}{fc_loss_known_exon}){
          $candidates{$gene_id}{fc_loss_known_exon} = $fc;
        }
        if ($fc < $candidates{$gene_id}{fc_loss_overall}){
          $candidates{$gene_id}{fc_loss_overall} = $fc;
        }
        unless(scalar(@{$candidates{$gene_id}{loss_known_exons}}) >= $feature_limit){
          push(@{$candidates{$gene_id}{loss_known_exons}}, $fid);
          $master_candidate_list{$id}{seq}="";
        }

      }elsif($type eq "Junction" && $supporting_ensembl_count > 0){
        $de_ref->{$fid}->{subtype} = "KnownJunction";
        if ($fc < $candidates{$gene_id}{fc_loss_known_junction}){
          $candidates{$gene_id}{fc_loss_known_junction} = $fc;
        }
        if ($fc < $candidates{$gene_id}{fc_loss_overall}){
          $candidates{$gene_id}{fc_loss_overall} = $fc;
        }
        unless(scalar(@{$candidates{$gene_id}{loss_known_junctions}}) >= $feature_limit){
          push(@{$candidates{$gene_id}{loss_known_junctions}}, $fid);
          $master_candidate_list{$id}{seq}="";
        }

      }elsif($type eq "Junction" && $supporting_ensembl_count == 0){
        $de_ref->{$fid}->{subtype} = "NovelJunction";
        if ($fc < $candidates{$gene_id}{fc_loss_novel_junction}){
          $candidates{$gene_id}{fc_loss_novel_junction} = $fc;
        }
        if ($fc < $candidates{$gene_id}{fc_loss_overall}){
          $candidates{$gene_id}{fc_loss_overall} = $fc;
        }
        unless(scalar(@{$candidates{$gene_id}{loss_novel_junctions}}) >= $feature_limit){
          push(@{$candidates{$gene_id}{loss_novel_junctions}}, $fid);
          $master_candidate_list{$id}{seq}="";
        }
      }
    }
  }

  return(\%candidates);
}


#######################################################################################################################################################################
#Grab the actual sequences for the candidate features                                                                                                                 #
#######################################################################################################################################################################
sub getCandidateSeqs{

  my %types;
  $types{'1'}{type} = "exonRegions";
  $types{'1'}{subtype} = "ExonRegion";
  $types{'1'}{file} = "exonRegions.fa.gz";
  $types{'2'}{type} = "exonJunctions";
  $types{'2'}{subtype} = "Junction";
  $types{'2'}{file} = "exonJunctions_"."$junction_seq_size"."mers.fa.gz";

  foreach my $c (sort {$a <=> $b} keys %types){
    my $type = $types{$c}{type};
    my $subtype = $types{$c}{subtype};
    print BLUE, "\n\tProcessing type: $type ($subtype)", RESET;

    my $file = $types{$c}{file};
    my $dir = "$annotation_dir"."$type/blastdb/";
    $dir = &checkDir('-dir'=>$dir, '-clear'=>"no");
    my $file_path = "$dir"."$file";

    open(SEQ, "zcat $file_path |") || die "\nCould not open annotation file: $file_path\n\n";
    my $match_found = 0;
    my $id;
    while(<SEQ>){
      chomp($_);
      my $line = $_;
      if ($_ =~ /\>(.*)/){
        $id = $1;
        if ($master_candidate_list{$id}){
          #print YELLOW, "\n$id", RESET;
          $match_found = 1;
        }
        next();
      }
      if ($match_found){
        $master_candidate_list{$id}{seq} = $_;
        #print YELLOW, "\n$line", RESET;
        $match_found = 0;
      }
    }
    close(SEQ);
  }

  my $seqs_found = keys %master_candidate_list;
  print BLUE, "\n\tFound $seqs_found SEQs", RESET;

  #Now go through each candidate sequence and determine the peptide sequence
  my $c = 0;
  while (my ($fid) = each %{$de_ref}){
    my $id = $feature_ref->{$fid}->{id};
    unless ($master_candidate_list{$id}){
      next();
    }

    $c++;
    if ($c == 1000){
      $| = 1; print BLUE, ".", RESET;  $| = 0;
      $c = 0;
    }
    my $seq = $master_candidate_list{$id}{seq};

    #print YELLOW, "\nFID: $fid\tID: $id\tTranslating seq:\n$seq", RESET;
    
    #Translate the feature sequence
    my $seq_obj = Bio::Seq->new('-seq'=>$seq, '-desc'=>'Sample Bio::Seq object', '-display_id' => 'something', '-accession_number' => 'accnum', '-alphabet' => 'dna');
    my @seqs = Bio::SeqUtils->translate_3frames($seq_obj);

    my @frames = ($seqs[0]->seq(), $seqs[1]->seq(), $seqs[2]->seq());

    #Figure out which frame is correct for this feature.  Do this by matching a test sequence from within it to the ORFs of the corresponding gene
    #First get the test peptide for each feature
    #For exons use the whole thing but inset by two peptides on each end
    my $inset = 1;
    my $type = $de_ref->{$fid}->{type};
    my $subtype = $de_ref->{$fid}->{subtype};
    my @subframes1;
    my @subframes2;
    if ($type eq "ExonRegion"){
      #Store a test peptide from the center of the exon
      foreach my $peptide (@frames){
        my $size = length($peptide);
        my $subframe = substr($peptide, $inset, $size-($inset*2));
        push(@subframes1, $subframe);
      }
    }elsif($type eq "Junction"){
      #Store a test peptide from the 5' side of the junction
      foreach my $peptide (@frames){
        my $size = length($peptide);
        my $half_size = sprintf("%.0f", $size/2);
        my $subframe = substr($peptide, $inset, $half_size-($inset*2));
        push(@subframes1, $subframe);
      }
      #Store a test peptide from the 3' end of the junction
      foreach my $peptide (@frames){
        my $size = length($peptide);
        my $half_size = sprintf("%.0f", $size/2);
        my $subframe = substr($peptide, $half_size+$inset, $half_size-($inset*2));
        push(@subframes2, $subframe);
      }
    }else{
      print RED, "\n\nType: $type not understood", RESET;
      exit();
    }
    #print YELLOW, "\n\n$id\t$fid\t$type\t$subtype\nFRA: @frames\nSUB1: @subframes1\nSUB2: @subframes2\nFrame Test: ", RESET;

    #Now test these test peptides against the gene's ORFs to look for a match
    my @gene_ids = @{$feature_ref->{$fid}->{gene_ids}};
    my $gene_id = $gene_ids[0];
    my $trans_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};
    my $matching_frame1 = 0;
    my $matching_frame2 = 0;
    my %matching_frames1;
    my %matching_frames2;
    my @orfs;
    foreach my $trans_id (sort keys %{$trans_ref}){
      if ($trans_ref->{$trans_id}->{cds_protein}){
        my $cds_protein = $trans_ref->{$trans_id}->{cds_protein};
        push(@orfs, $cds_protein);
      }
    }

    #Find frame matches for test peptide 1
    my $f = 0;
    foreach my $frame (@subframes1){
      $frame=~ s/\*/\\*/g;
      my $test = scalar(grep(/$frame/, @orfs));
      #print YELLOW, "$f = $test\t", RESET;
      if ($test){
        $matching_frame1=$f;
        $matching_frames1{$f}=1;
      }
      $f++;
    }
    my $matching_frame1_count = keys %matching_frames1;

    #Find frame matches for test peptide 2
    $f = 0;
    foreach my $frame (@subframes2){
      $frame=~ s/\*/\\*/g;
      my $test = scalar(grep(/$frame/, @orfs));
      #print YELLOW, "$f = $test\t", RESET;
      if ($test){
        $matching_frame2=$f;
        $matching_frames2{$f}=1;
      }
      $f++;
    }
    my $matching_frame2_count = keys %matching_frames2;

    #If the frame was successfully found get a new test peptide using this frame
    #Determine the number times this peptide occurs within the known ORFs of the gene (of how many) and within the entire ORFeome
    if ($matching_frame1_count == 1){
      my $target_protein = $frames[$matching_frame1];
  
      #Now get an amino acid of the target length from within the translated feature
      my $size = length($target_protein);
      my $center = sprintf("%.0f", $size/2);
      my $flank = sprintf("%.0f", $peptide_size/2);
      my $target_peptide = substr($target_protein, $center-$flank, $peptide_size);
      #print YELLOW, "\nSelecting peptide: $target_peptide", RESET;


      #Grab the cDNA sequence corresponding to this peptide
      my $cdna_size = (length($target_peptide))*3;
      my $cdna_center = sprintf("%.0f", ((length($seq))/2));
      my $cdna_flank = sprintf("%.0f", $cdna_size/2);
      my $target_seq = substr($seq, $cdna_center-$cdna_flank, $cdna_size);

      #How many ORFs does the source gene have?
      my $orf_count = scalar(@orfs);

      #How many of these ORFs contain the selected peptide?
      my $temp = $target_peptide;
      $temp=~ s/\*/\\*/g;
      my $within_gene_matches = scalar(grep(/$temp/, @orfs));

      #How many additional ORFs in the ORFeome contain the selected peptide?
      my $orfeome_matches = scalar(grep(/$temp/, @grand_orfs));
      my $outside_gene_matches = $orfeome_matches - $within_gene_matches;
      #print YELLOW, "\nThis peptide occurs in $within_gene_matches of $orf_count ORFs of THIS gene", RESET;
      #print YELLOW, "\nIt occurs in $outside_gene_matches ORFs from OTHER genes", RESET;

      #Was the frame maintained?
      my $frame_maintained = "no";
      if ($type eq "ExonRegion" || ($matching_frame2_count == 1 && $matching_frame1 == $matching_frame2)){
        $frame_maintained = "yes";
      }

      #Store these values for later
      $de_ref->{$fid}->{found_peptide} = 1;
      $de_ref->{$fid}->{target_peptide} = $target_peptide;
      $de_ref->{$fid}->{target_seq} = $target_seq;
      $de_ref->{$fid}->{within_gene_matches} = $within_gene_matches;
      $de_ref->{$fid}->{gene_orf_count} = $orf_count;
      $de_ref->{$fid}->{orfeome_matches} = $orfeome_matches;
      $de_ref->{$fid}->{outside_gene_matches} = $outside_gene_matches;
      $de_ref->{$fid}->{frame_maintained} = $frame_maintained;

    }else{
      $de_ref->{$fid}->{found_peptide} = 0;
    }

  }

  return();
}

#######################################################################################################################################################################
#Print a summary file (one for gains, one for losses) containing all candidate peptide features                                                                       #
#######################################################################################################################################################################
sub printSummaryFiles{
  my %args = @_;
  my $gain_file = $args{'-gain_file'};
  my $loss_file = $args{'-loss_file'};

  my %files;
  $files{'1'}{file} = $gain_file;
  $files{'1'}{sort_col} = "fc_gain_overall";
  $files{'1'}{known_exons_col} = "gain_known_exons";
  $files{'1'}{known_junctions_col} = "gain_known_junctions";
  $files{'1'}{novel_junctions_col} = "gain_novel_junctions";
  $files{'1'}{direction} = "GAIN";
  $files{'2'}{file} = $loss_file;
  $files{'2'}{sort_col} = "fc_loss_overall";
  $files{'2'}{known_exons_col} = "loss_known_exons";
  $files{'2'}{known_junctions_col} = "loss_known_junctions";
  $files{'2'}{novel_junctions_col} = "loss_novel_junctions";
  $files{'2'}{direction} = "LOSS";

  foreach my $f (sort {$a <=> $b} keys %files){
    my $file = $files{$f}{file};
    my $sort_col = $files{$f}{sort_col};
    my $known_exons_col = $files{$f}{known_exons_col};
    my $known_junctions_col = $files{$f}{known_junctions_col};
    my $novel_junctions_col = $files{$f}{novel_junctions_col};
    my $c = 0;

    print BLUE, "\n\tGenerating: $file", RESET;
    open (OUT, ">$file") || die "\n\nCould not open output file: $file\n\n";

    #Print header line
    print OUT "Grand Rank\tRank\tGene Name\tENSG\tFeature ID\tFeature Name\tFeature Type\tFold Change\tLibA Expression\tLibB Expression\tExons Skipped\tIs AE?\tExpressed Count (1 or 2)\tGene ORF count\tGene ORF matches\tORFeome matches\tFrame Maintained?\tPeptide Sequence\tcDNA Sequence\n";

    #Go through each candidate gene, and print all the peptide candidates for that gene
    my $gr = 0;
    foreach my $gene_id (sort {abs($candidates_ref->{$b}->{$sort_col}) <=> abs($candidates_ref->{$a}->{$sort_col})} keys %{$candidates_ref}){
      my $gene_name = $genes_ref->{$gene_id}->{gene_name};
      my $ensg = $genes_ref->{$gene_id}->{ensembl_g_id};

      #Get all the candidate FIDs for this gene
      my @fids;
      push (@fids, @{$candidates_ref->{$gene_id}->{$known_exons_col}});
      push (@fids, @{$candidates_ref->{$gene_id}->{$known_junctions_col}});
      push (@fids, @{$candidates_ref->{$gene_id}->{$novel_junctions_col}}); 

      my $fid_count = scalar(@fids);
      unless ($fid_count > 0){
        next();
      }
      
      #Order the candidates for this gene by their fold change
      my %fids;
      foreach my $fid (@fids){
        #Unless a peptide was found for this feature, skip it
        unless ($de_ref->{$fid}->{found_peptide}){
          next();
        }
        $fids{$fid}{fc} = $de_ref->{$fid}->{fold_change};
      }

      #Print them out according to this order
      my $r = 0;
      foreach my $fid (sort {abs($fids{$b}->{fc}) <=> abs($fids{$a}->{fc})} keys %fids){
        $gr++;
        $r++;
        my $fc = sprintf("%.2f", $fids{$fid}{fc});
        my $subtype = $de_ref->{$fid}->{subtype};
        my $feature_name = $de_ref->{$fid}->{seq_name};
        my $lib_a_exp = sprintf("%.2f", $de_ref->{$fid}->{a_norm});
        my $lib_b_exp = sprintf("%.2f", $de_ref->{$fid}->{b_norm});
        my $exons_skipped = $feature_ref->{$fid}->{exons_skipped};
        my $is_ae = $de_ref->{$fid}->{is_ae};
        my $expressed_count = $de_ref->{$fid}->{a_expressed} + $de_ref->{$fid}->{b_expressed};
        my $gene_orf_count = $de_ref->{$fid}->{gene_orf_count};
        my $within_gene_matches = $de_ref->{$fid}->{within_gene_matches};
        my $orfeome_matches = $de_ref->{$fid}->{orfeome_matches};
        my $target_peptide = $de_ref->{$fid}->{target_peptide};
        my $target_seq = $de_ref->{$fid}->{target_seq};
        my $frame_maintained = $de_ref->{$fid}->{frame_maintained};
        print OUT "$gr\t$r\t$gene_name\t$ensg\t$fid\t$feature_name\t$subtype\t$fc\t$lib_b_exp\t$lib_a_exp\t$exons_skipped\t$is_ae\t$expressed_count\t$gene_orf_count\t$within_gene_matches\t$orfeome_matches\t$frame_maintained\t$target_peptide\t$target_seq\n";
      }
    }
    close (OUT);
  }


  #Now go through the features again and produce the sorted html output tables for the following sets
  #1. Novel junctions only, gains only
  #2. Novel junctions only, losses only
  #3. All features, gains only
  #4. All features, losses only
  my %sets;
  $sets{'1'}{sort_col} = "fc_gain_novel_junction";
  $sets{'1'}{known_exons_col} = "gain_known_exons";
  $sets{'1'}{known_junctions_col} = "gain_known_junctions";
  $sets{'1'}{novel_junctions_col} = "gain_novel_junctions";
  $sets{'2'}{sort_col} = "fc_loss_novel_junction";
  $sets{'2'}{known_exons_col} = "loss_known_exons";
  $sets{'2'}{known_junctions_col} = "loss_known_junctions";
  $sets{'2'}{novel_junctions_col} = "loss_novel_junctions";
  $sets{'3'}{sort_col} = "fc_gain_overall";
  $sets{'3'}{known_exons_col} = "gain_known_exons";
  $sets{'3'}{known_junctions_col} = "gain_known_junctions";
  $sets{'3'}{novel_junctions_col} = "gain_novel_junctions";
  $sets{'4'}{sort_col} = "fc_loss_overall";
  $sets{'4'}{known_exons_col} = "loss_known_exons";
  $sets{'4'}{known_junctions_col} = "loss_known_junctions";
  $sets{'4'}{novel_junctions_col} = "loss_novel_junctions";

  #Add the table header for all tables
  my $header_row = "<TR><TD CLASS=\"Head3\">Grand<br>Rank</TD><TD CLASS=\"Head3\">Rank</TD><TD CLASS=\"Head3\">Gene<br>Name</TD><TD CLASS=\"Head3\">Feature<br>ID</TD><TD CLASS=\"Head3\">Feature<br>Name</TD><TD CLASS=\"Head3\">Feature<br>Type</TD><TD CLASS=\"Head3\">Fold<br>Change</TD><TD CLASS=\"Head3\">LibA<br>Level</TD><TD CLASS=\"Head3\">LibB<br>Level</TD><TD CLASS=\"Head3\">Exons<br>Skipped</TD><TD CLASS=\"Head3\">Is<br>AE?</TD><TD CLASS=\"Head3\">Expressed<br>Count<br>(1 or 2)</TD><TD CLASS=\"Head3\">Gene<br>ORF<br>Count</TD><TD CLASS=\"Head3\">Gene<br>ORF<br>matches</TD><TD CLASS=\"Head3\">ORFeome<br>matches</TD><TD CLASS=\"Head3\">Frame<br>Maintained?</TD><TD CLASS=\"Head3\">Peptide<br>Sequence</TD></TR>\n";
  $nj_gain_content .= $header_row;
  $nj_loss_content .= $header_row;
  $all_gain_content .= $header_row;
  $all_loss_content .= $header_row;

  foreach my $set (sort {$a <=> $b} keys %sets){
    my $sort_col = $sets{$set}{sort_col};
    my $known_exons_col = $sets{$set}{known_exons_col};
    my $known_junctions_col = $sets{$set}{known_junctions_col};
    my $novel_junctions_col = $sets{$set}{novel_junctions_col};

    #Go through each candidate gene, and print all the peptide candidates for that gene
    my $gr = 0;
    foreach my $gene_id (sort {abs($candidates_ref->{$b}->{$sort_col}) <=> abs($candidates_ref->{$a}->{$sort_col})} keys %{$candidates_ref}){
      my $gene_name = $genes_ref->{$gene_id}->{gene_name};
      my $ensg = $genes_ref->{$gene_id}->{ensembl_g_id};

      #Get all the candidate FIDs for this gene - limit to those needed for the current set
      my @fids;
      if ($set == 1 || $set == 2){
        push (@fids, @{$candidates_ref->{$gene_id}->{$novel_junctions_col}}); 
      }elsif ($set == 3 || $set == 4){
        push (@fids, @{$candidates_ref->{$gene_id}->{$known_exons_col}});
        push (@fids, @{$candidates_ref->{$gene_id}->{$known_junctions_col}});
        push (@fids, @{$candidates_ref->{$gene_id}->{$novel_junctions_col}}); 
      }
      my $fid_count = scalar(@fids);
      unless ($fid_count > 0){
        next();
      }
      
      #Order the candidates for this gene by their fold change
      my %fids;
      foreach my $fid (@fids){
        #Unless a peptide was found for this feature, skip it
        unless ($de_ref->{$fid}->{found_peptide}){
          next();
        }
        $fids{$fid}{fc} = $de_ref->{$fid}->{fold_change};
      }

      #Print them out according to this order
      my $r = 0;
      foreach my $fid (sort {abs($fids{$b}->{fc}) <=> abs($fids{$a}->{fc})} keys %fids){
        $gr++;
        $r++;
        #Only show the top_n results 
        if ($gr > $top_n){
          last();
        }

        my $fc = sprintf("%.2f", $fids{$fid}{fc});
        my $subtype = $de_ref->{$fid}->{subtype};
        my $feature_name = $de_ref->{$fid}->{seq_name};
        my $lib_a_exp = sprintf("%.2f", $de_ref->{$fid}->{a_norm});
        my $lib_b_exp = sprintf("%.2f", $de_ref->{$fid}->{b_norm});
        my $exons_skipped = $feature_ref->{$fid}->{exons_skipped};
        my $is_ae = $de_ref->{$fid}->{is_ae};
        my $expressed_count = $de_ref->{$fid}->{a_expressed} + $de_ref->{$fid}->{b_expressed};
        my $gene_orf_count = $de_ref->{$fid}->{gene_orf_count};
        my $within_gene_matches = $de_ref->{$fid}->{within_gene_matches};
        my $orfeome_matches = $de_ref->{$fid}->{orfeome_matches};
        my $target_peptide = $de_ref->{$fid}->{target_peptide};
        my $target_seq = $de_ref->{$fid}->{target_seq};
        my $gene_record_path = "genes/$genes_ref->{$gene_id}->{partition}/$ensg".".htm";
        my $gene_record_link = "<A HREF=\"$gene_record_path\">$gene_name</A>";
        my $frame_maintained = $de_ref->{$fid}->{frame_maintained};
        
        my $row = "<TR><TD CLASS=\"Head3\">$gr<TD CLASS=\"Head3\">$r</TD><TD CLASS=\"Head3\">$gene_record_link</TD><TD CLASS=\"Head3\">$fid</TD><TD CLASS=\"Head3\">$feature_name</TD><TD CLASS=\"Head3\">$subtype</TD><TD CLASS=\"Head3\">$fc</TD><TD CLASS=\"Head3\">$lib_b_exp</TD><TD CLASS=\"Head3\">$lib_a_exp</TD><TD CLASS=\"Head3\">$exons_skipped</TD><TD CLASS=\"Head3\">$is_ae</TD><TD CLASS=\"Head3\">$expressed_count</TD><TD CLASS=\"Head3\">$gene_orf_count</TD><TD CLASS=\"Head3\">$within_gene_matches</TD><TD CLASS=\"Head3\">$orfeome_matches</TD><TD CLASS=\"Head3\">$frame_maintained</TD><TD CLASS=\"Head3\">$target_peptide</TD></TR>\n";

        if ($set == 1){$nj_gain_content .= $row;}
        if ($set == 2){$nj_loss_content .= $row;}
        if ($set == 3){$all_gain_content .= $row;}
        if ($set == 4){$all_loss_content .= $row;}
      
      }
    }
  }

  return();
}

