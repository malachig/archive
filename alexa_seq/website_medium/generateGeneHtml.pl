#!/usr/bin/perl -w
#Written by Malachi Griffith and Obi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to generate a simple webpage to summarize expression, differential expression, splicing, etc. data for each individual gene
#All data will be loaded in, but will be processed on a gene by gene basis...

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use BerkeleyDB;

#Load the ALEXA modules
my $script_dir;
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    $script_dir = $1;
    push (@INC, $1);
  }
}
use website_medium::WEB qw(:all);
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);
use utilities::mapping qw(:all);

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $annotation_dir = '';
my $junction_seq_size = '';
my $boundary_seq_size = '';
my $entrez_annotations = '';
my $paths_file = '';  #File containing root paths to basic types of data (expression, DE and SI)
my $ensembl_version = '';
my $partition_file = '';
my $ucsc_build = '';
my $species = '';
my $track_url = '';
my $temp_dir = '';
my $test_genes  = '';
my $chr_filter = '';
my $search_page_url = '';
my $project_name = '';
my $species_name = '';
my $logfile = '';
my $alexa_home_path = '';
my $alexa_seq_path = '';
my $google_analytics_id = '';
my $test = '';

GetOptions('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'annotation_dir=s'=>\$annotation_dir, 
           'junction_seq_size=i'=>\$junction_seq_size, 'boundary_seq_size=i'=>\$boundary_seq_size, 'entrez_annotations=s'=>\$entrez_annotations,
           'paths_file=s'=>\$paths_file, 'ensembl_version=s'=>\$ensembl_version, 'partition_file=s'=>\$partition_file,
           'ucsc_build=s'=>\$ucsc_build, 'species=s'=>\$species, 'track_url=s'=>\$track_url, 'temp_dir=s'=>\$temp_dir, 'test_genes=s'=>\$test_genes, 'chr_filter=s'=>\$chr_filter,
           'search_page_url=s'=>\$search_page_url, 'project_name=s'=>\$project_name, 'species_name=s'=>\$species_name, 'logfile=s'=>\$logfile, 'alexa_home_path=s'=>\$alexa_home_path, 'alexa_seq_path=s'=>\$alexa_seq_path,
           'google_analytics_id=s'=>\$google_analytics_id, 'test=s'=>\$test);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script generates HTML pages to summarize expression, DE and SI results", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the base path to annotation files using:  --annotation_dir", RESET;
print GREEN, "\n\tSpecify the junction database sequence length using:  --junction_seq_size", RESET;
print GREEN, "\n\tSpecify the boundary database sequence length using:  --boundary_seq_size", RESET;
print GREEN, "\n\tSpecify a file containing root paths to data using: --paths_file", RESET;
print GREEN, "\n\tSpecify the version of EnsEMBL to get data for using: --ensembl_version", RESET;
print GREEN, "\n\tSpecify a file containing the coordinates of genome partitions used for the analysis using: --partition_file", RESET;
print GREEN, "\n\tSpecify the UCSC build to use for custom tracks using: --ucsc_build (e.g. hg18)", RESET;
print GREEN, "\n\tSpecify the species name using: --species (e.g. Homo sapiens)", RESET;
print GREEN, "\n\tSpecify the html path to UCSC track files using: --track_url", RESET;
print GREEN, "\n\tTo improve performance of file I/O, files will first be written to a temp_dir and then moved to the web dir:  --temp_dir", RESET;
print GREEN, "\n\tSpecify the URL to your Xapian-Omega search page using:  --search_page_url", RESET;
print GREEN, "\n\tSpecify a short project name using:  --project_name", RESET;
print GREEN, "\n\tSpecify the common name of the species using: --species_name", RESET;
print GREEN, "\n\tSpecify a log file using --logfile", RESET;
print GREEN, "\n\tSpecify the URL to your ALEXA home page using: --alexa_home_path", RESET;
print GREEN, "\n\tSpecify the URL to your ALEXA-Seq results page using: --alexa_seq_path", RESET;
print GREEN, "\n\tTo test on a small list of specific genes, use the following example: --test_genes='ENSG00000053747,ENSG00000114491'", RESET;
print GREEN, "\n\tTo limit to the genes of a specific genome region, use:  --chr_filter='3:16:121020102-126260001'", RESET;
print GREEN, "\n\tSpecify your Google Analytics ID using: --google_analytics_id", RESET;
print GREEN, "\n\tUse --test=1 to help testing R code", RESET;

print GREEN, "\n\nUsage: generateGeneHtml.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --annotation_dir=/projects/alexa/sequence_databases/hs_53_36o/  --junction_seq_size=62  --boundary_seq_size=62  --entrez_annotations=/projects/alexa/sequence_databases/hs_53_36o/Homo_sapiens.gene_info  --paths_file=/projects/malachig/alexa_seq/batch_jobs/5FU/5FU_Lib_Paths.txt  --ensembl_version=53  --ucsc_build=hg18  --species='Homo sapiens'  --partition_file=/projects/alexa/sequence_databases/hs_53_36o/Regions_250_Genes.txt  --track_url=http://www.alexaplatform.org/alexa_seq/5FU/ucsc/  --temp_dir=/projects/malachig/alexa_seq/temp/website/5FU/genes/  --search_page_url=http://www.bcgsc.ca/xapian-search/omega  --project_name='5-FU Resistance'  --species_name='Human'  --logfile=/projects/malachig/alexa_seq/logs/website/5FU/genes/generateGeneHtml_test_LOG.txt  --alexa_home_path=http://www.alexaplatform.org/index.htm  --alexa_seq_path=http://www.alexaplatform.org/alexa_seq/results.htm  --google_analytics_id=UA-xxxxxx-x\n\n", RESET;

#Check user supplied options
unless ($database && $server && $user && $password && $annotation_dir && $junction_seq_size && $boundary_seq_size && $entrez_annotations && $paths_file && $ensembl_version && $partition_file && $ucsc_build && $species && $track_url && $temp_dir && $search_page_url && $project_name && $species_name && $logfile && $alexa_home_path && $alexa_seq_path && $google_analytics_id){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}
$annotation_dir = &checkDir('-dir'=>$annotation_dir, '-clear'=>"no");

open(LOG, ">$logfile") || die "\n\nCould not open logfile: $logfile\n\n";

#Define the data types that will be processed in this script
my @types_list = qw (Gene Transcript ExonRegion Junction Boundary Intron ActiveIntronRegion SilentIntronRegion Intergenic ActiveIntergenicRegion SilentIntergenicRegion); 

my %types_list;
foreach my $type (@types_list){
  $types_list{$type}=1;
}
my $r_cmd;

#0.) Get annotations for genes, transcripts, and all seq types.
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
my $a_ref;

#Get list of genes to be processed
my @gene_id_list;
my %gene_id_list;

my $gene_data_file;
my $single_gene_data_file;
my $gene_model_file;
my ($region_number, $start_filter, $end_filter);
my $gene_transcripts_ref;

if ($test_genes){
  $temp_dir = &checkDir('-dir'=>$temp_dir, '-clear'=>"no");
  $gene_data_file = "$temp_dir"."GeneData_TEMP.txt";
  $single_gene_data_file = "$temp_dir"."single_gene_data_TEMP.txt";
  $gene_model_file = "$temp_dir"."gene_model_TEMP.txt";

  my @ensembl_g_ids = split(",", $test_genes);
  my %gene_ids = %{&getGeneIds ('-dbh'=>$alexa_dbh, '-ensembl_g_ids'=>\@ensembl_g_ids)};

  foreach my $ensg (keys %gene_ids){
    my $id = $gene_ids{$ensg}{alexa_gene_id};
    push (@gene_id_list, $id);
  }
  $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_id_list, '-sequence'=>"no");

}elsif($chr_filter){
  $temp_dir = &checkDir('-dir'=>$temp_dir, '-clear'=>"no");
  if ($chr_filter =~ /(.*)\:(\d+)\:(\d+)\-(\d+)/){
    $chr_filter = $1;
    $region_number = $2;
    $start_filter = $3;
    $end_filter = $4;
    unless ($end_filter > $start_filter){
      print RED, "\nStart of range must be smaller than end ($chr_filter)\n\n", RESET;
      print LOG "\nStart of range must be smaller than end ($chr_filter)\n\n";
      exit();
    }
  }else{
    print RED, "\nFormat of chr_filter not understood: $chr_filter (should be of the form:  Y:1:1-9999001)\n\n", RESET;
    print LOG "\nFormat of chr_filter not understood: $chr_filter (should be of the form:  Y:1:1-9999001)\n\n";
    exit();
  }

  $gene_data_file = "$temp_dir"."chr$chr_filter"."_"."$region_number"."_GeneData_TEMP.txt";
  $single_gene_data_file = "$temp_dir"."chr$chr_filter"."_"."$region_number"."_single_gene_data_TEMP.txt";
  $gene_model_file = "$temp_dir"."chr$chr_filter"."_"."$region_number"."_gene_model_TEMP.txt";

  #Get all genes for just this region
  @gene_id_list = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All', '-chromosome'=>$chr_filter, '-range'=>"$start_filter-$end_filter")};
  $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_id_list, '-sequence'=>"no");
 
}else{
  $temp_dir = &checkDir('-dir'=>$temp_dir, '-clear'=>"yes", '-force'=>"yes", '-recursive'=>"yes");
  $gene_data_file = "$temp_dir"."GeneData_TEMP.txt";
  $single_gene_data_file = "$temp_dir"."single_gene_data_TEMP.txt";
  $gene_model_file = "$temp_dir"."gene_model_TEMP.txt";
  @gene_id_list = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};
  $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_id_list, '-sequence'=>"no");
}

#Get annotations - only annotations corresponding to genes in the list supplied should be stored
foreach my $id (@gene_id_list){
  $gene_id_list{$id}=1;
}
$a_ref = &getAnnotations('-annotation_dir'=>$annotation_dir, '-dbh'=>$alexa_dbh, '-database'=>$database, '-types_list'=>\%types_list, '-gene_id_list'=>\@gene_id_list, '-junction_seq_size'=>$junction_seq_size, '-boundary_seq_size'=>$boundary_seq_size);
$alexa_dbh->disconnect();


#Import Entrez gene annotation info
#Entrez ID, EnsEMBL ID, Chromosome band, HUGO name
#Key on EnsEMBL ID
#Generate url to Entrez Database
my $ensg_to_entrez_ref = &importEntrezAnnotations('-entrez_annotations'=>$entrez_annotations);

#Create a hash that maps all features to their gene_ids
my $gene_seq_ref = &mapGenesToFeatures('-a_ref'=>$a_ref, '-types_list'=>\@types_list);


#1.) Get root paths to the data and import genome partition values
my $data_paths_ref = &getDataPaths('-file'=>$paths_file);
my $partitions_ref = &getPartitions('-file'=>$partition_file);

#Assign every annotated object to a partition
my $assigned_partitions = &assignPartitions('-a_ref'=>$a_ref, '-partitions_ref'=>$partitions_ref);
print BLUE, "\n\t\tAssigned $assigned_partitions records to a single genome partition", RESET;
print LOG "\n\t\tAssigned $assigned_partitions records to a single genome partition";

#Get paths refs
my $cutoff_paths_ref = $data_paths_ref->{'Expression'}->{paths};
my $expression_paths_ref = $data_paths_ref->{'Expression'}->{paths};
my $de_paths_ref = $data_paths_ref->{'DifferentialExpression'}->{paths};
my $ae_paths_ref = $data_paths_ref->{'DifferentialSplicing'}->{paths};

#Assign a cutoff value to each gene
my $cutoffs_ref = &assignCutoffs('-paths'=>$cutoff_paths_ref);


#1-B.) Get group IDs from Lib_Names file
my %lib_groups;
my $lib_names_file = $paths_file;
$lib_names_file =~ s/Paths/Names/;
open (LIB_NAMES, $lib_names_file) or die "can't fine $lib_names_file\n";
while(<LIB_NAMES>){
  my @line = split("\t", $_);
  my $lib_name = $line[0];
  my $group = $line[2];
  $lib_groups{$lib_name}=$group;
}
close (LIB_NAMES);


#2-A.) Import expression data - all known or expressed novel features should be imported for this summary
my $expression_data_ref = &getExpressionData('-paths'=>$expression_paths_ref, '-a_ref'=>$a_ref, '-ensembl_version'=>$ensembl_version, '-import_option'=>"all", '-types_list'=>\%types_list);

#Also import the alternative expression data
my $ae_data_ref = &getDifferentialSplicingData('-paths'=>$ae_paths_ref, '-a_ref'=>$a_ref, '-ensembl_version'=>$ensembl_version, '-import_option'=>"filtered", '-types_list'=>\%types_list);


my $type_count = keys %{$expression_data_ref};
print BLUE, "\n\nImported expression data for $type_count data types", RESET;
print LOG "\n\nImported expression data for $type_count data types";

my $genes_ref = $a_ref->{'Gene'};

print BLUE, "\n\nCreating gene-by-gene pages", RESET;
print LOG "\n\nCreating gene-by-gene pages";

my %expression_summary;
foreach my $gene_id (sort {$a <=> $b} @gene_id_list){
  print BLUE, "\n\tProcessing gene: $gene_id", RESET;
  print LOG "\n\tProcessing gene: $gene_id";

  my $partition = $genes_ref->{$gene_id}->{partition};

  my $gene_content = '';

  #3.) Basic info from ALEXA
  my $ensembl_g_id = $genes_ref->{$gene_id}->{ensembl_g_id};
  my $evidence = $genes_ref->{$gene_id}->{evidence};
  my $gene_name = $genes_ref->{$gene_id}->{gene_name};
  my $gene_type = $genes_ref->{$gene_id}->{gene_type};
  my $description = $genes_ref->{$gene_id}->{description};
  my $strand = $genes_ref->{$gene_id}->{chr_strand};
  if ($strand eq "1"){
    $strand = "+";
  }else{
    $strand = "-";
  }
  my $chromosome = $genes_ref->{$gene_id}->{chromosome};
  my $start = $genes_ref->{$gene_id}->{chr_start};
  my $end = $genes_ref->{$gene_id}->{chr_end};
  my $size = $genes_ref->{$gene_id}->{seq_length};

  #4.) Basic Entrez info
  my $band = "N/A";
  my $hugo = "N/A";
  my $entrez_link = "N/A";
  if ($ensg_to_entrez_ref->{$ensembl_g_id}){
    $band = join(" ", @{$ensg_to_entrez_ref->{$ensembl_g_id}->{band}});
    $hugo = join (" ", @{$ensg_to_entrez_ref->{$ensembl_g_id}->{hugo}});
    $entrez_link = join (" ", @{$ensg_to_entrez_ref->{$ensembl_g_id}->{url}});
  }

  #Create links to AspAlt Database
  my $link1 = "<A HREF=\"http://66.170.16.154/cgi-bin/aspalt/search.pl?select2="."$species"."&select3=Gene Symbol&&textfield1=$gene_name&search_next=1\" TARGET=\"_blank\" TITLE=\"Link to AspAlt\">$gene_name</A>";
  my $link2 = "<A HREF=\"http://66.170.16.154/cgi-bin/aspalt/search.pl?select2="."$species"."&select3=Ensembl_ID&&textfield1=$ensembl_g_id&search_next=1\" TARGET=\"_blank\" TITLE=\"Link to AspAlt\">$ensembl_g_id</A>";
  my $asp_alt_content = "AspAlt Database Records: $link1 | $link2";

  #Create link to EnsEMBL database
  my $species_formated = $species;
  $species_formated =~ s/ /\_/g;
  my $ensembl_link = "<A HREF=\"http://www.ensembl.org/$species_formated/Gene/Summary?g=$ensembl_g_id\" TITLE=\"Link to Ensembl gene record\">$ensembl_g_id</A>";

  #Store a list of libraries for later use
  my %library_list;
  foreach my $library (sort {$expression_paths_ref->{$a}->{line_order} <=> $expression_paths_ref->{$b}->{line_order}} keys %{$expression_paths_ref}){
    my $name = $expression_paths_ref->{$library}->{name};
    $library_list{$library}{name} = $name;
    $library_list{$library}{line_order} = $expression_paths_ref->{$library}->{line_order};
  }

  #5.) Extra info
  #5-A.) Generate links for each library - Also get the grand read count for each library for this gene
  my $links_content = &generateLinksContent('-gene_id'=>$gene_id);

  #5-B.) Get a count of all features - and create simple HTML object to summarize them
  my $feature_count_content = &generateFeatureCountContent('-a_ref'=>$a_ref, '-gene_id'=>$gene_id, '-gene_seq_ref'=>$gene_seq_ref, '-types_list'=>\@types_list);

  #5-C.) Create a simple HTML object to summarize transcript specific features
  my $trans_specific_content = &generateTranscriptSpecificContent('-a_ref'=>$a_ref, '-gene_id'=>$gene_id, '-gene_seq_ref'=>$gene_seq_ref, '-gene_name'=>$gene_name, '-ensembl_g_id'=>$ensembl_g_id);

  #7.) Expression summary
  #Create a table summarizing the expression level of all expressed features
  #This table should show all known features and any novel features that were expressed
  #It should be sorted by chromosome position (ascending chr_start for genes on positive strand and descending on chr_end for genes on -ve strand)
  #Display genes, then transcripts, then all other features
  #Display columns as follows:
  #Feature ID, Feature Name, Feature Type, Base Count (size), Sequence Supported (y/n), LibraryData {Expression level (X), Expressed (y/n),  Rank (within that feature type)}

  #7-A.) Build an expression object for this gene.  Key on data type and for library specific values store as hashes
  %expression_summary = (); #Summary of number of known expressed features (exon regions and exon-junctions) for each library
  my $gene_expression_ref = &buildExpressionObject('-gene_id'=>$gene_id, '-a_ref'=>$a_ref, '-gene_seq_ref'=>$gene_seq_ref, '-types_list'=>\@types_list, 
                                                   '-expression_data_ref'=>$expression_data_ref, '-ae_data_ref'=>$ae_data_ref, '-strand'=>$strand);

  #7-B.) Similarly, build a model object for this gene to allow a simple gene model to be drawn
  my @model_types_list = qw (ExonRegion Junction Intron Boundary ActiveIntronRegion);
  &buildGeneModel('-gene_id'=>$gene_id, '-a_ref'=>$a_ref, '-gene_seq_ref'=>$gene_seq_ref, '-types_list'=>\@model_types_list, 
                  '-output_file'=>$gene_model_file, '-ae_data_ref'=>$ae_data_ref);

  #Get the grand CDS start/end positions (if 'na', convert to 0)
  my $grand_cds_start = $gene_transcripts_ref->{$gene_id}->{grand_cds_start};
  my $grand_cds_end = $gene_transcripts_ref->{$gene_id}->{grand_cds_end};
  unless ($grand_cds_start =~ /\d+/ && $grand_cds_end =~ /\d+/){
    $grand_cds_start=0;
    $grand_cds_end=0;
  }

  #7-C.) Create a simple HTML object to summarize the number of known features and expressed known features for each library
  my $expression_summary_content = '';
  $expression_summary_content .= "\n<P CLASS=\"Indented12LR_s16_Bold\">Summary of known expressed features</P>\n";
  $expression_summary_content .= "\n<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[box6h]\" data-openimage=\"../../images/Minus_icon.gif\" data-closedimage=\"../../images/Plus_icon.gif\"><img src=\"../../images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $expression_summary_content .= "\n<div id=\"box6h\">";

  #Summary of features expressed above background
  $expression_summary_content .= "\n<P CLASS=\"Indented24LR_s16_Bold\">Expressed above background</P>\n";
  $expression_summary_content .= "<TABLE CLASS=\"Text2\">\n";
  $expression_summary_content .= "<TR><TD CLASS=\"Data8\">Library</TD><TD CLASS=\"Data8\">Exon Regions<BR>(Expressed/Total)</TD><TD CLASS=\"Data8\">Known Junctions<BR>(Expressed/Total)</TD></TR>\n";
  foreach my $library (sort {$library_list{$a}->{line_order} <=> $library_list{$b}->{line_order}} keys %library_list){
    my $name = $library_list{$library}{name};
    $expression_summary_content .= "<TR><TD CLASS=\"Data7\">$name:</TD><TD CLASS=\"Data7\">$expression_summary{$library}{expressed_exon_regions}/$expression_summary{$library}{exon_regions}</TD><TD CLASS=\"Data7\">$expression_summary{$library}{expressed_known_junctions}/$expression_summary{$library}{known_junctions}</TD></TR>\n";
  }
  $expression_summary_content .= "</TABLE><BR>\n";


  #Summary of features expressed above 0
  $expression_summary_content .= "\n<P CLASS=\"Indented24LR_s16_Bold\">Detected at any level above 0</P>\n";
  $expression_summary_content .= "<TABLE CLASS=\"Text2\">\n";
  $expression_summary_content .= "<TR><TD CLASS=\"Data8\">Library</TD><TD CLASS=\"Data8\">Exon Regions<BR>(Detected/Total)</TD><TD CLASS=\"Data8\">Known Junctions<BR>(Detected/Total)</TD></TR>\n";
  foreach my $library (sort {$library_list{$a}->{line_order} <=> $library_list{$b}->{line_order}} keys %library_list){
    my $name = $library_list{$library}{name};
    $expression_summary_content .= "<TR><TD CLASS=\"Data7\">$name:</TD><TD CLASS=\"Data7\">$expression_summary{$library}{detected_exon_regions}/$expression_summary{$library}{exon_regions}</TD><TD CLASS=\"Data7\">$expression_summary{$library}{detected_known_junctions}/$expression_summary{$library}{known_junctions}</TD></TR>\n";
  }
  $expression_summary_content .= "</TABLE><BR>\n</div><BR>\n";


  #7-D.) Generate an HTML table to summarize the expression data in more detail
  my $expression_table_content = &generateExpressionTable('-a_ref'=>$a_ref, '-gene_expression_ref'=>$gene_expression_ref, '-types_list'=>\@types_list, '-library_list'=>\%library_list, '-gene_name'=>$gene_name);

  &createGeneDataFile('-expression_data_ref'=>$expression_data_ref, '-library_list'=>\%library_list, '-temp_dir'=>$temp_dir);


  #Determine where the webpage for the current gene will be written to... define $web_path, $images_dir and $gene_data_dir
  my $partition_dir = "$temp_dir"."$partition"."/";
  unless(-e $partition_dir){
    system("mkdir $partition_dir");
  }
  my $images_dir = "$partition_dir"."images/";
  unless(-e $images_dir){
    system("mkdir $images_dir");
  }
  my $gene_data_dir = "$partition_dir"."data/";
  unless(-e $gene_data_dir){
    system("mkdir $gene_data_dir");
  }
  my $gene_data_path = "$gene_data_dir"."$ensembl_g_id".".txt";
  my $web_path = "$partition_dir"."$ensembl_g_id".".htm";

  my $expression_figures_content = &generateExpressionFigure('-a_ref'=>$a_ref, '-gene_expression_ref'=>$gene_expression_ref, '-types_list'=>\@types_list, '-library_list'=>\%library_list, 
                                                             '-images_dir'=>$images_dir, '-all_genes_data_file'=>$gene_data_file, '-single_gene_data_file'=>$single_gene_data_file,
                                                             '-gene_id'=>$gene_id, '-ensg_id'=>$ensembl_g_id, '-gene_name'=>$gene_name, '-cutoffs_ref'=>$cutoffs_ref,
                                                             '-grand_cds_start'=>$grand_cds_start, '-grand_cds_end'=>$grand_cds_end);
 


  #Gather the html content into a single object
  $gene_content = <<"EOF";
<!-- Start toggle All button section -->
<P CLASS=\"Indented12LR\"><a href="javascript:animatedcollapse.show(['box1s','box2h','box3h','box4h','box5h','box6h','box7h','box8h','leg1h','leg2h'])"><IMG SRC="../../images/Plus_icon.gif" CLASS="Pic_unbordered"></a> | <a href="javascript:animatedcollapse.hide(['box1s','box2h','box3h','box4h','box5h','box6h','box7h','box8h','leg1h','leg2h'])"><IMG SRC="../../images/Minus_icon.gif" CLASS="Pic_unbordered"></a> Show/Hide All Content</P>
<!-- End toggle All button section -->



<!-- Start gene stats section -->
<BR>
<P CLASS=\"Indented12LR_s16_Bold\">Basic Gene Stats</P>
<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[box1s]\" data-openimage=\"../../images/Minus_icon.gif\" data-closedimage=\"../../images/Plus_icon.gif\"><img src=\"../../images/Minus_icon.gif\" border=\"0\" /></a></P>
<div id=\"box1s\">
<P CLASS=\"Indented24LR_s16\">
Species: $species_name<BR>
Gene Name: '<i>$gene_name</i>' (HUGO: $hugo)<BR>
ALEXA Gene ID: $gene_id ($database)<BR>
EnsEMBL Gene ID: $ensembl_g_id<BR>
Entrez Gene Record(s): $entrez_link<BR>
Ensembl Gene Record: $ensembl_link<BR>
Evidence: $evidence<BR>
Gene Type: $gene_type<BR>
Location: chr$chromosome $start-$end ($strand): $band<BR>
Size (bp): $size<BR>
Description: $description<BR>
</P>
<BR>
</div><BR>
<!-- End gene stats section -->

<!-- Start of data links section -->
$links_content
<!-- End of data links section -->

<!-- Start of links to neighbourhood section for $gene_name ($ensembl_g_id) -->
<P CLASS=\"Indented12LR_s16_Bold\">Link to other genes in the same chromosome region as '<i>$gene_name</i>'</P>
<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[box3h]\" data-openimage=\"../../images/Minus_icon.gif\" data-closedimage=\"../../images/Plus_icon.gif\"><img src=\"../../images/Plus_icon.gif\" border=\"0\" /></a></P>
<div id=\"box3h\">
<P CLASS=\"Indented24LR_s16\"><A HREF=\"index.html\" TITLE=\"Link to nearby genes on the same chromosome\">$partition</A></P>
<BR>
</div><BR>
<!-- End of links to neighbourhood section for $gene_name ($ensembl_g_id) -->

<!-- Start of feature counts section for $gene_name ($ensembl_g_id) -->
$feature_count_content
<!--End of feature counts section -->

<!-- Start of transcript specific features section for $gene_name ($ensembl_g_id) -->
$trans_specific_content
<!-- End of feature counts section -->

<!-- Start of summary of expressed known features section for $gene_name ($ensembl_g_id) -->
$expression_summary_content 
<!-- End of summary of expressed known features section -->

<!-- Start of expression table section for $gene_name ($ensembl_g_id) -->
$expression_table_content
<!-- End of expression table section -->

<!-- Start of figures section for $gene_name ($ensembl_g_id) -->
$expression_figures_content
<!-- End of figures section -->
<BR>
<P CLASS=\"Indented12LR_s16_bold\">Download gene data file ($gene_name): <A HREF=\"data/$ensembl_g_id\.txt\">$ensembl_g_id\.txt</A></P><BR>
<BR>
EOF


  my $title = "Summary page for '$gene_name' ($ensembl_g_id) - Project: $project_name";
  my $meta_description = "Provides general gene info and summarizes expression data for all features of the gene $gene_name (aka $ensembl_g_id, $hugo): $description";
  my $meta_keywords = "$gene_name $ensembl_g_id, $hugo, $description";

  #Write out the page
  &writePage('-path'=>$web_path, '-title'=>$title, '-content'=>\$gene_content, '-css_path'=>"../../../ALEXA2.css", '-alexa_home_path'=>"$alexa_home_path", '-alexa_seq_home_path'=>"$alexa_seq_path", '-summary_path'=>"../../Summary.htm", '-genes_path'=>"../index.html", '-search_path'=>"$search_page_url", '-meta_description'=>$meta_description, '-meta_keywords'=>$meta_keywords, '-google_analytics'=>0, '-google_analytics_id'=>$google_analytics_id, '-collapse_div_script'=>'../../../animatedcollapse.js', '-jquery_min_script'=>'../../../jquery.min.js', '-div_count'=>9);

  #copy gene data file to web folder
  system ("cp $single_gene_data_file $gene_data_path");
  #system ("gzip $gene_data_path");
}

#Remove the temp files
if ($test){
  print YELLOW, "\n\nTest R script with the following R variables:\n\n$r_cmd\n\n", RESET;
}else{
  system ("rm -f $gene_data_file");
  system ("rm -f $single_gene_data_file");
  system ("rm -f $gene_model_file");
}

print MAGENTA, "\nWebsite files were written to $temp_dir\n\n", RESET;
print "SCRIPT COMPLETE\n\n";

print LOG "\nWebsite files were written to $temp_dir\n\n";
print LOG "SCRIPT COMPLETE\n\n";

close(LOG);
exit();


###########################################################################################################################
#importEntrezAnnotations                                                                                                  #
###########################################################################################################################
sub importEntrezAnnotations{
  my %args = @_;
  my $file = $args{'-entrez_annotations'};

  my %ensg_to_entrez;

  open (ENTREZ, "$file") || die "\nCould not open Entrez file: $file\n\n";

  while(<ENTREZ>){
    chomp($_);
    my @line = split("\t", $_);
    my $id = $line[1];
    my $entrez_name = $line[2];
    my $external_ids = $line[5];
    my $band = $line[7];
    my $hugo_name = $line[10];

    if ($external_ids =~ /Ensembl\:(\w+)/){
      my $ensembl_id = $1;

      #Craft link
      my $line_base = "http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=full_report&list_uids=";
      my $link_tmp = "$line_base"."$id";
      my $link = "<A HREF=\"$link_tmp\" TITLE=\"Link to Entrez gene record\">$entrez_name</A>";

      if ($ensg_to_entrez{$ensembl_id}){
        push(@{$ensg_to_entrez{$ensembl_id}{entrez_id}}, $id);
        push(@{$ensg_to_entrez{$ensembl_id}{band}}, $band);
        push(@{$ensg_to_entrez{$ensembl_id}{hugo}}, $hugo_name);
        push(@{$ensg_to_entrez{$ensembl_id}{url}}, $link);
      }else{
        my @ids;
        push(@ids, $id);
        my @bands;
        push(@bands, $band);
        my @hugos;
        push(@hugos, $hugo_name);
        my @links;
        push(@links, $link);
        $ensg_to_entrez{$ensembl_id}{entrez_id} = \@ids;
        $ensg_to_entrez{$ensembl_id}{band} = \@bands;
        $ensg_to_entrez{$ensembl_id}{hugo} = \@hugos;
        $ensg_to_entrez{$ensembl_id}{url} = \@links;
      }
    }
  }

  close(ENTREZ);

  return(\%ensg_to_entrez);
}


###########################################################################################################################
#mapGenesToFeatures                                                                                                       #
###########################################################################################################################
sub mapGenesToFeatures{
  my %args = @_;
  my $a_ref = $args{'-a_ref'};

  my %maps;

  print BLUE, "\n\nCreating gene to feature type mappings", RESET;
  print LOG "\n\nCreating gene to feature type mappings";

  my $c = 0;

  foreach my $data_type (sort {$a cmp $b} keys %{$a_ref}){
    print BLUE, "\n\tProcessing: $data_type", RESET;
    print LOG "\n\tProcessing: $data_type";

    my $seq_a_ref = $a_ref->{$data_type};

    foreach my $id (keys %{$seq_a_ref}){

      #Remember that certain types of features can be associated with multiple gene IDs (introns and intergenic regions)
      #Create the one-to-many associations so that each gene record will display all relevant introns/intergenics even if they are shared with other genes
      my $gene_id_list = $seq_a_ref->{$id}->{gene_id};

      unless ($gene_id_list){
        print RED, "\nInvalid gene string in mapGenesToFeatures()!", RESET;
        print LOG "\nInvalid gene string in mapGenesToFeatures()!";
        exit();
      }

      my @gene_ids = split(" ", $gene_id_list);

      foreach my $gene_id (@gene_ids){
        if ($maps{$data_type}{$gene_id}){
          $c++;
          push(@{$maps{$data_type}{$gene_id}}, $id)
        }else{
          $c++;
          my @ids;
          push(@ids, $id);
          $maps{$data_type}{$gene_id} = \@ids;
        }
      }
    }
  }
  print BLUE, "\n\n\tStored a total of $c gene-to-feature mappings", RESET;
  print LOG "\n\n\tStored a total of $c gene-to-feature mappings";
  return(\%maps);
}


###########################################################################################################################
#assign expression cutoff values for each gene                                                                            #
###########################################################################################################################
sub assignCutoffs{
  my %args = @_;
  my $paths_ref = $args{'-paths'};

  my %cutoffs;

  foreach my $library (sort {$a cmp $b} keys %{$paths_ref}){

    my $data_dir = $paths_ref->{$library}->{stats_path};
    unless ($data_dir =~ /.*\/$/){
      $data_dir = "$data_dir"."/";
    }
    my $cutoffs_file = "$data_dir"."$library"."_NORM1_average_coverage_cutoffs.txt";
    unless(-e $cutoffs_file){
      print RED, "\nCould not find expression cutoffs file for library: $library\n\t$cutoffs_file\n", RESET;
      print LOG "\nCould not find expression cutoffs file for library: $library\n\t$cutoffs_file\n";
      exit();
    }
    my $gene_cutoffs_ref = &importExpressionCutoffs ('-cutoffs_file'=>$cutoffs_file);

    foreach my $gene_id (keys %{$gene_cutoffs_ref}){
      $cutoffs{$gene_id}{$library} = $gene_cutoffs_ref->{$gene_id}->{cutoff};
    }
  }

  my $stored_cutoffs = keys %cutoffs;
  print BLUE, "\n\nStored $stored_cutoffs cutoff values", RESET;
  print LOG "\n\nStored $stored_cutoffs cutoff values";
  return(\%cutoffs);
}


###########################################################################################################################
#Generate links for each library - Also get the grand read count for each library for this gene
###########################################################################################################################
sub generateLinksContent{
  my %args = @_;
  my $gene_id = $args{'-gene_id'};
  my $links = '';
  $links .= "<P CLASS=\"Indented12LR_s16_Bold\">Data links for each library (displays expression data in UCSC Genome Browser)</P>\n";
  $links .= "<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[box2h]\" data-openimage=\"../../images/Minus_icon.gif\" data-closedimage=\"../../images/Plus_icon.gif\"><img src=\"../../images/Plus_icon.gif\" border=\"0\" /></a></P>\n";
  $links .= "<div id=\"box2h\">\n";
  $links .= "<P CLASS=\"Indented24LR_s16\">\n";

  my $gene_name = $genes_ref->{$gene_id}->{gene_name};
  my $chromosome = $genes_ref->{$gene_id}->{chromosome};
  my $partition = $genes_ref->{$gene_id}->{partition};
  my $start = $genes_ref->{$gene_id}->{chr_start};
  my $end = $genes_ref->{$gene_id}->{chr_end};
  my $display_start = $start-100;
  my $display_end = $end+100;

  #Create a link for each comparison for each library...
  foreach my $comparison (sort {$de_paths_ref->{$a}->{line_order} <=> $de_paths_ref->{$b}->{line_order}} keys %{$de_paths_ref}){
    $links .= "\nComparison: $comparison<BR>";
    my $libA = '';
    my $libB = '';
    if ($comparison =~ /(\w+)\_vs\_(\w+)/){
      $libA = $1;
      $libB = $2;
    }

    foreach my $library (sort {$expression_paths_ref->{$a}->{line_order} <=> $expression_paths_ref->{$b}->{line_order}} keys %{$expression_paths_ref}){
      unless ($library eq $libA || $library eq $libB){
        next();
      }
      my $current_lib = '';
      if ($library eq $libA){$current_lib = $libA;}
      if ($library eq $libB){$current_lib = $libB;}

      my $name = $expression_paths_ref->{$library}->{name};
      my $ucsc_link_clean = "http://genome.ucsc.edu/cgi-bin/hgTracks?db="."$ucsc_build"."&position=chr$chromosome:$display_start-$display_end&hgt.customText=$track_url"."DE/$library/$comparison/"."$partition"."_DE_merged.txt.gz&ctfile_"."$ucsc_build"."=";
      my $ucsc_link_persist = "http://genome.ucsc.edu/cgi-bin/hgTracks?db="."$ucsc_build"."&position=chr$chromosome:$display_start-$display_end&hgt.customText=$track_url"."DE/$library/$comparison/"."$partition"."_DE_merged.txt.gz";


      my $libraries_ref = $expression_data_ref->{'Gene'}->{libraries};
      my $data_ref = $libraries_ref->{$library}->{data};
      my $read_count = $data_ref->{$gene_id}->{read_count};
      $read_count = &commify($read_count);

      $links .= "\n$name ($current_lib) has $read_count total reads for '<i>$gene_name</i>'.  UCSC data links: (<A HREF=\"$ucsc_link_clean\" TARGET=\"_blank\" TITLE=\"Link to UCSC tracks - clears existing tracks\">C</A> | <A HREF=\"$ucsc_link_persist\" TARGET=\"_blank\" TITLE=\"Link to UCSC tracks - existing tracks will persist\">P</A>)<BR>";
    }
    $links .= "\n<BR>";
  }
  $links .= "</P><BR>\n</div><BR>\n";

  return($links);
}


###########################################################################################################################
#Get a count of all features - and create simple HTML object to summarize them                                            #
###########################################################################################################################
sub generateFeatureCountContent{
  my %args = @_;
  my $gene_id = $args{'-gene_id'};
  my $gene_seq_ref = $args{'-gene_seq_ref'};
  my @types_list = @{$args{'-types_list'}};
  my $a_ref = $args{'-a_ref'};

  my $feature_count = '';
  $feature_count .= "\n<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[box4h]\" data-openimage=\"../../images/Minus_icon.gif\" data-closedimage=\"../../images/Plus_icon.gif\"><img src=\"../../images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $feature_count .= "\n<div id=\"box4h\">";
  $feature_count .= "\n<P CLASS=\"Indented24LR_s16\">";

  my $grand_count = 0;
  foreach my $type (@types_list){
    my $seq_a_ref = $a_ref->{$type};

    my $count = 0;
    if ($gene_seq_ref->{$type}->{$gene_id}){
      my @feature_list = @{$gene_seq_ref->{$type}->{$gene_id}};
      $count = scalar(@feature_list);
      $grand_count += $count;
      $feature_count .= "\n$type: $count<BR>";

      #Deal with special cases for Junction and Boundary sequences
      if ($type eq "Junction" || $type eq "Boundary"){
        my $known_count = 0;
        my $novel_count = 0;
        foreach my $id (@feature_list){
          if ($seq_a_ref->{$id}->{known} == 1){
            $known_count++;
          }else{
            $novel_count++;
          }
        }
        my $known_type = "Known"."$type";
        my $novel_type = "Novel"."$type";
        $feature_count .= "\n&nbsp;&nbsp;&nbsp; $known_type: $known_count<BR>";
        $feature_count .= "\n&nbsp;&nbsp;&nbsp; $novel_type: $novel_count<BR>";
      } 
    }
  }
  $feature_count .= "</P><BR>";
  $feature_count .= "\n</div><BR>";
  $grand_count--;
  $feature_count = "\n<P CLASS=\"Indented12LR_s16_Bold\">Features defined for this gene: $grand_count</P>\n"."$feature_count";

  return($feature_count);
}


###########################################################################################################################
#Generate a summary of the transcript specific features of each transcript of this gene                                   #
###########################################################################################################################
sub generateTranscriptSpecificContent{
  my %args = @_;
  my $gene_id = $args{'-gene_id'};
  my $gene_name = $args{'-gene_name'};
  my $ensembl_g_id = $args{'-ensembl_g_id'};
  my $gene_seq_ref = $args{'-gene_seq_ref'};
  my $a_ref = $args{'-a_ref'};
  
  my $content = '';
  $content .= "<P CLASS=\"Indented12LR_s16_Bold\">Summary of transcript specific features for '<i>$gene_name</i>' ($ensembl_g_id)</P>\n";
  $content .= "\n<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[box5h]\" data-openimage=\"../../images/Minus_icon.gif\" data-closedimage=\"../../images/Plus_icon.gif\"><img src=\"../../images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $content .= "\n<div id=\"box5h\">";

  $content .= "<TABLE CLASS=\"Text2\">\n";

  my $type = "Transcript";
  my $seq_a_ref = $a_ref->{$type};

  if ($gene_seq_ref->{$type}->{$gene_id}){
    my @feature_list = @{$gene_seq_ref->{$type}->{$gene_id}};
    foreach my $id (@feature_list){

      #Get the list of transcript specific features
      my $feature_list = $seq_a_ref->{$id}->{feature_list};
      my @temp = split(",", $feature_list);
      my $formatted_list = join (", ", @temp);
      my $seq_name = $seq_a_ref->{$id}->{seq_name};
      $content .= "<TR><TD CLASS=\"Data6\">$seq_name:</TD><TD CLASS=\"Data5\">$formatted_list</TD></TR>\n";
    }
  }
  
  $content .= "</TABLE><BR>\n</div><BR>\n";

  return($content);
}


###########################################################################################################################
#Build an expression object keyed on data type - store library specific values as hashes                                  #
###########################################################################################################################
sub buildExpressionObject{
  my %args = @_;
  my $gene_id = $args{'-gene_id'};
  my $a_ref = $args{'-a_ref'};
  my $gene_seq_ref = $args{'-gene_seq_ref'};
  my @types_list = @{$args{'-types_list'}};
  my $expression_data_ref = $args{'-expression_data_ref'};
  my $ae_data_ref = $args{'-ae_data_ref'};
  my $strand = $args{'-strand'};

  my %expression;

  #Go through each feature of each data type and create an expression record for it if it is known or expressed in at least one library
  my $ec = 0;
  foreach my $type (@types_list){
    my %type_exp;
    my $seq_a_ref = $a_ref->{$type};

    my $libraries_ref = $expression_data_ref->{$type}->{libraries};

    #Get the list of feature ids for this type for this gene
    if ($gene_seq_ref->{$type}->{$gene_id}){
      my @ids = @{$gene_seq_ref->{$type}->{$gene_id}};

      foreach my $id (@ids){
        $ec++;

        #Make sure this element is 'known' or expressed in at least one library.  If not, delete it from the data object
        my $test1 = 0;

        #Get expression data for this id
        foreach my $library (keys %{$libraries_ref}){
          my $test2 = 0;
          my $test3 = 0;

          unless($expression_summary{$library}){
            $expression_summary{$library}{expressed_exon_regions} = 0;
            $expression_summary{$library}{detected_exon_regions} = 0;
            $expression_summary{$library}{exon_regions} = 0;
            $expression_summary{$library}{expressed_known_junctions} = 0;
            $expression_summary{$library}{detected_known_junctions} = 0;
            $expression_summary{$library}{known_junctions} = 0;
          }

          #print YELLOW, "\nType: $type\tID: $id\tLibrary: $library", RESET;
          my $data_ref = $libraries_ref->{$library}->{data};
          if ($data_ref->{$id}){
            #If this feature is defined in the expression object it must be either known or expressed
            my %tmp;
            my $norm_expression = $data_ref->{$id}->{norm_expression};
            my $norm_expression_log2 = $data_ref->{$id}->{norm_expression};
            unless($norm_expression eq "NA"){
              $norm_expression = sprintf("%.2f", $norm_expression);
             
              #Convert to log2 scale (first add 1 - to prevent NaN and -ve values)
              $norm_expression_log2 = (log($norm_expression+1))/(log(2));
              $norm_expression_log2 = sprintf("%.2f", $norm_expression_log2);
            }
            $tmp{norm_expression} = $norm_expression;
            $tmp{norm_expression_log2} = $norm_expression_log2;
            $tmp{expressed} =  $data_ref->{$id}->{expressed};

            $expression{$ec}{id} = $id;
            $expression{$ec}{seq_name} = $seq_a_ref->{$id}->{seq_name};
            $expression{$ec}{known} = $seq_a_ref->{$id}->{known};
            $expression{$ec}{data_type} = $type;
            $expression{$ec}{type_name} = $type;
            $expression{$ec}{is_ae} = 0;

            #Fix type to account for special cases: Junction and Boundary
            if ($type eq "Junction" || $type eq "Boundary"){
              if ($seq_a_ref->{$id}->{known} == 1){
                $expression{$ec}{type_name} = "Known"."$type";
              }else{
                $expression{$ec}{type_name} = "Novel"."$type";
              }
            }

            $expression{$ec}{chr_start} = $seq_a_ref->{$id}->{chr_start};
            $expression{$ec}{chr_end} = $seq_a_ref->{$id}->{chr_end};
            $expression{$ec}{base_count} = $seq_a_ref->{$id}->{base_count};

            #Make sure gene has lower coords than transcripts
            if ($type eq "Gene"){
              $expression{$ec}{chr_start}-=5;
              $expression{$ec}{chr_end}+=5;
            }

            if ($seq_a_ref->{$id}->{known} || $data_ref->{$id}->{expressed}){$test1 = 1;}
            if ($seq_a_ref->{$id}->{known} && $data_ref->{$id}->{expressed}){$test2 = 1;}
            
            unless($data_ref->{$id}->{norm_expression} eq "NA"){
              if ($seq_a_ref->{$id}->{known} && $data_ref->{$id}->{norm_expression} > 0){$test3 = 1;}
            }

            if ($seq_a_ref->{$id}->{known} && $type eq "ExonRegion"){$expression_summary{$library}{exon_regions}++;}
            if ($test2 && $type eq "ExonRegion"){$expression_summary{$library}{expressed_exon_regions}++;}
            if ($test3 && $type eq "ExonRegion"){$expression_summary{$library}{detected_exon_regions}++;}

            if ($seq_a_ref->{$id}->{known} && $type eq "Junction"){$expression_summary{$library}{known_junctions}++;}
            if ($test2 && $type eq "Junction"){$expression_summary{$library}{expressed_known_junctions}++;}
            if ($test3 && $type eq "Junction"){$expression_summary{$library}{detected_known_junctions}++;}

            if ($expression{$ec}{libraries}){
              my $libs_ref = $expression{$ec}{libraries};
              $libs_ref->{$library}->{data} = \%tmp;
            }else{
              my %libs;
              $libs{$library}{data} = \%tmp;
              $expression{$ec}{libraries} = \%libs;
            }
          }
        }
        if ($test1 eq "0" && $expression{$ec}){
          delete($expression{$ec});
        }
      }
    }
  }

  #Add an order to all features based on coordinates position and according to strand
  if ($strand eq "+"){
    my $order = 0;
    foreach my $ec (sort {$expression{$a}{chr_start} <=> $expression{$b}{chr_start} || $expression{$a}{chr_end} <=> $expression{$b}{chr_end}} keys %expression){
      $order++;
      $expression{$ec}{order} = $order;
    }
  }else{
    my $order = 0;
    foreach my $ec (sort {$expression{$b}{chr_end} <=> $expression{$a}{chr_end} || $expression{$b}{chr_start} <=> $expression{$a}{chr_start}} keys %expression){
      $order++;
      $expression{$ec}{order} = $order;
    }
  }

  #Check each feature to see if it was identified as alternatively expressed in at least one comparison - if so mark it
  foreach my $ec (keys %expression){
    my $id = $expression{$ec}{id};
    my $type = $expression{$ec}{data_type};

    foreach my $data_type_comp (keys %{$ae_data_ref}){
      my $comparison = $ae_data_ref->{$data_type_comp}->{comparison};
      my $data_type = $ae_data_ref->{$data_type_comp}->{data_type};
  
      if ($type eq $data_type){
        my $data_ref = $ae_data_ref->{$data_type_comp}->{data};
        if ($data_ref->{$id}){
          $expression{$ec}{is_ae} = 1;
        }
      }
    }
  }

  return(\%expression);
}


###########################################################################################################################
#Build an expression object keyed on data type - store library specific values as hashes                                  #
###########################################################################################################################
sub buildGeneModel{
  my %args = @_;
  my $gene_id = $args{'-gene_id'};
  my $a_ref = $args{'-a_ref'};
  my $gene_seq_ref = $args{'-gene_seq_ref'};
  my @types_list = @{$args{'-types_list'}};
  my $ae_data_ref = $args{'-ae_data_ref'};
  my $output_file = $args{'-output_file'};

  my %gene_model;

  #Go through each feature of each data type used in the gene model and store basic gene model info
  my $gm = 0;
  foreach my $type (@types_list){
    my $seq_a_ref = $a_ref->{$type};

    #Get the list of feature ids for this type for this gene
    if ($gene_seq_ref->{$type}->{$gene_id}){
      my @ids = @{$gene_seq_ref->{$type}->{$gene_id}};

      foreach my $id (@ids){
        $gm++;
        my $model_start = $seq_a_ref->{$id}->{model_start};
        my $model_end = $seq_a_ref->{$id}->{model_end};
        my $feature_name = $seq_a_ref->{$id}->{seq_name};
        my $skipped_exons = "NA";
        if ($type eq "Junction"){
           $skipped_exons = $seq_a_ref->{$id}->{exons_skipped};
        }

        #Shorten names
        if ($feature_name =~ /ER(\w+)/){
          $feature_name = $1;
        }elsif($feature_name =~ /E(\w+)\_E(\w+)/){
          $feature_name = "$1"."_"."$2";
        }elsif($feature_name =~ /I(\w+)/){
          $feature_name = $1;
        }elsif($feature_name =~ /E(\w+)/){
          $feature_name = $1;
        }

        $gene_model{$gm}{id}=$id;
        $gene_model{$gm}{start} = $model_start;
        $gene_model{$gm}{end} = $model_end;
        $gene_model{$gm}{name} = $feature_name;
        $gene_model{$gm}{type} = $type;
        $gene_model{$gm}{is_ae} = 0;
        $gene_model{$gm}{skipped_exons} = $skipped_exons;
        $gene_model{$gm}{order} = 0;
      }
    }
  }

  #If the exon has multiple pieces, keep the a, b, c, ... labeling, but otherwise remove it
  my %exon_counts;
  foreach my $gm (sort {$gene_model{$a}{start} <=> $gene_model{$b}{start} || $gene_model{$a}{end} <=> $gene_model{$b}{end}} keys %gene_model){
    my $type = $gene_model{$gm}{type};
    if ($type eq "ExonRegion"){
      my $feature_name = $gene_model{$gm}{name};
      if ($feature_name =~ /(\d+)\w+/){
        $exon_counts{$1}++;
      }
    }
  }

  foreach my $gm (sort {$gene_model{$a}{start} <=> $gene_model{$b}{start} || $gene_model{$a}{end} <=> $gene_model{$b}{end}} keys %gene_model){
    my $type = $gene_model{$gm}{type};
    if ($type eq "ExonRegion"){
      my $feature_name = $gene_model{$gm}{name};
      if ($feature_name =~ /(\d+)\w+/){
        my $exon_count = $1;
        if ($exon_counts{$exon_count} == 1){
          $gene_model{$gm}{name} = $exon_count;
        }
      }
    }
  }


  #Add an order to all features based on coordinate positions
  my $order = 0;
  foreach my $gm (sort {$gene_model{$a}{start} <=> $gene_model{$b}{start} || $gene_model{$a}{end} <=> $gene_model{$b}{end}} keys %gene_model){
    $order++;
    $gene_model{$gm}{order} = $order;
  }

  #Check each feature to see if it was identified as alternatively expressed in at least one comparison - if so mark it
  foreach my $gm (keys %gene_model){
    my $id = $gene_model{$gm}{id};
    my $type = $gene_model{$gm}{type};

    foreach my $data_type_comp (keys %{$ae_data_ref}){
      my $comparison = $ae_data_ref->{$data_type_comp}->{comparison};
      my $data_type = $ae_data_ref->{$data_type_comp}->{data_type};
  
      if ($type eq $data_type){
        my $data_ref = $ae_data_ref->{$data_type_comp}->{data};
        if ($data_ref->{$id}){
          $gene_model{$gm}{is_ae} = 1;
        }
      }
    }
  }

  #Write gene model to a temp file
  open(MODEL, ">$output_file") || die "\n\nCould not open temp model output file: $output_file\n\n";
  print MODEL "Start\tEnd\tName\tType\tIsAE\tSkippedExons\n";
  foreach my $gm (sort {$gene_model{$a}{order} <=> $gene_model{$b}{order}} keys %gene_model){
    print MODEL "$gene_model{$gm}{start}\t$gene_model{$gm}{end}\t$gene_model{$gm}{name}\t$gene_model{$gm}{type}\t$gene_model{$gm}{is_ae}\t$gene_model{$gm}{skipped_exons}\n";
  }
  close(MODEL);

  return();
}


###########################################################################################################################
#Generate an HTML table to summarize the expression data                                                                  #
###########################################################################################################################
sub generateExpressionTable{
  my %args = @_;
  my $a_ref = $args{'-a_ref'};
  my $gene_expression_ref = $args{'-gene_expression_ref'};
  my @types_list = @{$args{'-types_list'}};
  my %library_list = %{$args{'-library_list'}};
  my $gene_name = $args{'-gene_name'};

  #Display columns as follows:
  #Feature ID, Feature Name, Feature Type, Base Count (size), Sequence Supported (y/n), LibraryData {Expression level (X), Expressed (y/n),  Rank (within that feature type)}

  #Open table
  my $table_content = "\n<P CLASS=\"Indented12LR_s16_Bold\">Normalized expression data for each known or novel expressed feature of '<i>$gene_name</i>'</P>";
  $table_content .= "\n<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[box7h]\" data-openimage=\"../../images/Minus_icon.gif\" data-closedimage=\"../../images/Plus_icon.gif\"><img src=\"../../images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $table_content .= "\n<div id=\"box7h\">";

  $table_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[leg1h]\" data-openimage=\"../../images/QuestionOff_icon.gif\" data-closedimage=\"../../images/QuestionOn_icon.gif\"><img src=\"../../images/QuestionOn_icon.gif\" border=\"0\" /></a></P>";
  $table_content .= "\n<div id=\"leg1h\" style=\"width: 1100px;\">";

  $table_content .= "\n<P CLASS=\"Indented24LR_s16\"> The following table summarizes the expression values of all known features of the gene '<i>$gene_name</i>' in each library.  Known features (i.e. those that correspond to one or more EnsEMBL transcripts) are included regardless of their expression level. Novel features are only reported in this table if they are expressed above the level of intergenic background noise.  The 'FID' column reports the unique Feature ID for the gene, transcript, exon region, exon junction, etc.  The 'Name' column reports the Feature Name for the feature (e.g. E5a_E6a is the name of an exon junction corresponding to the connection of exons 5 and 6).  The 'Type' column reports the type of the feature and may be one of the following: Gene, Transcript, Exon Region, Known Junction, Novel Junction, Known Boundary, Novel Boundary, Silent or Active Intronic Regions, and Silent or Active Intergenic Regions (see our manuscript for more details).  The 'Base Count' column reports the number of nucleotide bases for the feature followed by the percent of these that do NOT correspond to repeats and the percent that are known to be protein coding in one or more EnsEMBL transcripts.  The 'Supporting ESTs/mRNAs' column reports the number of ESTs and mRNAs that when aligned to the human genome support the expression of a feature (i.e. the sequence alignment coordinates match an exon junction or boundaries of an exon, retained intron, etc.). The 'Conserved Species' column reports the number of 'other' species for which there was at least 1 supporting EST/mRNA alignment.  Finally, the remaining columns report the log2 expression value for each feature in each library. Bold values indicate expression above background.  Each expression value is also a hyperlink to a view of the feature's corresponding genomic region in the UCSC genome browser. Expression data will automatically be loaded as custom GFF and wiggle tracks.</P><BR>";
  $table_content .= "\n</div>";

  $table_content .= "\n<TABLE CLASS=\"Data2\">";

  #Print header row
  $table_content .= "\n  <TR><TD CLASS=\"Head3\">FID</TD><TD CLASS=\"Head3\">Name</TD><TD CLASS=\"Head3\">Type</TD><TD CLASS=\"Head3\">Base Count <BR> (% unmasked |<BR> % coding)</TD><TD CLASS=\"Head3\">Supporting<BR>$species_name<BR>ESTs/mRNAs</TD><TD CLASS=\"Head3\">Conserved<BR>Species</TD>";
  foreach my $library (sort {$library_list{$a}->{line_order} <=> $library_list{$b}->{line_order}} keys %library_list){
    my $name = $library_list{$library}{name};
    $table_content .= "<TD CLASS=\"Head3\">$name Log2<BR>Expression</TD>";
  }
  $table_content .= "\n  </TR>";

  #Print gene and transcipt feature rows - sorted
  #Print all other feature rows
  my @process_groups;
  push(@process_groups, "ExonRegion Junction Boundary Intron ActiveIntronRegion SilentIntronRegion Intergenic ActiveIntergenicRegion SilentIntergenicRegion");
  push(@process_groups, "Gene Transcript");

  foreach my $process_group (@process_groups){
    foreach my $ec (sort {$gene_expression_ref->{$a}->{order} <=> $gene_expression_ref->{$b}->{order}} keys %{$gene_expression_ref}){
      my $id = $gene_expression_ref->{$ec}->{id};
      my $seq_name = $gene_expression_ref->{$ec}->{seq_name};
      my $type = $gene_expression_ref->{$ec}->{data_type};
      my $type_name = $gene_expression_ref->{$ec}->{type_name};
      if ($process_group =~ /$type/){
        next();
      }
      my $seq_a_ref = $a_ref->{$type};
      my $fid = $seq_a_ref->{$id}->{fid};
      my $display_start = ($seq_a_ref->{$id}->{chr_start})-100;
      my $display_end = ($seq_a_ref->{$id}->{chr_end})+100;
      my $chromosome = $seq_a_ref->{$id}->{chromosome};
      if ($chromosome eq "MT"){$chromosome = "M";}
      my $ucsc_chromosome = "chr"."$chromosome";
      my $file_name_prefix = $seq_a_ref->{$id}->{partition};

      my $base_count = $seq_a_ref->{$id}->{base_count};
      my $unmasked_base_count = $seq_a_ref->{$id}->{unmasked_base_count};
      my $coding_base_count = $seq_a_ref->{$id}->{coding_base_count};
      my $unmasked_base_count_p = "N/A";
      my $coding_base_count_p = "N/A";
      if (($base_count =~ /\d+/) && ($unmasked_base_count =~ /\d+/) && ($coding_base_count =~ /\d+/)){
        $unmasked_base_count_p = sprintf("%.0f", (($unmasked_base_count/$base_count)*100));
        $coding_base_count_p = sprintf("%.0f", (($coding_base_count/$base_count)*100));
      }else{
        $base_count = "N/A";
      }

      my $supporting_seqs = $seq_a_ref->{$id}->{supporting_seqs};
      my $conserved_species_count = $seq_a_ref->{$id}->{conserved_species_count};

      my $class = "Data2";
      $table_content .= "\n  <TR><TD CLASS=\"Data2\">$fid</TD><TD CLASS=\"Data3\">$seq_name</TD><TD CLASS=\"Data1\">$type_name</TD><TD CLASS=\"Data2\">$base_count ($unmasked_base_count_p% | $coding_base_count_p%)</TD><TD CLASS=\"Data2\">$supporting_seqs</TD><TD CLASS=\"Data2\">$conserved_species_count</TD>";
      my $libs_ref = $gene_expression_ref->{$ec}->{libraries};
      foreach my $library (sort {$library_list{$a}->{line_order} <=> $library_list{$b}->{line_order}} keys %library_list){
        #Craft a link to the UCSC browser using the coordinates of the current sequence feature
        my $ucsc_link_clean = "http://genome.ucsc.edu/cgi-bin/hgTracks?db="."$ucsc_build"."&position=chr$chromosome:$display_start-$display_end&hgt.customText=$track_url"."$library/combined/"."$file_name_prefix"."_merged.txt.gz&ctfile_"."$ucsc_build=";
        my $ucsc_link_persist = "http://genome.ucsc.edu/cgi-bin/hgTracks?db="."$ucsc_build"."&position=chr$chromosome:$display_start-$display_end&hgt.customText=$track_url"."$library/combined/"."$file_name_prefix"."_merged.txt.gz";
        if ($libs_ref->{$library}){
          my $data_ref = $libs_ref->{$library}->{data};
          my $exp_data = "$data_ref->{norm_expression_log2} (<A HREF=\"$ucsc_link_clean\" TARGET=\"_blank\" TITLE=\"Link to UCSC tracks - clears existing tracks\">C</A> | <A HREF=\"$ucsc_link_persist\" TARGET=\"_blank\" TITLE=\"Link to UCSC tracks - existing tracks will persist\">P</A>)";
          unless($data_ref->{expressed} eq "NA"){
            if ($data_ref->{expressed} == 1){
              $class = "Data4";
            }
          }
          $table_content .= "<TD CLASS=\"$class\">$exp_data</TD>";
        }else{
          $table_content .= "<TD CLASS=\"$class\">N/A</TD>";
        }
      }
      $table_content .= "\n  </TR>";
    }
  }
  $table_content .= "\n</TABLE><BR>\n</div><BR>";
  return($table_content);
}


###########################################################################################################################
#Generate a data file containing gene expression values for all libraries                                                 #
###########################################################################################################################
sub createGeneDataFile{
  my %args = @_;
  my $expression_data_ref = $args{'-expression_data_ref'};
  my %library_list = %{$args{'-library_list'}};
  my $temp_dir = $args{'-temp_dir'};

  #Print out temp file to be imported into R
  #File. contains expression data for all expressed 'Gene' features - all non NA and non 0 values ...
  # - Single column with norm expression data values (one column per library)

  #Generate temp File1.
  my $libraries_ref = $expression_data_ref->{'Gene'}->{libraries};

  my %gene_data;
  foreach my $library (keys %{$libraries_ref}){
    $library_list{$library}{tmp} = '';
    my $data_ref = $libraries_ref->{$library}->{data};
    foreach my $id (keys %{$data_ref}){
      if ($data_ref->{$id}->{norm_expression} eq "NA"){
        next();
      }
      if ($data_ref->{$id}->{norm_expression} == 0){
        next();
      }
      $gene_data{$id}{$library} = $data_ref->{$id}->{norm_expression};
    }
  }
  #Now fill empty spaces with 'NA'
  foreach my $id(keys %gene_data){
    foreach my $library (keys %library_list){
      unless($gene_data{$id}{$library}){
        $gene_data{$id}{$library} = "NA";
      }
    }
  }
  open(OUT, ">$gene_data_file") || die "\nCould not open gene data output file: $gene_data_file\n\n";

  #Print the header
  my @header;
  push(@header, "GeneID");
  foreach my $library (sort {$library_list{$a}->{line_order} <=> $library_list{$b}->{line_order}} keys %library_list){
    push(@header, $library_list{$library}{name});
  }
  my $string = join("\t", @header);
  print OUT "$string\n";

  foreach my $id (keys %gene_data){
    my @data;
    push(@data, $id);
    foreach my $library (sort {$library_list{$a}->{line_order} <=> $library_list{$b}->{line_order}} keys %library_list){
      push(@data, $gene_data{$id}{$library});
    }
    my $string = join("\t", @data);
    print OUT "$string\n";
  }
  close(OUT);

  return();
}


###########################################################################################################################
#Generate a figure to display expression data                                                                             #
###########################################################################################################################
sub generateExpressionFigure{
  my %args = @_;
  my $a_ref = $args{'-a_ref'};
  my $gene_expression_ref = $args{'-gene_expression_ref'};
  my @types_list = @{$args{'-types_list'}};
  my %library_list = %{$args{'-library_list'}};
  my $images_dir = $args{'-images_dir'};
  my $all_genes_data_file = $args{'-all_genes_data_file'};
  my $single_gene_data_file = $args{'-single_gene_data_file'};
  my $gene_id = $args{'-gene_id'};
  my $ensg_id = $args{'-ensg_id'};
  my $gene_name = $args{'-gene_name'};
  my $cutoffs_ref = $args{'-cutoffs_ref'};
  my $grand_cds_start = $args{'-grand_cds_start'};
  my $grand_cds_end = $args{'-grand_cds_end'};

  #Print gene and transcipt feature rows - sorted
  #Print all other feature rows

  #File. contains expression data for all features of a single gene
  # - order, seq_name, seq_type, norm expression data value (one column per library)

  open (OUT, ">$single_gene_data_file") || die "\nCould not open output file: single_gene_data_file\n\n"; 

  #Get cutoff values, library names, and group ids
  my @library_names;
  my @intragenic_cutoffs;
  my @intergenic_cutoffs;
  my @group_ids;

  foreach my $library (sort {$library_list{$a}->{line_order} <=> $library_list{$b}->{line_order}} keys %library_list){
    push(@library_names, $library_list{$library}{name});
    push(@intragenic_cutoffs, $cutoffs_ref->{$gene_id}->{$library});
    push(@intergenic_cutoffs, $cutoffs_ref->{'0'}->{$library});
    push(@group_ids, $lib_groups{$library});
  }
  my $lib_string = join("\t", @library_names);
  my $lib_count = scalar(@library_names);

  #Print the header
  #Temp file to be fed to R will contain: display order, feature name, feature type, is AE?, expression values for each library
  print OUT "Order\tSeqName\tSeqType\tIsAE\t$lib_string\n";

  #Make sure the cutoff values were retrieved successfully
  unless($lib_count == scalar(@intragenic_cutoffs) && $lib_count == scalar(@intergenic_cutoffs)){
    print RED, "\n\nCould not find the expected number of cutoff values for gene: $gene_id\nIntragenic cutoffs: @intragenic_cutoffs\nIntergenic cutoffs: @intergenic_cutoffs\n\n", RESET;
    print LOG "\n\nCould not find the expected number of cutoff values for gene: $gene_id\nIntragenic cutoffs: @intragenic_cutoffs\nIntergenic cutoffs: @intergenic_cutoffs\n\n";
    exit();
  }

  my $intragenic_cutoffs_string = join(" ", @intragenic_cutoffs);
  my $intergenic_cutoffs_string = join(" ", @intergenic_cutoffs);
  my $groups_string = join(" ", @group_ids);

  #Generate temp File.
  my @process_groups;
  push(@process_groups, "ExonRegion Junction Boundary Intron ActiveIntronRegion SilentIntronRegion Intergenic ActiveIntergenicRegion SilentIntergenicRegion");
  push(@process_groups, "Gene Transcript");

  my $new_order = 0;
  my $expressed_data_count = 0;
  foreach my $process_group (@process_groups){

    foreach my $ec (sort {$gene_expression_ref->{$a}->{order} <=> $gene_expression_ref->{$b}->{order}} keys %{$gene_expression_ref}){
      my $type = $gene_expression_ref->{$ec}->{data_type};
      if ($process_group =~ /$type/){
        next();
      }

      my @data;
      $new_order++;

      my $seq_name = $gene_expression_ref->{$ec}->{seq_name};

      #Give transcripts shorter names - remove 'ENST' and padding 0's
      if ($type eq "Transcript"){
        if ($seq_name =~ /^[a-z]+[0]+(\d+)/i){
          $seq_name = $1;
        }
      }
      push(@data, $new_order);
      push(@data, $seq_name);
      push(@data, $gene_expression_ref->{$ec}->{type_name});
      push(@data, $gene_expression_ref->{$ec}->{is_ae});

      my $libs_ref = $gene_expression_ref->{$ec}->{libraries};
      foreach my $library (sort {$library_list{$a}->{line_order} <=> $library_list{$b}->{line_order}} keys %library_list){
        my $data_ref = $libs_ref->{$library}->{data};
        my $norm_expression = $data_ref->{norm_expression};
        if ($norm_expression eq "NA"){
          $norm_expression = 0;
        }
        push(@data, $norm_expression);

        if ($norm_expression > 0){
          $expressed_data_count++;
        }
      }
      my $data_string = join("\t", @data);
      print OUT "$data_string\n";
    }
  }
  close(OUT);

  my $image_desc = "Figures displaying gene and feature expression levels in each library";
  my $figure_content = "\n<P CLASS=\"Indented12LR_s16_Bold\">$image_desc</P>";
  $figure_content .= "\n<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[box8h]\" data-openimage=\"../../images/Minus_icon.gif\" data-closedimage=\"../../images/Plus_icon.gif\"><img src=\"../../images/Plus_icon.gif\" border=\"0\" /></a></P>";
  $figure_content .= "<div id=\"box8h\">";

  #Define the figure height and width (this should match what was used in the R script)
  my $image_width = 1000;
  my $image_height = (400*3);
  my $image_height2 = 600; #heatmap

  #Only generate a graph if some feature was expressed in some library
  if ($expressed_data_count > 0){

    #Using SVGs
    #Run R script
    my $r_script = "$script_dir"."/R_bin/geneFigure_SVG.2.R";

    #$r_cmd = "$r_script $all_genes_data_file $single_gene_data_file $gene_model_file \"$gene_name\" \"$ensg_id\" $gene_id $images_dir \"$intragenic_cutoffs_string\" \"$intergenic_cutoffs_string\" \"$groups_string\" $grand_cds_start $grand_cds_end";
    $r_cmd = "$r_script $all_genes_data_file $single_gene_data_file $gene_model_file \"$gene_name\" \"$ensg_id\" $gene_id $images_dir \"$intragenic_cutoffs_string\" \"$intergenic_cutoffs_string\" \"$groups_string\" $grand_cds_start $grand_cds_end 2>/dev/null";

    print BLUE, "\n\t$r_cmd\n", RESET;
    print LOG "\n\t$r_cmd\n";
    system($r_cmd);

    my $tmp_path = "images/$ensg_id".".svgz";
    my $tmp_path2 = "images/$ensg_id"."_hm.svgz";

    $figure_content .= "\n<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[leg2h]\" data-openimage=\"../../images/QuestionOff_icon.gif\" data-closedimage=\"../../images/QuestionOn_icon.gif\"><img src=\"../../images/QuestionOn_icon.gif\" border=\"0\" /></a></P>";
    $figure_content .= "\n<div id=\"leg2h\" style=\"width: 1000px;\">";
    $figure_content .= "\n<P CLASS=\"Indented24LR_s16\">The following figures illustrate the expression of the gene '<i>$gene_name</i>' and its component features for each library.  First, a simple gene model is depicted.  In this model, features that are alternatively expressed in at least one comparison are marked red (see below for the actual expression levels of each feature).  Next, the expression values are displayed for only the exons and known (or expressed novel) exon junctions.  This first display consists of a single plot with one colored line per library (as indicated in the legend).  Next, the expression of the gene relative to the distribution of all gene expression values is displayed as a histogram.  In these figures (one for each library) the expression of the current gene is indicated by a dotted red line.  Estimated cutoffs for background expression level corresponding to intergenic and intragenic noise are indicated by dotted black lines.  The bar plots following these histograms display the expression level of all individual features.  The color of the bars correspond to different feature types (enumerated as colored boxes in the legend).  As in the histograms, the estimated cutoffs level for intergenic and intragenic noise are indicated as dotted lines.  For genes with low expression values, these cutoffs converge to a single value.  Features with significant alternative expression values are highlighted yellow in the line plot, and marked with a pink asterix below.  *If you can not see the figure below, click <a href=\"../../../svg_help.htm\">here</a></P><BR>";
    $figure_content .= "\n</div>";
    $figure_content .= "\n<P CLASS=\"Indented24LR\"><OBJECT DATA=\"$tmp_path\" TYPE=\"image/svg+xml\" WIDTH=\"$image_width\" HEIGHT=\"$image_height\" TITLE=\"$image_desc\" CLASS=\"Pic_bordered\" STANDBY=\"$image_desc\"><EMBED SRC=\"$tmp_path\" TYPE=\"image/svg+xml\" WIDTH=\"$image_width\" HEIGHT=\"$image_height\" /></OBJECT></P>";
    $figure_content .= "\n<P CLASS=\"Indented24LR\"><OBJECT DATA=\"$tmp_path2\" TYPE=\"image/svg+xml\" WIDTH=\"$image_width\" HEIGHT=\"$image_height2\" TITLE=\"$image_desc\" CLASS=\"Pic_bordered\" STANDBY=\"$image_desc\"><EMBED SRC=\"$tmp_path2\" TYPE=\"image/svg+xml\" WIDTH=\"$image_width\" HEIGHT=\"$image_height2\" /></OBJECT></P>";

    $figure_content .= "\n<BR></div><BR>";


    #Compress the resulting file
    my $cmd_gz = "gzip -f -S z $images_dir"."$ensg_id".".svg"; 
    my $cmd_gz2 = "gzip -f -S z $images_dir"."$ensg_id"."_hm.svg"; 
    system($cmd_gz);
    system($cmd_gz2);

  }else{
    print BLUE, "\n\tInsufficient data to create figure - i.e. nothing was expressed above 0\n", RESET;
    print LOG "\n\tInsufficient data to create figure - i.e. nothing was expressed above 0\n";
    $figure_content .= "\n<P CLASS=\"Indented24LR_s16\">Insufficient data to create figures</P><BR>";
  }
  return($figure_content);
}


