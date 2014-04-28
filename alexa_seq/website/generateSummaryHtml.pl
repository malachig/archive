#!/usr/bin/perl -w
#Written by Malachi Griffith and Obi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to generate simple webpages that summarize 'Top' lists of differentially expression genes, transcripts, etc.

#DE = Differentially expressed
#DS = Differentially spliced

#Main page:
#Contains a brief introduction and links to the following summary pages
# - Included in this introduction should be basic stats about the databases of sequences used for this analysis (or maybe a seperate page for this...)
#Top Expressed Genes and Transcripts
#Top DE Genes and Transcripts
#Top DE Exons, Junctions, Introns and Intergenic Regions
#Top DS Transcripts, Exons, Junctions, Introns, etc.  (sorted by SI VALUE)
#Top DS Transcripts, Exons, Junctions, Introns, etc.  (sorted by RECIPROCITY VALUE)
#Top DS Transcripts, Exons, Junctions, Introns, etc.  (sorted by PERCENT CONTRIBUTION VALUE)

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

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $annotation_dir = '';
my $junction_seq_size = '';
my $boundary_seq_size = '';
my $paths_file = '';  #File containing root paths to basic types of data (expression, DE and SI)
my $ensembl_version = '';
my $web_root = '';
my $partition_file = '';
my $ucsc_build = '';
my $track_url = '';
my $search_page_url = '';
my $project_title = '';
my $project_name = '';
my $alexa_home_path = '';
my $alexa_seq_path = '';
my $google_analytics_id = '';
my $project_config_file = '';

GetOptions('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'annotation_dir=s'=>\$annotation_dir, 
           'junction_seq_size=i'=>\$junction_seq_size, 'boundary_seq_size=i'=>\$boundary_seq_size,
           'paths_file=s'=>\$paths_file, 'ensembl_version=s'=>\$ensembl_version, 'web_root=s'=>\$web_root, 'partition_file=s'=>\$partition_file,
           'ucsc_build=s'=>\$ucsc_build, 'track_url=s'=>\$track_url, 'search_page_url=s'=>\$search_page_url, 'project_title=s'=>\$project_title, 'project_config_file=s'=>\$project_config_file,
           'project_name=s'=>\$project_name, 'alexa_home_path=s'=>\$alexa_home_path, 'alexa_seq_path=s'=>\$alexa_seq_path, 'google_analytics_id=s'=>\$google_analytics_id);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script generates HTML pages to summarize expression, DE and SI results", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the base path to annotation files using:  --annotation_dir", RESET;
print GREEN, "\n\tSpecify the junction database sequence length using:  --junction_seq_size", RESET;
print GREEN, "\n\tSpecify the boundary database sequence length using:  --boundary_seq_size", RESET;
print GREEN, "\n\tSpecify a file containing root paths to data using: --paths_file", RESET;
print GREEN, "\n\tSpecify the version of EnsEMBL to get data for using: --ensembl_version", RESET;
print GREEN, "\n\tSpecify the complete path to the target web directory using: --web_root", RESET;
print GREEN, "\n\tSpecify a file containing the coordinates of genome partitions used for the analysis using: --partition_file", RESET;
print GREEN, "\n\tSpecify the UCSC build to use for custom tracks using: --ucsc_build (e.g. hg18)", RESET;
print GREEN, "\n\tSpecify the html path to UCSC track files using: --track_url", RESET;
print GREEN, "\n\tSpecify the URL to your Xapian-Omega search page using:  --search_page_url", RESET;
print GREEN, "\n\tSpecify a project title for the summary data table:  --project_title", RESET;
print GREEN, "\n\tSpecify a project name for the summary data table:  --project_name", RESET;
print GREEN, "\n\tSpecify the project config file for this project using: --project_config_file", RESET;
print GREEN, "\n\tSpecify the URL to your ALEXA home page using: --alexa_home_path", RESET;
print GREEN, "\n\tSpecify the URL to your ALEXA-Seq results page using: --alexa_seq_path", RESET;
print GREEN, "\n\tSpecify your Google Analytics ID using: --google_analytics_id", RESET;

print GREEN, "\n\nUsage: generateSummaryHtml.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --annotation_dir=/projects/malachig/sequence_databases/hs_53_36o/  --junction_seq_size=62  --boundary_seq_size=62  --paths_file=/projects/malachig/solexa/batch_jobs/library_paths.txt  --ensembl_version=53  --ucsc_build=hg18  --web_root=/projects/malachig/solexa/temp/  --partition_file=/projects/malachig/solexa/batch_jobs/EnsEMBL_53_Regions_50genes.txt  --track_url=http://www.alexaplatform.org/alexa_seq/5FU/ucsc/  --search_page_url=http://www.bcgsc.ca/xapian-search/omega  --project_title='Alternative expression analysis of 5-FU sensitive and resistant colorectal cancer cell lines'  --project_name=5FU  --project_config_file=/projects/malachig/alexa_seq/config_files/project_config_files/5FU/ALEXA_Seq_5FU.conf  --alexa_home_path=http://www.alexaplatform.org/index.htm  --alexa_seq_path=http://www.alexaplatform.org/alexa_seq/results.htm  --google_analytics_id=UA-xxxxxx-x\n\n", RESET;

#Check user supplied options
unless ($database && $server && $user && $password && $annotation_dir && $junction_seq_size && $boundary_seq_size && $paths_file && $ensembl_version && $web_root && $partition_file && $ucsc_build && $track_url && $search_page_url && $project_title && $project_config_file && $alexa_home_path && $alexa_seq_path && $google_analytics_id){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}
$web_root = &checkDir('-dir'=>$web_root, '-clear'=>"no");
my $web_images_dir = "$web_root"."images";
$web_images_dir = &checkDir('-dir'=>$web_images_dir, '-clear'=>"no");
$annotation_dir = &checkDir('-dir'=>$annotation_dir, '-clear'=>"no");

#1.) Get root paths to the data and import genome partition values
my $data_paths_ref = &getDataPaths('-file'=>$paths_file);
my $expression_paths_ref = $data_paths_ref->{'Expression'}->{paths};
my $de_paths_ref = $data_paths_ref->{'DifferentialExpression'}->{paths};
my $partitions_ref = &getPartitions('-file'=>$partition_file);
my $si_paths_ref = $data_paths_ref->{'DifferentialSplicing'}->{paths};

#Import gene level expression data. Get a list of genes that expressed above background in at least one condition
#NOTE:  This limits the total counts of exons, junctions, etc. expressed to only those the correspond to an expressed gene (in one or more libraries)
my %tmp;
$tmp{'Gene'}=1;
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
my $a_ref = &getAnnotations('-annotation_dir'=>$annotation_dir, '-dbh'=>$alexa_dbh, '-database'=>$database, '-types_list'=>\%tmp, '-light'=>1, '-junction_seq_size'=>$junction_seq_size, '-boundary_seq_size'=>$boundary_seq_size);
$alexa_dbh->disconnect();
my $expression_data_ref = &getExpressionData('-paths'=>$expression_paths_ref, '-a_ref'=>$a_ref, '-ensembl_version'=>$ensembl_version, '-import_option'=>"expressed", '-types_list'=>\%tmp);
my $libraries_ref = $expression_data_ref->{'Gene'}->{libraries};
my %gene_id_list;
foreach my $library (sort {$libraries_ref->{$a}->{line_order} <=> $libraries_ref->{$b}->{line_order}} keys %{$libraries_ref}){
  my $data_ref = $libraries_ref->{$library}->{data};
  while (my ($id) = each %{$data_ref}){
    if ($data_ref->{$id}->{expressed}){
      $gene_id_list{$id}=1;
    }
  }
}

#Define the data types that will be processed in this script
my @types_list = qw (Gene Transcript ExonRegion Junction KnownJunction NovelJunction Boundary KnownBoundary NovelBoundary Intron ActiveIntronRegion SilentIntronRegion Intergenic ActiveIntergenicRegion SilentIntergenicRegion); 
#my @types_list = qw (Gene ActiveIntergenicRegion); 

my $mem_message;
my %expression_data_info;
my $expression_data_ref_info = \%expression_data_info;
my %de_data_info;
my $de_data_ref_info = \%de_data_info;
my %si_data_info;
my $si_data_ref_info = \%si_data_info;
my %library_list_info;

foreach my $type (@types_list){

  my @types_list;
  push(@types_list, $type);

  my %types_list;
  foreach my $type (@types_list){
    $types_list{$type}=1;
  }

  #0.) Get annotations for genes, transcripts, and all seq types.
  my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
  %{$a_ref} = ();
  $a_ref = &getAnnotations('-annotation_dir'=>$annotation_dir, '-dbh'=>$alexa_dbh, '-database'=>$database, '-types_list'=>\%types_list, '-light'=>1, '-gene_id_list'=>\%gene_id_list, '-junction_seq_size'=>$junction_seq_size, '-boundary_seq_size'=>$boundary_seq_size);
  $alexa_dbh->disconnect();

  #Assign every annotated object to a genome partition (e.g. chr1_1)
  my $assigned_partitions = &assignPartitions('-a_ref'=>$a_ref, '-partitions_ref'=>$partitions_ref);
  print BLUE, "\n\t\tAssigned $assigned_partitions records to a single genome partition", RESET;

  #2-A.) Import expression data - only expressed features should be imported for this summary
  %{$expression_data_ref} = ();
  $expression_data_ref = &getExpressionData('-paths'=>$expression_paths_ref, '-a_ref'=>$a_ref, '-ensembl_version'=>$ensembl_version, '-import_option'=>"expressed", '-types_list'=>\%types_list);

  my $type_count = keys %{$expression_data_ref};
  print BLUE, "\n\nImported expression data for $type_count data types", RESET;

  #Hard code display order etc. for each data type
  &annotateDataTypes('-data_ref'=>$expression_data_ref);

  #Get a list of libraries, comparisons, etc. and their order according to the input file

  #2-B.) Now go through each data type and create a website to display the top N% (by expression rank) for it:
  print MAGENTA, "\n\nCreating EXPRESSION DATA website content for each data type (each page will summarize all libraries for that type)", RESET;
  my %library_list;
  foreach my $data_type (sort {$expression_data_ref->{$a}->{order} <=> $expression_data_ref->{$b}->{order}} keys %{$expression_data_ref}){

    print BLUE, "\n\nProcessing data type: $data_type", RESET;

    $libraries_ref = $expression_data_ref->{$data_type}->{libraries};

    my @library_names;
    my $lib_count = 0;
    foreach my $library (sort {$libraries_ref->{$a}->{line_order} <=> $libraries_ref->{$b}->{line_order}} keys %{$libraries_ref}){
      $lib_count++;
      #Make sure an images dir exists for this library
      my $images_path = "$web_root"."images/"."$library/";
      unless(-e $images_path && -d $images_path){
        my $cmd = "mkdir $images_path";
        system($cmd);
      }
      my $library_name = $expression_paths_ref->{$library}->{name};
      push(@library_names, $library_name);
      print BLUE, "\n\tProcessing library: $library", RESET;

      $library_list{$library}{name} = $library_name;
      $library_list{$library}{order} = $libraries_ref->{$library}->{line_order};

      $library_list_info{$library}{name} = $library_name;
      $library_list_info{$library}{order} = $libraries_ref->{$library}->{line_order};

      my $data_ref = $libraries_ref->{$library}->{data};

      #Rank the data
      my $assigned_ranks = &rankData('-data_ref'=>$data_ref, '-rank_variable'=>'norm_expression', '-descending'=>1, '-absolute'=>0);
      print BLUE, "\n\t\tAssigned ranks for $assigned_ranks records", RESET;
    }

    #Check for matrix file
    my $matrix_file = "data/matrix/Matrix_"."$data_type"."Expression_v$ensembl_version".".txt.gz";
    my $matrix_file_name = "$data_type"."Expression_v$ensembl_version".".txt.gz";
    my $matrix_file_path = "$web_root"."$matrix_file";
    if (-e $matrix_file_path){
      print BLUE, "\n\t\tFound matrix file: $matrix_file", RESET;
    }else{
      print YELLOW, "\n\t\tMatrix file missing: $matrix_file_path", RESET;
      $matrix_file="NA";
    }

    #Generate a content object for the data containing the top N% of the expressed features
    my $top_n_percent = 1;
    my $content_ref = &generateExpressionContent('-data_type'=>$data_type, '-data_ref'=>$libraries_ref, '-a_ref'=>$a_ref, '-top_n_percent'=>$top_n_percent, '-paths'=>$expression_paths_ref, 
                                                 '-matrix_file'=>$matrix_file, '-matrix_file_name'=>$matrix_file_name, '-web_images_dir'=>$web_images_dir, '-ucsc_build'=>$ucsc_build, '-track_url'=>$track_url);


    #Write a webpage for the data
    my $library_names = join (", ", @library_names);
    my $title = "Distributions of all expression values and Top $top_n_percent% of expressed $data_type features for libraries: $library_names";
    my $meta_description = "Provides a list of the top $top_n_percent% of expressed features of the type '$data_type' for the expression libraries '$library_names'";
    my $meta_keywords = "Expression, $data_type, $library_names";
    my $web_name = "Expression_"."$data_type".".htm";
    my $web_path = "$web_root"."$web_name";
    &writePage('-path'=>$web_path, '-title'=>$title, '-content'=>$content_ref, '-css_path'=>"../ALEXA2.css", '-alexa_home_path'=>"$alexa_home_path", '-alexa_seq_home_path'=>"$alexa_seq_path", '-summary_path'=>"Summary.htm", '-genes_path'=>"genes/index.html", '-search_path'=>$search_page_url, '-meta_description'=>$meta_description, '-meta_keywords'=>$meta_keywords, '-google_analytics'=>1, '-google_analytics_id'=>$google_analytics_id, '-collapse_div_script'=>'../animatedcollapse.js', '-jquery_min_script'=>'../jquery.min.js', '-div_count'=>$lib_count+2);

    $expression_data_ref->{$data_type}->{web_page_name} = $web_name;
    $expression_data_ref->{$data_type}->{web_page_path} = $web_path;

    #Store basic info from $expression_data_ref in a global object
    $expression_data_ref_info->{$data_type}->{web_page_name} = $web_name;
    $expression_data_ref_info->{$data_type}->{web_page_path} = $web_path;
    $expression_data_ref_info->{$data_type}->{level} = $expression_data_ref->{$data_type}->{level};
    $expression_data_ref_info->{$data_type}->{record_count} = $expression_data_ref->{$data_type}->{record_count};
  }

  #3-A.) Import differential expression data
  my $de_data_ref = &getDifferentialExpressionData('-paths'=>$de_paths_ref, '-a_ref'=>$a_ref, '-ensembl_version'=>$ensembl_version, '-import_option'=>"filtered", '-types_list'=>\%types_list);


  #3-B.) Now go through each data type and create a website to display the signficant DE events (sorted by differential expression rank) for it:
  print MAGENTA, "\n\nCreating DIFFERENTIAL EXPRESSION DATA website content for each data type (each page will summarize all libraries for that type)", RESET;
  foreach my $data_type (sort {$expression_data_ref->{$a}->{order} <=> $expression_data_ref->{$b}->{order}} keys %{$expression_data_ref}){
    print BLUE, "\n\nProcessing data type: $data_type", RESET;

   #Only process differential expression values if the data type was already defined for expression values and is defined in the DE object
    unless ($de_data_ref->{$data_type}){
      next();
    }
    my $comparisons_ref = $de_data_ref->{$data_type}->{comparisons};
    my $libraries_ref = $expression_data_ref->{$data_type}->{libraries};

    my @comparison_names;
    my $comp_count = 0;
    foreach my $comparison (sort {$comparisons_ref->{$a}->{line_order} <=> $comparisons_ref->{$b}->{line_order}} keys %{$comparisons_ref}){
      $comp_count++;

      #Make sure an images dir exists for this library
      my $images_path = "$web_root"."images/"."$comparison/";
      unless(-e $images_path && -d $images_path){
        my $cmd = "mkdir $images_path";
        system($cmd);
      }

      my $comparison_name = $de_paths_ref->{$comparison}->{name};
      push(@comparison_names, $comparison_name);
      print BLUE, "\n\tProcessing comparison: $comparison", RESET;

      my $data_ref = $comparisons_ref->{$comparison}->{data};

      #Rank the data
      my $assigned_ranks = &rankData('-data_ref'=>$data_ref, '-rank_variable'=>'fold_change', '-descending'=>1, '-absolute'=>1);
      print BLUE, "\n\t\tAssigned ranks for $assigned_ranks records", RESET;
    }

    my $content_ref = &generateDifferentialExpressionContent('-data_type'=>$data_type, '-data_ref'=>$comparisons_ref, '-a_ref'=>$a_ref, '-de_paths'=>$de_paths_ref, '-exp_paths'=>$expression_paths_ref,
                                                             '-web_images_dir'=>$web_images_dir, '-ucsc_build'=>$ucsc_build, '-track_url'=>$track_url);

    #Write a webpage for the data
    my $comparison_names = join (", ", @comparison_names);
    my $title = "Distributions and lists of all significant differential expression values of the type: $data_type for comparisons: $comparison_names (p-value < 0.05 after multiple testing correction and Fold-Change > 1.5)";
    my $meta_description = "Provides a list of all significant differentially expressed features of the type '$data_type' for the comparisons '$comparison_names'";
    my $meta_keywords = "Differential expression, $data_type, $comparison_names";
    my $web_name = "DE_"."$data_type".".htm";
    my $web_path = "$web_root"."$web_name";
    &writePage('-path'=>$web_path, '-title'=>$title, '-content'=>$content_ref, '-css_path'=>"../ALEXA2.css", '-alexa_home_path'=>"$alexa_home_path", '-alexa_seq_home_path'=>"$alexa_seq_path", '-summary_path'=>"Summary.htm", '-genes_path'=>"genes/index.html", '-search_path'=>$search_page_url, '-meta_description'=>$meta_description, '-meta_keywords'=>$meta_keywords, '-google_analytics'=>1, '-google_analytics_id'=>$google_analytics_id, '-collapse_div_script'=>'../animatedcollapse.js', '-jquery_min_script'=>'../jquery.min.js', '-div_count'=>(($comp_count*3)+2));

    $de_data_ref->{$data_type}->{web_page_name} = $web_name;
    $de_data_ref->{$data_type}->{web_page_path} = $web_path;
    
    #Store basic info from $de_data_ref in a global object
    $de_data_ref_info->{$data_type}->{web_page_name} = $web_name;
    $de_data_ref_info->{$data_type}->{web_page_path} = $web_path;
    $de_data_ref_info->{$data_type}->{record_count} = $de_data_ref->{$data_type}->{record_count};

  }

  #4-A.) Import differential splicing data
  my $si_data_ref = &getDifferentialSplicingData('-paths'=>$si_paths_ref, '-a_ref'=>$a_ref, '-ensembl_version'=>$ensembl_version, '-import_option'=>"filtered", '-types_list'=>\%types_list);

  #4-B.) Now go through each data type and create a website to display the significant differential splicing events (sorted by SI rank) for it:
  #      - Because of the complexity of the data, SI pages will only summarize a single comparison and data type
  #      - Different versions should be created, one for each sorting option
  print MAGENTA, "\n\nCreating DIFFERENTIAL SPLICING DATA website content for each data type (each page will summarize all libraries for that type)", RESET;
  foreach my $data_type_comp (sort {$si_data_ref->{$a}->{line_order} <=> $si_data_ref->{$b}->{line_order}} keys %{$si_data_ref}){

    print BLUE, "\n\nProcessing data type: $data_type_comp", RESET;
  
    my $data_ref = $si_data_ref->{$data_type_comp}->{data};
    my $comparison = $si_data_ref->{$data_type_comp}->{comparison};
    my $data_type = $si_data_ref->{$data_type_comp}->{data_type};

    #Make sure an images dir exists for this library
    my $images_path = "$web_root"."images/"."$comparison/";
    unless(-e $images_path && -d $images_path){
      my $cmd = "mkdir $images_path";
      system($cmd);
    }

    my $comparison_name = $si_paths_ref->{$comparison}->{name};
    print BLUE, "\n\tProcessing data_type comparison: $data_type_comp ($comparison_name)", RESET;

    #Rank using several different metrics for the same lists of events
    my @rank_terms = qw (si gene_fold_change seq_fold_change reciprocity percent_seq_log2_de);

    foreach my $rank_variable (@rank_terms){
      print BLUE, "\n\t\tRank data by: $rank_variable", RESET;
  
      #Rank the data
      my $assigned_ranks = &rankData('-data_ref'=>$data_ref, '-rank_variable'=>$rank_variable, '-descending'=>1, '-absolute'=>1);
      print BLUE, "\n\t\t\tAssigned ranks for $assigned_ranks records", RESET;

      my $content_ref = &generateDifferentialSplicingContent('-data_type_comp'=>$data_type_comp, '-si_data_ref'=>$si_data_ref, '-a_ref'=>$a_ref, '-si_paths'=>$si_paths_ref, '-exp_paths'=>$expression_paths_ref,
                                                             '-web_images_dir'=>$web_images_dir, '-ucsc_build'=>$ucsc_build, '-track_url'=>$track_url, '-rank_terms'=>\@rank_terms);


      #Write a webpage for the data
      my $title = "Alternative expression values of the type: $data_type for comparison: $comparison_name (abs(SI) value > 2)";
      my $meta_description = "Provides a list of all significant differentially spliced features of the type '$data_type' for the comparison '$comparison_name'";
      my $meta_keywords = "Differential splicing, splicing index, reciprocity, $data_type, $comparison_name";
  
      my $web_name = "SI_"."$data_type_comp"."_"."$rank_variable".".htm";
      my $web_path = "$web_root"."$web_name";
    &writePage('-path'=>$web_path, '-title'=>$title, '-content'=>$content_ref, '-css_path'=>"../ALEXA2.css", '-alexa_home_path'=>"$alexa_home_path", '-alexa_seq_home_path'=>"$alexa_seq_path", '-summary_path'=>"Summary.htm", '-genes_path'=>"genes/index.html",  '-search_path'=>$search_page_url, '-meta_description'=>$meta_description, '-meta_keywords'=>$meta_keywords, '-google_analytics'=>1, '-google_analytics_id'=>$google_analytics_id, '-collapse_div_script'=>'../animatedcollapse.js', '-jquery_min_script'=>'../jquery.min.js', '-div_count'=>5);
    }

    my $web_name = "SI_"."$data_type_comp"."_si.htm";
    my $web_path = "$web_root"."$web_name";

    $si_data_ref->{$data_type_comp}->{web_page_name} = $web_name;
    $si_data_ref->{$data_type_comp}->{web_page_path} = $web_path;

    #Store basic info from $si_data_ref in a global object
    $si_data_ref_info->{$data_type_comp}->{web_page_name} = $web_name;
    $si_data_ref_info->{$data_type_comp}->{web_page_path} = $web_path;
    $si_data_ref_info->{$data_type_comp}->{comparison} = $comparison;
    $si_data_ref_info->{$data_type_comp}->{data_type} = $data_type;
    $si_data_ref_info->{$data_type_comp}->{line_order} = $si_data_ref->{$data_type_comp}->{line_order};
    $si_data_ref_info->{$data_type_comp}->{record_count} = $si_data_ref->{$data_type_comp}->{record_count};
  }
  $mem_message = &memoryUsage();
  print YELLOW, "\n$mem_message", RESET;

}

#5.) Print out a main summary page with links to the other pages
print BLUE, "\n\nPrinting the final SUMMARY page with links to Expression, DE and SI data\n\n", RESET;
#Create a summary page which consists of a table of links to Expression, DE and SI data
#For expression and DE, the link will to be a page showing all libraries or comparisons
#For SI, a link will be required to a seperate page for each sample comparison

#Each row in the table will correspond to a single data type (Gene, Transcript, exon region, junction, etc.)
my $summary_content = '';
$summary_content = "$summary_content"."<TABLE CLASS=\"Data\">";

my $exp_icon = "images/EXP_icon.png";
my $de_icon = "images/DE_icon.png";
my $ae_icon = "images/AE_icon.png";

#Get SI comparisons
my $splicing_heads = '';
my %comp_list;
foreach my $data_type_comp (sort {$si_data_ref_info->{$a}->{line_order} <=> $si_data_ref_info->{$b}->{line_order}} keys %{$si_data_ref_info}){
  my $comparison = $si_data_ref_info->{$data_type_comp}->{comparison};
  my $comparison_name = $si_paths_ref->{$comparison}->{name};
  $si_data_ref_info->{$data_type_comp}->{comparison_name} = $comparison_name;
  $comp_list{$comparison_name}{line_order} = $si_data_ref_info->{$data_type_comp}->{line_order};
  $comp_list{$comparison_name}{comparison_id} = $comparison;
}
foreach my $comparison_name (sort {$comp_list{$a}{line_order} <=> $comp_list{$b}{line_order}} keys %comp_list){
  $splicing_heads = "$splicing_heads"."<TD CLASS = \"Head3\">Alternative Expression<BR>($comparison_name)</TD>";
}

#Create links to library summary pages
my $library_summary_links = ""; 
foreach my $library (sort {$library_list_info{$a}->{order} <=> $library_list_info{$b}->{order}} keys %library_list_info){
  my $library_name = $library_list_info{$library}{name};
  my $library_summary_file = "$library".".htm";
  $library_summary_links .=  "\n<P CLASS=\"Indented24LR_s16\"><a href=\"$library_summary_file\">$library</a> ($library_name)</P>";
}

#Check for results package (all results files, matrix files, etc)
my $results_package = "data/"."$project_name"."_ResultFiles.tar.gz";
my $results_package_name = "$project_name"."_ResultFiles.tar.gz";
my $results_package_path = "$web_root"."$results_package";
my $results_package_link = "<A HREF=\"$results_package\">$results_package_name</A>";
if (-e $results_package_path){
  print BLUE, "\n\t\tFound results package: $results_package", RESET;
}else{
  print YELLOW, "\n\t\tResults package missing: $results_package_path", RESET;
  $results_package_link="N/A";
}

#Generate links to each candidate gene list page (one per comparison)
my $candidate_list_links = '';
foreach my $comparison_name (sort {$comp_list{$a}{line_order} <=> $comp_list{$b}{line_order}} keys %comp_list){
  my $comparison_id = $comp_list{$comparison_name}{comparison_id};
  my $candidate_list_file = "$comparison_id".".htm";
  $candidate_list_links .= "\n<P CLASS=\"Indented24LR_s16\"><a href=\"$candidate_list_file\">$comparison_id</a> ($comparison_name)</P>";
}

#Generate links to each peptide list page (one per comparison)
my $peptide_list_links = '';
foreach my $comparison_name (sort {$comp_list{$a}{line_order} <=> $comp_list{$b}{line_order}} keys %comp_list){
  my $comparison_id = $comp_list{$comparison_name}{comparison_id};
  my $peptide_list_file = "$comparison_id"."_peptides.htm";
  $peptide_list_links .= "\n<P CLASS=\"Indented24LR_s16\"><a href=\"$peptide_list_file\">$comparison_id</a> ($comparison_name)</P>";
}

#Load the project config file and get basic info for this project
my $project_conf_ref = &loadConfigFile('-project_config_file'=>$project_config_file);
my $project_info_content = '';
$project_info_content .= "\n<P CLASS=\"Indented24LR_s16\">Project name: '$project_conf_ref->{PROJECT_NAME}'</P>";
$project_info_content .= "\n<P CLASS=\"Indented24LR_s16\">EnsEMBL version: $project_conf_ref->{ENSEMBL_VERSION}</P>";
$project_info_content .= "\n<P CLASS=\"Indented24LR_s16\">UCSC build: $project_conf_ref->{UCSC_BUILD}</P>";
$project_info_content .= "\n<P CLASS=\"Indented24LR_s16\">ALEXA-Seq database: <A HREF =\"../downloads.htm\">$project_conf_ref->{ALEXA_SEQ_DB}</A></P>";
$project_info_content .= "\n<P CLASS=\"Indented24LR_s16\">Species: $project_conf_ref->{SPECIES_NAME} ($project_conf_ref->{SPECIES_NAME_COMMON})</P>";
$project_info_content .= "\n<P CLASS=\"Indented24LR_s16\">Lanes of data analysed: $project_conf_ref->{LANE_COUNT}</P>";
$project_info_content .= "\n<P CLASS=\"Indented24LR_s16\">Library count: $project_conf_ref->{LIBRARY_COUNT}</P>";
$project_info_content .= "\n<P CLASS=\"Indented24LR_s16\">Comparisons defined: $project_conf_ref->{COMPARISON_COUNT}</P>";

#Header row
$summary_content = <<"EOF";

<P CLASS=\"Indented24LR_s16_bold\">Project information:</P>
$project_info_content
<BR>

<P CLASS=\"Indented24LR_s16_bold\">Library statistics (link to library summary):</P>
$library_summary_links
<BR>

<P CLASS=\"Indented24LR_s16_bold\">Download package of all results files: $results_package_link</P>
<BR>

<P CLASS=\"Indented24LR_s16_bold\">View candidate gene lists for each comparison (DE+AE genes, AE genes only, Gains only, Losses only):</P>
$candidate_list_links
<BR>

<P CLASS=\"Indented24LR_s16_bold\">View candidate peptide lists (peptides corresponding to: junction gains, junction losses, exon+junction gains, exon+junction losses):</P>
$peptide_list_links
<BR>


<P CLASS=\"Indented24LR_s16_bold\">Summary of expressed, differentially expressed and alternatively expressed features:</P>
<P CLASS=\"Indented24LR\"><a href=\"#\" rel=\"toggle[leg1s]\" data-openimage=\"images/QuestionOff_icon.gif\" data-closedimage=\"images/QuestionOn_icon.gif\"><img src=\"images/QuestionOff_icon.gif\" border=\"0\" /></a></P>
<div id=\"leg1s\" style=\"width: 1100px;\">
<P CLASS=\"Indented24LR_s16\">The following table summarizes the expression, differential expression and alternative expression of each feature type (exons, junctions, introns, etc.) for all EnsEMBL loci.  It does this by providing links to lists of: the most highly abundant expressed (EXP) features, all significant differentially expressed (DE) features and all significant alternatively expressed (AE) features. Click on an icon to view these lists for the desired feature type.  The number of features is enumerated beneath each icon.  Note that for EXP tables, only the top 1% of expressed features will be summarized.  Example 1: to see the most highly expressed exons in each library, click the 'EXP' icon in the 'ExonRegion' row of the table.  Example 2: to view a ranked list of significant differentially expressed transcripts, click the 'DE' icon in the 'Transcript row of the table.  Example 3: to view a ranked list of alternatively expressed, novel exon-exon junctions click the 'AE' icon in the 'NovelJunction' row of the table. For a detailed description of the feature types and how the ALEXA-Seq data viewer works, refer to our manuscript.  Note that AE and DE events are summarized individually for each feature type.  Since alternative expression is determined relative to gene level expression, the AE analysis is not applicable (N/A) for the gene and intergenic feature rows of this table.  To view data for a specific gene, use the 'SEARCH' link above.</P><BR>
</div>
<TABLE CLASS=\"Data2\">
  <TR>
    <TD CLASS=\"Head1\">Data Type</TD><TD CLASS = \"Head3\">Expression<BR>(All libraries)</TD><TD CLASS = \"Head3\">Differential Expression<BR>(All comparisons)</TD>$splicing_heads
  </TR>
EOF

foreach my $current_data_type (@types_list){

  my $level = $expression_data_ref_info->{$current_data_type}->{level};
  my $td_class = "Head1";
  if ($level == 2){
    $td_class = "Head2";
  }

  #Print the name for this row
  $summary_content .= "\n  <TR>\n    <TD CLASS=\"$td_class\">$current_data_type</TD>";

  #Expression data
  my $record = "<TD CLASS=\"Data2\">N/A</TD>";
  foreach my $data_type (sort {$a cmp $b} keys %{$expression_data_ref_info}){
    unless($data_type eq $current_data_type){
      next();
    }
    my $page_name = $expression_data_ref_info->{$data_type}->{web_page_name};
    my $record_count = &commify($expression_data_ref_info->{$data_type}->{record_count});
    my $img_desc = "Expression report for $data_type";
    $record = "<TD CLASS=\"Data2\"><A HREF=\"$page_name\"><IMG SRC=\"$exp_icon\" CLASS=\"Pic_unbordered\" ALT=\"$img_desc\" TITLE=\"$img_desc\"></A><BR>($record_count)</TD>";
  }
  $summary_content = "$summary_content"."$record";


  #Differential Expression data
  $record = "<TD CLASS=\"Data2\">N/A</TD>";
  foreach my $data_type (sort {$a cmp $b} keys %{$de_data_ref_info}){
    unless($data_type eq $current_data_type){
      next();
    }
    my $page_name = $de_data_ref_info->{$data_type}->{web_page_name};
    my $record_count = &commify($de_data_ref_info->{$data_type}->{record_count});
    my $img_desc = "Differential Expression report for $data_type";
    $record = "<TD CLASS=\"Data2\"><A HREF=\"$page_name\"><IMG SRC=\"$de_icon\" CLASS=\"Pic_unbordered\" ALT=\"$img_desc\" TITLE=\"$img_desc\"></A><BR>($record_count)</TD>";
  }
  $summary_content = "$summary_content"."$record";


  #Differential Splicing data
  $record = "<TD CLASS=\"Data2\">N/A</TD>";
  my $test = 0;
  foreach my $data_type_comp (sort {$si_data_ref_info->{$a}->{line_order} <=> $si_data_ref_info->{$b}->{line_order}} keys %{$si_data_ref_info}){
     my $data_type = $si_data_ref_info->{$data_type_comp}->{data_type};
     unless($data_type eq $current_data_type){
       next();
     }
     $test = 1;
     my $page_name = $si_data_ref_info->{$data_type_comp}->{web_page_name};
     my $record_count = &commify($si_data_ref_info->{$data_type_comp}->{record_count});
     my $img_desc = "Alternative Expression report for $data_type";
     $record = "<TD CLASS=\"Data2\"><A HREF=\"$page_name\"><IMG SRC=\"$ae_icon\" CLASS=\"Pic_unbordered\" ALT=\"$img_desc\" TITLE=\"$img_desc\"></A><BR>($record_count)</TD>";
     $summary_content = "$summary_content"."$record";
  }
  unless ($test){
    foreach my $comparison_name (sort {$comp_list{$a}{line_order} <=> $comp_list{$b}{line_order}} keys %comp_list){
      $summary_content = "$summary_content"."<TD CLASS=\"Data2\">N/A</TD>";
    }
  }

  #Finish the current row
  $summary_content = "$summary_content"."\n  </TR>";

}
#Finish the table
$summary_content = "$summary_content"."</TABLE><BR>";

#Generate the summary home page
my $meta_description = "Provides a central portal to lists of top expressed, differentially expressed and differentially spliced features for all genes";
my $meta_keywords = "Alternative expression analysis, alternative splicing, gene expression analysis, paired-end sequencing, Illumina, Solexa, Differential expression, Alternative Expression, differential splicing, alternative isoforms, transcript discovery";
my $web_path = "$web_root"."Summary".".htm";
&writePage('-path'=>$web_path, '-title'=>$project_title, '-content'=>\$summary_content, '-css_path'=>"../ALEXA2.css", '-alexa_home_path'=>"$alexa_home_path", '-alexa_seq_home_path'=>"$alexa_seq_path", '-summary_path'=>"Summary.htm", '-genes_path'=>"genes/index.html", '-search_path'=>$search_page_url, '-meta_description'=>$meta_description, '-meta_keywords'=>$meta_keywords, '-google_analytics'=>1, '-google_analytics_id'=>$google_analytics_id, '-collapse_div_script'=>'../animatedcollapse.js', '-jquery_min_script'=>'../jquery.min.js', '-div_count'=>12);

$mem_message = &memoryUsage();
print YELLOW, "\n\n$mem_message\n\n", RESET;

exit();


###################################################################################################################################################################
#Import project config file                                                                                                                                       #
###################################################################################################################################################################
sub loadConfigFile{
  my %args = @_;
  my $file = $args{'-project_config_file'};
  my %project_conf;
  open(CONF, "$file") || die "\n\nCould not open file: $file\n\n";
  my $lane_count = 0;
  my $library_count = 0;
  my $comparison_count = 0;
  my %test;
  while(<CONF>){
    chomp($_);
    #Skip comment lines
    if ($_ =~ /^\#/){
      next();
    }
    #Skip empty lines
    unless ($_ =~ /\w+|\d+/){
      next();
    }
    #Watch out for project specific configuration values
    if ($_ =~ /(.*)\=(.*)/){
      $project_conf{$1}=$2;
    }
    my @line = split(/ +/, $_);
    if ($line[0] =~ /^LANE/){
      $lane_count++;
    }
    if ($line[0] =~ /^LIBRARY/){
      $library_count++;
    }
    if ($line[0] =~ /^COMPARISON/){
      $comparison_count++;
    }
  }
  close(CONF);
  $project_conf{LANE_COUNT} = $lane_count;
  $project_conf{LIBRARY_COUNT} = $library_count;
  $project_conf{COMPARISON_COUNT} = $comparison_count;
  return(\%project_conf);
}

