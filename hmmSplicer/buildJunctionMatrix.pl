#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2010 Malachi Griffith

#This script builds a junction expression matrix and matching annotation matrix and stores them in a new analysis dir
#1.) For all the specified projects get the non-redundant list of annotations.
#  - Make sure an annotation file is available for each library
#  - Store the list of junctions in a hash (merge values across projects where appropriate, e.g. read count)
#2.) Get a list of library IDs and results files for each project considered
#     - Using this list get corresponding junction list files
#     - Build a hash of library IDs with a results file
#3.) Go through each junction list file and store the corresponding junction counts for the library
#    $junction{$jid}{$library}=$read_count
#    $junction{$jid}{TOTAL}+=$read_count
#    $junction{$jid}{LIB_COUNT}++
#    - Make sure the junction is defined in the list of annotated junctions
#    - Apply the filters (if any at this point). If a junction does not pass the filter, do not store it
#4.) Print the output junction expression matrix and matching annotation matrix to file
#    - For each junction, loop through each library in the list and print the counts to file.  For undefined libraries print a count of '0'
#    - At this point apply any filters that could not be applied earlier (e.g. Min library count and total read count)
#    - For each junction that gets printed, print corresponding annotation values

#For all of the steps above.  Allow for sorting of junctions.  This sorting may be applied during import of annotated junctions 
#Possible sorting options:
#Min. total read count
#Min. library count (# of libraries where this junction was observed)
#Min. score
#Min. left size & right size
#Anchored class list (e.g. DA, NDA, D, A, N)
#Max. exons skipped
#Splign validated

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;
use Benchmark;

my $analysis_dir = ''; #Base dir for junction analysis
my $results_dir = '';  #Base dir for output
my $project_list = ''; #Grand list of libaries to be considered - Name of this file will be used to create a results dir
my $matrix_desc = '';  #Description of comparison or group of projects being joined.

my $min_read_count = '';
my $min_lib_count = '';
my $min_score = '';
my $min_left_size = '';
my $min_right_size = '';
my $min_intron_size = '';
my $max_exons_skipped = '';
my $splign_validated = '';
my $library_list = '';

GetOptions ('analysis_dir=s'=>\$analysis_dir, 'results_dir=s'=>\$results_dir, 'project_list=s'=>\$project_list, 'matrix_desc=s'=>\$matrix_desc,
            'min_read_count=i'=>\$min_read_count, 'min_lib_count=i'=>\$min_lib_count, 'min_score=f'=>\$min_score, 'min_left_size=i'=>\$min_left_size, 'min_right_size=i'=>\$min_right_size,
            'min_intron_size=i'=>\$min_intron_size, 'max_exons_skipped=i'=>\$max_exons_skipped, 'splign_validated=s'=>\$splign_validated, 'library_list=s'=>\$library_list);

if ($analysis_dir && $results_dir && $project_list && $matrix_desc){
  #Parameters ok
}else{
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  print GREEN, "\n\nExample: buildJunctionMatrix  --analysis_dir=/projects/alexa2/hmmSplicer/  --project_list='SA_TN_Breast,HumanBodyMap,REMC,Tonsil'  --results_dir=/projects/alexa2/hmmSplicer/Summary/  --matrix_desc=TN-Breast_vs_Normals\n\n", RESET;
  print GREEN, "\n\tTo limit to a subset of library IDs within the projects specified provide a file using: --library_list\n", RESET;
  print GREEN, "\n\tOptional filters:", RESET;
  print GREEN, "\n\tMin. total read count.              --min_read_count (integer)", RESET;
  print GREEN, "\n\tMin. library count.                 --min_lib_count (integer)", RESET;
  print GREEN, "\n\tMin. hmmSplicer score.              --min_score (float)", RESET;
  print GREEN, "\n\tMin. size of left junction half.    --min_left_size (integer)", RESET;
  print GREEN, "\n\tMin. size of right junction half.   --min_right_size (integer)", RESET;
  print GREEN, "\n\tMin. size of intron.                --min_intron_size (integer)", RESET;
  print GREEN, "\n\tMax. exons skipped.                 --max_exons_skipped (integer)", RESET;
  print GREEN, "\n\tSplign validated.                   --splign_validated (0 or 1)", RESET;
  print GREEN, "\n\n", RESET;
  exit();
}
unless($min_read_count){$min_read_count = 0;}
unless($min_lib_count){$min_lib_count = 0;}
unless($min_score){$min_score = 0;}
unless($min_left_size){$min_left_size = 0;}
unless($min_right_size){$min_right_size = 0;}
unless($min_intron_size){$min_intron_size = 0;}
unless($max_exons_skipped){$max_exons_skipped = 1000000000000;}
unless($splign_validated){$splign_validated = 0;}

#Get list of projects
my @projects = split(",", $project_list);
my $project_count = scalar(@projects);

my $working_dir = &setupDir('-dir_name'=>$matrix_desc);

my %library_filter;
if ($library_list){
  unless (-e $library_list){
    print RED, "\n\nUser specified library list file not found: $library_list\n\n", RESET;
    exit();
  }
  open (LIBS, "$library_list") || die "\n\nCould not open library list file: $library_list\n\n";
  while(<LIBS>){
    chomp($_);
    my @line = split("\t", $_);
    $library_filter{$line[0]}=1;
  }
  close(LIBS);
}

#1.) For all the specified projects get the non-redundant list of annotations.
my $junc_anno = &importJunctionAnnotations('-projects'=>\@projects, '-analysis_dir'=>$analysis_dir);


#2.) Get a list of library IDs and results files for each project considered
my $class_header;
my %libraries = %{&getLibs('-base_dir'=>$analysis_dir, '-projects'=>\@projects)};


#3.) Go through each junction list file and store the corresponding junction counts for the library
my $junc_exp = &importJunctionCounts('-junc_anno'=>$junc_anno, '-libs'=>\%libraries);


#4.) For each junction, loop through each library in the list and print the counts to file.  For undefined libraries print a count of '0'
&writeMatrixFiles('-base_dir'=>$analysis_dir, '-junc_anno'=>$junc_anno, '-junc_exp'=>$junc_exp, '-results_dir'=>$results_dir, '-matrix_desc'=>$matrix_desc, '-libs'=>\%libraries);

#Print out the parameters used to build the matrix and store as a reference
unless($library_list){
  $library_list = "not specified";
}
my $param_file = "$working_dir"."MatrixParameters.txt";
open (PARAM, ">$param_file") || die "\n\nCould not open parameters file: $param_file\n\n";
print PARAM "analysis_dir= $analysis_dir\nresults_dir = results_dir\nproject_list = $project_list\nmatrix_desc = $matrix_desc\nmin_read_count = $min_read_count\nmin_lib_count = $min_lib_count\nmin_score = $min_score\nmin_left_size = $min_left_size\nmin_right_size = $min_right_size\nmin_intron_size = $min_intron_size\nmax_exons_skipped = $max_exons_skipped\nsplign_validated = $splign_validated\nlibrary_list = $library_list\n";
close(PARAM);

print "\n\n";

exit();


######################################################################################################################################
#Check input directory paths and setup the results dir                                                                               #
######################################################################################################################################
sub setupDir{
  my %args = @_;
  my $name = $args{'-dir_name'};

  #Basic checks of input
  unless ($analysis_dir =~ /.*\/$/){
    $analysis_dir .= "/";
  }
  unless (-e $analysis_dir && -d $analysis_dir){
    print RED, "\n\nCould not find analysis dir: $analysis_dir\n\n", RESET;
    exit();
  }
  unless ($results_dir =~ /.*\/$/){
    $results_dir .= "/";
  }
  unless (-e $results_dir && -d $results_dir){
    print RED, "\n\nCould not find results dir: $results_dir\n\n", RESET;
    exit();
  }

  my $working_dir = '';
  if ($project_count == 1){
    $working_dir = $results_dir;
  }else{
    $working_dir = "$results_dir"."$name/";
    unless (-e $working_dir){
      mkdir($working_dir);
    } 
  }

  print BLUE, "\n\nResulting matrix and support files will be written to:\n\t$working_dir", RESET;
  return($working_dir);
}



######################################################################################################################################
#For all the specified projects get the non-redundant list of annotations.                                                       #
######################################################################################################################################
sub importJunctionAnnotations{
  my %args = @_;
  my @projects = @{$args{'-projects'}};
  my $analysis_dir = $args{'-analysis_dir'};

  my %ja;

  print BLUE, "\n\nGetting merged junction annotations for projects: @projects", RESET;

  foreach my $project (@projects){
    print BLUE, "\n\tImporting $project ...", RESET;
    my $anno_file = "$analysis_dir"."$project/results/$project".".junction.list.anno.txt";
    unless (-e $anno_file){
      print RED, "\n\nCould not find project file for project: $project\n$anno_file\n\n", RESET;
      exit();
    }

    #Go through each annotation file and store values.  If multiple projects are being considered, merge the results
    #Read_count - cumulative read count
    #Score - highest observed
    #Left_Size - highest observed
    #Right_Size - highest observed

    open (ANNO, "$anno_file") || die "\n\nCould not open annotation file: $anno_file\n\n";
    my %columns;
    my $header = 1;
    while(<ANNO>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header){
        my $pos = 0;
        foreach my $name (@line){
          $columns{$name}{pos}=$pos;
          $pos++;
        }
        $header = 0;
        next();
      }
      my $jid = $line[$columns{'JID'}{pos}];

      if ($ja{$jid}){
        $ja{$jid}{read_count} += $line[$columns{'Read_Count'}{pos}];
        if ($line[$columns{'Score'}{pos}] > $ja{$jid}{score}){
          $ja{$jid}{score} = $line[$columns{'Score'}{pos}];
        }
        if ($line[$columns{'Left_Size'}{pos}] > $ja{$jid}{left_size}){
          $ja{$jid}{left_size} = $line[$columns{'Left_Size'}{pos}];
        }
        if ($line[$columns{'Right_Size'}{pos}] > $ja{$jid}{right_size}){
          $ja{$jid}{right_size} = $line[$columns{'Right_Size'}{pos}];
        }
      }else{
        $ja{$jid}{read_count} = $line[$columns{'Read_Count'}{pos}];
        $ja{$jid}{score} = $line[$columns{'Score'}{pos}];
        $ja{$jid}{left_size} = $line[$columns{'Left_Size'}{pos}];
        $ja{$jid}{right_size} = $line[$columns{'Right_Size'}{pos}];
        $ja{$jid}{intron_size} = $line[$columns{'Intron_Size'}{pos}];
        $ja{$jid}{splice_site} = $line[$columns{'Splice_Site'}{pos}];
        $ja{$jid}{anchored} = $line[$columns{'Anchored'}{pos}];
        $ja{$jid}{exons_skipped} = $line[$columns{'Exons_Skipped'}{pos}];
        $ja{$jid}{donors_skipped} = $line[$columns{'Donors_Skipped'}{pos}];        
        $ja{$jid}{acceptors_skipped} = $line[$columns{'Acceptors_Skipped'}{pos}];
        $ja{$jid}{splign_validated} = $line[$columns{'Splign_Validated'}{pos}];
        $ja{$jid}{splign_align_count} = $line[$columns{'Splign_Align_Count'}{pos}];
        $ja{$jid}{gene_name} = $line[$columns{'Gene_Name'}{pos}];
      }
    }
    close(ANNO);
  }
  return(\%ja);
}


######################################################################################################################################
#Get library IDs and results files associated with the projects                                                                      #
######################################################################################################################################
sub getLibs{
  my %args = @_;
  my $base_dir = $args{'-base_dir'};
  my @projects = @{$args{'-projects'}};

  print BLUE, "\n\nGetting library IDs and results files for projects: @projects", RESET;

  my %libs1; #Libs to consider
  my %libs2; #Libs where data was actually found
  foreach my $project (@projects){
    print BLUE, "\n\t$project", RESET;
    
    #Get the library list
    my $lib_list = "$base_dir"."$project/jobs/LibraryList.txt";
    print BLUE, "\n\t\t$lib_list", RESET;
     unless (-e $lib_list){
      print RED, "\n\nCould not find library list for $project in: $lib_list\n\n", RESET;
      exit();
    }
    open (LIB, "$lib_list") || die "\n\nCould not open library list: $lib_list\n\n";
    while(<LIB>){
      chomp($_);
      my @line = split("\t", $_);
      $libs1{$line[0]}{project}=$project;
    }
    close(LIB);

    #Get the library class entry if available
    my $lib_classes = "$base_dir"."$project/jobs/LibraryClasses.txt";
    print BLUE, "\n\t\t$lib_classes", RESET;
    unless (-e $lib_classes){
      print RED, "\n\nCould not find library classes for $project in: $lib_classes\n\n", RESET;
      exit();
    }
    open (CLASSES, "$lib_classes") || die "\n\nCould not open library classes: $lib_classes\n\n";
    my $header = 1;
    while(<CLASSES>){
      chomp($_);
      unless ($_ =~ /\w+/){
        next();
      }
      my @line = split("\t", $_);
      if ($header){
        $class_header = $_;
        $header = 0;
        next();
      }
      $libs1{$line[0]}{classes}=$_;
    }
    close(CLASSES);

  }
  my $lib1_count = keys %libs1;
  print BLUE, "\n\tFound $lib1_count libraries to consider", RESET;

  #Now determine the libraries that will actually be used
  foreach my $lib (sort keys %libs1){
    my $result_file = "$base_dir"."$libs1{$lib}{project}/results/$lib/$lib".".junction.list.txt";
    if (-e $result_file){
      if ($library_list){
        unless ($library_filter{$lib}){
          next();
        }
      }
      $libs2{$lib}{project} = $libs1{$lib}{project};
      $libs2{$lib}{classes} = $libs1{$lib}{classes};
      $libs2{$lib}{file} = $result_file;
    }
  }
  my $lib2_count = keys %libs2;
  print BLUE, "\n\tFound $lib2_count libraries with a results file", RESET;
  
  #Set the order for each library - order by project, and then within each project by library name
  my $order = 0;
  foreach my $project (@projects){
    foreach my $lib (sort keys %libs2){
      unless ($libs2{$lib}{project} eq $project){
        next();
      }
      $order++;
      $libs2{$lib}{order} = $order;
    }
  }
  return(\%libs2);
}


######################################################################################################################################
#Go through each junction list file and store the corresponding junction counts for the library                                      #
######################################################################################################################################
sub importJunctionCounts{
  my %args = @_;
  my $libs_ref = $args{'-libs'};
  my $junc_anno_ref = $args{'-junc_anno'};

  print BLUE, "\n\nImporting junction counts for each individual library", RESET;
  my %je;
  foreach my $lib (sort {$libs_ref->{$a}->{order} <=> $libs_ref->{$b}->{order}} keys %{$libs_ref}){
    my $lib_project = $libs_ref->{$lib}->{project};
    my $lib_file = $libs_ref->{$lib}->{file};
    print BLUE, "\n\tProcessing: $lib_project -> $lib", RESET;

    open (JUNC, "$lib_file") || die "\n\nCould not open junction read count file: $lib_file\n\n";
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
      my $rc = $line[$columns{'Read_Count'}{pos}];
        
      #Make sure the junction is defined in the list of annotated junctions
      unless($junc_anno->{$jid}){
        print RED, "\n\nEncountered a junction in a lib junction count file that was not in the annotation files - update annotations!!", RESET;
        exit();
      }

      #Apply filters (if any) at this point.  If the junction does not pass the filter, do not store it!
      unless($junc_anno_ref->{$jid}->{score} >= $min_score){next();}
      unless($junc_anno_ref->{$jid}->{left_size} >= $min_left_size){next();}
      unless($junc_anno_ref->{$jid}->{right_size} >= $min_right_size){next();}
      unless($junc_anno_ref->{$jid}->{intron_size} >= $min_intron_size){next();}
      if ($junc_anno_ref->{$jid}->{exons_skipped} =~ /\d+/){
        unless($junc_anno_ref->{$jid}->{exons_skipped} <= $max_exons_skipped){next();}
      }
      unless($junc_anno_ref->{$jid}->{splign_validated} >= $splign_validated){next();}


      #Store the junction, keyed on library
      $je{$jid}{$lib} = $rc;

      #Store the grand read count for this junction 
      $je{$jid}{TOTAL} += $rc;

      #Store the total libraries where this junction was observed
      $je{$jid}{LIB_COUNT}++;
    }
    close(JUNC);
  }
  return(\%je);
}


###########################################################################################################################################
#4.) For each junction, loop through each library in the list and print the counts to file.  For undefined libraries print a count of '0' #
###########################################################################################################################################
sub writeMatrixFiles{
  my %args = @_;
  my $junc_anno_ref = $args{'-junc_anno'};
  my $junc_exp_ref = $args{'-junc_exp'}; 
  my $results_dir = $args{'-results_dir'};
  my $matrix_desc = $args{'-matrix_desc'};
  my $libs_ref = $args{'-libs'};
  my $base_dir = $args{'-base_dir'};

  my $j_matrix_file = "$working_dir"."JunctionExpressionMatrix.tsv";
  my $new_lib_summary_file = "$working_dir"."LibrarySummary.tsv";
  my $new_lib_classes_file = "$working_dir"."LibraryClasses.tsv";
  print BLUE, "\n\nWriting matrix and library summary files:", RESET;
  print BLUE, "\n\t$j_matrix_file\n\t$new_lib_summary_file\n\t$new_lib_classes_file", RESET;

  open (J_OUT, ">$j_matrix_file") || die "\n\nCould not open output junction expression file: $j_matrix_file\n\n";

  my $header_string = "JID\tRead_Count\tScore\tLeft_Size\tRight_Size\tIntron_Size\tSplice_Site\tAnchored\tExons_Skipped\tDonors_Skipped\tAcceptors_Skipped\tGene_Name\t";
  foreach my $lib (sort {$libs_ref->{$a}->{order} <=> $libs_ref->{$b}->{order}} keys %{$libs_ref}){
    $header_string .= "$lib\t";
  }
  $header_string .= "Total\tLibrary_Count\n";
  print J_OUT "$header_string";

  #How many libraries are expected
  my $lib_count = keys %{$libs_ref};
  my $cols_expected = $lib_count+12+2;

  my $c = 0;
  print BLUE, "\n\nProcessing", RESET;
  foreach my $jid (sort keys %{$junc_exp_ref}){
    my $j_string = "$jid\t$junc_anno_ref->{$jid}->{read_count}\t$junc_anno_ref->{$jid}->{score}\t$junc_anno_ref->{$jid}->{left_size}\t$junc_anno_ref->{$jid}->{right_size}\t$junc_anno_ref->{$jid}->{intron_size}\t$junc_anno_ref->{$jid}->{splice_site}\t$junc_anno_ref->{$jid}->{anchored}\t$junc_anno_ref->{$jid}->{exons_skipped}\t$junc_anno_ref->{$jid}->{donors_skipped}\t$junc_anno_ref->{$jid}->{acceptors_skipped}\t$junc_anno_ref->{$jid}->{gene_name}\t";

    #Apply additional filters here (min.read.count & min.lib.count)
    unless($junc_exp_ref->{$jid}->{TOTAL} >= $min_read_count){next();}
    unless($junc_exp_ref->{$jid}->{LIB_COUNT} >= $min_lib_count){next();}

    $c++;
    if ($c == 10000){
      $| = 1; print BLUE, ".", RESET; $| = 0;
      $c = 0;
    }

    #For each junction, loop through each library (in order) and print the counts to file.  For undefined libraries print a count of '0'
    foreach my $lib (sort {$libs_ref->{$a}->{order} <=> $libs_ref->{$b}->{order}} keys %{$libs_ref}){
      if ($junc_exp_ref->{$jid}->{$lib}){
        $j_string .= "$junc_exp_ref->{$jid}->{$lib}\t";
      }else{
        $j_string .= "0\t";
      }
    }
    $j_string .= "$junc_exp_ref->{$jid}->{TOTAL}\t$junc_exp_ref->{$jid}->{LIB_COUNT}\n";

    my @cols = split("\t", $j_string);
    my $cols_found = scalar(@cols);
  
    unless ($cols_found == $cols_expected){
      print RED, "\n\nFound $cols_found columns to be written but $cols_expected were expected - aborting\n\n", RESET;
      exit();
    }

    print J_OUT "$j_string";

  }

  close(J_OUT);

  
  #Now write a summary file that includes only those libraries written to the matrix file
  my $header_line;
  my %lib_summary;
  my $order = 0;
  foreach my $project (@projects){
    my $lib_summary = "$base_dir"."$project/results/$project".".junction.library.summary.txt";
    my $header = 1;
    open (LIB_SUMMARY, "$lib_summary") || die "\n\nCould not find library summary file for project ($project): $lib_summary\n\n";
    while(<LIB_SUMMARY>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        $header_line = $_;
        $header = 0;
        next();
      }
      $order++;
      my $lib = $line[0];
      if ($libs_ref->{$lib}){
        $lib_summary{$lib}{line} = $_;
        $lib_summary{$lib}{order} = $order;
        $lib_summary{$lib}{project} = $project;
      }
    }
    close(LIB_SUMMARY);
  }
  open (NEW_LIB_SUMMARY, ">$new_lib_summary_file") || die "\n\nCould not open new library file: $new_lib_summary_file\n\n";
  print NEW_LIB_SUMMARY "$header_line"."Project\n";
  foreach my $lib (sort {$lib_summary{$a}->{order} <=> $lib_summary{$b}->{order}} keys %lib_summary){
    print NEW_LIB_SUMMARY "$lib_summary{$lib}{line}"."$lib_summary{$lib}{project}\n";
  }
  close(NEW_LIB_SUMMARY);

  #Now write a library classes file that includes only those libraries written to the matrix file
  open (NEW_LIB_CLASSES, ">$new_lib_classes_file") || die "\n\nCould not open new library classes file: $new_lib_classes_file\n\n";
  print NEW_LIB_CLASSES "$class_header\n";
  foreach my $lib (sort {$libs_ref->{$a}->{order} <=> $libs_ref->{$b}->{order}} keys %{$libs_ref}){
    my $class_string = $libs_ref->{$lib}->{classes};
    print NEW_LIB_CLASSES "$class_string\n";
  }
  close(NEW_LIB_CLASSES);

  return();
}


