#!/usr/bin/perl -w
#Written by Malachi Griffith and Obi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#This script takes an ALEXA-seq configuration file and Project specific configuration file as input and uses these to generate all necessary commands for ALEXA-Seq analysis
use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

#Load the ALEXA modules
my $script_dir;
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\.pl/){
    push (@INC, $1);
    $script_dir = $1;
  }
}
use utilities::utility qw(:all);

my $alexa_seq_config_file = '';
my $project_config_file = '';
my $commands_file = '';
my $cluster_commands = '';
my $clean = '';
my $demo = '';

GetOptions ('alexa_seq_config_file=s'=>\$alexa_seq_config_file, 'project_config_file=s'=>\$project_config_file, 'commands_file=s'=>\$commands_file, 'clean=s'=>\$clean, 'cluster_commands=s'=>\$cluster_commands, 'demo=s'=>\$demo);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script takes an ALEXA-Seq configuration file and Project specific configuration file as input and uses these to generate all necessary commands for ALEXA-Seq analysis", RESET;
print GREEN, "\n\tSpecify the complete path to your ALEXA-Seq configuration file using: --alexa_seq_config_file", RESET;
print GREEN, "\n\tSpecify the complete path to your Project specific configuration file using: --project_config_file", RESET;
print GREEN, "\n\tSpecify the file to which all commands will be written using: --commands_file", RESET;
print GREEN, "\n\n\tNOTE:", RESET;
print GREEN, "\n\tIf you have changed the LANE, LIBRARY or COMPARISON entries and wish to regenerate the commands cleanly, use: --clean=1", RESET;
print GREEN, "\n\tTo include GSC specific cluster submission commands (based on 'mqsub' wrapper for SGE queue submission manager) use: --cluster_commands=1", RESET;

print GREEN, "\n\nExample: createAnalysisCommands.pl  --alexa_seq_config_file=/home/malachig/svn/alexa_seq/config_files/ALEXA_Seq_PIPELINE.conf  --project_config_file=/home/malachig/svn/alexa_seq/config_files/project_config_files/LnCAP_AR_KnockIn/ALEXA_Seq_LnCAP_AR_KnockIn.conf  --commands_file=/home/malachig/svn/alexa_seq/config_files/project_config_files/LnCAP_AR_KnockIn/ALEXA_Seq_LnCAP_AR_KnockIn.commands\n\n", RESET;

unless ($alexa_seq_config_file && $project_config_file && $commands_file){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

unless (-e $alexa_seq_config_file && -e $project_config_file){
  print RED, "\n\nCould not verify existence of both config files ... check paths\n\n", RESET;
  exit();
}

unless($clean){
  $clean = 0;
}

#0.) Import config files
open(COMMANDS, ">$commands_file") || die "\n\nCould not open output commands file: $commands_file\n\n";
print MAGENTA, "\n\nImporting config files", RESET;
print MAGENTA, "\nPIPELINE: $alexa_seq_config_file";
print MAGENTA, "\nPROJECT:  $project_config_file";
print COMMANDS "#ALEXA-SEQ ANALYSIS - DETAILED INSTRUCTIONS";
print COMMANDS "\n#Using the following ALEXA-seq and Project specific config files for this analysis:";
print COMMANDS "\n#$alexa_seq_config_file\n#$project_config_file";

#Import ALEXA-Seq config file
my %alexa_seq_conf;
open(CONF, "$alexa_seq_config_file") || die "\n\nCould not open file: $alexa_seq_config_file\n\n";
while(<CONF>){
  chomp($_);
  #Skip comment lines
  if ($_ =~ /^\#/){
    next();
  }
  if ($_ =~ /(.*)\=(.*)/){
    $alexa_seq_conf{$1}=$2;
  }
}
close(CONF);

#Import project config file
my %project_conf;
my %tmp1; my %tmp2; my %tmp3; my %tmp4; my %tmp5;
$project_conf{LANE} = \%tmp1;
$project_conf{LIBRARY} = \%tmp2;
$project_conf{COMPARISON} = \%tmp3;
$project_conf{GROUP} = \%tmp4;
$project_conf{GROUP_COMPARISON} = \%tmp5;

open(CONF, "$project_config_file") || die "\n\nCould not open file: $project_config_file\n\n";
my $lane_count = 0;
my $library_count = 0;
my $comparison_count = 0;
my $group_count = 0;
my $group_comparison_count = 0;

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

  #Split on any whitespace characters
  my @line = split(/\s+/, $_);
  if ($line[0] =~ /^LANE/){
    $lane_count++;
    my $lanes_ref = $project_conf{LANE};
    $lanes_ref->{$lane_count}->{library_id} = $line[1];
    $lanes_ref->{$lane_count}->{flowcell_name} = $line[2];
    $lanes_ref->{$lane_count}->{lane_number} = $line[3];
    $lanes_ref->{$lane_count}->{data_path} = $line[4];
    $lanes_ref->{$lane_count}->{seq_file_type} = $line[5];
    $lanes_ref->{$lane_count}->{read_length} = $line[6];
    $lanes_ref->{$lane_count}->{read_trim} = $line[7];
    $lanes_ref->{$lane_count}->{max_n} = $line[8];
    $lanes_ref->{$lane_count}->{min_phred} = $line[9];
    $lanes_ref->{$lane_count}->{qual_type} = $line[10];

    my $lib_flowcell_lane = "$line[1]"."_"."$line[2]"."_"."$line[3]";
    $test{$lib_flowcell_lane}=1;
  }
  if ($line[0] =~ /^LIBRARY/){
    $library_count++;
    my $lib_ref = $project_conf{LIBRARY};
    my $lib_id = $line[1];
    $lib_ref->{$lib_id}->{library_count} = $library_count;
    $lib_ref->{$lib_id}->{library_name} = $line[2];
    $lib_ref->{$lib_id}->{group} = $line[3];
    unless ($line[3] == 1 || $line[3] == 2){
      print RED, "\n\nGroup value for LIBRARY entries must be '1' or '2'! Check configuration file\n\n", RESET;
      exit();
    }
  }
  if ($line[0] =~ /^COMPARISON/){
    $comparison_count++;
    my $comp_ref = $project_conf{COMPARISON};
    $comp_ref->{$comparison_count}->{libraryA_id} = $line[1];
    $comp_ref->{$comparison_count}->{libraryA_name} = $line[2];
    $comp_ref->{$comparison_count}->{libraryB_id} = $line[3];
    $comp_ref->{$comparison_count}->{libraryB_name} = $line[4];
    $comp_ref->{$comparison_count}->{comparison_id} = $line[5];
  }
  if ($line[0] =~ /^GROUP/){
    $group_count++;
    my $group_ref = $project_conf{GROUP};
    $group_ref->{$group_count}->{group_name} = $line[1];
    $group_ref->{$group_count}->{library_list} = $line[2];
  }
  if ($line[0] =~ /^GRP_COMPARISON/){
    $group_comparison_count++;
    my $group_comp_ref = $project_conf{GROUP_COMPARISON};
    $group_comp_ref->{$group_comparison_count}->{group_list} = $line[2];
    $group_comp_ref->{$group_comparison_count}->{group_comparison_id} = $line[1];
  }
}

close(CONF);
my $lanes_ref = $project_conf{LANE};
my $libraries_ref = $project_conf{LIBRARY};
my $comparisons_ref = $project_conf{COMPARISON};
my $groups_ref = $project_conf{GROUP};
my $group_comparisons_ref = $project_conf{GROUP_COMPARISON};

#Make sure that all flowcell-lane combinations are unique!
my $fcl_count = keys %test;
unless ($fcl_count == $lane_count){
  print RED, "\n\nThe number of unique flowcell lanes does not equal the number of lanes of data specified (each flowcell+lane combination must be unique!)\n\n", RESET;
  exit();
}

#print Dumper %project_conf;
print BLUE, "\n\nFound $lane_count lanes of data for $library_count libraries in $group_count groups with $comparison_count library and $group_comparison_count group comparison(s) to be performed", RESET;

#Perform basic checks of the config files
print MAGENTA, "\n\nChecking config files for problems", RESET;
my $warnings = &checkConfig('-alexa_seq_conf'=>\%alexa_seq_conf, '-project_conf'=>\%project_conf);

my $region_file = "$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/Regions_"."$project_conf{REGIONS_FILE_VERSION}"."_Genes.txt";

#0-A.) Install the EnsEMBL API specified in the config file (if necessary)
print MAGENTA, "\n\n0-A.) Install the EnsEMBL API specified in the config file (if necessary)", RESET;
print COMMANDS "\n\n\n#0-A.) Install the EnsEMBL API specified in the config file (if necessary)";
print COMMANDS "\n#Execute the following command to install the EnsEMBL API for version: $project_conf{ENSEMBL_VERSION}";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/installEnsemblAPI.pl  --install_dir=$script_dir/ensembl_api/  --ensembl_version=$project_conf{ENSEMBL_VERSION}";


#0-B.) Install the specified ALEXA-Seq annotation database (if neccessary)
print MAGENTA, "\n\n0-B.) Install the specified ALEXA-Seq annotation database (if neccessary)", RESET;
print COMMANDS "\n\n\n#0-B.) Install the specified ALEXA-Seq annotation database (if neccessary)";
print COMMANDS "\n#Execute the following command to check/install the ALEXA-Seq annotation database: $project_conf{ALEXA_SEQ_DB}";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/installAnnotationDb.pl  --annotation_dir=$alexa_seq_conf{SEQUENCE_DB_DIR}/  --db_build=$project_conf{ALEXA_SEQ_DB}  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER2}  --password=$alexa_seq_conf{ALEXA_PASSWORD2}";

#0-C.) Partition the annotation database (to improve memory and file I/O performance)
print MAGENTA, "\n\n#0-C.) Partition the annotation database (to improve memory and file I/O performance)", RESET;
print COMMANDS "\n\n\n#0-C.) Partition the annotation database (to improve memory and file I/O performance)";
print COMMANDS "\n#Execute the following command to partition the annotation database:";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/partitionAnnotationDb.pl  --region_file_version=$project_conf{REGIONS_FILE_VERSION}  --annotation_dir=$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/";

#1.) Create needed target directories
my $cmd1 = "$script_dir/process/createAnalysisDirs.pl  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --project=$project_conf{PROJECT_NAME}  --library_conf_file=$project_config_file  --ensembl_version=$project_conf{ENSEMBL_VERSION}";
print MAGENTA, "\n\n1.) Creating all necessary directories for this project - directories already present will be untouched!", RESET;
print COMMANDS "\n\n#The following commands should be run in order using the provided instructions";
print COMMANDS "\n\n\n#1.) Create needed target directories";
print COMMANDS "\n#This command should have already been run, but it is reproduced here for reference (no harm will be done be re-running it)";
print COMMANDS "\n$cmd1";

#This command actually needs to be executed as well as written to the commands file for reference
print BLUE, "\n\n$cmd1", RESET;
system($cmd1);

#2.) processRawSolexaReads.pl
print MAGENTA, "2.) Creating jobs to process raw reads", RESET;
print COMMANDS "\n\n\n#2.) Creating jobs to process raw reads\n#Concatenates *_seq.txt, *_qseq.txt or *_sequence.txt (fastq) files from source directory to create a combined master raw seq data file\n#These files are expected to be named as follows: s_1_1_0001_qseq.txt.bz2 (s_<Lane>_<Read-1-or-2>_Tile_<seq-or-qseq>.txt.<gz-or-bz2> OR s_<lane>_<Read-1-or-2>_sequence.txt.<gz-or-bz2>)\n#NOTE that the format of sequence data in the raw .seq file is as follows:\n#Lane  Tile X-coord    Y-Coord  Sequence\n#Bases which could not be resolved are represented by a '.'\n#File format descriptions for qseq and fastq files are available online\n#qseq: http://jumpgate.caltech.edu/wiki/QSeq\n#fastq: http://en.wikipedia.org/wiki/Fastq";

my $bash_file1 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/process/processRawSolexaReads.sh";
open(BASH, ">$bash_file1") || die "\n\nCould not open bash output file: $bash_file1\n\n";
foreach my $lane_count (sort {$a <=> $b} keys %{$lanes_ref}){
  my $source_dir = $lanes_ref->{$lane_count}->{data_path};
  my $flowcell_name = $lanes_ref->{$lane_count}->{flowcell_name};
  my $lane_number = $lanes_ref->{$lane_count}->{lane_number};
  my $library_id = $lanes_ref->{$lane_count}->{library_id};
  my $read_trim = $lanes_ref->{$lane_count}->{read_trim};
  my $seq_file_type = $lanes_ref->{$lane_count}->{seq_file_type};
  my $max_n = $lanes_ref->{$lane_count}->{max_n};
  my $min_phred = $lanes_ref->{$lane_count}->{min_phred};
  my $qual_type = $lanes_ref->{$lane_count}->{qual_type};
  my $read_length = $lanes_ref->{$lane_count}->{read_length};

  #Set low complexity cutoff to be 2/3 of read length
  my $lc_cutoff = sprintf("%.0f", ($read_length*(2/3)));

  my $cmd = "$script_dir/process/processRawSolexaReads.pl  --input_dir=$source_dir  --flowcell_name=$flowcell_name  --lane_number=$lane_number  --raw_read_dir=$alexa_seq_conf{ANALYSIS_DIR}/raw_seq_data/$library_id/  --read_record_dir=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/  --mdust_bin=$alexa_seq_conf{MDUST_BIN}  --temp_dir=$alexa_seq_conf{ANALYSIS_DIR}/temp/mdust/  --trim=$read_trim  --low_complexity_cutoff=$lc_cutoff  --n_cutoff=$max_n  --phred_cutoff=$min_phred  --file_type=$seq_file_type  --qual_type=$qual_type  --log_dir=$alexa_seq_conf{ANALYSIS_DIR}/logs/$library_id/  --force=yes";

  print BASH "$cmd\n";
}
close(BASH);
if ($cluster_commands){
  print COMMANDS "\n#Execute the following to submit parsing jobs to the cluster";
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}\ncd $alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/process/";
  print COMMANDS "\nmqsub  --file $bash_file1  --name prsr  --mkdir  --delay $alexa_seq_conf{VLONG_DELAY}  --qsub '-l mem_total=15G'";
}else{
  print COMMANDS "\n#The following shell script can either be run as a single job or the jobs contained can be submitted to a computer cluster";
  print COMMANDS "\nbash $bash_file1";
}

print COMMANDS "\n\n#STOP.  Take a look at the read records files created in the previous step for the libraries of these projects";
print COMMANDS "\n#Make sure the read IDs seem correct (i.e. FLOWCELL_LANE_TILE_X_Y)";


#3.) Gather basic info about each lane, library and comparison. Then get number of quality reads for each library and the average read length add this info to the appropriate records
my $project_lib_data_file = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$project_conf{PROJECT_NAME}_Lib_Data.txt";
my $project_lib_names_file = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$project_conf{PROJECT_NAME}_Lib_Names.txt";
my $project_lib_comparisons_file = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$project_conf{PROJECT_NAME}_Lib_Comparisons.txt";
my $project_lib_groups_file = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$project_conf{PROJECT_NAME}_Lib_Groups.txt";
my $project_group_comparisons_file = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$project_conf{PROJECT_NAME}_Group_Comparisons.txt";

print COMMANDS "\n\n\n#3.) Gathering basic info about each lane, library and comparison. Then get number of quality reads for each library and the average read length add this info to the appropriate records";
print COMMANDS "\n\n#Note: All of the individual steps of Step 3 and 4 can be run concurrently\n";
print MAGENTA, "\n\n3.) Gather basic info and creating jobs to determine the number of quality reads in each lane as well as the average read length", RESET;

#3-A.) Gathering basic info about each lane of data
print COMMANDS "\n#3-A.) Gathering basic info about each lane of data\n";
print MAGENTA, "\n\t3-A.) Gathering basic info about each lane of data", RESET;
print COMMANDS "#The following file stores this basic info with the following values:\n";
print COMMANDS "#library_id\tflowcell_name_lane\tflowcell_name_lane2\tflowcell_name\tlane_number\tread_length\tdouble_read_length\ttotal_quality_reads\tsource_dir\n";
print COMMANDS "#$project_lib_data_file\n";

#First create a tab delimited 'Project_Lib_Data.txt' file and store in the '/batch_jobs/Project/' dir
#Only if it is not already present
if ((! -e $project_lib_data_file) || ($clean == 1)){
  open (PLD, ">$project_lib_data_file") || die "\n\nCould not open project lib data file: $project_lib_data_file\n\n";
  foreach my $lane_count (sort {$a <=> $b} keys %{$lanes_ref}){
    my $source_dir = $lanes_ref->{$lane_count}->{data_path};
    my $flowcell_name = $lanes_ref->{$lane_count}->{flowcell_name};
    my $lane_number = $lanes_ref->{$lane_count}->{lane_number};
    my $flowcell_name_lane = "$flowcell_name"."_"."$lane_number";
    my $flowcell_name_lane2 = "$flowcell_name"."_Lane"."$lane_number";
    my $library_id = $lanes_ref->{$lane_count}->{library_id};
    my $read_length =  $lanes_ref->{$lane_count}->{read_length};
    my $double_read_length = $read_length*2;
    print PLD "$library_id\t$flowcell_name_lane\t$flowcell_name_lane2\t$flowcell_name\t$lane_number\t$read_length\t$double_read_length\tNA\t$source_dir\n";
  }
  close(PLD);
}

#3-B.) Gathering basic info about each library
#Next create a tab delimited 'Project_Lib_Names.txt' file and store in the '/batch_jobs/Project/' dir
#Only if it is not already present
print COMMANDS "\n#3-B.) Gathering basic info about each library\n";
print MAGENTA, "\n\t3-B.) Gathering basic info about each library", RESET;
print COMMANDS "#The following file stores this basic info with the following values:\n";
print COMMANDS "#library_id\tlibrary_name\tgroup_number\tjunction_correction_factor\taverage_read_length\n";
print COMMANDS "#$project_lib_names_file\n";
if ((! -e $project_lib_names_file) || ($clean == 1)){
  open (PLN, ">$project_lib_names_file") || die "\n\nCould not open project lib names file: $project_lib_names_file\n\n";
  foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
    my $library_name = $libraries_ref->{$library_id}->{library_name};
    my $group = $libraries_ref->{$library_id}->{group};
    print PLN "$library_id\t$library_name\t$group\t1\t0\n";
  }
  close(PLN);
}

#3-C.) Gathering basic info about each comparison
#Next create a tab delimited 'Project_Lib_Comparisons.txt' file and store in the '/batch_jobs/Project/' dir
#Only if it is not already present
print COMMANDS "\n#3-C.) Gathering basic info about each comparison\n";
print MAGENTA, "\n\t3-C.) Gathering basic info about each comparison", RESET;
print COMMANDS "#The following file stores this basic info with the following values:\n";
print COMMANDS "#library_A_id\tlibrary_A_name\tlibrary_B_id\tlibrary_B_name\tcomparison_name\tcomparison_id\n";
print COMMANDS "#$project_lib_comparisons_file\n";
if ((! -e $project_lib_comparisons_file) || ($clean == 1)){
  open (PLC, ">$project_lib_comparisons_file") || die "\n\nCould not open project lib comparisons file: $project_lib_comparisons_file\n\n";
  foreach my $comparison (sort {$a <=> $b} keys %{$comparisons_ref}){
    my $libA_id = $comparisons_ref->{$comparison}->{libraryA_id};
    my $libA_name = $comparisons_ref->{$comparison}->{libraryA_name};
    my $libB_id = $comparisons_ref->{$comparison}->{libraryB_id};
    my $libB_name = $comparisons_ref->{$comparison}->{libraryB_name};
    my $comparison_id = $comparisons_ref->{$comparison}->{comparison_id};
    my $comparison_name = $comparison_id;
    $comparison_name =~ s/_/ /g;
    $comparisons_ref->{$comparison}->{comparison_name} = $comparison_name;
    print PLC "$libA_id\t$libA_name\t$libB_id\t$libB_name\t$comparison_name\t$comparison_id\n";
  }
  close(PLC);
}

#3-D.) Gathering basic info about each group
#Next create a tab delimited 'Project_Lib_Groups.txt' file and store in the '/batch_jobs/Project/' dir
#Only if it is not already present
print COMMANDS "\n#3-D.) Gathering basic info about each group (if any defined)\n";
print MAGENTA, "\n\t3-D.) Gathering basic info about each group (if any defined)", RESET;
print COMMANDS "#The following file stores this basic info with the following values:\n";
print COMMANDS "#group_name\tlibrary_list\n";
print COMMANDS "#$project_lib_groups_file\n";
if ((! -e $project_lib_groups_file) || ($clean == 1)){
  open (PLG, ">$project_lib_groups_file") || die "\n\nCould not open project lib groups file: $project_lib_groups_file\n\n";
  foreach my $group (sort {$a <=> $b} keys %{$groups_ref}){
    my $group_name = $groups_ref->{$group}->{group_name};
    my $library_list = $groups_ref->{$group}->{library_list};
    print PLG "$group_name\t$library_list\n";
  }
  close(PLG);
}

#3-E.) Gathering basic info about each group comparison
#Next create a tab delimited 'Project_Group_Comparisons.txt' file and store in the '/batch_jobs/Project/' dir
#Only if it is not already present
print COMMANDS "\n#3-E.) Gathering basic info about each group comparison\n";
print MAGENTA, "\n\t3-E.) Gathering basic info about each group comparison", RESET;
print COMMANDS "#The following file stores this basic info with the following values:\n";
print COMMANDS "#group_comparison_name\tgroup_list\n";
print COMMANDS "#$project_group_comparisons_file\n";
if ((! -e $project_group_comparisons_file) || ($clean == 1)){
  open (PGC, ">$project_group_comparisons_file") || die "\n\nCould not open project group comparisons file: $project_group_comparisons_file\n\n";
  foreach my $group_comparison (sort {$a <=> $b} keys %{$group_comparisons_ref}){
    my $group_list = $group_comparisons_ref->{$group_comparison}->{group_list};
    my $group_comparison_id = $group_comparisons_ref->{$group_comparison}->{group_comparison_id};
    my $group_comparison_name = $group_comparison_id;
    $group_comparison_name =~ s/_/ /g;
    $group_comparisons_ref->{$group_comparison}->{group_comparison_name} = $group_comparison_name;
    print PGC "$group_list\t$group_comparison_name\t$group_comparison_id\n";
  }
  close(PGC);
}

#3-F.) Generate basic statistics for each library (quality read counts, average read length, etc.)
print COMMANDS "\n#3-F.) Generate basic statistics for each library\n";
print MAGENTA, "\n\t3-F.) Generate basic statistics for each library", RESET;
#Create jobs to get quality read counts
#Then merge these quality read counts into the Project_Lib_Data file
my $bash_file2 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/process/getQualityReadCounts.sh";
open(BASH, ">$bash_file2") || die "\n\nCould not open bash output file: $bash_file2\n\n";
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  #Print the quality read counts to a separate file
  my $result_file = "$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/Summary/qualityReadCounts.txt";
  print BASH "$script_dir/stats/getQualityReadCounts.pl  --read_record_dir=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/  --result_file=$result_file\n";
  print BASH "$script_dir/stats/mergeQualityReadCounts.pl  --lib_data_file=$project_lib_data_file  --quality_counts_file=$result_file\n"
}
close(BASH);
print COMMANDS "\n#Execute the following shell script to determine the number of quality read counts in each lane";
print COMMANDS "\nbash $bash_file2";


#Create jobs to get average read lengths for each library
#Then merge these quality read counts into the Project_Lib_Names file
my $bash_file4 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/process/getAverageReadLength.sh";
open(BASH, ">$bash_file4") || die "\n\nCould not open bash output file: $bash_file4\n\n";
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  #Print the quality read counts to a separate file
  my $result_file = "$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/Summary/averageReadLength.txt";
  print BASH "$script_dir/stats/getAverageReadLength.pl  --read_record_dir=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/  --result_file=$result_file  --library_id=$library_id\n";
  print BASH "$script_dir/stats/mergeAverageReadLength.pl  --lib_names_file=$project_lib_names_file  --average_length_file=$result_file\n"
}
close(BASH);
print COMMANDS "\n\n#Execute the following shell script to determine the average read length of the library";
print COMMANDS "\nbash $bash_file4";


#3-G.) Create Phred score CDF plots
print MAGENTA, "\n\t3-G.) Create Phred score CDF plots", RESET;

my $bash_file5c = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/stats/createPhredCDFplots.sh";
open(BASH, ">$bash_file5c") || die "\n\nCould not open bash output file: $bash_file5c\n\n";

#Gather flowcell_lanes for all libraries
my %library_lanes;
foreach my $lane_count (sort {$a <=> $b} keys %{$lanes_ref}){
  my $flowcell_name = $lanes_ref->{$lane_count}->{flowcell_name};
  my $lane_number = $lanes_ref->{$lane_count}->{lane_number};
  my $flowcell_name_lane = "$flowcell_name"."_"."$lane_number";
  my $library_id = $lanes_ref->{$lane_count}->{library_id};
  my $type = $lanes_ref->{$lane_count}->{seq_file_type};
  unless($type eq "seq"){
    $library_lanes{$library_id}{$flowcell_name_lane}++;
  }
}

foreach my $library (sort keys %library_lanes){
  my $phred_dir = "$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library/Summary/";
  my $result_dir = "$alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/$library/LibraryQuality/";
  my @lanes;
  foreach my $lane (sort keys %{$library_lanes{$library}}){
    push (@lanes, $lane);
  }
  print BASH "$script_dir/R_bin/summarizePhred.R $phred_dir $result_dir ",join(",",@lanes)," $library\n";
}
close BASH;

print COMMANDS "\n\n#Execute the following shell script to create Phred score CDF plots";
print COMMANDS "\nbash $bash_file5c";


#3-H.) Determine alignment parameters (lane-by-lane) to be used in the mapping and parsing jobs
print MAGENTA, "\n\t3-H.) Determine alignment parameters (lane-by-lane) to be used in the mapping and parsing jobs", RESET;

#Determine the word_size, min_bit_score, upper_bit_score, min_align_length etc. for each lane
#Note: most repeat elements have a min size of 70 bp so for these features assume READ_LENGTH does not exceed 70 for parameter calculations
#Note: for junctions/boundaries, assume READ_LENGTH does not exceed JUNCTION-BOUNDARY_SEQ_LENGTH for parameter calculations

#MISMATCHES_ALLOWED -> Provided in table file.  Allow 1 per 20 bp of read length.  (i.e. READ_LENGTH/20 ... then round up)
#READ_LENGTH        -> Supplied by user in config file
#MIN_OVERLAP        -> Supplied by user in config file (default is 6)
#N                  -> MIN_OVERLAP - 2
#MIN_ALIGN_LENGTH   -> READ_LENGTH-N (not used for junctions/boundaries)
#MIN_BIT_SCORE      -> For hits to transcripts, introns, intergenics: get from table file
#                   -> For repeats: get from table file but make sure it does not exceed 91.7 (corresponds to 2/3 of min length of repeat elements) 
#                   -> For junctions/boundaries: get from table file but does not exceed the PERFECT_BIT_SCORE corresponding to 1/2 JUNCTION-BOUNDARY_SEQ_LENGTH + MIN_OVERLAP (eg. (62/2)+6=37 -> bit score = 73.8)
#UPPER_BIT_SCORE    -> For hits to transcripts, introns, intergenics: get from table file
#                   -> For repeats: get from table file but make sure it does not exceed cutoff for a 70-mer (i.e. 107.0)
#                   -> For junctions/boundaries: get from table file but does not exceed the PERFECT_BIT_SCORE corresponding to 1/2 JUNCTION-BOUNDARY_SEQ_LENGTH + MIN_OVERLAP (eg. (62/2)+6=37 -> bit score = 73.8)
#WORD_SIZE          -> For hits to transcripts, introns, intergenics: (READ_LENGTH/ALLOWED_MISMATCHES+1) ... then round down.  No smaller than 11. 
#                   -> For repeats: Same as above BUT no larger than (70/ALLOWED_MISMATCHES+1) ... then round down
#                   -> For junctions/boundaries: Same as above BUT no larger than ((JUNCTION-BOUNDARY_SEQ_LENGTH)/ALLOWED_MISMATCHES+1) ... then round down 

#First import the table file containing read-length to parameter relationships
my $table_file = "$script_dir/utilities/bitScore_Vs_RL_Table.txt";
my %table;
open (TABLE, "$table_file") || die "\n\nCould not open table file: $table_file\n\n";
while(<TABLE>){
  chomp($_);
  if ($_ =~ /^\#/){
    next();
  }
  unless ($_ =~ /^\d+/){
    next();
  }
  my @line = split(/ +/, $_);
  $table{$line[0]}{allowed_mismatches} = $line[1];
  $table{$line[0]}{two_thirds_length} = $line[2];
  $table{$line[0]}{min_bit_score} = $line[3];
  $table{$line[0]}{upper_bit_score} = $line[4];
  $table{$line[0]}{perfect_bit_score} = $line[5];
}
close(TABLE);

my $min_overlap = $project_conf{MIN_OVERLAP};
my $junction_seq_length = $project_conf{JUNCTION_SEQ_LENGTH};
my $boundary_seq_length = $project_conf{BOUNDARY_SEQ_LENGTH};
my $n = $min_overlap-2;
foreach my $lane_count (sort {$a <=> $b} keys %{$lanes_ref}){
  my $flowcell_name = $lanes_ref->{$lane_count}->{flowcell_name};
  my $lane_number = $lanes_ref->{$lane_count}->{lane_number};
  my $flowcell_name_lane2 = "$flowcell_name"."_Lane"."$lane_number";

  my $min_repeat_length = 70;
  my $test_repeat_length = sprintf("%.0f", ($min_repeat_length*0.67));
  my $read_length = $lanes_ref->{$lane_count}->{read_length} - $lanes_ref->{$lane_count}->{read_trim};  # If read length is going to be trimmed, correct for this now!
  my $test_junction_length = sprintf("%.0f", (($junction_seq_length/2)+$min_overlap));
  my $test_boundary_length = sprintf("%.0f", (($boundary_seq_length/2)+$min_overlap));

  #Make sure all lengths are defined in the table file
  unless ($table{$read_length} && $table{$min_repeat_length} && $table{$test_repeat_length} && $table{$test_junction_length} && $table{$test_boundary_length}){
    print RED, "\n\nOne of the lengths being used to define parameters is not defined in the table file!!  - Aborting", RESET;
    print RED, "\n\nCheck you project config file ($project_config_file) and the table file ($table_file)", RESET;
    exit();
  }
  my $allowed_mismatches = $table{$read_length}{allowed_mismatches};
  my $allowed_mismatches_r = $table{$read_length}{allowed_mismatches};
  if ($allowed_mismatches_r > $table{$min_repeat_length}{allowed_mismatches}){$allowed_mismatches_r = $table{$min_repeat_length}{allowed_mismatches};}
  my $allowed_mismatches_j = $table{$read_length}{allowed_mismatches};
  if ($allowed_mismatches_j > $table{$junction_seq_length}{allowed_mismatches}){$allowed_mismatches_j = $table{$junction_seq_length}{allowed_mismatches};}
  my $allowed_mismatches_b = $table{$read_length}{allowed_mismatches};
  if ($allowed_mismatches_b > $table{$boundary_seq_length}{allowed_mismatches}){$allowed_mismatches_b = $table{$boundary_seq_length}{allowed_mismatches};}

  my $min_align_length = $read_length-$n;
  my $min_align_length_r = $min_align_length;
  if ($min_align_length_r > $min_repeat_length-$n){
    $min_align_length_r = $min_repeat_length-$n;
  }

  my $min_bit_score = $table{$read_length}{min_bit_score};
  my $min_bit_score_r = $table{$read_length}{min_bit_score};
  if ($min_bit_score_r > $table{$test_repeat_length}{min_bit_score}){$min_bit_score_r = $table{$test_repeat_length}{min_bit_score};}
  my $min_bit_score_j = $table{$read_length}{min_bit_score};
  my $min_bit_score_j_test = $table{$test_junction_length}{perfect_bit_score};
  if ($min_bit_score_j > $min_bit_score_j_test){$min_bit_score_j = $min_bit_score_j_test;}
  my $min_bit_score_b = $table{$read_length}{min_bit_score};
  my $min_bit_score_b_test = $table{$test_boundary_length}{perfect_bit_score};
  if ($min_bit_score_b > $min_bit_score_b_test){$min_bit_score_b = $min_bit_score_b_test;}

  my $upper_bit_score = $table{$read_length}{upper_bit_score};
  my $upper_bit_score_r = $table{$read_length}{upper_bit_score};
  if ($upper_bit_score_r > $table{$min_repeat_length}{upper_bit_score}){$upper_bit_score_r = $table{$min_repeat_length}{upper_bit_score};}
  my $upper_bit_score_j = $table{$read_length}{upper_bit_score};
  if ($upper_bit_score_j > $min_bit_score_j_test){$upper_bit_score_j = $min_bit_score_j_test;}
  my $upper_bit_score_b = $table{$read_length}{upper_bit_score};
  if ($upper_bit_score_b > $min_bit_score_b_test){$upper_bit_score_b = $min_bit_score_b_test;}

  my $word_size = sprintf("%.0f", ($read_length/($allowed_mismatches+1)));
  if ($word_size < 11){$word_size = 11;}

  my $word_size_r;
  if ($read_length <= $min_repeat_length){
    $word_size_r = sprintf("%.0f", ($read_length/($allowed_mismatches+1)));
  }else{
    $word_size_r = sprintf("%.0f", ($min_repeat_length/($allowed_mismatches_r+1)));
  }
  if ($word_size_r < 11){$word_size_r = 11;}

  my $word_size_j;
  if ($read_length <= $junction_seq_length){
    $word_size_j = sprintf("%.0f", ($read_length/($allowed_mismatches+1)));
  }else{
    $word_size_j = sprintf("%.0f", ($junction_seq_length/($allowed_mismatches_j+1)));
  }
  if ($word_size_j < 11){$word_size_j = 11;}

  my $word_size_b;
  if ($read_length <= $boundary_seq_length){
    $word_size_b = sprintf("%.0f", ($read_length/($allowed_mismatches+1)));
  }else{
    $word_size_b = sprintf("%.0f", ($boundary_seq_length/($allowed_mismatches_b+1)));
  }
  if ($word_size_b < 11){$word_size_b = 11;}

  #Store the lane-by-lane values for later
  $lanes_ref->{$lane_count}->{min_align_length} = $min_align_length;
  $lanes_ref->{$lane_count}->{min_align_length_r} = $min_align_length_r;
  $lanes_ref->{$lane_count}->{min_bit_score} = $min_bit_score;
  $lanes_ref->{$lane_count}->{min_bit_score_r} = $min_bit_score_r;
  $lanes_ref->{$lane_count}->{min_bit_score_j} = $min_bit_score_j;
  $lanes_ref->{$lane_count}->{min_bit_score_b} = $min_bit_score_b;
  $lanes_ref->{$lane_count}->{upper_bit_score} = $upper_bit_score;
  $lanes_ref->{$lane_count}->{upper_bit_score_r} = $upper_bit_score_r;
  $lanes_ref->{$lane_count}->{upper_bit_score_j} = $upper_bit_score_j;
  $lanes_ref->{$lane_count}->{upper_bit_score_b} = $upper_bit_score_b;
  $lanes_ref->{$lane_count}->{word_size} = $word_size;
  $lanes_ref->{$lane_count}->{word_size_r} = $word_size_r;
  $lanes_ref->{$lane_count}->{word_size_j} = $word_size_j;
  $lanes_ref->{$lane_count}->{word_size_b} = $word_size_b;

  #Summarize the parameters set for this lane
  #print YELLOW, "\n\n\t$flowcell_name_lane2\tRead Length=$read_length\tMin. Overlap = $min_overlap\tMin. Align. Length = $min_align_length\tParameters (general, repeats, junctions, boundaries)", RESET;
  #print YELLOW, "\n\t\tAllowed Mismatches: ($allowed_mismatches $allowed_mismatches_r $allowed_mismatches_j $allowed_mismatches_b)", RESET;
  #print YELLOW, "\n\t\tMin Bit Score: ($min_bit_score $min_bit_score_r $min_bit_score_j $min_bit_score_b)", RESET;
  #print YELLOW, "\n\t\tUpper Bit Score: ($upper_bit_score $upper_bit_score_r $upper_bit_score_j $upper_bit_score_b)", RESET;
  #print YELLOW, "\n\t\tWord Size: ($word_size $word_size_r $word_size_j $word_size_b)", RESET;
}


#4.) Create a fasta file of individual reads for all lanes of data
print COMMANDS "\n\n\n#4.) Create a fasta file of individual reads for all lanes of data";
print MAGENTA, "\n\n4.) Create a fasta file of individual reads for all lanes of data", RESET;

my $bash_file6 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/process/createSolexaReadFasta.sh";
open(BASH, ">$bash_file6") || die "\n\nCould not open bash output file: $bash_file6\n\n";
foreach my $lane_count (sort {$a <=> $b} keys %{$lanes_ref}){
  my $flowcell_name = $lanes_ref->{$lane_count}->{flowcell_name};
  my $lane_number = $lanes_ref->{$lane_count}->{lane_number};
  my $flowcell_name_lane2 = "$flowcell_name"."_Lane"."$lane_number";
  my $library_id = $lanes_ref->{$lane_count}->{library_id};
 
  my $cmd = "$script_dir/process/createSolexaReadFasta.pl  --read_file=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/$flowcell_name_lane2.txt.gz  --require_paired_reads=no  --exclusion_classes=\"Low_Quality Low_Complexity Duplicate\"  --fasta_file=$alexa_seq_conf{ANALYSIS_DIR}/fasta_seq_data/$library_id/$flowcell_name_lane2/$flowcell_name_lane2"."_QualityFiltered_Unpaired.fa  --log_file=$alexa_seq_conf{ANALYSIS_DIR}/logs/$library_id/$flowcell_name_lane2/createSolexaReadFasta_1_LOG.txt";
  print BASH "$cmd\n"
}
close(BASH);
if ($cluster_commands){
  print COMMANDS "\n#Execute the following to create a fasta file for each lane of data using the cluster";
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}\ncd $alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/process/";
  print COMMANDS "\nmqsub  --file $bash_file6  --name csrf  --mkdir  --delay $alexa_seq_conf{VLONG_DELAY}";
}else{
  print COMMANDS "\n#Execute the following shell script to create a fasta file for each lane of data";
  print COMMANDS "\nbash $bash_file6";
}

#5.) Create jobs to perform MAPPING of reads to various feature types
print MAGENTA, "\n\n5.) MAPPING", RESET;
print COMMANDS "\n\n\n#5.) MAPPING";
my %types;
$types{1}{name} = "repeats";
$types{1}{blast_db_name} = "repeats";
$types{1}{bwa_db_name} = "repeats.fa.gz";
$types{1}{abr} = "BR";
$types{2}{name} = "transcripts";
$types{2}{blast_db_name} = "transcripts";
$types{2}{bwa_db_name} = "transcripts.fa.gz";
$types{2}{abr} = "BT";
$types{3}{name} = "exonJunctions";
$types{3}{blast_db_name} = "exonJunctions_"."$project_conf{JUNCTION_SEQ_LENGTH}"."mers";
$types{3}{bwa_db_name} = "exonJunctions_"."$project_conf{JUNCTION_SEQ_LENGTH}"."mers.fa.gz";
$types{3}{abr} = "BEJ";
$types{4}{name} = "exonBoundaries";
$types{4}{blast_db_name} = "exonBoundaries_"."$project_conf{BOUNDARY_SEQ_LENGTH}"."mers";
$types{4}{bwa_db_name} = "exonBoundaries_"."$project_conf{BOUNDARY_SEQ_LENGTH}"."mers.fa.gz";
$types{4}{abr} = "BEB";
$types{5}{name} = "introns";
$types{5}{blast_db_name} = "introns";
$types{5}{bwa_db_name} = "introns.fa.gz";
$types{5}{abr} = "BI";
$types{6}{name} = "intergenics";
$types{6}{blast_db_name} = "intergenics";
$types{6}{bwa_db_name} = "intergenics.fa.gz";
$types{6}{abr} = "BIG";


#5-A.) Create jobs to perform mapping reads to various feature types
print COMMANDS "\n\n#5-A.) Create jobs to perform mapping reads to various feature types";
print MAGENTA, "\n\t5-A.) Create jobs to perform mapping reads to various feature types", RESET;
my $bash_file7 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/process/createMapBatch.sh";
open(BASH, ">$bash_file7") || die "\n\nCould not open bash output file: $bash_file7\n\n";
foreach my $tc (sort {$a <=> $b} keys %types){
  my $type_name = $types{$tc}{name};
  my $blast_db_name = $types{$tc}{blast_db_name};
  foreach my $lane_count (sort {$a <=> $b} keys %{$lanes_ref}){
    my $flowcell_name = $lanes_ref->{$lane_count}->{flowcell_name};
    my $library_id = $lanes_ref->{$lane_count}->{library_id};
    my $lane_number = $lanes_ref->{$lane_count}->{lane_number};
    my $flowcell_name_lane2 = "$flowcell_name"."_Lane"."$lane_number";

    my $min_bit_score = $lanes_ref->{$lane_count}->{min_bit_score};
    my $word_size = $lanes_ref->{$lane_count}->{word_size};
    my $overlap_options = '';

    if ($type_name =~ /repeats/){
      $min_bit_score = $lanes_ref->{$lane_count}->{min_bit_score_r};
      $word_size = $lanes_ref->{$lane_count}->{word_size_r};
    }
    if ($type_name =~ /exonJunctions/){
      $min_bit_score = $lanes_ref->{$lane_count}->{min_bit_score_j};
      $word_size = $lanes_ref->{$lane_count}->{word_size_j};
      $overlap_options = "--target_size=$project_conf{JUNCTION_SEQ_LENGTH}  --min_overlap=$project_conf{MIN_OVERLAP}";
    }
    if ($type_name =~ /exonBoundaries/){
      $min_bit_score = $lanes_ref->{$lane_count}->{min_bit_score_b};
      $word_size = $lanes_ref->{$lane_count}->{word_size_b};
      $overlap_options = "--target_size=$project_conf{BOUNDARY_SEQ_LENGTH}  --min_overlap=$project_conf{MIN_OVERLAP}";
    }

    if ($project_conf{ALIGNMENT_OPTION}==3){
      #Use BWA for all alignments
      my $bwa_db_name = $types{$tc}{bwa_db_name};
      $overlap_options = '';
      print BASH "$script_dir/process/createMapBatch.pl  --input_fasta_file=$alexa_seq_conf{ANALYSIS_DIR}/fasta_seq_data/$library_id/$flowcell_name_lane2/$flowcell_name_lane2"."_QualityFiltered_Unpaired.fa.gz  --working_dir=$alexa_seq_conf{ANALYSIS_DIR}/fasta_seq_data/$library_id/$flowcell_name_lane2/fasta_blocks/  --map_bin=$alexa_seq_conf{BWA_BIN_DIR}/bwa  --map_database=$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/$type_name/bwadb/$bwa_db_name  --word_size=$word_size  --min_bit_score=$min_bit_score  --map_filter_script=$script_dir/utilities/filterBwaStream.pl  --batch_file=$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$library_id/blast_vs_"."$type_name/$flowcell_name_lane2.sh  --map_results_dir=$alexa_seq_conf{ANALYSIS_DIR}/blast_results/$library_id/$flowcell_name_lane2/$type_name/  --job_size=$alexa_seq_conf{JOB_SIZE}  $overlap_options  --aligner=bwa\n";

    }elsif ($project_conf{ALIGNMENT_OPTION}==1 || $type_name =~ /repeats|exonJunctions|exonBoundaries/){
      #Use BLAST for all alignments
      print BASH "$script_dir/process/createMapBatch.pl  --input_fasta_file=$alexa_seq_conf{ANALYSIS_DIR}/fasta_seq_data/$library_id/$flowcell_name_lane2/$flowcell_name_lane2"."_QualityFiltered_Unpaired.fa.gz  --working_dir=$alexa_seq_conf{ANALYSIS_DIR}/fasta_seq_data/$library_id/$flowcell_name_lane2/fasta_blocks/  --map_bin=$alexa_seq_conf{BLAST_BIN_DIR}/blastall  --map_database=$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/$type_name/blastdb/$blast_db_name  --word_size=$word_size  --min_bit_score=$min_bit_score  --map_filter_script=$script_dir/utilities/filterBlastStream.pl  --batch_file=$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$library_id/blast_vs_"."$type_name/$flowcell_name_lane2.sh  --map_results_dir=$alexa_seq_conf{ANALYSIS_DIR}/blast_results/$library_id/$flowcell_name_lane2/$type_name/  --job_size=$alexa_seq_conf{JOB_SIZE}  $overlap_options  --aligner=blast\n";

    }elsif($project_conf{ALIGNMENT_OPTION}==2){
      #Use BWA for alignments except for repeats, junctions and boundaries
      my $bwa_db_name = $types{$tc}{bwa_db_name};
      $overlap_options = '';
      print BASH "$script_dir/process/createMapBatch.pl  --input_fasta_file=$alexa_seq_conf{ANALYSIS_DIR}/fasta_seq_data/$library_id/$flowcell_name_lane2/$flowcell_name_lane2"."_QualityFiltered_Unpaired.fa.gz  --working_dir=$alexa_seq_conf{ANALYSIS_DIR}/fasta_seq_data/$library_id/$flowcell_name_lane2/fasta_blocks/  --map_bin=$alexa_seq_conf{BWA_BIN_DIR}/bwa  --map_database=$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/$type_name/bwadb/$bwa_db_name  --word_size=$word_size  --min_bit_score=$min_bit_score  --map_filter_script=$script_dir/utilities/filterBwaStream.pl  --batch_file=$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$library_id/blast_vs_"."$type_name/$flowcell_name_lane2.sh  --map_results_dir=$alexa_seq_conf{ANALYSIS_DIR}/blast_results/$library_id/$flowcell_name_lane2/$type_name/  --job_size=$alexa_seq_conf{JOB_SIZE}  $overlap_options  --aligner=bwa\n";

    }else{
      print RED, "\n\nUnrecognized alignment option.  Check project config file...\n\n", RESET;
      exit();
    }
  }
}
print COMMANDS "\n#Execute the following shell script to mapping jobs for each lane of data";
print COMMANDS "\nbash $bash_file7";
close(BASH);

#5-B.) Join lane-by-lane jobs into one master batch for each feature type for each library
print COMMANDS "\n\n#5-B.) Join lane-by-lane jobs into one master batch for each feature type for each library";
print MAGENTA, "\n\t5-B.) Join lane-by-lane jobs into one master batch for each feature type for each library", RESET;
my $bash_file8 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/process/joinMapBatch.sh";
open(BASH, ">$bash_file8") || die "\n\nCould not open bash output file: $bash_file8\n\n";
foreach my $tc (sort {$a <=> $b} keys %types){
  my $type_name = $types{$tc}{name};
  foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
    my $joined_batch = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$library_id/blast_vs_"."$type_name/blast_vs_"."$type_name.sh";
    my $repair_batch = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$library_id/blast_vs_"."$type_name/blast_vs_"."$type_name"."_REPAIR.sh";
    print BASH "rm -f $joined_batch; rm -f $repair_batch; cat $alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$library_id/blast_vs_"."$type_name/*.sh > $joined_batch\n";
  }
}
close(BASH);
print COMMANDS "\n#Execute the following shell script to join the jobs";
print COMMANDS "\nbash $bash_file8";

#5-C.) Create submission cluster jobs, one for each feature type for each library
print COMMANDS "\n\n#5-C.) Create cluster submission jobs, one for each feature type for each library";
print MAGENTA, "\n\t5-C.) Create cluster submission jobs, one for each feature type for each library", RESET;
my $bash_file9 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/process/submitMapBatch.sh";
open(BASH, ">$bash_file9") || die "\n\nCould not open bash output file: $bash_file9\n\n";
foreach my $tc (sort {$a <=> $b} keys %types){
  my $type_name = $types{$tc}{name};
  my $type_abr = $types{$tc}{abr};
  foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
    if ($cluster_commands){
      print BASH "cd $alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$library_id/blast_vs_"."$type_name/\nmqsub --file blast_vs_"."$type_name.sh --name $type_abr"."_"."$library_id --mkdir --delay $alexa_seq_conf{MEDIUM_DELAY}\n";
    }else{
      print BASH "bash $alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$library_id/blast_vs_"."$type_name/blast_vs_"."$type_name.sh\n";
    }
  }
}
close(BASH);
if ($cluster_commands){
  print COMMANDS "\n#Log in to the cluster head node and execute the following shell script to submit your jobs";
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}\nbash $bash_file9\nexit";
}else{
  print COMMANDS "\n#Execute the following or submit the contained jobs to your cluster:";
  print COMMANDS "\nbash $bash_file9";
}

#5-D.) Create jobs to check cluster mapping jobs for completion and identify failed jobs that need to be repeated (if any)
print COMMANDS "\n\n#5-D.) Create jobs to check cluster mapping jobs for completion and identify failed jobs that need to be repeated (if any)";
print MAGENTA, "\n\t5-D.) Create jobs to check cluster mapping jobs for completion and identify failed jobs that need to be repeated (if any)", RESET;
my $bash_file10 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/process/checkMapBatch.sh";
open(BASH, ">$bash_file10") || die "\n\nCould not open bash output file: $bash_file10\n\n";
foreach my $tc (sort {$a <=> $b} keys %types){
  my $type_name = $types{$tc}{name};
  foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
    my $joined_batch = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$library_id/blast_vs_"."$type_name/blast_vs_"."$type_name.sh";
    my $repair_batch = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$library_id/blast_vs_"."$type_name/blast_vs_"."$type_name"."_REPAIR.sh";
    print BASH "$script_dir/process/checkClusterResults.pl  --input_file=$joined_batch  --output_file=$repair_batch\n";
  }
}
close(BASH);
print COMMANDS "\n#Execute the following shell script to check for completion of mapping jobs";
print COMMANDS "\nbash $bash_file10";
print COMMANDS "\n\n#If repeat jobs are necessary then submit these now";
print COMMANDS "\n#Use the following command to count the total number of jobs that need to be rerun (if none are needed, no files will be found):";
print COMMANDS "\nwc -l $alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/*/*/blast*_REPAIR.sh";

#6.) Create jobs to perform PARSING of reads assignments for each feature types
print MAGENTA, "\n\n6.) PARSING", RESET;
print COMMANDS "\n\n\n#6.) PARSING";
print COMMANDS "\n\n#Note: These steps must be run in order.  i.e. Don't start 6-B until 6-A is complete.";

#If the Project_Lib_Data file has already been generated, grab the total read counts from there
open(PLD, "$project_lib_data_file") || die "\n\nCould not open project lib data file\n\n";
my %lib_data;
while(<PLD>){
  chomp($_);
  my @line = split("\t", $_);
  $lib_data{$line[2]}{total_read_count} = $line[7];
  unless ($lib_data{$line[2]}{total_read_count} =~ /\d+/){
    $lib_data{$line[2]}{total_read_count} = 5000000;
  }
}  
close(PLD);

my $alexa_db_name = "ALEXA_"."$project_conf{ALEXA_SEQ_DB}";
my $parse_file_dir = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/process/";
my $expression_file_dir = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/expression/";
my $reads_analyzed_per_block = ($alexa_seq_conf{JOB_SIZE}/2);
my $parse_file;
my $parse_file_name;
my $total_lane_count = keys %{$lanes_ref};

#6-A.) REPEATS
print COMMANDS "\n\n#6-A.) REPEATS";
print COMMANDS "\n#Create jobs to check parse mapping results for repeat hits";
print MAGENTA, "\n\t6-A.) Create jobs to check parse mapping results for repeat hits", RESET;
$parse_file_name = "parseRepeats.sh";
$parse_file = "$parse_file_dir/$parse_file_name";
open (PARSE, ">$parse_file") || die "\n\nCould not open parse file: $parse_file\n\n";
foreach my $lane_count (sort {$a <=> $b} keys %{$lanes_ref}){
  my $flowcell_name = $lanes_ref->{$lane_count}->{flowcell_name};
  my $lane_number = $lanes_ref->{$lane_count}->{lane_number};
  my $flowcell_name_lane = "$flowcell_name"."_Lane"."$lane_number";
  my $library_id = $lanes_ref->{$lane_count}->{library_id};
  my $total_read_count = $lib_data{$flowcell_name_lane}{total_read_count};
  my $min_bit_score = $lanes_ref->{$lane_count}->{min_bit_score_r};
  my $upper_bit_score = $lanes_ref->{$lane_count}->{upper_bit_score_r};
  my $min_align_length = $lanes_ref->{$lane_count}->{min_align_length_r};

  print PARSE "source $alexa_seq_conf{SHELL_FILE}; $script_dir/expression/parseBlastResults_Repeats.pl  --read_records_infile=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/$flowcell_name_lane.txt.gz  --blast_results_dir=$alexa_seq_conf{ANALYSIS_DIR}/blast_results/$library_id/$flowcell_name_lane/repeats/  --reads_analyzed_per_block=$reads_analyzed_per_block  --working_dir=$alexa_seq_conf{ANALYSIS_DIR}/blast_results/$library_id/$flowcell_name_lane/temp/  --total_reads_analyzed=$total_read_count  --summary_outfile=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/Repeats/$flowcell_name_lane"."_Repeats.txt  --min_bit_score=$min_bit_score  --upper_bit_score=$upper_bit_score  --min_align_length=$min_align_length  --logfile=$alexa_seq_conf{ANALYSIS_DIR}/logs/$library_id/$flowcell_name_lane/parseBlastResults_Repeat_LOG.txt\n";
}
close (PARSE);
if ($cluster_commands){
  print COMMANDS "\n#Execute the following to submit parsing jobs to the cluster";
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}\ncd $parse_file_dir";
  print COMMANDS "\nmqsub  --file $parse_file_name  --name pr  --mkdir  --delay $alexa_seq_conf{VLONG_DELAY}";
}else{
  print COMMANDS "\n#Execute the following or submit the contained jobs to your cluster";
  print COMMANDS "\nbash $parse_file_dir"."$parse_file_name";
}
print COMMANDS "\n\n#Do not proceed to the next step until this step is complete!";
print COMMANDS "\n#Total number of parsing jobs that should be complete = $total_lane_count  (Use the following command to check for job completion):";
print COMMANDS "\n$script_dir/process/checkParsingJobs.pl  --project_name=$project_conf{PROJECT_NAME}  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --ensembl_version=$project_conf{ENSEMBL_VERSION}\n";


#6-B.) TRANSCRIPTS
print COMMANDS "\n\n#6-B.) TRANSCRIPTS";
print COMMANDS "\n#Create jobs to check parse mapping results for transcript hits";
print MAGENTA, "\n\t6-B.) Create jobs to check parse mapping results for transcript hits", RESET;
$parse_file_name = "parseTranscripts.sh";
$parse_file = "$parse_file_dir/$parse_file_name";
open (PARSE, ">$parse_file") || die "\n\nCould not open parse file: $parse_file\n\n";
foreach my $lane_count (sort {$a <=> $b} keys %{$lanes_ref}){
  my $flowcell_name = $lanes_ref->{$lane_count}->{flowcell_name};
  my $lane_number = $lanes_ref->{$lane_count}->{lane_number};
  my $flowcell_name_lane = "$flowcell_name"."_Lane"."$lane_number";
  my $library_id = $lanes_ref->{$lane_count}->{library_id};
  my $total_read_count = $lib_data{$flowcell_name_lane}{total_read_count};
  my $min_bit_score = $lanes_ref->{$lane_count}->{min_bit_score};
  my $upper_bit_score = $lanes_ref->{$lane_count}->{upper_bit_score};
  my $min_align_length = $lanes_ref->{$lane_count}->{min_align_length};

  print PARSE "source $alexa_seq_conf{SHELL_FILE}; $script_dir/expression/parseBlastResults_GeneExon.pl  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --working_dir=$alexa_seq_conf{ANALYSIS_DIR}/blast_results/$library_id/$flowcell_name_lane/temp/  --read_records_infile=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/$flowcell_name_lane.txt.gz  --blast_results_dir=$alexa_seq_conf{ANALYSIS_DIR}/blast_results/$library_id/$flowcell_name_lane/transcripts/  --reads_analyzed_per_block=$reads_analyzed_per_block  --total_reads_analyzed=$total_read_count  --summary_outfile=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/ENST_v$project_conf{ENSEMBL_VERSION}/$flowcell_name_lane"."_ENST_v$project_conf{ENSEMBL_VERSION}.txt  --filter_pseudogene_hits=no  --min_bit_score=$min_bit_score  --upper_bit_score=$upper_bit_score  --min_align_length=$min_align_length  --logfile=$alexa_seq_conf{ANALYSIS_DIR}/logs/$library_id/$flowcell_name_lane/parseBlastResults_GeneExon_LOG.txt\n";
}
close (PARSE);
if ($cluster_commands){
  print COMMANDS "\n#Execute the following to submit parsing jobs to the cluster";
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}\ncd $parse_file_dir";
  print COMMANDS "\nmqsub  --file $parse_file_name  --qsub '-l mem_total=7G'  --name pt  --mkdir  --delay $alexa_seq_conf{VLONG_DELAY}";
}else{
  print COMMANDS "\n#Execute the following or submit the contained jobs to your cluster";
  print COMMANDS "\nbash $parse_file_dir"."$parse_file_name";
}
print COMMANDS "\n\n#Do not proceed to the next step until this step is complete!";
print COMMANDS "\n#Total number of parsing jobs that should be complete = $total_lane_count  (Use the following command to check for job completion):";
print COMMANDS "\n$script_dir/process/checkParsingJobs.pl  --project_name=$project_conf{PROJECT_NAME}  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --ensembl_version=$project_conf{ENSEMBL_VERSION}\n";


#6-C.) JUNCTIONS
print COMMANDS "\n\n#6-C.) JUNCTIONS";
print COMMANDS "\n#Create jobs to check parse mapping results for exon junction hits";
print MAGENTA, "\n\t6-C.) Create jobs to check parse mapping results for exon junction hits", RESET;
$parse_file_name = "parseJunctions.sh";
$parse_file = "$parse_file_dir/$parse_file_name";
open (PARSE, ">$parse_file") || die "\n\nCould not open parse file: $parse_file\n\n";
foreach my $lane_count (sort {$a <=> $b} keys %{$lanes_ref}){
  my $flowcell_name = $lanes_ref->{$lane_count}->{flowcell_name};
  my $lane_number = $lanes_ref->{$lane_count}->{lane_number};
  my $flowcell_name_lane = "$flowcell_name"."_Lane"."$lane_number";
  my $library_id = $lanes_ref->{$lane_count}->{library_id};
  my $total_read_count = $lib_data{$flowcell_name_lane}{total_read_count};
  my $min_bit_score = $lanes_ref->{$lane_count}->{min_bit_score_j};
  my $upper_bit_score = $lanes_ref->{$lane_count}->{upper_bit_score_j};

  print PARSE "source $alexa_seq_conf{SHELL_FILE}; $script_dir/expression/parseBlastResults_JunctionBoundary.pl  --read_records_infile=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/$flowcell_name_lane.txt.gz  --blast_results_dir=$alexa_seq_conf{ANALYSIS_DIR}/blast_results/$library_id/$flowcell_name_lane/exonJunctions/  --reads_analyzed_per_block=$reads_analyzed_per_block  --total_reads_analyzed=$total_read_count  --working_dir=$alexa_seq_conf{ANALYSIS_DIR}/blast_results/$library_id/$flowcell_name_lane/temp/  --summary_outfile=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/Junctions_v$project_conf{ENSEMBL_VERSION}/$flowcell_name_lane"."_Junctions_v$project_conf{ENSEMBL_VERSION}.txt  --min_bit_score=$min_bit_score  --upper_bit_score=$upper_bit_score  --logfile=$alexa_seq_conf{ANALYSIS_DIR}/logs/$library_id/$flowcell_name_lane/parseBlastResults_Junction_LOG.txt  --seq_type=junction  --seq_size=$project_conf{JUNCTION_SEQ_LENGTH}\n";
}
close (PARSE);
if ($cluster_commands){
  print COMMANDS "\n#Execute the following to submit parsing jobs to the cluster";
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}\ncd $parse_file_dir";
  print COMMANDS "\nmqsub  --file $parse_file_name  --name pj  --mkdir  --delay $alexa_seq_conf{VLONG_DELAY}";
}else{
  print COMMANDS "\n#Execute the following or submit the contained jobs to your cluster";
  print COMMANDS "\nbash $parse_file_dir"."$parse_file_name";
}
print COMMANDS "\n\n#Do not proceed to the next step until this step is complete!";
print COMMANDS "\n#Total number of parsing jobs that should be complete = $total_lane_count  (Use the following command to check for job completion):";
print COMMANDS "\n$script_dir/process/checkParsingJobs.pl  --project_name=$project_conf{PROJECT_NAME}  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --ensembl_version=$project_conf{ENSEMBL_VERSION}\n";


#6-D.) BOUNDARIES
print COMMANDS "\n\n#6-D.) BOUNDARIES";
print COMMANDS "\n#Create jobs to check parse mapping results for alternative exon boundary hits";
print MAGENTA, "\n\t6-D.) Create jobs to check parse mapping results for alternative exon boundary hits", RESET;
$parse_file_name = "parseBoundaries.sh";
$parse_file = "$parse_file_dir/$parse_file_name";
open (PARSE, ">$parse_file") || die "\n\nCould not open parse file: $parse_file\n\n";
foreach my $lane_count (sort {$a <=> $b} keys %{$lanes_ref}){
  my $flowcell_name = $lanes_ref->{$lane_count}->{flowcell_name};
  my $lane_number = $lanes_ref->{$lane_count}->{lane_number};
  my $flowcell_name_lane = "$flowcell_name"."_Lane"."$lane_number";
  my $library_id = $lanes_ref->{$lane_count}->{library_id};
  my $total_read_count = $lib_data{$flowcell_name_lane}{total_read_count};
  my $min_bit_score = $lanes_ref->{$lane_count}->{min_bit_score_b};
  my $upper_bit_score = $lanes_ref->{$lane_count}->{upper_bit_score_b};

  print PARSE "source $alexa_seq_conf{SHELL_FILE}; $script_dir/expression/parseBlastResults_JunctionBoundary.pl  --read_records_infile=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/$flowcell_name_lane.txt.gz  --blast_results_dir=$alexa_seq_conf{ANALYSIS_DIR}/blast_results/$library_id/$flowcell_name_lane/exonBoundaries/  --reads_analyzed_per_block=$reads_analyzed_per_block  --total_reads_analyzed=$total_read_count  --working_dir=$alexa_seq_conf{ANALYSIS_DIR}/blast_results/$library_id/$flowcell_name_lane/temp/  --summary_outfile=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/Boundaries_v$project_conf{ENSEMBL_VERSION}/$flowcell_name_lane"."_Boundaries_v$project_conf{ENSEMBL_VERSION}.txt  --min_bit_score=$min_bit_score  --upper_bit_score=$upper_bit_score  --logfile=$alexa_seq_conf{ANALYSIS_DIR}/logs/$library_id/$flowcell_name_lane/parseBlastResults_Boundary_LOG.txt  --seq_type=boundary  --seq_size=$project_conf{BOUNDARY_SEQ_LENGTH}\n";
}
close (PARSE);
if ($cluster_commands){
  print COMMANDS "\n#Execute the following to submit parsing jobs to the cluster";
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}\ncd $parse_file_dir";
  print COMMANDS "\nmqsub  --file $parse_file_name  --name pb  --mkdir  --delay $alexa_seq_conf{VLONG_DELAY}";
}else{
  print COMMANDS "\n#Execute the following or submit the contained jobs to your cluster";
  print COMMANDS "\nbash $parse_file_dir"."$parse_file_name";
}
print COMMANDS "\n\n#Do not proceed to the next step until this step is complete!";
print COMMANDS "\n#Total number of parsing jobs that should be complete = $total_lane_count  (Use the following command to check for job completion):";
print COMMANDS "\n$script_dir/process/checkParsingJobs.pl  --project_name=$project_conf{PROJECT_NAME}  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --ensembl_version=$project_conf{ENSEMBL_VERSION}\n";


#6-E.) INTRONS
print COMMANDS "\n\n#6-E.) INTRONS";
print COMMANDS "\n#Create jobs to check parse mapping results for intron hits";
print MAGENTA, "\n\t6-E.) Create jobs to check parse mapping results for intron hits", RESET;
$parse_file_name = "parseIntrons.sh";
$parse_file = "$parse_file_dir/$parse_file_name";
open (PARSE, ">$parse_file") || die "\n\nCould not open parse file: $parse_file\n\n";
foreach my $lane_count (sort {$a <=> $b} keys %{$lanes_ref}){
  my $flowcell_name = $lanes_ref->{$lane_count}->{flowcell_name};
  my $lane_number = $lanes_ref->{$lane_count}->{lane_number};
  my $flowcell_name_lane = "$flowcell_name"."_Lane"."$lane_number";
  my $library_id = $lanes_ref->{$lane_count}->{library_id};
  my $total_read_count = $lib_data{$flowcell_name_lane}{total_read_count};
  my $min_bit_score = $lanes_ref->{$lane_count}->{min_bit_score};
  my $upper_bit_score = $lanes_ref->{$lane_count}->{upper_bit_score};
  my $min_align_length = $lanes_ref->{$lane_count}->{min_align_length};

  print PARSE "source $alexa_seq_conf{SHELL_FILE}; $script_dir/expression/parseBlastResults_Intron.pl  --intron_db=$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/introns/introns_annotated.txt.gz  --read_records_infile=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/$flowcell_name_lane.txt.gz  --blast_results_dir=$alexa_seq_conf{ANALYSIS_DIR}/blast_results/$library_id/$flowcell_name_lane/introns/  --reads_analyzed_per_block=$reads_analyzed_per_block  --working_dir=$alexa_seq_conf{ANALYSIS_DIR}/blast_results/$library_id/$flowcell_name_lane/temp/  --total_reads_analyzed=$total_read_count  --summary_outfile=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/Introns_v$project_conf{ENSEMBL_VERSION}/$flowcell_name_lane"."_Introns_v$project_conf{ENSEMBL_VERSION}.txt  --min_bit_score=$min_bit_score  --upper_bit_score=$upper_bit_score  --min_align_length=$min_align_length  --logfile=$alexa_seq_conf{ANALYSIS_DIR}/logs/$library_id/$flowcell_name_lane/parseBlastResults_Intron_LOG.txt\n";
}
close (PARSE);
if ($cluster_commands){
  print COMMANDS "\n#Execute the following to submit parsing jobs to the cluster";
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}\ncd $parse_file_dir";
  print COMMANDS "\nmqsub  --file $parse_file_name  --name pi  --mkdir  --delay $alexa_seq_conf{VLONG_DELAY}";
}else{
  print COMMANDS "\n#Execute the following or submit the contained jobs to your cluster";
  print COMMANDS "\nbash $parse_file_dir"."$parse_file_name";
}
print COMMANDS "\n\n#Do not proceed to the next step until this step is complete!";
print COMMANDS "\n#Total number of parsing jobs that should be complete = $total_lane_count  (Use the following command to check for job completion):";
print COMMANDS "\n$script_dir/process/checkParsingJobs.pl  --project_name=$project_conf{PROJECT_NAME}  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --ensembl_version=$project_conf{ENSEMBL_VERSION}\n";


#6-F.) INTERGENICS
print COMMANDS "\n\n#6-F.) INTERGENICS";
print COMMANDS "\n#Create jobs to check parse mapping results for intergenic hits";
print MAGENTA, "\n\t6-F.) Create jobs to check parse mapping results for intergenic hits", RESET;
$parse_file_name = "parseIntergenics.sh";
$parse_file = "$parse_file_dir/$parse_file_name";
open (PARSE, ">$parse_file") || die "\n\nCould not open parse file: $parse_file\n\n";
foreach my $lane_count (sort {$a <=> $b} keys %{$lanes_ref}){
  my $flowcell_name = $lanes_ref->{$lane_count}->{flowcell_name};
  my $lane_number = $lanes_ref->{$lane_count}->{lane_number};
  my $flowcell_name_lane = "$flowcell_name"."_Lane"."$lane_number";
  my $library_id = $lanes_ref->{$lane_count}->{library_id};
  my $total_read_count = $lib_data{$flowcell_name_lane}{total_read_count};
  my $min_bit_score = $lanes_ref->{$lane_count}->{min_bit_score};
  my $upper_bit_score = $lanes_ref->{$lane_count}->{upper_bit_score};
  my $min_align_length = $lanes_ref->{$lane_count}->{min_align_length};

  print PARSE "source $alexa_seq_conf{SHELL_FILE}; $script_dir/expression/parseBlastResults_Intergenic.pl  --intergenic_db=$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/intergenics/intergenics_annotated.txt.gz  --read_records_infile=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/$flowcell_name_lane.txt.gz  --blast_results_dir=$alexa_seq_conf{ANALYSIS_DIR}/blast_results/$library_id/$flowcell_name_lane/intergenics/  --reads_analyzed_per_block=$reads_analyzed_per_block  --working_dir=$alexa_seq_conf{ANALYSIS_DIR}/blast_results/$library_id/$flowcell_name_lane/temp/  --total_reads_analyzed=$total_read_count  --summary_outfile=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/Intergenics_v$project_conf{ENSEMBL_VERSION}/$flowcell_name_lane"."_Intergenics_v$project_conf{ENSEMBL_VERSION}.txt  --min_bit_score=$min_bit_score  --upper_bit_score=$upper_bit_score  --min_align_length=$min_align_length  --logfile=$alexa_seq_conf{ANALYSIS_DIR}/logs/$library_id/$flowcell_name_lane/parseBlastResults_Intergenic_LOG.txt\n";
}
close (PARSE);
if ($cluster_commands){
  print COMMANDS "\n#Execute the following to submit parsing jobs to the cluster";
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}\ncd $parse_file_dir";
  print COMMANDS "\nmqsub  --file $parse_file_name  --name pig  --mkdir  --delay $alexa_seq_conf{VLONG_DELAY}";
}else{
  print COMMANDS "\n#Execute the following or submit the contained jobs to your cluster";
  print COMMANDS "\nbash $parse_file_dir"."$parse_file_name";
}
print COMMANDS "\n\n#Do not proceed to the next step until this step is complete!";
print COMMANDS "\n#Total number of parsing jobs that should be complete = $total_lane_count  (Use the following command to check for job completion):";
print COMMANDS "\n$script_dir/process/checkParsingJobs.pl  --project_name=$project_conf{PROJECT_NAME}  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --ensembl_version=$project_conf{ENSEMBL_VERSION}\n";


#7.) Summarize read mapping and class assignment
#7-A.) Summarize reads assigned to each class for each library - all lanes combined
print COMMANDS "\n\n\n#SUMMARIZE READ MAPPING AND CLASS ASSIGNMENT";
print COMMANDS "\n\n#All individual steps of step 7 can be run concurrently\n";
print COMMANDS "\n#7-A.) Count reads assigned to each category/read class for each library by going through all read record files for that library";
print MAGENTA, "\n\n7-A.) Count reads assigned to each category/read class for each library by going through all read record files for that library", RESET;
my $bash_file11a = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/process/summarizeReadAssignments.sh";
open(BASH, ">$bash_file11a") || die "\n\nCould not open bash output file: $bash_file11a\n\n";
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  print BASH "$script_dir/bash/summarizeReadAssignments.sh $alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/ > $alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/$library_id/Expression_v$project_conf{ENSEMBL_VERSION}/ReadAssignmentStats.txt\n";
}
close(BASH);
print COMMANDS "\n#Execute the following shell script to summarize read assignments";
print COMMANDS "\nbash $bash_file11a";

#7-B.) Summarize mapped reads statistics lane-by-lane and for all lanes combined and determine the library size to be used for library size normalization
print COMMANDS "\n\n\n#7-B.) Summarize mapped reads statistics lane-by-lane and for all lanes combined and determine the library size to be used for library size normalization";
print MAGENTA, "\n\n7-B.) Summarize mapped reads statistics lane-by-lane and for all lanes combined and determine the library size to be used for library size normalization", RESET;
my $bash_file11b = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/process/calculateLibrarySize.sh";
open(BASH, ">$bash_file11b") || die "\n\nCould not open bash output file: $bash_file11b\n\n";
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  print BASH "$script_dir/expression/calculateLibrarySize.pl  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --library_id=$library_id\n";
}
close(BASH);
print COMMANDS "\n#Execute the following shell script to summarize read assignments";
print COMMANDS "\nbash $bash_file11b";

#7-C.) Summarize the relative position of mapped reads within transcript of increasing size: position bias
print COMMANDS "\n\n\n#7-C.) Summarize the relative position of mapped reads within transcript of increasing size: position bias";
print MAGENTA, "\n\n7-C.) Summarize the relative position of mapped reads within transcript of increasing size: position bias", RESET;
my $bash_file11c = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/process/positionBias.sh";
open(BASH, ">$bash_file11c") || die "\n\nCould not open bash output file: $bash_file11c\n\n";
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  print BASH "$script_dir/stats/positionBiasVsGeneSize.pl  --data_dir=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/ENST_v$project_conf{ENSEMBL_VERSION}/  --working_dir=$alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/$library_id/temp/  --results_dir=$alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/$library_id/ENST_v$project_conf{ENSEMBL_VERSION}/\n";
}
close(BASH);
print COMMANDS "\n#Execute the following shell script to summarize read assignments";
print COMMANDS "\nbash $bash_file11c";


#7-D.) Determine the tag redundancy (aka library complexity) of each library
print COMMANDS "\n\n\n#7-D.) Determine the tag redundancy (aka library complexity) of each library";
print MAGENTA, "\n\n7-D.) Determine the tag redundancy (aka library complexity) of each library", RESET;
my $bash_file5b = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/stats/estimateComplexity.sh";
open(BASH, ">$bash_file5b") || die "\n\nCould not open bash output file: $bash_file5b\n\n";
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  my $sample_size=1000000;
  if ($demo){$sample_size=100000;}
  print BASH "$script_dir/stats/estimateComplexity_SEQ.pl  --read_records_dir=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/  --col1=6  --col2=7  --temp_dir=$alexa_seq_conf{ANALYSIS_DIR}/temp/complexity/  --paired=yes  --outfile=$alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/$library_id/LibraryQuality/LibraryComplexity_SEQ.txt  --iterations=3  --header_present=yes  --block_size=$sample_size\n";

  print BASH "$script_dir/stats/estimateComplexity_MAP.pl  --read_records_dir=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --temp_dir=$alexa_seq_conf{ANALYSIS_DIR}/temp/complexity/  --outfile=$alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/$library_id/LibraryQuality/LibraryComplexity_MAP.txt  --iterations=3  --block_size=$sample_size\n";
}
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  print BASH "$script_dir/stats/summarizeComplexity.pl  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --library_id=$library_id  --script_dir=$script_dir/  --complexity_type=SEQ\n";
  print BASH "$script_dir/stats/summarizeComplexity.pl  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --library_id=$library_id  --script_dir=$script_dir/  --complexity_type=MAP\n";
}
close(BASH);
print COMMANDS "\n\n#Execute the following shell script to determine the tag redundancy of each library";
print COMMANDS "\nbash $bash_file5b";


#7-E.) Partition mapping results files into genome region partitions
print COMMANDS "\n\n\n#7-E.) Partition mapping results files into genome region partitions";
print MAGENTA, "\n\n7-E.) Partition mapping results files into genome region partitions", RESET;
my @class_list = qw (ENST NOVEL_JUNCTION NOVEL_BOUNDARY INTRON INTERGENIC);
my $bash_file11d = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/expression/partitionMapFiles.sh";
open(BASH, ">$bash_file11d") || die "\n\nCould not open bash output file: $bash_file11d\n\n";
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  foreach my $class (@class_list){
    my $seq_size = 1;
    if ($class eq "NOVEL_JUNCTION"){
      $seq_size=$project_conf{JUNCTION_SEQ_LENGTH};
    }
    if ($class eq "NOVEL_BOUNDARY"){
      $seq_size=$project_conf{BOUNDARY_SEQ_LENGTH};
    }
    print BASH "source $alexa_seq_conf{SHELL_FILE}; $script_dir/expression/partitionMapFiles.pl  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --annotation_dir=$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/  --seq_size=$seq_size  --library_id=$library_id  --genome_regions_file=$region_file  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --read_class=$class\n";
  }
}
close(BASH);
if($cluster_commands){
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}\ncd $expression_file_dir";
  print COMMANDS "\nmqsub  --file $bash_file11d  --name pmf  --mkdir  --delay $alexa_seq_conf{VLONG_DELAY}";
}else{
  print COMMANDS "\n#Execute the following or submit the contained jobs to your cluster";
  print COMMANDS "\nbash $bash_file11d";
}

print COMMANDS "\n\n#Check for successful completion of these jobs with the following command:";
print COMMANDS "\n$script_dir/process/checkMapPartitionJobs.pl  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --project_name=$project_conf{PROJECT_NAME}  --partition_file=$region_file";


#8.) EXPRESSION CALCULATIONS
print MAGENTA, "\n\n8.) EXPRESSION CALCULATIONS", RESET;
print COMMANDS "\n\n\n#8.) EXPRESSION CALCULATIONS";
print COMMANDS "\n\n#Note: Once steps 8-A and 8-B are complete, steps 8-C to 8-E can be run concurrently";

#8-A.) Generate batch jobs for generating expression values for all feature types:
print COMMANDS "\n\n#8-A.) Generate batch jobs for generating expression values for all feature types:";
print MAGENTA, "\n\t8-A.) Generate batch jobs for generating expression values for all feature types:", RESET;
my $bash_file12 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/expression/submitExpressionJobs.sh";
my $expression_cmd;
if ($cluster_commands){
  $expression_cmd = "$script_dir/process/createExpressionJobs.pl  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --ucsc_build=$project_conf{UCSC_BUILD}  --lib_names_file=$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$project_conf{PROJECT_NAME}"."_Lib_Names.txt  --script_dir=$script_dir/  --batch_dir=$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/expression/  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --annotation_dir=$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/  --regions_file_version=$project_conf{REGIONS_FILE_VERSION}  --ucsc_dir=$alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME}/ucsc/  --web_path=$alexa_seq_conf{WEB_URL}/$project_conf{PROJECT_NAME}/ucsc  --junction_size=$project_conf{JUNCTION_SEQ_LENGTH}  --boundary_size=$project_conf{BOUNDARY_SEQ_LENGTH}  --cluster_commands=1 > $bash_file12";
}else{
  $expression_cmd = "$script_dir/process/createExpressionJobs.pl  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --ucsc_build=$project_conf{UCSC_BUILD}  --lib_names_file=$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$project_conf{PROJECT_NAME}"."_Lib_Names.txt  --script_dir=$script_dir/  --batch_dir=$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/expression/  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --annotation_dir=$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/  --regions_file_version=$project_conf{REGIONS_FILE_VERSION}  --ucsc_dir=$alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME}/ucsc/  --web_path=$alexa_seq_conf{WEB_URL}/$project_conf{PROJECT_NAME}/ucsc  --junction_size=$project_conf{JUNCTION_SEQ_LENGTH}  --boundary_size=$project_conf{BOUNDARY_SEQ_LENGTH} > $bash_file12";
}
print BLUE "\n\n$expression_cmd\n", RESET;
system($expression_cmd);

if($cluster_commands){
  print COMMANDS "\n#Execute the following to submit expression calculation jobs to the cluster";
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}\ncd $alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/expression/";
  print COMMANDS "\nbash $bash_file12";
  print COMMANDS "\n\n#Do not proceed to the next step until all jobs have completed!";
}else{
  print COMMANDS "\n#Execute the following or submit the contained jobs to your cluster";
  print COMMANDS "\nbash $bash_file12";
}

#8-B.) Check all resulting expression jobs, join results files into final expression files, and clean up temp files:
print COMMANDS "\n\n#8-B.) Check all resulting expression jobs, join results files into final expression files, and clean up temp files:";
print MAGENTA, "\n\t8-B.) Check all resulting expression jobs, join results files into final expression files, and clean up temp files:", RESET;

my $lib_string = '';
my $lib_count = keys %{$libraries_ref};
my $lc = 0;
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  $lc++;
  if ($lc == $lib_count){
    $lib_string .= "$library_id";
  }else{
    $lib_string .= "$library_id,";
  }
}
print COMMANDS "\n$script_dir/process/checkExpressionJobs.pl  --source_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --batch_dir=$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/expression/  --libraries='$lib_string'  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --annotation_dir=$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/  --junction_seq_size=$project_conf{JUNCTION_SEQ_LENGTH}  --boundary_seq_size=$project_conf{BOUNDARY_SEQ_LENGTH}  --regions_file_version=$project_conf{REGIONS_FILE_VERSION}";

#8-C.) Calculate transcript specific expression values
print COMMANDS "\n\n#8-C.) Calculate transcript specific expression values (derived from exon and junction expression values)";
print MAGENTA, "\n\t8-C.) Calculate transcript specific expression values (derived from exon and junction expression values):", RESET;
my $bash_file13 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/expression/generateExpressionValues_Transcript.sh";
open(BASH, ">$bash_file13") || die "\n\nCould not open bash output file: $bash_file13\n\n";
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  print BASH "$script_dir/expression/generateExpressionValues_Transcript.pl  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --annotation_dir=$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/  --exon_region_file=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/ENST_v$project_conf{ENSEMBL_VERSION}/Summary/$library_id"."_ExonRegionExpression_v$project_conf{ENSEMBL_VERSION}.txt  --exon_junction_file=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/Junctions_v$project_conf{ENSEMBL_VERSION}/Summary/$library_id"."_JunctionExpression_v$project_conf{ENSEMBL_VERSION}.txt  --results_dir=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/Transcripts_v$project_conf{ENSEMBL_VERSION}/Summary/  --library_name=$library_id  --cutoffs_file=$alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/$library_id/Expression_v$project_conf{ENSEMBL_VERSION}/$library_id"."_NORM1_average_coverage_cutoffs.txt\n";
}
close(BASH);
print COMMANDS "\n#Execute the following shell script to calculate transcript specific expression values";
print COMMANDS "\nbash $bash_file13";

#8-D.) Compare exon region and junction/boundary expression values and determine the median difference (then use this to correct the expression values calculated above)
print COMMANDS "\n\n#8-D.) Compare exon region and junction/boundary expression values and determine the median difference (then use this to correct the expression values calculated above)";
print MAGENTA, "\n\t8-D.) Compare exon region and junction/boundary expression values and determine the median difference (then use this to correct the expression values calculated above):", RESET;
print COMMANDS "\n#NOTE: In the following analysis, many genes have canonical junctions for the calculation";
print COMMANDS "\n#However, for the exon boundaries, those with EnsEMBL support correspond to known alternative trancripts (not really analagous to the situation with junctions)";
print COMMANDS "\n#For this reason it is the exon junction derived correction factor will be used!";
print COMMANDS "\n#After determining the correction factor it will be added to the lib_names file:";
print COMMANDS "\n#$project_lib_names_file";
my $bash_file14 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/expression/compareExonAndJunctionExpression.sh";
open(BASH, ">$bash_file14") || die "\n\nCould not open bash output file: $bash_file14\n\n";
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  print BASH "$script_dir/process/compareExonAndJunctionExpression.pl  --exon_region_file=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/ENST_v$project_conf{ENSEMBL_VERSION}/Summary/$library_id"."_ExonRegionExpression_v$project_conf{ENSEMBL_VERSION}.txt  --exon_junction_file=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/Junctions_v$project_conf{ENSEMBL_VERSION}/Summary/$library_id"."_JunctionExpression_v$project_conf{ENSEMBL_VERSION}.txt  --exon_boundary_file=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/Boundaries_v$project_conf{ENSEMBL_VERSION}/Summary/$library_id"."_BoundaryExpression_v$project_conf{ENSEMBL_VERSION}.txt  --er_expression_cutoff=4.2  --library=$library_id  --library_name_file=$project_lib_names_file\n";
}
close(BASH);
print COMMANDS "\n#Execute the following shell script to calculate determine the exon junction and boundary correction factors";
print COMMANDS "\nbash $bash_file14";

#8-E.) Divide expression results for both JUNCTIONS and BOUNDARIES into those with or without EnsEMBL support (KNOWN or NOVEL)
print COMMANDS "\n\n#8-E.) Divide expression results for both JUNCTIONS and BOUNDARIES into those with or without EnsEMBL support (KNOWN or NOVEL)";
print MAGENTA, "\n\t8-E.) Divide expression results for both JUNCTIONS and BOUNDARIES into those with or without EnsEMBL support (KNOWN or NOVEL)", RESET;
my $bash_file15 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/expression/createNovelKnownFiles.sh";
open(BASH, ">$bash_file15") || die "\n\nCould not open bash output file: $bash_file15\n\n";
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  print BASH "\n$script_dir/bash/createNovelKnownFiles.sh $alexa_seq_conf{ANALYSIS_DIR}/read_records $library_id $project_conf{ENSEMBL_VERSION}";
}
close(BASH);
print COMMANDS "\n#Execute the following shell script to create known/novel versions of exon junction and boundary expression results";
print COMMANDS "\nbash $bash_file15";


#9.) SUMMARIZE EXPRESSION STATISTICS
print MAGENTA, "\n\n9.) SUMMARIZE EXPRESSION STATISTICS", RESET;
print COMMANDS "\n\n\n#9.) SUMMARIZE EXPRESSION STATISTICS";
print COMMANDS "\n\n#Note: All individual steps of step 9 can be run concurrently";

my @image_types_lc = qw (svg jpeg);
my @image_types_uc = qw (SVG JPEG);

#9-A.) Summarize expression values and generate expression cutoff values and other statistics for each library
print COMMANDS "\n\n#9-A.) Summarize expression values and generate expression cutoff values and other statistics for each library";
print MAGENTA, "\n\t9-A.) Summarize expression values and generate expression cutoff values and other statistics for each library:", RESET;
my $bash_file16 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/stats/summarizeExpressionValues.sh";
open(BASH, ">$bash_file16") || die "\n\nCould not open bash output file: $bash_file16\n\n";
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  foreach my $image_type (@image_types_lc){
    print BASH "$script_dir/stats/summarizeExpressionValues.pl  --library_name=$library_id  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --library_dir=$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/  --working_dir=$alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/$library_id/temp/  --results_dir=$alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/$library_id/Expression_v$project_conf{ENSEMBL_VERSION}/ --image_type=$image_type\n";
  }
}
close(BASH);
print COMMANDS "\n#Execute the following shell script to summarize expression values and generate cutoffs";
print COMMANDS "\nbash $bash_file16";

#9-B.) Summarize gene coverage statistics for each library
print COMMANDS "\n\n#9-B.) Summarize gene coverage statistics for each library";
print MAGENTA, "\n\t9-B.) Summarize gene coverage statistics for each library:", RESET;
my $bash_file17 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/stats/summarizeGeneExonCoverage.sh";
open(BASH, ">$bash_file17") || die "\n\nCould not open bash output file: $bash_file17\n\n";
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  my $library_name = $libraries_ref->{$library_id}->{library_name}; 
  foreach my $image_type (@image_types_uc){
    print BASH "$script_dir/R_bin/summarizeGeneExonCoverage_"."$image_type.R $alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/ENST_v$project_conf{ENSEMBL_VERSION}/Summary/$library_id"."_GeneExpression_v$project_conf{ENSEMBL_VERSION}.txt  $alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/ENST_v$project_conf{ENSEMBL_VERSION}/Summary/$library_id"."_ExonRegionExpression_v$project_conf{ENSEMBL_VERSION}.txt  $alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/$library_id/Generic/  $alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/$library_id/ENST_v$project_conf{ENSEMBL_VERSION}/  $project_conf{ENSEMBL_VERSION}  $library_name\n";
  }
}
close(BASH);
print COMMANDS "\n#Execute the following shell script to generate gene coverage statistics";
print COMMANDS "\nbash $bash_file17";


#9-C.) Summarize the number of features of each type are expressed above background
print COMMANDS "\n\n#9-C.) Summarize the number of features of each type are expressed above background";
print MAGENTA, "\n\t9-C.) Summarize the number of features of each type are expressed above background:", RESET;
my $bash_file18 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/stats/summarizeExpressedEvents.sh";
open(BASH, ">$bash_file18") || die "\n\nCould not open bash output file: $bash_file18\n\n";
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  print BASH "$script_dir/bash/summarizeExpressedEvents.sh $alexa_seq_conf{ANALYSIS_DIR} $project_conf{ENSEMBL_VERSION} $library_id\n";
}
close(BASH);
print COMMANDS "\n#Execute the following shell script to summarize features expressed above background";
print COMMANDS "\nbash $bash_file18";


#9-D.) Summarize the library expression skew
print COMMANDS "\n\n#9-D.) Summarize the library expression skew";
print MAGENTA, "\n\t9-D.) Summarize the library expression skew:", RESET;
my $bash_file18a = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/stats/summarizeLibrarySkew.sh";
open(BASH, ">$bash_file18a") || die "\n\nCould not open bash output file: $bash_file18\n\n";
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  print BASH "$script_dir/stats/summarizeLibrarySkew.pl  --library_name=$library_id  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/\n";
}
close(BASH);
print COMMANDS "\n#Execute the following shell script to summarize library expression skew";
print COMMANDS "\nbash $bash_file18a";


#9-E.)Update expression calculate jobs
print COMMANDS "\n\n\n#9-E.) Update expression calculation jobs";
print MAGENTA, "\n\n9-E.) Update expression calculation jobs:", RESET;
print COMMANDS "\n#IMPORTANT: To do this, run the command below";
my $bash_file18b = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/process/updateExpressionJobs.sh";
open(BASH, ">$bash_file18b") || die "\n\nCould not open bash output file: $bash_file18b\n\n";
print BASH "$expression_cmd";
close(BASH);
print COMMANDS "\nbash $bash_file18b";


#10.) Re-generate all expression values now that the cutoffs files and junction/boundary correction factors have been obtained
print COMMANDS "\n\n\n#10.) Re-generate all expression values now that the cutoffs files and junction/boundary correction factors have been obtained";
print MAGENTA, "\n\n10.) Re-generate all expression values now that the cutoffs files and junction/boundary correction factors have been obtained:", RESET;
print COMMANDS "\n#IMPORTANT: To do this, return to step #8 and repeat all steps back to this one before continuing";

#11.) Create expression matrix files for each feature. 
print COMMANDS "\n\n#11.) Create expression matrix files for each feature.";
print MAGENTA, "\n\n11.) Create expression matrix files for each feature.", RESET;

print COMMANDS "\n#Execute the following to create expression matrix files";
print COMMANDS "\n$script_dir/website/createMatrixFiles.pl  --project_name=$project_conf{PROJECT_NAME}  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --web_dir=$alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME}/data/"; 
print COMMANDS "\n\n#Note: you can start step 11 and then move on to the next step before it finishes";

#12.) Create merged UCSC track files combining all the track types created above
print COMMANDS "\n\n\n#12.) Create merged UCSC track files combining all the track types created above";
print MAGENTA, "\n\n12.) Create merged UCSC track files combining all the track types created above:", RESET;

#12-A.) First merge tracks for the various data types (exons, introns, etc.) into one file for each chromosome region
print COMMANDS "\n\n#12-A.) First merge tracks for the various data types (exons, introns, etc.) into one file for each chromosome region";
print MAGENTA, "\n\t12A.) First merge tracks for the various data types (exons, introns, etc.) into one file for each chromosome region:", RESET;
my $bash_file19 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/website/mergeUcscTrackFiles.sh";
open(BASH, ">$bash_file19") || die "\n\nCould not open bash output file: $bash_file19\n\n";
my $priority=80;
foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
  $priority++;
  my $group = $libraries_ref->{$library_id}->{group};
  my $library_name = $libraries_ref->{$library_id}->{library_name}; 
  print BASH "$script_dir/website/mergeUcscTrackFiles.pl  --input_dir=$alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME}/ucsc/$library_id/  --output_dir=$alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME}/ucsc/$library_id/combined/  --library_name=\"$library_name\"  --color_set=$group  --priority=$priority\n";
}
close(BASH);
if($cluster_commands){
  print COMMANDS "\n#Execute the following to submit wig merging jobs to the cluster";
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}\ncd $alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/website/";
  print COMMANDS "\nmqsub  --file $bash_file19  --qsub '-l mem_total=7G'  --name mutf  --mkdir  --delay $alexa_seq_conf{VLONG_DELAY}";
}else{
  print COMMANDS "\n#Execute the following shell script to merge UCSC tracks";
  print COMMANDS "\nbash $bash_file19";
}


#12-B.) Use the normalized expression wig track to generate a DE wig track.   Stabilize variance by adding 16 before taking log (DE = log2(Post+16) - log2(Pre+16))
print COMMANDS "\n\n#12-B.) Use the normalized expression wig track to generate Differential Expression UCSC wig tracks.   Stabilize variance by adding 16 before taking log (DE = log2(Post+16) - log2(Pre+16))";
print MAGENTA, "\n\t12-B.) Use the normalized expression wig track to generate Differential Expression UCSC wig tracks.", RESET;
print COMMANDS "\n#Desired comparisons  (the 'earlier' condition is listed first)";
print COMMANDS "\n#In the following calculations DE is (LibraryA - LibraryB)";
print COMMANDS "\n\n#Comparisons defined:";
my $bash_file20 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/website/mergeUcscTrackFiles_DE.sh";
open(BASH, ">$bash_file20") || die "\n\nCould not open bash output file: $bash_file20\n\n";
foreach my $comparison (sort {$a <=> $b} keys %{$comparisons_ref}){
  my $libA_id = $comparisons_ref->{$comparison}->{libraryA_id};
  my $libA_name = $comparisons_ref->{$comparison}->{libraryA_name};
  my $libB_id = $comparisons_ref->{$comparison}->{libraryB_id};
  my $libB_name = $comparisons_ref->{$comparison}->{libraryB_name};
  my $comparison_id = $comparisons_ref->{$comparison}->{comparison_id};
  my $comparison_name = $comparison_id;
  $comparison_name =~ s/_/ /g;
  my $comp_lib_id = "$libA_id"."_vs_"."$libB_id";
  print COMMANDS "\n#Library A ($libA_id aka $libA_name) - Library B ($libB_id aka $libB_name).  Comparison ($comp_lib_id aka $comparison_name)";
  print BASH "\n$script_dir/website/mergeUcscTrackFiles_DE.pl  --input_dir_A=$alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME}/ucsc/$libA_id/combined/  --input_dir_B=$alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME}/ucsc/$libB_id/combined/  --output_dir_A=$alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME}/ucsc/DE/$libA_id/$comp_lib_id/  --output_dir_B=$alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME}/ucsc/DE/$libB_id/$comp_lib_id/  --comparison_name=\"$comparison_name\"  --short_name=\"$comparison_id\"  --view_limit=\"-3.5:3.5\"";
}
close(BASH);
if($cluster_commands){
  print COMMANDS "\n#Execute the following to submit DE wig creation jobs to the cluster";
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}\ncd $alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/website/";
  print COMMANDS "\nmqsub  --file $bash_file20  --qsub '-l mem_total=7G'  --name mutfd  --mkdir  --delay $alexa_seq_conf{VLONG_DELAY}";
}else{
  print COMMANDS "\n#Execute the following shell script to create DE UCSC tracks";
  print COMMANDS "\nbash $bash_file20";
}


#13 and #14 Perform DE and AE pairwise and groupwise comparisons if any found
my $cmd;

#13A.) Calculate Differential Expression (DE) of Genes, Exons, Junctions, Boundaries, etc. for pairwise library comparisons
print COMMANDS "\n\n\n#Note: The following steps 13 and 14 can be run concurrently";
print COMMANDS "\n\n#13A.) Calculate Differential Expression (DE) of Genes, Exons, Junctions, Boundaries, etc. for pairwise library comparisons";
print MAGENTA, "\n\n13A.) Calculate Differential Expression (DE) of Genes, Exons, Junctions, Boundaries, etc. for pairwise library comparisons:", RESET;
if ($comparison_count>0){
  my $bash_file21 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/stats/calculateDifferentialExpression.sh";
  $cmd = "$script_dir/process/createComparisonJobs.pl  --comparison_type=DE  --comparison_file=$project_lib_comparisons_file  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --script_dir=$script_dir/stats/  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --project_name=$project_conf{PROJECT_NAME}  --outfile=$bash_file21  --silent=1";
  print BLUE "\n\n$cmd", RESET;
  system($cmd);
  print COMMANDS "\n#Execute the following shell script to calculate pairwise DE values";
  print COMMANDS "\nbash $bash_file21";
}else{
  print COMMANDS "\n#No pairwise library comparisons found";
  print MAGENTA, "\nNo pairwise library comparisons found", RESET;
}

#13B.) Calculate Differential Expression (DE) of Genes, Exons, Junctions, Boundaries, etc. for group comparisons
print COMMANDS "\n\n\n#13B.) Calculate Differential Expression (DE) of Genes, Exons, Junctions, Boundaries, etc. for group comparisons";
print MAGENTA, "\n\n13B.) Calculate Differential Expression (DE) of Genes, Exons, Junctions, Boundaries, etc. for group comparisons:", RESET;
if ($group_comparison_count>0){
  my $bash_file21b = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/stats/calculateGroupDifferentialExpression.sh";
  $cmd = "$script_dir/process/createGroupComparisonJobs.pl  --comparison_type=DE  --group_comparison_file=$project_group_comparisons_file  --group_members_file=$project_lib_groups_file  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --script_dir=$script_dir/stats/  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --project_name=$project_conf{PROJECT_NAME}  --outfile=$bash_file21b  --silent=1";
  print BLUE "\n\n$cmd", RESET;
  system($cmd);
  print COMMANDS "\n#Execute the following shell script to calculate groupwise DE values";
  print COMMANDS "\nbash $bash_file21b";
}else{
  print COMMANDS "\n#No groupwise library comparisons found";
  print MAGENTA, "\nNo groupwise library comparisons found", RESET;
}

#14A.) Calculate Alternative Expression (AE) of Genes, Exons, Junctions, Boundaries, etc. for pairwise library comparisons
print COMMANDS "\n\n\n#14A.) Calculate Alternative Expression (AE) of Genes, Exons, Junctions, Boundaries, etc. for pairwise library comparisons";
print MAGENTA, "\n\n14A.) Calculate Alternative Expression (AE) of Genes, Exons, Junctions, Boundaries, etc. for pairwise library comparisons:", RESET;
if ($comparison_count>0){
  my $bash_file22 = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/stats/calculateAlternativeExpression.sh";
  $cmd = "$script_dir/process/createComparisonJobs.pl  --comparison_type=SI  --comparison_file=$project_lib_comparisons_file  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --script_dir=$script_dir/stats/  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --project_name=$project_conf{PROJECT_NAME}  --outfile=$bash_file22  --silent=1";
  print BLUE "\n\n$cmd", RESET;
  system($cmd);
  print COMMANDS "\n#Execute the following shell script to calculate pairwise AE values";
  print COMMANDS "\nbash $bash_file22";
}else{
  print COMMANDS "\n#No pairwise library comparisons found";
  print MAGENTA, "\nNo pairwise library comparisons found", RESET;
}

#14B.) Calculate Alternative Expression (AE) of Genes, Exons, Junctions, Boundaries, etc. for group comparisons
print COMMANDS "\n\n\n#14B.) Calculate Alternative Expression (AE) of Genes, Exons, Junctions, Boundaries, etc. for group comparisons";
print MAGENTA, "\n\n14B.) Calculate Alternative Expression (AE) of Genes, Exons, Junctions, Boundaries, etc. for group comparisons:", RESET;
if ($group_comparison_count>0){
  my $bash_file22d = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/stats/calculateGroupAlternativeExpression.sh";
  $cmd = "$script_dir/process/createGroupComparisonJobs.pl  --comparison_type=SI  --group_comparison_file=$project_group_comparisons_file  --group_members_file=$project_lib_groups_file  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --script_dir=$script_dir/stats/  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --project_name=$project_conf{PROJECT_NAME}  --outfile=$bash_file22d  --silent=1";
  print BLUE "\n\n$cmd", RESET;
  system($cmd);
  print COMMANDS "\n#Execute the following shell script to calculate groupwise AE values";
  print COMMANDS "\nbash $bash_file22d";
}else{
  print COMMANDS "\n#No groupwise library comparisons found";
  print MAGENTA, "\nNo groupwise library comparisons found", RESET;
}


#15.) Populate ALEXA-Seq data viewer to display summary pages for Expression, DE and AE as well as gene-by-gene summaries and visualizations
print COMMANDS "\n\n\n#15.) Populate ALEXA-Seq data viewer to display summary pages for Expression, DE and SI as well as gene-by-gene summaries and visualizations";
print COMMANDS "\n\n#Note: Steps 15 A-G can be run concurrently";
print MAGENTA, "\n\n15.) Populate ALEXA-Seq data viewer to display summary pages for Expression, DE and SI as well as gene-by-gene summaries and visualizations", RESET;

my $project_library_paths_file = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/$project_conf{PROJECT_NAME}_Lib_Paths.txt";

#15-A.) Creating project library paths file ($project_library_paths_file) (unless it already exists)
print COMMANDS "\n\n#15-A.) Creating project library paths file ($project_library_paths_file) (unless it already exists)";
print MAGENTA, "\n\t15-A.) Creating project library paths file ($project_library_paths_file) (unless it already exists)", RESET;

#Make library paths file for this project:
if ((! -e $project_library_paths_file) || ($clean == 1)){
  open(PLP, ">$project_library_paths_file") || die "\n\nCould not open Project Library Paths file: $project_library_paths_file\n\n";
  print PLP "Library SampleName DataType DataPath StatsPath\n\n";
  print PLP "#Root directories for expression data and stats\n";
  foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
    my $library_name = $libraries_ref->{$library_id}->{library_name}; 
    my $data_path = "$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/";
    my $stats_path = "$alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/$library_id/Expression_v$project_conf{ENSEMBL_VERSION}/";
    print PLP "$library_id $library_name Expression $data_path $stats_path\n";
  }
  print PLP "\n#Root directories for differential expression data and stats\n";
  foreach my $comparison (sort {$a <=> $b} keys %{$comparisons_ref}){
    my $libA_id = $comparisons_ref->{$comparison}->{libraryA_id};
    my $libA_name = $comparisons_ref->{$comparison}->{libraryA_name};
    my $libB_id = $comparisons_ref->{$comparison}->{libraryB_id};
    my $libB_name = $comparisons_ref->{$comparison}->{libraryB_name};
    my $comparison_id = $comparisons_ref->{$comparison}->{comparison_id};
    my $comparison_name = $comparison_id;
    $comparison_name =~ s/_/ /g;
    my $comp1 = "$libA_id"."_vs_"."$libB_id";
    my $comp2 = "$libA_name"."_vs_"."$libB_name";
    my $data_path = "$alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/DE/$project_conf{PROJECT_NAME}/";
    my $stats_path = $data_path;
    print PLP "$comp1 $comp2 DifferentialExpression $data_path $stats_path\n";
  }
  print PLP "\n#Root directories for groupwise differential expression data and stats\n";
  foreach my $group_comparison (sort {$a <=> $b} keys %{$group_comparisons_ref}){
    my $group_comparison_id = $group_comparisons_ref->{$group_comparison}->{group_comparison_id};
    my $data_path = "$alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/DE/$project_conf{PROJECT_NAME}/";
    my $stats_path = $data_path;
    print PLP "$group_comparison_id $group_comparison_id GroupDifferentialExpression $data_path $stats_path\n";
  }
  print PLP "\n#Root directories for differential splicing data and stats\n";
  foreach my $comparison (sort {$a <=> $b} keys %{$comparisons_ref}){
    my $libA_id = $comparisons_ref->{$comparison}->{libraryA_id};
    my $libA_name = $comparisons_ref->{$comparison}->{libraryA_name};
    my $libB_id = $comparisons_ref->{$comparison}->{libraryB_id};
    my $libB_name = $comparisons_ref->{$comparison}->{libraryB_name};
    my $comparison_id = $comparisons_ref->{$comparison}->{comparison_id};
    my $comparison_name = $comparison_id;
    $comparison_name =~ s/_/ /g;
    my $comp1 = "$libA_id"."_vs_"."$libB_id";
    my $comp2 = "$libA_name"."_vs_"."$libB_name";
    my $data_path = "$alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/SI/$project_conf{PROJECT_NAME}/";
    my $stats_path = $data_path;
    print PLP "$comp1 $comp2 DifferentialSplicing $data_path $stats_path\n";
  }
  print PLP "\n#Root directories for groupwise differential splicing data and stats\n";
  foreach my $group_comparison (sort {$a <=> $b} keys %{$group_comparisons_ref}){
    my $group_comparison_id = $group_comparisons_ref->{$group_comparison}->{group_comparison_id};
    my $data_path = "$alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/SI/$project_conf{PROJECT_NAME}/";
    my $stats_path = $data_path;
    print PLP "$group_comparison_id $group_comparison_id GroupDifferentialSplicing $data_path $stats_path\n";
  }
  print PLP "\n#Root directories for quality data and stats\n";
  foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
    my $library_name = $libraries_ref->{$library_id}->{library_name}; 
    my $data_path = "$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/Summary/";
    my $stats_path = "$alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/$library_id/LibraryQuality/";
    print PLP "$library_id $library_name LibraryQuality $data_path $stats_path\n";
  }
  print PLP "\n#Root directories for transcript mapping data and stats\n";
  foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
    my $library_name = $libraries_ref->{$library_id}->{library_name}; 
    my $data_path = "$alexa_seq_conf{ANALYSIS_DIR}/read_records/$library_id/ENST_v$project_conf{ENSEMBL_VERSION}/";
    my $stats_path = "$alexa_seq_conf{ANALYSIS_DIR}/figures_and_stats/$library_id/ENST_v$project_conf{ENSEMBL_VERSION}/";
    print PLP "$library_id $library_name TranscriptMapping $data_path $stats_path\n";
  }
  print PLP "\n#Root directories for analysis commands and basic source data information\n";
  foreach my $library_id (sort {$libraries_ref->{$a}->{library_count} <=> $libraries_ref->{$b}->{library_count}} keys %{$libraries_ref}){
    my $library_name = $libraries_ref->{$library_id}->{library_name}; 
    my $data_path = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/";
    print PLP "$library_id $library_name AnalysisCommands $data_path $data_path\n";
  }

  close(PLP);
}

#15-B.) Package Expression, DE, AE and matrix summary files for download
print COMMANDS "\n\n#15-B.) Package Expression, DE, AE and matrix summary files for download";
print MAGENTA, "\n\t15-B.) Package Expression, DE, AE and matrix summary files for download", RESET;

print COMMANDS "\n#Execute the following to package all results files";
print COMMANDS "\n$script_dir/website/packageResultsFiles.pl  --project_name=$project_conf{PROJECT_NAME}  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --web_dir=$alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME}/data/"; 

#15-C.) Create grand summaries of Expression, Differential Expression and Alternative Expression results
print COMMANDS "\n\n#15-C.) Create grand summaries of Expression, Differential Expression and Alternative Expression results";
print MAGENTA, "\n\t15-C.) Create grand summaries of Expression, Differential Expression and Alternative Expression results", RESET;

print COMMANDS "\n#Execute the following to generate grand summary files";
print COMMANDS "\n$script_dir/$project_conf{WEBSITE_VERSION}/generateSummaryHtml.pl  --junction_seq_size=$project_conf{JUNCTION_SEQ_LENGTH}  --boundary_seq_size=$project_conf{BOUNDARY_SEQ_LENGTH}  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --annotation_dir=$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/  --paths_file=$project_library_paths_file  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --ucsc_build=$project_conf{UCSC_BUILD}  --web_root=$alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME}/  --partition_file=$region_file  --track_url=$alexa_seq_conf{WEB_URL}/$project_conf{PROJECT_NAME}/ucsc/  --search_page_url=$alexa_seq_conf{XAPIAN_URL}  --project_title=\"$project_conf{PROJECT_DESCRIPTION}\"  --project_name=$project_conf{PROJECT_NAME}  --project_config_file=$project_config_file  --alexa_home_path=$alexa_seq_conf{ALEXA_HOME_URL}  --alexa_seq_path=$alexa_seq_conf{ALEXA_SEQ_URL}  --google_analytics_id=$alexa_seq_conf{GOOGLE_ANALYTICS_ID}"; 

#15-D.) Create library summaries (statistics and graphics that describe the characteristics of each sequence library being compared)
print COMMANDS "\n\n#15-D.) Create library summaries (statistics and graphics that describe the characteristics of each sequence library being compared)";
print MAGENTA, "\n\t15-D.) Create library summaries (statistics and graphics that describe the characteristics of each sequence library being compared)", RESET;
print COMMANDS "\n#Execute the following to create a summary of each library";
print COMMANDS "\n$script_dir/website/generateLibrarySummaries.pl  --project_name=$project_conf{PROJECT_NAME}  --paths_file=$project_library_paths_file  --search_page_url=$alexa_seq_conf{XAPIAN_URL}  --target_web_dir=$alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME}/  --alexa_home_path=$alexa_seq_conf{ALEXA_HOME_URL}  --alexa_seq_path=$alexa_seq_conf{ALEXA_SEQ_URL}  --google_analytics_id=$alexa_seq_conf{GOOGLE_ANALYTICS_ID}";

#15-E.) Create candidate gene lists (AE+DE genes, AE genes only, Gains only, Losses only)
print COMMANDS "\n\n#15-E.) Create candidate gene lists (AE+DE genes, AE genes only, Gains only, Losses only)";
print MAGENTA, "\n\t15-E.) Create candidate gene lists (AE+DE genes, AE genes only, Gains only, Losses only)", RESET;
print COMMANDS "\n#Execute the following to create candidate lists for each comparison defined";
my $bash_file22b = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/website/generateCandidateLists.sh";
open(BASH, ">$bash_file22b") || die "\n\nCould not open bash output file: $bash_file22b\n\n";

foreach my $comparison (sort {$a <=> $b} keys %{$comparisons_ref}){
  my $comparison_name = $comparisons_ref->{$comparison}->{comparison_id};
  my $libA_id = $comparisons_ref->{$comparison}->{libraryA_id};
  my $libB_id = $comparisons_ref->{$comparison}->{libraryB_id};
  my $comp_lib_id = "$libA_id"."_vs_"."$libB_id";

  print BASH "$script_dir/website/generateCandidateLists.pl   --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --project=$project_conf{PROJECT_NAME}  --comparison_id=$comp_lib_id  --comparison_name=$comparison_name  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --annotation_dir=$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/  --junction_seq_size=$project_conf{JUNCTION_SEQ_LENGTH}  --boundary_seq_size=$project_conf{BOUNDARY_SEQ_LENGTH}  --partition_file=$region_file  --alexa_home_path=$alexa_seq_conf{ALEXA_HOME_URL}  --alexa_seq_path=$alexa_seq_conf{ALEXA_SEQ_URL}  --search_page_url=$alexa_seq_conf{XAPIAN_URL}  --google_analytics_id=$alexa_seq_conf{GOOGLE_ANALYTICS_ID}\n";
}
close(BASH);
print COMMANDS "\nbash $bash_file22b";


#15-F.) Create candidate peptide lists (junction gains, junction losses, exon+junction gains, exon+junction losses)
print COMMANDS "\n\n#15-F.) Create candidate peptide lists (junction gains, junction losses, exon+junction gains, exon+junction losses)";
print MAGENTA, "\n\t15-F.) Create candidate peptide lists (junction gains, junction losses, exon+junction gains, exon+junction losses)", RESET;
print COMMANDS "\n#Execute the following to create peptide lists for each comparison defined";
my $bash_file22c = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/website/generatePeptideLists.sh";
open(BASH, ">$bash_file22c") || die "\n\nCould not open bash output file: $bash_file22c\n\n";
foreach my $comparison (sort {$a <=> $b} keys %{$comparisons_ref}){
  my $comparison_name = $comparisons_ref->{$comparison}->{comparison_id};
  my $libA_id = $comparisons_ref->{$comparison}->{libraryA_id};
  my $libB_id = $comparisons_ref->{$comparison}->{libraryB_id};
  my $comp_lib_id = "$libA_id"."_vs_"."$libB_id";

  print BASH "$script_dir/website/generatePeptideLists.pl  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --peptide_size=15  --project=$project_conf{PROJECT_NAME}  --comparison_id=$comp_lib_id  --comparison_name=$comparison_name  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --annotation_dir=$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/  --junction_seq_size=$project_conf{JUNCTION_SEQ_LENGTH}  --partition_file=$region_file  --alexa_home_path=$alexa_seq_conf{ALEXA_HOME_URL}  --alexa_seq_path=$alexa_seq_conf{ALEXA_SEQ_URL}  --search_page_url=$alexa_seq_conf{XAPIAN_URL}  --google_analytics_id=$alexa_seq_conf{GOOGLE_ANALYTICS_ID}\n";

}
close(BASH);
print COMMANDS "\nbash $bash_file22c";


#15-G.) Create gene-by-gene summaries (statistics and graphics specific to individual genes)
print COMMANDS "\n\n#15-G.) Create gene-by-gene summaries (statistics and graphics specific to individual genes)";
print MAGENTA, "\n\t15-G.) Create gene-by-gene summaries (statistics and graphics specific to individual genes)", RESET;
my $bash_file23a = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/website/createGeneHtmlJobs.sh";
my $bash_file23b = "$alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/website/generateGeneHtml.sh";

my $create_jobs_script = "$script_dir/website/createGeneHtmlJobs.pl";
my $shell = "$alexa_seq_conf{SHELL_FILE}";
my $gene_script = "$script_dir/$project_conf{WEBSITE_VERSION}/generateGeneHtml.pl";
my $annotation_dir = "$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/";
my $entrez_annotations = "$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/$project_conf{SPECIES_NAME}.gene_info";
my $partition_file = $region_file;
my $partition_file_version = $project_conf{REGIONS_FILE_VERSION};
my $track_url = "$alexa_seq_conf{WEB_URL}/$project_conf{PROJECT_NAME}/ucsc/";
my $temp_dir = "$alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME}/genes/";
my $species_name_spaced = $project_conf{SPECIES_NAME};
$species_name_spaced =~ s/_/ /g;

open (BASH, ">$bash_file23a") || die "\n\nCould not open bash output file: $bash_file23a\n\n";
print BASH "\n$script_dir/website/createGeneHtmlJobs.pl  --shell=$alexa_seq_conf{SHELL_FILE}  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}  --gene_script=$gene_script  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --annotation_dir=$annotation_dir  --junction_seq_size=$project_conf{JUNCTION_SEQ_LENGTH}  --boundary_seq_size=$project_conf{BOUNDARY_SEQ_LENGTH}  --entrez_annotations=$entrez_annotations  --paths_file=$project_library_paths_file  --ensembl_version=$project_conf{ENSEMBL_VERSION}  --ucsc_build=$project_conf{UCSC_BUILD}  --species=\"$species_name_spaced\"  --species_name=$project_conf{SPECIES_NAME_COMMON}  --project_name=\"$project_conf{PROJECT_NAME}\"  --partition_file=$partition_file  --partition_file_version=$partition_file_version  --track_url=$track_url  --temp_dir=$temp_dir  --search_page_url=$alexa_seq_conf{XAPIAN_URL}  --alexa_home_path=$alexa_seq_conf{ALEXA_HOME_URL}  --alexa_seq_path=$alexa_seq_conf{ALEXA_SEQ_URL}  --google_analytics_id=$alexa_seq_conf{GOOGLE_ANALYTICS_ID}  --bash_file=$bash_file23b";
close(BASH);

print COMMANDS "\nbash $bash_file23a\n";

if ($cluster_commands){
  print COMMANDS "\n#Execute the following shell script to create gene-by-gene summaries.  Submit to cluster...";
  print COMMANDS "\ncd $alexa_seq_conf{ANALYSIS_DIR}/batch_jobs/$project_conf{PROJECT_NAME}/website/\nrm -fr GGH\nmqsub --file $bash_file23b  --name GGH  --mkdir  --delay $alexa_seq_conf{MEDIUM_DELAY}";
}else{
  print COMMANDS "\n#Execute the following or submit the contained jobs to your cluster";
  print COMMANDS "\nbash $bash_file23b";
}
print COMMANDS "\n\n#Check for completion of the correct number of jobs:";
print COMMANDS "\ngrep -v Chromosome $region_file | wc -l";
print COMMANDS "\ngrep \"SCRIPT COMPLETE\" $alexa_seq_conf{ANALYSIS_DIR}/logs/website/$project_conf{PROJECT_NAME}/genes/*LOG.txt | wc -l";


#15-H.) Create links between gene summaries and their neighbourhoods
print COMMANDS "\n\n#15-H.) Create links between gene summaries and their neighbourhoods";
print MAGENTA, "\n\t15-H.) Create links between gene summaries and their neighbourhoods", RESET;
print COMMANDS "\n#Execute the following to generate links";

print COMMANDS "\n$script_dir/website/generateLinksHtml.pl  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --annotation_dir=$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/  --junction_seq_size=$project_conf{JUNCTION_SEQ_LENGTH}  --boundary_seq_size=$project_conf{BOUNDARY_SEQ_LENGTH}  --partition_file=$alexa_seq_conf{SEQUENCE_DB_DIR}/$project_conf{ALEXA_SEQ_DB}/Regions_"."$project_conf{REGIONS_FILE_VERSION}"."_Genes.txt  --genes_dir=$alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME}/genes/  --search_page_url=$alexa_seq_conf{XAPIAN_URL}  --alexa_home_path=$alexa_seq_conf{ALEXA_HOME_URL}  --alexa_seq_path=$alexa_seq_conf{ALEXA_SEQ_URL}  --google_analytics_id=$alexa_seq_conf{GOOGLE_ANALYTICS_ID}"; 

#15-I.) Copy .png and .gff images, CSS files, java script files into the necessary directories ...
print COMMANDS "\n\n#15-I.) Copy .png and .gff images, CSS files, java script files into the necessary directories ...";
print MAGENTA, "\n\t15-I.) Copy .png and .gff images, CSS files, java script files into the necessary directories ...", RESET;
print COMMANDS "\n#Execute the following to copy the necessary files:";
print COMMANDS "\ncp -f $script_dir/website/web_files/*.css $script_dir/website/web_files/*.js $alexa_seq_conf{ANALYSIS_DIR}/temp/website/";
print COMMANDS "\ncp -f $script_dir/website/web_files/images/* $alexa_seq_conf{ANALYSIS_DIR}/temp/website/images/";
print COMMANDS "\ncp -f $script_dir/website/web_files/project/images/* $alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME}/images/";


#15-J.) Copy all web files from the temp dir into a public WWW directory (i.e. a directory served up by your web server)
print COMMANDS "\n\n#15-J.) Copy all web files from the temp dir into a public WWW directory (i.e. a directory served up by your web server)";
print MAGENTA, "\n\t15-J.) Copy all web files from the temp dir into a public WWW directory (i.e. a directory served up by your web server)", RESET;
print COMMANDS "\n#Execute the following to copy the necessary files:";

print COMMANDS "\ncp -f $script_dir/website/web_files/*.css $script_dir/website/web_files/*.js $alexa_seq_conf{WEB_DIR}/";
if ($demo){
  print COMMANDS "\nmv -f $alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME} $alexa_seq_conf{WEB_DIR}/";
  print COMMANDS "\n$cmd1";

}else{
  print COMMANDS "\ncp -fr $alexa_seq_conf{ANALYSIS_DIR}/temp/website/$project_conf{PROJECT_NAME} $alexa_seq_conf{WEB_DIR}/";
}
print COMMANDS "\nmkdir $alexa_seq_conf{WEB_DIR}/images/";
print COMMANDS "\ncp -fr $alexa_seq_conf{ANALYSIS_DIR}/temp/website/images/ $alexa_seq_conf{WEB_DIR}/";

#16.) Index all gene pages using the Xapain-Omega utility
print COMMANDS "\n\n\n#16.) Index all gene pages using the Xapain-Omega utility";
print MAGENTA, "\n\n16.) Index all gene pages using the Xapain-Omega utility", RESET;

print COMMANDS "\n\n#Make sure the paths in the conf file are correct:";
print COMMANDS "\n# $alexa_seq_conf{XAPIAN_DIR}/omega.conf";

print COMMANDS "\n\n#Update the search template page:";
print COMMANDS "\n# $alexa_seq_conf{XAPIAN_DIR}/searchdata/templates/query";

print COMMANDS "\n\n#Create the index files as follows (excluding various file types from indexing):";
print COMMANDS "\n#Delete the files in the data dir to create a fresh index or leave in place to add to the current index";
print COMMANDS "\ncd $alexa_seq_conf{XAPIAN_DIR}/searchdata/data/default/";
print COMMANDS "\n#rm -f *";
print COMMANDS "\n$alexa_seq_conf{XAPIAN_BIN} --duplicates ignore --db $alexa_seq_conf{XAPIAN_DIR}/searchdata/data/default/ $alexa_seq_conf{WEB_DIR}/$project_conf{PROJECT_NAME}/genes --url $alexa_seq_conf{WEB_URL}/$project_conf{PROJECT_NAME}/genes --mime-type doc: --mime-type jpg: --mime-type gz: --mime-type svg: --mime-type svgz: --mime-type txt: --mime-type xls: --mime-type gif: --mime-type html: ";

print COMMANDS "\n\n#Once the index is complete, test the search application by loading the following URL in your browser:";
print COMMANDS "\n#$alexa_seq_conf{XAPIAN_URL}";

print COMMANDS "\n\n#You can customize the look of this search page by modifying the Xapian-Omega template files provided by the developers";
print COMMANDS "\n#You can also start this customization by using the template files we modified.  These are stored here:";
print COMMANDS "\n#$script_dir/website/web_files/xapian_omega/templates/";

#17.) Clean-up the analysis temp files
print COMMANDS "\n\n\n#17.) Clean-up the analysis temp files";
print MAGENTA, "\n\n17.) Clean-up the analysis temp files", RESET;
print COMMANDS "\n#Execute the following to clean up temp files";
print COMMANDS "\n#IMPORTANT NOTE: Only do this when you are confidant that the analysis is complete and has run correctly (otherwise you will have to redo several steps)";

print COMMANDS "\n$script_dir/process/cleanUp.pl  --analysis_dir=$alexa_seq_conf{ANALYSIS_DIR}/  --project=$project_conf{PROJECT_NAME}  --library_conf_file=$project_config_file  --ensembl_version=$project_conf{ENSEMBL_VERSION}"; 
print COMMANDS "\n\n";
close(COMMANDS);


print YELLOW, "\n\nAll necessary commands to perform ALEXA-Seq analysis for the project '$project_conf{PROJECT_NAME}' have been written to following commands file:\n$commands_file\n\nFollow the instructions in this file to proceed...\n\n", RESET;

exit();

