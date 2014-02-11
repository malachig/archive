#!/usr/bin/perl -w

#Written by Malachi Griffith
#Take an input file of probe info (exon junction probes in this example), divide into into smaller blocks and store these in:
#  /home/user/alexa/ALEXA_version/specificity/exon_junction/probes
#Also output the probe sequences as a fasta file in:
#  /home/user/alexa/ALEXA_version/specificity/exon_junction/probe_fasta
#Then create blast commands to blast these against a database of mRNAs or ESTs. Resulting blast result files will be stored in:
#  /home/user/alexa/ALEXA_version/specificity/exon_junction/blast_results
#This script can also be used to blast probes against a complete genome (single contig for each chromosome)

#Finally create commands to run testProbeSpecificity.pl which will process the blast results files and store relevant info in
#new appended probe files in:
#In the case of negative probes the results are instead parsed by: summarizeBlastResults_NegativeProbes.pl

#    /home/user/alexa/ALEXA_version/specificity/exon_junction/parsed_results

#The two sets of commands for the cluster will be stored as bash scripts in:
#    /home/user/alexa/ALEXA_version/batch_scripts/specificity

#These two sets of commands can be run on the cluster
#Make sure the batch file is okay and run the following command:
#mqsub.py --file $batch_file.sh --name $job_name --mkdir
#This command will likely depend on how your cluster works

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

#ALEXA libraries
#When a script is initiated, use the full path of the script location at execution to add the perl module libraries to @INC
#This should allow this scripts to work regardless of the current working directory or the script location (where it was unpacked).
#The /utilities directory must remain in the same directory as this script but the entire code directory can be moved around
BEGIN {
  my $script_dir = &File::Basename::dirname($0);
  push (@INC, $script_dir);
}
use utilities::utility qw(:all);

#Initialize command line options
my $database = '';           #DB specification will be passed on to parsing script
my $server = '';             #DB specification will be passed on to parsing script
my $user = '';               #DB specification will be passed on to parsing script
my $password = '';           #DB specification will be passed on to parsing script
my $db_type = '';            #Database type ('EST', 'mRNA', or 'genomic')
my $master_probe_file = '';  #Input probe file
my $blast_bin = '';          #Full path to the blast executable to be used in the blast command - megablast or blastall
my $blast_database = '';     #Full path to blast database created with the same version of BLAST specied as blast_bin
my $word_size = '';          #Word search size (use a multiple of 4 for megablast (say 12), default for blastall is 11)
my $temp_dir = '';           #Directory where all other temp folders will be created
my $job_name = '';           #Use to create blast batch job file names
                             # - (eg. mRNA_vs_exonJunction -> blastall_mRNA_vs_exonJunction.sh + specificity_mRNA_vs_exonJunction.sh)
my $batch_dir = '';          #Directory to store batch job files
my $probe_parse_bin = '';    #Full path to script to be used to parse blast results
my $logfile = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 
	    'db_type=s'=>\$db_type, 'master_probe_file=s'=>\$master_probe_file, 'blast_bin=s'=>\$blast_bin, 'blast_database=s'=>\$blast_database,
	    'word_size=i'=>\$word_size, 'temp_dir=s'=>\$temp_dir, 'job_name=s'=>\$job_name, 'batch_dir=s'=>\$batch_dir,
	    'probe_parse_bin=s'=>\$probe_parse_bin, 'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script takes a probe file as input, breaks it into pieces and creates batch jobs for a cluster", RESET;
print GREEN, "\n\tTwo jobs are created.  One to conduct blast jobs, the other to parse and summarize the results of these", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password\n", RESET;
print GREEN, "\n\tSpecify the type of blast database being used using: --db_type ('EST', 'mRNA', 'genomic', 'enst', 'probe', or 'linker')", RESET;
print GREEN, "\n\tSpecify a tab-delimited input probe file using: --master_probe_file", RESET;
print GREEN, "\n\tSpecify the full path to the blast binary to be used (megablast or blastall) using: --blast_bin", RESET;
print GREEN, "\n\tSpecify the full path to the blast database to be used (created with same version of blast!) using: --blast_database", RESET;
print GREEN, "\n\tSpecify the word size to be used for searching using: --word_size", RESET;
print GREEN, "\n\t\tFor megablast try 12.  For blastall try 11.", RESET;
print GREEN, "\n\tSpecify the a temp directory for creating all temp probe, blast and results files using: --temp_dir", RESET;
print GREEN, "\n\tSpecify a name for this job using: --job_name (eg. mRNA_vs_exonJunction)", RESET;
print GREEN, "\n\tSpecify the directory for batch files to be stored using: --batch_dir", RESET;
print GREEN, "\n\tSpecify which script should be used to parse the results using: --probe_parse_bin", RESET;
print GREEN, "\n\t\tFor normal experimental probes use --probe_parse_bin=/home/user/alexa/Array_design/testProbeSpecificity_EST-mRNA.pl", RESET;
print GREEN, "\n\t\tFor negative probes (testing to make sure they do not match any expressed sequences at all) use testProbeSpecificity_Negative.pl";
print GREEN, "\n\tSpecify a logfile for output using: --logfile", RESET;
print GREEN, "\n\nExample: blastProbesVersus_Database.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --db_type=mRNA  --master_probe_file=/home/user/alexa/ALEXA_version/unfiltered_probes/exonJunctionProbes.txt  --blast_bin=/home/user/BioSw/BLAST2/blast2.2.15/bin/blastall  --blast_database=/home/user/alexa/ALEXA_version/data_sources/ucsc/blast_db/human_hg18_16Jan2007/human_mrna  --word_size=14  --temp_dir=/home/user/alexa/ALEXA_version/specificity/  --job_name=mRNA_vs_exonJunction  --batch_dir=/home/user/alexa/ALEXA_version/batch_scripts/specificity/exonJunction/  --probe_parse_bin=/home/user/alexa/Array_design/testProbeSpecificity_EST-mRNA.pl  --logfile=/home/user/alexa/ALEXA_version/logs/specificity/blastProbesVersus_mRNA_exonJunction_LOG.txt\n\n", RESET;

unless ($database && $server && $user && $password && $db_type && $master_probe_file && $blast_bin && $blast_database && $word_size && $temp_dir && $job_name && $batch_dir && $probe_parse_bin && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Open logfile for output
open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";
print LOG "\nUser Specified the following options:\ndatabase = $database\ndb_type = $db_type\nmaster_probe_file = $master_probe_file\nblast_bin = $blast_bin\nblast_database = $blast_database\nword_size = $word_size\ntemp_dir = $temp_dir\njob_name = $job_name\nbatch_dir = $batch_dir\nprobe_parse_bin = $probe_parse_bin\nlogfile = $logfile\n\n";

#Examples of blast binaries to use:
#/home/user/BioSw/BLAST2/blast2.2.14/bin/megablast
#/home/user/BioSw/BLAST2/blast2.2.14/bin/blastall

#1.) Check the blast binary path and add the program name to the blast command if blastall was chosen
unless (-e $blast_bin){
  print RED, "\nSpecified path to blast binary is not valid!!\n\n", RESET;
  close (LOG);
  exit();
}
if ($blast_bin =~ /blastall/){
  $blast_bin = "$blast_bin"." -p blastn";
}

#2.) Remind the user that the blast_database should match the type specified
if ($db_type eq "mRNA" || $db_type eq "EST" || $db_type eq "genomic" || $db_type eq "enst" || $db_type eq "probe" || $db_type eq "linker"){
  print YELLOW, "\nYou specified:\n\tdb_type = $db_type.\n\tblast_database = $blast_database\n\tALEXA_database = $database\n\tProbe_file = $master_probe_file\n\nDo these match each other (y/n)?? ", RESET;
  my $answer = <>;
  chomp($answer);
  unless($answer eq "y" || $answer eq "Y"){
    close (LOG);
    exit();
  }
}else{
  print RED, "\nMust specify a valid database type ('mRNA', 'EST', 'genomic', 'enst', 'probe', or 'linker')\n\n", RESET;
  close (LOG);
  exit();
}

#3.) Note: genomic option should only be used in conjunction with the NegativeControl parsing script
if ($db_type eq "genomic"){
  unless ($probe_parse_bin =~ /Negative\.pl/ || $probe_parse_bin =~ /testProbeSpecificity_Genomic\.pl/){
    print RED, "\nThe 'genomic' option should only be used in conjuction with testProbeSpecificity_Negative.pl or testProbeSpecificity_Genomic.pl!!\n\n", RESET;
    close (LOG);
    exit();
  }
}

if ($db_type eq "linker"){
  unless ($probe_parse_bin =~ /Negative\.pl/){
    print RED, "\nThe 'linker' option should only be used in conjuction with testProbeSpecificity_Negative.pl!!\n\n", RESET;
    close (LOG);
    exit();
  }
}

#4.) Define the number of probe sequences to write to each fasta file and probe file
my $block_size = 10000;

#5.) Remind user of the probe specificity test script that will be used to parse the blast results
print YELLOW, "\n\n$probe_parse_bin will be used to parse and summarize blast results.", RESET;
print YELLOW, "\nNote: testProbeSpecificity_Negative.pl should only be used for NegativeControl or Intron probes!!\n\n", RESET;

#6.) Test specified batch file directory and open output batch scripts
unless ($batch_dir =~ /.*\/$/){
  $batch_dir = "$batch_dir"."/";
}
unless (-e $batch_dir && -d $batch_dir){
  print RED, "\nSpecified batch directory: $batch_dir does not appear valid! Create this directory before proceeding\n\n", RESET;
  close (LOG);
  exit();
}
my $blast_batch = "$batch_dir"."blastall_"."$job_name".".sh";
my $specificity_batch = "$batch_dir"."specificity_"."$job_name".".sh";
open (BLAST_SH, ">$blast_batch") || die "\nCould not open blast.sh file: $blast_batch";
open (SPEC_SH, ">$specificity_batch") || die "\nCould not open specificity.sh file: $specificity_batch";

#7.) Prepare the working directory
#    Create $fasta_dir (probe_fasta), $probe_dir (probes), $results_dir (blast_results), and $parsed_results_dir (parsed_results)
unless ($temp_dir =~ /.*\/$/){
  $temp_dir = "$temp_dir"."/";
}
unless (-e $temp_dir && -d $temp_dir){
  print RED, "\nSpecified working directory: $temp_dir does not appear valid! Create this directory before proceeding\n\n", RESET;
  close (LOG);
  exit();
}
my $sub_dir = &createNewDir('-path'=>$temp_dir, '-new_dir_name'=>$job_name);
my $fasta_dir = &createNewDir('-path'=>$sub_dir, '-new_dir_name'=>"probe_fasta");
my $probe_dir = &createNewDir('-path'=>$sub_dir, '-new_dir_name'=>"probes");
my $results_dir = &createNewDir('-path'=>$sub_dir, '-new_dir_name'=>"blast_results");
my $parsed_results_dir = &createNewDir('-path'=>$sub_dir, '-new_dir_name'=>"parsed_results");

#8.) Conduct a preliminary test of the master input probe file 
#    Find out the number of lines, test file headers and determine location of probe sequences
open (PROBEFILE, "$master_probe_file") || die "\nCould not open probe file: $master_probe_file";
my $num_probes = 0;
my $probe_column = 0;
my $probe_id_column_found = 0;
my $sequence_column_found = 0;
my $first_line = 1;

while (<PROBEFILE>){
  my $line_record = $_;
  chomp ($line_record);

  if ($first_line == 1){
    my @line = split ("\t", $line_record);
    my $probe_column_count = 0;

    foreach my $value (@line){
      $probe_column_count++;
      if ($value =~ /^Probe_Count$/){
	$probe_id_column_found = 1;
      }
      if ($value =~ /^Sequence$/){
	$sequence_column_found = 1;
	$probe_column = $probe_column_count;
      }
    }
    $first_line = 0;

    unless ($probe_id_column_found == 1 && $sequence_column_found == 1){
      print RED, "\nFile: $master_probe_file does not appear to be a valid probe file (expected column missing) - exiting\n\n", RESET;
      close (PROBEFILE);
      close (LOG);
      exit();
    }
    next();
  }
  if ($line_record =~ /^\d+/){
    $num_probes++;
  }
}
print BLUE, "\nFound a total of $num_probes probes in the input probe file\n\n", RESET;
print LOG "\nFound a total of $num_probes probes in the input probe file\n\n";
close (PROBEFILE);

#9.) Generate output files and write batch commands to batch files
my $probe_count = 0;
my $total_probe_count = 0;
my $block_num = 1;

my ($fasta_file, $probe_sub_file, $blast_results_file, $parsed_results_file, $header_line);

#Open the master input probe file for processing
open (PROBEFILE, "$master_probe_file") || die "\nCould not open probe file: $master_probe_file";

#Open the first output fasta file
$fasta_file = "$fasta_dir"."probes_"."$block_num.fa";
open (FASTA_FILE, ">$fasta_file") || die "\nCould not open temp fasta file: $fasta_file";

#Open the first probe sub-file
$probe_sub_file = "$probe_dir"."probes_"."$block_num";
open (PROBE_SUB_FILE, ">$probe_sub_file") || die "\nCould not open probe sub file: $probe_sub_file";

#Figure out the name of the first output blast results file
$blast_results_file = "$results_dir"."blastResults_"."$db_type"."_$block_num";

#Figure out the name of the first parsed results file
$parsed_results_file = "$parsed_results_dir"."probes_specificity_"."$db_type"."_$block_num";

#Print out the BLAST command to be used
print BLAST_SH "$blast_bin -d $blast_database -i $fasta_file -o $blast_results_file -F F -W $word_size -m 8\n";

#Print out the PARSING command to be used - Tailor it to the parsing script specified
if ($probe_parse_bin =~ /Negative\.pl/){
  print SPEC_SH "$probe_parse_bin --probe_file $probe_sub_file --blast_file $blast_results_file --blast_type $db_type --outfile $parsed_results_file\n";
}elsif($probe_parse_bin =~ /testProbeSpecificity_Genomic\.pl/){
  print SPEC_SH "$probe_parse_bin  --probe_infile $probe_sub_file  --blast_results_infile $blast_results_file  --outfile $parsed_results_file\n";
}elsif ($probe_parse_bin =~ /EST\-mRNA\.pl/){
  print SPEC_SH "$probe_parse_bin --database $database --server $server --user $user --password $password --db_type $db_type --probe_infile $probe_sub_file --blast_results_infile $blast_results_file --outfile $parsed_results_file\n";
}elsif ($probe_parse_bin =~ /Ensembl\.pl/){
  print SPEC_SH "$probe_parse_bin --database $database --server $server --user $user --password $password --probe_infile $probe_sub_file --blast_results_infile $blast_results_file --outfile $parsed_results_file\n";
}elsif ($probe_parse_bin =~ /Probes\.pl/){
  print SPEC_SH "$probe_parse_bin --probe_infile $probe_sub_file --blast_results_infile $blast_results_file --outfile $parsed_results_file\n";
}else{
  print RED, "\nSpecified parsing script: $probe_parse_bin not understood!!\n\n", RESET;
  close (LOG);
  exit();
}

$first_line = 1;
while (<PROBEFILE>){

  #Grab the header line
  if ($first_line == 1){
    chomp ($_);
    $header_line = $_;
    $first_line = 0;
    #Print the header to the first probe sub file
    print PROBE_SUB_FILE "$header_line\n";
    next();
  }

  $probe_count++;
  $total_probe_count++;

  my $line_record = $_;
  chomp ($line_record);
  my @line = split ("\t", $_);
  my $probe_id = $line[0];
  my $probe_seq = $line[$probe_column-1];

  print FASTA_FILE ">$probe_id\n$probe_seq\n";
  print PROBE_SUB_FILE "$line_record\n";

  #Watch for end of file
  if ($total_probe_count == $num_probes){
    next();
  }

  #Once a full block is reached start on the next file
  if($probe_count == $block_size){

    #Close the current output files
    close (FASTA_FILE);
    close (PROBE_SUB_FILE);

    $block_num++;
    $fasta_file = "$fasta_dir"."probes_"."$block_num.fa";
    $probe_sub_file = "$probe_dir"."probes_"."$block_num";

    #Open the new output fasta file
    open (FASTA_FILE, ">$fasta_file") || die "\nCould not open temp fasta file: $fasta_file";

    #Open the new probe sub file
    open (PROBE_SUB_FILE, ">$probe_sub_file") || die "\nCould not open probe sub file: $probe_sub_file";

    #Print the header to this new probe sub file
    print PROBE_SUB_FILE "$header_line\n";

    #Figure out the name of the new output file
    $blast_results_file = "$results_dir"."blastResults_"."$db_type"."_$block_num";

    #Figure out the name of the new parsed results file
    $parsed_results_file = "$parsed_results_dir"."probes_specificity_"."$db_type"."_$block_num";

    #Print out the BLAST command to be used 
    print BLAST_SH "$blast_bin -d $blast_database -i $fasta_file -o $blast_results_file -F F -W $word_size -m 8\n";

    #Print out the PARSING command to be used - Tailor it to the parsing script specified
    if ($probe_parse_bin =~ /Negative\.pl/){
      print SPEC_SH "$probe_parse_bin --probe_file $probe_sub_file --blast_file $blast_results_file --blast_type $db_type --outfile $parsed_results_file\n";
    }elsif($probe_parse_bin =~ /testProbeSpecificity_Genomic\.pl/){
      print SPEC_SH "$probe_parse_bin  --probe_infile $probe_sub_file  --blast_results_infile $blast_results_file  --outfile $parsed_results_file\n";
    }elsif ($probe_parse_bin =~ /EST\-mRNA\.pl/){
      print SPEC_SH "$probe_parse_bin --database $database --server $server --user $user --password $password --db_type $db_type --probe_infile $probe_sub_file --blast_results_infile $blast_results_file --outfile $parsed_results_file\n";
    }elsif ($probe_parse_bin =~ /Ensembl\.pl/){
      print SPEC_SH "$probe_parse_bin --database $database --server $server --user $user --password $password --probe_infile $probe_sub_file --blast_results_infile $blast_results_file --outfile $parsed_results_file\n";
    }elsif ($probe_parse_bin =~ /Probes\.pl/){
      print SPEC_SH "$probe_parse_bin --probe_infile $probe_sub_file --blast_results_infile $blast_results_file --outfile $parsed_results_file\n";
    }else{
      print RED, "\nSpecified parsing script: $probe_parse_bin not understood!!\n\n", RESET;
      close (LOG);
      exit();
    }

    #Reset variables
    $probe_count = 0;
    $blast_results_file = '';
    $probe_sub_file = '';
    $parsed_results_file = '';
  }
}
close (FASTA_FILE);
close (PROBE_SUB_FILE);
close (PROBEFILE);
close (BLAST_SH);
close (SPEC_SH);

print BLUE, "\nPrinted a total of $total_probe_count probes to $block_num files in $sub_dir\n\n", RESET;
print LOG "\nPrinted a total of $total_probe_count probes to $block_num files in $sub_dir\n\n";

print BLUE, "Created the following batch scripts to conduct the specificity test for the supplied probe file:\n\t$blast_batch\n\t$specificity_batch\n\n", RESET;
print LOG "Created the following batch scripts to conduct the specificity test for the supplied probe file:\n\t$blast_batch\n\t$specificity_batch\n\n";

close (LOG);

exit();
