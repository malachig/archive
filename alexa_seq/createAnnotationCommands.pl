#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#This script takes an ALEXA-seq configuration file and Project specific configuration file as input and uses these to generate all necessary commands for creation of an ALEXA-Seq annotation database
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
my $annotation_config_file = '';
my $commands_file = '';
my $clean = '';
my $cluster_commands = '';

GetOptions ('alexa_seq_config_file=s'=>\$alexa_seq_config_file, 'annotation_config_file=s'=>\$annotation_config_file, 'commands_file=s'=>\$commands_file, 'clean=s'=>\$clean, 'cluster_commands=s'=>\$cluster_commands);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script takes an ALEXA-Seq configuration file and Project specific configuration file as input and uses these to generate all necessary commands for ALEXA-Seq analysis", RESET;
print GREEN, "\n\tSpecify the path to your ALEXA-Seq configuration file using: --alexa_seq_config_file", RESET;
print GREEN, "\n\tSpecify the path to your Annotation configuration file using: --annotation_config_file", RESET;
print GREEN, "\n\tSpecify the file to which all commands will be written using: --commands_file", RESET;

print GREEN, "\n\nExample: createAnnotationCommands.pl  --alexa_seq_config_file=/home/malachig/svn/alexa_seq/config_files/ALEXA_Seq_PIPELINE.conf  --annotation_config_file=/home/malachig/svn/alexa_seq/config_files/annotation/ALEXA_Seq_Annotation_hs_53_36o.conf  --commands_file=/home/malachig/svn/alexa_seq/config_files/annotation/ALEXA_Seq_Annotation_hs_53_36o.commands\n\n", RESET;

unless ($alexa_seq_config_file && $annotation_config_file && $commands_file){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

unless (-e $alexa_seq_config_file && -e $annotation_config_file){
  print RED, "\n\nCould not verify existence of both config files ... check paths\n\n", RESET;
  exit();
}

#0.) Import config files
open(COMMANDS, ">$commands_file") || die "\n\nCould not open output commands file: $commands_file\n\n";
print MAGENTA, "\n\nImporting config files", RESET;
print COMMANDS "#ALEXA-SEQ ANNOTATION - DETAILED INSTRUCTIONS";
print COMMANDS "\n#Using the following ALEXA-seq and Annotation config files for this analysis:";
print COMMANDS "\n#$alexa_seq_config_file\n#$annotation_config_file";

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

#Import annotation config file
my %annotation_conf;
open(CONF, "$annotation_config_file") || die "\n\nCould not open file: $annotation_config_file\n\n";
my $lane_count = 0;
my $library_count = 0;
my $comparison_count = 0;
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
  #Watch out for annotation configuration values
  if ($_ =~ /(.*)\=(.*)/){
    $annotation_conf{$1}=$2;
  }
}
close(CONF);


#Perform basic checks of the config files
print MAGENTA, "\n\nChecking config files for problems", RESET;
my $warnings = &checkConfig('-alexa_seq_conf'=>\%alexa_seq_conf, '-annotation_conf'=>\%annotation_conf);

my $alexa_db_name = "ALEXA_"."$annotation_conf{ALEXA_SEQ_DB}";
my $lc_name = lc($annotation_conf{SPECIES_NAME});
my $ensembl_code;
if ($annotation_conf{ALEXA_SEQ_DB} =~ /^\w{2}\_(.*)/){
  $ensembl_code = $1;
}
my $ensembl_db_name = "$lc_name"."_core_"."$ensembl_code";
my $annotation_path = "$alexa_seq_conf{SEQUENCE_DB_DIR}/$annotation_conf{ALEXA_SEQ_DB}";


#0.) Install the EnsEMBL API specified in the config file if necessary
print MAGENTA, "\n\n0.) Install the EnsEMBL API specified in the config file if necessary", RESET;
print COMMANDS "\n\n\n#0.) Install the EnsEMBL API specified in the config file if necessary";
print COMMANDS "\n#Execute the following command to install the EnsEMBL API for version: $annotation_conf{ENSEMBL_VERSION}";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/installEnsemblAPI.pl  --install_dir=$script_dir/ensembl_api/  --ensembl_version=$annotation_conf{ENSEMBL_VERSION}";

#1.) Install the EnsEMBL database if you have not already done so (grab the latest version or a specific version if you like).  Directions are provided here:
print MAGENTA, "\n\n1.) Install the EnsEMBL database if you have not already done so (grab the latest version or a specific version if you like).", RESET;
print COMMANDS "\n\n\n#1.) Install the EnsEMBL database if you have not already done so (grab the latest version or a specific version if you like).  Directions are provided here:";
print COMMANDS "\n#http://ensembl.org/info/docs/webcode/install/ensembl-data.html";
print COMMANDS "\n#If your centre maintains a local EnsEMBL database server, you can skip this step";
print COMMANDS "\n#According to your config files the target database should be: $ensembl_db_name";

#2.) Set up an ALEXA database for the desired version of the EnsEMBL database
print MAGENTA, "\n\n2.) Set up an ALEXA database for the desired version of the EnsEMBL database", RESET;
print COMMANDS "\n\n\n#2.) Set up an ALEXA database for the desired version of the EnsEMBL database";
print COMMANDS "\n#Log into a mysql server (with a user that can create databases) and perform the following commands:";
print COMMANDS "\nmysql -h $alexa_seq_conf{ALEXA_SERVER} -u $alexa_seq_conf{ALEXA_USER2} -p$alexa_seq_conf{ALEXA_PASSWORD2}";
print COMMANDS "\nDROP DATABASE $alexa_db_name;";
print COMMANDS "\nCREATE DATABASE $alexa_db_name;";
print COMMANDS "\nUSE $alexa_db_name;";
print COMMANDS "\nSOURCE $script_dir/sql/ALEXA_schema.sql;";
print COMMANDS "\nexit;";

#3.) Create needed target directories
my $cmd1 = "$script_dir/process/createAnnotationDirs.pl  --annotation_dir=$alexa_seq_conf{SEQUENCE_DB_DIR}/  --species_build=$annotation_conf{ALEXA_SEQ_DB}";

print MAGENTA, "\n\n3.) Create all necessary directories for this annotation ($annotation_conf{ALEXA_SEQ_DB}) - directories already present will be untouched!", RESET;
print COMMANDS "\n\n\n#3.) Create all necessary directories for this annotation ($annotation_conf{ALEXA_SEQ_DB}) - directories already present will be untouched!";
print COMMANDS "\n#This command should have already been run, but it is reproduced here for reference (no harm will be done be re-running it)";
print COMMANDS "\n$cmd1";

#This command actually needs to be executed as well as written to the commands file for reference
#print BLUE, "\n\n$cmd1", RESET;
#system($cmd1);



#4-A.) Get repeat elements for each species
print MAGENTA, "\n\n4-A.) Get repeat elements for each species (requires account at RepBase)", RESET;
print COMMANDS "\n\n\n#4-A.) Get repeat elements for each species";

print COMMANDS "\n# http://www.girinst.org/repbase/update/browse.php";
print COMMANDS "\n# Select RepeatClass 'All'.  Taxon '$annotation_conf{SPECIES_NAME}'.  And chose the '$annotation_conf{SPECIES_NAME} and Ancestral' download option";
print COMMANDS "\n# Download data in both fasta and EMBL format and store in the following dir:";
print COMMANDS "\n# $annotation_path/repeats/";
print COMMANDS "\n# Make sure the fasta file is named:  'repeats.fa'";
print COMMANDS "\n# Once you have downloaded these files from RepBase, execute the following:";
print COMMANDS "\n\ncd $annotation_path/repeats/\n";
print COMMANDS "\nperl -ne '\$_ =~ s/\\t/\|/g; print \"\$_\"' repeats.fa > blastdb/repeats.fa";
print COMMANDS "\nmv repeats.fa repeats.txt";
print COMMANDS "\ncd blastdb/";
print COMMANDS "\n$alexa_seq_conf{BLAST_BIN_DIR}/formatdb -i repeats.fa -t repeats -p F -o F -n repeats";
print COMMANDS "\ngzip repeats.fa";
if ($annotation_conf{ALIGNMENT_OPTION} == 2 || $annotation_conf{ALIGNMENT_OPTION} == 3){
  print COMMANDS "\n\n#Now create a bwa database";
  print COMMANDS "\ncp repeats.fa.gz ../bwadb/";
  print COMMANDS "\ncd ../bwadb/";
  print COMMANDS "\n$alexa_seq_conf{BWA_BIN_DIR}/bwa index repeats.fa.gz";
}


#4-B.) Get mRNA/EST alignment tables from UCSC
print MAGENTA, "\n\n4-B.) Get mRNA/EST alignment tables from UCSC", RESET;
print COMMANDS "\n\n\n#4-B.) Get mRNA/EST alignment tables from UCSC";
print COMMANDS "\ncd $annotation_path/mrna_est/";
print COMMANDS "\nwget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/$annotation_conf{UCSC_BUILD}/database/all_mrna.txt.gz' -O all_mrna.txt.gz";
print COMMANDS "\nwget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/$annotation_conf{UCSC_BUILD}/database/all_est.txt.gz' -O all_est.txt.gz";
print COMMANDS "\nwget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/$annotation_conf{UCSC_BUILD}/database/xenoMrna.txt.gz' -O xenoMrna.txt.gz";
print COMMANDS "\nwget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/$annotation_conf{UCSC_BUILD}/database/xenoEst.txt.gz' -O xenoEst.txt.gz";
print COMMANDS "\nwget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/$annotation_conf{UCSC_BUILD}/database/gbCdnaInfo.txt.gz' -O gbCdnaInfo.txt.gz";

#4-C.) Create a database of genbank-to-species mappings
print MAGENTA, "\n\n4-C.) Create a database of genbank-to-species mappings", RESET;
print COMMANDS "\n\n\n#4-C.) Create a database of genbank-to-species mappings";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/createGenBankIdToSpeciesDatabase.pl  --ucsc_align_dir=$annotation_path/mrna_est/   --out_dir=$annotation_path/mrna_est/";


#4-D.) Partition the UCSC alignment files by chromosome
print MAGENTA, "\n\n4-D.) Partition the UCSC alignment files by chromosome", RESET;
print COMMANDS "\n\n\n#4-D.) Partition the UCSC alignment files by chromosome";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/partitionUcscAlignmentFiles.pl  --ucsc_align_dir=$annotation_path/mrna_est/   --out_dir=$annotation_path/mrna_est/partitions/";


#4-E.) Get Entrez gene annotations for this species
print MAGENTA, "\n\n4-E.) Get Entrez gene annotations for this species", RESET;
print COMMANDS "\n\n\n#4-E.) Get Entrez gene annotations for this species";
print COMMANDS "\ncd $annotation_path/";
print COMMANDS "\nwget ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/$annotation_conf{ENTREZ_SPECIES_GROUP}/$annotation_conf{SPECIES_NAME}.gene_info.gz";
print COMMANDS "\ngunzip $annotation_conf{SPECIES_NAME}.gene_info.gz";


#5.) Populate the ALEXA database:
print MAGENTA, "\n\n5.) Populate the ALEXA database:", RESET;
print COMMANDS "\n\n\n#5.) Populate the ALEXA database:";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/getEnsemblGeneData.pl  --ensembl_api_version=$annotation_conf{ENSEMBL_VERSION}  --species=$annotation_conf{SPECIES_NAME_COMMON}  --ensembl_database=$ensembl_db_name  --ensembl_server=$alexa_seq_conf{ENSEMBL_SERVER}  --ensembl_user=$alexa_seq_conf{ENSEMBL_USER}  --ensembl_password=$alexa_seq_conf{ENSEMBL_PASSWORD}  --alexa_database=$alexa_db_name  --alexa_server=$alexa_seq_conf{ALEXA_SERVER}  --alexa_user=$alexa_seq_conf{ALEXA_USER2}  --alexa_password=$alexa_seq_conf{ALEXA_PASSWORD2}  --all_ids=1  --populate_database=yes  --logfile=$annotation_path/logs/getEnsemblGeneData_LOG.txt";

#6.) Get a list of chromosome names for this EnsEMBL build
print MAGENTA, "\n\n6.) Get a list of chromosome names for this EnsEMBL build", RESET;
print COMMANDS "\n\n\n#6.) Get a list of chromosome names for this EnsEMBL build:";
print COMMANDS "\necho 'SELECT DISTINCT chromosome from Gene;' | mysql -h $alexa_seq_conf{ALEXA_SERVER} -u $alexa_seq_conf{ALEXA_USER1} -p$alexa_seq_conf{ALEXA_PASSWORD1} -N $alexa_db_name > $annotation_path/ChromosomeNames.txt";
print COMMANDS "\ncd $annotation_path/";
print COMMANDS "\nls mrna_est/partitions/ | perl -ne 'chomp(\$_); if (\$_ =~ /(chr.*)\\_(est|mrna|xest|xmrna).txt.gz/){print \"\$1\\n\"}' | sort | uniq > ChromosomeNames_UCSC.txt";

print COMMANDS "\n\n#Note 1";
print COMMANDS "\n#*** Compare UCSC and EnsEMBL Chromosome names for each species - compensate for difference somehow  ***";
print COMMANDS "\n#*** Possibly by adding chromosome name re-maps to ALEXA_DB.pm ***";
print COMMANDS "\n\n#Note 2";
print COMMANDS "\n#*** Also for human it may be wise to eliminate haplotype chromsomes from the analysis";
print COMMANDS "\n#Information on these for the latest human build are here: http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/index.shtml";
print COMMANDS "\n#To eliminate these, add filter to getEnsemblGeneData.pl so that they are never imported";


#7-A.) Create jobs to break genomes into manageable pieces (blocks of say 250 genes - without causing breaks within genes)
print MAGENTA, "\n\n7-A.) Create jobs to break genomes into manageable pieces (blocks of say 250 genes - without causing breaks within genes)", RESET;
print COMMANDS "\n\n\n#7-A.) Create jobs to break genomes into manageable pieces (blocks of say 250 genes - without causing breaks within genes):";

print COMMANDS "\n\n#50 gene blocks";
print COMMANDS "\nperl -ne 'chomp(\$_); print \"source $alexa_seq_conf{SHELL_FILE}; $script_dir/alternativeExpressionDatabase/partitionEnsemblGenome.pl  --ensembl_api_version=$annotation_conf{ENSEMBL_VERSION}  --species=$annotation_conf{SPECIES_NAME_COMMON}  --ensembl_database=$ensembl_db_name  --ensembl_server=$alexa_seq_conf{ENSEMBL_SERVER}  --ensembl_user=$alexa_seq_conf{ENSEMBL_USER}  --ensembl_password=$alexa_seq_conf{ENSEMBL_PASSWORD}  --alexa_database=$alexa_db_name  --alexa_server=$alexa_seq_conf{ALEXA_SERVER}  --alexa_user=$alexa_seq_conf{ALEXA_USER1}  --alexa_password=$alexa_seq_conf{ALEXA_PASSWORD1}  --target_gene_count=50  --outfile=$annotation_path/jobs/temp/EnsEMBL_"."$annotation_conf{ENSEMBL_VERSION}"."_Regions_50genes_\$\_.txt  --chr_filter=\$\_\\n\"' $annotation_path/ChromosomeNames.txt > $annotation_path/jobs/partitionEnsemblGenome_50.sh";

print COMMANDS "\n\n#250 gene blocks";
print COMMANDS "\nperl -ne 'chomp(\$_); print \"source $alexa_seq_conf{SHELL_FILE}; $script_dir/alternativeExpressionDatabase/partitionEnsemblGenome.pl  --ensembl_api_version=$annotation_conf{ENSEMBL_VERSION}  --species=$annotation_conf{SPECIES_NAME_COMMON}  --ensembl_database=$ensembl_db_name  --ensembl_server=$alexa_seq_conf{ENSEMBL_SERVER}  --ensembl_user=$alexa_seq_conf{ENSEMBL_USER}  --ensembl_password=$alexa_seq_conf{ENSEMBL_PASSWORD}  --alexa_database=$alexa_db_name  --alexa_server=$alexa_seq_conf{ALEXA_SERVER}  --alexa_user=$alexa_seq_conf{ALEXA_USER1}  --alexa_password=$alexa_seq_conf{ALEXA_PASSWORD1}  --target_gene_count=250  --outfile=$annotation_path/jobs/temp/EnsEMBL_"."$annotation_conf{ENSEMBL_VERSION}"."_Regions_250genes_\$\_.txt  --chr_filter=\$\_\\n\"' $annotation_path/ChromosomeNames.txt > $annotation_path/jobs/partitionEnsemblGenome_250.sh";


#7-B.) Submit these jobs to a cluster
print MAGENTA, "\n\n7-B.) Submit these jobs to a cluster", RESET;
print COMMANDS "\n\n\n#7-B.) Submit these jobs to a cluster:";
if ($cluster_commands){
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}";
  print COMMANDS "\ncd $annotation_path/jobs/";
  print COMMANDS "\nmqsub --file partitionEnsemblGenome_50.sh  --name $annotation_conf{ALEXA_SEQ_DB}_PART_50 --mkdir  --delay 15";
  print COMMANDS "\nmqsub --file partitionEnsemblGenome_250.sh  --name $annotation_conf{ALEXA_SEQ_DB}_PART_250 --mkdir  --delay 15";
}else{
  print COMMANDS "\n$annotation_path/jobs/partitionEnsemblGenome_50.sh";
  print COMMANDS "\n$annotation_path/jobs/partitionEnsemblGenome_250.sh";
}

#7-C.) Check that all jobs completed successfully
print MAGENTA, "\n\n7-C.) Check that all jobs completed successfully", RESET;
print COMMANDS "\n\n\n#7-C.) Check that all jobs completed successfully:";
print COMMANDS "\n\nwc -l $annotation_path/jobs/partitionEnsemblGenome_50.sh";
print COMMANDS "\ngrep Chromosome $annotation_path/jobs/temp/*_50genes*.txt | wc -l";
print COMMANDS "\n\nwc -l $annotation_path/jobs/partitionEnsemblGenome_250.sh";
print COMMANDS "\ngrep Chromosome $annotation_path/jobs/temp/*_250genes*.txt | wc -l";


#7-D.) Concatenate the temp partition files into one file
print MAGENTA, "\n\n7-D.) Concatenate the temp partition files into one file", RESET;
print COMMANDS "\n\n\n#7-D.) Concatenate the temp partition files into one file:";
print COMMANDS "\ncd $annotation_path/";
print COMMANDS "\n\ngrep -h -v Chromosome jobs/temp/*_50genes*.txt > tmp.txt";
print COMMANDS "\nhead -q -n 1 jobs/temp/*_50genes*.txt | sort | uniq | cat - tmp.txt > Regions_50_Genes.txt";
print COMMANDS "\n\ngrep -h -v Chromosome jobs/temp/*_250genes*.txt > tmp.txt";
print COMMANDS "\nhead -q -n 1 jobs/temp/*_250genes*.txt | sort | uniq | cat - tmp.txt > Regions_250_Genes.txt";
print COMMANDS "\nrm -f tmp.txt";

#7-E.) Clean up the temp partition files
print MAGENTA, "\n\n7-E.) Clean up the temp partition files", RESET;
print COMMANDS "\n\n\n#7-E.) Clean up the temp partition files:";
print COMMANDS "\ncd $annotation_path/jobs/temp/";
print COMMANDS "\nrm -f *.txt";


#8.) Create a database file for genes 
print MAGENTA, "\n\n8.) Create a database file for genes", RESET;
print COMMANDS "\n\n\n#8.) Create a database file for genes:";
print COMMANDS "\nrm -f $annotation_path/genes/temp/*.txt";
print COMMANDS "\nperl -ne 'unless(\$_ =~ /Chromosome/){chomp(\$_); \@line=split(\"\\t\", \$_);  print \"$script_dir/alternativeExpressionDatabase/createGeneDatabase.pl  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --chr_filter=\\\"\$line[0]:\$line[1]:\$line[2]-\$line[3]\\\"  --outdir=$annotation_path/genes/temp/\\n\"}' $annotation_path/Regions_250_Genes.txt > $annotation_path/jobs/createGeneDatabase.sh";
print COMMANDS "\n\nbash $annotation_path/jobs/createGeneDatabase.sh";
print COMMANDS "\ncd $annotation_path/genes/";
print COMMANDS "\ngrep -h -v \"Gene_ID\" temp/*.txt | sort -n -k 1 > tmp.txt";
print COMMANDS "\nhead -q -n 1 temp/*.txt | sort | uniq | cat - tmp.txt > genes_annotated.txt";
print COMMANDS "\nrm -f tmp.txt";
print COMMANDS "\ngzip $annotation_path/genes/genes_annotated.txt";
print COMMANDS "\nrm -f temp/*.txt";

#9.) Create a database for Exon Regions 
print MAGENTA, "\n\n9.) Create a database for Exon Regions", RESET;
print COMMANDS "\n\n\n#9.) Create a database for Exon Regions:";
print COMMANDS "\nrm -f $annotation_path/exonRegions/temp/*.txt";
print COMMANDS "\nrm -f $annotation_path/logs/createExonRegionDatabase/*.txt";
print COMMANDS "\nperl -ne 'unless(\$_ =~ /Chromosome/){chomp(\$_); \@line=split(\"\\t\", \$_); print \"source $alexa_seq_conf{SHELL_FILE}; $script_dir/alternativeExpressionDatabase/createExonRegionDatabase.pl  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --chr_filter=\\\"\$line[0]:\$line[1]:\$line[2]-\$line[3]\\\"  --ucsc_align_dir=$annotation_path/mrna_est/partitions/  --genbank_mapfile=$annotation_path/mrna_est/GenBankToOrganism.btree  --wiggle=5  --outfile=$annotation_path/exonRegions/temp/exonRegions_annotated_\$line[0]_\$line[1].txt  --fasta_file=$annotation_path/exonRegions/blastdb/temp/exonRegions_\$line[0]_\$line[1].fa  --logfile=$annotation_path/logs/createExonRegionDatabase/createExonRegionDatabase_\$line[0]_\$line[1]_LOG.txt\\n\"}' $annotation_path/Regions_250_Genes.txt > $annotation_path/jobs/createExonRegionDatabase.sh";

print COMMANDS "\n\n#Then submit to cluster:";
if ($cluster_commands){
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}";
  print COMMANDS "\ncd $annotation_path/jobs/";
  print COMMANDS "\nmqsub --file createExonRegionDatabase.sh --name $annotation_conf{ALEXA_SEQ_DB}_CERD  --mkdir  --delay 5";
}else{
  print COMMANDS "\n$annotation_path/jobs/createExonRegionDatabase.sh";
}
print COMMANDS "\n\n#Then check success rate";
print COMMANDS "\ngrep -h \"SCRIPT COMPLETE\" $annotation_path/logs/createExonRegionDatabase/* | wc -l";
print COMMANDS "\ngrep -v Chromosome $annotation_path/Regions_250_Genes.txt | wc -l";
print COMMANDS "\ncd $annotation_path/exonRegions/";
print COMMANDS "\ngrep -h -v \"Gene_ID\" temp/*.txt > tmp.txt";
print COMMANDS "\nhead -q -n 1 temp/*.txt | sort | uniq | cat - tmp.txt > exonRegions_annotated.txt";
print COMMANDS "\nrm -f tmp.txt";
print COMMANDS "\ngzip exonRegions_annotated.txt";
print COMMANDS "\nrm -f temp/*.txt";
print COMMANDS "\ncd blastdb/";
print COMMANDS "\ncat temp/exonRegions*.fa > exonRegions.fa";
print COMMANDS "\n$alexa_seq_conf{BLAST_BIN_DIR}/formatdb -i exonRegions.fa -t exonRegions -p F -o T -n exonRegions";
print COMMANDS "\ngzip exonRegions.fa";
print COMMANDS "\nrm -f temp/*.fa";
if ($annotation_conf{ALIGNMENT_OPTION} == 2 || $annotation_conf{ALIGNMENT_OPTION} == 3){
  print COMMANDS "\n\n#Now create a bwa database";
  print COMMANDS "\ncp exonRegions.fa.gz ../bwadb/";
  print COMMANDS "\ncd ../bwadb/";
  print COMMANDS "\n$alexa_seq_conf{BWA_BIN_DIR}/bwa index exonRegions.fa.gz";
}

#10.) Create a database for Exon Junctions 
print MAGENTA, "\n\n10.) Create a database for Exon Junctions", RESET;
print COMMANDS "\n\n\n#10.) Create a database for Exon Junctions:";
print COMMANDS "\n\n#First extract the junctions";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/createExonJunctionDatabase.pl  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --target_length=$annotation_conf{JUNCTION_SEQ_LENGTH}  --outdir=$annotation_path/exonJunctions/  --logfile=$annotation_path/logs/createExonJunctionDatabase/createExonJunctionDatabase_LOG.txt";
print COMMANDS "\n\n#Then create jobs annotate them";
print COMMANDS "\ncd $annotation_path/";
print COMMANDS "\nrm -f exonJunctions/temp/*.txt";
print COMMANDS "\nrm -f logs/annotateExonJunctions/*.txt";
print COMMANDS "\nperl -ne 'unless(\$_ =~ /Chromosome/){chomp(\$_); \@line=split(\"\\t\", \$_); print \"source $alexa_seq_conf{SHELL_FILE}; $script_dir/alternativeExpressionDatabase/annotateExonJunctions.pl  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --chr_filter=\\\"\$line[0]:\$line[1]:\$line[2]-\$line[3]\\\"  --junction_file=$annotation_path/exonJunctions/exonJunctions_"."$annotation_conf{JUNCTION_SEQ_LENGTH}mers.txt  --ucsc_align_dir=$annotation_path/mrna_est/partitions/   --genbank_mapfile=$annotation_path/mrna_est/GenBankToOrganism.btree  --wiggle=3  --outfile=$annotation_path/exonJunctions/temp/exonJunctions_annotated_\$line[0]_\$line[1].txt  --logfile=$annotation_path/logs/annotateExonJunctions/annotateExonJunctions_\$line[0]_\$line[1]_LOG.txt\\n\"}' $annotation_path/Regions_250_Genes.txt > $annotation_path/jobs/annotateExonJunctions.sh";

print COMMANDS "\n\n#Then submit these jobs to the cluster";
if ($cluster_commands){
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}";
  print COMMANDS "\ncd $annotation_path/jobs/";
  print COMMANDS "\nmqsub --file annotateExonJunctions.sh --name $annotation_conf{ALEXA_SEQ_DB}_AEJ  --mkdir  --delay 5";
}else{
  print COMMANDS "\n$annotation_path/jobs/annotateExonJunctions.sh";
}
print COMMANDS "\n\n#Then check success rate";
print COMMANDS "\ngrep -h \"SCRIPT COMPLETE\" $annotation_path/logs/annotateExonJunctions/* | wc -l";
print COMMANDS "\ngrep -v Chromosome $annotation_path/Regions_250_Genes.txt | wc -l";
print COMMANDS "\n\n#Now join the annotation and fasta files together and create a blast database";
print COMMANDS "\ncd $annotation_path/exonJunctions/";
print COMMANDS "\ngrep -h -v Junction_ID temp/*.txt | sort -n > tmp.txt";
print COMMANDS "\nhead -q -n 1 temp/*.txt | sort | uniq | cat - tmp.txt > exonJunctions_"."$annotation_conf{JUNCTION_SEQ_LENGTH}mers_annotated.txt";
print COMMANDS "\nrm -f tmp.txt";
print COMMANDS "\ngzip exonJunctions_"."$annotation_conf{JUNCTION_SEQ_LENGTH}mers_annotated.txt";
print COMMANDS "\ncp exonJunctions_"."$annotation_conf{JUNCTION_SEQ_LENGTH}mers.fa blastdb/";
print COMMANDS "\ncd blastdb/";
print COMMANDS "\n$alexa_seq_conf{BLAST_BIN_DIR}/formatdb -i exonJunctions_"."$annotation_conf{JUNCTION_SEQ_LENGTH}mers.fa -t exonJunctions_"."$annotation_conf{JUNCTION_SEQ_LENGTH}mers -p F -o F -n exonJunctions_"."$annotation_conf{JUNCTION_SEQ_LENGTH}mers";
print COMMANDS "\ngzip exonJunctions_"."$annotation_conf{JUNCTION_SEQ_LENGTH}mers.fa";
if ($annotation_conf{ALIGNMENT_OPTION} == 2 || $annotation_conf{ALIGNMENT_OPTION} == 3){
  print COMMANDS "\n\n#Now create a bwa database";
  print COMMANDS "\ncp exonJunctions_"."$annotation_conf{JUNCTION_SEQ_LENGTH}mers.fa.gz ../bwadb/";
  print COMMANDS "\ncd ../bwadb/";
  print COMMANDS "\n$alexa_seq_conf{BWA_BIN_DIR}/bwa index exonJunctions_"."$annotation_conf{JUNCTION_SEQ_LENGTH}mers.fa.gz";
}
print COMMANDS "\n\n#If the above step seems to have worked then clean-up temp files...";
print COMMANDS "\ncd $annotation_path/exonJunctions/";
print COMMANDS "\nrm -f exonJunctions_"."$annotation_conf{JUNCTION_SEQ_LENGTH}mers.txt";
print COMMANDS "\nrm -f temp/*.txt";
print COMMANDS "\nrm -f exonJunctions_"."$annotation_conf{JUNCTION_SEQ_LENGTH}mers.fa";

#11.) Create a database for Exon Boundaries
print MAGENTA, "\n\n11.) Create a database for Exon Boundaries", RESET;
print COMMANDS "\n\n\n#11.) Create a database for Exon Boundaries:";
print COMMANDS "\n\n#First extract the boundaries";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/createExonBoundaryDatabase.pl  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --target_length=$annotation_conf{BOUNDARY_SEQ_LENGTH}  --outdir=$annotation_path/exonBoundaries/  --logfile=$annotation_path/logs/createExonBoundaryDatabase/createExonBoundaryDatabase_LOG.txt";
print COMMANDS "\n\n#Then create jobs annotate them";
print COMMANDS "\ncd $annotation_path/";
print COMMANDS "\nrm -f exonBoundaries/temp/*.txt";
print COMMANDS "\nrm -f logs/annotateExonBoundaries/*.txt";
print COMMANDS "\nperl -ne 'unless(\$_ =~ /Chromosome/){chomp(\$_); \@line=split(\"\\t\", \$_); print \"source $alexa_seq_conf{SHELL_FILE}; $script_dir/alternativeExpressionDatabase/annotateExonBoundaries.pl  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --chr_filter=\\\"\$line[0]:\$line[1]:\$line[2]-\$line[3]\\\"  --target_size=$annotation_conf{BOUNDARY_SEQ_LENGTH}  --boundary_file=$annotation_path/exonBoundaries/exonBoundaries_"."$annotation_conf{BOUNDARY_SEQ_LENGTH}mers.txt  --ucsc_align_dir=$annotation_path/mrna_est/partitions/   --genbank_mapfile=$annotation_path/mrna_est/GenBankToOrganism.btree  --wiggle=3  --outfile=$annotation_path/exonBoundaries/temp/exonBoundaries_annotated_\$line[0]_\$line[1].txt  --logfile=$annotation_path/logs/annotateExonBoundaries/annotateExonBoundaries_\$line[0]_\$line[1]_LOG.txt\\n\"}' $annotation_path/Regions_250_Genes.txt > $annotation_path/jobs/annotateExonBoundaries.sh";

print COMMANDS "\n\n#Then submit these jobs to the cluster";
if ($cluster_commands){
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}";
  print COMMANDS "\ncd $annotation_path/jobs/";
  print COMMANDS "\nmqsub --file annotateExonBoundaries.sh --name $annotation_conf{ALEXA_SEQ_DB}_AEB  --mkdir  --delay 5";
}else{
  print COMMANDS "\n$annotation_path/jobs/annotateExonBoundaries.sh";
}
print COMMANDS "\n\n#Then check success rate";
print COMMANDS "\ngrep -h \"SCRIPT COMPLETE\" $annotation_path/logs/annotateExonBoundaries/* | wc -l";
print COMMANDS "\ngrep -v Chromosome $annotation_path/Regions_250_Genes.txt | wc -l";
print COMMANDS "\n\n#Now join the annotation and fasta files together and create a blast database";
print COMMANDS "\ncd $annotation_path/exonBoundaries/";
print COMMANDS "\ngrep -h -v Boundary_ID temp/*.txt | sort -n > tmp.txt";
print COMMANDS "\nhead -q -n 1 temp/*.txt | sort | uniq | cat - tmp.txt > exonBoundaries_"."$annotation_conf{BOUNDARY_SEQ_LENGTH}mers_annotated.txt";
print COMMANDS "\nrm -f tmp.txt";
print COMMANDS "\ngzip exonBoundaries_"."$annotation_conf{BOUNDARY_SEQ_LENGTH}mers_annotated.txt";
print COMMANDS "\ncp exonBoundaries_"."$annotation_conf{BOUNDARY_SEQ_LENGTH}mers.fa blastdb/";
print COMMANDS "\ncd blastdb/";
print COMMANDS "\n$alexa_seq_conf{BLAST_BIN_DIR}/formatdb -i exonBoundaries_"."$annotation_conf{BOUNDARY_SEQ_LENGTH}mers.fa -t exonBoundaries_"."$annotation_conf{BOUNDARY_SEQ_LENGTH}mers -p F -o F -n exonBoundaries_"."$annotation_conf{BOUNDARY_SEQ_LENGTH}mers";
print COMMANDS "\ngzip exonBoundaries_"."$annotation_conf{BOUNDARY_SEQ_LENGTH}mers.fa";
if ($annotation_conf{ALIGNMENT_OPTION} == 2 || $annotation_conf{ALIGNMENT_OPTION} == 3){
  print COMMANDS "\n\n#Now create a bwa database";
  print COMMANDS "\ncp exonBoundaries_"."$annotation_conf{BOUNDARY_SEQ_LENGTH}mers.fa.gz ../bwadb/";
  print COMMANDS "\ncd ../bwadb/";
  print COMMANDS "\n$alexa_seq_conf{BWA_BIN_DIR}/bwa index exonBoundaries_"."$annotation_conf{BOUNDARY_SEQ_LENGTH}mers.fa.gz";
}
print COMMANDS "\n\n#If the above step seems to have worked then clean-up temp files...";
print COMMANDS "\ncd $annotation_path/exonBoundaries/";
print COMMANDS "\nrm -f exonBoundaries_"."$annotation_conf{BOUNDARY_SEQ_LENGTH}mers.txt";
print COMMANDS "\nrm -f temp/*.txt";
print COMMANDS "\nrm -f exonBoundaries_"."$annotation_conf{BOUNDARY_SEQ_LENGTH}mers.fa";

#12.) Create a database of all introns
print MAGENTA, "\n\n12.) Create a database of all introns", RESET;
print COMMANDS "\n\n\n#12.) Create a database of all introns:";
print COMMANDS "\n\n#First create the jobs";
print COMMANDS "\nperl -ne 'chomp (\$_); \@line =split(\"\\t\", \$_); unless (\$_ =~ /Chromosome/){print \"source $alexa_seq_conf{SHELL_FILE}; $script_dir/alternativeExpressionDatabase/createIntronDatabase.pl  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --min_size=42  --chr_filter=\\\"\$line[0]:\$line[1]:\$line[2]-\$line[3]\\\"  --ucsc_align_dir=$annotation_path/mrna_est/partitions/  --genbank_mapfile=$annotation_path/mrna_est/GenBankToOrganism.btree  --wiggle=5  --outdir=$annotation_path/introns/temp/  --logfile=$annotation_path/logs/createIntronDatabase/createIntronDatabase_\$line[0]_\$line[1]_LOG.txt\\n\"}' $annotation_path/Regions_50_Genes.txt > $annotation_path/jobs/createIntronDatabase.sh";
print COMMANDS "\nrm -f $annotation_path/logs/createIntronDatabase/*.txt";
print COMMANDS "\nrm -f $annotation_path/introns/temp/*.fa";
print COMMANDS "\nrm -f $annotation_path/introns/temp/*.txt";
print COMMANDS "\n\n#Then submit these jobs to the cluster";
if ($cluster_commands){
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}";
  print COMMANDS "\ncd $annotation_path/jobs/";
  print COMMANDS "\nmqsub  --file createIntronDatabase.sh  --name $annotation_conf{ALEXA_SEQ_DB}_CID  --mkdir  --delay 15";
}else{
  print COMMANDS "\n$annotation_path/jobs/createIntronDatabase.sh";
}
print COMMANDS "\n\n#Then check success rate";
print COMMANDS "\ngrep -h \"SCRIPT COMPLETE\" $annotation_path/logs/createIntronDatabase/* | wc -l";
print COMMANDS "\ngrep -v Chromosome $annotation_path/Regions_250_Genes.txt | wc -l";
print COMMANDS "\ngrep -v Chromosome $annotation_path/Regions_50_Genes.txt | wc -l";
print COMMANDS "\n\n#Retrieve jobs that need to be rerun - get list of jobs that completed according to their log files - and then filter the batch file with these";
print COMMANDS "\ncd $annotation_path/";
print COMMANDS "\ngrep \"SCRIPT COMPLETE\" logs/createIntronDatabase/* | perl -ne 'if (\$_ =~ /(create.*\\.txt)/){print \"\$1\\n\"}' > logs/tmp.txt";
print COMMANDS "\ngrep -v -f logs/tmp.txt jobs/createIntronDatabase.sh > jobs/createIntronDatabase_REPAIR.sh";
print COMMANDS "\nrm -f logs/tmp.txt";
print COMMANDS "\n\n#Now join the annotation and fasta files together and create a blast database";
print COMMANDS "\ncd $annotation_path/introns/";
print COMMANDS "\ngrep -h -v \"Intron_ID\" temp/*introns_annotated.txt > tmp.txt";
print COMMANDS "\nhead -q -n 1 temp/*_introns_annotated.txt | sort | uniq | cat - tmp.txt > introns_annotated.txt";
print COMMANDS "\ngzip introns_annotated.txt";
print COMMANDS "\ngrep -h -v \"Active_Region_ID\" temp/*intronRegionsActive.txt > tmp.txt";
print COMMANDS "\nhead -q -n 1 temp/*_intronRegionsActive.txt | sort| uniq | cat - tmp.txt > activeIntronRegions.txt";
print COMMANDS "\ngzip activeIntronRegions.txt";
print COMMANDS "\ngrep -h -v \"Silent_Region_ID\" temp/*intronRegionsSilent.txt > tmp.txt";
print COMMANDS "\nhead -q -n 1 temp/*_intronRegionsSilent.txt | sort | uniq | cat - tmp.txt > silentIntronRegions.txt";
print COMMANDS "\ngzip silentIntronRegions.txt";
print COMMANDS "\ncat temp/*.fa > blastdb/introns.fa";
print COMMANDS "\nrm -f tmp.txt";
print COMMANDS "\ncd blastdb/";
print COMMANDS "\n$alexa_seq_conf{BLAST_BIN_DIR}/formatdb -i introns.fa -t introns -p F -o F -n introns";
print COMMANDS "\ngzip introns.fa";

if ($annotation_conf{ALIGNMENT_OPTION} == 2 || $annotation_conf{ALIGNMENT_OPTION} == 3){
  print COMMANDS "\n\n#Now create a bwa database";
  print COMMANDS "\ncp introns.fa.gz ../bwadb/";
  print COMMANDS "\ncd ../bwadb/";
  print COMMANDS "\n$alexa_seq_conf{BWA_BIN_DIR}/bwa index introns.fa.gz";
}
print COMMANDS "\n\n#If the above step seems to have worked then clean-up temp files...";
print COMMANDS "\ncd $annotation_path/introns/";
print COMMANDS "\nrm -f temp/*.fa";
print COMMANDS "\nrm -f temp/*_annotated.txt";
print COMMANDS "\nrm -f temp/*_intronRegionsActive.txt";
print COMMANDS "\nrm -f temp/*_intronRegionsSilent.txt";


#13.) Create a database of all intergenic regions
print MAGENTA, "\n\n13.) Create a database of all intergenic regions", RESET;
print COMMANDS "\n\n\n#13.) Create a database of all intergenic regions:";
print COMMANDS "\n\n#First create the jobs";
print COMMANDS "\nperl -ne 'chomp (\$_); \@line=split(\"\\t\", \$_); unless (\$_ =~ /Chromosome/){print \"source $alexa_seq_conf{SHELL_FILE}; $script_dir/alternativeExpressionDatabase/createIntergenicDatabase.pl  --ensembl_api_version=$annotation_conf{ENSEMBL_VERSION}  --species=\\\"$annotation_conf{SPECIES_NAME}\\\"  --ensembl_database=$ensembl_db_name  --ensembl_server=$alexa_seq_conf{ENSEMBL_SERVER}  --ensembl_user=$alexa_seq_conf{ENSEMBL_USER}  --ensembl_password=$alexa_seq_conf{ENSEMBL_PASSWORD}  --alexa_database=$alexa_db_name  --alexa_server=$alexa_seq_conf{ALEXA_SERVER}  --alexa_user=$alexa_seq_conf{ALEXA_USER1}  --alexa_password=$alexa_seq_conf{ALEXA_PASSWORD1}  --min_size=42  --chr_filter=\\\"\$line[0]:\$line[1]:\$line[2]-\$line[3]\\\"  --ucsc_align_dir=$annotation_path/mrna_est/partitions/  --genbank_mapfile=$annotation_path/mrna_est/GenBankToOrganism.btree  --wiggle=5  --outdir=$annotation_path/intergenics/temp/  --logfile=$annotation_path//logs/createIntergenicDatabase/createIntergenicDatabase_\$line[0]_\$line[1]_LOG.txt\\n\"}' $annotation_path/Regions_50_Genes.txt > $annotation_path/jobs/createIntergenicDatabase.sh";
print COMMANDS "\n\nrm -f $annotation_path/logs/createIntergenicDatabase/*.txt";
print COMMANDS "\nrm -f $annotation_path/intergenics/temp/*.fa";
print COMMANDS "\nrm -f $annotation_path/intergenics/temp/*Active.txt";
print COMMANDS "\nrm -f $annotation_path/intergenics/temp/*Silent.txt";
print COMMANDS "\nrm -f $annotation_path/intergenics/temp/*.txt";
print COMMANDS "\n\n#Then submit these jobs to the cluster";
if ($cluster_commands){
  print COMMANDS "\nssh $alexa_seq_conf{CLUSTER_HEAD_NODE}";
  print COMMANDS "\ncd $annotation_path/jobs/";
  print COMMANDS "\nmqsub  --file createIntergenicDatabase.sh  --name $annotation_conf{ALEXA_SEQ_DB}_CIGD  --mkdir  --delay 5";
}else{
  print COMMANDS "\n$annotation_path/jobs/createIntergenicDatabase.sh";
}
print COMMANDS "\n\n#Then check success rate";
print COMMANDS "\ngrep -h \"SCRIPT COMPLETE\" $annotation_path/logs/createIntergenicDatabase/* | wc -l";
print COMMANDS "\ngrep -v Chromosome $annotation_path/Regions_250_Genes.txt | wc -l";
print COMMANDS "\ngrep -v Chromosome $annotation_path/Regions_50_Genes.txt | wc -l";
print COMMANDS "\n\n#Retrieve jobs that need to be rerun - get list of jobs that completed according to their log files - and then filter the batch file with these";
print COMMANDS "\ncd $annotation_path/";
print COMMANDS "\ngrep \"SCRIPT COMPLETE\" logs/createIntergenicDatabase/* | perl -ne 'if (\$_ =~ /(create.*\\.txt)/){print \"\$1\\n\"}' > logs/tmp.txt";
print COMMANDS "\ngrep -v -f logs/tmp.txt jobs/createIntergenicDatabase.sh > jobs/createIntergenicDatabase_REPAIR.sh";
print COMMANDS "\nrm -f logs/tmp.txt";
print COMMANDS "\n\n#Now join the annotation and fasta files together and create a blast database";
print COMMANDS "\ncd $annotation_path/intergenics/";
print COMMANDS "\ngrep -h -v \"Intergenic_ID\" temp/*intergenics_annotated.txt > tmp.txt";
print COMMANDS "\nhead -q -n 1 temp/*_annotated.txt | sort | uniq | cat - tmp.txt > intergenics_annotated.txt";
print COMMANDS "\ngzip intergenics_annotated.txt";
print COMMANDS "\ngrep -h -v \"Active_Region_ID\" temp/*intergenicRegionsActive.txt > tmp.txt";
print COMMANDS "\nhead -q -n 1 temp/*_intergenicRegionsActive.txt | sort | uniq | cat - tmp.txt > activeIntergenicRegions.txt";
print COMMANDS "\ngzip activeIntergenicRegions.txt";
print COMMANDS "\ngrep -h -v \"Silent_Region_ID\" temp/*intergenicRegionsSilent.txt > tmp.txt";
print COMMANDS "\nhead -q -n 1 temp/*_intergenicRegionsSilent.txt | sort | uniq | cat - tmp.txt > silentIntergenicRegions.txt";
print COMMANDS "\ngzip silentIntergenicRegions.txt";
print COMMANDS "\ncat temp/*.fa > blastdb/intergenics.fa";
print COMMANDS "\nrm -f tmp.txt";
print COMMANDS "\ncd blastdb/";
print COMMANDS "\n$alexa_seq_conf{BLAST_BIN_DIR}/formatdb -i intergenics.fa -t intergenics -p F -o F -n intergenics";
print COMMANDS "\ngzip intergenics.fa";
if ($annotation_conf{ALIGNMENT_OPTION} == 2 || $annotation_conf{ALIGNMENT_OPTION} == 3){
  print COMMANDS "\n\n#Now create a bwa database";
  print COMMANDS "\ncp intergenics.fa.gz ../bwadb/";
  print COMMANDS "\ncd ../bwadb/";
  print COMMANDS "\n$alexa_seq_conf{BWA_BIN_DIR}/bwa index intergenics.fa.gz";
}
print COMMANDS "\n\n#If the above step seems to have worked then clean-up temp files...";
print COMMANDS "\ncd $annotation_path/intergenics/";
print COMMANDS "\nrm -f temp/*.fa";
print COMMANDS "\nrm -f temp/*_annotated.txt";
print COMMANDS "\nrm -f temp/*_intergenicRegionsActive.txt";
print COMMANDS "\nrm -f temp/*_intergenicRegionsSilent.txt";


#14.) Create a database for transcript specific annotation records
print MAGENTA, "\n\n14.) Create a database for transcript specific annotation records", RESET;
print COMMANDS "\n\n\n#14.) Create a database for transcript specific annotation records:";
print COMMANDS "\n\n$script_dir/alternativeExpressionDatabase/createTranscriptDatabase.pl  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --annotation_dir=$annotation_path/  --junction_seq_size=$annotation_conf{JUNCTION_SEQ_LENGTH}  --outdir=$annotation_path/transcripts/";
print COMMANDS "\n\ngzip $annotation_path/transcripts/transcripts_annotated.txt";
print COMMANDS "\n\n#Create a blast database for all full-length transcript sequences";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/createEnsemblTranscriptFasta.pl  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER1}  --password=$alexa_seq_conf{ALEXA_PASSWORD1}  --outfile=$annotation_path/transcripts/blastdb/transcripts.fa  --logfile=$annotation_path/logs/createEnsemblTranscriptsFasta_LOG.txt";
print COMMANDS "\n\ncd $annotation_path/transcripts/blastdb/";
print COMMANDS "\n$alexa_seq_conf{BLAST_BIN_DIR}/formatdb -i transcripts.fa -t transcripts -p F -o F -n transcripts";
print COMMANDS "\ngzip $annotation_path/transcripts/blastdb/transcripts.fa";
if ($annotation_conf{ALIGNMENT_OPTION} == 2 || $annotation_conf{ALIGNMENT_OPTION} == 3){
  print COMMANDS "\n\n#Now create a bwa database";
  print COMMANDS "\ncp transcripts.fa.gz ../bwadb/";
  print COMMANDS "\ncd ../bwadb/";
  print COMMANDS "\n$alexa_seq_conf{BWA_BIN_DIR}/bwa index transcripts.fa.gz";
}

#15.) For all feature types, append a unique FID (an id that is unique across all feature types)
print MAGENTA, "\n\n15.) For all feature types, append a unique FID (an id that is unique across all feature types)", RESET;
print COMMANDS "\n\n\n#15.) For all feature types, append a unique FID (an id that is unique across all feature types):";
print COMMANDS "\n#*** All previous steps should be done first ! ***";
print COMMANDS "\n#Use a prefix for each feature type:";
print COMMANDS "\n# G = genes\n# T = transcripts\n# ER = exon regions\n# EJ = exon junctions\n# EB = exon boundaries\n# IN = introns\n# AIN = active intron region\n# SIN = silent intron region\n# IG = intergenic region\n# AIG = active intergenic region\n# SIG = silent intergenic region\n";
print COMMANDS "\ncd $annotation_path/";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/appendFeatureId.pl  --infile=genes/genes_annotated.txt.gz  --outfile=genes/tmp.txt --prefix='G'";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/appendFeatureId.pl  --infile=transcripts/transcripts_annotated.txt.gz  --outfile=transcripts/tmp.txt --prefix='T'";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/appendFeatureId.pl  --infile=exonRegions/exonRegions_annotated.txt.gz  --outfile=exonRegions/tmp.txt --prefix='ER'";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/appendFeatureId.pl  --infile=exonJunctions/exonJunctions_"."$annotation_conf{JUNCTION_SEQ_LENGTH}mers_annotated.txt.gz  --outfile=exonJunctions/tmp.txt --prefix='EJ'";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/appendFeatureId.pl  --infile=exonBoundaries/exonBoundaries_"."$annotation_conf{BOUNDARY_SEQ_LENGTH}mers_annotated.txt.gz  --outfile=exonBoundaries/tmp.txt --prefix='EB'";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/appendFeatureId.pl  --infile=introns/introns_annotated.txt.gz  --outfile=introns/tmp.txt --prefix='IN'";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/appendFeatureId.pl  --infile=introns/activeIntronRegions.txt.gz  --outfile=introns/tmp.txt --prefix='AIN'";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/appendFeatureId.pl  --infile=introns/silentIntronRegions.txt.gz  --outfile=introns/tmp.txt --prefix='SIN'";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/appendFeatureId.pl  --infile=intergenics/intergenics_annotated.txt.gz  --outfile=intergenics/tmp.txt --prefix='IG'";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/appendFeatureId.pl  --infile=intergenics/activeIntergenicRegions.txt.gz  --outfile=intergenics/tmp.txt --prefix='AIG'";
print COMMANDS "\n$script_dir/alternativeExpressionDatabase/appendFeatureId.pl  --infile=intergenics/silentIntergenicRegions.txt.gz  --outfile=intergenics/tmp.txt --prefix='SIG'";


#16.) Create a backup of the ALEXA database
print MAGENTA, "\n\n16.) Create a backup of the ALEXA database", RESET;
print COMMANDS "\n\n\n#16.) Create a backup of the ALEXA database:";
print COMMANDS "\n\nssh $alexa_seq_conf{ALEXA_SERVER}";
print COMMANDS "\n$script_dir/sql/backupAlexaDb.pl  --database=$alexa_db_name  --server=$alexa_seq_conf{ALEXA_SERVER}  --user=$alexa_seq_conf{ALEXA_USER2}  --password=$alexa_seq_conf{ALEXA_PASSWORD2}  --mysql_dump_bin=/usr/bin/mysqldump  --working_dir=$annotation_path/alexa_db/  --logfile=$annotation_path/logs/backAlexaDb_LOG.txt";


#17.) Create a backup of the ALEXA database
print MAGENTA, "\n\n17.) Create a summary statistics file for the build", RESET;
print COMMANDS "\n\n\n#17.) Create a summary statistics file for the build:";
print COMMANDS "\n#Summarize disk usage of the entire database package (compressed)";
print COMMANDS "\n#Summarize the grand total number of features";
print COMMANDS "\n#Summarize the number of each type of feature";
print COMMANDS "\n#If applicable, summarize the number of features with EnsEMBL support, mRNA support, EST support, conservation support";
print COMMANDS "\n\n$script_dir/bash/summarizeAnnotationDatabase.sh $annotation_conf{ALEXA_SEQ_DB} $annotation_conf{ENSEMBL_VERSION} $annotation_conf{JUNCTION_SEQ_LENGTH} $alexa_seq_conf{SEQUENCE_DB_DIR}";


#18.) Package and compress DB files
print MAGENTA, "\n\n18.) Package and compress DB files", RESET;
print COMMANDS "\n\n\n#18.) Package and compress DB files:";
print COMMANDS "\ncd $alexa_seq_conf{SEQUENCE_DB_DIR}/";
print COMMANDS "\ntar -cf $annotation_conf{ALEXA_SEQ_DB}.tar $annotation_conf{ALEXA_SEQ_DB}/*.txt $annotation_conf{ALEXA_SEQ_DB}/*.gene_info $annotation_conf{ALEXA_SEQ_DB}/alexa_db $annotation_conf{ALEXA_SEQ_DB}/exonBoundaries/ $annotation_conf{ALEXA_SEQ_DB}/exonJunctions/ $annotation_conf{ALEXA_SEQ_DB}/exonRegions/ $annotation_conf{ALEXA_SEQ_DB}/genes/ $annotation_conf{ALEXA_SEQ_DB}/intergenics/ $annotation_conf{ALEXA_SEQ_DB}/introns/ $annotation_conf{ALEXA_SEQ_DB}/transcripts/ $annotation_conf{ALEXA_SEQ_DB}/repeats/";
print COMMANDS "\ngzip --best $annotation_conf{ALEXA_SEQ_DB}.tar";

print YELLOW, "\n\nAll necessary commands to create an ALEXA-Seq annotation DB for '$annotation_conf{ALEXA_SEQ_DB}' have been written to following commands file:\n$commands_file\n\nFollow the instructions in this file to proceed...\n\n", RESET;

print COMMANDS "\n\n";
print "\n\n";
close(COMMANDS);

exit();
