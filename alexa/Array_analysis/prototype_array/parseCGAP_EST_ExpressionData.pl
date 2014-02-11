#!/usr/local/bin/perl58 -w
#Written by Malachi Griffith
#This takes a series of input files containing library, gene and expression info for various CGAP EST libraries
#It will then create an output file where each line contains:
#EntrezGene/LocusID, ExpressionCount in Condition1 (ex Cerebellum), ExpressionCount in Condition2 (ex Prostate cancer)
#All neccessary files come from ftp://ftp1.nci.nih.gov/pub/CGAP/

#To do this will involve the following steps
#1.) Import Unigene/CGAP library info
#2.) Import Gene info.  Get Unigene to Entrez Gene maps
#3.) Import expression counts for each gene for each library
#4.) Summarize the number of counts for the libraries matching certain criteria

use strict;
use Data::Dumper;
use Getopt::Long;

#Initialize command line options
my $library_data = '';
my $gene_data = '';
my $expression_data = '';
my $out_file = '';

GetOptions ('library_data=s'=>\$library_data,'gene_data=s'=>\$gene_data, 'expression_data=s'=>\$expression_data,
	    'out_file=s'=>\$out_file);

#Provide instruction to the user
print "\n\nUsage:";
print "\n\tSpecify the file containing CGAP library data : --library_data";
print "\n\tSpecify the file containing CGAP/Unigene gene data using: --gene_data";
print "\n\tSpecify the file containing EST Expression data using: --expression_data";
print "\n\tSpecify the output file using: --out_file";
print "\n\nExample: parseCGAP_EST_ExpressionData.pl --library_data  --exon_probes=filteredExonProbes.txt --negative_probes=filteredNegativeControlProbes.txt --nimblegen_data=All_pair.txt --out=dataSummary.txt\n\n";

#Make sure all options were specified
unless ($junction_probes && $exon_probes && $negative_probes && $nimblegen_data_file && $out_file){
  print "\nOptions missing!\n\n";
  exit();
}
