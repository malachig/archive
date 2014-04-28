#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Get all ensembl exon coords.  These are not quite the same as 'exon regions' that are used throughout the ALEXA-Seq pipeline
use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use BerkeleyDB;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $outfile = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'outfile=s'=>\$outfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the output file using:  --outfile", RESET;
print GREEN, "\n\nExample: getExonCoords.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --outfile=ensembl_exons_hs_53_36o.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $outfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};
my $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");
my $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

#Process each gene
foreach my $gene_id (@gene_ids){

  if ($genes_ref->{$gene_id}->{chromosome} eq "MT"){$genes_ref->{$gene_id}->{chromosome} = "M";}

  #Get a list of non-redundant exon for this gene
  my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};
  my $strand = $genes_ref->{$gene_id}->{chr_strand};
  my $chr = $genes_ref->{$gene_id}->{chromosome};
  my $gene_name = $genes_ref->{$gene_id}->{gene_name};
  my $ensg = $genes_ref->{$gene_id}->{ensembl_g_id};

  if ($strand eq "1"){
    $strand = "+";
  }else{
    $strand = "-";
  }
  my %reference_exons;
  foreach my $trans_id (keys %{$transcripts_ref}){
    my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

    foreach my $trans_exon_id (keys %{$exons_ref}){
      my $trans_exon_start = $exons_ref->{$trans_exon_id}->{exon_start};
      my $trans_exon_end = $exons_ref->{$trans_exon_id}->{exon_end};

      #Check each of the reference exons to see if one of these is the same, otherwise add it to the list
      my $redundant_exon = 0;
      foreach my $ref_exon_id (keys %reference_exons){
	if ($trans_exon_start == $reference_exons{$ref_exon_id}{start} && $trans_exon_end == $reference_exons{$ref_exon_id}{end}){
	  $redundant_exon = 1;
	}
      }
      #Unless the current transcript exon was found to be redundant, add it to the list
      unless ($redundant_exon == 1){
        my $coords_ref = &convertGeneCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$trans_exon_start, '-end_pos'=>$trans_exon_end, '-ordered'=>"yes");
	$reference_exons{$trans_exon_id}{start} = $coords_ref->{$gene_id}->{chr_start};
	$reference_exons{$trans_exon_id}{end} = $coords_ref->{$gene_id}->{chr_end};
      }
    }
  }

  foreach my $e (sort {$reference_exons{$a}->{start} <=> $reference_exons{$a}->{end}} keys %reference_exons){
    my $start = $reference_exons{$e}{start};
    my $end = $reference_exons{$e}{end};
    print OUT "chr$chr:$start-$end\t$strand\t$gene_name\t$ensg\n";
  }
}
$alexa_dbh->disconnect();
close(OUT);
exit();

