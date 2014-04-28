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
my $junction_db = '';
my $junction_seq_file = '';
my $outfile = '';
my $ensembl_version = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'junction_seq_file=s'=>\$junction_seq_file, 'junction_db=s'=>\$junction_db, 'outfile=s'=>\$outfile, 'ensembl_version=s'=>\$ensembl_version);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the junction DB to analyze using:  --junction_db", RESET;
print GREEN, "\n\tSpecify the junction sequence file using:  --junction_seq_file", RESET;
print GREEN, "\n\tSpecify an output file using: --outfile", RESET;
print GREEN, "\n\tSpecify the ensembl version using:  --ensembl_version", RESET;
print GREEN, "\n\nExample: identifyJunctionFrames.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --junction_db=exonJunctions_60mers_annotated.txt.gz  --junction_seq_file=exonJunctions_60mers.fa.gz  --outfile=exonJunctions_60mers_annotated.txt2  --ensembl_version=53\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $junction_db && $junction_seq_file && $outfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Load the specified ensembl API and BioPerl
&loadEnsemblApi('-api'=>$ensembl_version);
require Bio::Perl;

my $mem_message;
$|=1;

my $debug_limit = 1000;

#Get basic gene info
print BLUE, "\n\nGetting basic gene info:", RESET;
my $genes_ref;
my $gene_transcripts_ref;
my @grand_orfs;
&getBasicGeneInfo();
$mem_message = &memoryUsage();
print YELLOW, "\n$mem_message", RESET;

#Grab the actual sequences for the junctions
print BLUE, "\n\nRetrieving junction sequences:", RESET;
my $j_seq_ref = &getJunctionSeqs('-j_seq_file'=>$junction_seq_file);
$mem_message = &memoryUsage();
print YELLOW, "\n$mem_message", RESET;

#Now go through the annotation file
print BLUE, "\n\nSearch for reading frames and maintenance of known frame:\n", RESET;
open (ANN, "zcat $junction_db |") || die "\n\nCould not open annotation file: $junction_db\n\n";
open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
my $header = 1;
my $c = 0;
while(<ANN>){
  $c+=2;
  if ($c >= 200000){
    print BLUE, ".", RESET;
    $c = 0;
  }
  chomp($_);
  my @line = split("\t", $_);
  if ($header == 1){
    print OUT "$_\tFrame\tFrameMaintained\n";
    $header = 0;

    #If the file already has frame info, abort
    if ($_ =~ /Frame/){
      print YELLOW, "\n\nFile already has frame columns - aborting", RESET;
      close(ANN);
      close(OUT);
      system ("rm -f $outfile");
      exit();
    }
    next();
  }
  my $id = $line[0];
  my $gene_id = $line[1];
  my $seq = $j_seq_ref->{$id}->{seq};

  #Translate the feature sequence
  my $seq_obj = Bio::Seq->new('-seq'=>$seq, '-desc'=>'Sample Bio::Seq object', '-display_id' => 'something', '-accession_number' => 'accnum', '-alphabet' => 'dna');
  my @seqs = Bio::SeqUtils->translate_3frames($seq_obj);

  my @frames = ($seqs[0]->seq(), $seqs[1]->seq(), $seqs[2]->seq());

  #Figure out which frame is correct for this feature.  Do this by matching a test sequence from within it to the ORFs of the corresponding gene
  #First get the test peptide for each junction
  my $inset = 1;
  my @subframes1;
  my @subframes2;
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

  #Now test these test peptides against the gene's ORFs to look for a match
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
  my $frame_maintained = "NA";
  my $peptide_match;
  if ($matching_frame1_count >= 1){
    $frame_maintained = "no";
    if ($matching_frame1_count == 1 && $matching_frame2_count == 1 && $matching_frame1 == $matching_frame2){
      $frame_maintained = "yes";
    }
    $peptide_match = $frames[$matching_frame1];
  }else{
    $matching_frame1="NA";
    $peptide_match = "NA";
  }
  print OUT "$_\t$matching_frame1\t$frame_maintained\n";
  #print "$_\t$matching_frame1\t$frame_maintained\t$peptide_match\n";
}

close(ANN);
close(OUT);

$|=0;

my $cmd_gzip = "gzip -f $outfile";
print BLUE, "\n\nCompressing output file:\n$cmd_gzip\n\n", RESET;
system($cmd_gzip);


exit();


#######################################################################################################################################################################
#Get basic gene info from ALEXA DB                                                                                                                                    #
#######################################################################################################################################################################
sub getBasicGeneInfo{

  my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
  my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};
  my $storable_name = "$database"."_AllGenes_GeneInfo_WithSeq.storable";
  $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"yes", '-storable'=>$storable_name);

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

  return();
}


#######################################################################################################################################################################
#Grab the actual junction sequences                                                                                                                                   #
#######################################################################################################################################################################
sub getJunctionSeqs{
  my %args = @_;
  my $file = $args{'-j_seq_file'};

  my %j; 
  open(SEQ, "zcat $file |") || die "\nCould not open junction seq file file: $file\n\n";
  my $match_found = 0;
  my $id;
  my $c = 0;
  while(<SEQ>){
    $c++;
#    if ($c >= $debug_limit){
#      last();
#    }

    chomp($_);
    my $line = $_;
    if ($_ =~ /\>(.*)/){
      $id = $1;
      $match_found = 1;
      next();
    }
    if ($match_found){
      $j{$id}{seq} = $_;
      $match_found = 0;
    }
  }
  close(SEQ);
  
  return(\%j);
}



