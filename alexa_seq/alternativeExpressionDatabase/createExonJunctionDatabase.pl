#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to generate all possible exon-exon junction sequences for each ensembl gene
#The user can specify a probe length (which must be an even number for simplicity)
#Sequences will be centered on the junction

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

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
my $target_length = '';
my $verbose = '';
my $outdir = '';
my $logfile = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'target_length=i'=>\$target_length, 'verbose=s'=>\$verbose, 'outdir=s'=>\$outdir,
	    'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the target probe length (must be a multiple of 2) using: --target_length (eg. 66)", RESET;
print GREEN, "\n\t\tIf you want verbose output, use: --verbose=yes", RESET;
print GREEN, "\n\tSpecify a path for output files using: --outdir", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of the design process using: --logfile", RESET;
print GREEN, "\n\nExample: createExonJunctionDatabase.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --target_length=62   --outdir=/projects/malachig/sequence_databases/hs_53_36o/exonJunctions/  --logfile=/projects/malachig/sequence_databases/hs_53_36o/logs/createExonJunctionDatabase/createExonJunctionDatabase_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $target_length && $outdir && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

my $current_junction_id = 0;

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#Print out the parameters supplied by the user to the logfile for future reference
print LOG "\nUser Specified the following options:\ndatabase = $database\ntarget_length = $target_length\noutdir = $outdir\nlogfile = $logfile\n\n";

my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};

my $target_gene_count = @gene_ids;
print BLUE, "\nFound $target_gene_count genes\n\n", RESET;
print LOG "\nFound $target_gene_count genes\n\n";

my $gene_count = 0;
my $total_successful_junctions = 0;
my $total_possible_junctions = 0;

#Check output dir and set output file names
$outdir = &checkDir('-dir'=>$outdir, '-clear'=>"no");
my $db_file = "$outdir"."exonJunctions_"."$target_length"."mers.txt";
my $fasta_file = "$outdir"."exonJunctions_"."$target_length"."mers.fa";

#Open the junction output files
print BLUE, "\nAll data will be written to $db_file and $fasta_file\n\n", RESET;
print LOG "\nAll data will be written to $db_file and $fasta_file\n\n";
open (DB, ">$db_file") || die "\nCould not open database output file: $db_file";
open (FASTA, ">$fasta_file") || die "\nCould not open database output file: $fasta_file";

#Print out the header line for the output file
print DB "Junction_ID\tGene_ID\tEnsEMBL_Gene_ID\tGene_Name\tChromosome\tStrand\tUnit1_start\tUnit1_end\tUnit2_start\tUnit2_end\tUnit1_start_chr\tUnit1_end_chr\tUnit2_start_chr\tUnit2_end_chr\tBase_Count\tUnMasked_Base_Count\tCoding_Base_Count\n";

#Get the gene sequence and other gene info for all genes
print BLUE, "\nGetting gene sequence data", RESET;
print LOG "\nGetting gene sequence data";
my $storable_name = "$database"."_AllGenes_GeneInfo_WithSeq.storable";
my $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"yes", '-storable'=>$storable_name);

#Get the transcripts for these genes to allow comparison to actual known transcripts
$storable_name = "$database"."_AllGenes_TranscriptInfo_WithSeq.storable";
my $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"yes", '-storable'=>$storable_name);

my $masked_gene_ref = &getMaskedGene ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);


#Build transcript sequences for all genes
foreach my $gene_id (keys %{$gene_transcripts_ref}){
  my $trans_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};
  my %trans_seqs;
  foreach my $trans_id (keys %{$trans_ref}){
    my $exons_ref = $trans_ref->{$trans_id}->{exons};

    my $seq = '';
    foreach my $exon_id (sort {$exons_ref->{$a}->{exon_start} <=> $exons_ref->{$b}->{exon_start}} keys %{$exons_ref}){
      $seq = "$seq"."$exons_ref->{$exon_id}->{sequence}";
    }
    $trans_seqs{$trans_id} = $seq;
  }
  $gene_transcripts_ref->{$gene_id}->{trans_seqs} = \%trans_seqs;
}


#Get the exons for all genes
print BLUE, "\nGetting exon coordinate data", RESET;
print LOG "\nGetting exon coordinate data";
my $gene_exons_ref = &getExons ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"yes");

#Get the total theoretical junctions possible for all genes
print BLUE, "\nGetting theoretical junction counts data", RESET;
print LOG "\nGetting theoretical junction counts data";
my $junction_counts_ref = &junctionProbeCombinations('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);

foreach my $gene_id (@gene_ids){
  $gene_count++;

  my $successful_junctions = 0;

  print CYAN, "\n\n\n*****************************************************************************", RESET;
  print CYAN, "\nGENE: $gene_count\tALEXA: $gene_id\tENSG: $genes_ref->{$gene_id}->{ensembl_g_id}", RESET;
  print LOG "\n\n\n*****************************************************************************";
  print LOG "\nGENE: $gene_count\tALEXA: $gene_id\tENSG: $genes_ref->{$gene_id}->{ensembl_g_id}";

  my $gene_seq = $genes_ref->{$gene_id}->{sequence};
  my $trans_seq_ref = $gene_transcripts_ref->{$gene_id}->{trans_seqs};

  my $masked_gene_seq = $masked_gene_ref->{$gene_id}->{sequence};
  my $protein_bases_ref = &getProteinBases ('-gene_id'=>$gene_id, '-genes_ref'=>$genes_ref, '-gene_transcripts_ref'=>$gene_transcripts_ref);


  #Get all exons for this gene
  my $exons_ref = $gene_exons_ref->{$gene_id}->{exons};
  my $exon_count = keys %{$exons_ref};

  #Go through each exon and find the total number of unique end positions
  #Keep track of the exon IDs associated with each start and end positions (store as an array to keep track of redundant start/stop positions)
  my %start_positions;
  my %end_positions;
  foreach my $exon_id (sort keys %{$exons_ref}){

    #Unique start positions and associated exon IDs
    my $exon_start = $exons_ref->{$exon_id}->{exon_start};
    if ($start_positions{$exon_start}{exon_ids}){
      my @tmp = @{$start_positions{$exon_start}{exon_ids}};
      push (@tmp, $exon_id);
      $start_positions{$exon_start}{exon_ids} = \@tmp;
    }else{
      my @tmp;
      push (@tmp, $exon_id);
      $start_positions{$exon_start}{exon_ids} = \@tmp;
    }

    #Unique end positions and associated exon IDs
    my $exon_end = $exons_ref->{$exon_id}->{exon_end};
    if ($end_positions{$exon_end}{exon_ids}){
      my @tmp = @{$end_positions{$exon_end}{exon_ids}};
      push (@tmp, $exon_id);
      $end_positions{$exon_end}{exon_ids} = \@tmp;
    }else{
      my @tmp;
      push (@tmp, $exon_id);
      $end_positions{$exon_end}{exon_ids} = \@tmp;
    }
  }

  #If the gene has only one exon - skip
  if ($exon_count == 1){
    print CYAN, "\n\tGene has only a single exon - skipping", RESET;
    next();
  }

  #NOTE: Some exons have the same start or end position!
  #Junctions should be tracked primarily with respect to their gene coordinates as specified in the GeneProbe table
  #We will also attempt to track which exons are associated with each junction

  #Find every valid NON-Redundant exon-exon pair, consider every possible permutation and check for valid connections by start/end position
  #i.e. Find the number of valid connections between each end position and all start positions
  foreach my $exon1_end_pos (sort {$a <=> $b} keys %end_positions){

    #For each end position, compare to all unique start positions of exons, any that are greater represent valid connections
    foreach my $exon2_start_pos (sort {$a <=> $b} keys %start_positions){

      #Consider whether this combination of exon1 and exon2 is valid
      if ($exon2_start_pos > $exon1_end_pos){
	my %junctions;
	my $junction_count = 0;

	#Determine the start and end position of the junction halves based on the length of the junction required
	#First 1/2 of junction
	my $unit1_start_pos = ($exon1_end_pos-($target_length/2));
	my $unit1_length = ($target_length/2);
	my $unit1_end_pos = $exon1_end_pos;
	my $unit1_seq = substr($gene_seq, $unit1_start_pos, $unit1_length);
	my $seq1_length = length($unit1_seq);

        unless($unit1_seq){
          my $test_length = length($gene_seq);
          print RED, "\nCould not retrieve sequence 1 (Gene: $gene_id $genes_ref->{$gene_id}->{gene_name}) with substr(SEQ, $unit1_start_pos, $unit1_length).  SEQ lenth is $test_length\n\n", RESET;
          exit();
        }

	#Second 1/2 of junction
	my $unit2_start_pos = $exon2_start_pos;
	my $unit2_length = ($target_length/2);
	my $unit2_end_pos = ($unit2_start_pos + $unit2_length);
	my $unit2_seq = substr($gene_seq, $unit2_start_pos-1, $unit2_length);

        unless($unit2_seq){
          my $test_length = length($gene_seq);
          print RED, "\nCould not retrieve sequence 2 (Gene: $gene_id $genes_ref->{$gene_id}->{gene_name}) with substr(SEQ, $unit2_start_pos-1, $unit2_length).  SEQ lenth is $test_length\n\n", RESET;
          exit();
        }

	my $seq2_length = length($unit2_seq);
	
	my $junction_seq = "$unit1_seq"."$unit2_seq";
	my $junction_seq_length = length($junction_seq);

	#Still need to consider the length of each exon.  There may be several exons that end/start at a particular position
	#But not all of them are neccessarily long enough to allow a junction to be designed ...
	#This length issue is a problem which is difficult to fix completely...
	
	#The exon_end_pos and exon_start_pos combination used for this junction are associated with at least two exons
	#Get all the exons involved in this specific junction and check if the resulting unit_1_start position is within the range of the exon1 options
	#Similarly, check if the resulting unit_2_end position is within the range of the exon2 possibilities
	#If the junction exceeds the boundaries of all possible exons at each side it must be rejected, otherwise, the exons for which
	#it is valid will be noted
	my @exon1_ids = @{$end_positions{$exon1_end_pos}{exon_ids}};
	my @exon2_ids = @{$start_positions{$exon2_start_pos}{exon_ids}};

	my @valid_exon1_ids;
	foreach my $exon1_id (@exon1_ids){
	  if ($unit1_start_pos <= $exons_ref->{$exon1_id}->{exon_start}){
	    if ($verbose eq "yes"){
	      print YELLOW, "\n\tJunction exceeds 5' boundary of exon: $exon1_id\tJUNCTION = $junction_seq", RESET;
	    }
	    print LOG "\n\tJunction exceeds 5' boundary of exon: $exon1_id\tJUNCTION = $junction_seq";
	  }else{
	    push (@valid_exon1_ids, $exon1_id);
	  }
	}
	my @valid_exon2_ids;
	foreach my $exon2_id (@exon2_ids){
	  if ($unit2_end_pos >= $exons_ref->{$exon2_id}->{exon_end}){
	    if ($verbose eq "yes"){
	      print YELLOW, "\n\tJunction exceeds 3' boundary of exon: $exon2_id\tJUNCTION = $junction_seq", RESET;
	    }
	      print LOG "\n\tJunction exceeds 3' boundary of exon: $exon2_id\tJUNCTION = $junction_seq";
	  }else{
	    push (@valid_exon2_ids, $exon2_id);
	  }
	}
	my $valid_exon1_ids = @valid_exon1_ids;
	my $valid_exon2_ids = @valid_exon2_ids;



	#Deal with junctions that fail due to short exons on one or both ends
	#Make sure the junction is valid for at least one of the pairs of exons with redundant start/stop positions
	#If one or both ends of the junction sequence exceeded an exon boundary try building the junction via a different method
	#When you have these small exons, the combinatorial possibilities become ridiculuous
	#TRY TO CAPTURE AT LEAST ANY OF THE JUNCTIONS THAT ARE KNOWN TO OCCUR IN ENSEMBL TRANSCRIPTS!!
	#Basically the exon sequence will be expanded into flanking exons
	#BUT the coordinates will still be reported for only the part involving directly adjacent exons
	# - so the coordinates will indicate this the junction is really shorter than it is...

	if ($valid_exon1_ids == 0 || $valid_exon2_ids == 0){

	  #Find the exons from this gene that have: $exon1_end_pos
	  #Find the exons from this gene that have: $exon2_start_pos
	  my @exon1_list;
	  foreach my $exon_id (sort keys %{$exons_ref}){
	    if ($exon1_end_pos == $exons_ref->{$exon_id}->{exon_end}){
	      push (@exon1_list, $exon_id);
	    }
	  }
	  my @exon2_list;
	  foreach my $exon_id (sort keys %{$exons_ref}){
	    if ($exon2_start_pos == $exons_ref->{$exon_id}->{exon_start}){
	      push (@exon2_list, $exon_id);
	    }
	  }

	  #Build the pairwise possible sequence combinations for this junction
	  my %seqs;
	  my $flank = ($target_length+10)/2;
	  foreach my $e1_id (@exon1_list){
	    foreach my $e2_id (@exon2_list){
	      my $sequence = "$exons_ref->{$e1_id}->{sequence}"."$exons_ref->{$e2_id}->{sequence}";
	      my $e1_length = length($exons_ref->{$e1_id}->{sequence});
	      my $e2_length = length($exons_ref->{$e2_id}->{sequence});

	      #Go through the known transcript sequences of this gene and see if these combos actual occur in real transcripts
	      foreach my $trans_id (keys %{$trans_seq_ref}){
		my $trans_seq = $trans_seq_ref->{$trans_id};

		if ($trans_seq =~ /(\w{$flank}$sequence\w{$flank})/){

		  my $seq_with_flank = $1;

		  my $start = ($flank + $e1_length)-($target_length/2);
		  my $length = $target_length;
		  my $new_junction = substr($seq_with_flank, $start, $length);

		  $seqs{$new_junction}{tmp} = '';
		}
	      }
	    }
	  }

	  #Deal with the outcome of this attempt
	  my $new_junctions_found = keys %seqs;
	  if ($new_junctions_found == 0){
	    #If none of these junctions matched known transcripts, skip
	    next();
	  }elsif($new_junctions_found == 1){
	    #A single replacement known junction was found!  Replace the current junction seq with this seq
	    foreach my $nj (keys %seqs){
	      $junction_seq = $nj;
	      if ($verbose eq "yes"){
		print YELLOW, "\n\tFound a single known new junction:\n\t$junction_seq", RESET;
	      }
	    }
	  }else{
	    #Deal with multiple, known, new junctions found??
	    if ($verbose eq "yes"){
	      print YELLOW, "\n\tFound MULTIPLE known new junction:\n", RESET;
	      print Dumper %seqs;
	    }
	    next();
	  }

	  if(($unit1_start_pos < 0) || ($unit2_end_pos >= length($gene_seq))){
	    print RED, "\n\tUnexpected junction coords!!\n\n", RESET;
	    print RED, "\nUnit1_start=$unit1_start_pos\tUnit1_end=$unit1_end_pos", RESET;
	    print RED, "\nUnit2_start=$unit2_start_pos\tUnit2_end=$unit2_end_pos", RESET;
	    next();
	  }

	}

	#Sanity check of resulting junction length
	my $result_length = length($junction_seq);
	unless ($result_length == $target_length){
	  print RED, "\n\tUnexpected junction length ($result_length)!!\n\n", RESET;
	  print RED, "\nJUNCTION: $junction_seq", RESET;
	  print RED, "\nUNIT1: $unit1_seq\tUNIT2: $unit2_seq", RESET;
	  print RED, "\nExon1: @exon1_ids\tExon2: @exon2_ids", RESET;
	  print RED, "\nUnit1_start=$unit1_start_pos\tUnit1_end=$unit1_end_pos", RESET;
	  print RED, "\nUnit2_start=$unit2_start_pos\tUnit2_end=$unit2_end_pos", RESET;
	  close (LOG);
	  $alexa_dbh->disconnect();
	  exit();
	}
	
	#NOTE: Some exons actually contain N's from the underlying genomic sequence in ensembl!  
	#For simplicity, junctions that incorporate these unknown bases should be skipped!

	#Check for presence of genomic N's and other non-valid letters such as Ambiguiety codes
	if($junction_seq =~ /[^ATCG]/){
	  if ($verbose eq "yes"){
	    print YELLOW, "\n\tThis junction: $junction_seq contains an invalid bases such as an N or R", RESET;
	  }
	  print LOG "\n\tThis junction: $junction_seq contains an invalid bases such as an N or R";
	  next();
	}

	$junction_count++;

	$junctions{$junction_count}{sequence} = $junction_seq;
	$junctions{$junction_count}{junction_seq_length} = $junction_seq_length;
	$junctions{$junction_count}{exon1_end} = $exon1_end_pos;
	$junctions{$junction_count}{exon2_start} = $exon2_start_pos;
	$junctions{$junction_count}{unit1_start} = $unit1_start_pos;
	$junctions{$junction_count}{unit1_end} = $unit1_end_pos;
	$junctions{$junction_count}{unit1_seq} = $unit1_seq;
	$junctions{$junction_count}{unit1_length} = $unit1_length;
	$junctions{$junction_count}{unit2_start} = $unit2_start_pos;
	$junctions{$junction_count}{unit2_end} = $unit2_end_pos;
	$junctions{$junction_count}{unit2_seq} = $unit2_seq;
	$junctions{$junction_count}{unit2_length} = $unit2_length;

        #Calculate chromosome coordinates
        my $coords1_ref = &convertGeneCoordinates('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$unit1_start_pos, '-end_pos'=>$unit1_end_pos, '-ordered'=>"yes");
        my $coords2_ref = &convertGeneCoordinates('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$unit2_start_pos, '-end_pos'=>$unit2_end_pos, '-ordered'=>"yes");

        #Make sure all 4 coords are sorted from smallest to largest, regardless of strand from which they came 
        my @temp = ($coords1_ref->{$gene_id}->{chr_start}, $coords1_ref->{$gene_id}->{chr_end}, $coords2_ref->{$gene_id}->{chr_start}, $coords2_ref->{$gene_id}->{chr_end});
        my @temp_sorted = sort {$a <=> $b} @temp;

        #Get the base_count, unmasked_base_count, and protein_coding count of each junction

        #Masked Seq
        my $start1 = $junctions{$junction_count}{unit1_start};
        my $start2 = $junctions{$junction_count}{unit2_start};
        my $size = $target_length/2;
        my $masked_exon_region_seq1 = substr ($masked_gene_seq, $start1, $size);
        my $masked_exon_region_seq2 = substr ($masked_gene_seq, $start2, $size);
        my $masked_exon_region_seq = "$masked_exon_region_seq1"."$masked_exon_region_seq2";
        my $masked_n_count = ($masked_exon_region_seq =~ tr/N/N/);
        my $unmasked_base_count = $target_length - $masked_n_count;
  
        #Protein coding bases
        my $chr_start1 = $temp_sorted[0];
        my $chr_end1 = $temp_sorted[1];
        my $chr_start2 = $temp_sorted[2];
        my $chr_end2 = $temp_sorted[3];
        my $coding_base_count1 = &getProteinBaseCount ('-protein_bases_ref'=>$protein_bases_ref, '-chr_start'=>$chr_start1, '-chr_end'=>$chr_end1-1);
        my $coding_base_count2 = &getProteinBaseCount ('-protein_bases_ref'=>$protein_bases_ref, '-chr_start'=>$chr_start2, '-chr_end'=>$chr_end2-1);
        my $coding_base_count = $coding_base_count1+$coding_base_count2;

	#Print junction info to Output junction file
	$successful_junctions++;
	$total_successful_junctions++;
	$current_junction_id++;

	my $j_count = $junction_count;

	#Junction_ID\tGene_ID\tEnsEMBL_Gene_ID\tSequence\tChromosome\tStrand\tUnit1_start\tUnit1_end\tUnit2_start\tUnit2_end\tUnit1_start_chr\tUnit1_end_chr\tUnit2_start_chr\tUnit2_end_chr
	print DB "$current_junction_id\t$gene_id\t$genes_ref->{$gene_id}->{ensembl_g_id}\t$genes_ref->{$gene_id}->{gene_name}\t$genes_ref->{$gene_id}->{chromosome}\t$genes_ref->{$gene_id}->{chr_strand}\t$junctions{$j_count}{unit1_start}\t$junctions{$j_count}{unit1_end}\t$junctions{$j_count}{unit2_start}\t$junctions{$j_count}{unit2_end}\t$temp_sorted[0]\t$temp_sorted[1]\t$temp_sorted[2]\t$temp_sorted[3]\t$target_length\t$unmasked_base_count\t$coding_base_count\n";

        print FASTA ">$current_junction_id\n$junctions{$j_count}{sequence}\n";

      }#Valid junction loop
    }#Junction start position loop (unit_2_start)
  }#Junction end position loop (unit_1_end)

  my $possible_junctions = $junction_counts_ref->{$gene_id}->{exon_exon};

  $total_possible_junctions += $possible_junctions;

  print CYAN, "\n\nSUMMARY for Gene ID: $gene_id", RESET;
  print CYAN, "\nNumber Exons = $exon_count\tPossible Exon-Exon Junctions = $possible_junctions\tSuccessful Junctions = $successful_junctions", RESET;
  print LOG "\n\nSUMMARY for Gene ID: $gene_id";
  print LOG "\nNumber Exons = $exon_count\tPossible Exon-Exon Junctions = $possible_junctions\tSuccessful Junctions = $successful_junctions";

}#Gene loop

print CYAN, "\n\nTotal Possible Junctions = $total_possible_junctions\nTotal Successful Junctions = $total_successful_junctions\n\n", RESET;
print LOG "\n\nTotal Possible Junctions = $total_possible_junctions\nTotal Successful Junctions = $total_successful_junctions\n\n";

close (DB);
close (FASTA);
close (LOG);

#Close database connection
$alexa_dbh->disconnect();

exit();

