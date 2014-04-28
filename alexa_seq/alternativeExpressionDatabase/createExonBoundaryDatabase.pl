#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to identify all possible exon-intron boundary sequences for each ensembl gene
#The user can specify a sequence length, which must be an even number for simplicity. The resulting probes will be centered on the boundary.

#Note: Intron-exon boundaries can potentially be used to detect intron inclusion events.  For overlapping exons, the case is not so simple.
#In this case the 'Intron-exon' probe is actually spanning an alternate splice site at the 3' or 5' end of the exon.
#In fact any time expression seems to be present for an intron-exon probe it may be an indication that an alternate splice site is being used
#In other words, it doesn't neccessarily mean that an intron was included, it may also mean that the exon is simply longer than expected because
#an alternate 5' or 3' splice site was used.

#Sequences will be centered on the boundary

#To allow for testing of this script, the user can also specify a single target Ensembl gene ID or design probes for all Ensembl Genes

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
	    'target_length=i'=>\$target_length, 'verbose=s'=>\$verbose, 'outdir=s'=>\$outdir, 'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the target boundary sequence length (must be a multiple of 2) using: --target_length (eg. 62)", RESET;
print GREEN, "\n\t\tIf you want verbose output, use: --verbose=yes", RESET;
print GREEN, "\n\tSpecify a path for output files using: --outdir", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of the design process using: --logfile", RESET;
print GREEN, "\n\nExample: createExonBoundaryDatabase.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --target_length=62  --outdir=/projects/malachig/sequence_databases/hs_53_36o/exonBoundaries/  --logfile=/projects/malachig/sequence_databases/hs_53_36o/logs/createExonBoundaryDatabase/createExonBoundaryDatabase_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $target_length && $outdir && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#Print out the parameters supplied by the user to the logfile for future reference
print LOG "\nUser Specified the following options:\ndatabase = $database\ntarget_length = $target_length\noutdir = $outdir\nlogfile = $logfile\n\n";

my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};
my $gene_count = 0;
my $total_successful_boundaries = 0;
my $total_possible_boundaries = 0;
my $current_boundary_id = 0;

#Check output dir and set output file names
$outdir = &checkDir('-dir'=>$outdir, '-clear'=>"no");
my $db_file = "$outdir"."exonBoundaries_"."$target_length"."mers.txt";
my $fasta_file = "$outdir"."exonBoundaries_"."$target_length"."mers.fa";


#Open the boundary output files
print BLUE, "\nAll data will be written to $db_file and $fasta_file\n\n", RESET;
print LOG "\nAll data will be written to $db_file and $fasta_file\n\n";
open (DB, ">$db_file") || die "\nCould not open database output file: $db_file";
open (FASTA, ">$fasta_file") || die "\nCould not open database output file: $fasta_file";


#Print out the header line for the output file
#The probe type will be Inton-Exon or Exon-Intron
print DB "Boundary_ID\tGene_ID\tEnsEMBL_Gene_ID\tGene_Name\tChromosome\tStrand\tBoundary_Type\tUnit1_start\tUnit1_end\tUnit1_start_chr\tUnit1_end_chr\tBase_Count\tUnMasked_Base_Count\tCoding_Base_Count\n";

#Get the gene sequence and other gene info for all genes
print BLUE, "\nGetting gene sequence data", RESET;
print LOG "\nGetting gene sequence data";
my $storable_gene = "$database"."_AllGenes_GeneInfo_WithSeq.storable";
my $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids,'-sequence'=>"yes", '-storable'=>$storable_gene);

#Get the transcripts for these genes to allow comparison to actual known transcripts
my $storable_name = "$database"."_AllGenes_TranscriptInfo_NoSeq.storable";
my $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no", '-storable'=>$storable_name);

my $masked_gene_ref = &getMaskedGene ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);



#Get the exons for all genes
print BLUE, "\nGetting exon coordinate data", RESET;
print LOG "\nGetting exon coordinate data";
my $storable_exon = "$database"."_AllExons_WithSeq.storable";
my $gene_exons_ref = &getExons ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"yes", '-storable'=>$storable_exon);

#Get the total theoretical probes possible for all genes
print BLUE, "\nGetting theoretical boundary counts data", RESET;
print LOG "\nGetting theoretical boundary counts data";
my $boundary_counts_ref = &junctionProbeCombinations('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);

#Process genes
foreach my $gene_id (@gene_ids){
  $gene_count++;

  #Keep track of exons, successful and unsuccessful boundaries for this gene
  my $successful_boundaries = 0;

  #Get the start and end coordinates of the entire gene
  my $gene_start = $genes_ref->{$gene_id}->{gene_start};
  my $gene_end = $genes_ref->{$gene_id}->{gene_end};
  my $ensembl_g_id = $genes_ref->{$gene_id}->{ensembl_g_id};
  my $gene_name = $genes_ref->{$gene_id}->{gene_name};
  my $chromosome = $genes_ref->{$gene_id}->{chromosome};
  my $strand = $genes_ref->{$gene_id}->{chr_strand};

  my $masked_gene_seq = $masked_gene_ref->{$gene_id}->{sequence};
  my $protein_bases_ref = &getProteinBases ('-gene_id'=>$gene_id, '-genes_ref'=>$genes_ref, '-gene_transcripts_ref'=>$gene_transcripts_ref);

  print CYAN, "\n\n\n**************************************************************************************************************", RESET;
  print CYAN, "\nGENE: $gene_count\tALEXA: $gene_id\tENSG: $genes_ref->{$gene_id}->{ensembl_g_id}\tGene Coords: ($gene_start - $gene_end)", RESET;
  print LOG "\n\n\n**************************************************************************************************************";
  print LOG "\nGENE: $gene_count\tALEXA: $gene_id\tENSG: $genes_ref->{$gene_id}->{ensembl_g_id}\tGene Coords: ($gene_start - $gene_end)";

  my $gene_seq = $genes_ref->{$gene_id}->{sequence};

  #1.) Get all exons for this gene
  my $exons_ref = $gene_exons_ref->{$gene_id}->{exons};
  my $exon_count = keys %{$exons_ref};

  #2.) Create a hash to track which exons correspond to which exons
  my %junction_maps_5p;
  my %junction_maps_3p;
  foreach my $exon_id (sort keys %{$exons_ref}){
    my $start = $exons_ref->{$exon_id}->{exon_start};
    my $end = $exons_ref->{$exon_id}->{exon_end};

    #Exon start positions (5p intron-exon junctions)
    if ($junction_maps_5p{$start}{exon_ids}){
      my @tmp = @{$junction_maps_5p{$start}{exon_ids}};
      push (@tmp, $exon_id);
      $junction_maps_5p{$start}{exon_ids} = \@tmp;
    }else{
      my @tmp;
      push (@tmp, $exon_id);
      $junction_maps_5p{$start}{exon_ids} = \@tmp;
    }
    #Exon end positions (3p exon-intron junctions)
    if ($junction_maps_3p{$end}{exon_ids}){
      my @tmp = @{$junction_maps_3p{$end}{exon_ids}};
      push (@tmp, $exon_id);
      $junction_maps_3p{$end}{exon_ids} = \@tmp;
    }else{
      my @tmp;
      push (@tmp, $exon_id);
      $junction_maps_3p{$end}{exon_ids} = \@tmp;
    }
  }

  my @exon_array;
  foreach my $exon_id (sort {$exons_ref->{$a}->{exon_start} <=> $exons_ref->{$b}->{exon_start}} keys %{$exons_ref}){
    push (@exon_array, $exon_id);
  }

  #If the gene has only one exon - skip
  if ($exon_count == 1){
    print CYAN, "\n\tGene has only a single exon - skipping", RESET;
    next();
  }

  #4.) Go through each exon sequentially and get the exon-intron junction sequences
  my %boundary_coords_5p;  #Used to keep track of the 5 prime exon ends targeted
  my %boundary_coords_3p;  #Used to keep track of the 3 prime exon ends targeted

  #PROCESS EXONS
  for (my $i = 0; $i < $exon_count; $i++){

    my %boundaries_5p;
    my %boundaries_3p;
    my $boundary_count_5p = 0;
    my $boundary_count_3p = 0;

    my $exon_id = $exon_array[$i];
    if ($verbose eq "yes"){
      print YELLOW, "\n\nProcess EXON: $exon_id\tSTART: $exons_ref->{$exon_id}->{exon_start}\tEND: $exons_ref->{$exon_id}->{exon_end}", RESET;
    }
    print LOG "\n\nProcess EXON: $exon_id\tSTART: $exons_ref->{$exon_id}->{exon_start}\tEND: $exons_ref->{$exon_id}->{exon_end}";

    my $exon_seq = $exons_ref->{$exon_id}->{sequence};

    #If the exon being considered is to short to allow for the specified target length, skip it
    my $exon_length = length($exon_seq);

    if ($exon_length < ($target_length/2)){
      if ($verbose eq "yes"){
	print YELLOW, "\n\t\tExon: $exon_id  Exon is too short ($exon_length bp) to allow a sequence of the specified length", RESET;
      }
      print LOG "\n\t\tExon: $exon_id  Exon is too short ($exon_length bp) to allow a sequence of the specified length";
      next();
    }

    #Define the exon sequence
    my $gene_seq = $genes_ref->{$gene_id}->{sequence};
    my $gene_length = length($gene_seq);
    my $exon_start = $exons_ref->{$exon_id}->{exon_start};
    my $exon_end = $exons_ref->{$exon_id}->{exon_end};

    ######################################################
    #4-A.) Get the 5-PRIME INTRON-EXON boundary sequence #
    ######################################################

    #Avoid processing the first exon (as there is no intron sequence defined upstream of it)
    unless ($exon_start == $gene_start){

      my $intron_start_5p = ($exon_start - ($target_length/2))-1;
      my $intron_end_5p = $exon_start-1;
      my $intron_seq_5p = substr ($gene_seq, $intron_start_5p, $target_length);
      my $boundary_5p_length = length($intron_seq_5p);
      my $test_seq_5p = $intron_seq_5p;

      #If the proposed 5-prime intron start is before the beginning of the gene, skip it (always the case for the first exon)
      if ($intron_start_5p <= 0){
        if ($verbose eq "yes"){
	  print YELLOW, "\n\t\tExon: $exon_id  5-PRIME  Not enough intron sequence 5 prime of this exon", RESET;
	}
	print LOG "\n\t\tExon: $exon_id  5-PRIME  Not enough intron sequence 5 prime of this exon";
	$intron_seq_5p = '';
      }

      #NOTE: Some exons actually contain N's from the underlying genomic sequence in ensembl!
      #For simplicity, boundaries that incorporate these unknown bases should be skipped!
      #Check for presence of genomic N's and other non-valid letters such as Ambiguiety codes
      if($test_seq_5p =~ /[^ATCG]/){
	if ($verbose eq "yes"){
	  print YELLOW, "\n\t\tExon: $exon_id 5-PRIME  Ensembl region for this sequence contains N's or other invalid bases!", RESET;
	}
	print LOG "\n\t\tExon: $exon_id 5-PRIME  Ensembl region for this sequence contains N's or other invalid bases!";
	$intron_seq_5p = '';
      }

      #Keep track of the 5-prime start position used to avoid repeats
      if ($boundary_coords_5p{$intron_start_5p}){
	if ($verbose eq "yes"){
	  print YELLOW, "\n\t\tExon: $exon_id  5-PRIME  Already have this intron-exon junction from a previous exon", RESET;
	}
	print LOG "\n\t\tExon: $exon_id  5-PRIME  Already have this intron-exon junction from a previous exon";
      }elsif($intron_seq_5p){
	#Store successful intron-exon sequence
	$boundary_coords_5p{$intron_start_5p}{junk}='';

	my $exon_start_pos = $exon_start;

	#Sanity check
	unless (length($intron_seq_5p) == $target_length){
	  print RED, "\nExon: $exon_id  5-PRIME  Incorrect sequence length\n\n", RESET;
	  close (LOG);
	  $alexa_dbh->disconnect();
	  exit();
	}

        my $seq_end_5p = ($exon_start+($target_length/2));

	#Store info for this Intron-Exon sequence
        $boundary_count_5p++;
	$boundaries_5p{$boundary_count_5p}{sequence} = $intron_seq_5p;
	$boundaries_5p{$boundary_count_5p}{exon_start} = $exon_start;
	$boundaries_5p{$boundary_count_5p}{exon_end} = $exon_end;
	$boundaries_5p{$boundary_count_5p}{unit1_start} = $intron_start_5p;
	$boundaries_5p{$boundary_count_5p}{unit1_end} = $seq_end_5p;
     
      }
    }

    ########################################################
    #4-B.) Get the 3-PRIME (EXON-INTRON) boundary sequence #
    ########################################################

    #Avoid processing the last exon (as there is no intron sequence defined downstream of it)
    unless ($exon_end == $gene_end){

      my $exon_start_3p = ($exon_end - ($target_length/2)); #Unit1_start
      my $intron_end_3p = $exon_start_3p+$target_length+1;      #Unit2_end
      my $intron_seq_3p = substr ($gene_seq, $exon_start_3p, $target_length);
      my $boundary_3p_length = length($intron_seq_3p);
      my $test_seq_3p = $intron_seq_3p;

      #If the proposed 3-prime intron end is past the end of the gene, skip it (always the case for the last exon)
      if ($intron_end_3p > $gene_length){
	if ($verbose eq "yes"){
	  print YELLOW, "\n\t\tExon: $exon_id  3-PRIME  Not enough intron sequence 3 prime of this exon", RESET;
	}
	print LOG "\n\t\tExon: $exon_id  3-PRIME  Not enough intron sequence 3 prime of this exon";
	$intron_seq_3p = '';
      }

      #NOTE: Some exons actually contain N's from the underlying genomic sequence in ensembl!
      #For simplicity, boundaries that incorporate these unknown bases should be skipped!
      #Check for presence of genomic N's and other non-valid letters such as Ambiguiety codes
      if($test_seq_3p =~ /[^ATCG]/){
	if ($verbose eq "yes"){
	  print YELLOW, "\n\t\tExon: $exon_id 3-PRIME  Ensembl region for this sequence contains N's or other invalid bases!", RESET;
	}
	print LOG "\n\t\tExon: $exon_id 3-PRIME  Ensembl region for this sequence contains N's or other invalid bases!";
	$intron_seq_3p = '';
      }

      #Keep track of the 5-prime and 3-prime start positions used to avoid repeats
      if ($boundary_coords_3p{$exon_start_3p}){
	if ($verbose eq "yes"){
	  print YELLOW, "\n\t\tExon: $exon_id  3-PRIME  Already have this intron-exon junction from a previous exon", RESET;
	}
	print LOG "\n\t\tExon: $exon_id  3-PRIME  Already have this intron-exon junction from a previous exon";

      }elsif($intron_seq_3p){
	#Store successful intron-exon sequence
	$boundary_coords_3p{$exon_start_3p}{junk}='';

	my $exon_end_3p = $exon_end;     #Unit1_end
	my $intron_start_3p = $exon_end+1; #Unit2_start

	#Sanity check
	unless (length($intron_seq_3p) == $target_length){
	  print RED, "\nExon: $exon_id  3-PRIME  Incorrect sequence length\n\n", RESET;
	  close (LOG);
	  $alexa_dbh->disconnect();
	  exit();
	}

	#Store info for this Exon-Intron sequence
        $boundary_count_3p++;

	$boundaries_3p{$boundary_count_3p}{sequence} = $intron_seq_3p;
	$boundaries_3p{$boundary_count_3p}{exon_start} = $exon_start;
	$boundaries_3p{$boundary_count_3p}{exon_end} = $exon_end;
	$boundaries_3p{$boundary_count_3p}{unit1_start} = $exon_start_3p;
	$boundaries_3p{$boundary_count_3p}{unit1_end} = $intron_end_3p;

      }
    }

    #5.) Now that sequences have been generated for both ends of this exon - print them out

    #Intron-Exon
    #First make sure a sequence was actually found
    if ($boundary_count_5p == 1){
      $current_boundary_id++;

      #Print boundary info to output file
      $successful_boundaries++;
      $total_successful_boundaries++;
      my $type = "IE";

      #Get chromosome coordinates for this boundary sequence
      my $coords_ref = &convertGeneCoordinates('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$boundaries_5p{1}{unit1_start}, '-end_pos'=>$boundaries_5p{1}{unit1_end}, '-ordered'=>"yes");

      #Get the base_count, unmasked_base_count, and protein_coding count of each junction
      #Masked Seq
      my $base_count = $target_length;
      my $start = $boundaries_5p{1}{unit1_start};
      my $masked_region_seq = substr ($masked_gene_seq, $start, $target_length);
      my $masked_n_count = ($masked_region_seq =~ tr/N/N/);
      my $unmasked_base_count = $target_length - $masked_n_count;
  
      #Protein coding bases
      my $chr_start = $coords_ref->{$gene_id}->{chr_start};
      my $chr_end = $coords_ref->{$gene_id}->{chr_end};
      my $coding_base_count = &getProteinBaseCount ('-protein_bases_ref'=>$protein_bases_ref, '-chr_start'=>$chr_start+1, '-chr_end'=>$chr_end-1);

      #Print info out to exonBoundary database file: 
      print DB "$current_boundary_id\t$gene_id\t$ensembl_g_id\t$gene_name\t$chromosome\t$strand\t$type\t$boundaries_5p{1}{unit1_start}\t$boundaries_5p{1}{unit1_end}\t$coords_ref->{$gene_id}->{chr_start}\t$coords_ref->{$gene_id}->{chr_end}\t$base_count\t$unmasked_base_count\t$coding_base_count\n";

      #Print ID and sequence to a fasta file
      print FASTA ">$current_boundary_id\n$boundaries_5p{1}{sequence}\n";

    }

    #Exon-Intron
    #First make sure a sequence was actually found
    if ($boundary_count_3p == 1){
      $current_boundary_id++;

      #Print boundary info to output file
      $successful_boundaries++;
      $total_successful_boundaries++;
      my $type = "EI";

      #Get chromosome coordinates for this boundary sequence
      my $coords_ref = &convertGeneCoordinates('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$boundaries_3p{1}{unit1_start}, '-end_pos'=>$boundaries_3p{1}{unit1_end}, '-ordered'=>"yes");

      #Get the base_count, unmasked_base_count, and protein_coding count of each junction
      #Masked Seq
      my $base_count = $target_length;
      my $start = $boundaries_3p{1}{unit1_start};
      my $masked_region_seq = substr ($masked_gene_seq, $start, $target_length);
      my $masked_n_count = ($masked_region_seq =~ tr/N/N/);
      my $unmasked_base_count = $target_length - $masked_n_count;
  
      #Protein coding bases
      my $chr_start = $coords_ref->{$gene_id}->{chr_start};
      my $chr_end = $coords_ref->{$gene_id}->{chr_end};
      my $coding_base_count = &getProteinBaseCount ('-protein_bases_ref'=>$protein_bases_ref, '-chr_start'=>$chr_start+1, '-chr_end'=>$chr_end-1);

      #Print info out to exonBoundary database file: 
      print DB "$current_boundary_id\t$gene_id\t$ensembl_g_id\t$gene_name\t$chromosome\t$strand\t$type\t$boundaries_3p{1}{unit1_start}\t$boundaries_3p{1}{unit1_end}\t$coords_ref->{$gene_id}->{chr_start}\t$coords_ref->{$gene_id}->{chr_end}\t$base_count\t$unmasked_base_count\t$coding_base_count\n";

      #Print ID and sequence to a fasta file
      print FASTA ">$current_boundary_id\n$boundaries_3p{1}{sequence}\n";

    }

  }#Process Exons Loop

  #Calculate the total theoretical boundaries possible for this gene
  my $possible_boundaries = ($boundary_counts_ref->{$gene_id}->{intron_exon});
  $total_possible_boundaries += $possible_boundaries;

  print CYAN, "\n\nSUMMARY for Gene ID: $gene_id", RESET;
  print CYAN, "\nNumber Exons = $exon_count\tPossible Exon-Intron Boundaries = $possible_boundaries\tSuccessful Boundaries = $successful_boundaries", RESET;
  print LOG "\n\nSUMMARY for Gene ID: $gene_id";
  print LOG "\nNumber Exons = $exon_count\tPossible Exon-Intron Boundaries = $possible_boundaries\tSuccessful Boundaries = $successful_boundaries";

}#Gene Loop

print CYAN, "\n\nFinal Summary of all intron junction boundaries examined", RESET;
print CYAN, "\nTotal Possible Boundaries = $total_possible_boundaries\nTotal Successful Boundaries = $total_successful_boundaries\n\n", RESET;
print LOG "\n\nFinal Summary of all intron junction boundaries examined";
print LOG "\nTotal Possible Boundaries = $total_possible_boundaries\nTotal Successful Boundaries = $total_successful_boundaries\n\n";

close (DB);
close(FASTA);
close (LOG);

#Close database connection
$alexa_dbh->disconnect();

exit();
