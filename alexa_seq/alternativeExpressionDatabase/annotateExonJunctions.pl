#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to analyze a file of exon junction sequences and determine for each how many exons are skipped 
#in the theoretical exon-exon join it is attempting to interrogate
#Specifically it will look at each junction record in an input file and get the Exon IDs for Unit1 and Unit2.
#It will then determine the number of exons that are completely skipped by connecting Unit1 and Unit2 (a theoretical splice event)
#To do this, the script will require a hash of all genes and their exons (with coordinates).
#This information will then be added onto the end of each exon junction record.
#By looking at the splicing events currently contained in Ensembl I have determined that 95% of all exon-skipping events involve a maximum
#of 5 skipped exons.


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
use utilities::ALEXA_DB qw(:all);
use utilities::utility qw(:all);
use utilities::ucsc qw(:all);

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $chr_filter = '';
my $junction_file = '';
my $ucsc_align_dir = '';
my $genbank_mapfile = '';
my $wiggle = '';
my $outfile = '';
my $logfile = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'chr_filter=s'=>\$chr_filter, 
            'wiggle=i'=>\$wiggle, 'ucsc_align_dir=s'=>\$ucsc_align_dir, 'genbank_mapfile=s'=>\$genbank_mapfile,
            'junction_file=s'=>\$junction_file, 'outfile=s'=>\$outfile, 'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script determines the number of exons skipped by each exon-exon junction sequence specified in an input file", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password\n", RESET;
print GREEN, "\n\tSpecify a genomic region to analyze using: --chr_filter", RESET;
print GREEN, "\n\tSpecify a tab-delimited input exon-exon junction file using: --junction_file", RESET;
print GREEN, "\n\t\tIt is assumed the the 1st column contains a unique ID", RESET;
print GREEN, "\n\tSpecify a directory containing UCSC mRNA/EST/xmRNA and xEST files using: --ucsc_align_dir", RESET;
print GREEN, "\n\tSpecify a berkeley DB containing GenBank-To-Species mappings using:  --genbank_mapfile", RESET;
print GREEN, "\n\tJunctions for xmRNA and xEST will be allowed to 'wiggle' by the amount specied by:  --wiggle", RESET;
print GREEN, "\n\tSpecify the name of the resulting outfile using: --outfile", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of this script: --logfile", RESET;
print GREEN, "\n\nExample: annotateExonJunctions.pl  --database=ALEXA_hs_53_36o  --server=server_name  --user=username  --password=pwd  --chr_filter='3:16:121020102-126260001'  --junction_file=/projects/malachig/sequence_databases/hs_53_36o/ensembl_exonJunctions_hs_53_36o/exonJunctions_62mers.txt  --ucsc_align_dir=/projects/malachig/sequence_databases/hs_53_36o/mrna_est/partitions/   --genbank_mapfile=/projects/malachig/sequence_databases/hs_53_36o/mrna_est/GenBankToOrganism.btree  --wiggle=3  --outfile=/projects/malachig/sequence_databases/hs_53_36o/ensembl_exonJunctions_hs_53_36o/temp/exonJunctions_annotated_3_16.txt  --logfile=/projects/malachig/sequence_databases/hs_53_36o/logs/annotateExonJunctions/annotateExonJunctions_3_16_LOG.txt\n\n", RESET;

unless ($database && $server && $user && $password && $chr_filter && $junction_file && $ucsc_align_dir && $genbank_mapfile && ($wiggle =~ /\d+/) && $outfile && $logfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

my $region_number;
my $start_filter;
my $end_filter;
if ($chr_filter =~ /(.*)\:(\d+)\:(\d+)\-(\d+)/){
  $chr_filter = $1;
  $region_number = $2;
  $start_filter = $3;
  $end_filter = $4;
  unless ($end_filter > $start_filter){
    print RED, "\nStart of range must be smaller than end ($chr_filter)\n\n", RESET;
    exit();
  }
}else{
  print RED, "\nFormat of chr_filter not understood: $chr_filter (should be of the form:  Y:1:1-9999001)\n\n", RESET;
  exit();
}

#Check the temp dir before proceeding.  Make sure it is empty
$ucsc_align_dir = &checkDir('-dir'=>$ucsc_align_dir, '-clear'=>"no");

#Get the UCSC alignment files (one each for target species mRNA and EST and one each for all other species mRNA and EST)
my $tmp_chr = "chr$chr_filter";
if ($tmp_chr eq "chrMT"){$tmp_chr = "chrM";}

my $ucsc_mrna_table = "$ucsc_align_dir"."$tmp_chr"."_mrna.txt.gz";
my $ucsc_est_table = "$ucsc_align_dir"."$tmp_chr"."_est.txt.gz";
my $ucsc_xmrna_table = "$ucsc_align_dir"."$tmp_chr"."_xmrna.txt.gz";
my $ucsc_xest_table = "$ucsc_align_dir"."$tmp_chr"."_xest.txt.gz";

#Load the GenBank-To-Species mapfile
my %gb_org_map;
tie(%gb_org_map, 'BerkeleyDB::Btree', -Cachesize =>256000000, -Filename=> $genbank_mapfile , -Flags => DB_RDONLY) or die "can't open file $genbank_mapfile: $! $BerkeleyDB::Error\n";


#Global variable storing letters for labeling of Acceptors and Donors
#my @letters = ('aa' .. 'zz');
my @letters = ('a' .. 'z');

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#Print out the parameters supplied by the user to the logfile for future reference
print LOG "\nUser Specified the following options:\ndatabase = $database\njunction_file = $junction_file\noutfile = $outfile\nlogfile = $logfile\n\n";

#1.) First get all the junction IDs from the input file, as well as their Unit1 and Unit2 exon IDs
#Remember that in the case of overlapping exons, there may be multiple unit1 and Unit2 exons
my %junctions;
my %master_gene_list;

print BLUE, "\n1.) Parsing exon-exon junction file: $junction_file for coordinate data ($chr_filter:$region_number:$start_filter-$end_filter)\n\n", RESET;
print LOG "\n1.) Parsing exon-exon junction file: $junction_file for coordinate data ($chr_filter:$region_number:$start_filter-$end_filter)\n\n";

open (JUNCTIONS, "$junction_file") || die "\nCould not open input junction file: $junction_file\n\n";

my $header = 1;
my %columns;
while (<JUNCTIONS>){
  chomp($_);
  my @line = split("\t", $_);

  if ($header == 1){
    my $column_count = 0;

    foreach my $column (@line){
      $columns{$column}{column_pos} = $column_count;
      $column_count++;
    }

    #Check for critical columns and their names
    unless ($columns{'Gene_ID'} && $columns{'Unit1_end'} && $columns{'Unit2_start'}){
      print RED, "\nCritical column missing or named incorrectly in junction file: $junction_file, check input file\n\n", RESET;
      close (LOG);
      exit();
    }

    $header = 0;
    next();
  }

  #Skip this junction unless it is within the target region
  my $chr = $line[$columns{'Chromosome'}{column_pos}];
  my @tmp = ($line[$columns{'Unit1_start_chr'}{column_pos}], $line[$columns{'Unit1_end_chr'}{column_pos}], $line[$columns{'Unit2_start_chr'}{column_pos}], $line[$columns{'Unit2_end_chr'}{column_pos}]);
  my @tmp_sort = sort {$a <=> $b} @tmp;
  my $lower = $tmp_sort[0];
  my $upper = $tmp_sort[3];

  unless ($chr eq $chr_filter && ($lower >= $start_filter && $lower <= $end_filter) && ($upper >= $start_filter && $upper <= $end_filter)){
    next();
  }

  my $junction_id = $line[0];
  $master_gene_list{$line[$columns{'Gene_ID'}{column_pos}]}{tmp} = '';

  $junctions{$junction_id}{gene_id} = $line[$columns{'Gene_ID'}{column_pos}];
  $junctions{$junction_id}{chromosome} = $line[$columns{'Chromosome'}{column_pos}];
  $junctions{$junction_id}{strand} = $line[$columns{'Strand'}{column_pos}];

  $junctions{$junction_id}{unit1_start} = $line[$columns{'Unit1_start'}{column_pos}];
  $junctions{$junction_id}{unit1_end} = $line[$columns{'Unit1_end'}{column_pos}];
  $junctions{$junction_id}{unit2_start} = $line[$columns{'Unit2_start'}{column_pos}];
  $junctions{$junction_id}{unit2_end} = $line[$columns{'Unit2_end'}{column_pos}];

  $junctions{$junction_id}{unit1_start_chr} = $line[$columns{'Unit1_start_chr'}{column_pos}];
  $junctions{$junction_id}{unit1_end_chr} = $line[$columns{'Unit1_end_chr'}{column_pos}];
  $junctions{$junction_id}{unit2_start_chr} = $line[$columns{'Unit2_start_chr'}{column_pos}];
  $junctions{$junction_id}{unit2_end_chr} = $line[$columns{'Unit2_end_chr'}{column_pos}];

  #Format corrections to deal with UCSC
  if ($junctions{$junction_id}{strand} eq "-1"){
    $junctions{$junction_id}{strand} = "-";
  }else{
    $junctions{$junction_id}{strand} = "+";
  }
  $junctions{$junction_id}{chromosome} = "chr"."$junctions{$junction_id}{chromosome}";

}
close (JUNCTIONS);


#2.) Get gene info from ALEXA database
#Only use those genes actually found in the input junction file
my @gene_ids = keys %master_gene_list;

my $gene_count = @gene_ids;
my $total_junction_count = keys %junctions;

print BLUE, "\nFound $total_junction_count exon-exon_junctions for a total of $gene_count genes\n\n", RESET;
print LOG "\nFound $total_junction_count exon-exon junctions for a total of $gene_count genes\n\n";

print BLUE, "\n2.) Extracting Gene and Exon data from ALEXA\n\n", RESET;
print LOG "\n2.) Extracting Gene and Exon data from ALEXA\n\n";

my $genes_ref;
my $gene_transcripts_ref;
my $gene_exon_content_ref;
if ($gene_count > 0){
  my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
  $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");
  $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");
  $gene_exon_content_ref = &getExonContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);
  $alexa_dbh->disconnect();
}else{
  #Some regions may not have any multi-exon genes and therefore no junctions
  #If there are no genes found, it is pointless to continue.
  print "\nSince no junctions were found in this region ...\nSCRIPT COMPLETE\n";
  print LOG "\nSince no junctions were found in this region ...\nSCRIPT COMPLETE\n";
  close(LOG);
  exit();
}

print BLUE, "\n\n2-A.) Generating reference exon objects\n\n", RESET;
print LOG "\n\n2-A.) Generating reference exon objects\n\n";

my $ref_exon_counter = 0;
foreach my $gene_id (@gene_ids){

  my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};
  my $trans_count = keys %{$transcripts_ref};

  $gene_transcripts_ref->{$gene_id}->{trans_count} = $trans_count;

  #2-A.) Assemble a reference set of exons (a superset of all non-redundant exons)
  my %reference_exons;
  foreach my $trans_id (keys %{$transcripts_ref}){
    my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

    foreach my $trans_exon_id (keys %{$exons_ref}){
      my $trans_exon_start = $exons_ref->{$trans_exon_id}->{exon_start};
      my $trans_exon_end = $exons_ref->{$trans_exon_id}->{exon_end};

      #Check each of the reference exons to see if one of these is the same, otherwise add it to the list
      my $redundant_exon = 0;
      foreach my $ref_exon_id (keys %reference_exons){
	if ($trans_exon_start == $reference_exons{$ref_exon_id}{exon_start} && $trans_exon_end == $reference_exons{$ref_exon_id}{exon_end}){
	  $redundant_exon = 1;
	}
      }
      #Unless the current transcript exon was found to be redundant, add it to the list
      unless ($redundant_exon == 1){
        $ref_exon_counter++;
	$reference_exons{$ref_exon_counter}{exon_start} = $trans_exon_start;
	$reference_exons{$ref_exon_counter}{exon_end} = $trans_exon_end;
      }
    }
  }

  #2-B.) Get arrays to represent the reference exons
  my $ref_exon_count = keys %reference_exons;

  my @reference_exon_starts;
  my @reference_exon_ends;
  foreach my $ref_exon_id (sort {$reference_exons{$a}->{exon_start} <=> $reference_exons{$b}->{exon_start}} keys %reference_exons){
    push (@reference_exon_starts, $reference_exons{$ref_exon_id}{exon_start});
    push (@reference_exon_ends, $reference_exons{$ref_exon_id}{exon_end});
  }

  unless($ref_exon_count == 1){
    $gene_transcripts_ref->{$gene_id}->{ref_exon_starts} = \@reference_exon_starts;
    $gene_transcripts_ref->{$gene_id}->{ref_exon_ends} = \@reference_exon_ends;
    $gene_transcripts_ref->{$gene_id}->{ref_exon_count} = $ref_exon_count;
    $gene_transcripts_ref->{$gene_id}->{ref_exons} = \%reference_exons;
  }

  #2-C.) Build a list of non-redundant exon-exon junctions coordinates actually used by transcripts of this gene
  my %junctions; #Non-redundant list of known exon-exon junctions for all transcripts of the current gene
  my %starts;    #Non-redundant list of known start positions for all transcripts of the current gene
  my %ends;      #Non-redundant list of known end positions for all transcripts of the current gene
  foreach my $trans_id (keys %{$transcripts_ref}){
    my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

    my @starts;
    my @ends;

    #Build arrays of the exon starts and ends for this transcript
    foreach my $trans_exon_id (sort {$exons_ref->{$a}->{exon_start} <=> $exons_ref->{$b}->{exon_start}} keys %{$exons_ref}){
      my $trans_exon_start = $exons_ref->{$trans_exon_id}->{exon_start};
      my $trans_exon_end = $exons_ref->{$trans_exon_id}->{exon_end};
      push (@starts, $trans_exon_start);
      push (@ends, $trans_exon_end);

      $starts{$trans_exon_start}{label} = "na";
      $ends{$trans_exon_end}{label} = "na";
    }
    my $exon_count = scalar(@starts);

    for (my $i = 0; $i < $exon_count-1; $i++){

      my $junction = "$ends[$i]"."_"."$starts[$i+1]";

      if ($junctions{$junction}){
	$junctions{$junction}{trans_count}++;
	push(@{$junctions{$junction}{transcripts}}, $trans_id)
      }else{
	$junctions{$junction}{trans_count}= 1;
        my @temp;
	push(@temp, $trans_id);
        $junctions{$junction}{transcripts}  = \@temp;
      }
    }
  }
  $gene_transcripts_ref->{$gene_id}->{junctions} = \%junctions;

  #Now assign 'exon order' to each exon content block
  my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};
  my $exon_count = keys %{$exon_content_ref};

  my $exon_order = 0;

  foreach my $exon_id (sort {$exon_content_ref->{$a}->{start} <=> $exon_content_ref->{$b}->{start}} keys %{$exon_content_ref}){
    my $ec_start = $exon_content_ref->{$exon_id}->{start};
    my $ec_end = $exon_content_ref->{$exon_id}->{end};

    #Note that since we are dealing with GENE coordinates which are always stored 5' to 3', the strand does no matter at this stage
    $exon_order++;
    $exon_content_ref->{$exon_id}->{exon_order} = $exon_order;

    #For this exon content block, identify the exon starts and ends that are within it:
    #Then label the donors (exon ends) and acceptor (exon starts) from A to Z (5' to 3')

    my %ec_starts; #Starts that are actually within the current exon content block
    my %ec_ends;   #Ends that are actually within the current exon content block

    foreach my $start (keys %starts){
      if ($start >= $ec_start && $start <= $ec_end){
	$ec_starts{$start}{tmp} = "na";
      }
    }
    foreach my $end (keys %ends){
      if ($end >= $ec_start && $end <= $ec_end){
	$ec_ends{$end}{tmp} = "na";
      }
    }

    #Now apply the labels to each Acceptor or Donor
    my @labels = @letters;
    foreach my $ec_start (sort {$a <=> $b} keys %ec_starts){
      $starts{$ec_start}{label} = shift (@labels);
    }
    @labels = @letters;
    foreach my $ec_end (sort {$a <=> $b} keys %ec_ends){
      $ends{$ec_end}{label} = shift (@labels);
    }
  }

  $gene_exon_content_ref->{$gene_id}->{acceptors} = \%starts;
  $gene_exon_content_ref->{$gene_id}->{donors} = \%ends;
}

print BLUE, "\n2-B.) Getting EXON CONTENT of each gene\n\n", RESET;
print LOG "\n2-B.) Getting EXON CONTENT of each gene\n\n";
my $counter = 0;
foreach my $gene_id (@gene_ids){

  my $chromosome = $genes_ref->{$gene_id}->{chromosome};
  my $chr_strand = $genes_ref->{$gene_id}->{chr_strand};
  my $chr_start = $genes_ref->{$gene_id}->{chr_start};
  my $chr_end = $genes_ref->{$gene_id}->{chr_end};
  my $gene_start = $genes_ref->{$gene_id}->{gene_start};
  my $gene_end = $genes_ref->{$gene_id}->{gene_end};

  $counter++;
  if ($counter == 100){
    $counter = 0;
    $| = 1; print BLUE, ".", RESET; $| = 0;
  }

  #Calculate the size of each transcript by adding up the size of its exons
  my $size = 0;
  my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};

  $genes_ref->{$gene_id}->{exon_content_size} = 0; #Total number of bases covered by exons of this gene

  foreach my $exon_id (keys %{$exon_content_ref}){

    my $start = $exon_content_ref->{$exon_id}->{start};
    my $end = $exon_content_ref->{$exon_id}->{end};

    my $coords_ref = &convertGeneCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$start, '-end_pos'=>$end, '-ordered'=>"yes");

    $exon_content_ref->{$exon_id}->{chr_start} = $coords_ref->{$gene_id}->{chr_start};
    $exon_content_ref->{$exon_id}->{chr_end} = $coords_ref->{$gene_id}->{chr_end};
    $exon_content_ref->{$exon_id}->{size} = ($coords_ref->{$gene_id}->{chr_end} - $coords_ref->{$gene_id}->{chr_start})+1;
  }
}


#3.) Import all mRNA and EST observed exon-exon junctions
#    - NOTE that these databases are massive
#    - Organize into hashes according to chromosome and strand
#    - Use berkley DB to store these data structures temporarily and clean-up at the end of the script
my $mrna_junctions = &importUcscAlignments('-import_type'=>"junction", '-align_file'=>$ucsc_mrna_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"mrna", '-filter'=>"0", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);
my $est_junctions = &importUcscAlignments('-import_type'=>"junction",'-align_file'=>$ucsc_est_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"est", '-filter'=>"1", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);
my $xmrna_junctions = &importUcscAlignments('-import_type'=>"junction",'-align_file'=>$ucsc_xmrna_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"xmrna", '-filter'=>"2", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);
my $xest_junctions = &importUcscAlignments('-import_type'=>"junction",'-align_file'=>$ucsc_xest_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"xest", '-filter'=>"2", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);

#print Dumper $mrna_junctions;
#foreach my $chr_strand (keys %{$mrna_junctions}){
#  my $ref = $mrna_junctions->{$chr_strand};
#  foreach my $exon (sort keys %{$ref}){
#    print YELLOW, "\n$exon", RESET;
#  }
#}


#4.) Go through each exon-exon junction and get the number of exons skipped for its theoretical target
#    - Also determine whether it corresponds to a known EnsEMBL transcript
#    - Also determine whether is corresponds to a known EST or mRNA mapped to the same strand by UCSC
my $total_known_junction_count = 0;
my $unknown_chromosome_records = 0;
my $mrna_supported_count = 0;
my $est_supported_count = 0;
my $xmrna_supported_count = 0;
my $xest_supported_count = 0;

print BLUE, "\n\nDetermining the number of exon skips for each exon junction - and whether is corresponds to a known EnsEMBL transcript or not\n\n", RESET;
print LOG "\n\nDetermining the number of exon skips for each exon junction - and whether is corresponds to a known EnsEMBL transcript or not\n\n";

print BLUE, "\n\nAt this time also annotate each junction according to the exons they involve (e.g. exon1_exon2 or exon1_exon3)\n\n", RESET;
print LOG "\n\nnAt this time also annotate each junction according to the exons they involve (e.g. exon1_exon2 or exon1_exon3)\n\n";

my %failed_chrs;

$counter = 0;
foreach my $junction_id (sort {$a <=> $b} keys %junctions){

  $counter++;
  if ($counter == 10000){
    print BLUE ".";
    $counter = 0;
  }

  my $gene_id = $junctions{$junction_id}{gene_id};

  my $exon1_start = $junctions{$junction_id}{unit1_start};
  my $exon1_end = $junctions{$junction_id}{unit1_end};
  my $exon2_start = $junctions{$junction_id}{unit2_start};
  my $exon2_end = $junctions{$junction_id}{unit2_end};

  my $exon1_start_chr = $junctions{$junction_id}{unit1_start_chr};
  my $exon1_end_chr = $junctions{$junction_id}{unit1_end_chr};
  my $exon2_start_chr = $junctions{$junction_id}{unit2_start_chr};
  my $exon2_end_chr = $junctions{$junction_id}{unit2_end_chr};

  my $exon1_end_chr_sorted = $junctions{$junction_id}{unit1_end_chr};
  my $exon2_start_chr_sorted = $junctions{$junction_id}{unit2_start_chr};


  #Sanity check 
  unless ($exon2_start > $exon1_end){
    print RED, "\nExon1 end and Exon2 start positions do not make sense!\n\n", RESET;
    close (LOG);
    exit();
  }
  my $skipped_exons = 0;

  #Compare the coordinate of this junction to the reference exons and identify the number of skipped exons
  my @ref_exon_starts = @{$gene_transcripts_ref->{$gene_id}->{ref_exon_starts}};
  my @ref_exon_ends = @{$gene_transcripts_ref->{$gene_id}->{ref_exon_ends}};
  my $ref_exon_count = $gene_transcripts_ref->{$gene_id}->{ref_exon_count};

  my $found_exon1 = 0;

  #First identify the last exon which corresponds to the exon1_end position
  for (my $i = 0; $i < $ref_exon_count-1; $i++){

    #Look for the exon which corresponds to the exon2_start position
    if ($exon2_start >= $ref_exon_starts[$i] && $exon2_start <= $ref_exon_ends[$i]){
      last();
    }

    #Look for the exon which corresponds to the exon1_end position
    if ($exon1_end >= $ref_exon_starts[$i] && $exon1_end <= $ref_exon_ends[$i]){
      $found_exon1 = 1;
      next();
    }

    #If exon1 has been found in the ordered reference set but exon2 has not been found yet, this is a skipped exon
    #The boundaries of the current reference exon should also be completely between the two exons targetted
    if ($found_exon1 == 1 && $exon1_end < $ref_exon_starts[$i] && $exon2_start > $ref_exon_starts[$i] && $exon1_end < $ref_exon_ends[$i] && $exon2_start > $ref_exon_ends[$i]){
      $skipped_exons++;
      next();
    }
  }
  $junctions{$junction_id}{skipped_exons} = $skipped_exons;


  #Check if the exon-exon junction coordinates of this exon-exon junction correspond to a known EnsEMBL transcript
  my $current_junction =  "$exon1_end"."_"."$exon2_start";

  my $junctions_ref = $gene_transcripts_ref->{$gene_id}->{junctions};

  $junctions{$junction_id}{specific_trans_id} = "na";

  if ($junctions_ref->{$current_junction}){
    $junctions{$junction_id}{known_junction} = $junctions_ref->{$current_junction}->{trans_count};
    $total_known_junction_count++;

    #Also note whether it is specific to a single transcript 
    # - Always the case when there is only one transcript, but rarely the case when multiple transcripts exist
    # - Record the transcript ID for these cases, 
    # - Denote each transcript ID as 'S' or 'M' according to whether they correspond to a Single or Multiple transcript gene

    my $trans_id_modifier;
    if ($gene_transcripts_ref->{$gene_id}->{trans_count} == 1){
      $trans_id_modifier = "S";
    }else{
      $trans_id_modifier = "M";
    }	    

    if ($junctions_ref->{$current_junction}->{trans_count} == 1){
      my $trans_id = @{$junctions_ref->{$current_junction}->{transcripts}}[0];
      $junctions{$junction_id}{specific_trans_id} = "$trans_id_modifier"."_"."$trans_id";
    }
  }else{
    $junctions{$junction_id}{known_junction} = 0;
  }

  #Test the current junction against junctions from mRNAs/ESTs/xmRNAs/xESTs
  #Allow the junction to wiggle on either side by the amount specified by the user
  #Since seqs from other species rarely align perfectly, good support for an exon skipping event can be missed if you require the junction to match exactly to the base
  $junctions{$junction_id}{mrna_supported_junction} = 0;
  $junctions{$junction_id}{est_supported_junction} = 0;
  $junctions{$junction_id}{xmrna_supported_junction} = 0;
  $junctions{$junction_id}{xest_supported_junction} = 0;
  $junctions{$junction_id}{conserved_species_count} = 0;
  my $current_chr_strand = "$junctions{$junction_id}{chromosome}"."_"."$junctions{$junction_id}{strand}";
  &testJunctionSupport('-junction_id'=>$junction_id, '-junction_start'=>$exon1_end_chr_sorted, '-junction_end'=>$exon2_start_chr_sorted, '-chr_strand'=>$current_chr_strand, '-wiggle'=>$wiggle);


  #Annotate each exon-exon junction to identify which exons are involved
  #Use the 'exon_order' values assigned to each exon content block
  #Label will be according to the block of overlapping exons (if applicable)
  #In the case of overlapping exons, alternate acceptor, or donor positions will be labeled

  my $strand = $genes_ref->{$gene_id}->{chr_strand};
  my $exon1_order_value = "na";
  my $exon2_order_value = "na";

  my $exon1_found = 0;
  my $exon2_found = 0;

  my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};

  foreach my $exon_id (sort {$exon_content_ref->{$a}->{start} <=> $exon_content_ref->{$b}->{start}} keys %{$exon_content_ref}){

    my $ref_exon_start = $exon_content_ref->{$exon_id}->{start};
    my $ref_exon_end = $exon_content_ref->{$exon_id}->{end};
    my $ref_exon_order = $exon_content_ref->{$exon_id}->{exon_order};

    #Is the current junction (end of unit1 portion) contained within this exon content block?
    if ($exon1_end >= $ref_exon_start && $exon1_end <= $ref_exon_end){
      $exon1_found++;
      $exon1_order_value = $ref_exon_order;
    }

    #Is the current junction (start of unit2 portion) contained with this exon content block?
    if ($exon2_start >= $ref_exon_start && $exon2_start <= $ref_exon_end){
      $exon2_found++;
      $exon2_order_value = $ref_exon_order;
    }
  }

  if ($exon1_found > 1 || $exon2_found > 1){
    #Since we are comparing exon junction coordinates to exon-content blocks, we should always find only a single target exon
    print RED,  "\n\tAmbigous junction to exon mapping!  Junction: $junction_id\tExon1's Found = $exon1_found\tExon2's Found = $exon2_found", RESET;
    exit();

  }elsif ($exon1_found == 0 || $exon2_found == 0){
    #Every junction coordinate should correspond to at least one exon content block
    print RED,  "\n\tUnassiged junction to exon mapping!  Junction: $junction_id", RESET;
    exit();

  }else{

    #print YELLOW, "\n\tJunction: $junction_id\tStrand: $strand\tWill be annotated as: $exon1_order_value"."_"."$exon2_order_value", RESET;

  }

  #Assign the exon_names value for this junction

  my $acceptors_ref = $gene_exon_content_ref->{$gene_id}->{acceptors};
  my $donors_ref = $gene_exon_content_ref->{$gene_id}->{donors};

  my $acceptor_label = $acceptors_ref->{$exon2_start}->{label};
  my $donor_label = $donors_ref->{$exon1_end}->{label};

  unless (defined($exon1_order_value)){$exon1_order_value=0;}
  unless (defined($exon2_order_value)){$exon2_order_value=0;}
  unless (defined($donor_label)){$donor_label = "x";}
  unless (defined($acceptor_label)){$acceptor_label = "x";}
  $junctions{$junction_id}{exon_names} = "E"."$exon1_order_value"."$donor_label"."_E"."$exon2_order_value"."$acceptor_label";

}

print BLUE, "\n\n$total_known_junction_count junctions of $total_junction_count correspond to known EnsEMBL transcripts\n\n", RESET;
print LOG "\n\n$total_known_junction_count junctions of $total_junction_count correspond to known EnsEMBL transcripts\n\n";
print BLUE, "\n\n$mrna_supported_count junctions of $total_junction_count are supported by one or more mRNAs\n\n", RESET;
print LOG "\n\n$mrna_supported_count junctions of $total_junction_count are supported by one or more mRNAs\n\n", RESET;
print BLUE, "\n\n$est_supported_count junctions of $total_junction_count are supported by one or more ESTs\n\n", RESET;
print LOG "\n\n$est_supported_count junctions of $total_junction_count are supported by one or more ESTs\n\n", RESET;
print BLUE, "\n\n$xmrna_supported_count junctions of $total_junction_count are supported by one or more xmRNAs\n\n", RESET;
print LOG "\n\n$xmrna_supported_count junctions of $total_junction_count are supported by one or more xmRNAs\n\n", RESET;
print BLUE, "\n\n$xest_supported_count junctions of $total_junction_count are supported by one or more xESTs\n\n", RESET;
print LOG "\n\n$xest_supported_count junctions of $total_junction_count are supported by one or more xESTs\n\n", RESET;
print BLUE, "\n\n$unknown_chromosome_records junctions of $total_junction_count had an unrecognized chromosome name\n\n", RESET;
print LOG "\n\n$unknown_chromosome_records junctions of $total_junction_count had an unrecognized chromosome name\n\n";


#6.) Finally print out the results
print BLUE, "\nPrinting output to file: $outfile\n\n", RESET;
print LOG "\nPrinting output to file: $outfile\n\n";

open (JUNCTIONS, "$junction_file") || die "\nCould not open input junction file: $junction_file\n\n";
open (OUTFILE, ">$outfile") || die "\nCould not open output junction file: $outfile\n\n";

$header = 1;
while (<JUNCTIONS>){
  chomp($_);
  if ($header == 1){
    $header = 0;
    print OUTFILE "$_\tExons_Skipped\tSupporting_EnsEMBL_Count\tSupporting_mRNA_Count\tSupporting_EST_Count\tSupporting_xmRNA_Count\tSupporting_xEST_Count\tConserved_Species_Count\tSeq_Name\tSpecific_Trans_ID\n";
    next();
  }

  my @line = split("\t", $_);
  my $junction_id = $line[0];
 
  #Only print out records for junctions within the current region
  unless($junctions{$junction_id}){
    next();
  }

  my $gene_id = $junctions{$junction_id}{gene_id};

  unless ($junctions{$junction_id}){
    print RED, "\nProbe_id not found in hash!\n\n", RESET;
    close (LOG);
    exit();
  }
  chomp ($_);
  print OUTFILE "$_\t$junctions{$junction_id}{skipped_exons}\t$junctions{$junction_id}{known_junction}\t$junctions{$junction_id}{mrna_supported_junction}\t$junctions{$junction_id}{est_supported_junction}\t$junctions{$junction_id}{xmrna_supported_junction}\t$junctions{$junction_id}{xest_supported_junction}\t$junctions{$junction_id}{conserved_species_count}\t$junctions{$junction_id}{exon_names}\t$junctions{$junction_id}{specific_trans_id}\n";
}
print "SCRIPT COMPLETE";
print LOG "SCRIPT COMPLETE";

close (JUNCTIONS);
close (OUTFILE);
close (LOG);

untie(%gb_org_map);

exit();


##############################################################################################################################3
#Test support for a junction against mRNA/EST/xmRNA/xEST
##############################################################################################################################3
sub testJunctionSupport{
  my %args = @_;
  my $junction_id = $args{'-junction_id'};
  my $current_chr_strand = $args{'-chr_strand'};
  my $junction_start = $args{'-junction_start'};
  my $junction_end = $args{'-junction_end'};
  my $wiggle = $args{'-wiggle'};

  my $current_junction_chr = "$junction_start"."_"."$junction_end";
  #print YELLOW, "\nDEBUG: $junction_id\t$current_chr_strand\t$junction_start\t$junction_end\t$wiggle\t$current_junction_chr", RESET;

  my $mrna_supported_test = 0;
  my $est_supported_test = 0;
  my $xmrna_supported_test = 0;
  my $xest_supported_test = 0;

  #TARGET SPECIES mRNAS - Check if the exon-exon junction coordinates of this exon-exon junction correspond to a known mRNA sequence from the target species
  if ($mrna_junctions->{$current_chr_strand}){
    my $mrna_junctions_ref = $mrna_junctions->{$current_chr_strand};

    if ($mrna_junctions_ref->{$current_junction_chr}){
      $junctions{$junction_id}{mrna_supported_junction} = $mrna_junctions_ref->{$current_junction_chr}->{c};
      $mrna_supported_test = 1;
    }
  }else{
    $unknown_chromosome_records++;
    $junctions{$junction_id}{mrna_supported_junction} = "na";

    if ($failed_chrs{$junctions{$junction_id}{chromosome}}){
      $failed_chrs{$junctions{$junction_id}{chromosome}}{count}++;
    }else{
      $failed_chrs{$junctions{$junction_id}{chromosome}}{count} = 1;
    }
  }
  #TARGET SPECIES ESTS - Check if the exon-exon junction coordinates of this exon-exon junction correspond to a known EST sequence from the target species
  if ($est_junctions->{$current_chr_strand}){
    my $est_junctions_ref = $est_junctions->{$current_chr_strand};

    if ($est_junctions_ref->{$current_junction_chr}){
      $junctions{$junction_id}{est_supported_junction} = $est_junctions_ref->{$current_junction_chr}->{c};
      $est_supported_test = 1;
    }
  }else{
    $junctions{$junction_id}{est_supported_junction} = "na";
  }

  #Use the amount of 'wiggle' specified by the user to generate a list of junctions
  my @s;
  my @e;
  push(@s, $junction_start);
  push(@e, $junction_end);
  for (my $i = 1; $i <= $wiggle; $i++){
    push(@s, $junction_start+$i);
    push(@s, $junction_start-$i);
    push(@e, $junction_end+$i);
    push(@e, $junction_end-$i);
  }
  my @current_junction_list;
  foreach my $start (@s){
    foreach my $end (@e){
      my $current_junction = "$start"."_"."$end";
      push (@current_junction_list, $current_junction);
    }
  }

  #Keep a list of all species IDs encountered for the current junction_id
  my %species_list;

  #Only apply the wiggle to xEST and xmRNAs!
  foreach my $current_junction_chr (@current_junction_list){

    #OTHER SPECIES MRNAS - Check if the exon-exon junction coordinates of this exon-exon junction correspond to a known mRNA sequence from any other species
    if ($xmrna_junctions->{$current_chr_strand}){
      my $xmrna_junctions_ref = $xmrna_junctions->{$current_chr_strand};

      if ($xmrna_junctions_ref->{$current_junction_chr}){
        $junctions{$junction_id}{xmrna_supported_junction} += $xmrna_junctions_ref->{$current_junction_chr}->{c};
        $xmrna_supported_test = 1;
        my $s_ref = $xmrna_junctions_ref->{$current_junction_chr}->{s};
        foreach my $s (keys %{$s_ref}){
          $species_list{$s}='';
        }
      }
    }else{
      $junctions{$junction_id}{xmrna_supported_junction} = "na";
    }
    #OTHER SPECIES ESTS - Check if the exon-exon junction coordinates of this exon-exon junction correspond to a known EST sequence from any other species
    if ($xest_junctions->{$current_chr_strand}){
      my $xest_junctions_ref = $xest_junctions->{$current_chr_strand};

      if ($xest_junctions_ref->{$current_junction_chr}){
        $junctions{$junction_id}{xest_supported_junction} += $xest_junctions_ref->{$current_junction_chr}->{c};
        $xest_supported_test = 1;
        my $s_ref = $xest_junctions_ref->{$current_junction_chr}->{s};
        foreach my $s (keys %{$s_ref}){
          $species_list{$s}='';
        }
      }
    }else{
      $junctions{$junction_id}{xest_supported_junction} = "na";
    }
  }

  if ($mrna_supported_test == 1){$mrna_supported_count++;}
  if ($est_supported_test == 1){$est_supported_count++;}
  if ($xmrna_supported_test == 1){$xmrna_supported_count++;}
  if ($xest_supported_test == 1){$xest_supported_count++;}

  #Count up all the species encountered in junction matches from xmrnas and xests
  my $conserved_species_count = keys %species_list;
  $junctions{$junction_id}{conserved_species_count} = $conserved_species_count;

  return();
}




