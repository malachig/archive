#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to create a database of intron sequences for the purposes of read mapping
#These sequences should be repeat masked
#An 'intron' in this context means the sequence between exon content blocks (which are created by merging overlapping exons of all transcripts of a single gene locus)
#Gene models are EnsEMBL transcripts/genes taken from an ALEXA database

#A chr_filter option is used to process one chromosome at a time for performance reasons

#1-A.) Get gene models from ALEXA for all genes - get completed masked sequence for each gene locus
#1-B.) Get exon content blocks for each gene (use existing method in ALEXA_DB.pm)
#1-C.) Using the exon content blocks, extract intron content coordinates (use existing method in ALEXA_DB.pm)
#    - One of the caveats of this approach is that it will miss EnsEMBL known intron retentions
#    - However, reads which map to these 'introns' will captured when mapping to the known transcript with the intron retention
#    - These intron known intron retension can also be identified by looking at known exon boundary extension sequences 
#      (i.e. the boundaries that correspond to the exons of the transcript which does not retain the intron)
#    - Finally, each known intron retension should have already been catalogued previously as an 'Exon Region'.
#1-D.) Get chromosome coordinates for all exon and intron regions

#2.) Go through UCSC mRNA AND EST databases and determine the coverage of each intronic region.
#    - Using this coverage map, identify active and silent regions within each intron
#    - For ESTs, apply basic quality filters, only used spliced ESTs, etc.
#    - Use this EST and mRNA data to define the regions of each intron that have evidence for expression versus those that do not...
#    - i.e. identify the following two sets of regions: IntronRegions_Active vs IntronRegions_Silent

#3.) Print output files:

#3-A.) Print out intron sequence records to a DB file
#    - Intron_ID, Gene_ID, EnsEMBL_Gene_ID, Intron_Name, Chromosome, Strand, Unit1_start, Unit1_end, Unit1_start_chr, Unit1_end_chr, Length, UnMasked_Base_Count, Active_Base_Count, Silent_Base_Count, MultipleGenes, Gene_ID_list 

#3-B.) Print out a fasta file of the intron sequences
#    i.e. ">$intron_id\n$intron_sequence\n"

#3-C.) Print out DB files for 'IntronRegion_Active' and 'IntronRegion_Silent'
#      - Active_Region_ID, Intron_ID, Gene_ID, EnsEMBL_Gene_ID, Active_Region_Name, Chromosome, Strand, Unit1_start, Unit1_end, Unit1_start_chr, Unit1_end_chr, Length, UnMasked_Base_Count
#      - Silent_Region_ID, Intron_ID, Gene_ID, EnsEMBL_Gene_ID, Silent_Region_Name, Chromosome, Strand, Unit1_start, Unit1_end, Unit1_start_chr, Unit1_end_chr, Length, UnMasked_Base_Count


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
use utilities::ucsc qw(:all);

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $min_size = '';
my $chr_filter = '';
my $ucsc_align_dir = '';
my $genbank_mapfile = '';
my $wiggle = '';
my $outdir = '';
my $logfile = '';
my $verbose = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'min_size=i'=>\$min_size,
            'chr_filter=s'=>\$chr_filter, 'ucsc_align_dir=s'=>\$ucsc_align_dir, 'genbank_mapfile=s'=>\$genbank_mapfile, 'wiggle=i'=>\$wiggle,
            'outdir=s'=>\$outdir, 'logfile=s'=>\$logfile, 'verbose=s'=>\$verbose);


#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the minimum intron size (introns with less unmasked bases than this will be skipped) using:  --min_size", RESET;
print GREEN, "\n\tSpecify a single chromosome region to be processed using: --chr_filter", RESET;
print GREEN, "\n\tSpecify the path to be used for output files using: --outdir", RESET;
print GREEN, "\n\tSpecify a directory containing UCSC mRNA/EST/xmRNA and xEST files using: --ucsc_align_dir", RESET;
print GREEN, "\n\tSpecify a berkeley DB containing GenBank-To-Species mappings using:  --genbank_mapfile", RESET;
print GREEN, "\n\tJunctions for xmRNA and xEST will be allowed to 'wiggle' by the amount specied by:  --wiggle", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of the design process using: --logfile", RESET;
print GREEN, "\n\tIf you want verbose output, use: --verbose=yes", RESET;

print GREEN, "\n\nExample: createIntronDatabase.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --min_size=42  --chr_filter='3:16:121020102-126260001'  --ucsc_align_dir=/projects/malachig/sequence_databases/hs_53_36o/mrna_est/partitions/  --genbank_mapfile=/projects/malachig/sequence_databases/hs_53_36o/mrna_est/GenBankToOrganism.btree  --wiggle=5  --outdir=/projects/malachig/sequence_databases/hs_53_36o/ensembl_introns_hs_53_36o/temp/  --logfile=/projects/malachig/sequence_databases/hs_53_36o/logs/createIntronDatabase/createIntronDatabase_3_16_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $min_size && $chr_filter && $ucsc_align_dir && $genbank_mapfile && ($wiggle =~ /\d+/) && $outdir && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
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

if ($chr_filter eq "MT"){$chr_filter = "M";}
my $tmp_chr = "chr$chr_filter";

$ucsc_align_dir = &checkDir('-dir'=>$ucsc_align_dir, '-clear'=>"no");
my $ucsc_mrna_table = "$ucsc_align_dir"."$tmp_chr"."_mrna.txt.gz";
my $ucsc_est_table = "$ucsc_align_dir"."$tmp_chr"."_est.txt.gz";
my $ucsc_xmrna_table = "$ucsc_align_dir"."$tmp_chr"."_xmrna.txt.gz";
my $ucsc_xest_table = "$ucsc_align_dir"."$tmp_chr"."_xest.txt.gz";

#Load the GenBank-To-Species mapfile
my %gb_org_map;
tie(%gb_org_map, 'BerkeleyDB::Btree', -Cachesize =>256000000, -Filename=> $genbank_mapfile , -Flags => DB_RDONLY) or die "can't open file $genbank_mapfile: $! $BerkeleyDB::Error\n";

#Check output dir
$outdir = &checkDir('-dir'=>$outdir, '-clear'=>"no");

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#Global counters
my $grand_base_count = 0;
my $grand_intron_seq_overlaps = 0;
my $grand_cumulative_coverage_bases = 0;
my $grand_unique_covered_bases = 0;
my $grand_active_region_count = 0;
my $grand_silent_region_count = 0;
my $grand_small_introns = 0;

#1.) Get gene models from ALEXA for all genes - get completed masked sequence for each gene locus

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
my $genes_ref;
my $genes_ref_masked;
my $gene_intron_content_ref;
my $current_intron_id = 0;
my $current_active_region_id = 0;
my $current_silent_region_id = 0;

$| = 1; print BLUE, "\n\n1.) Getting basic gene info for chr: $chr_filter:$start_filter-$end_filter", RESET; $| = 0;
print LOG "\n\n1.) Getting basic gene info for chr: $chr_filter:$start_filter-$end_filter";
my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All', '-chromosome'=>$chr_filter, '-range'=>"$start_filter-$end_filter")};
my %introns;
&getBasicGeneInfo ('-gene_ids'=> \@gene_ids);

#Close database connection
$alexa_dbh->disconnect();

#For both EST and mRNA seq files, write a berkley DB containing only those records corresponding to the current chromosome
my $mrna_ref = &parseSeqFile('-infile'=>$ucsc_mrna_table);
my $est_ref = &parseSeqFile('-infile'=>$ucsc_est_table);
my $mrna_record_count = keys %{$mrna_ref};
my $est_record_count = keys%{$est_ref};

#Use a different function to get the coordinates of mrna, est, xmrna and xest exons
my $mrna_exons = &importUcscAlignments('-import_type'=>"exon", '-align_file'=>$ucsc_mrna_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"mrna", '-filter'=>"0", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);
my $est_exons = &importUcscAlignments('-import_type'=>"exon", '-align_file'=>$ucsc_est_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"est", '-filter'=>"1", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);
my $xmrna_exons = &importUcscAlignments('-import_type'=>"exon", '-align_file'=>$ucsc_xmrna_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"xmrna", '-filter'=>"2", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);
my $xest_exons = &importUcscAlignments('-import_type'=>"exon", '-align_file'=>$ucsc_xest_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"xest", '-filter'=>"2", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);

#Now process intron-by-intron
my %intron_coverage;

my $count_down = keys %introns;

$| = 1; print BLUE, "\n\n2.) Initializing coverage hash, adding coverage from mRNAs/ESTs, determining 'Active' and 'Silent regions.  Intron-By-Intron\n", RESET; $| = 0;
print LOG "\n\n2.) Initializing coverage hash, adding coverage from mRNAs/ESTs, determining 'Active' and 'Silent regions.  Intron-By-Intron\n";

foreach my $intron_id (sort {$introns{$a}->{start_chr} <=> $introns{$b}->{start_chr}} keys %introns){
  $count_down--;

  #Skip failing introns and intergenic regions!!
  unless (($introns{$intron_id}{pass} == 1) && $introns{$intron_id}{gene_id} =~ /\d+/){
    next();
  }

  my $strand = $introns{$intron_id}{strand};
  $| = 1; print BLUE, "\n\n($count_down) Processing Non-Exonic Region $intron_id\tchr$chr_filter:$introns{$intron_id}{start_chr}-$introns{$intron_id}{end_chr} (strand: $strand) (Gene IDs: @{$introns{$intron_id}{gene_id_list}})", RESET; $| = 0;
  print LOG "\n\n($count_down) Processing Non-Exonic Region $intron_id\tchr$chr_filter:$introns{$intron_id}{start_chr}-$introns{$intron_id}{end_chr} (strand: $strand) (Gene IDs: @{$introns{$intron_id}{gene_id_list}})";
 

  #2.) Go through UCSC mRNA AND EST databases and determine the coverage of each intronic region.
  %intron_coverage = ();

  #2-A.) Initialize a coverage hash for all intron base positions (keyed as geneID_position)
  &initializeCoverage('-intron_id'=>$intron_id);

  #2-B.) Add coverage from ESTs/mRNAs to this hash
  if ($mrna_record_count > 0){
    &addCoverage('-alignment_object'=>$mrna_ref, '-intron_id'=>$intron_id, '-type'=>'mRNA');
  }
  if ($est_record_count > 0){
    &addCoverage('-alignment_object'=>$est_ref, '-intron_id'=>$intron_id, '-type'=>'EST');
  }

  #2-C.) Extract the coordinates corresponding to 'Active' and 'Silent' regions within each intron
  &identifyActiveSilentRegions('-intron_id'=>$intron_id);

  #2-D.) Determine the conservation support
  &determineSequenceSupport();

}

#3.) Print output files
$| = 1; print BLUE, "\n\n\n3.) Printing output files to: $outdir\n", RESET; $| = 0;
print LOG "\n\n\n3.) Printing output files to: $outdir\n";
&printOutputFiles();

untie(%gb_org_map);

print BLUE, "\n\nFINAL SUMMARY", RESET;
print BLUE, "\nIntron base count = $grand_base_count", RESET;
print BLUE, "\nIntron-To-Sequence (EST/mRNA) Overlaps = $grand_intron_seq_overlaps", RESET;
print BLUE, "\nCumulative EST/mRNA coverage bases = $grand_cumulative_coverage_bases", RESET;
print BLUE, "\nUnique EST/mRNA covered bases = $grand_unique_covered_bases", RESET;
print BLUE, "\nActive intron region count (EST/mRNA evidence) = $grand_active_region_count", RESET;
print BLUE, "\nSilent region count (EST/mRNA evidence) = $grand_silent_region_count", RESET;
print BLUE, "\nNumber of introns skipped because they were too small = $grand_small_introns\n\n", RESET;
print BLUE, "\nSCRIPT COMPLETE\n\n", RESET;

print LOG "\n\nFINAL SUMMARY";
print LOG "\nIntron base count = $grand_base_count";
print LOG "\nIntron-To-Sequence (EST/mRNA) Overlaps = $grand_intron_seq_overlaps";
print LOG "\nCumulative EST/mRNA coverage bases = $grand_cumulative_coverage_bases";
print LOG "\nUnique EST/mRNA covered bases = $grand_unique_covered_bases";
print LOG "\nActive intron region count (EST/mRNA evidence) = $grand_active_region_count";
print LOG "\nSilent region count (EST/mRNA evidence) = $grand_silent_region_count";
print LOG "\nNumber of introns skipped because they were too small = $grand_small_introns\n\n";
print LOG "\nSCRIPT COMPLETE\n\n";

close(LOG);
exit();


################################################################################################
#Get gene models, masked gene sequence, intron content, etc.                                   #
################################################################################################
sub getBasicGeneInfo{
  my %args = @_;
  my @gene_ids = @{$args{'-gene_ids'}};

  #Get basic info about each gene
  #filter out genes that are not within the specified region
  $| = 1; print BLUE, "\nGetting gene models", RESET; $| = 0; 
  print LOG "\nGetting gene models";
  $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

  #Get masked sequence for each gene
  $| = 1; print BLUE, "\nGetting masked gene sequences", RESET; $| = 0; 
  print LOG "\nGetting masked gene sequences";
  $genes_ref_masked = &getMaskedGene ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);

  #Get exon content for each gene
  $| = 1; print BLUE, "\nGetting exon content", RESET; $| = 0; 
  print LOG "\nGetting exon content";
  my $gene_exon_content_ref = &getExonContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);

  #Get intron content for each gene (only used to determine names of introns. - intron regions are defined at chromosome level using exon content) 
  $| = 1; print BLUE, "\nGetting intron content", RESET; $| = 0; 
  print LOG "\nGetting intron content";
  my $gene_intron_content_ref = &getIntronContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);

  #Get chromosome coordinates for exon content and use these to define a chromosome wide exon-content hash
  $| = 1; print BLUE, "\nGetting exon content chromosome coordinates and create chromosome-wide coverage hash", RESET; $| = 0; 
  print LOG "\nGetting exon content chromosome coordinates and create chromosome-wide coverage hash";

  my %exonic_coverage_chr;
  my $file_name = "$outdir"."$chr_filter"."_"."$region_number".".exonic_coverage.btree";
  my $rm_cmd = "rm -f $file_name";
  system($rm_cmd);
  print YELLOW, "\nCreating binary tree file: $file_name", RESET;
  tie(%exonic_coverage_chr, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $file_name , -Flags => DB_CREATE) or die "can't open file $file_name: $! $BerkeleyDB::Error\n";
  my $h_ref = \%exonic_coverage_chr;

  foreach my $gene_id (keys %{$genes_ref}){
    my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};

    foreach my $exon_id (keys %{$exon_content_ref}){

      my $start = $exon_content_ref->{$exon_id}->{start};
      my $end = $exon_content_ref->{$exon_id}->{end};
      my $coords_ref = &convertGeneCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$start, '-end_pos'=>$end, '-ordered'=>"yes");

      my $chr_start = $coords_ref->{$gene_id}->{chr_start};
      my $chr_end = $coords_ref->{$gene_id}->{chr_end};

      $exon_content_ref->{$exon_id}->{chr_start} = $chr_start;
      $exon_content_ref->{$exon_id}->{chr_end} =  $chr_end;

      #Use these to population chromosome wide hash of exonic positions
      for (my $i = $chr_start; $i <= $chr_end; $i++){
	if ($h_ref->{$i}){
          $h_ref->{$i}++;
	}else{
	  $h_ref->{$i} = 1;
	}
      }
    }
  }

  #Now identify the non-redundant intron regions for this chromosome
  #To do this first identify the exonic content coordinates, and then use these to extract the intron coordinates
  my %intron_regions;
  my $previous_pos;
  my $first = 1;
  my @exon_starts;
  my @exon_ends;

  #Get the start stop positions
  foreach my $pos (sort {$a <=> $b} keys %{$h_ref}){
    if ($first == 1){
      push(@exon_starts, $pos);
      $first = 0;
      $previous_pos = $pos;
      next();
    }
    if ($pos > ($previous_pos + 1)){
      push(@exon_ends, $previous_pos);
      push(@exon_starts, $pos);
    }
    $previous_pos = $pos;
  }
  #If the final end position did not get added, add it now
  my $starts = @exon_starts;
  my $ends = @exon_ends;
  if ($starts == $ends){
    #do nothing
  }elsif($starts == ($ends + 1)){
    push (@exon_ends, $previous_pos);
  }else{
    print "\nIncorrect number of start/end coordinates!!\n\n";
    exit();
  }

  #The exonic base coverage hash for the entire chromosome should no longer be needed:
  untie(%exonic_coverage_chr);
  system($rm_cmd);


  #Now build a final list of intronic regions - perform initial filters on these before adding them in the list
  #NOTE: This list will initially contain intergenic regions.  These will be flagged later
  my $intron_count = 0;
  for (my $i = 0; $i < $starts-1; $i++){
    my $intron_start = ($exon_ends[$i])+1;
    my $intron_end = ($exon_starts[$i+1])-1;

    #Determine the size of each intron
    my $size = ($intron_end-$intron_start)+1;

    $intron_count++;
    $introns{$intron_count}{start_chr} = $intron_start;
    $introns{$intron_count}{end_chr} = $intron_end;
    $introns{$intron_count}{size} = $size;
    $introns{$intron_count}{pass} = 1;   
    $introns{$intron_count}{overlaping_gene_count} = 0;         #Number of genes that have some overlap with the region
    $introns{$intron_count}{contained_within_gene_count} = 0;   #Number of genes that completely flank the intronic region
    $introns{$intron_count}{gene_id} = "NA";                    #Will be NA for intergenic regions, for cases where the intron is contained within multiple genes, chose one arbitrarily (will be used to retrieve sequence)
    my @tmp;
    $introns{$intron_count}{gene_id_list} = \@tmp;

    #Skip introns that are too small, or have too few unmasked bases to allow mapping of reads???
    if ($size < $min_size){
      $grand_small_introns++;
      $introns{$intron_count}{pass} = 0;
    }
  }

  my $chr_intron_count = keys %introns;

  $| = 1; print BLUE, "\n\nFound a total of $chr_intron_count introns for this chromosome region: $chr_filter:$region_number:$start_filter-$end_filter", RESET; $| = 0; 
  print LOG "\n\nFound a total of $chr_intron_count introns for this chromosome region: $chr_filter:$region_number:$start_filter-$end_filter";

  #Now go through each intron region and determine:
  #1.) How many genes it overlaps (should be 0 for intergenic regions - mark these)
  #2.) How many genes completely flank it - assign one of these to allow retrieval of sequence
  $| = 1; print BLUE, "\n\nDetermining overlap between introns and genes", RESET; $| = 0; 
  print LOG "\n\nDetermining overlap between introns and genes";
 
  foreach my $i_id (keys %introns){
    my $i_start = $introns{$i_id}{start_chr};
    my $i_end = $introns{$i_id}{end_chr};

    my $overlap_count = 0;
    my $contained_count = 0;

    foreach my $gene_id (keys %{$genes_ref}){
      my $g_start = $genes_ref->{$gene_id}->{chr_start};
      my $g_end = $genes_ref->{$gene_id}->{chr_end};
      my $temp = $g_start;
      if ($g_start > $g_end){
        $temp = $g_start;
        $g_start = $g_end;
        $g_end = $temp;
      }

      #Check for overlap of the intron with each gene (or cases where the intron completely flanks a gene)
      if (($i_start >= $g_start && $i_start <= $g_end) || ($i_end >= $g_start && $i_end <= $g_end) || ($i_start < $g_start && $i_start > $g_end)){
        $overlap_count++;
      }

      #Check for cases where the intron is completely contained within a gene
      if ($i_start >= $g_start && $i_end <= $g_end){
        $contained_count++;
        $introns{$i_id}{gene_id} = $gene_id;   #Single gene ID to be used for retrieving sequence
        push(@{$introns{$i_id}{gene_id_list}}, $gene_id); #If multiple gene IDs are valid (introns of overlapping genes), keep this list
      }
    }
    $introns{$i_id}{overlaping_gene_count} = $overlap_count;
    $introns{$i_id}{contained_within_gene_count} = $contained_count;
  }

  #1.) Now go through each intron defined, and if it was found to be within a known gene, convert the chromosome coordinates to gene coordinates
  #2.) Now go through each intron and use the gene_id to get the masked/unmasked base count
  #3.) Now, if the intron is completely flanked by only one gene, determine it's name (e.g. I1, I2, etc.) using the gene coordinates determined
  $| = 1; print BLUE, "\n\nGetting gene coordinates, # unmasked bases, and intron names", RESET; $| = 0; 
  print LOG "\n\nGetting gene coordinates, # unmasked bases, and intron names";

  foreach my $i_id (keys %introns){
    $introns{$i_id}{start} = "NA";
    $introns{$i_id}{end} = "NA";
    $introns{$i_id}{strand} = "NA";
    $introns{$i_id}{multi_gene} = "NA";
    $introns{$i_id}{unmasked_base_count} = 0;
    $introns{$i_id}{active_base_count} = 0;
    $introns{$i_id}{silent_base_count} = 0;
    $introns{$i_id}{name} = "Ix";

    my $gene_id = $introns{$i_id}{gene_id};
    if ($gene_id =~ /\d+/){
      
      my $chr_start = $introns{$i_id}{start_chr};
      my $chr_end = $introns{$i_id}{end_chr};
      my $coords_ref = &convertGenomicCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$chr_start, '-end_pos'=>$chr_end, '-ordered'=>"yes");
       
      my $start = $coords_ref->{$gene_id}->{gene_start};
      my $end = $coords_ref->{$gene_id}->{gene_end};
      $introns{$i_id}{start} = $start;
      $introns{$i_id}{end} = $end;
      $introns{$i_id}{strand} = $genes_ref->{$gene_id}->{chr_strand};

      my $size = $introns{$i_id}{size};
      my $masked_gene_seq = $genes_ref_masked->{$gene_id}->{sequence};
      my $masked_region_seq = substr ($masked_gene_seq, $start, $size);
      my $masked_n_count = ($masked_region_seq =~ tr/N/N/);
      my $unmasked_count = $size - $masked_n_count;

      $introns{$i_id}{unmasked_base_count} = $unmasked_count;

      #Fail introns with too few unmasked bases
      if ($unmasked_count < $min_size){
        $grand_small_introns++;
        $introns{$i_id}{pass} = 0;
      }

      if ($introns{$i_id}{contained_within_gene_count} eq "1"){
        my $intron_content_ref = $gene_intron_content_ref->{$gene_id}->{intron_content};
        my $pos = 0;
        foreach my $ic (sort {$intron_content_ref->{$a}->{start} <=> $intron_content_ref->{$b}->{start}} keys %{$intron_content_ref}){
          $pos++;
          my $ic_start = $intron_content_ref->{$ic}->{start};
          my $ic_end = $intron_content_ref->{$ic}->{end};
          if ($start >= $ic_start && $end <= $ic_end){
            $introns{$i_id}{name} = "I$pos";
            last();
          }
        }
        $introns{$i_id}{multi_gene} = "No";
      }else{
        #If this intron is contained within multiple genes, the gene-level start/end and strand are not really valid
        #However the strand and start/end of the intron relative to a gene are still NEEDED to allow coordinate conversion later!!!
        #Both the reference gene used for this purpose as well as the complete list of overlapping genes will be stored
        $introns{$i_id}{multi_gene} = "Yes";
      }
    }
  }

  #Clear un-needed hashes to save memory
  %{$gene_exon_content_ref} = ();
  %{$gene_intron_content_ref} = ();

  return();
}


################################################################################################
#Parse an EST/mRNA file for only the current chromosome                                        #
################################################################################################
sub parseSeqFile{
  my %args = @_;
  my $infile = $args{'-infile'};

  my %seqs;
  my $record_count = 0;

  $| = 1; print BLUE, "\nGrabbing UCSC alignments that match the current chr from: $infile", RESET; $| = 0; 
  print LOG "\nGrabbing UCSC alignments that match the current chr from: $infile"; 

  unless (-e $infile){
    print YELLOW, "\n\nmRNA/EST file not found: $infile - returning 0 alignments for comparison", RESET;
    return(\%seqs);
  }

  open (ALIGN, "zcat $infile |") || die "\nCould not open mrna/est alignment file: $infile\n\n";

  while (<ALIGN>){
    chomp($_);
    my @line = split("\t", $_);
    my $matches = $line[1];
    my $mismatches = $line[2];
    my $n_count = $line[4];
    my $query_insertions = $line[5];

    my $strand = $line[9];
    my $chr = $line[14];
    my @exon_sizes = split(",", $line[19]);
    my @starts = split(",", $line[21]);
    my $exon_count = scalar(@starts);

    my $chr_test = "chr$chr_filter";

    my @adjusted_starts;
    my @ends;
    my @tmp_sizes = @exon_sizes;
    foreach my $start (@starts){
      my $exon_size = shift (@tmp_sizes);
      my $end = $start + $exon_size;

      push (@ends, $end);
      push (@adjusted_starts, $start+1);
    }

    #Unless the chromosome of the current EST/mRNA is the same as the current chromosome being processed, skip this record
    #Make sure the alignment has some overlap with the current target region
    my @tmp = (@adjusted_starts, @ends);
    my @tmp_sort = sort {$a <=> $b} @tmp;
    my $lower = $tmp_sort[0];
    my $upper = $tmp_sort[scalar(@tmp_sort)-1];
    my $chr_tmp = "chr"."$chr_filter";
    unless (($chr eq $chr_test) && ($lower >= $start_filter && $lower <= $end_filter) && ($upper >= $start_filter && $upper <= $end_filter)){
      next();
    }


    #QUALITY FILTER OF EST/mRNA SEQUENCES
    #Apply quality filters to EST/mRNA sequences...

    #if ($exon_count == 1){next();}  #Skip single exon alignments

    unless($matches > 300 && $mismatches < 5 && $n_count < 2 && $query_insertions < 5){  #Skip short alignments or those containing too many mismatches or insertions
      next();
    }

    #Store records
    $record_count++;

    $seqs{$record_count}{strand} = $strand;
    $seqs{$record_count}{starts} = \@starts;
    $seqs{$record_count}{sizes} = \@exon_sizes;
  }

  close(ALIGN);

  $| = 1; print BLUE, "\nFound a total of $record_count alignments for this chromosome, stored as a hash\n\n", RESET; $| = 0; 
  print LOG "\nFound a total of $record_count alignments for this chromosome, stored as a hash\n\n"; 

  return(\%seqs);
}


################################################################################################
#Initialize coverage hashes for the intron content of each intron                              #
################################################################################################
sub initializeCoverage{
  my %args = @_;
  my $intron_id = $args{'-intron_id'};

  my $base_count = 0;

  #Construct this hash as keyed on "IntronID_ChromosomePosition".  This make it easy to use a Berkeley DB to store this data structure if memory usage become a problem
  #For this coverage hash, inset the coordinates of each intron at each end to ensure that we are not considering the last or first base of real exons as introns

  my $chr_start = ($introns{$intron_id}{start_chr})+1;
  my $chr_end = ($introns{$intron_id}{end_chr})-1;

  for (my $i = $chr_start; $i <= $chr_end; $i++){
    my $intron_pos_id = "$intron_id"."_"."$i";
    $intron_coverage{$intron_pos_id} = 0;
    $base_count++;
    $grand_base_count++;
  }

  $| = 1; print BLUE, "\n\tBaseCount: $base_count\t", RESET; $| = 0; 
  print LOG "\n\tBaseCount: $base_count\t", RESET; 
  
  return();
}


################################################################################################
#Add coverage from mRNA/EST alignments to intron content block defined for EnsEMBL genes       #
################################################################################################
sub addCoverage{
  my %args = @_;
  my $seqs_ref = $args{'-alignment_object'};
  my $intron_id = $args{'-intron_id'};
  my $type = $args{'-type'};

  my $intron_seq_overlaps = 0;
  my $cumulative_coverage_bases = 0;
  my $unique_covered_bases = 0;

  foreach my $seq (keys %{$seqs_ref}){
    my $strand = $seqs_ref->{$seq}->{strand};
    my @exon_sizes = @{$seqs_ref->{$seq}->{sizes}};
    my @starts = @{$seqs_ref->{$seq}->{starts}};
    my @adjusted_starts;

    #For mRNAs/ESTs overlap on either strand will be considered as expression activity for those regions!  Since we dont have strand info for our reads, this is safer.

    my $chr_test = "chr$chr_filter";

    my @ends;
    foreach my $start (@starts){
      my $exon_size = shift (@exon_sizes);
      my $end = $start + $exon_size;
      push (@ends, $end);
      push (@adjusted_starts, $start+1);
    }
    my %seq_exons_tmp;
    my @tmp_ends = @ends;
    my $seq_exon_count = 0;

    foreach my $start (@adjusted_starts){
      $seq_exon_count++;
      my $end = shift(@tmp_ends);
      $seq_exons_tmp{$seq_exon_count}{start} = $start;
      $seq_exons_tmp{$seq_exon_count}{end} = $end;
    }

    #Get the outer coordinates of the EST/mRNA
    my $seq_start = $adjusted_starts[0];
    my $seq_end = $ends[(scalar(@ends)-1)];

    #Look for overlap.
    my $intron_start_chr = $introns{$intron_id}{start_chr};
    my $intron_end_chr = $introns{$intron_id}{end_chr};

    #Determine if there is overlap between this EST/mRNA and the current intron
    my $overlap = 0;
    if ($seq_start > $intron_start_chr && $seq_start < $intron_end_chr){  #Start of seq is within gene
      $overlap = 1;
    }
    if ($seq_end > $intron_start_chr && $seq_end < $intron_end_chr){      #End of seq is within gene
      $overlap = 1;
    }
    if ($seq_start <= $intron_start_chr && $seq_end >= $intron_end_chr){  #Seq completely flanks gene
      $overlap = 1;
    }

    #If overlap is present, add whatever coverage corresponds to previously identified intron coordinates
    if ($overlap == 1){
      $grand_intron_seq_overlaps++;
      $intron_seq_overlaps++;

      #Go through each 'exon' from the EST/mRNA sequence and add coverage where the coordinates correspond to intron content blocks
      foreach my $seq_exon (keys %seq_exons_tmp){

        my $seq_exon_start = $seq_exons_tmp{$seq_exon}{start};
        my $seq_exon_end = $seq_exons_tmp{$seq_exon}{end};

        for (my $i = $seq_exon_start; $i <= $seq_exon_end; $i++){
          my $intron_pos_id = "$intron_id"."_"."$i";

          #If this base position corresponds to an intron base position, increment it
          if (defined($intron_coverage{$intron_pos_id})){
            $grand_cumulative_coverage_bases++;
            $cumulative_coverage_bases++;
            $intron_coverage{$intron_pos_id}++;
          }
        }
      }
    }
  }

  #How many unique intronic base positions were actually covered by 1 or mRNA/ESTs
  foreach my $intron_pos_id (keys %intron_coverage){
    if ($intron_coverage{$intron_pos_id} > 0){
      $grand_unique_covered_bases++;      
      $unique_covered_bases++;
    }
  }
  $| = 1; print BLUE, "\n\t$type Coverage.\tIntron-Seq-Overlaps: $intron_seq_overlaps\tCumulativeCoverageBases: $cumulative_coverage_bases\tUniqueCoveredBases: $unique_covered_bases\t", RESET; $| = 0;
  print LOG "\n\t$type Coverage.\tIntron-Seq-Overlaps: $intron_seq_overlaps\tCumulativeCoverageBases: $cumulative_coverage_bases\tUniqueCoveredBases: $unique_covered_bases\t";

  return();
}

################################################################################################
#Extract the coordinates corresponding to 'Active' and 'Silent' regions within each intron  #
################################################################################################
sub identifyActiveSilentRegions{
  my %args = @_;
  my $intron_id = $args{'-intron_id'};

  my @active_coords;
  my @silent_coords;

  #Now get all the intronic positions specifically for this gene
  foreach my $intron_pos_id (keys %intron_coverage){
    my $tmp_intron_id;
    my $pos;

    if ($intron_pos_id =~ /(\d+)\_(\d+)/){
      $tmp_intron_id = $1;
      $pos = $2;
    }else{
      print RED, "\n\nSomething wrong in intron_position hash...\n", RESET;
      exit();
    }

    if ($intron_coverage{$intron_pos_id} > 0){
      push(@active_coords, $pos);
    }else{
      push(@silent_coords, $pos);
    }
  }

  #Now use the arrays of 'Active' and 'Silent' coords to identify the coordinates of 'Active' and 'Silent' regions
  my %active_regions;
  my %silent_regions;

  my @active_coords_sort = sort{$a<=>$b}(@active_coords);
  my @silent_coords_sort = sort{$a<=>$b}(@silent_coords);
  @silent_coords = ();
  @active_coords = ();

  #A.) Active regions first
  my $previous_pos;
  my $first = 1;
  my @active_starts;
  my @active_ends;

  #Get the start stop positions
  foreach my $pos (@active_coords_sort){
    if ($first == 1){
      push(@active_starts, $pos);
      $first = 0;
      $previous_pos = $pos;
      next();
    }
    if ($pos > ($previous_pos + 1)){
      push(@active_ends, $previous_pos);
      push(@active_starts, $pos);
    }
    $previous_pos = $pos;
  }
  #If the final end position did not get added, add it now
  my $starts = @active_starts;
  my $ends = @active_ends;
  if ($starts == $ends){
    #do nothing
  }elsif($starts == ($ends + 1)){
    push (@active_ends, $previous_pos);
  }else{
    print "\nIncorrect number of start/end coordinates!!\n\n";
    exit();
  }

  #Now build an ordered hash of start/end positions 
  my $region_count = 0;
  foreach my $start (@active_starts){
    $region_count++;
    my $end = shift @active_ends;

    #Skip very small regions (say less than 5 bp) 
    my $size = (abs($end-$start))+1;
    unless ($size <= 5){
      $current_active_region_id++;
      $active_regions{$region_count}{active_region_id} = "$chr_filter"."_"."$region_number"."_"."$current_active_region_id";
      $active_regions{$region_count}{chr_start} = $start;
      $active_regions{$region_count}{chr_end} = $end;
      $active_regions{$region_count}{size} = $size;
    }
  }

  #B.) Silent regions next
  my $previous_pos2;
  $first = 1;
  my @silent_starts;
  my @silent_ends;

  #Get the start stop positions
  foreach my $pos (@silent_coords_sort){
    if ($first == 1){
      push(@silent_starts, $pos);
      $first = 0;
      $previous_pos2 = $pos;
      next();
    }
    if ($pos > ($previous_pos2 + 1)){
      push(@silent_ends, $previous_pos2);
      push(@silent_starts, $pos);
    }
    $previous_pos2 = $pos;
  }
  #If the final end position did not get added, add it now
  $starts = @silent_starts;
  $ends = @silent_ends;
  if ($starts == $ends){
    #do nothing
  }elsif($starts == ($ends + 1)){
    push (@silent_ends, $previous_pos2);
  }else{
    print RED, "\nIncorrect number of start/end coordinates!!\n\n", RESET;
    exit();
  }

  #Now build an ordered hash of start/end positions 
  $region_count = 0;
  foreach my $start (@silent_starts){
    $region_count++;
    my $end = shift @silent_ends;

    #Skip very small regions (say less than 5 bp) 
    my $size = (abs($end-$start))+1;
    unless ($size <= 5){
      $current_silent_region_id++;
      $silent_regions{$region_count}{silent_region_id} = "$chr_filter"."_"."$region_number"."_"."$current_silent_region_id";
      $silent_regions{$region_count}{chr_start} = $start;
      $silent_regions{$region_count}{chr_end} = $end;
      $silent_regions{$region_count}{size} = $size;
    }
  }


  #Now go through the active and silent regions
  #i.) Calculate their gene-level coordinates from their chromosome coordinates
  my $gene_id = $introns{$intron_id}{gene_id};

  foreach my $ar (keys %active_regions){
    if ($gene_id =~ /\d+/){
      my $start = $active_regions{$ar}{chr_start};
      my $end = $active_regions{$ar}{chr_end};
      my $coords_ref = &convertGenomicCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$start, '-end_pos'=>$end, '-ordered'=>"yes");
      $active_regions{$ar}{start} = $coords_ref->{$gene_id}->{gene_start};
      $active_regions{$ar}{end} = $coords_ref->{$gene_id}->{gene_end};
    }else{
      $active_regions{$ar}{start} = "NA";
      $active_regions{$ar}{end} = "NA";
    }
  }

  foreach my $sr (keys %silent_regions){
    if ($gene_id =~ /\d+/){
      my $start = $silent_regions{$sr}{chr_start};
      my $end = $silent_regions{$sr}{chr_end};
      my $coords_ref = &convertGenomicCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$start, '-end_pos'=>$end, '-ordered'=>"yes");
      $silent_regions{$sr}{start} = $coords_ref->{$gene_id}->{gene_start};
      $silent_regions{$sr}{end} = $coords_ref->{$gene_id}->{gene_end};
    }else{
      $silent_regions{$sr}{start} = "NA";
      $silent_regions{$sr}{end} = "NA";
    }
  }


  #ii.) Determine which intron they belong to and assign them a name according to their order within that intron
  #     - e.g. I2_SR1 (intron2-silentRegion1), I2_SR2, I2_AR1 (intron2-activeRegion1)
  #     - If an intron corresponds to multiple genes, the ID will look like: Ix_SR1, Ix_SR2
  #iii.) For each intron, add up the active and silent bases

  #A.) First for introns that are contained within a single gene... use gene based coordinates
  if ($introns{$intron_id}{contained_within_gene_count} eq "1"){

    my $active_base_count = 0;
    my $silent_base_count = 0;
    my $intron_ar_count = 0;
    my $intron_sr_count = 0;

    foreach my $ar (sort {$active_regions{$a}->{start} <=> $active_regions{$b}->{start}} keys %active_regions){
      $active_base_count += $active_regions{$ar}{size};
      $active_regions{$ar}{intron_id} = $intron_id;
      $intron_ar_count++;
      my $ar_name = "$introns{$intron_id}{name}"."_AR"."$intron_ar_count";
      $active_regions{$ar}{name} = $ar_name;
    }

    foreach my $sr (sort {$silent_regions{$a}->{start} <=> $silent_regions{$b}->{start}} keys %silent_regions){
      $silent_base_count += ($silent_regions{$sr}{size});
      $silent_regions{$sr}{intron_id} = $intron_id;
      $intron_sr_count++;
      my $sr_name = "$introns{$intron_id}{name}"."_SR"."$intron_sr_count";
      $silent_regions{$sr}{name} = $sr_name;
    }

    $introns{$intron_id}{active_base_count} = $active_base_count;
    $introns{$intron_id}{silent_base_count} = $silent_base_count+2; 

  }else{
    #B.) Deal with introns that have overlaping genes or that are actually intergenic - simply use chromosome coordinates for these (instead of gene coordinates)
    my $active_base_count = 0;
    my $silent_base_count = 0;
    my $intron_ar_count = 0;
    my $intron_sr_count = 0;

    foreach my $ar (sort {$active_regions{$a}->{chr_start} <=> $active_regions{$b}->{chr_start}} keys %active_regions){
      $active_base_count += $active_regions{$ar}{size};
      $active_regions{$ar}{intron_id} = $intron_id;
      $intron_ar_count++;
      my $ar_name = "$introns{$intron_id}{name}"."_AR"."$intron_ar_count";
      $active_regions{$ar}{name} = $ar_name;
    }

    foreach my $sr (sort {$silent_regions{$a}->{chr_start} <=> $silent_regions{$b}->{chr_start}} keys %silent_regions){
      $silent_base_count += ($silent_regions{$sr}{size});
      $silent_regions{$sr}{intron_id} = $intron_id;
      $intron_sr_count++;
      my $sr_name = "$introns{$intron_id}{name}"."_SR"."$intron_sr_count";
      $silent_regions{$sr}{name} = $sr_name;
    }

    $introns{$intron_id}{active_base_count} = $active_base_count;
    $introns{$intron_id}{silent_base_count} = $silent_base_count+2; 
  }

  $introns{$intron_id}{active_regions} = \%active_regions;
  $introns{$intron_id}{silent_regions} = \%silent_regions;
  my $active_region_count = keys %active_regions;
  my $silent_region_count = keys %silent_regions;
  $introns{$intron_id}{active_region_count} = $active_region_count;
  $introns{$intron_id}{silent_region_count} = $silent_region_count;

  $grand_active_region_count += $active_region_count;
  $grand_silent_region_count += $silent_region_count;

  #NOTE: Since very small active/silent regions (< 5 bp) are skipped, the total silent and active region base counts may not always add up to the total intron region size 

  $| = 1; print BLUE, "\n\tActiveRegions: $active_region_count\tSilentRegions: $silent_region_count\t", RESET; $| = 0;
  print LOG "\n\tActiveRegions: $active_region_count\tSilentRegions: $silent_region_count\t";

  return();
}


################################################################################################
#Determine mRNA/EST sequence support for each exon region                                      #
################################################################################################
sub determineSequenceSupport{
  my %args = @_;

  #Go through each intron, activeIntronRegion and silentIntronRegion and determine if is completely contained within an xmRNA or xEST exon
  #Test both strands for each region...

  #A.) First initialize the variables
  foreach my $i_id (sort {$introns{$a}->{start_chr} <=> $introns{$b}->{start_chr}} keys %introns){
    #Skip failing introns and intergenic regions!!
    unless (($introns{$i_id}{pass} == 1) && $introns{$i_id}{gene_id} =~ /\d+/){
      next();
    }
    $introns{$i_id}{supporting_mrna_count} = 0;
    $introns{$i_id}{supporting_est_count} = 0;
    $introns{$i_id}{supporting_xmrna_count} = 0;
    $introns{$i_id}{supporting_xest_count} = 0;
    $introns{$i_id}{conserved_species_count} = 0;
    my $active_regions_ref = $introns{$i_id}{active_regions};
    my $silent_regions_ref = $introns{$i_id}{silent_regions};
    #Active regions first
    foreach my $ar (sort {$active_regions_ref->{$a}->{chr_start} <=> $active_regions_ref->{$b}->{chr_start}} keys %{$active_regions_ref}){
      $active_regions_ref->{$ar}->{supporting_mrna_count} = 0;
      $active_regions_ref->{$ar}->{supporting_est_count} = 0;
      $active_regions_ref->{$ar}->{supporting_xmrna_count} = 0;
      $active_regions_ref->{$ar}->{supporting_xest_count} = 0;
      $active_regions_ref->{$ar}->{conserved_species_count} = 0;
    }
    #Silent regions next
    foreach my $sr (sort {$silent_regions_ref->{$a}->{chr_start} <=> $silent_regions_ref->{$b}->{chr_start}} keys %{$silent_regions_ref}){
      $silent_regions_ref->{$sr}->{supporting_mrna_count} = 0;
      $silent_regions_ref->{$sr}->{supporting_est_count} = 0;
      $silent_regions_ref->{$sr}->{supporting_xmrna_count} = 0;
      $silent_regions_ref->{$sr}->{supporting_xest_count} = 0;
      $silent_regions_ref->{$sr}->{conserved_species_count} = 0;
    }
  } 

  #B.) Next test for overlaping xmRNA and xEST coordinates on both strands for each region
  my %seqs;
  $seqs{mrna} = $mrna_exons;
  $seqs{est} = $est_exons;
  $seqs{xmrna} = $xmrna_exons;
  $seqs{xest} = $xest_exons;

  foreach my $i_id (sort {$introns{$a}->{start_chr} <=> $introns{$b}->{start_chr}} keys %introns){
    #Skip failing introns and intergenic regions!!
    unless (($introns{$i_id}{pass} == 1) && $introns{$i_id}{gene_id} =~ /\d+/){
      next();
    }

    #If this intron corresponds to a single gene only, then only look for overlap on the matching strand, otherwise test both strands 
    my @chr_strands;
    if ($introns{$i_id}{contained_within_gene_count} eq "1"){
      my $strand = $introns{$i_id}{strand};
      if ($strand eq "1"){
        $strand = "+";
      }else{
        $strand = "-";
      }
      my $chr_strand = "chr$chr_filter"."_"."$strand";
      push(@chr_strands, $chr_strand);
    }else{
      my $chr_strand1 = "chr$chr_filter"."_"."+";
      my $chr_strand2 = "chr$chr_filter"."_"."-";
      push(@chr_strands, $chr_strand1);
      push(@chr_strands, $chr_strand2);
    }

    #Go through each seq type (mrna, est, xmrna, and xest) and test for overlaps with that type
    my $i_start_chr = $introns{$i_id}{start_chr};
    my $i_end_chr = $introns{$i_id}{end_chr};
  
    &testSeqOverlap('-id'=>$i_id, '-regions_ref'=>\%introns, '-strands'=>\@chr_strands, '-seqs_ref'=>\%seqs, '-start_chr'=>$i_start_chr, '-end_chr'=>$i_end_chr);

    my $active_regions_ref = $introns{$i_id}{active_regions};
    my $silent_regions_ref = $introns{$i_id}{silent_regions};
    
    #Active regions first
    foreach my $ar (sort {$active_regions_ref->{$a}->{chr_start} <=> $active_regions_ref->{$b}->{chr_start}} keys %{$active_regions_ref}){
      my $ar_start_chr = $active_regions_ref->{$ar}->{chr_start};
      my $ar_end_chr = $active_regions_ref->{$ar}->{chr_end};
      &testSeqOverlap('-id'=>$ar, '-regions_ref'=>$active_regions_ref, '-strands'=>\@chr_strands, '-seqs_ref'=>\%seqs, '-start_chr'=>$ar_start_chr, '-end_chr'=>$ar_end_chr);
    }
    #Silent regions next
    foreach my $sr (sort {$silent_regions_ref->{$a}->{chr_start} <=> $silent_regions_ref->{$b}->{chr_start}} keys %{$silent_regions_ref}){
      my $sr_start_chr = $silent_regions_ref->{$sr}->{chr_start};
      my $sr_end_chr = $silent_regions_ref->{$sr}->{chr_end};
      &testSeqOverlap('-id'=>$sr, '-regions_ref'=>$silent_regions_ref, '-strands'=>\@chr_strands, '-seqs_ref'=>\%seqs, '-start_chr'=>$sr_start_chr, '-end_chr'=>$sr_end_chr);
    }
  } 

  return();
}

########################################
#Helper function for previous function #
########################################
sub testSeqOverlap{
  my %args = @_;
  my $id = $args{'-id'};
  my $regions_ref = $args{'-regions_ref'};
  my @chr_strands = @{$args{'-strands'}};
  my $seqs_ref = $args{'-seqs_ref'};
  my $start_chr = $args{'-start_chr'};
  my $end_chr = $args{'-end_chr'};

  #Keep a list of all species IDs encountered for the current junction_id
  my %species_list;

  #Go through each seq type (mrna, est, xmrna, xest)
  foreach my $type (keys %{$seqs_ref}){
    my $data_ref = $seqs_ref->{$type};
    my $field_name = "supporting_"."$type"."_count";

    #Make sure same some seq exons were defined 
    my $data_found = 0;
    foreach my $chr_strand (@chr_strands){
      if ($data_ref->{$chr_strand}){
        $data_found = 1;
      }
    }
    unless($data_found){
      $regions_ref->{$id}->{$field_name} = "na";
      next();
    }

    #For OTHER SPECIES - allow some wiggle room in the coordinates of the boundary by making it a bit smaller
    my $start_test = $start_chr;
    my $end_test = $end_chr;
    my $x_flag = 0;
    if ($type eq "xmrna" || $type eq "xest"){
      $start_test += $wiggle;
      $end_test -= $wiggle;
      $x_flag = 1;
    }

    #Test for overlap on both strands
    foreach my $chr_strand (@chr_strands){
      if ($data_ref->{$chr_strand}){
        my $seq_exons_ref = $data_ref->{$chr_strand};

        foreach my $seq_exon (keys %{$seq_exons_ref}){

          my $string = $seq_exon;
          if ($string =~ /(\d+)\_(\d+)/){
            my $seq_exon_start = $1;
            my $seq_exon_end = $2;

            if (($start_test >= $seq_exon_start) && ($start_test <= $seq_exon_end) && ($end_test >= $seq_exon_start) && ($end_test <= $seq_exon_end)){
              $regions_ref->{$id}->{$field_name} += ($seq_exons_ref->{$seq_exon}->{c});

              if ($x_flag){
                my $s_ref = $seq_exons_ref->{$seq_exon}->{s};
                foreach my $s (keys %{$s_ref}){
                  $species_list{$s}='';
                }
              }
            }
          }
        }
      }
    }

    #Count up all the species encountered in boundary matches from xmrnas and xests
    my $conserved_species_count = keys %species_list;
    $regions_ref->{$id}->{conserved_species_count} = $conserved_species_count;
  }
  return();
}


################################################################################################
#Print four output files                                                                       #
################################################################################################
sub printOutputFiles{
  my %args = @_;
  
  my $intron_db_file = "$outdir"."chr$chr_filter"."_"."$region_number"."_introns_annotated.txt";
  my $fasta_file = "$outdir"."chr$chr_filter"."_"."$region_number"."_introns.fa";
  my $intron_regions_active_file = "$outdir"."chr$chr_filter"."_"."$region_number"."_intronRegionsActive.txt";
  my $intron_regions_silent_file = "$outdir"."chr$chr_filter"."_"."$region_number"."_intronRegionsSilent.txt";

  open (INTRON_DB, ">$intron_db_file") || die "\nCould not open intron db file: $intron_db_file\n\n";
  open (FASTA, ">$fasta_file") || die "\nCould not open fasta file: $fasta_file\n\n";
  open (ACTIVE, ">$intron_regions_active_file") || die "\nCould not open active intron regions file: $intron_regions_active_file\n\n";
  open (SILENT, ">$intron_regions_silent_file") || die "\nCould not open silent intron regions file: $intron_regions_silent_file\n\n";

  print INTRON_DB "Intron_ID\tGene_ID\tEnsEMBL_Gene_ID\tGene_Name\tSeq_Name\tChromosome\tStrand\tUnit1_start\tUnit1_end\tUnit1_start_chr\tUnit1_end_chr\tBase_Count\tUnMasked_Base_Count\tCoding_Base_Count\tActive_Region_Count\tActive_Base_Count\tSilent_Region_Count\tSilent_Base_Count\tMultiple_Genes\tGene_ID_List\tSupporting_mRNA_Count\tSupporting_EST_Count\tSupporting_xmRNA_Count\tSupporting_xEST_Count\tConserved_Species_Count\n";

  print ACTIVE "Active_Region_ID\tIntron_ID\tGene_ID\tEnsEMBL_Gene_ID\tGene_Name\tSeq_Name\tChromosome\tStrand\tUnit1_start\tUnit1_end\tUnit1_start_chr\tUnit1_end_chr\tBase_Count\tUnMasked_Base_Count\tCoding_Base_Count\tSupporting_mRNA_Count\tSupporting_EST_Count\tSupporting_xmRNA_Count\tSupporting_xEST_Count\tConserved_Species_Count\n";
  print SILENT "Silent_Region_ID\tIntron_ID\tGene_ID\tEnsEMBL_Gene_ID\tGene_Name\tSeq_Name\tChromosome\tStrand\tUnit1_start\tUnit1_end\tUnit1_start_chr\tUnit1_end_chr\tBase_Count\tUnMasked_Base_Count\tCoding_Base_Count\tSupporting_mRNA_Count\tSupporting_EST_Count\tSupporting_xmRNA_Count\tSupporting_xEST_Count\tConserved_Species_Count\n";

  foreach my $i_id (sort {$introns{$a}->{start_chr} <=> $introns{$b}->{start_chr}} keys %introns){

    #Skip failing introns and intergenic regions!!
    unless (($introns{$i_id}{pass} == 1) && $introns{$i_id}{gene_id} =~ /\d+/){
      next();
    }

    my $ensembl_g_id;
    my $gene_name;
    my $gene_id = $introns{$i_id}{gene_id};
    my $masked_gene_seq = $genes_ref_masked->{$gene_id}->{sequence};
    my $masked_intron_seq = substr ($masked_gene_seq, $introns{$i_id}{start}, $introns{$i_id}{size});
    my $strand = $introns{$i_id}{strand};
    
    if ($introns{$i_id}{contained_within_gene_count} eq "1"){
      $ensembl_g_id = $genes_ref->{$gene_id}->{ensembl_g_id};
      $gene_name = $genes_ref->{$gene_id}->{gene_name};
    }else{
      $ensembl_g_id = "NA";
      $gene_name = "NA";
      $strand = "NA";
    }

    #3-A.) Print out intron sequence records to a DB file
    #    - Intron_ID, Gene_ID, EnsEMBL_Gene_ID, Gene_Name, Intron_Name, Chromosome, Strand, Unit1_start, Unit1_end, Unit1_start_chr, Unit1_end_chr, Length, UnMasked_Base_Count, Active_Base_Count, Silent_Base_Count, MultiGene, Gene_ID_list
    my $chr_intron_id = "$chr_filter"."_"."$region_number"."_"."$i_id";
    print INTRON_DB "$chr_intron_id\t$gene_id\t$ensembl_g_id\t$gene_name\t$introns{$i_id}{name}\t$chr_filter\t$introns{$i_id}{strand}\t$introns{$i_id}{start}\t$introns{$i_id}{end}\t$introns{$i_id}{start_chr}\t$introns{$i_id}->{end_chr}\t$introns{$i_id}{size}\t$introns{$i_id}{unmasked_base_count}\t0\t$introns{$i_id}{active_region_count}\t$introns{$i_id}{active_base_count}\t$introns{$i_id}{silent_region_count}\t$introns{$i_id}{silent_base_count}\t$introns{$i_id}{multi_gene}\t@{$introns{$i_id}{gene_id_list}}\t$introns{$i_id}{supporting_mrna_count}\t$introns{$i_id}{supporting_est_count}\t$introns{$i_id}{supporting_xmrna_count}\t$introns{$i_id}{supporting_xest_count}\t$introns{$i_id}{conserved_species_count}\n";

    #3-B.) Print out a fasta file of the intron sequences
    #    i.e. ">$intron_id\n$intron_sequence\n"
    print FASTA ">$chr_intron_id\n$masked_intron_seq\n";


    #3-C.) Print out DB files for 'IntronRegion_Active' and 'IntronRegion_Silent'
    #Active_Region_ID, Intron_ID, Gene_ID, EnsEMBL_Gene_ID, Gene_Name, Active_Region_Name, Chromosome, Strand, Unit1_start, Unit1_end, Unit1_start_chr, Unit1_end_chr, Length, UnMasked_Base_Count
    #Silent_Region_ID, Intron_ID, Gene_ID, EnsEMBL_Gene_ID, Gene_Name, Silent_Region_Name, Chromosome, Strand, Unit1_start, Unit1_end, Unit1_start_chr, Unit1_end_chr, Length, UnMasked_Base_Count
    my $active_regions_ref = $introns{$i_id}{active_regions};
    my $silent_regions_ref = $introns{$i_id}{silent_regions};

    #Active regions first
    foreach my $ar (sort {$active_regions_ref->{$a}->{chr_start} <=> $active_regions_ref->{$b}->{chr_start}} keys %{$active_regions_ref}){

      #Determine the number of unmasked bases within this active region...
      my $start = $active_regions_ref->{$ar}->{start};
      my $size = $active_regions_ref->{$ar}->{size};
      my $masked_region_seq = substr ($masked_gene_seq, $start, $size);
      my $masked_n_count = ($masked_region_seq =~ tr/N/N/);
      my $unmasked_count = $size - $masked_n_count;

      print ACTIVE "$active_regions_ref->{$ar}->{active_region_id}\t$chr_intron_id\t@{$introns{$i_id}{gene_id_list}}\t$ensembl_g_id\t$gene_name\t$active_regions_ref->{$ar}->{name}\t$chr_filter\t$strand\t$active_regions_ref->{$ar}->{start}\t$active_regions_ref->{$ar}->{end}\t$active_regions_ref->{$ar}->{chr_start}\t$active_regions_ref->{$ar}->{chr_end}\t$size\t$unmasked_count\t0\t$active_regions_ref->{$ar}->{supporting_mrna_count}\t$active_regions_ref->{$ar}->{supporting_est_count}\t$active_regions_ref->{$ar}->{supporting_xmrna_count}\t$active_regions_ref->{$ar}->{supporting_xest_count}\t$active_regions_ref->{$ar}->{conserved_species_count}\n";
    }

    #Silent regions next
    foreach my $sr (sort {$silent_regions_ref->{$a}->{chr_start} <=> $silent_regions_ref->{$b}->{chr_start}} keys %{$silent_regions_ref}){

      #Determine the number of unmasked bases within this silent region...
      my $start = $silent_regions_ref->{$sr}->{start};
      my $size = $silent_regions_ref->{$sr}->{size};
      my $masked_region_seq = substr ($masked_gene_seq, $start, $size);
      my $masked_n_count = ($masked_region_seq =~ tr/N/N/);
      my $unmasked_count = $size - $masked_n_count;

      print SILENT "$silent_regions_ref->{$sr}->{silent_region_id}\t$chr_intron_id\t@{$introns{$i_id}{gene_id_list}}\t$ensembl_g_id\t$gene_name\t$silent_regions_ref->{$sr}->{name}\t$chr_filter\t$strand\t$silent_regions_ref->{$sr}->{start}\t$silent_regions_ref->{$sr}->{end}\t$silent_regions_ref->{$sr}->{chr_start}\t$silent_regions_ref->{$sr}->{chr_end}\t$size\t$unmasked_count\t0\t$silent_regions_ref->{$sr}->{supporting_mrna_count}\t$silent_regions_ref->{$sr}->{supporting_est_count}\t$silent_regions_ref->{$sr}->{supporting_xmrna_count}\t$silent_regions_ref->{$sr}->{supporting_xest_count}\t$silent_regions_ref->{$sr}->{conserved_species_count}\n";
    }
  } 

  close(INTRON_DB);
  close(FASTA);
  close(ACTIVE);
  close(SILENT);
  
  return();
}
