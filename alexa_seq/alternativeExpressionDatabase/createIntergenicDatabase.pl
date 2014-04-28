#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to create a database of intergenic sequences for the purposes of read mapping
#These sequences should be repeat masked
#An 'intergenic region' in this context means the sequence between EnsEMBL genes (which are created by merging overlapping genes for an entire chromosome)
#Gene models are EnsEMBL transcripts/genes taken from an ALEXA database

#A chr_filter option is used to process one chromosome at a time for performance reasons

#1-A.) Get the complete masked chromosome sequence from EnsEMBL - or simply an EnsEMBL API object that allows their retrieval later via chromosome coordinates
#1-B.) Get gene models from ALEXA for all genes
#1-C.) Using the outer chromosome coordinates of all genes mask out genic regions in the chromosome
#1-D.) Get chromosome coordinates for the intergenic regions that remain for this chromosome

#2.) Go through UCSC mRNA AND EST databases and determine the coverage of each intergenic region.
#    - Using this coverage map, identify active and silent regions within each intergenic region
#    - For ESTs, apply basic quality filters, only used spliced ESTs, etc.
#    - Use this EST and mRNA data to define the regions of each intergenic region that have evidence for expression versus those that do not...
#    - i.e. identify the following two sets of regions: IntergenicRegions_Active vs IntergenicRegions_Silent

#3.) Print output files:

#3-A.) Print out intergenic sequence records to a DB file
#    - Intergenic_ID, Chromosome, Strand, Unit1_start, Unit1_end, Unit1_start_chr, Unit1_end_chr, Length, UnMasked_Base_Count, Active_Base_Count, Silent_Base_Count, Upstream_Gene_ID, Downstream_Gene_ID

#3-B.) Print out a fasta file of the intergenic sequences
#    i.e. ">$intergenic_id\n$intergenic_sequence\n"

#NOTE: Remember that in ALEXA, transcripts, exons, introns etc. are reverse complememented if they are on the negative strand (because source gene sequences are stored in ALEXA this way)
#      - However, these INTERGENIC sequences are all taken from complete chromosome sequences 

#3-C.) Print out DB files for 'IntergenicRegion_Active' and 'IntergenicRegion_Silent'
#      - Active_Region_ID, Intergenic_ID, Gene_ID, EnsEMBL_Gene_ID, Active_Region_Name, Chromosome, Strand, Unit1_start_chr, Unit1_end_chr, Length, UnMasked_Base_Count
#      - Silent_Region_ID, Intergenic_ID, Gene_ID, EnsEMBL_Gene_ID, Silent_Region_Name, Chromosome, Strand, Unit1_start_chr, Unit1_end_chr, Length, UnMasked_Base_Count


use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use Benchmark;
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
my $ensembl_api_version = ''; #Version of EnsEMBL to use
my $species = '';
my $ensembl_database = '';
my $ensembl_server = '';
my $ensembl_user = '';
my $ensembl_password = '';
my $alexa_database = '';
my $alexa_server = '';
my $alexa_user = '';
my $alexa_password = '';
my $min_size = '';
my $chr_filter = '';
my $ucsc_align_dir = '';
my $genbank_mapfile = '';
my $wiggle = '';
my $outdir = '';
my $logfile = '';
my $verbose = '';

GetOptions ('ensembl_api_version=s'=>\$ensembl_api_version, 'species=s'=>\$species,
            'ensembl_database=s'=>\$ensembl_database, 'ensembl_server=s'=>\$ensembl_server, 'ensembl_user=s'=>\$ensembl_user, 'ensembl_password=s'=>\$ensembl_password, 
            'alexa_database=s'=>\$alexa_database, 'alexa_server=s'=>\$alexa_server, 'alexa_user=s'=>\$alexa_user, 'alexa_password=s'=>\$alexa_password, 
            'min_size=i'=>\$min_size, 'chr_filter=s'=>\$chr_filter, 'ucsc_align_dir=s'=>\$ucsc_align_dir, 'genbank_mapfile=s'=>\$genbank_mapfile, 'wiggle=i'=>\$wiggle,
            'outdir=s'=>\$outdir, 'logfile=s'=>\$logfile, 'verbose=s'=>\$verbose);


#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the correct EnsEMBL API version using: --ensembl_api_version (41, 42, etc.)", RESET;
print GREEN, "\n\tSpecify the SOURCE Ensembl Database, Server, User and Password using: --ensembl_database  --ensembl_server  --ensembl_user  and  --ensembl_password", RESET;
print GREEN, "\n\tThe species using: --species (e.g. --species=Human or --species='Homo sapiens')", RESET;
print GREEN, "\n\t\tMake sure the species you supply matches the EnsEMBL database you supply!!", RESET;
print GREEN, "\n\tSpecify the TARGET Database and Server to query using: --alexa_database and --alexa_server", RESET;
print GREEN, "\n\tSpecify the User and Password for access using: --alexa_user and --alexa_password\n", RESET;
print GREEN, "\n\tSpecify the minimum intergenic region size (regions with less unmasked bases than this will be skipped) using:  --min_size", RESET;
print GREEN, "\n\tSpecify a single chromosome to be processed using: --chr_filter", RESET;
print GREEN, "\n\tSpecify the path to be used for output files using: --outdir", RESET;
print GREEN, "\n\tSpecify a directory containing UCSC mRNA/EST/xmRNA and xEST files using: --ucsc_align_dir", RESET;
print GREEN, "\n\tSpecify a berkeley DB containing GenBank-To-Species mappings using:  --genbank_mapfile", RESET;
print GREEN, "\n\tJunctions for xmRNA and xEST will be allowed to 'wiggle' by the amount specied by:  --wiggle", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of the design process using: --logfile", RESET;
print GREEN, "\n\tIf you want verbose output, use: --verbose=yes", RESET;

print GREEN, "\n\nExample: createIntergenicDatabase.pl  --ensembl_api_version=53  --species=Human  --ensembl_database=homo_sapiens_core_53_36o  --ensembl_server=ensembl01.bcgsc.ca  --ensembl_user=ensembl  --ensembl_password=ensembl  --alexa_database=ALEXA_hs_53_36o  --alexa_server=jango.bcgsc.ca  --alexa_user=viewer  --alexa_password=viewer  --min_size=42  --chr_filter='3:16:121020102-126260001'  --ucsc_align_dir=/projects/malachig/sequence_databases/hs_53_36o/mrna_est/partitions/  --genbank_mapfile=/projects/malachig/sequence_databases/hs_53_36o/mrna_est/GenBankToOrganism.btree  --wiggle=5  --outdir=/projects/malachig/sequence_databases/hs_53_36o/intergenics/temp/  --logfile=/projects/malachig/sequence_databases/hs_53_36o/logs/createIntergenicDatabase/createIntergenicDatabase_3_16_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($ensembl_api_version && $species && $ensembl_database && $ensembl_server && $ensembl_user && defined($ensembl_password) && $alexa_database && $alexa_server && $alexa_user && $alexa_password && $min_size && $chr_filter && $ucsc_align_dir && $genbank_mapfile && ($wiggle =~ /\d+/) && $outdir && $logfile){
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

my $tmp_chr;
if ($chr_filter eq "MT"){
  $tmp_chr = "chrM";
}else{
  $tmp_chr = "chr$chr_filter";
}
$ucsc_align_dir = &checkDir('-dir'=>$ucsc_align_dir, '-clear'=>"no");
my $ucsc_mrna_table = "$ucsc_align_dir"."$tmp_chr"."_mrna.txt.gz";
my $ucsc_est_table = "$ucsc_align_dir"."$tmp_chr"."_est.txt.gz";
my $ucsc_xmrna_table = "$ucsc_align_dir"."$tmp_chr"."_xmrna.txt.gz";
my $ucsc_xest_table = "$ucsc_align_dir"."$tmp_chr"."_xest.txt.gz";

#Load the GenBank-To-Species mapfile
my %gb_org_map;
tie(%gb_org_map, 'BerkeleyDB::Btree', -Cachesize =>256000000, -Filename=> $genbank_mapfile , -Flags => DB_RDONLY) or die "can't open file $genbank_mapfile: $! $BerkeleyDB::Error\n";


#Global time and process variables
my $t1 = new Benchmark;
my $pid = $$;

#Load the specified ensembl API
&loadEnsemblApi('-api'=>$ensembl_api_version);

#Check output dir
$outdir = &checkDir('-dir'=>$outdir, '-clear'=>"no");

#Set cache size for Berkley DB objects
my $cachesize = 25600000;


open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#Global counters
my $grand_base_count = 0;
my $grand_intergenic_seq_overlaps = 0;
my $grand_cumulative_coverage_bases = 0;
my $grand_unique_covered_bases = 0;
my $grand_active_region_count = 0;
my $grand_silent_region_count = 0;
my $grand_small_intergenics = 0;
my $grand_masked_intergenics = 0;

##CONNECT TO ENSEMBL SERVER
my $registry = "Bio::EnsEMBL::Registry";
$registry->load_registry_from_db(-host =>$ensembl_server, -user =>$ensembl_user, -pass=>$ensembl_password, -db_version =>$ensembl_api_version);

#Get a connection to the local Ensembl CORE database
my $ensembl_core_api = Bio::EnsEMBL::Registry->get_DBAdaptor($species, "core");

#Get a slice for the entire chromosome
my $slice_adaptor = $registry->get_adaptor($species, 'Core', 'Slice');

#If retrieving a slice fails at the chromsome level - try as a supercontig (neccessary for genes on NT_XXXXXX objects)
#Get the length of the chromosome
#my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr_filter, $start_filter, $end_filter);
#unless($slice){
#  $slice = $slice_adaptor->fetch_by_region('supercontig', $chr_filter, $start_filter, $end_filter);
#}
my $slice = $slice_adaptor->fetch_by_region('toplevel', $chr_filter, $start_filter, $end_filter);

my $masked_seq_object = $slice->get_repeatmasked_seq();
my $chr_length = length($masked_seq_object->seq());

if ($chr_filter eq "MT"){$chr_filter = "M";}

print BLUE, "\n\nSuccessfully retrieved a masked seq for chr: $chr_filter:$start_filter-$end_filter of length $chr_length\n\n", RESET;
&processUpdate();

#my $seq_file_name = "$outdir"."$chr_filter".".chr_masked_seq.btree";
#my $seq_rm_cmd = "rm -f $seq_file_name";
#system($seq_rm_cmd);
my %chr_seq;
$chr_seq{$chr_filter} = $masked_seq_object->seq();
&processUpdate();


#1.) Get gene models from ALEXA for all genes - get completed masked sequence for each gene locus

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$alexa_database, '-server'=>$alexa_server, '-user'=>$alexa_user, '-password'=>$alexa_password);
my $genes_ref;
my $current_intergenic_id = 0;
my $current_active_region_id = 0;
my $current_silent_region_id = 0;

my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All', '-chromosome'=>$chr_filter, '-range'=>"$start_filter-$end_filter")};

my %intergenics;
$| = 1; print BLUE, "\n\n1.) Getting basic gene info for chr: $chr_filter:$start_filter-$end_filter", RESET; $| = 0;
print LOG "\n\n1.) Getting basic gene info for chr: $chr_filter:$start_filter-$end_filter";
&getBasicGeneInfo ('-gene_ids'=> \@gene_ids);

#Close database connection
$alexa_dbh->disconnect();

#For both EST and mRNA seq files, write a berkley DB containing only those records corresponding to the current chromosome
my $mrna_ref;
my $est_ref;
$mrna_ref = &parseSeqFile('-infile'=>$ucsc_mrna_table);
&processUpdate();
$est_ref = &parseSeqFile('-infile'=>$ucsc_est_table);
&processUpdate();
my $mrna_record_count = keys %{$mrna_ref};
my $est_record_count = keys%{$est_ref};

#Use a different function to get the coordinates of mrna, est, xmrna and xest exons
my $mrna_exons = &importUcscAlignments('-import_type'=>"exon", '-align_file'=>$ucsc_mrna_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"mrna", '-filter'=>"0", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);
my $est_exons = &importUcscAlignments('-import_type'=>"exon", '-align_file'=>$ucsc_est_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"est", '-filter'=>"1", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);
my $xmrna_exons = &importUcscAlignments('-import_type'=>"exon", '-align_file'=>$ucsc_xmrna_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"xmrna", '-filter'=>"2", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);
my $xest_exons = &importUcscAlignments('-import_type'=>"exon", '-align_file'=>$ucsc_xest_table, '-gb_org_map'=>\%gb_org_map, '-seq_type'=>"xest", '-filter'=>"2", '-chr_filter'=>$chr_filter, '-start_filter'=>$start_filter, '-end_filter'=>$end_filter);

#Create berkley DB objects to store data on a intergenic-by-intergenic basis
my %intergenic_coverage;
my $cov_file_name = "$outdir"."$chr_filter"."_"."$region_number".".intergenic_coverage.btree";
my $cov_rm_cmd = "rm -f $cov_file_name";
system($cov_rm_cmd);
tie(%intergenic_coverage, 'BerkeleyDB::Btree', -Cachesize => $cachesize, -Filename=> $cov_file_name , -Flags => DB_CREATE) or die "can't open file $cov_file_name: $! $BerkeleyDB::Error\n";

my %active_coords;
my $ac_file_name = "$outdir"."$chr_filter"."_"."$region_number".".active_coords.btree";
my $ac_rm_cmd = "rm -f $ac_file_name";
system($ac_rm_cmd);
tie(%active_coords, 'BerkeleyDB::Btree', -Cachesize => $cachesize, -Filename=> $ac_file_name , -Flags => DB_CREATE) or die "can't open file $ac_file_name: $! $BerkeleyDB::Error\n";

my %silent_coords;
my $sc_file_name = "$outdir"."$chr_filter"."_"."$region_number".".silent_coords.btree";
my $sc_rm_cmd = "rm -f $sc_file_name";
system($sc_rm_cmd);
tie(%silent_coords, 'BerkeleyDB::Btree', -Cachesize => $cachesize, -Filename=> $sc_file_name , -Flags => DB_CREATE) or die "can't open file $sc_file_name: $! $BerkeleyDB::Error\n";


#Now process intergenic-by-intergenic
my $count_down = keys %intergenics;

$| = 1; print BLUE, "\n\n2.) Initializing coverage hash, adding coverage from mRNAs/ESTs, determining 'Active' and 'Silent regions.  Intergenic-By-Intergenic\n", RESET; $| = 0;
print LOG "\n\n2.) Initializing coverage hash, adding coverage from mRNAs/ESTs, determining 'Active' and 'Silent regions.  Intergenic-By-Intergenic\n";

foreach my $i_id (sort {$intergenics{$a}->{start_chr} <=> $intergenics{$b}->{start_chr}} keys %intergenics){
  $count_down--;

  my $strand = $intergenics{$i_id}{strand};
  $| = 1; print BLUE, "\n\n($count_down) Processing Intergenic Region $i_id\tSize $intergenics{$i_id}{size} bp ($intergenics{$i_id}{unmasked_base_count} bp unmasked)\tchr$chr_filter:$intergenics{$i_id}{start_chr}-$intergenics{$i_id}{end_chr} (strand: $strand)", RESET; $| = 0;
  print LOG "\n\n($count_down) Processing Intergenic Region $i_id\tSize $intergenics{$i_id}{size} bp ($intergenics{$i_id}{unmasked_base_count} bp unmasked)\tchr$chr_filter:$intergenics{$i_id}{start_chr}-$intergenics{$i_id}{end_chr} (strand: $strand)";
 

  #Skip failing intergenics regions!!
  unless (($intergenics{$i_id}{pass} == 1)){
    print YELLOW, "\n\tFailed - skipping", RESET;
    print LOG "\n\tFailed - skipping";
    next();
  }

  #2.) Go through UCSC mRNA AND EST databases and determine the coverage of each intergenic region.
  %intergenic_coverage = ();

  #2-A.) Initialize a coverage hash for all intergenic base positions (keyed as geneID_position)
  &initializeCoverage('-i_id'=>$i_id);  #&processUpdate('-display_option'=>2);

  #2-B.) Add coverage from ESTs/mRNAs to this hash
  if ($mrna_record_count > 0){
    &addCoverage('-alignment_object'=>$mrna_ref, '-i_id'=>$i_id, '-type'=>'mRNA');  #&processUpdate('-display_option'=>2);
  }
  if ($est_record_count > 0){
    &addCoverage('-alignment_object'=>$est_ref, '-i_id'=>$i_id, '-type'=>'EST');  #&processUpdate('-display_option'=>2);
  }

  #2-C.) Extract the coordinates corresponding to 'Active' and 'Silent' regions within each intergenic region
  &identifyActiveSilentRegions('-i_id'=>$i_id);
  &processUpdate('-display_option'=>1);

  #2-D.) Determine the conservation support
  &determineSequenceSupport();

}

#3.) Print output files
$| = 1; print BLUE, "\n\n\n3.) Printing output files to: $outdir\n", RESET; $| = 0;
print LOG "\n\n\n3.) Printing output files to: $outdir\n";
&printOutputFiles();

#Untie and delete berkley DB files
untie(%intergenic_coverage);
system($cov_rm_cmd);
untie(%active_coords);
system($ac_rm_cmd);
untie(%silent_coords);
system($sc_rm_cmd);

untie(%gb_org_map);

my $intergenic_count = keys %intergenics;

print BLUE, "\n\nFINAL SUMMARY - chr$chr_filter", RESET;
print BLUE, "\nIntergenic regions found = $intergenic_count", RESET;
print BLUE, "\nIntron base count = $grand_base_count", RESET;
print BLUE, "\nIntron-To-Sequence (EST/mRNA) Overlaps = $grand_intergenic_seq_overlaps", RESET;
print BLUE, "\nCumulative EST/mRNA coverage bases = $grand_cumulative_coverage_bases", RESET;
print BLUE, "\nUnique EST/mRNA covered bases = $grand_unique_covered_bases", RESET;
print BLUE, "\nActive region count (EST/mRNA evidence) = $grand_active_region_count", RESET;
print BLUE, "\nSilent region count (EST/mRNA evidence) = $grand_silent_region_count", RESET;
print BLUE, "\nNumber of intergenics skipped because they were too small = $grand_small_intergenics", RESET;
print BLUE, "\nNumber of intergenics skipped because they had too few unmasked bases = $grand_masked_intergenics\n\n", RESET;
print BLUE, "\nSCRIPT COMPLETE\n\n", RESET;

print LOG "\n\nFINAL SUMMARY - chr$chr_filter";
print LOG "\nIntergenic regions found = $intergenic_count";
print LOG "\nIntron base count = $grand_base_count";
print LOG "\nIntron-To-Sequence (EST/mRNA) Overlaps = $grand_intergenic_seq_overlaps";
print LOG "\nCumulative EST/mRNA coverage bases = $grand_cumulative_coverage_bases";
print LOG "\nUnique EST/mRNA covered bases = $grand_unique_covered_bases";
print LOG "\nActive region count (EST/mRNA evidence) = $grand_active_region_count";
print LOG "\nSilent region count (EST/mRNA evidence) = $grand_silent_region_count";
print LOG "\nNumber of intergenics skipped because they were too small = $grand_small_intergenics";
print LOG "\nNumber of intergenics skipped because they had too few unmasked bases = $grand_masked_intergenics\n\n";
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
  $| = 1; print BLUE, "\n\nGetting gene models", RESET; $| = 0; 
  print LOG "\nGetting gene models";
  &processUpdate();
  $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

  #Identify all genic regions by merging overlaping genes (using there gene coordinates)
  $| = 1; print BLUE, "\n\nUsing chromosome coordinates to determine chromosome-region genic regions", RESET; $| = 0; 
  print LOG "\nUsing chromosome coordinates to determine chromosome-region genic regions";
  &processUpdate();
  my %genic_regions;
  my $gr_count = 0;

  #First make sure the chr_start is smaller than chr_end
  foreach my $g (keys %{$genes_ref}){
    my $tmp = $genes_ref->{$g}->{chr_start};
    if ($genes_ref->{$g}->{chr_start} > $genes_ref->{$g}->{chr_end}){
      $genes_ref->{$g}->{chr_start} = $genes_ref->{$g}->{chr_end};
      $genes_ref->{$g}->{chr_end} = $tmp;
    }
  }

  my $count = 0;
  foreach my $g (sort {$genes_ref->{$a}->{chr_start} <=> $genes_ref->{$b}->{chr_start}} keys %{$genes_ref}){
    my $g_start = $genes_ref->{$g}->{chr_start};
    my $g_end = $genes_ref->{$g}->{chr_end};

    $count++;
    if ($count == 10){
      $count = 0;
      $| = 1; print BLUE, ".", RESET;  $| = 0;
    }

    #For this gene, go through all genic regions thus far and look for overlaps
    my $overlap = 0;
    foreach my $gr (sort {$genic_regions{$a}->{chr_start} <=> $genic_regions{$b}->{chr_start}} keys %genic_regions){
      my $gr_start = $genic_regions{$gr}{chr_start};
      my $gr_end = $genic_regions{$gr}{chr_end};

      #Test overlap
      if (($g_start >= $gr_start && $g_start <= $gr_end) || ($g_end >= $gr_start && $g_end <= $gr_end) || ($g_start <= $gr_start && $g_end >= $gr_end)){
        $overlap = 1;
      }

      if ($overlap == 1){
        #Merge this gene into an existing genic region and get the adjusted coordinates
        my @coords = ($g_start, $g_end, $gr_start, $gr_end);
        my @coords_sort = sort {$a <=> $b} @coords;
        $genic_regions{$gr_count}{chr_start} = $coords_sort[0];
        $genic_regions{$gr_count}{chr_end} = $coords_sort[3];
        last();
      }
    }

    if ($overlap == 0){
      #Enter a new genic region
      $gr_count++;
      $genic_regions{$gr_count}{chr_start} = $g_start;
      $genic_regions{$gr_count}{chr_end} = $g_end;
    }
  }

  #Now go through the entire length of the current chromosome region in 1 million bp blocks
  #If there is no genic region within 500k bases in either direction, create a dummy genic region of 2 bp to prevent extremely large intergenic regions from occuring
  #Otherwise massive regions occur (particular at the beginning and end of chromosomes)
  my $block_size = 1000000;
  
  for (my $i = $start_filter+$block_size; $i < $chr_length-1; $i+=$block_size){

    #Is there a gene within $block_size/2 from the current $i
    my $i_start = ($i - ($block_size/2));
    my $i_end = ($i + ($block_size/2));
    my $overlap = 0;
    foreach my $gr (sort {$genic_regions{$a}->{chr_start} <=> $genic_regions{$b}->{chr_start}} keys %genic_regions){
      my $gr_start = $genic_regions{$gr}{chr_start};
      my $gr_end = $genic_regions{$gr}{chr_end};

      if (($i_start >= $gr_start && $i_start <= $gr_end) || ($i_end >= $gr_start && $i_end <= $gr_end) || ($i_start <= $gr_start && $i_end >= $gr_end)){
        $overlap = 1;
        last();
      }
    }
    if ($overlap == 0){
      my $dummy_start = $i;
      my $dummy_end = $i+1;
      print YELLOW, "\n\tCreating a dummy genic region with coords: $dummy_start - $dummy_end", RESET;
      $gr_count++;
      $genic_regions{$gr_count}{chr_start} = $dummy_start;
      $genic_regions{$gr_count}{chr_end} = $dummy_end;
    }
  }

  my @genic_starts;
  my @genic_ends;
  foreach my $gr (sort {$genic_regions{$a}->{chr_start} <=> $genic_regions{$b}->{chr_start}} keys %genic_regions){
    push(@genic_starts, $genic_regions{$gr}{chr_start});
    push(@genic_ends, $genic_regions{$gr}{chr_end});
  }

  #Artificially create a genic region at the beginning and end of the chromosome
  unshift (@genic_starts, $start_filter-1);
  unshift (@genic_ends, $start_filter);
  push (@genic_starts, $end_filter-1);
  push (@genic_ends, $end_filter);

  #print YELLOW, "\n\nSTARTS: @genic_starts\n\nENDS: @genic_ends\n\n", RESET;

  my $starts = scalar(@genic_starts);
  my $ends = @genic_ends;
  if ($starts == $ends){
    #do nothing
  }else{
    print "\nIncorrect number of start/end coordinates!!\n\n";
    exit();
  }


  #Now identify the non-redundant intergenic regions for this chromosome region
  #To do this first identify the genic content coordinates, and then use these to extract the intergenic coordinates
  $| = 1; print BLUE, "\n\nGetting the start/end coordinates of intergenic regions using this chromosome-region coverage hash", RESET; $| = 0; 
  print LOG "\nGetting the start/end coordinates of intergenic regions using this chromosome-region coverage hash";
  &processUpdate();

  #Now build a final list of intergenic regions - perform initial filters on these before adding them in the list
  my $intergenic_count = 0;
  for (my $i = 0; $i < $starts-1; $i++){
    my $intergenic_start = ($genic_ends[$i])+1;
    my $intergenic_end = ($genic_starts[$i+1])-1;

    #Determine the size of each intergenic region
    my $size = ($intergenic_end-$intergenic_start)+1;

    $intergenic_count++;

    #print YELLOW, "\n\tIG: $intergenic_count\tStart: $intergenic_start\tEnd: $intergenic_end", RESET;

    $intergenics{$intergenic_count}{start_chr} = $intergenic_start;
    $intergenics{$intergenic_count}{end_chr} = $intergenic_end;
    $intergenics{$intergenic_count}{strand} = "1";  #Strand is always +ve for these
    $intergenics{$intergenic_count}{name} = "IG$intergenic_count";
    $intergenics{$intergenic_count}{size} = $size;
    $intergenics{$intergenic_count}{pass} = 1;
    $intergenics{$intergenic_count}{unmasked_base_count} = 0;
    $intergenics{$intergenic_count}{active_base_count} = 0;
    $intergenics{$intergenic_count}{silent_base_count} = 0;
    $intergenics{$intergenic_count}{upstream_gene} = "na";
    $intergenics{$intergenic_count}{downstream_gene} = "na";
  }

  my $chr_intergenic_count = keys %intergenics;

  $| = 1; print BLUE, "\n\nFound a total of $chr_intergenic_count intergenic regions for this chromosome region", RESET; $| = 0; 
  print LOG "\n\nFound a total of $chr_intergenic_count intergenic regions for this chromosome region";

  #Now go through each intergenic region and determine:
  #1.) The gene closest in the upstream direction
  #2.) The gene closest in the downstream direction
  $| = 1; print BLUE, "\n\nDetermining upstream and downstream genes", RESET; $| = 0; 
  print LOG "\n\nDetermining upstream and downstream genes";
  &processUpdate();

  $count = 0;
  foreach my $i_id (sort {$intergenics{$a}->{start_chr} <=> $intergenics{$b}->{start_chr}} keys %intergenics){
    $count++;
    if ($count == 10){
      $count = 0;
      $| = 1; print BLUE, ".", RESET;  $| = 0;
    }

    my $i_start = $intergenics{$i_id}{start_chr};
    my $i_end = $intergenics{$i_id}{end_chr};

    my $upstream_gene_id = "na";
    my $downstream_gene_id = "na";

    my $upstream_distance = 1000000000000000000000;
    my $downstream_distance = 1000000000000000000000;
    foreach my $gene_id (keys %{$genes_ref}){
      my $g_start = $genes_ref->{$gene_id}->{chr_start};
      my $g_end = $genes_ref->{$gene_id}->{chr_end};
      my $temp = $g_start;
      if ($g_start > $g_end){
        $temp = $g_start;
        $g_start = $g_end;
        $g_end = $temp;
      }

      #print YELLOW, "\nDEBUG: i_id: $i_id ($i_start - $i_end)\tg_id ($g_start - $g_end): $gene_id\tup_dis: $upstream_distance\tdown_dis: $downstream_distance", RESET;

      #Upstream
      my $distance = $i_start - $g_end;
      #print YELLOW, "\n\tDistance_U: $distance", RESET;
      if ($distance < $upstream_distance && $distance >= 0){
        $upstream_distance = $distance;
        $upstream_gene_id = $gene_id;
      }
      #Downstream
      $distance = $g_start - $i_end;
      #print YELLOW, "\n\tDistance_D: $distance", RESET;
      if ($distance < $downstream_distance && $distance >= 0){
        $downstream_distance = $distance;
        $downstream_gene_id = $gene_id;
      }
    }
    $intergenics{$i_id}{upstream_gene} = $upstream_gene_id;
    $intergenics{$i_id}{downstream_gene} = $downstream_gene_id;
  }

  #Now go through each intergenic regions and get the masked/unmasked base count
  $| = 1; print BLUE, "\n\nGetting # unmasked bases", RESET; $| = 0; 
  print LOG "\n\nGetting # unmasked bases";
  &processUpdate();

  $count = 0;
  foreach my $i_id (sort {$intergenics{$a}->{start_chr} <=> $intergenics{$b}->{start_chr}} keys %intergenics){
    $count++;
    if ($count == 10){
      $count = 0;
      $| = 1; print BLUE, ".", RESET;  $| = 0;
    }

    my $i_start = $intergenics{$i_id}{start_chr};
    my $i_end = $intergenics{$i_id}{end_chr};
    my $size = $intergenics{$i_id}{size};

    my $corrected_start = ($i_start - $start_filter)+1;

    #print YELLOW, "\nGet seq for: region = ($i_start - $i_end)\tcorrected_start = $corrected_start\tsize = $size", RESET;
    my $masked_n_count = 0;
    my $unmasked_count  = 0;
    if ($size > 0){
      $masked_n_count = (substr($chr_seq{$chr_filter}, $corrected_start, $size) =~ tr/N/N/);
      $unmasked_count = $size - $masked_n_count;
    }

    $intergenics{$i_id}{unmasked_base_count} = $unmasked_count;

    #print YELLOW, "\nI_ID: $i_id\tsize: $size\tunmasked_count: $unmasked_count", RESET;

    #Fail intergenic regions that are too small, or have too few unmasked bases to allow mapping of reads
    if ($size < $min_size){
      $grand_small_intergenics++;
      $intergenics{$i_id}{pass} = 0;
    }
    if ($size >= $min_size && $unmasked_count < $min_size){
      $grand_masked_intergenics++;
      $intergenics{$i_id}{pass} = 0;
    }
  }
  
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

  $| = 1; print BLUE, "\n\nGrabbing UCSC alignments that match the current chr from: $infile", RESET; $| = 0; 
  print LOG "\n\nGrabbing UCSC alignments that match the current chr from: $infile"; 

  unless (-e $infile){
    print YELLOW, "\n\nmRNA/EST file not found: $infile - returning 0 alignments for comparison", RESET;
    return(\%seqs);
  }

  open (ALIGN, "zcat $infile |") || die "\nCould not open mrna/est alignment file: $infile\n\n";

  while (<ALIGN>){
    chomp($_);
    #Get columns 10,15,20,22
    my @line = split("\t", $_);
    my $matches = $line[1];
    my $mismatches = $line[2];
    my $n_count = $line[4];
    my $query_insertions = $line[5];
    my $strand = $line[9];
    #my $acc = $line[10];
    my $chr = $line[14];

    my @exon_sizes = split(",", $line[19]);
    my @starts = split(",", $line[21]);
    my $exon_count = scalar(@starts);

    my $chr_test = "chr$chr_filter";

    #Unless the chromosome of the current EST/mRNA is the same as the current chromosome being processed, skip this record
    unless ($chr eq $chr_test){
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

    #$seqs{$record_count}{acc} = $acc;
    #$seqs{$record_count}{strand} = $strand;
    $seqs{$record_count}{starts} = \@starts;
    $seqs{$record_count}{sizes} = \@exon_sizes;
  }

  close(ALIGN);

  $| = 1; print BLUE, "\nFound a total of $record_count alignments for this chromosome, stored as a hash\n", RESET; $| = 0; 
  print LOG "\nFound a total of $record_count alignments for this chromosome, stored as a hash\n"; 

  return(\%seqs);
}

################################################################################################
#Initialize coverage hashes for each intergenic region                                         #
################################################################################################
sub initializeCoverage{
  my %args = @_;
  my $i_id = $args{'-i_id'};

  my $base_count = 0;

  #Construct this hash as keyed on "IntergenicID_ChromosomePosition".  This makes it easy to use a Berkeley DB to store this data structure if memory usage become a problem
  #For this coverage hash, inset the coordinates of each intergenic region at each end to ensure that we are not considering the last or first base of intergenic regions 

  my $chr_start = ($intergenics{$i_id}{start_chr})+1;
  my $chr_end = ($intergenics{$i_id}{end_chr})-1;

  for (my $i = $chr_start; $i <= $chr_end; $i++){
    my $i_pos_id = "$i_id"."_"."$i";
    $intergenic_coverage{$i} = 0;
    $base_count++;
    $grand_base_count++;
  }

  $| = 1; print BLUE, "\n\tBaseCount: $base_count\t", RESET; $| = 0; 
  print LOG "\n\tBaseCount: $base_count\t", RESET; 
  
  return();
}


################################################################################################
#Add coverage from mRNA/EST alignments to intergenic content block defined for EnsEMBL genes   #
################################################################################################
sub addCoverage{
  my %args = @_;
  my $seqs_ref = $args{'-alignment_object'};
  my $i_id = $args{'-i_id'};
  my $type = $args{'-type'};

  my $intergenic_seq_overlaps = 0;
  my $cumulative_coverage_bases = 0;
  my $unique_covered_bases = 0;

  foreach my $seq (keys %{$seqs_ref}){
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
    my $i_start_chr = $intergenics{$i_id}{start_chr};
    my $i_end_chr = $intergenics{$i_id}{end_chr};

    #Determine if there is overlap between this EST/mRNA and the current intergenic region
    my $overlap = 0;
    if ($seq_start > $i_start_chr && $seq_start < $i_end_chr){  #Start of seq is within intergenic
      $overlap = 1;
    }
    if ($seq_end > $i_start_chr && $seq_end < $i_end_chr){      #End of seq is within intergenic
      $overlap = 1;
    }
    if ($seq_start <= $i_start_chr && $seq_end >= $i_end_chr){  #Seq completely flanks intergenic
      $overlap = 1;
    }

    #If overlap is present, add whatever coverage corresponds to previously identified intergenic coordinates
    if ($overlap == 1){
      $grand_intergenic_seq_overlaps++;
      $intergenic_seq_overlaps++;

      #Go through each 'exon' from the EST/mRNA sequence and add coverage where the coordinates correspond to intergenic content blocks
      foreach my $seq_exon (keys %seq_exons_tmp){

        my $seq_exon_start = $seq_exons_tmp{$seq_exon}{start};
        my $seq_exon_end = $seq_exons_tmp{$seq_exon}{end};

        #print YELLOW, "\nSEQ_EXON: chr$chr_filter:$seq_exon_start-$seq_exon_end", RESET;

        for (my $i = $seq_exon_start; $i <= $seq_exon_end; $i++){

          #If this base position corresponds to an intergenic base position, increment it
          if (defined($intergenic_coverage{$i})){
            $grand_cumulative_coverage_bases++;
            $cumulative_coverage_bases++;
            $intergenic_coverage{$i}++;
          }
        }
      }
    }
  }

  #How many unique intergenic base positions were actually covered by 1 or mRNA/ESTs
  foreach my $i (keys %intergenic_coverage){
    if ($intergenic_coverage{$i} > 0){
      $grand_unique_covered_bases++;      
      $unique_covered_bases++;
    }
  }
  $| = 1; print BLUE, "\n\t$type Coverage.\tIntergenic-Seq-Overlaps: $intergenic_seq_overlaps\tCumulativeCoverageBases: $cumulative_coverage_bases\tUniqueCoveredBases: $unique_covered_bases\t", RESET; $| = 0;
  print LOG "\n\t$type Coverage.\tIntergenic-Seq-Overlaps: $intergenic_seq_overlaps\tCumulativeCoverageBases: $cumulative_coverage_bases\tUniqueCoveredBases: $unique_covered_bases\t";

  return();
}

########################################################################################################
#Extract the coordinates corresponding to 'Active' and 'Silent' regions within each intergenic region  #
########################################################################################################
sub identifyActiveSilentRegions{
  my %args = @_;
  my $i_id = $args{'-i_id'};

  #Now get all the intergenic positions specifically for this intergenic region
  foreach my $i (keys %intergenic_coverage){

    if ($intergenic_coverage{$i} > 0){
      $active_coords{$i}=1;
    }else{
      $silent_coords{$i}=1;
    }
  }

  #Now use the arrays of 'Active' and 'Silent' coords to identify the coordinates of 'Active' and 'Silent' regions
  my %active_regions;
  my %silent_regions;

  #A.) Active regions first
  my $previous_pos;
  my $first = 1;
  my @active_starts;
  my @active_ends;

  #Get the start stop positions
  foreach my $pos (sort {$a <=> $b} keys %active_coords){
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
  %active_coords = ();

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
  foreach my $pos (sort {$a<=>$b} keys %silent_coords){
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
  %silent_coords = ();

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

  #ii.) Assign each ACTIVE/SILENT region a name according to their order within the intergenic region
  #     - e.g. IG1_SR1 (intergenicRegion1-silentRegion1), IG2_SR2, IG2_AR1 (intergenicRegion2-activeRegion1)
  #iii.) For each intergenic, add up the active and silent bases
  my $active_base_count = 0;
  my $silent_base_count = 0;
  my $intergenic_ar_count = 0;
  my $intergenic_sr_count = 0;

  foreach my $ar (sort {$active_regions{$a}->{chr_start} <=> $active_regions{$b}->{chr_start}} keys %active_regions){
    $active_base_count += $active_regions{$ar}{size};
    $active_regions{$ar}{intergenic_id} = $i_id;
    $intergenic_ar_count++;
    my $ar_name = "$intergenics{$i_id}{name}"."_AR"."$intergenic_ar_count";
    $active_regions{$ar}{name} = $ar_name;
  }

  foreach my $sr (sort {$silent_regions{$a}->{chr_start} <=> $silent_regions{$b}->{chr_start}} keys %silent_regions){
    $silent_base_count += ($silent_regions{$sr}{size});
    $silent_regions{$sr}{intergenic_id} = $i_id;
    $intergenic_sr_count++;
    my $sr_name = "$intergenics{$i_id}{name}"."_SR"."$intergenic_sr_count";
    $silent_regions{$sr}{name} = $sr_name;
  }

  $intergenics{$i_id}{active_base_count} = $active_base_count;
  $intergenics{$i_id}{silent_base_count} = $silent_base_count+2; 

  $intergenics{$i_id}{active_regions} = \%active_regions;
  $intergenics{$i_id}{silent_regions} = \%silent_regions;
  my $active_region_count = keys %active_regions;
  my $silent_region_count = keys %silent_regions;
  $intergenics{$i_id}{active_region_count} = $active_region_count;
  $intergenics{$i_id}{silent_region_count} = $silent_region_count;

  $grand_active_region_count += $active_region_count;
  $grand_silent_region_count += $silent_region_count;

  #NOTE: Since very small active/silent regions (< 5 bp) are skipped, the total silent and active region base counts may not always add up to the total intergenic region size 
  $| = 1; print BLUE, "\n\tActiveRegions: $active_region_count\tSilentRegions: $silent_region_count\t", RESET; $| = 0;
  print LOG "\n\tActiveRegions: $active_region_count\tSilentRegions: $silent_region_count\t";

  return();
}


################################################################################################
#Determine mRNA/EST sequence support for each exon region                                      
################################################################################################
sub determineSequenceSupport{
  my %args = @_;

  #Go through each intron, activeIntronRegion and silentIntronRegion and determine if is completely contained within an xmRNA or xEST exon
  #Test both strands for each region...

  #A.) First initialize the variables
  foreach my $i_id (sort {$intergenics{$a}->{start_chr} <=> $intergenics{$b}->{start_chr}} keys %intergenics){
    #Skip failing intergenics and intergenic regions!!
    unless (($intergenics{$i_id}{pass} == 1)){
      next();
    }
    $intergenics{$i_id}{supporting_mrna_count} = 0;
    $intergenics{$i_id}{supporting_est_count} = 0;
    $intergenics{$i_id}{supporting_xmrna_count} = 0;
    $intergenics{$i_id}{supporting_xest_count} = 0;
    $intergenics{$i_id}{conserved_species_count} = 0;
    my $active_regions_ref = $intergenics{$i_id}{active_regions};
    my $silent_regions_ref = $intergenics{$i_id}{silent_regions};
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

  foreach my $i_id (sort {$intergenics{$a}->{start_chr} <=> $intergenics{$b}->{start_chr}} keys %intergenics){
    #Skip failing intergenics and intergenic regions!!
    unless (($intergenics{$i_id}{pass} == 1)){
      next();
    }

    #test both strands 
    my $chr_strand1 = "chr$chr_filter"."_"."+";
    my $chr_strand2 = "chr$chr_filter"."_"."-";
    my @chr_strands;
    push(@chr_strands, $chr_strand1);
    push(@chr_strands, $chr_strand2);

    #Go through each seq type (mrna, est, xmrna, and xest) and test for overlaps with that type
    my $i_start_chr = $intergenics{$i_id}{start_chr};
    my $i_end_chr = $intergenics{$i_id}{end_chr};
  
    &testSeqOverlap('-id'=>$i_id, '-regions_ref'=>\%intergenics, '-strands'=>\@chr_strands, '-seqs_ref'=>\%seqs, '-start_chr'=>$i_start_chr, '-end_chr'=>$i_end_chr);

    my $active_regions_ref = $intergenics{$i_id}{active_regions};
    my $silent_regions_ref = $intergenics{$i_id}{silent_regions};
    
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
  
  my $intergenic_db_file = "$outdir"."chr$chr_filter"."_"."$region_number"."_intergenics_annotated.txt";
  my $fasta_file = "$outdir"."chr$chr_filter"."_"."$region_number"."_intergenics.fa";
  my $intergenic_regions_active_file = "$outdir"."chr$chr_filter"."_"."$region_number"."_intergenicRegionsActive.txt";
  my $intergenic_regions_silent_file = "$outdir"."chr$chr_filter"."_"."$region_number"."_intergenicRegionsSilent.txt";

  open (INTERGENIC_DB, ">$intergenic_db_file") || die "\nCould not open intergenic db file: $intergenic_db_file\n\n";
  open (FASTA, ">$fasta_file") || die "\nCould not open fasta file: $fasta_file\n\n";
  open (ACTIVE, ">$intergenic_regions_active_file") || die "\nCould not open active intergenic regions file: $intergenic_regions_active_file\n\n";
  open (SILENT, ">$intergenic_regions_silent_file") || die "\nCould not open silent intergenic regions file: $intergenic_regions_silent_file\n\n";

  print INTERGENIC_DB "Intergenic_ID\tSeq_Name\tChromosome\tStrand\tUnit1_start_chr\tUnit1_end_chr\tBase_Count\tUnMasked_Base_Count\tCoding_Base_Count\tActive_Region_Count\tActive_Base_Count\tSilent_Region_Count\tSilent_Base_Count\tUpstream_Gene_ID\tDownstream_Gene_ID\tSupporting_mRNA_Count\tSupporting_EST_Count\tSupporting_xmRNA_Count\tSupporting_xEST_Count\tConserved_Species_Count\n";

  print ACTIVE "Active_Region_ID\tIntergenic_ID\tSeq_Name\tChromosome\tStrand\tUnit1_start_chr\tUnit1_end_chr\tBase_Count\tUnMasked_Base_Count\tCoding_Base_Count\tSupporting_mRNA_Count\tSupporting_EST_Count\tSupporting_xmRNA_Count\tSupporting_xEST_Count\tConserved_Species_Count\n";
  print SILENT "Silent_Region_ID\tIntergenic_ID\tSeq_Name\tChromosome\tStrand\tUnit1_start_chr\tUnit1_end_chr\tBase_Count\tUnMasked_Base_Count\tCoding_Base_Count\tSupporting_mRNA_Count\tSupporting_EST_Count\tSupporting_xmRNA_Count\tSupporting_xEST_Count\tConserved_Species_Count\n";

  my $count = 0;
  foreach my $i_id (sort {$intergenics{$a}->{start_chr} <=> $intergenics{$b}->{start_chr}} keys %intergenics){
    $count++;
    if ($count == 10){
      $count = 0;
      $| = 1; print BLUE, ".", RESET;  $| = 0;
    }

    my $strand = $intergenics{$i_id}{strand};

    #Skip failing intergenic regions!!
    unless ($intergenics{$i_id}{pass} == 1){
      next();
    }

    #Get the masked sequence

    my $start = $intergenics{$i_id}{start_chr};
    my $size = $intergenics{$i_id}{size};

    #3-A.) Print out intergenic sequence records to a DB file
    my $chr_intergenic_id = "$chr_filter"."_"."$region_number"."_"."$i_id";

    print INTERGENIC_DB "$chr_intergenic_id\t$intergenics{$i_id}{name}\t$chr_filter\t$intergenics{$i_id}{strand}\t$intergenics{$i_id}{start_chr}\t$intergenics{$i_id}->{end_chr}\t$intergenics{$i_id}{size}\t$intergenics{$i_id}{unmasked_base_count}\t0\t$intergenics{$i_id}{active_region_count}\t$intergenics{$i_id}{active_base_count}\t$intergenics{$i_id}{silent_region_count}\t$intergenics{$i_id}{silent_base_count}\t$intergenics{$i_id}{upstream_gene}\t$intergenics{$i_id}{downstream_gene}\t$intergenics{$i_id}{supporting_mrna_count}\t$intergenics{$i_id}{supporting_est_count}\t$intergenics{$i_id}{supporting_xmrna_count}\t$intergenics{$i_id}{supporting_xest_count}\t$intergenics{$i_id}{conserved_species_count}\n";

    #3-B.) Print out a fasta file of the intergenic sequences
    #    i.e. ">$intergenic_id\n$intergenic_sequence\n"
    my $corrected_start = ($start - $start_filter)+1;
    print FASTA ">$chr_intergenic_id\n", substr($chr_seq{$chr_filter}, $corrected_start, $size), "\n";


    #3-C.) Print out DB files for 'IntergenicRegion_Active' and 'IntergenicRegion_Silent'
    my $active_regions_ref = $intergenics{$i_id}{active_regions};
    my $silent_regions_ref = $intergenics{$i_id}{silent_regions};

    #Active regions first
    foreach my $ar (sort {$active_regions_ref->{$a}->{chr_start} <=> $active_regions_ref->{$b}->{chr_start}} keys %{$active_regions_ref}){

      #Determine the number of unmasked bases with this active region...
      my $start = $active_regions_ref->{$ar}->{chr_start};
      my $size = $active_regions_ref->{$ar}->{size};
      my $corrected_start = ($start - $start_filter)+1;
      my $masked_n_count = (substr($chr_seq{$chr_filter}, $corrected_start, $size) =~ tr/N/N/);
      my $unmasked_count = $size - $masked_n_count;

      print ACTIVE "$active_regions_ref->{$ar}->{active_region_id}\t$chr_intergenic_id\t$active_regions_ref->{$ar}->{name}\t$chr_filter\t$strand\t$active_regions_ref->{$ar}->{chr_start}\t$active_regions_ref->{$ar}->{chr_end}\t$size\t$unmasked_count\t0\t$active_regions_ref->{$ar}->{supporting_mrna_count}\t$active_regions_ref->{$ar}->{supporting_est_count}\t$active_regions_ref->{$ar}->{supporting_xmrna_count}\t$active_regions_ref->{$ar}->{supporting_xest_count}\t$active_regions_ref->{$ar}->{conserved_species_count}\n";
    }

    #Silent regions next
    foreach my $sr (sort {$silent_regions_ref->{$a}->{chr_start} <=> $silent_regions_ref->{$b}->{chr_start}} keys %{$silent_regions_ref}){

      #Determine the number of unmasked bases with this active region...
      my $start = $silent_regions_ref->{$sr}->{chr_start};
      my $corrected_start = ($start - $start_filter)+1;
      my $size = $silent_regions_ref->{$sr}->{size};
      my $masked_n_count = (substr($chr_seq{$chr_filter}, $corrected_start, $size) =~ tr/N/N/);
      my $unmasked_count = $size - $masked_n_count;

      print SILENT "$silent_regions_ref->{$sr}->{silent_region_id}\t$chr_intergenic_id\t$silent_regions_ref->{$sr}->{name}\t$chr_filter\t$strand\t$silent_regions_ref->{$sr}->{chr_start}\t$silent_regions_ref->{$sr}->{chr_end}\t$size\t$unmasked_count\t0\t$silent_regions_ref->{$sr}->{supporting_mrna_count}\t$silent_regions_ref->{$sr}->{supporting_est_count}\t$silent_regions_ref->{$sr}->{supporting_xmrna_count}\t$silent_regions_ref->{$sr}->{supporting_xest_count}\t$silent_regions_ref->{$sr}->{conserved_species_count}\n";
    }
  } 

  close(INTERGENIC_DB);
  close(FASTA);
  close(ACTIVE);
  close(SILENT);
  
  return();
}


##################################################################################################
#Display compute time and memory usage                                                           #
##################################################################################################
sub processUpdate{
  my %args = @_;
  my $display_option = $args{'-display_option'};

  unless($display_option){
    $display_option = 0;
  }

  #Get memory usage in k and percent 
  my $ps_cmd = "ps -p $pid -o \"rss pmem\"";
  my $result = `$ps_cmd`;
  #print Dumper $result;
  my $memory_k = "Unknown";
  my $memory_m = "Unknown";
  my $memory_g = "Unknown";
  my $pmem = "Unknown";
  if ($result =~ /(\d+)\s+(\d+\.\d+)/m){
    $memory_k = $1;
    $memory_m = $memory_k/1024;
    $memory_g = sprintf("%.2f", ($memory_m/1024));
    $pmem = $2;
  } 

  #Get time since last check
  my $t2 = new Benchmark;
  my $td = timediff($t2, $t1);
  my $time_string = timestr($td);
  my $seconds = "Unknown";
  if ($time_string =~ /(\d+)\s+wallclock/){
    $seconds = $1;
  }

  $| = 1;
  if ($display_option == 1){
    print CYAN, "\n\tPROCESS: Time since last update: $seconds seconds\tCurrent memory usage: $memory_g Gb ($pmem%)\n", RESET;
  }elsif($display_option == 2){
    print CYAN, "\tCurrent memory usage: $memory_g Gb ($pmem%)", RESET;
  }else{
    print CYAN, "\nPROCESS: Time since last update: $seconds seconds\tCurrent memory usage: $memory_g Gb ($pmem%)\n", RESET;
  }
  $| = 0;
  $t1 = $t2;
  return();
}

