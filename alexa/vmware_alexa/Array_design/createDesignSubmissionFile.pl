#!/usr/bin/perl -w
#Written by Malachi Griffith

#The purpose of this script is to generate a design submission file for a custom microarray.
#This file is essentially a list of probe sequences with unique IDs.
#A custom array manufacturer (NimbleGen, Agilent, etc.) would use this to create a physical microarray representing the design file

#To use this script you must first have created or downloaded a complete ALEXA database for your species/genome of interest
#The script works by accepting a list of EnsEMBL gene IDs and getting the probes corresponding to each of these from the ALEXA database
#An design file is then generated with the following information: (probe id, target gene id, probe sequence, source database)
#An annotation file is also generated to describe each of the genes targeted and a log file records the creation steps.

#STEPS
#1.) Import file containing list of genes to be targeted (specified by EnsEMBL Gene ID)
#    - Check all gene IDs to make sure they exist in the specified database
#    - 'Invalid Gene ID. - (not in database at all - wrong species or EnsEMBL version?)

#2.) Check filtered sets available
#    - If more than one filtered set has been defined in the database, ask the user which to use
#    - Define the filtered probe set ID to be used for following steps

#3.) Test output directory specified by the user
#    - This directory will be used to store the resulting submission file, annotation file, log file, success summary, etc.
#    - Create a sub-directory within this directory to place a probe file for creation of custom UCSC tracks
#    - Open all output files

#4.)  Go through each valid gene (process them in batches of 100) and do the following:

#4a.)  Get gene info (number exons, number transcripts, etc. that will be needed for all steps)
#      - Also get the filtered probes for this gene (or unfiltered is specified by the user)

#4b.)  Check if the gene has too many exons (more than the max specified by the user) - skip gene if neccessary
#      - Set the 'class' of such genes as 'Too Many Exons' - Valid gene ID but with more exons than allowed by the user
#      - Set the 'status' of such genes as 'Skipped'

#4c.)  Check for genes with 0 probes - A valid gene ID but 0 successful probes (could be pseudogene, predicted gene, etc.) - skip gene if neccessary
#      - Set the 'class' of such genes as 'No Probes'
#      - Set the 'status' of such genes as 'Skipped'

#4d.)  Determine which probes will actually be placed on the array.
#      - For exon-junction probesets, skip those that exceed the exon-skip limit specified by the user
#      - Select probes for each probeset according to the specified probeset size
#        - Select 'best' exon/intron region probes - Based on Tm, length and minimizing overlap
#        - Select 'best' exon junction/boundary probes - Based on Tm


#4e.)  - Count the resulting set of passing selected probes for each gene
#      - Determine the gene probe coverage based on the probes selected (and see if it meets the user specifed cutoff) - skip gene if neccessary
#      - Check for genes that dont meet the probe coverage cutoff - Valid gene ID but with less than the desired 'percent probe coverage'
#      - Set the 'class' of such genes as 'Insufficient Probe Coverage'
#      - Set the 'status' of such genes as 'Skipped'

#4f.) If user specified that only EST/mRNA supported alternative junction probes are allowed, apply this filter now
#     - Identify probes with a near perfect match (95% of length) to an EST, mRNA or EnsEMBL transcript from the target locus
#     - All other probes will be set to a status of 'Fail'
#     - Note that this must be done AFTER checking for gene probe coverage (otherwise it will affect that calculation)

#4g.)  Set the final status and class of genes making it this far
#      - Set the 'class' of such genes as 'Approved' - no problems

#4h.)  Go through each gene in the block and starting adding probes to the array
#      - Check the number of probes already placed on the array and compare to the user specified array capacity
#      - If these probes can be added while still allowing enough room to add the user specified number of negative controls, then add them
#      - Now set the 'status' of such genes as 'Targeted' and continue
#      - If there is not enough space to accomodate this gene set the 'status' as 'No space' and skip this gene.

#4i.)  Write to output files
#      - Write probes for this gene to gene submission file (probe id, target gene id, probe sequence, source database)
#      - Write a record for this gene to the annotation file
#      - Write out probe records to a file suitable for generating custom UCSC tracks:
#        - (Probe_Count,Gene_ID,Probe_Type,Unit1_start_chr,Unit1_end_chr,Unit2_start_chr,Unit2_end_chr
#      - Write an entry for this gene in a design success/failure file:
#        - Contains a line for each Gene ID in the input file (same order)
#        - For each gene ID provide the following info:
#          - Input ID, ALEXA gene id, gene_type, gene_evidence, gene_name, exon_count, total_probe_count, percent_probe_coverage
#          - Specify the 'class' defined above: 
#          - Also specify a status: 'Targeted', 'Skipped', 'No space'

#5.)  Once all genes have been processed, fill remaining space on the array with Negative Control (NC) probes.
#     - Use filtered NC probes and select them to uniformly cover the range of Tm and length of the experimental probes


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
use utilities::ALEXA_DB qw(:all);

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $target_tm = '';
my $target_length = '';
my $max_length_variance = '';
my $probeset_size = '';
my $all_genes = '';
my $gene_list = '';
my $filtered_probes = '';
my $exon_skip_limit = '';
my $known_events_only = '';
my $array_capacity = '';
my $min_gene_probe_coverage = '';
my $max_exons = '';
my $min_nc_probes = '';
my $out_dir = '';
my $verbose = '';

GetOptions ('database=s'=>\$database, 'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'target_tm=f'=>\$target_tm, 'target_length=i'=>\$target_length, 'max_length_variance=i'=>\$max_length_variance, 'probeset_size=i'=>\$probeset_size,
	    'all_genes=s'=>\$all_genes, 'gene_list=s'=>\$gene_list,
	    'filtered_probes=s'=>\$filtered_probes, 'exon_skip_limit=i'=>\$exon_skip_limit,
	    'known_events_only=s'=>\$known_events_only, 'array_capacity=i'=>\$array_capacity,
	    'min_gene_probe_coverage=f'=>\$min_gene_probe_coverage, 'max_exons=i'=>\$max_exons,
	    'min_nc_probes=i'=>\$min_nc_probes, 'out_dir=s'=>\$out_dir, 'verbose=s'=>\$verbose);

#Provide instruction to the user
print GREEN, "\n\nThis script generates an array design submission file for all genes in an ALEXA database or genes specified in an input file\n\n", RESET;
print GREEN, "\nNOTE: NimbleGen's internal control probes will not be imported", RESET;

print GREEN, "\nParameters:", RESET;
print GREEN, "\n\tSpecify the target database and server using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password using: --user and --password", RESET;
print GREEN, "\n\tSpecify the target tm, length and length variance used for this database design using: ",RESET;
print GREEN, "\n\t\t--target_tm (e.g. 67.0), --target_length (e.g. 36) and --max_length_variance (e.g. 8 or whatever was used to generate exon probes)", RESET;
print GREEN, "\n\t\tMake sure these are the correct values for this design (for pre-computed designs, refer to the ALEXA website for this info)", RESET;
print GREEN, "\n\tSpecify the desired number of probes per probeset (exon, intron or junction) using: --probeset_size (e.g. 3)", RESET;
print GREEN, "\n\tSpecify what genes are to be targeted:", RESET;
print GREEN, "\n\t\tTo target all genes use: --all_genes=yes", RESET;
print GREEN, "\n\t\tOtherwise provide an input file with one EnsEMBL gene ID per line (no header) using: --gene_list", RESET;
print GREEN, "\n\tTo allow all probes regardless of quality use: --filtered_probes=no (not recommended)", RESET;
print GREEN, "\n\tSpecify the maximum number of exons skipped to be allowed for exon-junction probes using: --exon_skip_limit (0 for no limit)", RESET;
print GREEN, "\n\tTo limit the array to only probes supported by a 95% matching EST, mRNA or EnsEMBL transcript use: --known_events_only=yes", RESET;
print GREEN, "\n\tSpecify the number of spots available on your custom array platform using: --array_capacity (e.g. 385000, 0 for no limit)", RESET;
print GREEN, "\n\tSpecify the minumum percent gene probe coverage (actual probes/ideal number) using: --min_gene_probe_coverage (0 to ignore)", RESET;
print GREEN, "\n\tSpecify the maximum number of exons for a single gene using: --max_exons (0 to ignore)", RESET;
print GREEN, "\n\tSpecify the minimum number of negative control probes you want on the final array using: --min_nc_probes (0 for no minimum)", RESET;
print GREEN, "\n\t\tNegative control probes will be used to fill empty space on your array after targeting as many genes as possible", RESET;
print GREEN, "\n\tSpecify the output directory for all result files using: --out_dir", RESET;

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tcreateDesignSubmissionFile.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --target_tm=67.0  --target_length=36  --max_length_variance=8  --probeset_size=3  --gene_list=target_genes.txt  --filtered_probes=yes  --exon_skip_limit=0  --known_events_only=yes  --array_capacity=385000  --min_gene_probe_coverage=75.0  --max_exons=40  --min_nc_probes=4000  --out_dir=test1_hs_35_35h\n\n", RESET;

#Make sure all required parameters were specified and do basic checks of the format of each parameter
unless ($database && $server && $user && $password && $target_tm && ($target_length =~ /\d+/) && ($max_length_variance =~ /\d+/) && ($probeset_size =~ /\d+/) && ($filtered_probes =~ /^yes$|^no$/i) && ($exon_skip_limit =~ /\d+/) && ($known_events_only =~ /^yes$|^no$/i) && ($array_capacity =~ /\d+/) && ($min_gene_probe_coverage =~ /\d+/) && ($max_exons =~ /\d+/) && ($min_nc_probes =~ /\d+/) && $out_dir){
  print RED, "\nRequired input parameter(s) missing or incorrect format!\n\n", RESET;
  exit();
}
unless (($all_genes =~ /^yes$/i) || $gene_list){
  print RED, "\nYou must specify an input gene list with --gene_list or use all genes with --all_genes=yes\n\n", RESET;
  exit();
}

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);


#1.) Import file containing list of genes to be targeted (specified by EnsEMBL Gene ID)
#    - Check all gene IDs to make sure they exist in the specified database
my @target_genes;
if ($all_genes =~ /^yes$/i){
  #Note that while I specify All genes here, they may not have all had probes selected (most likely pseudogenes have been excluded at the least)
  @target_genes = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};
}else{
  #Open the input file, test that it is valid, get the EnsEMBL gene IDs present and query the ALEXA database to see if they are valid
  @target_genes = @{&importGeneList('-dbh'=>$alexa_dbh, '-gene_list'=>$gene_list)};
}


#2.) Check filtered sets available
#    - If more than one filtered set has been defined in the database, ask the user which to use
#    - Define the filtered probe set ID to be used for following steps
my $filtered_set;
unless ($filtered_probes =~ /^no$/){
  $filtered_set = &selectFilteredSet ('-dbh'=>$alexa_dbh);
}


#3.) Test output directory specified by the user
#    - This directory will be used to store the resulting submission file, annotation file, log file, success summary, etc.
#    - Create a sub-directory within this directory to place a probe file for creation of custom UCSC tracks
#    - Open all output files
unless ($out_dir =~ /.*\/$/){
  $out_dir = "$out_dir"."/";
}
my $working_dir = &createNewDir('-path'=>$out_dir, '-new_dir_name'=>"probes");
my $submission_file = "$out_dir"."Design_Probe_Submission.txt";
my $annotation_file = "$out_dir"."Design_Gene_Annotation.txt";
my $design_success_file = "$out_dir"."Design_Success_Log.txt";
my $probe_record_file = "$working_dir"."ProbeRecordFile.txt";

#Open output files
open (SUBMISSION, ">$submission_file") || die "\nCould not open submission file for writing: $submission_file\n\n";
open (ANNOTATION, ">$annotation_file") || die "\nCould not open annotation file for writing: $annotation_file\n\n";
open (DESIGN_SUCCESS, ">$design_success_file") || die "\nCould not open design success log file for writing: $design_success_file\n\n";
open (PROBES, ">$probe_record_file") || die "\nCould not open probe record file for writing: $probe_record_file\n\n";

#Print headers
print DESIGN_SUCCESS "EnsEMBL_Gene_ID\tALEXA_Gene_ID\tGene_Name\tGene_Type\tEvidence\tTranscript_Count\tExon_Count\tFinal_Probe_Count\tGene_Probe_Coverage\tDesign_Class\tDesign_Status\n";

print ANNOTATION "EnsEMBL_Gene_ID\tALEXA_Gene_ID\tGene_Source\tGene_Name\tGene_Type\tEvidence\tTranscript_Count\tExon_Count\tFinal_Probe_Count\tGene_Probe_Coverage\tChromosome\tStrand\tChr_Start\tChr_End\n";

print SUBMISSION "ALEXA_Probe_ID\tALEXA_ProbeSet_ID\tEnsEMBL_Gene_ID\tALEXA_Gene_ID\tGene_Source\tSequence\tProbe_Type\n";

print PROBES "Probe_Count\tProbeSet_ID\tGene_ID\tProbe_Type\tSequence\tProbe_Length\tProbe_Tm\tChromosome\tStrand\tUnit1_start_chr\tUnit1_end_chr\tUnit2_start_chr\tUnit2_end_chr\n";

#4.)  Go through each valid gene (process blocks of genes at a time) and do the following:
my $gene_count = @target_genes;

print BLUE, "\nBegin processing $gene_count genes is batches of 100\n\n", RESET;

my $genes_processed = 0;
my $total_junction_probes = 0;
my $total_known_junction_probes = 0;
my $total_exon_junction_probes = 0;
my $total_exons_skipped_limit_probes = 0;
my $total_junction_probes2 = 0;
my $total_junction_probes_size_limit = 0;
my $exon_skip_exceptions = 0;

my $space_limit = $array_capacity - $min_nc_probes;
my $current_array_usage = 0;

my $too_many_exons_count = 0;
my $insufficient_coverage_count = 0;
my $targeted_count = 0;
my $no_space_count = 0;

my @sub_gene_list;
my %gene_data; my $gene_data_ref = \%gene_data;

foreach my $gene_id (@target_genes){
  $genes_processed++;
  push (@sub_gene_list, $gene_id);

  #If the sub gene list contains 100 genes or the last gene has been reached submit the query and reset variables
  my $current_gene_count = @sub_gene_list;

  if ($current_gene_count == 100 || $genes_processed == $gene_count){

    #PROCESS A SET OF 100 GENES

    #Get gene info for these genes to be used later

    #4a.)  Get gene info (number exons, number transcripts, etc. that will be needed for all steps)
    #      - Also get the filtered probes for this gene (or unfiltered is specified by the user)
    #      - Also calculate the number of possible exon regions and junctions for each gene
    $gene_data_ref = &getGeneData('-dbh'=>$alexa_dbh, '-gene_ids'=>\@sub_gene_list, '-exon_skip_limit'=>$exon_skip_limit, '-max_length_variance'=>$max_length_variance);


    #4b.)  Check if each gene has too many exons (more than the max specified by the user) - skip gene if neccessary
    #      - Set the 'class' of such genes as 'Too Many Exons' - Valid gene ID but with more exons than allowed by the user
    #      - Set the 'status' of such genes as 'Skipped'
    foreach my $gene_id (@sub_gene_list){
      #print YELLOW, "\nDEBUG $sub_count: gene: $gene_id max_exons: $max_exons vs exon_count: $gene_data_ref->{$gene_id}->{exon_count}\n", RESET;
      if ($max_exons > 0 && $gene_data_ref->{$gene_id}->{exon_count} > $max_exons){
	$gene_data_ref->{$gene_id}->{class} = "Too Many Exons";
	$gene_data_ref->{$gene_id}->{status} = "Skipped";
      }
    }

    #4c.)  Check for genes with 0 probes - A valid gene ID but 0 successful probes (could be pseudogene, predicted gene, etc.) - skip gene if neccessary
    #      - Set the 'class' of such genes as 'No Probes'
    #      - Set the 'status' of such genes as 'Skipped'
    foreach my $gene_id (@sub_gene_list){

      #If this gene has already been skipped, dont bother processing if for this step
      if ($gene_data_ref->{$gene_id}->{status} eq "Skipped"){
	next();
      }

      if ($gene_data_ref->{$gene_id}->{probe_count} == 0){
	$gene_data_ref->{$gene_id}->{class} = "No Probes";
	$gene_data_ref->{$gene_id}->{status} = "Skipped";
      }
    }

    #4d.)  Determine which probes will actually be placed on the array.
    #      - For exon-junction probesets, skip those that exceed the exon-skip limit specified by the user
    #      - Select probes for each probeset according to the specified probeset size
    #        - Select 'best' exon/intron region probes - Based on Tm, length and minimizing overlap
    #        - Select 'best' exon junction/boundary probes - Based on Tm


    #Skip exon-exon probes that correspond to more than the exon_skip_limit
    #The status of these probes will be changed to 'Fail'
    if ($exon_skip_limit > 0){
      &testExonSkipLimit('-gene_object'=>$gene_data_ref);
    }

    #Select 'best' exon/intron region probes - Based on Tm, length and minimizing overlap
    &selectBestRegionProbes('-gene_object'=>$gene_data_ref, '-probeset_size'=>$probeset_size, '-target_tm'=>$target_tm, '-target_length'=>$target_length);

    #Select 'best' exon junction/boundary probes - Based on Tm
    &selectBestJunctionProbes('-gene_object'=>$gene_data_ref, '-probeset_size'=>$probeset_size, '-target_tm'=>$target_tm);

    #4e.) Determine the gene probe coverage based on the probes selected (and see if it meets the user specifed cutoff) - skip gene if neccessary
    #      - Check for genes that dont meet the probe coverage cutoff - Valid gene ID but with less than the desired 'percent probe coverage'
    #      - Set the 'class' of such genes as 'Insufficient Probe Coverage'
    #      - Set the 'status' of such genes as 'Skipped'
    #      - Note that the exon_skip_limit has already been taken into account in determining the theoretical number of junctions and status of probes
      &determineProbeCoverage('-dbh'=>$alexa_dbh, '-gene_object'=>$gene_data_ref,
			      '-probeset_size'=>$probeset_size, '-min_gene_probe_coverage'=>$min_gene_probe_coverage);


    #4f.) If user specified that only EST/mRNA supported alternative junction probes are allowed, apply this filter now
    #Identify probes with a near perfect match (95% of length) to an EST, mRNA or EnsEMBL transcript from the target locus
    #All other probes will be set to a status of 'Fail'
    #Note that this must be done AFTER checking for gene probe coverage (otherwise it will affect that calculation)
    #Only apply to exon-exon, exon-intron and intron-exon probes (not intron!)
    &identifySequenceSupportedProbes('-gene_object'=>$gene_data_ref);

    #Count the resulting set of passing selected probes for each gene
    foreach my $gene_id (@sub_gene_list){
      #If this gene has already been skipped, dont bother processing if for this step
      if ($gene_data_ref->{$gene_id}->{status} eq "Skipped"){
	next();
      }
      $gene_data_ref->{$gene_id}->{final_probe_count} = 0;
      my $probesets_ref = $gene_data_ref->{$gene_id}->{probesets};
      foreach my $probeset_id (keys %{$probesets_ref}){
	my $probes_ref = $probesets_ref->{$probeset_id}->{probes};

	my $probes_selected = 0;

	foreach my $probe_id (sort {$probes_ref->{$a}->{tm_diff} <=> $probes_ref->{$b}->{tm_diff}} keys %{$probes_ref}){

	  #Go through the junction probes again and for those that exceeded the exon_skip_limit but had sequence support, allow them back in
	  #But still observe the probeset size limit
	  $probes_selected++;
	  if ($probes_ref->{$probe_id}->{probe_type} eq "Exon-Exon"){
	    if ($probes_ref->{$probe_id}->{sequence_support} =~/^yes$/i){
	      if ($exon_skip_limit > 0){
		if ($probes_ref->{$probe_id}->{exons_skipped} > $exon_skip_limit){
		  unless ($probes_selected > $probeset_size){
		    #Change it back to 'Pass' status if neccessary
		    $exon_skip_exceptions++;
		    $probes_ref->{$probe_id}->{status} = "Pass";
		  }
		}
	      }
	    }
	  }

	  if ($probes_ref->{$probe_id}->{status} eq "Pass"){
	    $gene_data_ref->{$gene_id}->{final_probe_count}++;
	  }
	}
      }
    }

    #4g.)  Set the final status and class of genes making it this far
    #      - Set the 'class' of such genes as 'Approved' - no problems
    foreach my $gene_id (@sub_gene_list){
      if ($gene_data_ref->{$gene_id}->{class} eq "Unknown"){
	$gene_data_ref->{$gene_id}->{class} = "Approved";
      }
      if ($verbose eq "yes"){
	print YELLOW, "\nGene: $gene_id\t$gene_data_ref->{$gene_id}->{status}\t$gene_data_ref->{$gene_id}->{class}\tExons: $gene_data_ref->{$gene_id}->{exon_count}\tCoverage: $gene_data_ref->{$gene_id}->{gene_probe_coverage}", RESET;
      }
    }

    #4h.)  Go through each gene in the block and starting adding probes to the array
    #      - Check the number of probes already placed on the array and compare to the user specified array capacity
    #      - If these probes can be added while still allowing enough room to add the user specified number of negative controls, then add them
    #      - Now set the 'status' of such genes as 'Targeted' and continue
    #      - If there is not enough space to accomodate this gene set the 'status' as 'No Space' and skip this gene.

    #4i.)  As the array is being filled, write out results to the output files
    #      - Write probes for this gene to gene submission file (probe_id, probeset_id, Ensembl_Gene_ID, ALEXA_Gene_ID, source database, probe sequence, probe_type)
    #      - Write a record for this gene to the annotation file
    #      - Write out probe records to a file suitable for generating custom UCSC tracks:
    #        - (Probe_Count,ProbeSet_ID,Gene_ID,Probe_Type,Sequence, Tm, Chr,Strand,Unit1_start_chr,Unit1_end_chr,Unit2_start_chr,Unit2_end_chr)
    #      - Write an entry for this gene in a design success/failure file:
    #        - Contains a line for each Gene ID in the input file (same order)
    #        - For each gene ID provide the following info:
    #          - Input ID, ALEXA gene id, gene_name, gene_type, gene_evidence, transcript_count, exon_count, total_probe_count, percent_probe_coverage
    #          - Specify the 'class' defined above: 
    #          - Also specify a status: 'Targeted', 'Skipped', 'No Space'

    foreach my $gene_id (@sub_gene_list){
      #If this gene has already been skipped, dont bother processing if for this step
      if ($gene_data_ref->{$gene_id}->{status} eq "Skipped"){

	#Write an entry for this skipped gene in the design log
	print DESIGN_SUCCESS "$gene_data_ref->{$gene_id}->{ensembl_g_id}\t$gene_id\t$gene_data_ref->{$gene_id}->{gene_name}\t$gene_data_ref->{$gene_id}->{gene_type}\t$gene_data_ref->{$gene_id}->{evidence}\t$gene_data_ref->{$gene_id}->{trans_count}\t$gene_data_ref->{$gene_id}->{exon_count}\tna\t$gene_data_ref->{$gene_id}->{gene_probe_coverage}\t$gene_data_ref->{$gene_id}->{class}\t$gene_data_ref->{$gene_id}->{status}\n";

	next();
      }
      my $final_probe_count = $gene_data_ref->{$gene_id}->{final_probe_count};

      #See if the probes of this gene will fit in the remaining space on the array
      if (($current_array_usage + $final_probe_count) < $space_limit || $array_capacity eq "0"){

	#Mark the gene as targted and write an entry in the design log
	$gene_data_ref->{$gene_id}->{status} = "Targeted";
	print DESIGN_SUCCESS "$gene_data_ref->{$gene_id}->{ensembl_g_id}\t$gene_id\t$gene_data_ref->{$gene_id}->{gene_name}\t$gene_data_ref->{$gene_id}->{gene_type}\t$gene_data_ref->{$gene_id}->{evidence}\t$gene_data_ref->{$gene_id}->{trans_count}\t$gene_data_ref->{$gene_id}->{exon_count}\t$gene_data_ref->{$gene_id}->{final_probe_count}\t$gene_data_ref->{$gene_id}->{gene_probe_coverage}\t$gene_data_ref->{$gene_id}->{class}\t$gene_data_ref->{$gene_id}->{status}\n";


	#Print annotation entry for this gene
	print ANNOTATION "$gene_data_ref->{$gene_id}->{ensembl_g_id}\t$gene_id\t$gene_data_ref->{$gene_id}->{source}\t$gene_data_ref->{$gene_id}->{gene_name}\t$gene_data_ref->{$gene_id}->{gene_type}\t$gene_data_ref->{$gene_id}->{evidence}\t$gene_data_ref->{$gene_id}->{trans_count}\t$gene_data_ref->{$gene_id}->{exon_count}\t$gene_data_ref->{$gene_id}->{final_probe_count}\t$gene_data_ref->{$gene_id}->{gene_probe_coverage}\tchr$gene_data_ref->{$gene_id}->{chromosome}\t$gene_data_ref->{$gene_id}->{chr_strand}\t$gene_data_ref->{$gene_id}->{chr_start}\t$gene_data_ref->{$gene_id}->{chr_end}\n";


	#Add probes for this gene to the array, print out an annotation for the record 
	my $probesets_ref = $gene_data_ref->{$gene_id}->{probesets};
	foreach my $probeset_id (sort {$probesets_ref->{$a} <=> $probesets_ref->{$b}} keys %{$probesets_ref}){
	  my $probes_ref = $probesets_ref->{$probeset_id}->{probes};
	  foreach my $probe_id (sort {$probes_ref->{$a} <=> $probes_ref->{$b}} keys %{$probes_ref}){
	    if ($probes_ref->{$probe_id}->{status} eq "Pass"){

	      #Print probe to submission file
	      print SUBMISSION "ALEXA_P_$probe_id\tALEXA_PS_$probeset_id\t$gene_data_ref->{$gene_id}->{ensembl_g_id}\t$gene_id\t$gene_data_ref->{$gene_id}->{source}\t$probes_ref->{$probe_id}->{sequence}\t$probes_ref->{$probe_id}->{probe_type}\n";

	      #Print probe record to probe record file
	      #(Probe_Count,ProbeSet_ID,Gene_ID,Probe_Type,Sequence, Tm, Chr,Strand,Unit1_start_chr,Unit1_end_chr,Unit2_start_chr,Unit2_end_chr)
	      print PROBES "$probe_id\t$probeset_id\t$gene_id\t$probes_ref->{$probe_id}->{probe_type}\t$probes_ref->{$probe_id}->{sequence}\t$probes_ref->{$probe_id}->{probe_length}\t$probes_ref->{$probe_id}->{tm_celsius}\tchr$gene_data_ref->{$gene_id}->{chromosome}\t$gene_data_ref->{$gene_id}->{chr_strand}\t$probes_ref->{$probe_id}->{unit1_start_chr}\t$probes_ref->{$probe_id}->{unit1_end_chr}\t$probes_ref->{$probe_id}->{unit2_start_chr}\t$probes_ref->{$probe_id}->{unit2_end_chr}\n";

	    }
	  }
	}

	$current_array_usage += $final_probe_count;
      }else{
	#No space for the probes of this gene!

	#Mark the gene as 'No Space' and write an entry in the design log
	$gene_data_ref->{$gene_id}->{status} = "No Space";
	print DESIGN_SUCCESS "$gene_data_ref->{$gene_id}->{ensembl_g_id}\t$gene_id\t$gene_data_ref->{$gene_id}->{gene_name}\t$gene_data_ref->{$gene_id}->{gene_type}\t$gene_data_ref->{$gene_id}->{evidence}\t$gene_data_ref->{$gene_id}->{trans_count}\t$gene_data_ref->{$gene_id}->{exon_count}\tna\t$gene_data_ref->{$gene_id}->{gene_probe_coverage}\t$gene_data_ref->{$gene_id}->{class}\t$gene_data_ref->{$gene_id}->{status}\n";



      }
    }

    foreach my $gene_id (@sub_gene_list){
      if ($gene_data_ref->{$gene_id}->{class} eq "Too Many Exons"){
	$too_many_exons_count++;
      }
      if ($gene_data_ref->{$gene_id}->{class} eq "Insufficient Probe Coverage"){
	$insufficient_coverage_count++;
      }
      if ($gene_data_ref->{$gene_id}->{status} eq "Targeted"){
	$targeted_count++;
      }
      if ($gene_data_ref->{$gene_id}->{status} eq "No Space"){
	$no_space_count++;
      }
    }

    #Reset variables for batch of genes
    $| = 1;
    print ".";
    $| = 0;

    @sub_gene_list = ();
    %{$gene_data_ref} = ();
  }

}


#5.)  Once all genes have been processed, fill remaining space on the array with Negative Control (NC) probes.
#     - Use filtered NC probes and select them to uniformly cover the range of Tm and length of the experimental probes

#If the user specified unlimited capacity, use the min number of NC probes requested.  Otherwise fill space remaining on the array
my $desired_nc_probes;
if ($array_capacity eq "0"){
  $desired_nc_probes = $min_nc_probes;
}else{
  $desired_nc_probes = $array_capacity - $current_array_usage;
}

if ($desired_nc_probes == 0){
  print YELLOW, "\nWARNING - 0 negative probes will be selected! -- Not recommended\n\n", RESET;
}else{
  &selectNegativeControlProbes('-dbh'=>$alexa_dbh, '-nc_probe_count'=>$desired_nc_probes, '-filtered'=>$filtered_probes, '-filtered_set'=>$filtered_set);
}



#Print out the total number of junction probes that had sequence support
#This includes those from all genes that pass the max_exons test and have at least one probe (i.e. some genes that were later skipped)
#This is purely for interests sake
if ($known_events_only =~ /^yes$/i){
  print BLUE, "\n\n$total_known_junction_probes / $total_junction_probes junction probes had sequence support\n\n", RESET;
}
if ($exon_skip_limit > 0){
  print BLUE, "\n$total_exons_skipped_limit_probes / $total_exon_junction_probes probes exceed the exon skip limit\n\n", RESET;
  print BLUE, "\n$exon_skip_exceptions exon-exon junction probes were allowed despite the 'exon skip limit' because they had sequence support\n\n", RESET;
}

print BLUE, "\n$total_junction_probes_size_limit / $total_junction_probes2 probes exceed the probeset size limit\n\n", RESET;


#Summarize the number of genes added to the array the class of those not added to the array
print BLUE "\n\n$too_many_exons_count genes had too many exons and were skipped", RESET;
print BLUE "\n$insufficient_coverage_count genes had insufficient gene probe coverage and were skipped", RESET;
print BLUE "\n$targeted_count genes were approved and successfully targeted", RESET;
print BLUE "\n$no_space_count genes were approved but could not be targeted because the specified array capacity ran out\n\n", RESET;

close (SUBMISSION);
close (ANNOTATION);
close (DESIGN_SUCCESS);
close (PROBES);

#Close database connection
$alexa_dbh->disconnect();

exit();


########################################################################################################################
#Import the gene list specified by the user, get ALEXA IDs for each EnsEMBL gene ID specified                          #
#Check that each EnsEMBL gene ID is valid and report invalid ones to the user                                          #
#If the input list contains invalid IDs, the user should fix their list before continuing                              #
########################################################################################################################
sub importGeneList{
  my %args = @_;
  my $dbh = $args{'-dbh'};
  my $gene_list = $args{'-gene_list'};

  print BLUE, "\n\nParsing input genes list", RESET;

  my %input_genes;
  my %valid_genes;
  my %invalid_genes;

  #Get the gene ID from each line
  #This script assumes only one ID perline and this ID consists of any non-space characters
  open (GENES, "$gene_list") || die "\nCould not open the specified gene list: $gene_list\n\n";
  my $line_count = 0;
  while (<GENES>){
    if ($_ =~ /(\S+)/){
      $line_count++;
      $input_genes{$1}{order} = $line_count;
    }else{
      print RED, "\nLine format of input file not understood! LINE: $_\n", RESET;
      exit();
    }
  }
  close (GENES);

  my @input_genes = keys %input_genes;
  %valid_genes = %{&getGeneIds ('-dbh'=>$dbh, '-ensembl_g_ids'=>\@input_genes)};

  my $input_gene_count = keys %input_genes;
  my $valid_gene_count = keys %valid_genes;

  print BLUE, "\n\tFound $input_gene_count genes in put gene list and $valid_gene_count of these were found in the ALEXA database", RESET;

  unless ($input_gene_count == $valid_gene_count){

    print YELLOW, "\n\nThe following genes in the input list were not found in the ALEXA database specified", RESET;
    foreach my $input_gene (keys %input_genes){
      unless ($valid_genes{$input_gene}){
	print YELLOW, "\n\t$input_gene", RESET;
      }
    }
    print YELLOW, "\n\nThese IDs are either invalid or belong to a different database version or species than the one you specified", RESET;
    print YELLOW, "\nFix this issue and try again\n\n", RESET;
    exit();
  }

  print "\n\n";

  my @final_genes;

  #Get the alexa gene ID for each of these genes - keep the order the same as the input gene list!
  foreach my $input_gene (sort {$input_genes{$a}->{order} <=> $input_genes{$b}->{order}} keys %input_genes){
    push (@final_genes, $valid_genes{$input_gene}{alexa_gene_id});
  }

  return(\@final_genes);
}

###########################################################################################################################
#4a.)  Get gene info (number exons, number transcripts, etc. that will be needed for all steps)                           #
###########################################################################################################################
sub getGeneData{
  my %args= @_;
  my $dbh = $args{'-dbh'};
  my $gene_ids_ref = $args{'-gene_ids'};
  my $exon_skip_limit = $args{'-exon_skip_limit'};
  my $max_length_variance = $args{'-max_length_variance'};

  #A.) Create a gene data object by first getting basic gene info
  my $gene_info_ref = &getGeneInfo ('-dbh'=>$dbh, '-gene_ids'=>$gene_ids_ref, '-sequence'=>"no", '-silent'=>"yes");

  #B.) Get all exons for these genes
  my $gene_exons_ref = &getExons ('-dbh'=>$dbh, '-gene_ids'=>$gene_ids_ref, '-sequence'=>"no", '-silent'=>"yes");

  #C.) Get additional information such as the number of unique exons and add this to the gene data object
  my $gene_transcripts_ref = &getTranscripts('-dbh'=>$dbh, '-gene_ids'=>$gene_ids_ref, '-sequence'=>"no", '-silent'=>"yes");

  #D.) Get the filtered or unfiltered probes for these genes
  my $gene_probes_ref;
  if ($filtered_probes =~ /^no$/i){
    $gene_probes_ref = &getGeneProbes('-dbh'=>$dbh, '-gene_ids'=>$gene_ids_ref, '-filtered'=>"no", '-silent'=>"yes");
  }else{
    $gene_probes_ref = &getGeneProbes('-dbh'=>$dbh, '-gene_ids'=>$gene_ids_ref, '-filtered'=>"yes", '-filter_set'=>$filtered_set, '-silent'=>"yes");
  }

  #Now add all the additional data to the main gene data info object
  foreach my $gene_id (keys %{$gene_info_ref}){
    my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};

    #E.) Number of transcripts for this gene
    my $trans_count = keys %{$transcripts_ref};
    $gene_info_ref->{$gene_id}->{trans_count} = $trans_count;

    #F.) Number of non-redundant exons for this gene
    my %nr_exons;
    foreach my $trans_id (keys %{$transcripts_ref}){
      my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

      foreach my $exon_id (keys %{$exons_ref}){
	my $start_stop = "$exons_ref->{$exon_id}->{exon_start}"."."."$exons_ref->{$exon_id}->{exon_end}";
	$nr_exons{$start_stop}{tmp} = 'na';
      }
    }
    my $exon_count = keys %nr_exons;
    $gene_info_ref->{$gene_id}->{exon_count} = $exon_count;

    #Add gene probes data to this object
    $gene_info_ref->{$gene_id}->{probe_count} = $gene_probes_ref->{$gene_id}->{probe_count};
    $gene_info_ref->{$gene_id}->{probesets} = $gene_probes_ref->{$gene_id}->{probesets};

    #G.) Initialize the 'class' of each gene
    $gene_info_ref->{$gene_id}->{class} = 'Unknown';  #'Too Many Exons', 'No Probes', 'Insufficient Probe Coverage', 'Approved'
    $gene_info_ref->{$gene_id}->{status} = 'Unknown'; #'Skipped', 'Targeted', 'No Space' (would have been targeted but for the lack of space)

    #H.) Determine the number of exon regions that would be targeted by the design
    #First define clusters of overlapping exons, then define exon regions of these clusters
    my $exons_ref = $gene_exons_ref->{$gene_id}->{exons};
    my $cluster_count = 0;
    my %clusters;

  EXON: foreach my $exon_id (sort {$exons_ref->{$a}->{exon_start} <=> $exons_ref->{$b}->{exon_start}} keys %{$exons_ref}){

      #Get the start and end positions of the current exon
      my $exon_start = $exons_ref->{$exon_id}->{exon_start};
      my $exon_end = $exons_ref->{$exon_id}->{exon_end};

      #Go through each cluster and see if the current exon overlaps with one of the exons already present
      foreach my $cluster_id (sort keys %clusters){
	my @tmp_exons = @{$clusters{$cluster_id}{exons}};
	foreach my $exon_id_test (@tmp_exons){

	  my $exon_start_test = $exons_ref->{$exon_id_test}->{exon_start};
	  my $exon_end_test = $exons_ref->{$exon_id_test}->{exon_end};

	  #See if the start of the current exon is within the range of the test exon
	  if ($exon_start >= $exon_start_test && $exon_start <= $exon_end_test){
	    push (@tmp_exons, $exon_id);
	    $clusters{$cluster_id}{exons} = \@tmp_exons;
	    next(EXON);
	  }
	  #See if the end of the current exon is within the range of the test exon
	  if ($exon_end >= $exon_start_test && $exon_end <= $exon_end_test){
	    push (@tmp_exons, $exon_id);
	    $clusters{$cluster_id}{exons} = \@tmp_exons;
	    next(EXON);
	  }
	  #See if the current exon completely flanks the test exon - if so it should be added to the cluster
	  if ($exon_start <= $exon_start_test && $exon_end >= $exon_end_test){
	    push (@tmp_exons, $exon_id);
	    $clusters{$cluster_id}{exons} = \@tmp_exons;
	    next(EXON);
	  }
	}
      }
      #If the current exon was not added to any of the current clusters - create a new cluster
      $cluster_count++;
      my @tmp_exons;
      push (@tmp_exons, $exon_id);
      $clusters{$cluster_count}{exons} = \@tmp_exons;
    }

    #Now using the exon clusters identified here, determine the number of exon regions exactly as was done when generating the exon probes
    my $total_exon_region_count = 0;
    foreach my $cluster (sort {$a <=> $b} keys %clusters){
      my $cluster_exons_aref = $clusters{$cluster}{exons};
      my $cluster_size = @{$cluster_exons_aref};

      if ($cluster_size == 1){
	#If cluster consists of only a single exon - there is only one exon region
	$total_exon_region_count++;

      }else{
	#If the cluster consists of multiple exons, determine the number of regions to be targeted
	my $exon_regions_ref = &defineExonRegions('-exon_object'=>$exons_ref, '-exon_ids'=>$cluster_exons_aref, '-max_length_variance'=>$max_length_variance);
	my $exon_region_count = keys %{$exon_regions_ref};
	$total_exon_region_count += $exon_region_count;
      }
    }
    #Add the info gathered above to the gene object
    $gene_info_ref->{$gene_id}->{exon_region_count} = $total_exon_region_count;

  }

  #I.) Get the theoretical number of exon-exon and exon-intron junction probes that would ideally be designed for each gene
  my @multi_exon_genes;
  foreach my $gene_id (keys %{$gene_info_ref}){

    $gene_info_ref->{$gene_id}->{exon_junction_count} = 0;
    $gene_info_ref->{$gene_id}->{exon_boundary_count} = 0;

    unless ($gene_info_ref->{$gene_id}->{exon_count} == 1){
      push (@multi_exon_genes, $gene_id);
    }
  }
  my $probe_counts_ref = &junctionProbeCombinations ('-dbh'=>$dbh, '-gene_ids'=>\@multi_exon_genes, '-exon_skip_limit'=>$exon_skip_limit, '-silent'=>"yes");

  foreach my $gene_id (@multi_exon_genes){

    $gene_info_ref->{$gene_id}->{exon_junction_count} = $probe_counts_ref->{$gene_id}->{exon_exon};
    $gene_info_ref->{$gene_id}->{exon_boundary_count} = $probe_counts_ref->{$gene_id}->{intron_exon};
  }

  return($gene_info_ref);
}


################################################################################################
#For exon clusters that have been derived from overlapping exons, identify informative regions #
#of these exons to use for probe design.                                                       #
#For each of these regions, make note of which of the exons from the cluster of overlapping    #
#exons are involved in each of the defined regions                                             #
################################################################################################
sub defineExonRegions{
  my %args = @_;
  my $exon_object_ref = $args{'-exon_object'};
  my $exon_ids_aref = $args{'-exon_ids'};
  my $max_length_variance = $args{'-max_length_variance'};

  #The exons in @exon_ids comprise a single cluster of overlapping exons
  #First try to identify regions of these exons that are unique to each exon

  #Compile a non-redundant list of start/end positions
  #Note that if two exons are in a cluster and have exactly the same start/end positions (which does happen in Ensembl!) then when
  #converted to a non-redundant list they will resolve to a single region which is good.
  my %junctions;
  my @junctions;
  foreach my $exon_id (@{$exon_ids_aref}){
    #Use a hash to compile a non-redundant list of start/end positions
    my $exon_start = $exon_object_ref->{$exon_id}->{exon_start};
    my $exon_end = $exon_object_ref->{$exon_id}->{exon_end};
    $junctions{$exon_start}{tmp}='na';
    $junctions{$exon_end}{tmp}='na';
 }
  #Create a sorted array of these junctions
  foreach my $junct (sort {$a <=> $b} keys %junctions){
    push (@junctions, $junct);
  }

  #Now consider the regions between the junctions and identify the exons associated with each
  my $number_junctions = @junctions;
  my %regions;
  my $region_count;

  for (my $i = 0; $i < $number_junctions-1; $i++){
    my $region_start = ($junctions[$i])+1;
    my $region_end = ($junctions[$i+1])-1;

    #Confirm that the region selected is actually valid and larger than the required probe size
    if ($region_end <= $region_start){
      next();
    }
    my $region_size = ($region_end - $region_start);
    if ($region_size <= ($target_length + $max_length_variance)){
      next();
    }

    #Check each exon to see which overlap within this region at either end - or flank it completely
    my @exons_within_region;
    foreach my $exon_id (@{$exon_ids_aref}){
      my $exon_start = $exon_object_ref->{$exon_id}->{exon_start};
      my $exon_end = $exon_object_ref->{$exon_id}->{exon_end};

      #Is the start position of the selected region within this exon?
      if ($region_start >= $exon_start && $region_start <= $exon_end){
	push (@exons_within_region, $exon_id);
	next();
      }
      #Is the end position of the selected region within this exon?
      if ($region_end >= $exon_start && $region_end <= $exon_end){
	push (@exons_within_region, $exon_id);
	next();
      }
      #Does the selected region completely flank this exon?
      if ($region_start <= $exon_start && $region_end >= $exon_end){
	push (@exons_within_region, $exon_id);
	next();
      }

    }
    $region_count++;
    $regions{$region_count}{region_start} = $region_start;
    $regions{$region_count}{region_end} = $region_end;
    $regions{$region_count}{exons} = \@exons_within_region;
  }

  #Keep track of the exons that are successfully covered by the regions defined
  #Check each region defined to see if it is completely within an exon.
  #For every exon, I want a region that is completely within the boundaries of the exon.  For those exons where this is not true,
  #define a region within the exon and note which other exons it covers - for this evaluation, if a probe covers any amount of sequence
  #of another exon it will be noted.  This is a conservative approach, if there is any ambiguity at all regarding which exons a probe covers,
  #it must be noted.  Nevertheless I want at least one probe that is completely within each exon.
  my %exons_covered;
  foreach my $region (sort {$a <=> $b} keys %regions){

    my $region_start = $regions{$region}{region_start};
    my $region_end = $regions{$region}{region_end};

    #Check this region against the initial list of exons for this exon cluster
    foreach my $exon_id (@{$exon_ids_aref}){
      my $exon_start = $exon_object_ref->{$exon_id}->{exon_start};
      my $exon_end = $exon_object_ref->{$exon_id}->{exon_end};

      #Which exons completely cover this region?
      if ($region_start >= $exon_start && $region_end <= $exon_end){
	$exons_covered{$exon_id}{tmp} = 'na';
      }
    }
  }
  #Now go through the original list of exons for this exon cluster and see which have not been successfully covered
  foreach my $exon_id (@{$exon_ids_aref}){
    unless ($exons_covered{$exon_id}{tmp}){

      #Still need a probe region for this exon!
      my $region_start = ($exon_object_ref->{$exon_id}->{exon_start})+1;
      my $region_end = ($exon_object_ref->{$exon_id}->{exon_end})-1;

      #Make sure this region fits the required probe length
      my $region_size = $region_end - $region_start;
      if ($region_size <= ($target_length + $max_length_variance)){
	next();
      }

      #Define the region and note the exons that it covers.
      my @exons_within_region;
      foreach my $exon_id (@{$exon_ids_aref}){
	my $exon_start = $exon_object_ref->{$exon_id}->{exon_start};
	my $exon_end = $exon_object_ref->{$exon_id}->{exon_end};

	#Is the start position of the selected region within this exon?
	if ($region_start >= $exon_start && $region_start <= $exon_end){
	  push (@exons_within_region, $exon_id);
	  next();
	}
	#Is the end position of the selected region within this exon?
	if ($region_end >= $exon_start && $region_end <= $exon_end){
	  push (@exons_within_region, $exon_id);
	  next();
	}
      }
      $region_count++;
      $regions{$region_count}{region_start} = $region_start;
      $regions{$region_count}{region_end} = $region_end;
      $regions{$region_count}{exons} = \@exons_within_region;
    }
  }

  return(\%regions);
}


##########################################################################################################################
#Skip exon-exon probes that correspond to more than the exon_skip_limit                                                  #
#The status of these probes will be changed to 'Fail'                                                                    #
##########################################################################################################################
sub testExonSkipLimit{
  my %args = @_;
  my $gene_object_ref = $args{'-gene_object'};

  foreach my $gene_id (keys %{$gene_object_ref}){

    #If this gene has already been skipped, dont bother processing if for this step
    if ($gene_object_ref->{$gene_id}->{status} eq "Skipped"){
      next();
    }

    my $probesets_ref = $gene_object_ref->{$gene_id}->{probesets};

    foreach my $probeset_id (keys %{$probesets_ref}){
      my $probes_ref = $probesets_ref->{$probeset_id}->{probes};

      foreach my $probe_id (keys %{$probes_ref}){
	$total_exon_junction_probes++;
	my $probe_type = $probes_ref->{$probe_id}->{probe_type};

	if ($probe_type eq "Exon-Exon"){

	  if ($probes_ref->{$probe_id}->{exons_skipped} > $exon_skip_limit){
	    $probes_ref->{$probe_id}->{status} = "Fail";
	    $total_exons_skipped_limit_probes++;
	  }
	}
      }
    }
  }

  return();
}


####################################################################################################################################
#Select best exon or intron region probe from those that passed previous tests - based on the Tm,target-length, and postion    #
####################################################################################################################################
sub selectBestRegionProbes{
  my %args = @_;
  my $gene_object_ref = $args{'-gene_object'};
  my $probeset_size = $args{'-probeset_size'};
  my $target_tm = $args{'-target_tm'};
  my $target_length = $args{'-target_length'};



  foreach my $gene_id (keys %{$gene_object_ref}){

    #If this gene has already been skipped, dont bother processing if for this step
    if ($gene_object_ref->{$gene_id}->{status} eq "Skipped"){
      next();
    }

    my $probesets_ref = $gene_object_ref->{$gene_id}->{probesets};

    #Go through each probeset (if it corresponds to an Exon or Intron)
    foreach my $probeset_id (keys %{$probesets_ref}){

      #A.) Unless this is actually an exon or intron probeset, skip it.
      unless ($probesets_ref->{$probeset_id}->{probe_type} eq "Exon" || $probesets_ref->{$probeset_id}->{probe_type} eq "Intron"){
	next();
      }

      #B.) Now go through each of the eligable (passed) probes for this probeset and calculate a score based on the Tm and length of the probe
      my %candidate_probes;
      my $probes_ref = $probesets_ref->{$probeset_id}->{probes};

      foreach my $probe_id (keys %{$probes_ref}){
	if ($probes_ref->{$probe_id}->{status} eq "Pass"){

	  my $tm = $probes_ref->{$probe_id}->{tm_celsius};
	  my $tm_diff = abs($tm - $target_tm);

	  my $length = $probes_ref->{$probe_id}->{probe_length};
	  my $length_diff = abs($length - $target_length);

	  #Calculate the score as follows so that Tm gives up to a value of 6 times the Tm range of probes and length is taken at par
	  #This means that if exon probes were allowed to vary by +/- 3.5 degrees and +/- 8 bp in length the total score possible is 21 + 8 = 29
	  #A perfect score of 0 would occur for probes with the target length and tm.  Slightly more weight is given to the influence of Tm than length
	  #This should be okay considering that the Tm of exon probes is only allowed to vary by 3.5 degrees anyway.
	  my $quality_score = ($tm_diff * 6) + $length_diff;

	  $candidate_probes{$probe_id}{quality_score} = $quality_score;
	}
      }

      if ($verbose eq "yes"){
	print YELLOW, "\n\nProbeset: $probeset_id", RESET;
	print YELLOW, "\nCandidate probes sorted by quality score:", RESET;

	foreach my $cp_id (sort {$candidate_probes{$a}->{quality_score} <=> $candidate_probes{$b}->{quality_score}} keys %candidate_probes){
	  print YELLOW, "\n$cp_id\tTm: $probes_ref->{$cp_id}->{tm_celsius}\tLength: $probes_ref->{$cp_id}->{probe_length}\tscore: $candidate_probes{$cp_id}{quality_score}\tStart:$probes_ref->{$cp_id}->{unit1_start}\tEnd: $probes_ref->{$cp_id}->{unit1_end}", RESET;
	}
      }

      #C.) Now go the list of candidate probes as many times as the number of desired probes in the final probeset
      #    Select the best 'N' probes to fill this requirement.  Use Tm and probe length to score probes, and attempt to avoid overlaping probes
      my %best_region_probes;
    SELECT:for (my $i = 0; $i < $probeset_size; $i++){
	
	#Go through the candidate probes starting with the one with the best score
	foreach my $candidate_probe_id (sort {$candidate_probes{$a}->{quality_score} <=> $candidate_probes{$b}->{quality_score}} keys %candidate_probes){

	  #If this probe is already part of the list of those chosen, continue
	  if ($best_region_probes{$candidate_probe_id}){
	    next();
	  }

	  #Go through each of the probes already selected and see if this current probe can be added
	  my $overlap = 0;
	  foreach my $best_probe_id (keys %best_region_probes){

	    #Test overlap - 5' end of candidate probe
	    if (($probes_ref->{$candidate_probe_id}->{unit1_start} >= $probes_ref->{$best_probe_id}->{unit1_start}) && ($probes_ref->{$candidate_probe_id}->{unit1_start} <= $probes_ref->{$best_probe_id}->{unit1_end})){
	      $overlap = 1;
	    }
	    #Test overlap - 3' end of candidate probe
	    if (($probes_ref->{$candidate_probe_id}->{unit1_end} >= $probes_ref->{$best_probe_id}->{unit1_start}) && ($probes_ref->{$candidate_probe_id}->{unit1_end} <= $probes_ref->{$best_probe_id}->{unit1_end})){
	      $overlap = 1;
	    }
	    #Test overlap - complete flank
	    if (($probes_ref->{$candidate_probe_id}->{unit1_start} <= $probes_ref->{$best_probe_id}->{unit1_start}) && ($probes_ref->{$candidate_probe_id}->{unit1_end} >= $probes_ref->{$best_probe_id}->{unit1_end})){
	      $overlap = 1;
	    }
	  }

	  #If no overlap was found, add this probe to the best_region_probes list and continue to the next selection
	  unless ($overlap == 1){
	    $best_region_probes{$candidate_probe_id}{quality_score} = $candidate_probes{$candidate_probe_id}{quality_score};
	    next SELECT;
	  }
	}

	#If all candidate probes were considered, but a suitable probe was not found (because of overlap to already selected probes):
	# - Go through the candidates again and simply choose the next best one without worrying about overlap
	foreach my $candidate_probe_id (sort {$candidate_probes{$a}->{quality_score} <=> $candidate_probes{$b}->{quality_score}} keys %candidate_probes){

	  #If this probe is already part of the list of those chosen, continue
	  if ($best_region_probes{$candidate_probe_id}){
	    next();
	  }
	  $best_region_probes{$candidate_probe_id}{quality_score} = $candidate_probes{$candidate_probe_id}{quality_score};
	  next SELECT;
	}
      }

      #DEBUG
      if ($verbose eq "yes"){
	print YELLOW, "\n\nBest region probes actually selected probes sorted by quality score:", RESET;

	foreach my $brp_id (sort {$best_region_probes{$a}->{quality_score} <=> $best_region_probes{$b}->{quality_score}} keys %best_region_probes){
	  print YELLOW, "\n$brp_id\tTm: $probes_ref->{$brp_id}->{tm_celsius}\tLength: $probes_ref->{$brp_id}->{probe_length}\tscore: $best_region_probes{$brp_id}{quality_score}\tStart:$probes_ref->{$brp_id}->{unit1_start}\tEnd: $probes_ref->{$brp_id}->{unit1_end}", RESET;
	}
      }

      #F.) Finally, go through all the probes of this probeset and mark all probes as 'Fail' unless they are one of the chosen 'Best' probes
      foreach my $probe_id (keys %{$probes_ref}){
	unless ($best_region_probes{$probe_id}){
	  $probes_ref->{$probe_id}->{status} = "Fail";
	}
      }
    }
  }

  return ();
}



##########################################################################################################################
#Select 'best' exon junction/boundary probes - Based on Tm                                                               #
##########################################################################################################################
sub selectBestJunctionProbes{
  my %args = @_;
  my $gene_object_ref = $args{'-gene_object'};
  my $probeset_size = $args{'-probeset_size'};
  my $target_tm = $args{'-target_tm'};

  foreach my $gene_id (keys %{$gene_object_ref}){

    #If this gene has already been skipped, dont bother processing if for this step
    if ($gene_object_ref->{$gene_id}->{status} eq "Skipped"){
      next();
    }

    my $probesets_ref = $gene_object_ref->{$gene_id}->{probesets};

    foreach my $probeset_id (keys %{$probesets_ref}){
      my $probes_ref = $probesets_ref->{$probeset_id}->{probes};

      #Get the difference between the probe Tm and target Tm for every probe
      foreach my $probe_id (keys %{$probes_ref}){
	my $tm = $probes_ref->{$probe_id}->{tm_celsius};
	my $tm_diff = abs($tm - $target_tm);
	$probes_ref->{$probe_id}->{tm_diff} = $tm_diff;
      }

      my $probe_type = $probesets_ref->{$probeset_id}->{probe_type};

      #Only process exon junction and boundary probes here
      if ($probe_type eq "Exon-Exon" || $probe_type eq "Exon-Intron" || $probe_type eq "Intron-Exon"){
	next();
      }

      #Now go through the probes again, sort by Tm diff, and once the probeset size is reached mark them as 'Fail'
      my $probes_selected = 0;
      foreach my $probe_id (sort {$probes_ref->{$a}->{tm_diff} <=> $probes_ref->{$b}->{tm_diff}} keys %{$probes_ref}){
	$probes_selected++;
	$total_junction_probes2++;

	if ($probes_selected > $probeset_size){
	  $probes_ref->{$probe_id}->{status} = "Fail";
	  $total_junction_probes_size_limit++;
	}
      }
    }
  }
  return();
}


##########################################################################################################################
#determineProbeCoverage
##########################################################################################################################
sub determineProbeCoverage{
  my %args = @_;
  my $dbh = $args{'-dbh'};
  my $gene_object_ref = $args{'-gene_object'};
  my $probeset_size = $args{'-probeset_size'};
  my $min_gene_probe_coverage = $args{'-min_gene_probe_coverage'};

  foreach my $gene_id (keys %{$gene_object_ref}){

    my $probesets_ref = $gene_object_ref->{$gene_id}->{probesets};

    my $passing_e_probes = 0;
    my $passing_ej_probes = 0;
    my $passing_eb_probes = 0;

    foreach my $probeset_id (keys %{$probesets_ref}){
      my $probes_ref = $probesets_ref->{$probeset_id}->{probes};

      my $probe_type = $probesets_ref->{$probeset_id}->{probe_type};

      #Dont count intron probes towards gene coverage
      if ($probe_type eq "Intron"){
	next();
      }

      foreach my $probe_id (keys %{$probes_ref}){
	if ($probes_ref->{$probe_id}->{status} eq "Fail"){
	  next();
	}

	if ($probe_type eq "Exon"){
	  $passing_e_probes++;
	}elsif ($probe_type eq "Exon-Exon"){
	  $passing_ej_probes++;
	}elsif ($probe_type eq "Exon-Intron" || $probe_type eq "Intron-Exon"){
	  $passing_eb_probes++;
	}
      }
    }

    my $exon_region_count = ($gene_object_ref->{$gene_id}->{exon_region_count})*$probeset_size;
    my $exon_junction_count = ($gene_object_ref->{$gene_id}->{exon_junction_count})*$probeset_size;
    my $exon_boundary_count = ($gene_object_ref->{$gene_id}->{exon_boundary_count})*$probeset_size;

    if ($verbose eq "yes"){
      print YELLOW, "\n\nGene: $gene_id\tProbeset_size: $probeset_size\tExons: $gene_object_ref->{$gene_id}->{exon_count}", RESET;
      print YELLOW, "\n\tPassing Probes\tE: $passing_e_probes\tEJ: $passing_ej_probes\tEB: $passing_eb_probes", RESET;
      print YELLOW, "\n\tTheoretical\tE: $exon_region_count\tEJ: $exon_junction_count\tEB: $exon_boundary_count", RESET;
    }

    my $actual_count = $passing_e_probes + $passing_ej_probes + $passing_eb_probes;
    my $theoretical_count = $exon_region_count + $exon_junction_count + $exon_boundary_count;
    my $gene_coverage = ($actual_count/$theoretical_count)*100;

    $gene_object_ref->{$gene_id}->{gene_probe_coverage} = $gene_coverage;

    if ($gene_coverage < $min_gene_probe_coverage && $gene_object_ref->{$gene_id}->{status} eq "Unknown"){
      $gene_object_ref->{$gene_id}->{class} = "Insufficient Probe Coverage";
      $gene_object_ref->{$gene_id}->{status} = "Skipped";
    }

  }
  return();
}


###########################################################################################################################
#Identify probes with a near perfect match (95% of length) to an EST, mRNA or EnsEMBL transcript from the target locus
#Only apply to exon-exon, exon-intron and intron-exon probes
###########################################################################################################################
sub identifySequenceSupportedProbes{
  my %args = @_;
  my $gene_object_ref = $args{'-gene_object'};

  foreach my $gene_id (keys %{$gene_object_ref}){

    #If this gene has already been skipped, dont bother processing if for this step
    if ($gene_object_ref->{$gene_id}->{status} eq "Skipped"){
      next();
    }

    my $probesets_ref = $gene_object_ref->{$gene_id}->{probesets};

    foreach my $probeset_id (keys %{$probesets_ref}){
      my $probes_ref = $probesets_ref->{$probeset_id}->{probes};

      foreach my $probe_id (keys %{$probes_ref}){

	#For each probe see if there is a supporting sequence of the required length or longer
	my $probe_type = $probes_ref->{$probe_id}->{probe_type};
	
	if ($probe_type eq "Exon-Exon" || $probe_type eq "Exon-Intron" || $probe_type eq "Intron-Exon"){
	  my $probe_length = $probes_ref->{$probe_id}->{probe_length};
	  my $est_thl = $probes_ref->{$probe_id}->{est_thl};
	  my $mrna_thl = $probes_ref->{$probe_id}->{mrna_thl};
	  my $enst_thl = $probes_ref->{$probe_id}->{enst_thl};
	  my $sequence_support = 0;

	  if ($est_thl =~ /\d+/){
	    my $hit_percent = ($est_thl/$probe_length)*100;
	    if ($hit_percent >= 95){$sequence_support = 1;}
	  }
	  if ($mrna_thl =~ /\d+/){
	    my $hit_percent = ($mrna_thl/$probe_length)*100;
	    if ($hit_percent >= 95){$sequence_support = 1;}
	  }
	  if ($enst_thl =~ /\d+/){
	    my $hit_percent = ($enst_thl/$probe_length)*100;
	    if ($hit_percent >= 95){$sequence_support = 1;}
	  }

	  #Now fail the probe if it didnt have a supporting sequence (all probe start with a passing status)
	  if ($sequence_support == 1){
	    if ($known_events_only =~ /^yes$/i){
	      $total_junction_probes++;
	      $total_known_junction_probes++;
	    }
	    $probes_ref->{$probe_id}->{sequence_support} = "Yes";
	  }else{
	    if ($known_events_only =~ /^yes$/i){
	      $total_junction_probes++;
	      $probes_ref->{$probe_id}->{status} = "Fail";
	    }
	  }
	}
      }
    }
  }

  return();
}



#############################################################################################################################
#Select negative control probes foreach bin so that they uniformly cover all Tm and Probe-Lengths                       #
#############################################################################################################################
sub selectNegativeControlProbes{
  my %args = @_;
  my $dbh = $args{'-dbh'};
  my $target_nc_probe_count = $args{'-nc_probe_count'};
  my $filtered = $args{'-filtered'};
  my $filtered_set = $args{'-filtered_set'};

  print BLUE, "\n\nSelecting the desired number ($target_nc_probe_count) of negative control probes, drawing evenly from passing probes of each Tm/Length bin\n", RESET;

  #Get hash of all NC probes from the database
  my $probes_ref;
  if ($filtered =~ /^yes$/i){
    print BLUE, "\n\tGetting all filtered Negative Control probes from ALEXA - may take a while\n", RESET;
    $probes_ref = &getNcProbes ('-dbh'=>$dbh, '-filtered'=>"yes", '-filtered_set'=>$filtered_set);
  }else{
    print BLUE, "\n\tGetting all Unfiltered Negative Control probes from ALEXA -may take a while", RESET;
    print YELLOW" \n\tNOTE: Using unfiltered probes is not recommended ...", RESET;
    $probes_ref = &getNcProbes ('-dbh'=>$dbh, '-filtered'=>"no");
  }


  #Organize all probes according to their probe length and Tm (to 0.1 decimal places)
  my %probe_bins;

  my $potential_candidates = 0;
  foreach my $probe_id (keys %{$probes_ref}){

    $potential_candidates++;
    my $length = $probes_ref->{$probe_id}->{probe_length};
    my $probe_tm = $probes_ref->{$probe_id}->{tm_celsius};
    my $tm = sprintf("%.1f", $probe_tm);

    if ($probe_bins{$length}){
      #Length has been seen before
      my $probe_length_bin_ref = $probe_bins{$length}{probe_length_bins};

      if ($probe_length_bin_ref->{$tm}){
	#Tm has been observed before
	my $probe_tm_bin_ref = $probe_length_bin_ref->{$tm}->{tm_probes};
	$probe_tm_bin_ref->{$probe_id}->{tm} = $probe_tm;

      }else{
	#First probe of this Tm
	my %probe_tm_bin;
	$probe_tm_bin{$probe_id}{tm} = $probe_tm;
	$probe_length_bin_ref->{$tm}->{tm_probes} = \%probe_tm_bin;
      }
    }else{
      #First probe of this length and Tm
      my %probe_tm_bin;
      $probe_tm_bin{$probe_id}{tm} = $probe_tm;

      my %probe_length_bin;
      $probe_length_bin{$tm}{tm_probes} = \%probe_tm_bin;

      $probe_bins{$length}{probe_length_bins} = \%probe_length_bin;
    }
  }

  #Make sure enough passing NC probes were found to accomodate the desired target number
  unless ($potential_candidates > $target_nc_probe_count){
    print RED, "\nNot enough passing negative control probes ($potential_candidates) were found to accomodate the desired target: $target_nc_probe_count", RESET;
    exit();
  }

  #Now go through each length-bin -> tm-bin -> probe_list  and select one probe from each list.
  #Repeat entire process until desired number are chosen
  my %final_nc_probes;
  my $final_nc_probe_count = 0;

  my $iter = 0;
  #Until desired probe count is reached
  ITER:while ($final_nc_probe_count < $target_nc_probe_count){
      $iter++;

      my $length_bin_count = keys %probe_bins;
      if ($verbose eq "yes"){
	print BLUE, "\nIter: $iter.  Going through $length_bin_count length bins of probes", RESET;
      }

      #Go through each length bin
      LENGTH:foreach my $length (sort {$a <=> $b} keys %probe_bins){

	  my $tm_bin_ref = $probe_bins{$length}{probe_length_bins};
	  my $bin_count = keys %{$tm_bin_ref};

	  if ($verbose eq "yes"){
	    print BLUE, "\n\tGoing through $bin_count tm bins of probes", RESET;
	  }

	  #Go through each Tm bin of this length
	TM:foreach my $tm (sort {$a <=> $b} keys %{$tm_bin_ref}){

	    my $tm_probes_ref = $tm_bin_ref->{$tm}->{tm_probes};

	    #Select 1 probe from this bin (make sure it hasnt already been added)
	    foreach my $probe_id (sort {$tm_probes_ref->{$a}->{tm} <=> $tm_probes_ref->{$b}->{tm}} keys %{$tm_probes_ref}){

	      unless($final_nc_probes{$probe_id}){
		$final_nc_probe_count++;
		$final_nc_probes{$probe_id}{tm_celsius} = $tm_probes_ref->{$probe_id}->{tm};
		$final_nc_probes{$probe_id}{nc_probe_count} = $final_nc_probe_count;

		#See if the desired number of probes has been reached
		if ($final_nc_probe_count == $target_nc_probe_count){
		  last LENGTH;
		}else{
		  #otherwise just go to the next Tm bin of probes
		  next TM;
		}
	      }
	    }
	  }
	}
    }

  my $success_count = keys %final_nc_probes;
  print BLUE, "\n\tSelected $success_count negative control probes, drawing evenly from passing probes of each Tm/Length bin\n", RESET;

  #Go through all the negative control probes print them out (i.e. add them to the array) if they were chosen above
  foreach my $probe_id (sort {$a <=> $b} keys %{$probes_ref}){
    unless ($final_nc_probes{$probe_id}){
      next();
    }
    my $type = "Control-Negative";

    #Print probe to submission file
    print SUBMISSION "ALEXA_P_$probe_id\tALEXA_PS_$probes_ref->{$probe_id}->{probe_set_id}\tna\tna\tna\t$probes_ref->{$probe_id}->{sequence}\t$type\n";

    #Print probe record to probe record file
    print PROBES "$probe_id\t$probes_ref->{$probe_id}->{probe_set_id}\tna\t$type\t$probes_ref->{$probe_id}->{sequence}\t$probes_ref->{$probe_id}->{probe_length}\t$probes_ref->{$probe_id}->{tm_celsius}\tna\tna\tna\tna\tna\tna\n";



  }

  return();
}
