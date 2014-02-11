#!/usr/bin/perl -w
#Written by Malachi Griffith

#This script takes an Affy probeset-to-ensembl_gene map file as input attempts to refine the mappings so that only core probesets remain
#The basic mapping was done by simply checking each probeset against the boundaries of the gene locus.  To be considered a 'CORE' probeset,
#it must lie completely within the boundaries of a known Ensembl exon within the target gene.

#1.) Parse the input mapfile and get probeset-to-ensembl-gene mappings
#NOTE: For now use the coordinates for the entire probeset.  This means that all probes in the probe set must be completely
#      bounded by an Ensembl exon.  Probeset data is provided for both NCBI build 34 and 35.  Individual probe coordinates are only provided for build 34
#     - To use the probe coordinates will therefore require converting the coordinates.

#2.)  Get a hash of all ALEXA genes and their gene coordinates from the ALEXA database
#     - Get the exon coordinates for each gene (eliminate duplicate exons)

#2a.) For each probeset, get the coordinates (chr, strand, start, stop).  Make note of the number of probes in each probeset
#2b.) Attempt to map each probeset to an Ensembl exon
#2b.) Some probesets may map to multiple exons (in the case of overlaping exon with different start/end coords)
#     - As long as a probeset maps completely within at least one known Ensembl exon, allow it into the core set
#3.)  Create an output file with the following data:
#     - Probe_set_ID, EnsemblGeneID, ProbeCount
#     - A seperate script: createEnsemblMetaProbesetFile.pl will be used to convert this file to a meta-probeset file to be used with ExACT


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
my $map_file = '';
my $out_file = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'map_file=s'=>\$map_file, 'out_file=s'=>\$out_file);

#Provide instruction to the user
print BLUE, "\n\nUsage:", RESET;
print BLUE, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print BLUE, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print BLUE, "\n\tSpecify the input probeset-to-ensembl map file using: --map_file", RESET;
print BLUE, "\n\tSpecify the output file using: --out_file", RESET;
print BLUE, "\n\nExample: identifyCoreEnsemblProbesets.pl  --database=ALEXA_hs_35_35h  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --map_file=probesets_mapped_to_ensembl_hs_35_35h.txt  --out_file=CORE_probesets_mapped_to_ensembl_hs_35_35h.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $map_file && $out_file){
  print RED, "\nOptions missing!\n\n", RESET;
  exit();
}

#1.) Parse input map file containing probeset-to-Ensembl genes mappings
#    For each probeset, get the coordinates (chr, strand, start, stop) for each of the probes.  Make note of the number of probes in each probeset
#    - To save memory, immediately check the coordinates and only store those that correspond to EnsEMBL exons
my %gene_probeset_map;
&parseMapFile('-map_file'=>$map_file);

#2.) Get additional info for all ALEXA genes found in the mapfile
#    - Just need exon coordinates for non-redundant exons of each gene

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

#For each gene, go through all the probes of that gene and get genomic coordinates
&getGeneExonData('-dbh'=>$alexa_dbh,'-gene_probeset_map'=>\%gene_probeset_map);

#Close database connection
$alexa_dbh->disconnect();

#print Dumper %gene_probeset_map;

#3.) Actually do the mapping required and create an output file with the following data:
#     - Probe_set_ID, alexa_gene_id, EnsemblGeneID, ProbeCount
&mapProbesetsToExons('-gene_probeset_map'=>\%gene_probeset_map, '-out_file'=>$out_file);

exit();


#######################################################################################################################
#1.) Parse input map file containing probeset-to-Ensembl genes mappings                                               #
#######################################################################################################################
sub parseMapFile{
  my %args = @_;
  my $map_file = $args{'-map_file'};

  #Populate hash: %gene_probeset_map;

  my $record_count = 0;
  my $gene_count = 0;

  print GREEN, "\nBegin importing Gene-to-Probeset mappings from: $map_file\n\n", RESET;

  open (MAPFILE, "$map_file") || die "\nCould not open mapfile: $map_file\n\n";

  while(<MAPFILE>){
    my $line = $_;
    chomp($line);

    my @line = split ("\t", $line);

    if ($line[0] =~/^\d+/){
      $record_count++;

      my $probeset_id = $line[0];
      my $alexa_gene_id = $line[1];
      my $ensembl_gene_id = $line[2];
      my $transcript_cluster_id = $line[3];
      my $chromosome = $line[4];
      my $strand = $line[5];
      my $start = $line[6];
      my $stop = $line[7];
      my $probe_count = $line[8];

      if ($gene_probeset_map{$alexa_gene_id}){
	
	#Grab a copy of the hash reference stored in the gene_probe_map hash
	#Use this reference to directly add to the probeset list of this gene
	my $probesets_ref = $gene_probeset_map{$alexa_gene_id}{probesets};
	
	$probesets_ref->{$probeset_id}->{alexa_gene_id} = $alexa_gene_id;
	$probesets_ref->{$probeset_id}->{ensembl_gene_id} = $ensembl_gene_id;
	$probesets_ref->{$probeset_id}->{transcript_cluster_id} = $transcript_cluster_id;
	$probesets_ref->{$probeset_id}->{chromosome} = $chromosome;
	$probesets_ref->{$probeset_id}->{strand} = $strand;
	$probesets_ref->{$probeset_id}->{start} = $start;
	$probesets_ref->{$probeset_id}->{stop} = $stop;
	$probesets_ref->{$probeset_id}->{probe_count} = $probe_count;
	
      }else{
	#Gene has not been seen before
	$gene_count++;
	my %probesets;
	$probesets{$probeset_id}{alexa_gene_id} = $alexa_gene_id;
	$probesets{$probeset_id}{ensembl_gene_id} = $ensembl_gene_id;
	$probesets{$probeset_id}{transcript_cluster_id} = $transcript_cluster_id;
	$probesets{$probeset_id}{chromosome} = $chromosome;
	$probesets{$probeset_id}{strand} = $strand;
	$probesets{$probeset_id}{start} = $start;
	$probesets{$probeset_id}{stop} = $stop;
	$probesets{$probeset_id}{probe_count} = $probe_count;
	$gene_probeset_map{$alexa_gene_id}{probesets} = \%probesets;
      }
    }
  }
  close (MAPFILE);
  print GREEN, "\nProcessed: $record_count probeset records corresponding to $gene_count genes\n\n", RESET;

  return();
}

#######################################################################################################################
#2.) Get additional info for all ALEXA genes found in the mapfile                                                     #
#######################################################################################################################
sub getGeneExonData{
  my %args = @_;
  my $dbh = $args{'-dbh'};
  my $gene_probeset_map_ref = $args{'-gene_probeset_map'};

  print GREEN, "\nGetting exon data for each gene from ALEXA\n\n", RESET;

  my @gene_ids = keys %{$gene_probeset_map_ref};
  my $gene_exons_ref = &getExons ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

  foreach my $gene_id (keys %{$gene_probeset_map_ref}){
    #Get the UNIQUE exons for a particular gene ID (all exons from all transcripts)
    my $exons = $gene_exons_ref->{$gene_id}->{exons};

    #Assemble a reference set of exons (a superset of all non-redundant exons)
    my %reference_exons;

    foreach my $exon_id (sort {$a <=> $b} keys %{$exons}){
      my $exon_start = $exons->{$exon_id}->{exon_start};
      my $exon_end = $exons->{$exon_id}->{exon_end};

      #Check each of the reference exons to see if one of these is the same, otherwise add it to the list
      my $redundant_exon = 0;
      foreach my $ref_exon_id (keys %reference_exons){
	if ($exon_start == $reference_exons{$ref_exon_id}{exon_start} && $exon_end == $reference_exons{$ref_exon_id}{exon_end}){
	  $redundant_exon = 1;
	}
      }
      #Unless the current transcript exon was found to be redundant, add it to the list
      unless ($redundant_exon == 1){

	#First get the corresponding chromosome level coordinates
	my %coords = %{&convertGeneCoordinates ('-dbh'=>$dbh, '-gene_id'=>$gene_id, '-start_pos'=>$exon_start, '-end_pos'=>$exon_end)};
	my $chr_start = $coords{$gene_id}{chr_start};
	my $chr_end = $coords{$gene_id}{chr_end};

	$reference_exons{$exon_id}{exon_start} = $chr_start;
	$reference_exons{$exon_id}{exon_end} = $chr_end;
      }
    }
    $gene_probeset_map_ref->{$gene_id}->{exons} = \%reference_exons;
  }

  return();
}


#######################################################################################################################
#3.) Actually do the mapping required and create an output file with the following data:                              #
#    - Probe_set_ID, alexa_gene_id, EnsemblGeneID, ProbeCount                                                         #
#######################################################################################################################
sub mapProbesetsToExons{
  my %args = @_;
  my $gene_probeset_map_ref = $args{'-gene_probeset_map'};
  my $out_file = $args{'-out_file'};

  my $probe_record_count = 0;
  my $core_probes_found = 0;
  my $gene_count = 0;

  print GREEN, "\nAttempting to map each probeset to an Ensembl exon from its gene\n\n", RESET;

  open (OUTFILE, ">$out_file") || die "\nCould not open output file: $out_file\n\n";

  print OUTFILE "probe_set_id\talexa_gene_id\tensembl_gene_id\ttranscript_cluster_id\tchromosome\tstrand\tstart\tstop\tprobe_count\n";


  foreach my $gene_id (sort {$a <=> $b} keys %{$gene_probeset_map_ref}){
    $gene_count++;

    #Get the probesets and exons for this gene (by reference only!)
    my $probesets_ref = $gene_probeset_map_ref->{$gene_id}->{probesets};
    my $exons_ref = $gene_probeset_map_ref->{$gene_id}->{exons};

    #Now go through each probeset and compare its coordinates to all the exons
    PROBESET:foreach my $probeset_id (sort {$a <=> $b} keys %{$probesets_ref}){
	$probe_record_count++;

	foreach my $exon (sort {$a <=> $b} keys %{$exons_ref}){

	  #First make sure the probes are on the right strand.
	  #Since this script is starting with probesets mapped unambiguously to Ensembl genes, we dont need to worry about strand
	  #The exons will be on the same strand as the gene each probeset was mapped to (obviously)
	  my $p_start = $probesets_ref->{$probeset_id}->{start};
	  my $p_end = $probesets_ref->{$probeset_id}->{stop};
	  my $e_start = $exons_ref->{$exon}->{exon_start};
	  my $e_end = $exons_ref->{$exon}->{exon_end};

	  if ($p_start >= $e_start && $p_start <= $e_end && $p_end >= $e_start && $p_end <= $e_end){
	    $core_probes_found++;

	    #Found a matching exon for this probeset
	    print BLUE, "\n$probe_record_count: $probeset_id ($probesets_ref->{$probeset_id}->{chromosome}:$p_start-$p_end $probesets_ref->{$probeset_id}->{strand}) maps to an exon of $gene_id", RESET;

	    print OUTFILE "$probeset_id\t$gene_id\t$probesets_ref->{$probeset_id}->{ensembl_gene_id}\t$probesets_ref->{$probeset_id}{transcript_cluster_id}\t$probesets_ref->{$probeset_id}->{chromosome}\t$probesets_ref->{$probeset_id}->{strand}\t$p_start\t$p_end\t$probesets_ref->{$probeset_id}->{probe_count}\n";



	    #Only one exon match is required, skip to the next probeset
	    next PROBESET;
	}
      }
    }
  }
  print GREEN, "\nProcessed $probe_record_count corresponding to $gene_count Ensembl genes\n\n", RESET;
  print GREEN, "\nFound $core_probes_found core probes that map within an Ensembl exon\n\n", RESET;

  close (OUTFILE);

  return();
}
