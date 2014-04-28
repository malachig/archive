#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Create an annotation database for known 'transcripts' that can be profiled by unique exon regions and exon-exon junctions

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $annotation_dir = '';
my $junction_seq_size = '';
my $outdir = '';

GetOptions('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'annotation_dir=s'=>\$annotation_dir, 'junction_seq_size=i'=>\$junction_seq_size, 'outdir=s'=>\$outdir);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password\n", RESET;
print GREEN, "\n\tSpecify the path to a directory containing exon region and exon junction annotation files using: --annotation_dir", RESET;
print GREEN, "\n\tSpecify the junction database sequence length using:  --junction_seq_size", RESET;
print GREEN, "\n\tSpecify an output directory for the annotation file using:  --outdir", RESET;

print GREEN, "\n\nExample: createTranscriptDatabase.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --annotation_dir=/projects/malachig/sequence_databases/hs_53_36o/  --junction_seq_size=62  --outdir=/projects/malachig/sequence_databases/hs_53_36o/transcripts/\n\n", RESET;

#Check user supplied options
unless ($database && $server && $user && $password && $annotation_dir && $junction_seq_size && $outdir){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}
$annotation_dir = &checkDir('-dir'=>$annotation_dir, '-clear'=>"no");
$outdir = &checkDir('-dir'=>$outdir, '-clear'=>"no");

my $exon_region_file = "$annotation_dir"."exonRegions/exonRegions_annotated.txt.gz";
my $exon_junction_file = "$annotation_dir"."exonJunctions/exonJunctions_"."$junction_seq_size"."mers_annotated.txt.gz";

unless((-e $exon_region_file) && (-e $exon_junction_file)){
  print RED, "\nOne of the required input annotation files was not found!\n\t$exon_region_file\n\t$exon_junction_file\n\n", RESET;
  exit();
}


#1.) Get all the exon region and exon-junction expression values that correspond to specific trancsripts
my $er_ref = &getAnnotationData('-file'=>$exon_region_file);
my $j_ref = &getAnnotationData('-file'=>$exon_junction_file);

#2.) Import gene and transcript info from ALEXA 
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};
my $g_storable_name = "$database"."_AllGenes_GeneInfo_NoSeq.storable";
my $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no", '-storable'=>$g_storable_name);

my $trans_storable = "$database"."_AllGenes_TranscriptInfo_NoSeq.storable";
my $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no", '-storable'=>$trans_storable);

$alexa_dbh->disconnect();

#Initialize value in the transcripts objects
foreach my $g (sort {$a <=> $b} keys %{$genes_ref}){
  my $trans_ref = $gene_transcripts_ref->{$g}->{transcripts};

  foreach my $t (sort {$a <=> $b} keys %{$trans_ref}){
    $trans_ref->{$t}->{base_count} = 0;
    $trans_ref->{$t}->{unmasked_base_count} = 0;
    $trans_ref->{$t}->{coding_base_count} = 0;
    $trans_ref->{$t}->{specific_seq_count} = 0;
    $trans_ref->{$t}->{specific_exon_regions} = 0;
    $trans_ref->{$t}->{specific_exon_junctions} = 0;
    my %tmp;
    $trans_ref->{$t}->{feature_list} = \%tmp;
  }
}

#3-A.) First store the transcript specific sequences within the appropriate transcript object
foreach my $er (keys %{$er_ref}){
  my $gene_id = $er_ref->{$er}->{gene_id};
  my $trans_id = $er_ref->{$er}->{trans_id};
  my $base_count = $er_ref->{$er}->{base_count};
  my $unmasked_base_count = $er_ref->{$er}->{unmasked_base_count};
  my $coding_base_count = $er_ref->{$er}->{coding_base_count};
  my $seq_name = $er_ref->{$er}->{seq_name};
  my $pos = $er_ref->{$er}->{unit1_start};

  my $trans_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};
  $trans_ref->{$trans_id}->{specific_seq_count}++;
  $trans_ref->{$trans_id}->{specific_exon_regions}++;
  $trans_ref->{$trans_id}->{base_count} += $base_count;
  $trans_ref->{$trans_id}->{unmasked_base_count} += $unmasked_base_count;
  $trans_ref->{$trans_id}->{coding_base_count} += $coding_base_count;
  my $feature_list_ref = $trans_ref->{$trans_id}->{feature_list};
  $feature_list_ref->{$seq_name}->{pos} = $pos;
}

foreach my $j (keys %{$j_ref}){
  my $gene_id = $j_ref->{$j}->{gene_id};
  my $trans_id = $j_ref->{$j}->{trans_id};
  my $base_count = $j_ref->{$j}->{base_count};
  my $unmasked_base_count = $j_ref->{$j}->{unmasked_base_count};
  my $coding_base_count = $j_ref->{$j}->{coding_base_count};
  my $seq_name = $j_ref->{$j}->{seq_name};
  my $pos = $j_ref->{$j}->{unit1_start};

  my $trans_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};
  $trans_ref->{$trans_id}->{specific_seq_count}++;
  $trans_ref->{$trans_id}->{specific_exon_junctions}++;
  $trans_ref->{$trans_id}->{base_count} += $base_count;
  $trans_ref->{$trans_id}->{unmasked_base_count} += $unmasked_base_count;
  $trans_ref->{$trans_id}->{coding_base_count} += $coding_base_count;
  my $feature_list_ref = $trans_ref->{$trans_id}->{feature_list};
  $feature_list_ref->{$seq_name}->{pos} = $pos;
}

#3-B.) Now go through each transcript of each gene, and set values to 'NA' if no specific sequences were found
foreach my $g (sort {$a <=> $b} keys %{$genes_ref}){
  my $trans_ref = $gene_transcripts_ref->{$g}->{transcripts};

  my $trans_count = keys %{$trans_ref};

  foreach my $t (sort {$a <=> $b} keys %{$trans_ref}){
    if ($trans_ref->{$t}->{specific_seq_count} == 0){
      #Initialize variables for transcripts that did not have transcript specific elements to use...
      $trans_ref->{$t}->{base_count} = "NA";
      $trans_ref->{$t}->{unmasked_base_count} = "NA";
      $trans_ref->{$t}->{coding_base_count} = "NA";
      $trans_ref->{$t}->{specific_exon_regions} = "NA";
      $trans_ref->{$t}->{specific_exon_junctions} = "NA";
      $trans_ref->{$t}->{feature_list_string} = "NA";
    }else{
      my $feature_list_ref = $trans_ref->{$t}->{feature_list};
      my $feature_count = keys %{$feature_list_ref};
      if ($feature_count > 0){
        my @feature_list;
        foreach my $feature (sort {$feature_list_ref->{$a}->{pos} <=> $feature_list_ref->{$b}->{pos}} keys %{$feature_list_ref}){
          push(@feature_list, $feature);
        }
        $trans_ref->{$t}->{feature_list_string} = join(",", @feature_list);
      }
    }
  }
}

#Print out the summary
#Trans_ID, EnsEMBL_Trans_ID, Gene_ID, EnsEMBL_Gene_ID, Gene_Name, Chromosome, Strand, Unit1_start_chr, Unit1_end_chr, Base_Count, UnMasked_Base_Count, Coding_Base_Count, Specific_Exon_Region_Count, Specific_Junction_Count
my $transcript_output_file = "$outdir"."transcripts_annotated".".txt";
open (TRANS, ">$transcript_output_file") || die "\nCould not open transcript output file: $transcript_output_file\n\n";
print TRANS "Trans_ID\tEnsEMBL_Trans_ID\tGene_ID\tEnsEMBL_Gene_ID\tGene_Name\tChromosome\tStrand\tTranscript_Size\tUnit1_start_chr\tUnit1_end_chr\tBase_Count\tUnMasked_Base_Count\tCoding_Base_Count\tSpecific_Exon_Region_Count\tSpecific_Junction_Count\tSeq_Name\tFeature_List\n";

foreach my $g (sort {$a <=> $b} keys %{$genes_ref}){
  my $trans_ref = $gene_transcripts_ref->{$g}->{transcripts};

  foreach my $t (sort {$a <=> $b} keys %{$trans_ref}){

    #Get chromosome coordinates for this transcript 
    my $start = $trans_ref->{$t}->{transcript_start};
    my $end = $trans_ref->{$t}->{transcript_end};
    my $coords_ref = &convertGeneCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$g, '-start_pos'=>$start, '-end_pos'=>$end, '-ordered'=>"yes");
    my $chr_start = $coords_ref->{$g}->{chr_start};
    my $chr_end = $coords_ref->{$g}->{chr_end};

    print TRANS "$t\t$trans_ref->{$t}->{ensembl_t_id}\t$g\t$genes_ref->{$g}->{ensembl_g_id}\t$genes_ref->{$g}->{gene_name}\t$genes_ref->{$g}->{chromosome}\t$genes_ref->{$g}->{chr_strand}\t$trans_ref->{$t}->{transcript_size}\t$chr_start\t$chr_end\t$trans_ref->{$t}->{base_count}\t$trans_ref->{$t}->{unmasked_base_count}\t$trans_ref->{$t}->{coding_base_count}\t$trans_ref->{$t}->{specific_exon_regions}\t$trans_ref->{$t}->{specific_exon_junctions}\t$trans_ref->{$t}->{ensembl_t_id}\t$trans_ref->{$t}->{feature_list_string}\n";

  }
}
close (TRANS);

#Summarize the total memory usage at close (since Perl doesnt usually release memory ... this should be the max used by the script):
my $pid = $$;
my $ps_query = `ps -p $pid -o pmem,rss`;
my @process_info = split ("\n", $ps_query);
my $memory_usage = '';
my $memory_usage_p = '';
if ($process_info[1] =~ /(\S+)\s+(\S+)/){
  $memory_usage_p = $1;
  $memory_usage = $2;
}
my $memory_usage_m = sprintf("%.1f", ($memory_usage/1024));
print YELLOW, "\n\nMemory usage at end of script: $memory_usage_m Mb ($memory_usage_p%)\n\n", RESET; 

exit();


######################################################################################################################################################
#parse input annotation file and store                                                                                                               #
######################################################################################################################################################
sub getAnnotationData{
  my %args = @_;
  my $file = $args{'-file'};
  
  my %data;

  my %columns;
  my $header = 1;

  my $data_count = 0;

  print BLUE, "\n\nImporting data from: $file", RESET;

  open (FILE, "zcat $file |") || die "\nCould not open input file: $file";
  while(<FILE>){
    my $line = $_;
    chomp($line);
    my @line = split("\t", $line);

    if ($header == 1){
      $header = 0;
      my $pos = 0;
      foreach my $head (@line){
        $columns{$head}{column_position} = $pos;
        $pos++;
      }
      next();
    }

    my $seq_id = $line[0];
    my $base_count = $line[$columns{'Base_Count'}{column_position}];
    my $unmasked_base_count = $line[$columns{'UnMasked_Base_Count'}{column_position}];
    my $coding_base_count = $line[$columns{'Coding_Base_Count'}{column_position}];
    my $gene_id = $line[$columns{'Gene_ID'}{column_position}];
    my $specific_trans_id = $line[$columns{'Specific_Trans_ID'}{column_position}];
    my $seq_name = $line[$columns{'Seq_Name'}{column_position}];
    my $unit1_start = $line[$columns{'Unit1_start'}{column_position}];

    $data_count++;

    if ($specific_trans_id =~ /(\w+)\_(\d+)/){
      my $type = $1;
      my $trans_id = $2;
      $data{$seq_id}{gene_id} = $gene_id;
      $data{$seq_id}{type} = $type;
      $data{$seq_id}{trans_id} = $trans_id;
      $data{$seq_id}{base_count} = $base_count;
      $data{$seq_id}{unmasked_base_count} = $unmasked_base_count;
      $data{$seq_id}{coding_base_count} = $coding_base_count;
      $data{$seq_id}{seq_name} = $seq_name;
      $data{$seq_id}{unit1_start} = $unit1_start;
    }else{
      next();
    }
  }
  close(FILE);

  my $seq_count = keys %data;
  print BLUE, "\nScanned $data_count sequences and found $seq_count transcript specific sequences (exon regions or junctions)\n\n", RESET;

  return(\%data);
}


