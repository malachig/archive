#!/usr/bin/perl

#Full command line command:
#/home/bcho/Perl/GeneJunction/identifyExpressedJunctions.pl --analysis_dir=/projects/alexa2/alexa_seq/ --library_id=HS04391 --gene_name=CCNK --target_ensg_id=ENSG00000090061 --ensembl_version=53 --alexa_db=ALEXA_hs_53_36o --unit_size=10 --min_count=1 --track_color=255,0,0 --results_dir=/home/bcho/Temp/ --test=1

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;

my $script_dir;
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\.pl/){
    push (@INC, $1);
    $script_dir = $1;
  }
}

my $analysis_dir = '';
my $library_id = '';
my $target_ensg_id = '';
my $ensembl_version = '';
my $alexa_db = '';
my $size = '';
my $min_count = '';
my $track_color = '';
my $results_dir = '';
my $test = '';

my $gene_name;
my $chr_start;
my $chr_end;
my $chr_strand;
my $chr;

    #Need ensembl_version parameter (e.g. 53)
    #Need alexa_db parameter (e.g. ALEXA_hs_53_36o)
    #Change '--test_id' to '--target_ensg_id'
    #Change '--lib_id' to 'track name'
    #Change '--color' to '--track_color'
    #Add '--test' option which is the number of files to be processed (this will be passed to getReadID and getReadSeq and use to limit the files processed)
    #Also change error messages to report the script they are coming from
    #Add explanation of each option to the 'unless()' block below
    #Add a '-results_dir' parameter.  Go to that dir, create a new dir named '$library_id', enter that dir, create a new dir named '$gene_name', enter that dir, and write all files there
    #mkdir($results_dir)
    #chdir($results_dir)


GetOptions( 'analysis_dir=s' => \$analysis_dir,
            'library_id=s' => \$library_id,
            'gene_name=s' => \$gene_name,
            'target_ensg_id=s' => \$target_ensg_id,
            'ensembl_version=s' => \$ensembl_version,
            'alexa_db=s' => \$alexa_db,
            'unit_size=i' => \$size,
            'min_count=i' => \$min_count,
            'track_color=s' => \$track_color,
            'results_dir=s' => \$results_dir,
            'test=i' => \$test
          );


unless ($analysis_dir && $library_id && $gene_name && $target_ensg_id && $ensembl_version && $alexa_db && $size && $min_count && $track_color && $results_dir)
{
  print MAGENTA, "Please specify the options, using the following format:\n";
  print "--analysis_dir=/projects/alexa2/alexa_seq/ --library_id=HS04391 --gene_name=CCNK --target_ensg_id=ENSG00000090061 --ensembl_version=53 --alexa_db=ALEXA_hs_53_36o --unit_size=10 --min_count=10 --track_color=255,0,0 --results_dir=/home/bcho/Temp/ --test=1\n\n";

  print "analysis_dir: directory where the libraries are stored\n";
  print "library_id: ID of the library\n";
  print "gene_name: name of the target gene\n";
  print "target_ensg_id: ID of the target ensemble gene\n";
  print "ensembl_version: version of ensembl\n";
  print "alexa_db: alexa database\n";
  print "unit_size: minimum number of bases matched bordering either side of the junction\n";
  print "min_count: minimum count for the number of junctions found\n";
  print "track_color: track color to be displayed in browser\n";
  print "results_dir: directory where the resulting 5 files will be written\n";
  print "test: number of files inside library to test (test=0 or omit option for all files inside library)\n";
  
  print "\nidentifyExpressedJunctions.pl terminated\n\n", RESET;
  exit();
}

unless($test){
  $test=0;
}

unless ($results_dir =~ /\/$/)
{
  #add a slash (/) at the end of the dir string
  $results_dir = "$results_dir"."\/";
}


#Creating and moving to results_dir
my $library_id_dir = "$results_dir"."$library_id/";
my $gene_name_dir = "$results_dir"."$library_id/$gene_name/";
unless (-e $results_dir && -d $results_dir){mkdir($results_dir);}
unless (-e $library_id_dir && -d $library_id_dir){mkdir($library_id_dir);}
unless (-e $gene_name_dir && -d $gene_name_dir){mkdir($gene_name_dir);}


#Goal:
#1.) Identify all the reads that map to a specific gene (or list of genes) in the first set of files
#    - the target gene can be associated with R1, R2 or both
#2.) Store the read IDs for these reads 
#3.) Using these read IDs, retrieve the R1 and R2 sequences from the read records files
#4.) Store them as a fasta file
my $read_seq_file = "$gene_name_dir"."readSeq.fa";
my $getReadID_cmd = "$script_dir/getReadID.pl  --analysis_dir=$analysis_dir  --library_id=$library_id  --gene_name=$gene_name  --ensembl_version=$ensembl_version  --test=$test  |  $script_dir/getReadSeq.pl  --analysis_dir=$analysis_dir  --library_id=$library_id  --test=$test > $read_seq_file";
print BLUE, "\n\nGetting read IDs mapped to $gene_name from ALEXA-Seq ENST mapped reads files ...", RESET;
print YELLOW, "\n$getReadID_cmd", RESET;
system($getReadID_cmd);

#5.) Retrieve the corresponding gene sequence (like the CAV32 genomic sequence we have been mapping reads to)
my $getGeneSeq_cmd = "$script_dir/getGeneSeq.pl  --target_ensg_id=$target_ensg_id  --alexa_db=$alexa_db  --working_dir=$gene_name_dir";
print BLUE, "\n\nGetting complete genome sequence for $gene_name from ALEXA-DB $alexa_db ...", RESET;
print YELLOW, "\n$getGeneSeq_cmd", RESET;
system($getGeneSeq_cmd);

my $gene_info_file = "$gene_name_dir"."gene_info.txt";
open (GENE_INFO, "$gene_info_file");
my $n = 1;
while(my $line = <GENE_INFO>)
{
  #read in the gene information bits from getGeneSeq.pl
  chomp $line;
  if ($n == 1)
  {
    $gene_name = $line;
    $n++; next();
  }
  if ($n == 2)
  {
    $chr_start = $line;
    $n++; next();
  }
  if ($n == 3)
  {
    $chr_end = $line;
    $n++; next();
  }
  if ($n == 4)
  {
    $chr_strand = $line;
    $n++;
    if ($chr_strand == 1)
    {
      $chr_strand = '+';
    }
    else
    {
      $chr_strand = '-';
    }
    next();
  }
  if ($n == 5)
  {
    $chr = $line;
  }
}
close (GENE_INFO);


#6.) Use spidey to map all reads in the fasta file to the genomic sequence obtained in the previous step
#7.) Summarize the junctions
my $summary_file = "$gene_name_dir"."summary.txt";
my $gene_seq_file = "$gene_name_dir"."$gene_name"."_seq.fa";
my $spidey_cmd = "$script_dir/spidey.linux -i $gene_seq_file -m $read_seq_file -p 1 2>/dev/null  |  $script_dir/junction.pl  --size=$size | $script_dir/convertCoords.pl  --strand=$chr_strand  --chr_start=$chr_start  --chr_end=$chr_end > $summary_file";
print BLUE, "\n\nUsing spidey to align reads to gene sequence ...", RESET;
print YELLOW, "\n$spidey_cmd", RESET;
system($spidey_cmd);

#8.) Produce a GFF file
my $track_name = "$gene_name"."_"."$library_id"; 
my $gff_file = "$gene_name_dir"."$track_name.gff";
my $makeGFF_cmd = "cat $summary_file | $script_dir/makeGFF.pl  --strand=$chr_strand  --chr=$chr  --min_count=$min_count  --track_name=$track_name  --track_color=$track_color > $gff_file";
print BLUE, "\n\nGenerating a GFF file to display the results ...", RESET;
print YELLOW, "\n$makeGFF_cmd", RESET;
system($makeGFF_cmd);

print MAGENTA, "\n\nAll results stored here: $gene_name_dir", RESET;

print "\n\nSCRIPT COMPLETE\n\n";

exit();


