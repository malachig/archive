#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use File::Copy;

my $reference_fasta = "/storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa";
my $sample = '';
my $alignfile = '';
my $outdir = '';

GetOptions ('sample=s'=>\$sample, 'alignfile=s'=>\$alignfile, 'outdir=s'=>\$outdir);

#Make sure parameters were provided
unless($sample && $alignfile && $outdir){
  print STDERR "\n\nRequired parameters missing\n\n";
  print STDERR "Usage:  organizeBams.pl --sample='BrMET008.1.Tumor.DNA' --alignfile='/storage1/fs1/gpdunn/Active/Project_0001_Clinical/gc2609/gbm_antigen_dunn/model_data/4aeb78cb9b434b72852ffe748361657d/build6c3ff930322e4155957d4f19e40b5526/results/tumor.cram' --outdir='/scratch1/fs1/mgriffit/gbm_submission'\n\n";
  exit();
}

#Make sure reference_fasta, alignfile and outdir exist
unless (-e $outdir && -d $outdir){
  print STDERR "\n\nOutdir ($outdir) could not be found or is not a directory\n\n";
  exit();
}
unless (-e $alignfile){
  print STDERR "\n\nAlignment file ($alignfile) could not be found\n\n";
  exit();
}
unless (-e $reference_fasta){
  print STDERR "\n\nReference fasta file ($reference_fasta) could not be found\n\n";
  exit();
}

#Make sure outdir path has a trailing '/'
unless ($outdir =~ /\/$/){
  $outdir .= "/";
}

#1. Create a short description like: Exome capture sequence data for DNA isolated from region X of tumor sample Y
my $sample_base_name = ''; #GBM070
my $sample_type = ''; #Tumor/Normal 
my $nucleic_acid_type = ''; #DNA/RNA
my $region = ''; #1,2,3,4
my $data_type = ''; #Exome/RNA-Seq

if ($sample =~ /^(\S+)\./){
  $sample_base_name = $1;
}
if ($sample =~ /tumor/i){$sample_type = "Tumor";}
if ($sample =~ /normal/i){$sample_type = "Normal";}
if ($sample =~ /dna/i){$nucleic_acid_type = "DNA"; $data_type = "Exome";}
if ($sample =~ /rna/i){$nucleic_acid_type = "RNA"; $data_type = "RNA-Seq";}
if ($sample =~ /(\d+)\.tumor/i){$region = $1;}
if ($sample =~ /normal/i){$region = "N/A";}

my $description = $data_type . " data generated for " . $nucleic_acid_type . " isolated from region " . $region . " of " . $sample_type . " sample " .  "$sample";

#Make sure all elements were resolved
unless ($sample_base_name && $sample_type && $nucleic_acid_type && $region && $data_type){
  print STDERR "\n\nFailed to resolve information from sample_name: $sample\nsample_base_name = $sample_base_name\nsample_type = $sample_type\nnucleic_acid_type = $nucleic_acid_type\nregion = $region\ndata_type = $data_type\n\n";
  exit();
}
print "Data Description:\n$description\n\n";

#2. Create a design description for exome and RNA-seq data
my $design_description = '';
if ($data_type eq "Exome"){
  $design_description = "Exome hybrid capture data was created using the IDT exome reagent: xGen Exome Research Panel v2";
}elsif($data_type eq "RNA-Seq"){
  $design_description = "RNA-seq data was generated from total RNA subjected to ribosomal RNA reduction followed by a strand-specific library construction approach";
}
print "Design Description:\n$design_description\n\n";

#3. If the alignment is in CRAM format, convert to BAM, otherwise copy it to the outdir.
#   Either way name it using the sample name
my $new_filename = $sample . ".bam";
my $new_path = $outdir . $new_filename;

#If the new file is already in the new location use it:
if (-e $new_path){
  print "\nUsing existing bam file already in the target location:\n$new_path\n\n";
}else{
  if($alignfile =~ /.*\/(\S+)\.bam$/){
    print "\ncp $alignfile $new_path\n\n";
    copy($alignfile, $new_path);

  }elsif($alignfile =~ /.*\/(\S+)\.cram$/){
    #Convert cram to BAM
    my $convert_command = "samtools view -b -T $reference_fasta -o $new_path $alignfile";
    print "\n$convert_command\n\n"; 
    system($convert_command);

  }else{
    print STDERR "\n\nFile extension not recognized for file:\n$alignfile\n\n";
    exit();
  }
}

#4. Calculate MD5sum for the BAM file and parse the result

#If the md5 is already calculated use that file
my $md5_outfile = $outdir . $sample . ".md5sum";
my $md5 = '';
if (-e $md5_outfile){
  print "\nUsing existing bam file already in the target location:\n$md5_outfile\n\n";
}else{
  my $md5_command = "md5sum $new_path > $md5_outfile";
  print "\n$md5_command\n\n";
  system($md5_command);

  open(my $fh, '<:encoding(UTF-8)', $md5_outfile) or die "Could not open file '$md5_outfile' $!";
  while (my $row = <$fh>) {
    chomp $row;
    if ($row =~ /^(\S+)/){
      $md5 = $1;
    }else{
      print STDERR "\n\nCould not extract md5 result from file: $md5_outfile";
      exit();
    }
    print "$row\n";
  }
  close $fh;
}

#5. Output a record line with each element needed for the dbGaP submission form and store in an output file:
#Sample, Data Description, Design Description, filename, filepath, MD5_checksum
my $record_file = $outdir . $sample . ".dbgap_record.txt";
my $record_line = "$sample\t$description\t$design_description\t$new_filename\t$new_path\t$md5";
print "\n$record_line\n\n";
open(FH, '>', $record_file) or die $!;
print FH "$record_line\n";
close(FH);

print "\nBAM FILE PROCESSING FOR $sample COMPLETED\n\n";

exit;
