#!/usr/bin/perl -w
#Written by Malachi Griffith

#Parses all ace files for assemblies in a parent directory and generates a fasta files for the contigs from the most recent Ace file of each assembly

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

use lib '/home/malachig/svn/clone_assembly/';
use lib "/home/malachig/perl/bioperl-1.4";    #Bioperl - used by ACE3_GSC
use utilities::ACE3_GSC;

#Initialize command line options
my $assembly_dir = '';
my $fasta_file = '';

GetOptions ('assembly_dir=s'=>\$assembly_dir, 'fasta_file=s'=>\$fasta_file);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the top directory containing clone assembly directories: --assembly_dir", RESET;
print GREEN, "\n\nExample:  parseAceFileToFasta.pl  --assembly_dir=/home/malachig/AlternativeSplicing/UMPS_Clone_Assemblies/  --fasta_file=CloneFasta.txt\n\n", RESET;

#Make sure all neccessary parameters were specified
unless ($assembly_dir && $fasta_file){
  print RED, "\nRequired input parameter(s) missing or incorrect format!\n\n", RESET;
  exit();
}

#Get clone assembly directories with a specified parent directory - make note of the most recent ACE file (by build number)
my $clone_assembly_dirs_ref = &getCloneAssemblyDirs('-assembly_dir'=>$assembly_dir);

#parse the ace files collected in the previous step and write a single fasta file
&generateFastaFile('-clone_assemblies'=>$clone_assembly_dirs_ref, '-fasta_file'=>$fasta_file);


#print Dumper $clone_assembly_dirs_ref;





exit();



####################################################################################################################################
#Get clone assembly directories with a specified parent directory - grab the most recent ACE file (by build number)
####################################################################################################################################
sub getCloneAssemblyDirs{
  my %args = @_;
  my $assembly_dir = $args{'-assembly_dir'};

  my %dirs;

  #Add trailing '/' to directory paths if they were forgotten
  unless ($assembly_dir =~ /.*\/$/){
    $assembly_dir = "$assembly_dir"."/";
  }

  #Check the validity of required directories specified by the user
  unless (-d $assembly_dir && -e $assembly_dir){
    print RED, "\n\nInput assembly directory does not appear to be a valid directory:\n\tassembly_dir = $assembly_dir\n", RESET;
    exit();
  }


  print BLUE, "\n\nChecking specified input directory for clone assemblies ...\n", RESET;

  opendir(DIRHANDLE, "$assembly_dir") || die "\nCannot open directory: $assembly_dir\n\n";
  my @files = readdir(DIRHANDLE);

  my $count = 0;
  foreach my $file (@files){

    my $path = "$assembly_dir"."$file"."/";

    unless (-d $path){
      print YELLOW, "\n\t$path is not a clone assembly directory - skipping", RESET;
      next();
    }
    if ($file eq "." || $file eq ".."){
      print YELLOW, "\n\t$path is not a clone assembly directory - skipping", RESET;
      next();
    }

    #Make sure this path contains a edit_dir!
    my $edit_dir_path = "$path"."edit_dir"."/";

    unless (-e $edit_dir_path && -d $edit_dir_path){
      print YELLOW, "\n\t$edit_dir_path does not contain an edit_dir - skipping", RESET;
      next();
    }

    $count++;
    $dirs{$count}{clone_name} = $file;
    $dirs{$count}{assembly_path} = $path;
    $dirs{$count}{edit_dir_path} = $edit_dir_path;

    #Get list of builds (ace files) from assembly directory
    my @ace_files = `ls $edit_dir_path*.ace.*`;
    chomp @ace_files;

    my %ace_files;

    my $current_build = 0;

    #Use regex to grab build numbers and verify ace file names
    foreach my $ace_file (@ace_files){
      if ($ace_file =~ /.*\.fasta\.screen\.ace\.(\d+)$/){
	my $build = $1;
	if ($build > $current_build){
	  $current_build = $build;
	}
	$ace_files{$build}{file_path} = $ace_file;
      }
    }

    $dirs{$count}{ace_files} = \%ace_files;
    $dirs{$count}{current_build} = $current_build;
  }
  print BLUE, "\nFound $count clone assembly directories\n\n", RESET;

  return(\%dirs);
}


####################################################################################################################################
#parse the ace files collected in the previous step and write a single fasta file
####################################################################################################################################
sub generateFastaFile{
  my %args = @_;
  my $fasta_file = $args{'-fasta_file'};
  my $clone_assembly_obj = $args{'-clone_assemblies'};

  print BLUE, "\n\nWriting clone sequence from most recent ace files to: $fasta_file\n\n", RESET;

  open (FASTA, ">$fasta_file") || die "\nCould not open fasta output file: $fasta_file\n\n";

  foreach my $clone_count (sort {$a <=> $b} keys %{$clone_assembly_obj}){

    my $current_build = $clone_assembly_obj->{$clone_count}->{current_build};
    my $clone_name = $clone_assembly_obj->{$clone_count}->{clone_name};

    print YELLOW, "\n$clone_count\tClone: $clone_name\tBuild: $current_build", RESET;

    my $acefiles_ref = $clone_assembly_obj->{$clone_count}->{ace_files};
    my $acefile = $acefiles_ref->{$current_build}->{file_path};

    #Set parameters for the AC3 ace object
    my @params = ('acefile' => "$acefile", 'ace_verbose' => 1);

    #Instantiate and AC3 ace object
    my $aceObj = Bio::Tools::ACE3->new(@params);

    #Get an array with the contig names in order to use accessor methods 
    my @con_names_ref = @{$aceObj->get_contigs()};

    #Get the total Number of Contigs in this acefile
    my $total_contigs = $aceObj->total_number_contigs();

    if ($total_contigs == 0){
      print RED, "\nAce File has zero contigs\n\n", RESET;
      return (0);
    }elsif ($total_contigs > 1){
      print YELLOW, "\nClone $clone_name has more than one contig!\n", RESET;
    }

    foreach my $contig_number (@con_names_ref){
      my $num_reads = $aceObj->number_reads($contig_number);

      #Get contig length
      my $length = $aceObj->contig_length($contig_number);

      #Get sequence object (pads removed).
      my $consensus_seq_obj = $aceObj->consensus_Seq_Obj($contig_number);
      my $contig_sequence = $consensus_seq_obj->seq();

      #Get raw sequence consensus scores as a string
      my $qualities = $aceObj->quality($contig_number);

      #Get raw sequence (with pads)
      my $padded_sequence = $aceObj->consensus($contig_number);

      #Get the average phred for this contig
      my $average_phred = $aceObj->contig_average_phred($contig_number);

      #Get the error_rate (per 10 kb) for this contig
      my $error_rate = $aceObj->contig_error_rate($contig_number);

      #Get the lowest phred score for this contig (ie. the phred value of the worst base)
      my $lowest_score = $aceObj->contig_lowest_score($contig_number);

      #Clean up contig sequence to remove 'x's
      my $clean_seq = $contig_sequence;
      $clean_seq =~ s/x//g;
      $clean_seq =~ s/X//g;

      my $new_length = length($clean_seq);

      print FASTA ">$clone_name (CONTIG=$contig_number  READS=$num_reads  LENGTH=$new_length  AVG_PHRED=$average_phred  ERROR_RATE=$error_rate  LOWEST_SCORE=$lowest_score)\n$clean_seq\n";
    }

  }

  close (FASTA);

  print "\n\n";

  return();
}
