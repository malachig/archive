#!/usr/bin/perl -w
# Malachi Griffith

#Overview
#Identify mutations in individual libraries associated with splicing events in those libraries
#1.) Look for mutation in library X where there is a corresponding exon skip event that flanks the mutation position
#    - These are junctions of class ('DA', 'NDA', 'D', or 'A')
#    - They must also be an exon skip of at least 1
#2.) Look for 'D' or 'A' junctions that are within X bases of a mutation position
#    - Ideally these should be anchored to an exon that is normally connected to the Donor or Acceptor at the mutated position
#3.) Create a BED file with all mutant positions for visualization purposes
#    - Manually check mutations associated with a junction identified in (1) or (2) against the junctions identified as expressed in that library...

#Notes
#One mutation may be contained within multiple exon-skipping events... 
#It is perhaps unwise to exclude known exon skipping junctions because known transcript annotations may contain aberrant exon skipping events
#The known vs. novel status should be noted though, as this will influence interpretation of the results

#Output
#Mutation_Pos, Gene Name, Junction ID, Junction class, Number exons skipped, Correlation type (contained vs. proximal), distance, number of libs with mutation, number of libs with junction, list of libs with the mutation, list of libs with the junction, list of KNOWN junctions associated with the mutation (i.e. junction for which it is actually a splice site mutation)

#Steps
#1.) Get libary to library name mappings
#2.) Get a list of file pairs to be considered.   The mutation file and junction for a single library
#3.) Get the junction annotations for the whole project
#4.) For each file pair, store the list of mutation positions and junctions observed
#    - Only store junctions that are 'DA', 'NDA', 'D', or 'A'
#5.) Identify the known junction(s) each mutation is associated with.  Use this association to infer the strand that is relevant to the splice site mutation
#    - Make note of whether the splice site mutation is occuring at a Donor or Acceptor site
#    - A splice site mutation at Donor will be on the left side for +ve strand junctions and on the right side for -ve strand junctions
#    - A splice site mutation at Acceptor will be on the right side for +ve strand junctions and on the left side for -ve strand junctions

#6.) Look for overlaps between all the junctions stored for the library and all the mutation positions of the library
#    - Use BEDTools to identify overlaps?
#7.) Store the overlaps (key on mutation_pos+junction_pos)
#    - Keep a count of how many libraries have the same overlap, also how many libraries have the mutation, and how many libraries have the junction
#8.) Print out a summary file for all comparisons with the Output values listed above

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;
use Benchmark;

my $ref_junctions_file = '';
my $mutant_dir = '';
my $bedtools_bin_dir = '';
my $analysis_dir = '';
my $project = '';
my $outdir = '';
my $max_distance = '';

GetOptions ('ref_junctions_file=s'=>\$ref_junctions_file, 'mutant_dir=s'=>\$mutant_dir, 'bedtools_bin_dir=s'=>\$bedtools_bin_dir, 'analysis_dir=s'=>\$analysis_dir, 'project=s'=>\$project, 'outdir=s'=>\$outdir,
            'max_distance=i'=>\$max_distance);

if ($ref_junctions_file && $mutant_dir && $bedtools_bin_dir && $analysis_dir && $project && $outdir && $max_distance){
  #Parameters ok
}else{
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  print GREEN, "\nSpecify a file containing reference junctions using: --ref_junctions", RESET;
  print GREEN, "\nSpecify the directory containing library by library splice site mutation files using:  --mutant_dir", RESET;
  print GREEN, "\nSpecify the path to the BEDTools binary dir using: --bedtools_bin_dir", RESET;
  print GREEN, "\nSpecify the base analysis dir using: --analysis_dir", RESET;
  print GREEN, "\nSpecify the project name using: --project", RESET;
  print GREEN, "\nSpecify the number of bases on either side of the mutation position to scan for alternative donor/acceptor site usage using:  --max_distance", RESET;
  print GREEN, "\nSpecify the name of the output file using: --outdir", RESET;
  print GREEN, "\n\nExample: correlateSNV_vs_Splice.pl  --ref_junctions=/projects/alexa2/hmmSplicer/ReferenceAnnotations/hg18/ALL.junc  --mutant_dir=/projects/rgoyaprj2/projects/breast_cancer/WGA/EXCAP/bb_0.2/  --bedtools_bin_dir=/home/malachig/tools/BEDTools-Version-2.10.1/bin/  --analysis_dir=/projects/alexa2/hmmSplicer/  --project=SA_TN_Breast  --outdir=/projects/alexa2/hmmSplicer/SA_TN_Breast/summary/spliceSiteMutations/  --max_distance=50000\n\n", RESET;
  exit();
}
unless (-e $ref_junctions_file){
  print RED, "\n\nCould not find reference junctions file: $ref_junctions_file\n\n", RESET;
  exit();
}
unless ($mutant_dir =~ /\/$/){
  $mutant_dir .= "/";
}
unless (-e $mutant_dir && -d $mutant_dir){
  print RED, "\n\nMutant dir: $mutant_dir not found\n\n", RESET;
  exit();
}
unless ($analysis_dir =~ /\/$/){
  $analysis_dir .= "/";
}
unless (-e $analysis_dir && -d $analysis_dir){
  print RED, "\n\nAnalysis dir: $analysis_dir not found\n\n", RESET;
  exit();
}
unless ($outdir =~ /\/$/){
  $outdir .= "/";
}
unless (-e $outdir && -d $outdir){
  print RED, "\n\nAnalysis dir: $outdir not found\n\n", RESET;
  exit();
}
unless ($bedtools_bin_dir =~ /\/$/){
  $bedtools_bin_dir .= "/";
}
unless (-e $bedtools_bin_dir && -d $bedtools_bin_dir){
  print RED, "\n\nBEDTools binary dir: $bedtools_bin_dir not found\n\n", RESET;
  exit();
}


#Get the library to library name mappings from the classes file
#Also gather file pairs at this time (one mutation file and junction file per library)
my $class_file = "$analysis_dir"."$project/jobs/LibraryClasses.txt";
my $libs_ref = &getLibInfo('-lib_file'=>$class_file);


#Store the junction annotations for this project
#Only store junctions of type 'DA', 'NDA', 'D' and 'A'

my @projects;
push(@projects, $project);
my $junc_anno = &importJunctionAnnotations('-projects'=>\@projects, '-analysis_dir'=>$analysis_dir);


#Store the mutations and junctions for each library
#Only store junctions that were defined in the annotated list
my %lib_mutations;
my %lib_junctions;
&importLibraryData('-libs'=>$libs_ref);
#print Dumper %lib_mutations;

#Determine the known junction associated with each mutation and use this info to infer the strand
my $anchored_mutations = &getMutatedJunctions('-infile'=>$ref_junctions_file, '-mutations'=>\%lib_mutations);

#Determine the overlap between all stored junction events and splice site mutations
#Only consider the mutations that were successfully anchored in the previous step...
#The junction must be anchored to at least one of the known junctions (on the same strand) associated with the mutation identified in the previous step
#Possible outcomes for expressed/observed junctions compared to mutation positions:
#1.) No overlap (flanking or proximal) - simply ignore these junctions
#2.) The overlaping junctions found may be one of the known junctions already found, in this case, note the libraries and read counts
#3.) A putative alternate junction on the same strand and within k bases of the mutation or flanking it completely... Store these as alternative junctions associated with the mutation.  Note the libraries and read counts
&getAlternateJunctions('-junctions'=>\%lib_junctions);
#print Dumper $anchored_mutations;

#Print output files.  Create three output files
#1.) A main results file
#2.) A short list that has been filtered (on the mutated vs. non-mutated ratio for example)
#3.) A file to feed into R for figure generation (mutation id, gene name, mutated lib ids, known junction ids, alternative junction ids)
&printOutputFiles('-dir'=>$outdir, '-data'=>$anchored_mutations);

print "\n\n";

exit();

##################################################################################################################################################
#Get lib ids and names                                                                                                                           #
##################################################################################################################################################
sub getLibInfo{
  my %args = @_;
  my $lib_file = $args{'-lib_file'};

  print BLUE, "\n\nGathering mutation and detected junction files for each library in: $lib_file", RESET;

  #Get the lib ids and names
  my %libs1;
  open(LIBS, "$lib_file") || die "\n\nCould not open lib info file: $lib_file\n\n";
  my $header = 1;
  while(<LIBS>){
    if ($header == 1){
      $header = 0;
      next();
    }
    chomp($_);
    my @line = split("\t", $_);
    $libs1{$line[0]}{name}=$line[1];
  }
  close(LIBS);
  #print Dumper %libs1;

  #Check for mutation and junction files and return the subset of libs that has both
  my %libs2;
  foreach my $lib (keys %libs1){
    my $name = $libs1{$lib}{name};

    #Check for a mutation file using either the lib_id or name
    my $mf1 = "$mutant_dir"."$lib".".somatics.tsv";
    my $mf2 = "$mutant_dir"."$name".".somatics.tsv";
    my $jf = "$analysis_dir"."$project/results/$lib/$lib".".junction.list.txt";

    if (-e $mf1 && -e $jf){
      $libs2{$lib}{name} = $name;
      $libs2{$lib}{mutation_file} = $mf1;
      $libs2{$lib}{junction_file} = $jf;
    }elsif(-e $mf2 && -e $jf){
      $libs2{$lib}{name} = $name;
      $libs2{$lib}{mutation_file} = $mf2;
      $libs2{$lib}{junction_file} = $jf;
    }else{
      #print YELLOW, "\n\tDid not find at two valid files:\n\t$mf1\n\t$mf2\n\t$jf", RESET;
      next();
    }
  }
  #print Dumper %libs2;
  my $fc = keys %libs2;
  print BLUE, "\n\tFound $fc libraries with matching mutation and junction files for comparison", RESET;

  return(\%libs2);
}


######################################################################################################################################
#For all the specified projects get the non-redundant list of annotations.                                                       #
######################################################################################################################################
sub importJunctionAnnotations{
  my %args = @_;
  my @projects = @{$args{'-projects'}};
  my $analysis_dir = $args{'-analysis_dir'};

  my %ja;

  print BLUE, "\n\nGetting merged junction annotations for projects: @projects", RESET;

  my $jc = 0;
  foreach my $project (@projects){
    print BLUE, "\n\tImporting $project ...", RESET;
    my $anno_file = "$analysis_dir"."$project/results/$project".".junction.list.anno.txt";
    unless (-e $anno_file){
      print RED, "\n\nCould not find project file for project: $project\n$anno_file\n\n", RESET;
      exit();
    }

    #Go through each annotation file and store values.  If multiple projects are being considered, merge the results
    #Read_count - cumulative read count
    #Score - highest observed
    #Left_Size - highest observed
    #Right_Size - highest observed

    open (ANNO, "$anno_file") || die "\n\nCould not open annotation file: $anno_file\n\n";
    my %columns;
    my $header = 1;
    while(<ANNO>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header){
        my $pos = 0;
        foreach my $name (@line){
          $columns{$name}{pos}=$pos;
          $pos++;
        }
        $header = 0;
        next();
      }
      my $jid = $line[$columns{'JID'}{pos}];
      my $anchored = $line[$columns{'Anchored'}{pos}];
      my $exons_skipped = $line[$columns{'Exons_Skipped'}{pos}];

      unless ($anchored =~ /^DA$|^NDA$|^D$|^A$/){
        next();
      }

      if ($ja{$jid}){
        $ja{$jid}{read_count} += $line[$columns{'Read_Count'}{pos}];
        if ($line[$columns{'Score'}{pos}] > $ja{$jid}{score}){
          $ja{$jid}{score} = $line[$columns{'Score'}{pos}];
        }
        if ($line[$columns{'Left_Size'}{pos}] > $ja{$jid}{left_size}){
          $ja{$jid}{left_size} = $line[$columns{'Left_Size'}{pos}];
        }
        if ($line[$columns{'Right_Size'}{pos}] > $ja{$jid}{right_size}){
          $ja{$jid}{right_size} = $line[$columns{'Right_Size'}{pos}];
        }
      }else{
        $jc++;
        $ja{$jid}{read_count} = $line[$columns{'Read_Count'}{pos}];
        $ja{$jid}{score} = $line[$columns{'Score'}{pos}];
        $ja{$jid}{left_size} = $line[$columns{'Left_Size'}{pos}];
        $ja{$jid}{right_size} = $line[$columns{'Right_Size'}{pos}];
        $ja{$jid}{intron_size} = $line[$columns{'Intron_Size'}{pos}];
        $ja{$jid}{splice_site} = $line[$columns{'Splice_Site'}{pos}];
        $ja{$jid}{anchored} = $line[$columns{'Anchored'}{pos}];
        $ja{$jid}{exons_skipped} = $line[$columns{'Exons_Skipped'}{pos}];
        $ja{$jid}{donors_skipped} = $line[$columns{'Donors_Skipped'}{pos}];        
        $ja{$jid}{acceptors_skipped} = $line[$columns{'Acceptors_Skipped'}{pos}];
        $ja{$jid}{splign_validated} = $line[$columns{'Splign_Validated'}{pos}];
        $ja{$jid}{splign_align_count} = $line[$columns{'Splign_Align_Count'}{pos}];
        $ja{$jid}{gene_name} = $line[$columns{'Gene_Name'}{pos}];
      }
    }
    close(ANNO);
  }
  print BLUE, "\n\tStored $jc junctions that are anchored", RESET;
  return(\%ja);
}


######################################################################################################################################
#Import library by library mutations and exon skip junctions                                                                         #
######################################################################################################################################
sub importLibraryData{
  my %args = @_;
  my $libs_ref = $args{'-libs'};

  print BLUE, "\n\nImporting library by library junctions and mutations", RESET;

  foreach my $lib_id (keys %{$libs_ref}){
    my $lib_name = $libs_ref->{$lib_id}->{name};
    print BLUE, "\n\tProcessing: $lib_id ($lib_name)", RESET;
    my $j_file = $libs_ref->{$lib_id}->{junction_file};
    my $m_file = $libs_ref->{$lib_id}->{mutation_file};

    #Junctions first [chrX:100-200(-)]
    open (JUN, "$j_file") || die "\n\nCould not open junction file: $j_file\n\n";
    my $j_header = 1;
    while(<JUN>){
      if ($j_header){
        $j_header = 0;
        next();
      }
      my @line = split("\t", $_);
      my $jid = $line[0];
      
      #Unless this junction was stored according to the filter skip it
      unless ($junc_anno->{$jid}){
        next();
      }
      $lib_junctions{$jid}{$lib_id}=1;
    }
    close(JUN);
  
    #Mutations now [X:1000]
    open (MUT, "$m_file") || die "\n\nCould not open mutation file: $m_file\n\n";
    while(<MUT>){
      chomp($_);
      my @line = split("\t", $_);
      my $mid = $line[0];
      #Strip of the 'chr' if present
      if ($mid =~ /^chr(.*)/){
        $mid = $1;
      }
      $lib_mutations{$mid}{$lib_id}=1;
    }
    close(MUT);
  }

  my $jcount = keys %lib_junctions;
  my $mcount = keys %lib_mutations;
  my $lcount = keys %{$libs_ref};
  print BLUE, "\nStored $jcount distinct junctions and $mcount distinct mutations for $lcount libraries", RESET;
  
  return();
}


######################################################################################################################################
#Determine the known junction associated with each mutation and use this info to infer the strand                                    #
#How to deal with splice site mutations that correspond to a known canonical junction plus an exon skipping junction?
#Select the one with the smallest intron size?
######################################################################################################################################
sub getMutatedJunctions{
  my %args = @_;
  my $infile = $args{'-infile'};
  my $mutations_ref = $args{'-mutations'};

  print BLUE, "\n\nAssociating each mutation with known junctions (inferring donor/acceptor and strand)", RESET;

  my %mut;

  
  #First find all the overlaps between known junctions and mutation positions
  my $m_temp = "temp_mutations.bed";
  my $j_temp = "temp_ref_junctions.bed";
  my $o_temp = "temp_ref_overlaps.bed";
  open (MBED, ">$m_temp") || die "\n\nCould not open temp mutations BED file: $m_temp\n\n";
  open (JBED, ">$j_temp") || die "\n\nCould not open temp junctions BED file: $j_temp\n\n";

  #Print temp BED file of known junctions
  my %ref_junctions;
  open (REF_JUNC, "$infile") || die "\n\nCould not open reference junctions file\n\n";
  while(<REF_JUNC>){
    chomp($_);
    my @line = split("\t", $_);
    my $jid = $line[0];
    my $chr = $line[1];
    my $left = $line[2];
    my $right = $line[3];
    my $strand = $line[4];
    my $t_count = $line[5];
    my $gene_name = $line[6];
    my $trans_list = $line[7];
    $ref_junctions{$jid}{chr} = $chr;
    $ref_junctions{$jid}{left} = $left;
    $ref_junctions{$jid}{right} = $right;
    $ref_junctions{$jid}{strand} = $strand;
    $ref_junctions{$jid}{t_count} = $t_count;
    $ref_junctions{$jid}{gene_name} = $gene_name;
    $ref_junctions{$jid}{trans_list} = $trans_list;

    print JBED "$chr\t$left\t$right\t$jid\t.\t$strand\n";
  }
  close(REF_JUNC);

  #Print temp BED file of observed mutations (X:1000)
  foreach my $mut (keys %{$mutations_ref}){
    if ($mut =~ /(\w+)\:(\d+)/){
      my $chr = $1;
      my $left = $2-1;
      my $right = $2+1;

      #Fake strand for now
      print MBED "chr$chr\t$left\t$right\t$mut\t.\t+\n";
    }else{
      print RED, "\n\nCould not understand mutation record format: $mut\n\n", RESET;
      exit();
    }
  }

  close(MBED);
  close(JBED);

  #Determine overlap of mutations and all known junctions - some overlaps may be missed because my known junctions do not include ACEview as used by Rodrigo
  #Use bedtools
  my $bed_cmd1 = "$bedtools_bin_dir"."intersectBed -a $j_temp -b $m_temp -wa -wb > $o_temp";
  print BLUE, "\n\t$bed_cmd1", RESET;
  system($bed_cmd1);

  #No go through each overlapping junction/mutation pair.  Find ones that correspond to an actual splice site mutation.
  #For these, note the junction, the strand and whether the mutation is at the donor or acceptor splice site
  #The junction ID will be contained in column 4, the mutation ID in column 10
  open (OVERLAP, "$o_temp") || die "\n\nCould not open overlaps file: $o_temp\n\n";
  while(<OVERLAP>){
    chomp($_);
    my @line = split("\t", $_);
    my $jid = $line[3];
    my $strand = $line[5];
    my $mid = $line[9];

    #Get the mutation position
    my $chr;
    my $pos;
    if ($mid =~ /(\w+)\:(\d+)/){
      $chr = $1;
      $pos = $2;
    }else{
      print RED, "\n\nCould not understand mutation record format: $mid\n\n", RESET;
      exit();
    }

    #Get the junction donor/acceptor splice site position coords
    my ($left, $left1, $left2, $right, $right1, $right2, $intron_size);
    if ($jid =~ /chr\w+\:(\d+)\-(\d+)\(\S\)/){
      $left = $1;
      $right = $2;
      $left1 = $left+1;
      $left2 = $left+2;
      $right1 = $right-2;
      $right2 = $right-1;
      $intron_size = ($right-$left)+1;
    }else{
      print RED, "\n\nCould not interpret junction ID: $jid\n\n", RESET;
      exit();
    }

    #Does the mutation position actually correspond to the donor or acceptor of this junction?
    unless ($pos == $left1 || $pos == $left2 || $pos == $right1 || $pos == $right2){
      #print YELLOW, "\n\nMutation does not correspond to a splice site of the overlaping junction:\nmid:$mid\tjid:$jid\tpos:$pos\tleft1:$left1\tleft2:$left2\tright1:$right1\tright2:$right2", RESET;
      next();
    }
    
    my $da_status = '';
    if ($strand eq "+"){
      if ($pos == $left1 || $pos == $left2){
        $da_status = "Donor";
      }else{
        $da_status = "Acceptor";
      }
    }else{
      if ($pos == $left1 || $pos == $left2){
        $da_status = "Acceptor";
      }else{
        $da_status = "Donor";
      }
    }

    #Get a list of libraries containing this mutation and attach this list to the anchored mutations object
    my %lib_ids = %{$mutations_ref->{$mid}};
    my %lib_names;
    foreach my $lib_id (keys %lib_ids){
      my $lib_name = $libs_ref->{$lib_id}->{name};
      $lib_names{$lib_name}=1;
    }
    my @tmp1 = sort keys %lib_ids;
    my $mutation_lib_count = scalar(@tmp1);
    my $mutation_lib_list = join(",", @tmp1);
    
    my @tmp2 = sort keys %lib_names;
    my $mutation_lib_name_list = join(",", @tmp2);
    
    my @genes = split(",", $ref_junctions{$jid}{gene_name});

    #print BLUE, "\n\nMutation corresponds to a splice site of the overlaping junction:\nmid:$mid\tjid:$jid\tpos:$pos\tleft1:$left1\tleft2:$left2\tright1:$right1\tright2:$right2", RESET;
    if ($mut{$mid}){
      my $j_ref = $mut{$mid}{known_junctions};
      $j_ref->{$jid}->{strand} = $strand;
      $j_ref->{$jid}->{donor_acceptor} = $da_status;
      $j_ref->{$jid}->{left} = $left;
      $j_ref->{$jid}->{right} = $right;
      $j_ref->{$jid}->{intron_size} = $intron_size;
      $j_ref->{$jid}->{gene_name} = $ref_junctions{$jid}{gene_name};
      
      #Store a grand collection of gene IDS
      my $genes_ref = $mut{$mid}{gene_names};
      foreach my $gene (@genes){
        unless ($gene eq "n/a"){
          $genes_ref->{$gene}=1;
        }
      }
    }else{
      my %tmp;
      my %tmp2;
      $tmp{$jid}{strand} = $strand;
      $tmp{$jid}{donor_acceptor} = $da_status;
      $tmp{$jid}{left} = $left;
      $tmp{$jid}{right} = $right;
      $tmp{$jid}{intron_size} = $intron_size;
      $tmp{$jid}{gene_name} = $ref_junctions{$jid}{gene_name};
      $mut{$mid}{donor_acceptor} = $da_status;
      $mut{$mid}{known_junctions} = \%tmp;
      $mut{$mid}{alternate_junctions} = \%tmp2;
      $mut{$mid}{strand} = $strand;
      $mut{$mid}{chromosome} = $chr;
      $mut{$mid}{mutation_position} = $pos;
      $mut{$mid}{mutation_lib_list} = $mutation_lib_list;
      $mut{$mid}{mutation_lib_name_list} = $mutation_lib_name_list;
      $mut{$mid}{mutation_lib_count} = $mutation_lib_count;

      #Store a grand collection of gene IDS
      my %genes;
      foreach my $gene (@genes){
        unless ($gene eq "n/a"){
          $genes{$gene}=1;
        }
      }
      $mut{$mid}{gene_names} = \%genes;
    }
  }
  close (OVERLAP);

  #Note that some mutations correspond to multiple junctions that share the same donor/acceptor site (exon skips of different boundaries on the other end)
  #If they all share the same acceptor/donor site, they should all be affected by a mutation... (i.e. reduce usage and possible increased usage of some other donor/acceptor site nearby)


  #Clean up the temp files
  my $rm_cmd = "rm -f $m_temp $j_temp $o_temp";
  print BLUE, "\n\t$rm_cmd", RESET;
  system($rm_cmd);

  return(\%mut);
}


######################################################################################################################################
#Determine the overlap between all anchored splice site mutations and all stored expressed junctions                                 #
######################################################################################################################################
sub getAlternateJunctions{
  my %args = @_;
  my $junctions_ref = $args{'-junctions'};

  #Determine the overlap between all stored junction events and splice site mutations
  #Only consider the mutations that were successfully anchored in the previous step...
  #The junction must be anchored to at least one of the known junctions (on the same strand) associated with the mutation identified in the previous step
  #Possible outcomes for expressed/observed junctions compared to mutation positions:
  #1.) No overlap (flanking or proximal) - simply ignore these junctions
  #2.) The overlaping junctions found may be one of the known junctions already found, in this case, note the libraries and read counts
  #3.) A putative alternate junction on the same strand and within k bases of the mutation or flanking it completely... Store these as alternative junctions associated with the mutation.  Note the libraries and read counts

  print BLUE, "\n\nDetermining overlap between all observed junctions and all mutation positions", RESET;
  #my %lib_mutations;
  #my %lib_junctions;

  my $m_temp = "temp_mutations.bed";
  my $j_temp = "temp_obs_junctions.bed";
  my $o_temp = "temp_obs_overlaps.bed";
  open (MBED, ">$m_temp") || die "\n\nCould not open temp mutations BED file: $m_temp\n\n";
  open (JBED, ">$j_temp") || die "\n\nCould not open temp junctions BED file: $j_temp\n\n";

  #Temp BED file for anchored mutations
  foreach my $mid (keys %{$anchored_mutations}){
    if ($mid =~ /(\w+)\:(\d+)/){
      my $chr = $1;
      my $left = $2-$max_distance;
      my $right = $2+$max_distance;
      my $strand = $anchored_mutations->{$mid}->{strand};
      print MBED "chr$chr\t$left\t$right\t$mid\t.\t$strand\n";
    }else{
      print RED, "\n\nCould not understand mutation record format: $mid\n\n", RESET;
      exit();
    }
  }
  #Temp BED file for all observed junctions
  foreach my $jid (keys %lib_junctions){
    if ($jid =~ /chr(\w+)\:(\d+)\-(\d+)\((\S)\)/){
      my $chr = $1;
      my $left = $2;
      my $right = $3;
      my $strand = $4;
      print JBED "chr$chr\t$left\t$right\t$jid\t.\t$strand\n";
    }else{
      print RED, "\n\nCould not interpret junction ID: $jid\n\n", RESET;
      exit();
    }
  }
  close(MBED);
  close(JBED);

  #Determine overlap of anchored mutations and all expressed junctions (observed in at least one library)
  #Use bedtools
  my $bed_cmd1 = "$bedtools_bin_dir"."intersectBed -a $j_temp -b $m_temp -wa -wb -s > $o_temp";
  print BLUE, "\n\t$bed_cmd1", RESET;
  system($bed_cmd1);

  open (OVERLAP, "$o_temp") || die "\n\nCould not open overlaps file: $o_temp\n\n";
  while(<OVERLAP>){
    chomp($_);
    my @line = split("\t", $_);
    my $jid = $line[3];
    my $mid = $line[9];
    my $known_junctions = $anchored_mutations->{$mid}->{known_junctions};
    my $mutation_position = $anchored_mutations->{$mid}->{mutation_position};
    my $mutation_lib_list = $anchored_mutations->{$mid}->{mutation_lib_list};

    my ($chr, $left, $right, $strand, $intron_size);
    if ($jid =~ /chr(\w+)\:(\d+)\-(\d+)\((\S)\)/){
      $chr = $1;
      $left = $2;
      $right = $3;
      $strand = $4;
      $intron_size = ($right-$left)+1;      
    }else{
      print RED, "\n\nCould not interpret junction ID: $jid\n\n", RESET;
      exit();
    }

    #Three possibilities: 'known', 'flanking', 'proximal'
    my $junction_overlap_type = '';

    if ($known_junctions->{$jid}){
      #Is this junction a known junction (exact match) anchored to the mutation?
      $junction_overlap_type = "known";
    }else{
      #Is this junction a new junction, that is anchored to at least one of the known junctions?
      #It has to be anchored at the correct side.  i.e. if the donor was mutated, look for novel connections using the acceptor (and a different donor)
      my $anchored = 0;
      foreach my $kjid (keys %{$known_junctions}){
        my $kstrand = $known_junctions->{$kjid}->{strand};
        my $kleft = $known_junctions->{$kjid}->{left};
        my $kright = $known_junctions->{$kjid}->{right};
        my $da_status = $known_junctions->{$kjid}->{donor_acceptor};

        if ($kstrand eq "+" && $da_status eq "Donor" && $kright == $right){
          $anchored = 1;
        }elsif($kstrand eq "+" && $da_status eq "Acceptor" && $kleft == $left){
          $anchored = 1;          
        }elsif($kstrand eq "-" && $da_status eq "Acceptor" && $kright == $right){
          $anchored = 1;                    
        }elsif($kstrand eq "-" && $da_status eq "Donor" && $kleft == $left){
          $anchored = 1;                    
        }else{
          
        }
      }
      unless ($anchored){
        next();
      }
      #The novel junction is anchored, but how does it relate to the mutation?  flanking or proximal
      if ($left <= $mutation_position && $right >= $mutation_position){
        $junction_overlap_type = "flanking";
      }else{
        $junction_overlap_type = "proximal";
      }
    }

    #Is the alternate junction expressed in at least one of the mutation containing libraries? (an alternate junction only in libs without the mutation is a red herring...)
    my %lib_ids = %{$junctions_ref->{$jid}};
    my @tmp1 = sort keys %lib_ids;
    my $junction_libs_list = join(",", @tmp1);

    my %lib_names;
    foreach my $lib_id (keys %lib_ids){
      my $lib_name = $libs_ref->{$lib_id}->{name};
      $lib_names{$lib_name}=1;
    }
    my @tmp2 = sort keys %lib_names;
    my $junction_libs_names_list = join(",", @tmp2);

    my $matching_libs = 0;
    my $total_libs = scalar(@tmp1);
    foreach my $lib (@tmp1){
      if ($mutation_lib_list =~ /$lib/){
        $matching_libs++;
      }
    }
    #Store the alternate junction(s) - known junctions are already stored seperately
    if (($junction_overlap_type eq "flanking" || $junction_overlap_type eq "proximal") && $matching_libs > 0){
      my $lib_ratio = $total_libs/$matching_libs;
      my $alternate_junctions = $anchored_mutations->{$mid}->{alternate_junctions};
      $alternate_junctions->{$jid}->{left} = $left;
      $alternate_junctions->{$jid}->{right} = $right;
      $alternate_junctions->{$jid}->{intron_size} = $intron_size;
      $alternate_junctions->{$jid}->{junction_overlap_type} = $junction_overlap_type;
      $alternate_junctions->{$jid}->{matching_libs} = $matching_libs;
      $alternate_junctions->{$jid}->{total_libs} = $total_libs;
      $alternate_junctions->{$jid}->{junction_libs_list} = $junction_libs_list;
      $alternate_junctions->{$jid}->{junction_libs_names_list} = $junction_libs_names_list;
      $alternate_junctions->{$jid}->{lib_ratio} = $lib_ratio;
      $alternate_junctions->{$jid}->{exons_skipped} = $junc_anno->{$jid}->{exons_skipped};
      $alternate_junctions->{$jid}->{anchored} = $junc_anno->{$jid}->{anchored};
    }
  }
  close(OVERLAP);

  #Clean up the temp files
  my $rm_cmd = "rm -f $m_temp $j_temp $o_temp";
  print BLUE, "\n\t$rm_cmd", RESET;
  system($rm_cmd);

  return();
}


######################################################################################################################################
#Print output files.  Create three output files
######################################################################################################################################
sub printOutputFiles{ 
  my %args = @_;
  my $dir = $args{'-dir'};
  my $mut = $args{'-data'};

  my $outfile1 = "$outdir"."SSMutations_AND_GainedJuncs.tsv";
  my $outfile2 = "$outdir"."SSMutations_AND_GainedJuncs.Shortlist.tsv";
  my $outfile3 = "$outdir"."SSMutations_Junctions_IDs.tsv";
  my $outfile4 = "$outdir"."SSMutations_JunctionCoords.tsv";

  print BLUE, "\n\nPrinting results to output files:\n\t$outfile1\n\t$outfile2\n\t$outfile3\n\t$outfile4", RESET;

  #1.) A main results file
  #Mutation_Pos, Gene Name, Junction ID, Junction class, Number exons skipped, Correlation type (contained vs. proximal), distance, number of libs with mutation, number of libs with junction, list of libs with the mutation, list of libs with the junction, list of KNOWN junctions associated with the mutation (i.e. junction for which it is actually a splice site mutation)
  #2.) A short list that has been filtered (on the mutated vs. non-mutated ratio for example)
  open (OUT1, ">$outfile1") || die "\n\nCould not open outfile1: $outfile1\n\n";
  open (OUT2, ">$outfile2") || die "\n\nCould not open outfile2: $outfile2\n\n";
  my $header = "mcount\tmid\tgene_string\tjid\tdonor_vs_acceptor\tanchored\texons_skipped\tjunction_overlap_type\tintron_size\tmutation_lib_count\tmutation_lib_list\tmatching_libs\ttotal_libs\tlib_ratio";
  print OUT1 "$header\n";
  print OUT2 "$header\tjunction_libs_list\n";

  my $mcount = 0;
  my %tmp;
  foreach my $mid (sort keys %{$mut}){

    my $mutation_lib_name_list = $mut->{$mid}->{mutation_lib_name_list};
    my $mutation_lib_count = $mut->{$mid}->{mutation_lib_count};
    my $donor_acceptor = $mut->{$mid}->{donor_acceptor};
    my $strand = $mut->{$mid}->{strand};
    my $alternate_junctions = $mut->{$mid}->{alternate_junctions};
    my $genes_ref = $mut->{$mid}->{gene_names};
    my @gene_list = keys %{$genes_ref};
    my $gene_string = join(",", @gene_list);
    unless ($gene_string){
      $gene_string = "n/a";
    }
    $mut->{$mid}->{gene_string} = $gene_string;
    my $alt_junc_count = keys %{$alternate_junctions};
    unless ($alt_junc_count > 0){
      next();
    }
    $mcount++;
    foreach my $jid (sort {$alternate_junctions->{$a}->{intron_size} <=> $alternate_junctions->{$b}->{intron_size}} keys %{$alternate_junctions}){
      my $matching_libs = $alternate_junctions->{$jid}->{matching_libs};
      my $total_libs = $alternate_junctions->{$jid}->{total_libs};
      my $junction_libs_names_list = $alternate_junctions->{$jid}->{junction_libs_names_list};
      my $lib_ratio = $alternate_junctions->{$jid}->{lib_ratio};
      my $intron_size = $alternate_junctions->{$jid}->{intron_size};
      my $junction_overlap_type = $alternate_junctions->{$jid}->{junction_overlap_type};
      my $exons_skipped = $alternate_junctions->{$jid}->{exons_skipped};
      my $anchored = $alternate_junctions->{$jid}->{anchored};

      my $string = "$mid\t$gene_string\t$jid\t$donor_acceptor\t$anchored\t$exons_skipped\t$junction_overlap_type\t$intron_size\t$mutation_lib_count\t$mutation_lib_name_list\t$matching_libs\t$total_libs\t$lib_ratio";
      print OUT1 "$mcount\t$string\n";

      #Print out a short list limiting based on the library ratio.  i.e. if an alternative junction associated with a mutation is observed in lots of other libraries that do NOT have the mutation, these will get filtered out
      #e.g. 1 library with alt junction and mutation, but 2 libraries with the alt junction and no mutation (3/1 = lib ratio of 3).
      if ($lib_ratio <= 3){
        $tmp{$mid}=1;
        my $mcount = keys %tmp;
        print OUT2 "$mcount\t$string\t$junction_libs_names_list\n";
      }
    }
  }
  close(OUT1);
  close(OUT2);

  #3.) A file to feed into R for figure generation (mutation id, donor_acceptor, gene name, chr, strand, min_start, max_end, mutated lib ids, known junction ids, alternative junction ids)
  open (OUT3, ">$outfile3") || die "\n\nCould not open outfile3: $outfile3\n\n";
  print OUT3 "mid\tpos\tdonor_acceptor\tgene_string\tchr\tstrand\tmin_start\tmax_end\tmutation_lib_list\tmutation_lib_name_list\tknown_junction_ids_string\talternate_junction_ids_string\n";
  foreach my $mid (sort keys %{$mut}){
    my $gene_string = $mut->{$mid}->{gene_string};
    my $mutation_lib_list = $mut->{$mid}->{mutation_lib_list};
    my $mutation_lib_name_list = $mut->{$mid}->{mutation_lib_name_list};    
    my $donor_acceptor = $mut->{$mid}->{donor_acceptor};
    my $strand = $mut->{$mid}->{strand};
    my $mutation_position = $mut->{$mid}->{mutation_position};
    my $known_junctions = $mut->{$mid}->{known_junctions};
    my @known_junction_ids = keys %{$known_junctions};
    my $known_junction_ids_string = join(",", @known_junction_ids);
    unless($known_junction_ids_string){
      $known_junction_ids_string = "NA";
    }

    my $alternate_junctions = $mut->{$mid}->{alternate_junctions};
    my @alternate_junction_ids = keys %{$alternate_junctions};
    my $alternate_junction_ids_string = join(",", @alternate_junction_ids);
    unless($alternate_junction_ids_string){
      $alternate_junction_ids_string = "NA";
    }

    #Determine the outer coordinates of the junctions of interest
    my @coords;
    my $chr;
    foreach my $j (@known_junction_ids){
      if ($j =~ /(chr\w+)\:(\d+)\-(\d+)\((\S)\)/){
        $chr=$1;
        my $start=$2;
        my $end=$3;
        push(@coords, $start);
        push(@coords, $end);
      }
    }
    foreach my $j (@alternate_junction_ids){
      if ($j =~ /(chr\w+)\:(\d+)\-(\d+)\((\S)\)/){
        $chr=$1;
        my $start=$2;
        my $end=$3;
        push(@coords, $start);
        push(@coords, $end);
      }
    }
    my @coords_sort = sort {$a <=> $b} @coords;
    my $min_start = $coords_sort[0];
    my $max_end = $coords_sort[(scalar(@coords_sort))-1];

    print OUT3 "$mid\t$mutation_position\t$donor_acceptor\t$gene_string\t$chr\t$strand\t$min_start\t$max_end\t$mutation_lib_list\t$mutation_lib_name_list\t$known_junction_ids_string\t$alternate_junction_ids_string\n";
  }
  close(OUT3);

  #Print out a file with junction coordinates for convenience in building gene model diagrams in R
  my %junc_list;
  foreach my $mid (sort keys %{$mut}){
    my $known_junctions = $mut->{$mid}->{known_junctions};
    my @known_junction_ids = keys %{$known_junctions};
    my $alternate_junctions = $mut->{$mid}->{alternate_junctions};
    my @alternate_junction_ids = keys %{$alternate_junctions};
    my @juncs = (@known_junction_ids,@alternate_junction_ids);
    foreach my $j (@juncs){
      if ($j =~ /(chr\w+)\:(\d+)\-(\d+)\((\S)\)/){
        my $chr=$1;
        my $start=$2;
        my $end=$3;
        my $strand=$4;
        $junc_list{$j}{chr}=$chr;
        $junc_list{$j}{strand}=$strand;
        $junc_list{$j}{start}=$start;
        $junc_list{$j}{end}=$end;
        if ($alternate_junctions->{$j}){
          $junc_list{$j}{type}="Alternate";
        }else{
          $junc_list{$j}{type}="Canonical";
        }
      }
    }
  }
  open (OUT4, ">$outfile4") || die "\n\nCould not open outfile4: $outfile4\n\n";
  print OUT4 "jid\tchr\tstrand\tstart\tend\ttype\n";
  foreach my $j (sort keys %junc_list){
    print OUT4 "$j\t$junc_list{$j}{chr}\t$junc_list{$j}{strand}\t$junc_list{$j}{start}\t$junc_list{$j}{end}\t$junc_list{$j}{type}\n";
  }
  close(OUT4);

  return();
}
