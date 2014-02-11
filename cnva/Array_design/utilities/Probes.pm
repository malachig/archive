=head1 NAME

Probes.pm - Basic methods for scoring probes or choosing optimal probes

=head1 SYNOPSIS

use Probes qw(:all);

=head2 NOTE

Currently located in '~/Array_design/utilities'

=head2 RECENT CHANGES

None.  Last modified 09 January 2007

=head1 DESCRIPTION

Generic utility for dealing with probe files and probe sequences

=head1 EXAMPLES

use lib './';

use utilities::Probes qw(:all);

=head1 SEE ALSO

None

=head1 BUGS

Contact author via email

=head1 AUTHOR

Written by Malachi Griffith (malachig@bcgsc.ca)

=head1 ACKNOWLEDGEMENTS

University of British Columbia Graduate Studies

Michael Smith Foundation for Health Research

Natural Sciences and Engineering Research Council

Genome British Columbia

=head1 AFFLIATIONS

Malachi Griffith is supervised by Marco A. Marra

Genome Sciences Centre, BC Cancer Research Centre, BC Cancer Agency, UBC Faculty of Medicine - Medical Genetics

=head1 SUBROUTINES

=cut

package utilities::Probes;
require Exporter;

@ISA = qw(Exporter);
@EXPORT = qw();

@EXPORT_OK = qw(&tmCalc &tmConverter &pairFoldCalc &getCurrentProbeIds);

%EXPORT_TAGS = (all => [qw(&tmCalc &tmConverter &pairFoldCalc &getCurrentProbeIds)]);

use strict;
use Data::Dumper;

=head2 tmCalc

=over 3

=item Function:

Calculate the Tm of a probe using the Nearest Neighbour methods of Breslauer et.al, 1986, 

and using the values from Sugimoto et. al, 1996

Basically reproduces the calculation at http://www.basic.nwu.edu/biotools/oligocalc.html#helpthermo

=item Return:

Tm value of probe in degrees Kelvin

=item Args:

'-sequence' => probe sequence string
'-silent' => 1 or 0

=item Example(s):

my $temp_k = &tmCalc('-sequence'=>$probe_seq, '-silent' => 0);

=back

=cut

#######################
#tmCalc               #
#######################
sub tmCalc{
  my %args = @_;

  my $seq = $args{'-sequence'};
  my $silent = $args{'-silent'};

  #Create a hash to store the deltaH and deltaS values for every possible 2 bp combination (AA/TT, AT/TA, etc.)
  my %thermo;

  #DeltaH values are calories/mol - Taken from Sugimoto, et al., NAR, 1996 - As done in Pan et al., Molecular Cell 2004.
  $thermo{deltaH}{aa} = -8000;
  $thermo{deltaH}{tt} = -8000;
  $thermo{deltaH}{at} = -5600;
  $thermo{deltaH}{ta} = -6600;
  $thermo{deltaH}{ca} = -8200;  $thermo{deltaH}{tg} = -8200;
  $thermo{deltaH}{ct} = -6600;  $thermo{deltaH}{ag} = -6600;
  $thermo{deltaH}{ga} = -8800;  $thermo{deltaH}{tc} = -8800;
  $thermo{deltaH}{gt} = -9400;  $thermo{deltaH}{ac} = -9400;
  $thermo{deltaH}{cg} = -11800;
  $thermo{deltaH}{gc} = -10500;
  $thermo{deltaH}{gg} = -10900;
  $thermo{deltaH}{cc} = -10900;

  #DeltaS values are calories/mol/K
  $thermo{deltaS}{aa} = -21.9;
  $thermo{deltaS}{tt} = -21.9;
  $thermo{deltaS}{at} = -15.2;
  $thermo{deltaS}{ta} = -18.4;
  $thermo{deltaS}{ca} = -21.0;  $thermo{deltaS}{tg} = -21.0;
  $thermo{deltaS}{ct} = -16.4;  $thermo{deltaS}{ag} = -16.4;
  $thermo{deltaS}{ga} = -23.5;  $thermo{deltaS}{tc} = -23.5;
  $thermo{deltaS}{gt} = -25.5;  $thermo{deltaS}{ac} = -25.5;
  $thermo{deltaS}{cg} = -29.0;
  $thermo{deltaS}{gc} = -26.4;
  $thermo{deltaS}{gg} = -28.4;
  $thermo{deltaS}{cc} = -28.4;

  #Initialize counters for each of the possible adjacent nucleotide pairs
  my %neighbour_counts;
  $neighbour_counts{aa} = 0;
  $neighbour_counts{tt} = 0;
  $neighbour_counts{at} = 0;
  $neighbour_counts{ta} = 0;
  $neighbour_counts{ca} = 0;  $neighbour_counts{tg} = 0;
  $neighbour_counts{ct} = 0;  $neighbour_counts{ag} = 0;
  $neighbour_counts{ga} = 0;  $neighbour_counts{tc} = 0;
  $neighbour_counts{gt} = 0;  $neighbour_counts{ac} = 0;
  $neighbour_counts{cg} = 0;
  $neighbour_counts{gc} = 0;
  $neighbour_counts{gg} = 0;
  $neighbour_counts{cc} = 0;

  my $seq_length = length($seq);

  unless($seq_length > 1){
    print "\nProbe sequence has length of 1 or less\n\n";
    exit();
  }

  #Go through the probe sequence and get all the nucleotide neighbours 
  for (my $i = 0; $i <= $seq_length-2; $i++){

    my $basepair1 = substr($seq, $i, 1);
    my $basepair2 = substr($seq, $i+1, 1);
    my $couple = "$basepair1"."$basepair2";
    my $lower_couple = lc($couple);

    #Count the neighbors by incrementing the appropriate counter
    $neighbour_counts{$lower_couple}++;
  }

  my $h_nuc = 5000; #Calories per mol
  $h_nuc = 3400; #According to Sugimoto - at least thats what this website says 

  #Reproduce the calculation at http://www.basic.nwu.edu/biotools/oligocalc.html#helpthermo
  #NOTE:  ** This calculation will give you the Nearest Neighbor Value shown on this website.**

  #The nearest neighbor and thermodynamic calculations are done essentially as described by Breslauer et al., (1986) Proc. Nat. Acad. Sci. 83:3746-50
  #but using the values published by Sugimoto et al., (1996) Nucl. Acids Res. 24:4501-4505 (Abstract). 
  #assumes that the sequences are not symmetric and contain at least one G or C. The minimum length for the query sequence is 8.

  #To find the deltaH multiply the deltaH value for each couple by the number of occurences of each pair, add all these up and add the initiation value 'h_nuc'
  #in Calories/mol
  my $deltaH = - (($neighbour_counts{aa} * $thermo{deltaH}{aa}) + ($neighbour_counts{tt} * $thermo{deltaH}{tt}) +
		  ($neighbour_counts{at} * $thermo{deltaH}{at}) + ($neighbour_counts{ta} * $thermo{deltaH}{ta}) +
		  ($neighbour_counts{ca} * $thermo{deltaH}{ca}) + ($neighbour_counts{tg} * $thermo{deltaH}{tg}) +
		  ($neighbour_counts{ct} * $thermo{deltaH}{ct}) + ($neighbour_counts{ag} * $thermo{deltaH}{ag}) +
		  ($neighbour_counts{ga} * $thermo{deltaH}{ga}) + ($neighbour_counts{tc} * $thermo{deltaH}{tc}) +
		  ($neighbour_counts{gt} * $thermo{deltaH}{gt}) + ($neighbour_counts{ac} * $thermo{deltaH}{ac}) +
		  ($neighbour_counts{cg} * $thermo{deltaH}{cg}) + ($neighbour_counts{gc} * $thermo{deltaH}{gc}) +
		  ($neighbour_counts{gg} * $thermo{deltaH}{gg}) + ($neighbour_counts{cc} * $thermo{deltaH}{cc}));



  #in calories/mol/K
  my $deltaS = - (($neighbour_counts{aa} * $thermo{deltaS}{aa}) + ($neighbour_counts{tt} * $thermo{deltaS}{tt}) +
		  ($neighbour_counts{at} * $thermo{deltaS}{at}) + ($neighbour_counts{ta} * $thermo{deltaS}{ta}) +
		  ($neighbour_counts{ca} * $thermo{deltaS}{ca}) + ($neighbour_counts{tg} * $thermo{deltaS}{tg}) +
		  ($neighbour_counts{ct} * $thermo{deltaS}{ct}) + ($neighbour_counts{ag} * $thermo{deltaS}{ag}) +
		  ($neighbour_counts{ga} * $thermo{deltaS}{ga}) + ($neighbour_counts{tc} * $thermo{deltaS}{tc}) +
		  ($neighbour_counts{gt} * $thermo{deltaS}{gt}) + ($neighbour_counts{ac} * $thermo{deltaS}{ac}) +
		  ($neighbour_counts{cg} * $thermo{deltaS}{cg}) + ($neighbour_counts{gc} * $thermo{deltaS}{gc}) +
		  ($neighbour_counts{gg} * $thermo{deltaS}{gg}) + ($neighbour_counts{cc} * $thermo{deltaS}{cc}));

  #This formula corresponds to that found at the website:
  #http://www.basic.nwu.edu/biotools/oligocalc.html#helpthermo

  #Assume a salt concentration of 50 mM (milliMolar)
  my $na_salt_conc = (50*1e-3);

  #Assume a primer concentration of 50 nM (nanoMolar)
  my $primer_conc = (50*1e-9);

  #The Gas constant R = 1.987 cal/mol * K
  #Assumes annealing occurs at pH of 7.0
  #Tm assumes the sequences are not symmetric and contain at least one G or C
  #The oligo should be at least 8 nucleotides long for accurate calculations
  #Calculations use 50nM primer conc, and 50 mM salt conc.

  my $tm_k = ($deltaH - $h_nuc)/($deltaS + (1.987 * log(1/$primer_conc))) + (16.6*(log($na_salt_conc)/log(10)));

  unless ($silent == 1){
    print "\nDeltaH = $deltaH\tDeltaS = $deltaS\tTm_K = $tm_k";
  }

  return($tm_k);

}


=head2 tmConverter

=over 3

=item Function:

Convert between celcius and kelvin melting temperature values

=item Return:

Tm as float

=item Args:

'-tm' => $probe_Tm
'-scale' => 'Celcius' or 'Kelvin'

=item Example(s):

my $Tm_convert = tmConverter ('-tm'=>$probe_Tm, '-scale'=>'Kelvin');

=back

=cut


#######################
#tmConverter          #
#######################
sub tmConverter{

  my %args = @_;
  my $tm = $args{'-tm'};
  my $scale = $args{'-scale'};

  my $tm_convert;

  unless ($tm){
    print "\nRequired parameter for tmConverter missing\n\n";
    exit();
  }

  unless ($scale eq "Celcius" || $scale eq "Kelvin"){
    print "\nTemperature scale specification not understood by tmConverter\n\n";
    exit();
  }

  #If the temperature provide was in degrees Kelvin:
  if ($scale eq "Kelvin"){
    $tm_convert = $tm-273.15;
  }

  #If the temperature provided was in degrees celcius:
  if ($scale eq "Celcius"){
    $tm_convert = 273.15+$tm;
  }
  return($tm_convert)
}


=head2 pairFoldCalc

=over 3

=item Function:

Uses Simfold/PairFold from RNAsoft to calculate the minimum energy of the fold of an RNA sequence to itself (within probe folding)

or alternately it can find the minimum energy of two copies of the same interacting with each other

The smaller the value, the more stable the fold (i.e. probes with large -ves are bad!)

Simfold considers folding within a sequence.  If only one sequence is provided and PairFold is specified a self-self comparison is conducted

If a second sequence is provided and PairFold is specified, they are folded against each other

If only a single sequence is provided to PairFold will not execute

See /home/user/pairFold/MultiRNAFold-1.1 for details of this software

=item Return:

Variable containing requested score

=item Args:

'-sequence1' => $sequence    #String variable of probe sequence
'-sequence2' => $sequence2   #String variable of comparison sequence 
'-program' => $program       #'simfold' or 'pairfold'
'-silent' => 1 or 0          #Whether values will be printed out or just silently returned
'-bin_dir' => $bin_dir       #Specify the location of the pairfold/simfold bins - allows you to copy these to /tmp for cluster jobs
                             #to reduce the number of system calls over the network/filer

=item Example(s):

my $simfold_score = &pairFoldCalc('-sequence1'=>$sequence, '-program'=>$program, '-silent'=>0, '-bin_dir'=>$bin_dir);

=back

=cut


#######################
#pairFoldCalc         #
#######################
sub pairFoldCalc{

  my %args = @_;
  my $sequence1 = $args{'-sequence1'};
  my $sequence2 = $args{'-sequence2'};
  my $program = $args{'-program'};
  my $silent = $args{'-silent'};
  my $bin_dir = $args{'-bin_dir'};

  my $score;

  unless ($sequence1 && $program){
    print "\nRequired parameter for pairFoldCalc() missing\n\n";
    exit();
  }
  unless ($silent == 0 || $silent == 1){
    print "\nRequired parameter '-silent' for pairFoldCalc() missing\n\n";
    exit();
  }
  unless ($bin_dir){
    print "\nRequired parameter '-bin_dir' for pairFoldCalc() missing\n\n";
    exit();
  }

  #Confirm format of program specified
  unless ($program eq "simfold" || $program eq "pairfold"){
    print "\nSpecified program: $program for pairFoldCalc does not seem valid\n\n";
    exit();
  }

  #Define the program command
  my $cmd;
  if ($program eq "simfold"){
    $cmd = "$bin_dir"."$program"." $sequence1";
  }elsif ($program eq "pairfold" && $sequence2){
    $cmd = "$bin_dir"."$program"." $sequence1 $sequence2";
  }else{
    $cmd = "$bin_dir"."$program"." $sequence1 $sequence1";
  }

  my $result = `$cmd`;

  #print "\nDEBUG\n$result\nDEBUG\n";

  my @lines = split ("\n", $result);

  foreach my $line (@lines){

    if ($line =~ /\s+(\d+\.\d+)/){ #Positive Number
      $score = $1;
    }elsif ($line =~ /\s+(\-\d+\.\d+)/){ #Negative Number
      $score = $1;
    }elsif ($line =~ /\s+(\d+)/){ #Zero or simple integer
      $score = $1;
    }
  }

  unless ($score){
    print "\nLine format not understood!\n\n";
    print "\n$result";
    exit();
  }

  unless ($silent == 1){
    print "\n$program Score = $score\n";
  }

  return($score)
}


=head2 getCurrentProbeIds

=over 3

=item Function:

Examine a user specified directory containing probe IDs

Determine the current maximum Probe_ID and Probeset_ID

If the specified directory is empty, start with 1

In either case prompt the user to confirm the choice

If the 'force' option is specified, the max IDs are used without asking the user

=item Return:

Array (current_max_probe_ID, current_max_probeset_ID)

=item Args:

'-probe_dir' => $probe_dir    #Full path to directory containing probe files

=item Example(s):

my @current_ids = &getCurrentProbeIds('-probe_dir'=>$probe_dir);

=back

=cut


#######################
#getCurrentProbeIds   #
#######################
sub getCurrentProbeIds{
  my %args = @_;
  my $probe_dir = $args{'-probe_dir'};
  my $force = $args{'-force'};

  unless ($probe_dir){
    print "\nRequired parameter for getCurrentProbeIds() missing\n\n";
    exit();
  }

  unless ($probe_dir =~ /.*\/$/){
    $probe_dir = "$probe_dir"."/";
  }

  #Make sure the specified directory is valid
  unless (-e $probe_dir && -d $probe_dir){
    print "\nSpecified directory does not appear to be valid: $probe_dir\n\n";
    exit();
  }

  #Get files from this directory
  print "\nSearching $probe_dir for probe files";

  my @probe_files;
  opendir(DIRHANDLE, "$probe_dir") || die "\nCannot open directory: $probe_dir\n\n";
  my @test_files = readdir(DIRHANDLE);

  foreach my $test_file (@test_files){
    my $file_path = "$probe_dir"."$test_file";

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      print "\n\t$file_path  is a directory - skipping";
      next();
    }
    push(@probe_files, $file_path);
  }

  #If no files are found, use 0 as the current max probe and probeset ids
  my $file_count = @probe_files;
  if ($file_count == 0){

    #If the force option was used
    if ($force){
      if ($force eq "yes"){
	print "\n\nDirectory appears empty and force option was selected. Will start counting probe and probeset IDs from 1!";
	my @ids = (0,0);
	return(@ids);
      }
    }

    print "\n\nDirectory appears empty. Start counting probe and probeset IDs from 1 (y/n)? ";

    my $answer = <>;
    chomp($answer);
    if ($answer eq "y" || $answer eq "Y"){
      #Start counting IDs from scratch
      my @ids = (0,0);
      return(@ids);

    }else{
      #Despite the lack of existing probe files, the user wishes to specify the IDs to use
      print "\nDesired starting Probe_id? : ";
      my $probe_id = <>;
      chomp($probe_id);

      unless($probe_id =~ /^\d+/){
	print "\nNot a valid Probe ID!  Use a simple integer!\n\n";
	exit();
      }

      print "\nDesired starting ProbeSet_id? : ";
      my $probeset_id = <>;
      chomp($probeset_id);

      unless($probeset_id =~ /^\d+/){
	print "\nNot a valid ProbeSet ID!  Use a simple integer!\n\n";
	exit();
      }
      my @ids;
      push (@ids, $probe_id-1);
      push (@ids, $probeset_id-1);
      return(@ids);
    }
  }

  #Some files were found
  print "\nFound files:";
  foreach my $file (@probe_files){
    chomp($file);
    print "\n\t$file";
  }
  print "\n\n";

  #Go through each file, check if it is a valid probe file, if so determine the max probe_id and probeset_id
  my $grand_max_probe_id = 0;
  my $grand_max_probeset_id = 0;

  my %file_summary;

 FILE:foreach my $file (@probe_files){
    chomp($file);

    open (INFILE, "$file") || die "\nCould not open input file: $file\n\n";

    my $first_line = 1;
    my $first_data = 1;
    my $max_probe_id = 0;
    my $min_probe_id = 0;
    my $max_probeset_id = 0;
    my $min_probeset_id = 0;

    while (<INFILE>){
      chomp($_);

      #Process the header line and determine the column which contains the probe_tm values
      if ($first_line == 1){
        my @header = split("\t", $_);

	#If this file is a valid probe file it should have 'Probe_Count' and 'ProbeSet_ID' as the first two column
	unless ($header[0] eq "Probe_Count" && $header[1] eq "ProbeSet_ID"){
	  print "\nFile: $file does not appear to be a valid probe_file (unexpected or missing headers) - skipping\n\n";
	  next FILE;
	}
	$first_line = 0;
	next();
      }
      my @line = split("\t", $_);

      if ($first_data == 1){
	$max_probe_id = $line[0];
	$min_probe_id = $line[0];
	$max_probeset_id = $line[1];
	$min_probeset_id = $line[1];

	$first_data = 0;
	next();
      }

      if ($line[0] > $max_probe_id){$max_probe_id = $line[0];}
      if ($line[0] < $min_probe_id){$min_probe_id = $line[0];}
      if ($line[1] > $max_probeset_id){$max_probeset_id = $line[1];}
      if ($line[1] < $min_probeset_id){$min_probeset_id = $line[1];}
    }

    $file_summary{$file}{min_probe_id} = $min_probe_id;
    $file_summary{$file}{max_probe_id} = $max_probe_id;
    $file_summary{$file}{min_probeset_id} = $min_probeset_id;
    $file_summary{$file}{max_probeset_id} = $max_probeset_id;

    if ($max_probe_id > $grand_max_probe_id){$grand_max_probe_id = $max_probe_id;}
    if ($max_probeset_id > $grand_max_probeset_id){$grand_max_probeset_id = $max_probeset_id;}
  }

  foreach my $file (sort {$file_summary{$a}->{min_probe_id} <=> $file_summary{$b}->{min_probe_id}} keys %file_summary){
    print "\nFile: $file\n\tProbe_ids: ($file_summary{$file}{min_probe_id} - $file_summary{$file}{max_probe_id})\n\tProbeset_ids: ($file_summary{$file}{min_probeset_id} - $file_summary{$file}{max_probeset_id})\n\n";
  }

  #If the force option was used
  if ($force){
    if ($force eq "yes"){
      print "\n\nWill start counting probe and probeset IDs from current maxes! ($grand_max_probe_id and $grand_max_probeset_id)";
      my @ids = ($grand_max_probe_id,$grand_max_probeset_id);
      return(@ids);
    }
  }

  print "\nCurrent max probe_id and probesets ids are: $grand_max_probe_id and $grand_max_probeset_id\n\n";
  print "\nContinue counting from these values (y/n)? ";

  my $answer = <>;
  chomp($answer);
  if ($answer eq "y" || $answer eq "Y"){
    #Use the current max IDs
    my @ids = ($grand_max_probe_id,$grand_max_probeset_id);
    return(@ids);

  }else{
    #Despite existing files with probes, the user wishes to specify the IDs to use
    print "\nDesired starting Probe_id? : ";
    my $probe_id = <>;
    chomp($probe_id);

    unless($probe_id =~ /^\d+/){
      print "\nNot a valid Probe ID!  Use a simple integer!\n\n";
      exit();
    }

    print "\nDesired starting ProbeSet_id? : ";
    my $probeset_id = <>;
    chomp($probeset_id);

    unless($probeset_id =~ /^\d+/){
      print "\nNot a valid ProbeSet ID!  Use a simple integer!\n\n";
      exit();
    }
    my @ids;
    push (@ids, $probe_id-1);
    push (@ids, $probeset_id-1);
    return(@ids);
  }
  return();
}


1;


