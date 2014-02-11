#!/usr/bin/perl5.8.5
#!/usr/local/bin/perl58 -w
#Written by Malachi Griffith

#This script takes a microarray data file as input and calculates AS splicing index values by a variety of methods
#The input file should contain, normalized, background corrected intensity values
#The user will indicate the columns of intensity values for each of the two conditions to be compared
#The results will be appended to the input file and dumped to a new results file

#For many of the calculations made below, the position and type of the probe must be considered.  This will require retrieval of data from
#the corresponding ALEXA database for the Array design used to generate the input data.

#Note: 'In' = Inclusion, 'Ex' = Exclusion

use strict;
use Data::Dumper;
use Getopt::Long;
use Math::Complex;

use lib '/home/malachig/AlternativeSplicing/perl_bin';
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);
use Statistics::Descriptive;

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $probe_file = '';
my $out_file = '';
my $column_i1 = '';  #Sample 1 column of intensity values (eg. LnCAP, BGC2 values)
my $sample1 = '';    #Name of the first sample
my $column_i2 = '';  #Sample 2 column of intensity values (eg. Brain, BGC2 values)
my $sample2 = '';    #Name of the second sample

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'probe_file=s'=>\$probe_file, 'out_file=s'=>\$out_file, 'column_i1=i'=>\$column_i1, 'sample1=s'=>\$sample1,
	    'column_i2=i'=>\$column_i2, 'sample2=s'=>\$sample2);

#Provide instruction to the user
print "\n\nUsage:";
print "\n\tSpecify the database and server to query using: --database and --server";
print "\n\tSpecify the user and password for access using: --user and --password";
print "\n\tSpecify the input probe microarray data file using: --probe_file";
print "\n\tSpecify the column contain intensity values for sample 1 (eg LnCAP, BGC2 values) using: --column_i1";
print "\n\tSpecify the name of this first sample (eg LnCAP) using: --sample1";
print "\n\tSpecify the column contain intensity values for sample 2 (eg Brain, BGC2 values) using: --column_i2";
print "\n\tSpecify the name of this second sample (eg Brain) using: --sample2";
print "\n\tSpecify the output file (resulting calculations will be appended to the input file): --out_file";
print "\n\nExample: calculate_AS_spliceIndexValues.pl --database ALEXA_hs_31_35d --server jango.bcgsc.ca --user malachig --password pwd --probe_file Probes_LnBr_BGC1-2_DE_24.txt --column_i1=22 --sample1 LnCAP --column_i2=23 --sample2 Brain --out_file test.out\n\n";

#Make sure all options were specified
unless ($database && $server && $user && $password && $column_i1 && $column_i2 && $probe_file && $out_file){
  print "\nOptions missing!\n\n";
  exit();
}

#1.) Parse input probe intensity data file
my %probes;
my %genes;
my $header_line = '';

&parseDataFiles('-probe_file'=>$probe_file, '-column_1'=>$column_i1, '-column_2'=>$column_i2);
my $probes_found = keys %probes;
my $genes_found = keys %genes;
print "\n\nRetrieved intensity data for $probes_found probes, corresponding to $genes_found genes\n";

#2.) Get gene coordinates for each probe and the intron/exon structure of each gene

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = connectDB($database, $server, $user, $password);

#For each gene, go through all the probes of that gene and get genomic coordinates
&getProbeCoordinates('-dbh'=>$alexa_dbh);

#Close database connection
$alexa_dbh->disconnect();

#3.) Calculate inclusion/exclusion index values for each exon skipping probe.  These represent differential expression of an AS event
#    The values also take into account the level of gene expression (normalized to gene level expression in each tissue)
#    Inclusion - Refers to an exon (or consecutive series of exons) being included in the transcript in tissue 1 vs 2
#    Exclusion - Refers to the exon (or consecutive series of exons) being excluded in the transcript in tissue 1 vs 2
#    Both of these values will be anchored to the exon-skipping probe as a reference (will be reported on this probe line)
#    - The actual score takes into account many other probes (exon and canonical exon-exon probes within and outside the putative affected region)
&calculate_Skip_InEx_IndexValues();

#4.) Create output file
&printOutputFile('-out_file'=>$out_file);





exit();


#######################################################################################################################
#1.) Parse input probe and gene data files                                                                            #
#######################################################################################################################
sub parseDataFiles{
  my %args = @_;

  my $probe_file = $args{'-probe_file'};
  my $column1 = $args{'-column_1'} - 1; #Subtract 1 from the specified column to get the array index
  my $column2 = $args{'-column_2'} - 1; #Subtract 1 from the specified column to get the array index

  #Open input file
  open (PROBE, "$probe_file") || die "\nCould not open probe file: $probe_file\n\n";

  #Build probes object
  my $header = 1;
  my @headers;
  my %headers;
  my $probe_count = 0;
  my $line_count = 0;

  print "\n\nProcessing input probe file: $probe_file";

  while (<PROBE>){
    chomp($_);

    #Get the header positions and names
    if ($header == 1){
      $header_line = $_;
      $header = 0;
      @headers = split ("\t", $header_line);
      my $header_count = 0;

      my $i1_name = $headers[$column1];
      my $i2_name = $headers[$column2];

      print "\n\nUsing the intensity values: $i1_name versus $i2_name !!!";

      foreach my $h (@headers){
	$headers{$h}{pos} = $header_count;
	$header_count++;
      }
      #Make sure basic required columns were found
      unless ($headers{'Probe_type'} && $headers{'Exons_skipped'} && $headers{'Tm'}){
	print "\nCritical field (Probe_type or Exons_skipped or Tm) missing\n\n";
	exit();
      }
      next();
    }
    $probe_count++;
    $line_count++;

    my $line_record = $_;
    my @line = split ("\t", $_);

    my $probe_id;
    if ($line[0] =~ /(\d+)/){
      $probe_id = $1;
    }

    my $gene_id;
    my $id = $line[$headers{'AlexaGene_ID'}{pos}];
    if ($id =~ /(\d+)/){
      $gene_id = $1;
    }elsif($id eq "NA"){
      #Skip probes that do not have a gene ID such as negative control probes
      next();
    }else{
      print "\nGene ID format: $id not undestood!\n\n";
      exit();
    }

    #Populate probes hash
    $probes{$probe_id}{AlexaGene_ID} = $gene_id;
    $probes{$probe_id}{Probe_type} = $line[$headers{'Probe_type'}{pos}];
    $probes{$probe_id}{Exons_skipped} = $line[$headers{'Exons_skipped'}{pos}];
    $probes{$probe_id}{Tm} = $line[$headers{'Tm'}{pos}];
    $probes{$probe_id}{i1} = $line[$column1];
    $probes{$probe_id}{i2} = $line[$column2];

    #NOTE: *** Make sure intensity values are at least 1 (do not allow ratios 0-0.99)
    #Since background corrected values are expressed as a ratio of signal over background of Tm matched random probes, fractions do not make sense
    #Ratios less than one will cause bad things to happen in the calculations below
    if ($probes{$probe_id}{i1} < 1){
      $probes{$probe_id}{i1} = 1.0;
    }
    if ($probes{$probe_id}{i2} < 1){
      $probes{$probe_id}{i2} = 1.0;
    }

    #Store complete line info to allow printout of file later
    $probes{$probe_id}{line_count} = $line_count;
    $probes{$probe_id}{line_record} = $line_record;

  }
  close PROBE;

  #Go through all probes and associate with gene records
  print "\n\nAssociating probe and gene info\n";
  foreach my $probe_id (keys %probes){
    my $gene_id = $probes{$probe_id}{AlexaGene_ID};

    my %gene_probes;
    if ($genes{$gene_id}{probes}){
      %gene_probes = %{$genes{$gene_id}{probes}};
    }

    $gene_probes{$probe_id}{tmp} = '';

    $genes{$gene_id}{probes} = \%gene_probes;
  }

  return();
}


#######################################################################################################################
#2.) Get gene coordinates for each probe                                                                              #
#######################################################################################################################
sub getProbeCoordinates{
  my %args = @_;

  my $dbh = $args{'-dbh'};

  print "\nGetting gene coordinates for each probe\n\n";
  foreach my $gene_id (sort {$a <=> $b} keys %genes){

    print "\nProcessing Gene: $gene_id";

    #Get the probe IDs for this gene
    my %gene_probes = %{$genes{$gene_id}{probes}};

    #Get probe coordinates
    foreach my $probe_id (keys %gene_probes){
      my %probe_info = %{&getProbeInfo ('-dbh'=>$dbh, '-probe_count_id'=>$probe_id)};

      $probes{$probe_id}{unit1_start} = $probe_info{unit1_start};
      $probes{$probe_id}{unit1_end} = $probe_info{unit1_end};
      $probes{$probe_id}{unit2_start} = $probe_info{unit2_start};
      $probes{$probe_id}{unit2_end} = $probe_info{unit2_end};
    }

    #Get Gene coordinates
    my %gene_info = %{getGeneInfo ('-dbh'=>$dbh, '-gene_id'=>$gene_id)};
    $genes{$gene_id}{chr_start} = $gene_info{$gene_id}{chr_start};
    $genes{$gene_id}{chr_end} = $gene_info{$gene_id}{chr_end};
    $genes{$gene_id}{chromosome} = $gene_info{$gene_id}{chromosome};

    my $strand = $gene_info{$gene_id}{chr_strand};
    if ($strand eq "1"){
      $genes{$gene_id}{chr_strand} = "+";
    }else{
      $genes{$gene_id}{chr_strand} = "-";
    }
  }

  return();
}

#######################################################################################################################################
#3.) Calculate inclusion/exclusion index values for each exon skipping probe.  These represent differential expression of an AS event #
#######################################################################################################################################
sub calculate_Skip_InEx_IndexValues{

  print "\nCalculating Exon-skip Inclusion/Exclusion values\n\n";

  #Go through each gene and calculate inclusion/exclusion index values for eligable probes (exon-exon probes that skip 1 or more exons)
  foreach my $gene_id (sort {$a <=> $b} keys %genes){

    print "\nCalculating for gene: $gene_id";

    #Get the probes for this gene
    my %gene_probes = %{$genes{$gene_id}{probes}};

    foreach my $probe_id (sort keys %gene_probes){
      $probes{$probe_id}{exon_skip_in} = 'na';          #ratio of means method
      $probes{$probe_id}{exon_skip_ex} = 'na';          #ratio of means method
      $probes{$probe_id}{exon_skip_in_ratios} = 'na';   #mean of ratios method
      $probes{$probe_id}{exon_skip_ex_ratios} = 'na';   #mean of ratios method
      $probes{$probe_id}{exon_skip_ncount} = 'na';

      #Unless this probe is an exon-exon probe that skips at least one exon, skip it
      unless ($probes{$probe_id}{Probe_type} eq "Exon-Exon"){
	next();
      }
      unless ($probes{$probe_id}{Exons_skipped} > 0){
	next();
      }

      my $ref_unit1_end = $probes{$probe_id}{unit1_end};
      my $ref_unit2_start = $probes{$probe_id}{unit2_start};

      #1.) For each exon skipping probe, identify probe sets to use in the calculation of inclusion and exclusion index values
      #    Get the corresponding intensity values for each tissue
      my @s1_inclusion;
      my @s1_exclusion;
      my @s1_gene_expression;

      my @s2_inclusion;
      my @s2_exclusion;
      my @s2_gene_expression;

      my @inclusion_ratios;
      my @exclusion_ratios;
      my @gene_expression_ratios;

      push (@s1_exclusion, $probes{$probe_id}{i1});
      push (@s2_exclusion, $probes{$probe_id}{i2});
      push (@exclusion_ratios, ($probes{$probe_id}{i1} / $probes{$probe_id}{i2}));

      #Go through each probe for this gene and compare to the current reference probe to identify the neccessary probes for calculation
      foreach my $probe_comp_id (sort keys %gene_probes){
	my $exons_skipped = $probes{$probe_comp_id}{Exons_skipped};
	my $probe_type = $probes{$probe_comp_id}{Probe_type};

	#Skip probes that are not Exon or Canonical junction probes
	unless ($probe_type eq "Exon" || $probe_type eq "Exon-Exon"){
	  next();
	}
	if ($probe_type eq "Exon-Exon"){
	  unless ($exons_skipped == 0){
	    next();
	  }
	}

	#A.) Probes that indicate inclusion of the exon (or series of exons) - ('In', numerators)
	#    -> @s_inclusion
	#    - Exon and canonical junction probes within the affected region
	#    - Exon probes that are entirely within the unit1_end and unit2_start of the current exon-skip probe
	#    - Canonical junction probes that are entirely within the unit1_end and unit2_start of the current exon-skip probe
	#    - OR canonical junction probes that have the same unit1_end or unit2_start as the current exon-skip probe

	#exon probes
	if ($probe_type eq "Exon"){
	  if ($probes{$probe_comp_id}{unit1_start} >= $ref_unit1_end && $probes{$probe_comp_id}{unit1_end} <= $ref_unit2_start){
	    push (@s1_inclusion, $probes{$probe_comp_id}{i1});
	    push (@s2_inclusion, $probes{$probe_comp_id}{i2});
	    push (@inclusion_ratios, ($probes{$probe_comp_id}{i1} / $probes{$probe_comp_id}{i2}));
	    next();
	  }
	}
	#canonical probes
	if ($probe_type eq "Exon-Exon" && $exons_skipped == 0){
	  #entirely within region of interest
	  if ($probes{$probe_comp_id}{unit1_start} >= $ref_unit1_end && $probes{$probe_comp_id}{unit2_end} <= $ref_unit2_start){
	    push (@s1_inclusion, $probes{$probe_comp_id}{i1});
	    push (@s2_inclusion, $probes{$probe_comp_id}{i2});
	    push (@inclusion_ratios, ($probes{$probe_comp_id}{i1} / $probes{$probe_comp_id}{i2}));
	    next();
	  }
	  #have same unit1_end or unit2_start as current exon-skip probe
	  if ($probes{$probe_comp_id}{unit1_end} == $ref_unit1_end || $probes{$probe_comp_id}{unit2_start} == $ref_unit2_start){
	    push (@s1_inclusion, $probes{$probe_comp_id}{i1});
	    push (@s2_inclusion, $probes{$probe_comp_id}{i2});
	    push (@inclusion_ratios, ($probes{$probe_comp_id}{i1} / $probes{$probe_comp_id}{i2}));
	    next();
	  }
	}

	#B.) Probes to estimate the gene level expression (outside the AS region of interest) - ('In' and 'Ex' denominators)
	#    -> @s_gene_expression
	#    - Exon and canonical junction probes outside the affected region
	#    - Exon and canonical junction probes that are entirely outside of the unit1_end and unit2_start region (no overlap)

	#exon probes
	if ($probe_type eq "Exon"){
	  #entirely outside the region of interest
	  if ($probes{$probe_comp_id}{unit1_start} > $ref_unit2_start || $probes{$probe_comp_id}{unit1_end} < $ref_unit1_end){
	    push (@s1_gene_expression, $probes{$probe_comp_id}{i1});
	    push (@s2_gene_expression, $probes{$probe_comp_id}{i2});
	    push (@gene_expression_ratios, ($probes{$probe_comp_id}{i1} / $probes{$probe_comp_id}{i2}));
	    next();
	  }
	}
	#canonical probes
	if ($probe_type eq "Exon-Exon" && $exons_skipped == 0){
	  #entirely outside the region of interest
	  if ($probes{$probe_comp_id}{unit1_start} > $ref_unit2_start || $probes{$probe_comp_id}{unit2_end} < $ref_unit1_end){
	    push (@s1_gene_expression, $probes{$probe_comp_id}{i1});
	    push (@s2_gene_expression, $probes{$probe_comp_id}{i2});
	    push (@gene_expression_ratios, ($probes{$probe_comp_id}{i1} / $probes{$probe_comp_id}{i2}));
	    next();
	  }
	}

	#C.) Probe that indicates exclusion of the exon ('Ex', numerator)
	#    -> @s_exclusion
	#    - current exon skip probe being considered
	#    - this info was gathered above
      }#Comparison probe loop

      #2.) Determine whether enough probes of each category (A,B,C) were found.  If not, the 'In' and 'Ex' values will be set to 'na'
      #    - (A) For category A, we expect at least three probes for the simplest case (single exon) skipped.  Require at least 2.
      #    - (B) For category B, we expect at least two probes to be outside the region of interest for the simplest case (single exon skip
      #          in a gene with only three exons).  Require at least 2.
      #    - (C) Only one probe is expected and required for category C.

      #Since the @s1 and @s2 arrays contain values for the same probes, only need to count 1.
      my $inclusion_count = @s1_inclusion;
      my $gene_expression_count = @s1_gene_expression;
      my $exclusion_count = @s1_exclusion;

      #Make note of the number of probe data points that will be used for the following calculation
      $probes{$probe_id}{exon_skip_ncount} = $inclusion_count + $gene_expression_count + $exclusion_count;

      #In insufficient values were found for the calculation, skip this exon-skip probe
      unless ($inclusion_count >= 2 && $gene_expression_count >= 2 && $exclusion_count == 1){
	next();
      }

      #3.) Actually calculate the In/Ex values for the current reference probe and associate this value with the probe record
      #    Note: The following calculation finds the mean for each category, calculates a gene expression normalized value, and then
      #          calculates a ratio between the two samples (ratio of means)
      #    Could also calculate ratios between the two samples for every probe, calculate the mean of these ratios for each category
      #    and then calculate a gene expression normalized value... (mean of ratios). Li et al. 2006 uses this approach
      #    - not sure which of these makes more sense.

      #3-A.) ratio of means
      my $s1_in_stat = Statistics::Descriptive::Full->new();
      $s1_in_stat->add_data(@s1_inclusion);
      my $s1_in_mean = $s1_in_stat->mean();

      my $s1_ge_stat = Statistics::Descriptive::Full->new();
      $s1_ge_stat->add_data(@s1_gene_expression);
      my $s1_ge_mean = $s1_ge_stat->mean();

      my $s2_in_stat = Statistics::Descriptive::Full->new();
      $s2_in_stat->add_data(@s2_inclusion);
      my $s2_in_mean = $s2_in_stat->mean();

      my $s2_ge_stat = Statistics::Descriptive::Full->new();
      $s2_ge_stat->add_data(@s2_gene_expression);
      my $s2_ge_mean = $s2_ge_stat->mean();

      my $in_s1_vs_s2 = logn(($s1_in_mean/$s1_ge_mean),2) - logn(($s2_in_mean/$s2_ge_mean),2);
      my $ex_s1_vs_s2 = logn(($s1_exclusion[0]/$s1_ge_mean),2) - logn(($s2_exclusion[0]/$s2_ge_mean),2);
      my $in_ex_diff = abs($in_s1_vs_s2 - $ex_s1_vs_s2);

      print "\n\tProbe: $probe_id";
      print "\n\tRatio of Means.\tIn: $in_s1_vs_s2\tEx: $ex_s1_vs_s2\tDiff: $in_ex_diff\tn = $probes{$probe_id}{exon_skip_ncount}";

      $probes{$probe_id}{exon_skip_in} = $in_s1_vs_s2;          #ratio of means method
      $probes{$probe_id}{exon_skip_ex} = $ex_s1_vs_s2;          #ratio of means method
      $probes{$probe_id}{exon_skip_in_ex_diff} = $in_ex_diff;   #ratio of means method

      #3-B.) mean of ratios
      my $s12_in_stat = Statistics::Descriptive::Full->new();
      $s12_in_stat->add_data(@inclusion_ratios);
      my $s12_in_mean = $s12_in_stat->mean();

      my $s12_ge_stat = Statistics::Descriptive::Full->new();
      $s12_ge_stat->add_data(@gene_expression_ratios);
      my $s12_ge_mean = $s12_ge_stat->mean();

      my $in_s1_vs_s2_ratios = logn($s12_in_mean,2) - logn($s12_ge_mean,2);
      my $ex_s1_vs_s2_ratios = logn($exclusion_ratios[0],2) - logn($s12_ge_mean,2);
      my $in_ex_diff_ratios = abs($in_s1_vs_s2_ratios - $ex_s1_vs_s2_ratios);

      print "\n\tMean of ratios.\tIn: $in_s1_vs_s2_ratios\tEx: $ex_s1_vs_s2_ratios\tDiff: $in_ex_diff_ratios\tn = $probes{$probe_id}{exon_skip_ncount}";


      $probes{$probe_id}{exon_skip_in_ratios} = $in_s1_vs_s2_ratios;        #mean of ratios method
      $probes{$probe_id}{exon_skip_ex_ratios} = $ex_s1_vs_s2_ratios;        #mean of ratios method
      $probes{$probe_id}{exon_skip_in_ex_diff_ratios} = $in_ex_diff_ratios; #mean of ratios method

    }#Reference exon-skip probe looop

  }#Gene loop

  return();
}


###########################################################################################################################
#4.) Create output file                                                                                                   #
###########################################################################################################################
sub printOutputFile{
  my %args = @_;
  my $out_file = $args{'-out_file'};

  print "\n\nPrinting output file: $out_file\n\n";

  open (OUT, ">$out_file") || die "\nCould not open output file: $out_file\n\n";

  print OUT "$header_line\texon_skip_ncount\texon_skip_in\texon_skip_ex\texon_skip_in_ex_diff\texon_skip_in_ratios\texon_skip_ex_ratios\texon_skip_in_ex_diff_ratios\n";

  foreach my $probe_id (sort {$probes{$a}->{line_count} <=> $probes{$b}->{line_count}} keys %probes){
    print OUT "$probes{$probe_id}{line_record}\t$probes{$probe_id}{exon_skip_ncount}\t$probes{$probe_id}{exon_skip_in}\t$probes{$probe_id}{exon_skip_ex}\t$probes{$probe_id}{exon_skip_in_ex_diff}\t$probes{$probe_id}{exon_skip_in_ratios}\t$probes{$probe_id}{exon_skip_ex_ratios}\t$probes{$probe_id}{exon_skip_in_ex_diff_ratios}\n";
  }

  close (OUT);

  return();
}





