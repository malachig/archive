#!/usr/bin/perl -w
#Written by Malachi Griffith
#This script takes a tab-delimited input file containing microarray data and creates UCSC tracks to represent
#expression values for ALEXA probes and genes.  The first column should be an ALEXA probe ID.
#Other columns contain various raw and background corrected probe intenstity values, as well as DE estimates, etc.

#Process:
#1.) Parse input file and get probe IDs and expression values
#    - Example files: Probes_LnBr_BGC1-2_DE_24.txt, Probes_LnBr_BGC1-2_DE_60.txt, GeneStats_24.txt, GeneStats_60.txt
#    - create a probe hash, with the various scores
#    - create a gene hash, with the various scores
#2.) Use ALEXA to get genomic coordinates for each probe
#    - for each probe ID, get the gene coordinates for the probe position and convert to genomic coordinates
#    - for each gene, display the exon content for the gene (will be used to display the estimate of gene expression at the gene level)
#3.) Create UCSC tracks to display the following data
#    A) BGC1 ratios for each tissue
#    B) BGC2 ratios for each tissue
#    C) DE for each probe (tissue A versus tissue B)
#    D) DE at the gene level (tissue A versus tissue B)
#    E) DE probes that are 2-fold or greater


use strict;
use Data::Dumper;
use Getopt::Long;

#Always use the seqdev folder as the library root and then specify the packages to use.
use lib '/usr/local/ulib/beta/seqdev';
use utilities::utility qw(:all);

use lib '/home/malachig/AlternativeSplicing/perl_bin';
use utilities::ALEXA_DB qw(:all);

#Initialize command line options
my $probe_file = '';
my $gene_file = '';
my $hybe_type = '';
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $out_dir = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'probe_file=s'=>\$probe_file, 'gene_file=s'=>\$gene_file, 'hybe_type=s'=>\$hybe_type, 'out_dir=s'=>\$out_dir);

#Provide instruction to the user
print "\n\nUsage:";
print "\n\tSpecify the database and server to query using: --database and --server";
print "\n\tSpecify the user and password for access using: --user and --password";
print "\n\tSpecify the input probe microarray data file using: --probe_file";
print "\n\tSpecify the input gene summary data file using: --gene_file";
print "\n\tSpecify the hybe type code for labeling purposes (eg. 24 or 60) using: --hybe_type";
print "\n\tSpecify the output dir where UCSC track will be written using: --out_dir";
print "\n\nExample: createArrayDataUCSC_tracks.pl --database ALEXA_hs_31_35d --server jango.bcgsc.ca --user malachig --password pwd --probe_file Probes_LnBr_BGC1-2_DE_24.txt --gene_file GeneStats_24.txt --hybe_type 24 --out_dir /home/malachig/AlternativeSplicing/Array_analysis/ucsc_tracks\n\n";

#Make sure all options were specified
unless ($database && $server && $user && $password && $probe_file && $gene_file && $hybe_type && $out_dir){
  print "\nOptions missing!\n\n";
  exit();
}

#1.) Parse input probe and gene data files
my %probes;
my %genes;
&parseDataFiles('-probe_file'=>$probe_file, '-gene_file'=>$gene_file);

#2.) Get genomic coordinates for each probe and the gene/exon structure

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = connectDB($database, $server, $user, $password);

#For each gene, go through all the probes of that gene and get genomic coordinates
&getGenomicCoordinates('-dbh'=>$alexa_dbh);

#Close database connection
$alexa_dbh->disconnect();

#3.) Create UCSC tracks to display the following data
&generateUCSC_tracks('-ucsc_dir'=>$out_dir);


exit();


#######################################################################################################################
#Parse input probe and gene data files                                                                                #
#######################################################################################################################
sub parseDataFiles{
  my %args = @_;

  my $probe_file = $args{'-probe_file'};
  my $gene_file = $args{'-gene_file'};

  #Open input file
  open (PROBE, "$probe_file") || die "\nCould not open probe file: $probe_file\n\n";

  #Build probes object
  my $header_line = '';
  my $header = 1;
  my @headers;
  my %headers;
  my $probe_count = 0;

  print "\nProcessing input probe file: $probe_file\n\n";

  while (<PROBE>){
    chomp($_);

    #Get the header positions and names
    if ($header == 1){
      $header_line = $_;
      $header = 0;
      @headers = split ("\t", $header_line);
      my $header_count = 0;
      foreach my $h (@headers){
	$headers{$h}{pos} = $header_count;
	$header_count++;
      }
      next();
    }
    $probe_count++;

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
    my $field_name;
    $probes{$probe_id}{AlexaGene_ID} = $gene_id;
    $probes{$probe_id}{Probe_type} = $line[$headers{'Probe_type'}{pos}];
    $probes{$probe_id}{Tm} = $line[$headers{'Tm'}{pos}];
    $field_name = "LnCAP_"."$hybe_type"."_bgc_ratio";
    $probes{$probe_id}{LnCAP_bgc_ratio} = $line[$headers{$field_name}{pos}];
    $field_name = "Brain_"."$hybe_type"."_bgc_ratio";
    $probes{$probe_id}{Brain_bgc_ratio} = $line[$headers{$field_name}{pos}];
    $field_name = "LnCAP_"."$hybe_type"."_bgc2_ratio";
    $probes{$probe_id}{LnCAP_bgc2_ratio} = $line[$headers{$field_name}{pos}];
    $field_name = "Brain_"."$hybe_type"."_bgc2_ratio";
    $probes{$probe_id}{Brain_bgc2_ratio} = $line[$headers{$field_name}{pos}];
    $field_name = "DE_"."$hybe_type"."_bgc";
    $probes{$probe_id}{DE_bgc} = $line[$headers{$field_name}{pos}];
    $field_name = "DE_"."$hybe_type"."_bgc2";
    $probes{$probe_id}{DE_bgc2} = $line[$headers{$field_name}{pos}];

    #Create abreviations for each probe type
    my $probe_type_abr = '';
    if ($probes{$probe_id}{Probe_type} eq "Exon"){
      $probe_type_abr = "E";
    }elsif($probes{$probe_id}{Probe_type} eq "Exon-Exon"){
      $probe_type_abr = "EE"
    }elsif($probes{$probe_id}{Probe_type} eq "Exon-Intron"){
      $probe_type_abr = "EI";
    }elsif($probes{$probe_id}{Probe_type} eq "Intron-Exon"){
      $probe_type_abr = "IE";
    }
    $probes{$probe_id}{Probe_type_abr} = $probe_type_abr;
  }
  close PROBE;

  open (GENE, "$gene_file") || die "\nCould not open gene file: $gene_file\n\n";

  #Build probes object
  $header_line = '';
  $header = 1;
  my @gene_headers;
  my %gene_headers;
  my $gene_count = 0;

  print "\nProcessing input gene file: $gene_file\n\n";

  while (<GENE>){
    chomp($_);

    #Get the header positions and names
    if ($header == 1){
      $header_line = $_;
      $header = 0;
      @gene_headers = split ("\t", $header_line);
      my $header_count = 0;
      foreach my $h (@gene_headers){
	$gene_headers{$h}{pos} = $header_count;
	$header_count++;
      }
      next();
    }
    $gene_count++;

    my @line = split ("\t", $_);

    my $gene_id;
    if ($line[0] =~ /(\d+)/){
      $gene_id = $1;
    }

    #Populate probes hash
    $genes{$gene_id}{Entrez_ID} = $line[$gene_headers{'Entrez ID'}{pos}];
    $genes{$gene_id}{Mean_BGC2_Ratio_LnCAP} = $line[$gene_headers{'Mean BGC2 Ratio LnCAP'}{pos}];
    $genes{$gene_id}{Mean_BGC2_Ratio_Brain} = $line[$gene_headers{'Mean BGC2 Ratio Brain'}{pos}];
    $genes{$gene_id}{Mean_DE_BGC2_Ratio} = $line[$gene_headers{'Mean DE BGC2 Ratio'}{pos}];

  }
  close GENE;

  #Go through all probes and associate with gene records
  print "\nAssociating probe and gene info\n\n";
  foreach my $probe_id (keys %probes){
    my $gene_id = $probes{$probe_id}{AlexaGene_ID};

    my %gene_probes;
    if ($genes{$gene_id}{probes}){
      %gene_probes = %{$genes{$gene_id}{probes}};
    }

    $gene_probes{$probe_id}{Probe_type} = $probes{$probe_id}{Probe_type};
    $gene_probes{$probe_id}{Tm} = $probes{$probe_id}{Tm};
    $gene_probes{$probe_id}{LnCAP_bgc_ratio} = $probes{$probe_id}{LnCAP_bgc_ratio};
    $gene_probes{$probe_id}{Brain_bgc_ratio} = $probes{$probe_id}{Brain_bgc_ratio};
    $gene_probes{$probe_id}{LnCAP_bgc2_ratio} = $probes{$probe_id}{LnCAP_bgc2_ratio};
    $gene_probes{$probe_id}{Brain_bgc2_ratio} = $probes{$probe_id}{Brain_bgc2_ratio};
    $gene_probes{$probe_id}{DE_bgc} = $probes{$probe_id}{DE_bgc};
    $gene_probes{$probe_id}{DE_bgc2} = $probes{$probe_id}{DE_bgc2};
    $gene_probes{$probe_id}{Probe_type_abr} = $probes{$probe_id}{Probe_type_abr};
    $genes{$gene_id}{probes} = \%gene_probes;
  }

  return();
}


#######################################################################################################################
#Get genomic coordinates for each probe                                                                               #
#######################################################################################################################
sub getGenomicCoordinates{
  my %args = @_;

  my $dbh = $args{'-dbh'};

  print "\nGetting genomic coordinates for each probe\n\n";
  foreach my $gene_id (sort {$a <=> $b} keys %genes){

    print "\nProcessing Gene: $gene_id";

    my %gene_probes = %{$genes{$gene_id}{probes}};

    #Get probe coordinates
    foreach my $probe_id (keys %gene_probes){
      my %probe = %{&getProbeInfo ('-dbh'=>$dbh, '-probe_count_id'=>$probe_id)};

      #Every probe should have unit1 coordinates
      my %coords_unit1 = %{&convertGeneCoordinates ('-dbh'=>$dbh, '-gene_id'=>$gene_id,
						    '-start_pos'=>$probe{unit1_start},
						    '-end_pos'=>$probe{unit1_end})};

      $gene_probes{$probe_id}{unit1_start} = $coords_unit1{$gene_id}{chr_start};
      $gene_probes{$probe_id}{unit1_end} = $coords_unit1{$gene_id}{chr_end};
      $gene_probes{$probe_id}{strand} = $coords_unit1{$gene_id}{strand};

      #Only exon-exon probes will have unit2 coordinates
      if ($probe{unit2_start} =~ /\d+/ && $probe{unit2_end} =~ /\d+/){
	my %coords_unit2 = %{&convertGeneCoordinates ('-dbh'=>$dbh, '-gene_id'=>$gene_id,
						    '-start_pos'=>$probe{unit2_start},
						    '-end_pos'=>$probe{unit2_end})};

	$gene_probes{$probe_id}{unit2_start} = $coords_unit2{$gene_id}{chr_start};
	$gene_probes{$probe_id}{unit2_end} = $coords_unit2{$gene_id}{chr_end};
      }else{
	$gene_probes{$probe_id}{unit2_start} = $probe{unit2_start};
	$gene_probes{$probe_id}{unit2_end} = $probe{unit2_end};
      }
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

    #Get exon content coordinates - Genomic coordinates representing all exon content for this gene from any transcript
    my %exonContent = %{&getExonContent ('-dbh'=>$dbh, '-gene_id'=>$gene_id)};

    foreach my $ec (keys %exonContent){
      my %coords = %{&convertGeneCoordinates ('-dbh'=>$dbh, '-gene_id'=>$gene_id, '-start_pos'=>$exonContent{$ec}{start}, 
					      '-end_pos'=>$exonContent{$ec}{end})};
      $exonContent{$ec}{chr_start} = $coords{$gene_id}{chr_start};
      $exonContent{$ec}{chr_end} = $coords{$gene_id}{chr_end};
    }

    $genes{$gene_id}{exon_content} = \%exonContent;
    $genes{$gene_id}{probes} = \%gene_probes;

    #DEBUG
    #print Dumper %gene_probes;
    #return();

  }

  return();

}


#################################################################################################################################
#3.) Generate custom track files for the UCSC browser                                                                           #
#################################################################################################################################
sub generateUCSC_tracks{
  my %args = @_;
  my $ucsc_dir = $args{'-ucsc_dir'};

  #First step:
  #Get the max observed value for each number to be displayed so that a score can be calculated on a scale of 0-1000
  #If the values are ratios with a minimum value of 1, get a score on a scale of 0-1000 by:
  #  taking the difference from 1 then divide by the max value's difference from 1 and multiplying by 1000, then add a very small value to prevent scores of 0
  #  Note that DE values are already on a scale of 0 to X so this is not neccesary.

  #Gene Values
  my ($max_mean_bgc2_ratio_lncap,$max_mean_bgc2_ratio_brain,$max_mean_DE_bgc2_ratio) = (0,0,0);

  #Probe Values
  my ($min_tm, $max_tm, $max_lncap_bgc_ratio, $max_brain_bgc_ratio, $max_lncap_bgc2_ratio, $max_brain_bgc2_ratio, $max_de_bgc, $max_de_bgc2) = (100,0,0,0,0,0,0,0);

  foreach my $gene_id (sort {$a <=> $b} keys %genes){
    if ($genes{$gene_id}{Mean_BGC2_Ratio_LnCAP} > $max_mean_bgc2_ratio_lncap){
      $max_mean_bgc2_ratio_lncap = $genes{$gene_id}{Mean_BGC2_Ratio_LnCAP};
    }
    if ($genes{$gene_id}{Mean_BGC2_Ratio_Brain} > $max_mean_bgc2_ratio_brain){
      $max_mean_bgc2_ratio_brain = $genes{$gene_id}{Mean_BGC2_Ratio_Brain};
    }
    if (abs($genes{$gene_id}{Mean_DE_BGC2_Ratio}) > $max_mean_DE_bgc2_ratio){
      $max_mean_DE_bgc2_ratio = abs($genes{$gene_id}{Mean_DE_BGC2_Ratio});
    }
    my %probe_obj = %{$genes{$gene_id}{probes}};

    foreach my $probe_id (keys %probe_obj){

      if ($probe_obj{$probe_id}{Tm} < $min_tm){
	$min_tm = $probe_obj{$probe_id}{Tm};
      }
      if ($probe_obj{$probe_id}{Tm} > $max_tm){
	$max_tm = $probe_obj{$probe_id}{Tm};
      }
      if ($probe_obj{$probe_id}{LnCAP_bgc_ratio} > $max_lncap_bgc_ratio){
	$max_lncap_bgc_ratio = $probe_obj{$probe_id}{LnCAP_bgc_ratio};
      }
      if ($probe_obj{$probe_id}{Brain_bgc_ratio} > $max_brain_bgc_ratio){
	$max_brain_bgc_ratio = $probe_obj{$probe_id}{Brain_bgc_ratio};
      }
      if ($probe_obj{$probe_id}{LnCAP_bgc2_ratio} > $max_lncap_bgc2_ratio){
	$max_lncap_bgc2_ratio = $probe_obj{$probe_id}{LnCAP_bgc2_ratio};
      }
      if ($probe_obj{$probe_id}{Brain_bgc2_ratio} > $max_brain_bgc2_ratio){
	$max_brain_bgc2_ratio = $probe_obj{$probe_id}{Brain_bgc2_ratio};
      }
      if ($probe_obj{$probe_id}{DE_bgc} > $max_de_bgc){
	$max_de_bgc = $probe_obj{$probe_id}{DE_bgc};
      }
      if ($probe_obj{$probe_id}{DE_bgc2} > $max_de_bgc2){
	$max_de_bgc2 = $probe_obj{$probe_id}{DE_bgc2};
      }
    }
  }

  #COLORS
  #LnCAP or LnCAP Over-expressed 100,50,0 (Shades of brown)
  #Brain or Brain Over-expressed 0,0,0 (Shades of black)
  #Splicing Index scores 0,60,120 (Shades of blue)

  foreach my $gene_id (sort {$a <=> $b} keys %genes){

    #Open a UCSC track file for this gene
    my $outfile = "$ucsc_dir/"."ALEXA_"."$gene_id".".txt";

    print "\nCreating: $outfile";

    open (UCSC, ">$outfile") || die "\nCould not open test file: $outfile\n\n";

    #Set default browser settings
    print UCSC "#Browser line";
    print UCSC "\nbrowser position chr$genes{$gene_id}{chromosome}:$genes{$gene_id}{chr_start}-$genes{$gene_id}{chr_end}";  #Default start pos
    print UCSC "\nbrowser dense exoniphy intronEst est";
    print UCSC "\nbrowser full knownGene refGene ensGene";
    print UCSC "\nbrowser pack multiz8way";


    #A.) Gene Level Expression and Differential Expression Tracks

    #A-1.) LnCAP
    print UCSC "\n\n#Gene Level Expression LnCAP $hybe_type";
    print UCSC "\ntrack name=mean_BGC2_ratio_LnCAP description=\"Mean BGC2-Ratio LnCAP RNA $hybe_type-mer (max = $max_mean_bgc2_ratio_lncap)\" color=100,50,0 useScore=1 visibility=2";

    #Data lines: seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup
    print UCSC "\n#BEGIN DATA - Mean_BGC2_Ratio_LnCAP";

    #Get the value and calculate a score
    my $mean_bgc2_ratio_lncap = $genes{$gene_id}{Mean_BGC2_Ratio_LnCAP};

    my $feature_name = "$gene_id"."_"."($mean_bgc2_ratio_lncap)";
    my $score = ((($mean_bgc2_ratio_lncap-1)/($max_mean_bgc2_ratio_lncap-1))*1000)+0.0000000001;

    #Go through each segment of this gene's exon content
    my %exonContent = %{$genes{$gene_id}{exon_content}};
    foreach my $ec (sort {$a <=> $b} keys %exonContent){
      print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tGene\t$exonContent{$ec}{chr_start}\t$exonContent{$ec}{chr_end}\t$score\t$genes{$gene_id}{chr_strand}\t.\t$feature_name";
    }

    #A-2.) Brain
    print UCSC "\n\n#Gene Level Expression Brain $hybe_type";
    print UCSC "\ntrack name=mean_BGC2_ratio_Brain description=\"Mean BGC2-Ratio Brain RNA $hybe_type-mer (max = $max_mean_bgc2_ratio_brain)\" color=0,0,0 useScore=1 visibility=2";

    print UCSC "\n#BEGIN DATA - Mean_BGC2_Ratio_Brain";

    #Get the value and calculate a score
    my $mean_bgc2_ratio_brain = $genes{$gene_id}{Mean_BGC2_Ratio_Brain};

    $feature_name = "$gene_id"."_"."($mean_bgc2_ratio_brain)";
    $score = ((($mean_bgc2_ratio_brain-1)/($max_mean_bgc2_ratio_brain-1))*1000)+0.0000000001;

    #Go through each segment of this gene's exon content
    foreach my $ec (sort {$a <=> $b} keys %exonContent){
      print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tGene\t$exonContent{$ec}{chr_start}\t$exonContent{$ec}{chr_end}\t$score\t$genes{$gene_id}{chr_strand}\t.\t$feature_name";
    }

    #A-3.) DE LnCAP vs. Brain
    #Divide into two categories (up in LnCAP and up in Brain)
    #UP in LnCAP
    if ($genes{$gene_id}{Mean_DE_BGC2_Ratio} >= 0){
      print UCSC "\n\n#Gene Level Differential Expression - OverExpressed in LnCAP $hybe_type";
      print UCSC "\ntrack name=mean_DE_Over_LnCAP description=\"Mean Log2 BGC2-Ratio OverExpressed in LnCAP RNA $hybe_type-mer (max = $max_mean_DE_bgc2_ratio)\" color=100,50,0 useScore=1 visibility=2";

      #Data lines: seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup
      print UCSC "\n#BEGIN DATA - Mean_DE_BGC2_Ratio OverExpressed in LnCAP";

      #Get the value and calculate a score
      my $mean_DE_bgc2_ratio_over_lncap = abs($genes{$gene_id}{Mean_DE_BGC2_Ratio});

      $feature_name = "$gene_id"."_"."($genes{$gene_id}{Mean_DE_BGC2_Ratio})";
      $score = ($mean_DE_bgc2_ratio_over_lncap/$max_mean_DE_bgc2_ratio)*1000;

      #Go through each segment of this gene's exon content
      foreach my $ec (sort {$a <=> $b} keys %exonContent){
	print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tGene\t$exonContent{$ec}{chr_start}\t$exonContent{$ec}{chr_end}\t$score\t$genes{$gene_id}{chr_strand}\t.\t$feature_name";
      }
    }else{
    #UP in Brain

      print UCSC "\n\n#Gene Level Differential Expression - OverExpressed in Brain $hybe_type";
      print UCSC "\ntrack name=mean_DE_Over_Brain description=\"Mean Log2 BGC2-Ratio OverExpressed in Brain RNA $hybe_type-mer (max = $max_mean_DE_bgc2_ratio)\" color=0,0,0 useScore=1 visibility=2";

      #Data lines: seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup
      print UCSC "\n#BEGIN DATA - Mean_DE_BGC2_Ratio OverExpressed in Brain";

      #Get the value and calculate a score
      my $mean_DE_bgc2_ratio_over_brain = abs($genes{$gene_id}{Mean_DE_BGC2_Ratio});

      $feature_name = "$gene_id"."_"."($genes{$gene_id}{Mean_DE_BGC2_Ratio})";
      $score = ($mean_DE_bgc2_ratio_over_brain/$max_mean_DE_bgc2_ratio)*1000;

      #Go through each segment of this gene's exon content
      foreach my $ec (sort {$a <=> $b} keys %exonContent){
	print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tGene\t$exonContent{$ec}{chr_start}\t$exonContent{$ec}{chr_end}\t$score\t$genes{$gene_id}{chr_strand}\t.\t$feature_name";
      }
    }


    #B Probe Level Expression and Differential Expression Tracks
    my %probe_obj = %{$genes{$gene_id}{probes}};

    #B-1.) Probe Tm's
    print UCSC "\n\n#Probe Tms";
    print UCSC "\ntrack name=probe_tm description=\"Probe Tm (min = $min_tm  max = $max_tm)\" color=0,60,120 useScore=1 visibility=2";

    #Data lines: seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup
    print UCSC "\n#BEGIN DATA - Probe Tms";

    foreach my $probe_id (sort {$a <=> $b} keys %probe_obj){

      #Get the value and calculate a score
      my $probe_tm = $probe_obj{$probe_id}{Tm};
      my $x = sprintf("%.1f", $probe_tm);
      my $feature_name = "$probe_obj{$probe_id}{Probe_type_abr}"."$probe_id"."_"."($x)";
      my $score = ((($probe_tm-$min_tm)/($max_tm-$min_tm))*1000)+0.0000000001;
      print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit1_start}\t$probe_obj{$probe_id}{unit1_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      if ($probe_obj{$probe_id}{unit2_start} =~ /\d+/ && $probe_obj{$probe_id}{unit2_end} =~ /\d+/){
	print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit2_start}\t$probe_obj{$probe_id}{unit2_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      }
    }

    #B-2.) LnCAP_bgc_ratio
    print UCSC "\n\n#Probe Level Expression - LnCAP BGC Ratio";
    my $max_v = sprintf("%.3f", $max_lncap_bgc_ratio);
    print UCSC "\ntrack name=lnCAP_bgc_ratio description=\"BGC-Ratio LnCAP RNA $hybe_type-mer (max = $max_v)\" color=100,50,0 useScore=1 visibility=2";

    #Data lines: seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup
    print UCSC "\n#BEGIN DATA - BGC_Ratio_LnCAP";

    my $gene_probe_count = 0;
    foreach my $probe_id (sort {$a <=> $b} keys %probe_obj){
      $gene_probe_count++;

      #Get the value and calculate a score
      my $lncap_bgc_ratio = $probe_obj{$probe_id}{LnCAP_bgc_ratio};
      my $x = sprintf("%.1f", $lncap_bgc_ratio);
      my $feature_name = "$probe_obj{$probe_id}{Probe_type_abr}"."_$gene_probe_count"."_"."($x)";
      my $score = ((($lncap_bgc_ratio)/($max_lncap_bgc_ratio))*1000);

      print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit1_start}\t$probe_obj{$probe_id}{unit1_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      if ($probe_obj{$probe_id}{unit2_start} =~ /\d+/ && $probe_obj{$probe_id}{unit2_end} =~ /\d+/){
	print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit2_start}\t$probe_obj{$probe_id}{unit2_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      }
    }

    #B-3.) Brain_bgc_ratio
    print UCSC "\n\n#Probe Level Expression - Brain BGC Ratio";
    $max_v = sprintf("%.3f", $max_brain_bgc_ratio);
    print UCSC "\ntrack name=Brain_bgc_ratio description=\"BGC-Ratio Brain RNA $hybe_type-mer (max = $max_v)\" color=0,0,0 useScore=1 visibility=2";

    #Data lines: seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup
    print UCSC "\n#BEGIN DATA - BGC_Ratio_Brain";

    $gene_probe_count = 0;
    foreach my $probe_id (sort {$a <=> $b} keys %probe_obj){
      $gene_probe_count++;
      #Get the value and calculate a score
      my $brain_bgc_ratio = $probe_obj{$probe_id}{Brain_bgc_ratio};
      my $x = sprintf("%.1f", $brain_bgc_ratio);
      my $feature_name = "$probe_obj{$probe_id}{Probe_type_abr}"."_$gene_probe_count"."_"."($x)";
      my $score = ((($brain_bgc_ratio)/($max_brain_bgc_ratio))*1000);

      print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit1_start}\t$probe_obj{$probe_id}{unit1_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      if ($probe_obj{$probe_id}{unit2_start} =~ /\d+/ && $probe_obj{$probe_id}{unit2_end} =~ /\d+/){
	print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit2_start}\t$probe_obj{$probe_id}{unit2_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      }
    }

    #B-4.) LnCAP_bgc2_ratio
    print UCSC "\n\n#Probe Level Expression - LnCAP BGC2 Ratio";
    $max_v = sprintf("%.3f", $max_lncap_bgc2_ratio);
    print UCSC "\ntrack name=LnCAP_bgc2_ratio description=\"BGC2-Ratio LnCAP RNA $hybe_type-mer (max = $max_v)\" color=100,50,0 useScore=1 visibility=2";

    #Data lines: seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup
    print UCSC "\n#BEGIN DATA - BGC2_Ratio_LnCAP";

    $gene_probe_count = 0;
    foreach my $probe_id (sort {$a <=> $b} keys %probe_obj){
      $gene_probe_count++;

      #Get the value and calculate a score
      my $lncap_bgc2_ratio = $probe_obj{$probe_id}{LnCAP_bgc2_ratio};
      my $x = sprintf("%.1f", $lncap_bgc2_ratio);
      my $feature_name = "$probe_obj{$probe_id}{Probe_type_abr}"."_$gene_probe_count"."_"."($x)";
      my $score = ((($lncap_bgc2_ratio)/($max_lncap_bgc2_ratio))*1000);

      print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit1_start}\t$probe_obj{$probe_id}{unit1_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      if ($probe_obj{$probe_id}{unit2_start} =~ /\d+/ && $probe_obj{$probe_id}{unit2_end} =~ /\d+/){
	print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit2_start}\t$probe_obj{$probe_id}{unit2_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      }
    }

    #B-5.) Brain_bgc2_ratio
    print UCSC "\n\n#Probe Level Expression - Brain BGC2 Ratio";
    $max_v = sprintf("%.3f", $max_brain_bgc2_ratio);
    print UCSC "\ntrack name=Brain_bgc2_ratio description=\"BGC2-Ratio Brain RNA $hybe_type-mer (max = $max_v)\" color=0,0,0 useScore=1 visibility=2";

    #Data lines: seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup
    print UCSC "\n#BEGIN DATA - BGC2_Ratio_Brain";

    $gene_probe_count = 0;
    foreach my $probe_id (sort {$a <=> $b} keys %probe_obj){
      $gene_probe_count++;

      #Get the value and calculate a score
      my $brain_bgc2_ratio = $probe_obj{$probe_id}{Brain_bgc2_ratio};
      my $x = sprintf("%.1f", $brain_bgc2_ratio);
      my $feature_name = "$probe_obj{$probe_id}{Probe_type_abr}"."_$gene_probe_count"."_"."($x)";
      my $score = ((($brain_bgc2_ratio)/($max_brain_bgc2_ratio))*1000);

      print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit1_start}\t$probe_obj{$probe_id}{unit1_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      if ($probe_obj{$probe_id}{unit2_start} =~ /\d+/ && $probe_obj{$probe_id}{unit2_end} =~ /\d+/){
	print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit2_start}\t$probe_obj{$probe_id}{unit2_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      }
    }

    #B-6-I.) Probe DE BGC-Ratios LnCAP versus Brain - Overexpressed in LnCAP
    print UCSC "\n\n#Probe Level Differential Expression - LnCAP versus Brain BGC Ratios - Overexpressed in LnCAP";
    $max_v = sprintf("%.3f", $max_de_bgc);
    print UCSC "\ntrack name=DE_bgc_ratio_LnCAP_OverExpressed description=\"DE Log2 BGC-Ratio LnCAP Over-Expressed $hybe_type-mer (max = $max_v)\" color=100,50,0 useScore=1 visibility=2";

    #Data lines: seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup
    print UCSC "\n#BEGIN DATA - BGC Ratio LnCAP Overexpressed";

    $gene_probe_count = 0;
    foreach my $probe_id (sort {$a <=> $b} keys %probe_obj){
      $gene_probe_count++;

      #Get the value and calculate a score
      my $de_bgc_ratio = $probe_obj{$probe_id}{DE_bgc};

      if ($de_bgc_ratio < 0){
	next();
      }

      my $x = sprintf("%.1f", $de_bgc_ratio);
      my $feature_name = "$probe_obj{$probe_id}{Probe_type_abr}"."_$gene_probe_count"."_"."($x)";
      my $score = ((abs($de_bgc_ratio)/$max_de_bgc)*1000)+0.0000000001;

      print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit1_start}\t$probe_obj{$probe_id}{unit1_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      if ($probe_obj{$probe_id}{unit2_start} =~ /\d+/ && $probe_obj{$probe_id}{unit2_end} =~ /\d+/){
	print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit2_start}\t$probe_obj{$probe_id}{unit2_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      }
    }

    #B-6-II.) Probe DE BGC-Ratios LnCAP versus Brain - Overexpressed in Brain
    print UCSC "\n\n#Probe Level Differential Expression - LnCAP versus Brain BGC Ratios - Overexpressed in Brain";
    $max_v = sprintf("%.3f", $max_de_bgc);
    print UCSC "\ntrack name=DE_bgc_ratio_Brain_OverExpressed description=\"DE Log2 BGC-Ratio Brain Over-Expressed $hybe_type-mer (max = $max_v)\" color=0,0,0 useScore=1 visibility=2";

    #Data lines: seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup
    print UCSC "\n#BEGIN DATA - BGC Ratio Brain Overexpressed";

    $gene_probe_count = 0;
    foreach my $probe_id (sort {$a <=> $b} keys %probe_obj){
      $gene_probe_count++;
      #Get the value and calculate a score
      my $de_bgc_ratio = $probe_obj{$probe_id}{DE_bgc};

      if ($de_bgc_ratio >= 0){
	next();
      }

      my $x = sprintf("%.1f", $de_bgc_ratio);
      my $feature_name = "$probe_obj{$probe_id}{Probe_type_abr}"."_$gene_probe_count"."_"."($x)";
      my $score = ((abs($de_bgc_ratio)/$max_de_bgc)*1000)+0.0000000001;

      print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit1_start}\t$probe_obj{$probe_id}{unit1_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      if ($probe_obj{$probe_id}{unit2_start} =~ /\d+/ && $probe_obj{$probe_id}{unit2_end} =~ /\d+/){
	print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit2_start}\t$probe_obj{$probe_id}{unit2_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      }
    }


    #B-7-I.) Probe DE BGC2-Ratios LnCAP versus Brain - Overexpressed in LnCAP
    print UCSC "\n\n#Probe Level Differential Expression - LnCAP versus Brain BGC2 Ratios - Overexpressed in LnCAP";
    $max_v = sprintf("%.3f", $max_de_bgc2);
    print UCSC "\ntrack name=DE_bgc2_ratio_LnCAP_OverExpressed description=\"DE Log2 BGC2-Ratio LnCAP Over-Expressed $hybe_type-mer (max = $max_v)\" color=100,50,0 useScore=1 visibility=2";

    #Data lines: seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup
    print UCSC "\n#BEGIN DATA - BGC2 Ratio LnCAP Overexpressed";

    $gene_probe_count = 0;
    foreach my $probe_id (sort {$a <=> $b} keys %probe_obj){
      $gene_probe_count++;

      #Get the value and calculate a score
      my $de_bgc2_ratio = $probe_obj{$probe_id}{DE_bgc2};

      if ($de_bgc2_ratio < 0){
	next();
      }

      my $x = sprintf("%.1f", $de_bgc2_ratio);
      my $feature_name = "$probe_obj{$probe_id}{Probe_type_abr}"."_$gene_probe_count"."_"."($x)";
      my $score = ((abs($de_bgc2_ratio)/$max_de_bgc2)*1000)+0.0000000001;

      print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit1_start}\t$probe_obj{$probe_id}{unit1_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      if ($probe_obj{$probe_id}{unit2_start} =~ /\d+/ && $probe_obj{$probe_id}{unit2_end} =~ /\d+/){
	print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit2_start}\t$probe_obj{$probe_id}{unit2_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      }
    }


    #B-7-II.) SIGNIFICANT Probe DE BGC2-Ratios LnCAP versus Brain - Overexpressed in LnCAP
    my $sig_level = 1.5;
    print UCSC "\n\n#Probe Level Differential Expression - LnCAP versus Brain BGC2 Ratios - Overexpressed in LnCAP ($sig_level FOLD)";
    $max_v = sprintf("%.3f", $max_de_bgc2);
    print UCSC "\ntrack name=DE_bgc2_ratio_LnCAP_OverExpressed_SIG description=\"DE Log2 BGC2-Ratio LnCAP Over-Expressed $hybe_type-mer (max = $max_v) ($sig_level FOLD)\" color=100,50,0 useScore=1 visibility=2";

    #Data lines: seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup
    print UCSC "\n#BEGIN DATA - BGC2 Ratio LnCAP Overexpressed ($sig_level FOLD)";

    $gene_probe_count = 0;
    foreach my $probe_id (sort {$a <=> $b} keys %probe_obj){
      $gene_probe_count++;

      #Get the value and calculate a score
      my $de_bgc2_ratio = $probe_obj{$probe_id}{DE_bgc2};

      if ($de_bgc2_ratio < $sig_level){
	next();
      }

      my $x = sprintf("%.1f", $de_bgc2_ratio);
      my $feature_name = "$probe_obj{$probe_id}{Probe_type_abr}"."_$gene_probe_count"."_"."($x)";
      my $score = ((abs($de_bgc2_ratio)/$max_de_bgc2)*1000)+0.0000000001;

      print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit1_start}\t$probe_obj{$probe_id}{unit1_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      if ($probe_obj{$probe_id}{unit2_start} =~ /\d+/ && $probe_obj{$probe_id}{unit2_end} =~ /\d+/){
	print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit2_start}\t$probe_obj{$probe_id}{unit2_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      }
    }

    #B-7-III.) Probe DE BGC2-Ratios LnCAP versus Brain - Overexpressed in Brain
    print UCSC "\n\n#Probe Level Differential Expression - LnCAP versus Brain BGC2 Ratios - Overexpressed in Brain";
    $max_v = sprintf("%.3f", $max_de_bgc2);
    print UCSC "\ntrack name=DE_bgc2_ratio_Brain_OverExpressed description=\"DE Log2 BGC2-Ratio Brain Over-Expressed $hybe_type-mer (max = $max_v)\" color=0,0,0 useScore=1 visibility=2";

    #Data lines: seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup
    print UCSC "\n#BEGIN DATA - BGC2 Ratio Brain Overexpressed";

    $gene_probe_count = 0;
    foreach my $probe_id (sort {$a <=> $b} keys %probe_obj){
      $gene_probe_count++;

      #Get the value and calculate a score
      my $de_bgc2_ratio = $probe_obj{$probe_id}{DE_bgc2};

      if ($de_bgc2_ratio >= 0){
	next();
      }

      my $x = sprintf("%.1f", $de_bgc2_ratio);
      my $feature_name = "$probe_obj{$probe_id}{Probe_type_abr}"."_$gene_probe_count"."_"."($x)";
      my $score = ((abs($de_bgc2_ratio)/$max_de_bgc2)*1000)+0.0000000001;

      print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit1_start}\t$probe_obj{$probe_id}{unit1_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      if ($probe_obj{$probe_id}{unit2_start} =~ /\d+/ && $probe_obj{$probe_id}{unit2_end} =~ /\d+/){
	print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit2_start}\t$probe_obj{$probe_id}{unit2_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      }
    }

    #B-7-IV.) SIGNIFICANT Probe DE BGC2-Ratios LnCAP versus Brain - Overexpressed in Brain
    $sig_level = -1.5;
    print UCSC "\n\n#Probe Level Differential Expression - LnCAP versus Brain BGC2 Ratios - Overexpressed in Brain ($sig_level FOLD)";
    $max_v = sprintf("%.3f", $max_de_bgc2);
    print UCSC "\ntrack name=DE_bgc2_ratio_Brain_OverExpressed_SIG description=\"DE Log2 BGC2-Ratio Brain Over-Expressed $hybe_type-mer (max = $max_v) ($sig_level FOLD)\" color=0,0,0 useScore=1 visibility=2";

    #Data lines: seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tgroup
    print UCSC "\n#BEGIN DATA - BGC2 Ratio Brain Overexpressed ($sig_level FOLD)";

    $gene_probe_count = 0;
    foreach my $probe_id (sort {$a <=> $b} keys %probe_obj){
      $gene_probe_count++;

      #Get the value and calculate a score
      my $de_bgc2_ratio = $probe_obj{$probe_id}{DE_bgc2};

      if ($de_bgc2_ratio > $sig_level){
	next();
      }

      my $x = sprintf("%.1f", $de_bgc2_ratio);
      my $feature_name = "$probe_obj{$probe_id}{Probe_type_abr}"."_$gene_probe_count"."_"."($x)";
      my $score = ((abs($de_bgc2_ratio)/$max_de_bgc2)*1000)+0.0000000001;

      print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit1_start}\t$probe_obj{$probe_id}{unit1_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      if ($probe_obj{$probe_id}{unit2_start} =~ /\d+/ && $probe_obj{$probe_id}{unit2_end} =~ /\d+/){
	print UCSC "\nchr$genes{$gene_id}{chromosome}\tALEXA\tProbe\t$probe_obj{$probe_id}{unit2_start}\t$probe_obj{$probe_id}{unit2_end}\t$score\t$probe_obj{$probe_id}{strand}\t.\t$feature_name";
      }
    }

    close (UCSC);

    #DEBUG
    #return();

  }
  return();
}


