#!/usr/bin/perl -w
# Malachi Griffith

#Use after using R to create a SVGs for each splice site mutations and corresponding exon-exon junction expression
#Create a webpage that summarizes the mutations
#First display the number of each class of event:
#1.) Reciprocal; 2.) Gain of alternative junction; 3.) Loss of canonical junction; 4.) Gain of canonical junction 5.) Unassigned
#Next create a ranked table for each of these categories
#Each table should contain:
#Mutation ID, gene, library list, acceptor/donor, display coordinates, sort value, links to SVGs (both clickable and non-clickable versions)

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;
use Benchmark;

my $infile = '';
my $outdir = '';
GetOptions ('infile=s'=>\$infile, 'outdir=s'=>\$outdir);

if ($infile && $outdir){
  #Parameters ok
}else{
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  print GREEN, "\nSpecify the name of the input file using: --infile", RESET;
  print GREEN, "\nSpecify the name of the output file using: --outdir", RESET;
  print GREEN, "\n\nExample: summarizeSpliceSiteMutations.pl  --infile=/projects/alexa2/hmmSplicer/SA_TN_Breast/summary/spliceSiteMutations/figures/RankedList.tsv  --outdir=/projects/alexa2/hmmSplicer/SA_TN_Breast/summary/spliceSiteMutations/figures/\n\n", RESET;
  exit();
}

unless(-e $infile){
  print RED, "\n\nNo such input file: $infile", RESET;
  exit();
}
unless ($outdir =~ /\/$/){
  $outdir .= "/";
}
unless (-e $outdir && -d $outdir){
  print RED, "\n\nNo such output dir: $outdir", RESET;
  exit();
}

my $title = "Splice related somatic mutations and their affect on RNA expression (exon-exon junction usage)";

#Import the input file
my $mut_ref = &importMutations();
#print Dumper $mut_ref;

#Get counts for each of the five classes;
my $reciprocal_count = 0;
my $gain_alter_count = 0;
my $loss_canon_count = 0;
my $gain_canon_count = 0;
my $unassigned_count = 0;
foreach my $mid (keys %{$mut_ref}){
  my $class = $mut_ref->{$mid}->{class};
  if ($class eq "RECIPROCAL"){$reciprocal_count++;}
  if ($class eq "GAIN_ALTER"){$gain_alter_count++;}
  if ($class eq "LOSS_CANON"){$loss_canon_count++;}
  if ($class eq "GAIN_CANON"){$gain_canon_count++;}
  if ($class eq "UNCHANGED"){$unassigned_count++;}
}

print "\n\nreciprocal_count = $reciprocal_count\tgain_alter_count = $gain_alter_count\tloss_canon_count = $loss_canon_count\tgain_canon_count = $gain_canon_count\tunassigned_count = $unassigned_count\n\n";


#Process each of the five classes and create an html table using the appropriate sort field
my $table1_content = &createTable('-class'=>"RECIPROCAL", '-sort_column'=>"max_recip_diff");
my $table2_content = &createTable('-class'=>"LOSS_CANON", '-sort_column'=>"max_loss_canon");
my $table3_content = &createTable('-class'=>"GAIN_ALTER", '-sort_column'=>"max_gain_alter");
my $table4_content = &createTable('-class'=>"GAIN_CANON", '-sort_column'=>"max_gain_canon");
my $table5_content = &createTable('-class'=>"UNCHANGED", '-sort_column'=>"avg_median_diff");


#Copy neccessary html files to the outdir
system("cp -f /home/malachig/svn/alexa_seq/website/web_files/ALEXA2.css .");
system("cp -f /home/malachig/svn/alexa_seq/website/web_files/jquery.min.js .");
system("cp -f /home/malachig/svn/alexa_seq/website/web_files/animatedcollapse.js .");
system("cp -f /home/malachig/svn/alexa_seq/website/web_files/images/Plus_icon.gif .");
system("cp -f /home/malachig/svn/alexa_seq/website/web_files/images/Minus_icon.gif .");

#Write out the html output file
my $outfile = "$outdir"."index.htm";
open (HTML, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";

print HTML <<"EOF";
<HTML>
<HEAD>
<TITLE>$title</TITLE>
<link rel="stylesheet" type="text/css" href="ALEXA2.css">

<!-- Expandable area code -->
<script type="text/javascript" src="jquery.min.js"></script>
<script type="text/javascript" src="animatedcollapse.js"></script>
<script type="text/javascript">
animatedcollapse.addDiv('table1', 'fade=1, speed=500, hide=1 persist=1')
animatedcollapse.addDiv('table2', 'fade=1, speed=500, hide=1 persist=1')
animatedcollapse.addDiv('table3', 'fade=1, speed=500, hide=1 persist=1')
animatedcollapse.addDiv('table4', 'fade=1, speed=500, hide=1 persist=1')
animatedcollapse.addDiv('table5', 'fade=1, speed=500, hide=1 persist=1')
animatedcollapse.ontoggle=function(\$, divobj, state){}
animatedcollapse.init()
</script>

</HEAD>
<BODY>

<!-- Summary -->
<BR>
<P CLASS=\"Indented12LR_s19_bold\">$title</P><BR>
<P CLASS=\"Indented12LR_s16_bold\">Summary of splicing or expression changes putatively associated with somatic mutations</P>
<P CLASS=\"Indented12LR_s16\">
Reciprocal = $reciprocal_count (Loss of a canonical junction and corresponding gain of an alternative junction)<BR>
Loss canonical = $loss_canon_count (Loss of a canonical junction)<BR>
Gain alternative = $gain_alter_count (Gain of an alternative junction only)<BR>
Gain canonical = $gain_canon_count (Gain of a canonical junction)<BR>
Unassigned = $unassigned_count (No change - unassigned to the categories above)<BR>
</P>
<BR>

<!-- Table1 -->
<P CLASS=\"Indented12LR_s16_Bold\">Reciprocal = $reciprocal_count (Loss of a canonical junction and corresponding gain of an alternative junction)</P>
<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[table1]\" data-openimage=\"Minus_icon.gif\" data-closedimage=\"Plus_icon.gif\"><img src=\"Plus_icon.gif\" border=\"0\" /></a></P>
<div id=\"table1\">
$table1_content
"</div><BR>

<!-- Table2 -->
<P CLASS=\"Indented12LR_s16_Bold\">Loss canonical = $loss_canon_count (Loss of a canonical junction)</P>
<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[table2]\" data-openimage=\"Minus_icon.gif\" data-closedimage=\"Plus_icon.gif\"><img src=\"Plus_icon.gif\" border=\"0\" /></a></P>
<div id=\"table2\">
$table2_content
"</div><BR>

<!-- Table3 -->
<P CLASS=\"Indented12LR_s16_Bold\">Gain alternative = $gain_alter_count (Gain of an alternative junction only)</P>
<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[table3]\" data-openimage=\"Minus_icon.gif\" data-closedimage=\"Plus_icon.gif\"><img src=\"Plus_icon.gif\" border=\"0\" /></a></P>
<div id=\"table3\">
$table3_content
"</div><BR>

<!-- Table4 -->
<P CLASS=\"Indented12LR_s16_Bold\">Gain canonical = $gain_canon_count (Gain of a canonical junction)</P>
<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[table4]\" data-openimage=\"Minus_icon.gif\" data-closedimage=\"Plus_icon.gif\"><img src=\"Plus_icon.gif\" border=\"0\" /></a></P>
<div id=\"table4\">
$table4_content
"</div><BR>

<!-- Table5 -->
<P CLASS=\"Indented12LR_s16_Bold\">Unassigned = $unassigned_count (No change - unassigned to the categories above)</P>
<P CLASS=\"Indented12LR\"><a href=\"#\" rel=\"toggle[table5]\" data-openimage=\"Minus_icon.gif\" data-closedimage=\"Plus_icon.gif\"><img src=\"Plus_icon.gif\" border=\"0\" /></a></P>
<div id=\"table5\">
$table5_content
"</div><BR>

</BODY>
</HMTL>

EOF






close(HTML);

exit();


###########################################################################################################################
#Import mutations from input file                                                                                         #
###########################################################################################################################
sub importMutations{
  my %mut;
  open (MUT, "$infile") || die "\n\nCould not open file: $infile\n\n";
  my $header = 1;
  my %columns;
  while(<MUT>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header ==1){
      my $c = 0;
      foreach my $head (@line){
        $columns{$head}{position} = $c;
        $c++;
      }
      $header = 0;
      next();
    }
    my $mid = $line[0];
    $mut{$mid}{pos} = $line[$columns{'pos'}{position}];
    $mut{$mid}{donor_acceptor} = $line[$columns{'donor_acceptor'}{position}];
    $mut{$mid}{gene_string} = $line[$columns{'gene_string'}{position}];
    $mut{$mid}{chr} = $line[$columns{'chr'}{position}];
    $mut{$mid}{strand} = $line[$columns{'strand'}{position}];
    $mut{$mid}{min_start} = $line[$columns{'min_start'}{position}];
    $mut{$mid}{max_end} = $line[$columns{'max_end'}{position}];
    $mut{$mid}{mutation_lib_list} = $line[$columns{'mutation_lib_list'}{position}];
    $mut{$mid}{mutation_lib_name_list} = $line[$columns{'mutation_lib_name_list'}{position}];
    $mut{$mid}{known_junction_ids_string} = $line[$columns{'known_junction_ids_string'}{position}];
    $mut{$mid}{alternate_junction_ids_string} = $line[$columns{'alternate_junction_ids_string'}{position}];
    $mut{$mid}{cum_median_diff} = $line[$columns{'cum_median_diff'}{position}];
    $mut{$mid}{avg_median_diff} = $line[$columns{'avg_median_diff'}{position}];
    $mut{$mid}{max_gain_canon} = $line[$columns{'max_gain_canon'}{position}];
    $mut{$mid}{max_loss_canon} = $line[$columns{'max_loss_canon'}{position}];
    $mut{$mid}{max_gain_alter} = $line[$columns{'max_gain_alter'}{position}];
    $mut{$mid}{max_loss_alter} = $line[$columns{'max_loss_alter'}{position}];
    $mut{$mid}{max_recip_diff} = $line[$columns{'max_recip_diff'}{position}];
    $mut{$mid}{class} = $line[$columns{'class'}{position}];
    if ($mut{$mid}{cum_median_diff} eq "NA"){$mut{$mid}{cum_median_diff} = 0;}
    if ($mut{$mid}{avg_median_diff} eq "NA"){$mut{$mid}{avg_median_diff} = 0;}
    if ($mut{$mid}{max_gain_canon} eq "NA"){$mut{$mid}{max_gain_canon} = 0;}
    if ($mut{$mid}{max_loss_canon} eq "NA"){$mut{$mid}{max_loss_canon} = 0;}
    if ($mut{$mid}{max_gain_alter} eq "NA"){$mut{$mid}{max_gain_alter} = 0;}
    if ($mut{$mid}{max_loss_alter} eq "NA"){$mut{$mid}{max_loss_alter} = 0;}
    if ($mut{$mid}{max_recip_diff} eq "NA"){$mut{$mid}{max_recip_diff} = 0;}

  }
  close(MUT);
  return(\%mut);
}


###########################################################################################################################
#Process each of the five classes and create an html table using the appropriate sort field                               #
###########################################################################################################################
sub createTable{
  my %args = @_;
  my $target_class = $args{'-class'};
  my $sort_column = $args{'-sort_column'};

  my $content = '';

  #Open table
  $content .= "\n<TABLE CLASS=\"Data\">";
  
  #Table header
  #print "Mutation\tGene\tLibrary List\tDonor/Acceptor\tCoords\tSort value ($sort_column)\tLink to SVG\n";
  $content .= "\n<TR> <TD CLASS=\"Head3\">Mutation</TD> <TD CLASS=\"Head3\">Gene</TD> <TD CLASS=\"Head3\">Library List</TD> <TD CLASS=\"Head3\">Donor/Acceptor</TD> <TD CLASS=\"Head3\">Coords</TD> <TD CLASS=\"Head3\">Sort value ($sort_column)</TD> <TD CLASS=\"Head3\">Link to Pretty SVG</TD> <TD CLASS=\"Head3\">Link to Active SVG</TD></TR>";

  foreach my $mid (sort {abs($mut_ref->{$b}->{$sort_column}) <=> abs($mut_ref->{$a}->{$sort_column})} keys %{$mut_ref}){
    my $sort_val = $mut_ref->{$mid}->{$sort_column};
    my $sort_val_f = sprintf("%.3f", $sort_val);
    my $class = $mut_ref->{$mid}->{class};
    unless ($class eq $target_class){next();}

    #Mutation ID, gene, library list, acceptor/donor, display coordinates, sort value, links to SVGs (both clickable and non-clickable versions)
    my $mutation_lib_name_list = $mut_ref->{$mid}->{mutation_lib_name_list};
    my $gene_string = $mut_ref->{$mid}->{gene_string};
    my $donor_acceptor = $mut_ref->{$mid}->{donor_acceptor};
    my $chr = $mut_ref->{$mid}->{chr};
    my $mut_pos = $mut_ref->{$mid}->{pos};
    my $strand = $mut_ref->{$mid}->{strand};
    my $min_start = $mut_ref->{$mid}->{min_start};
    my $max_end = $mut_ref->{$mid}->{max_end};
    my $display_coord = "$chr:$min_start-$max_end($strand)";
    my $svg_file_name1 = "$chr"."_"."$mut_pos".".svg";
    my $svg_file_name2 = "$chr"."_"."$mut_pos".".active.svg";
    my $svg_link1 = "<A HREF=\"$svg_file_name1\" TITLE=\"$mid Pretty SVG\">$mid</A>"; #for new window add: TARGET=\"_blank\" into A tag def
    my $svg_link2 = "<A HREF=\"$svg_file_name2\" TITLE=\"$mid Active SVG\">$mid</A>"; #for new window add: TARGET=\"_blank\" into A tag def

    #print "$mid\t$gene_string\t$mutation_lib_name_list\t$donor_acceptor\t$display_coord\t$sort_val_f\tSVG LINK\n";
    $content .= "<TR> <TD CLASS=\"Data1\">$mid</TD> <TD CLASS=\"Data1\">$gene_string</TD> <TD CLASS=\"Data1\">$mutation_lib_name_list</TD> <TD CLASS=\"Data1\">$donor_acceptor</TD> <TD CLASS=\"Data1\">$display_coord</TD> <TD CLASS=\"Data1\">$sort_val_f</TD> <TD CLASS=\"Data1\">$svg_link1</TD> <TD CLASS=\"Data1\">$svg_link2</TD></TR>";
  }
 
  #Close table
  $content .= "\n</TABLE>";
  

  return($content);
}





