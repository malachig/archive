#!/usr/bin/perl


#
#This program extracts all the junctions from a spidey output
#
#By: Brian Cho, BCGSC
#Supervisor: Malachi Griffith
#



#Convert qseq file into a fasta file (.fa):
#cat RA0001.qseq.txt | cut -f 9 | perl -ne 'unless ($_ =~ /\./){$c++; print ">$c\n$_"}' > RA0001.fa


#Example spidey command:
#/home/bcho/tools/spidey.linux -i RatH-Genomic.fa -m RatH-RefcDNA.fa -p 1

#To pipe spidey output to junction.pl (this perl script):
#/home/bcho/tools/spidey.linux -i RatH-Genomic.fa -m RatH-RefcDNA.fa -p 1 | /home/bcho/Perl/spidey/junction.pl --size=10

#command line script
#/home/bcho/tools/spidey.linux -i RatH-Genomic_RN4_plus.fa -m test2.fa -p 1 | /home/bcho/Perl/spidey/junction.pl --size=10 | /home/bcho/Perl/spidey/convertCoords.pl --strand='-' --chr_start=14606372 --chr_end=14694051

#merge files and send to makeGFF.pl
#/home/bcho/Perl/spidey/merge.pl --dir=/home/bcho/Cav32_Project/Rat_Illumina_Data/cluster_test/RA0001/spidey/ | /home/bcho/Perl/spidey/makeGFF.pl --strand='-' --chr=10 --min_count=10 --lib_id=RA0001 --color=255,0,0

use strict;
use Getopt::Long;
use Data::Dumper;
use Term::ANSIColor qw(:constants);

my @junctionID;
my @genomicStart;
my @genomicEnd;
my @donor;
my @acceptor;
my @unit1_size;
my @unit2_size;


my $n = 0;
my $totalExons = 0;
my $strand;
my $min_align_size;

my %junctions;




GetOptions( 'size=i' => \$min_align_size);
unless ($min_align_size =~ /\d+/)
{
  print MAGENTA, "Please specify the minimum alignment size of exon using the following format:\n";
  print "--size=10\n";
  print "junction.pl terminated\n\n", RESET;
  exit();
}



while((my $line = <STDIN>))
{
  if ($line =~ /^\-\-SPIDEY/ || eof)
  {
    #Reached end of block of data from spidey (corresponds to one read/cDNA mapped)
      #The first iteration of the "--SPIDEY" line will not have any data stored in the arrays --> not a big deal
    #Stores all arrays into %junctions hash and then resets all arrays
    for (my $m = 0; $m < $totalExons-1; $m++)
    {
      if ($unit1_size[$m] >= $min_align_size && $unit2_size[$m] >= $min_align_size)
      {
      

        if ($junctions{"$genomicStart[$m]\:$genomicEnd[$m]"})
        {
          #If the junction has been defined already, the count (number of observations) is incremented
          $junctions{"$genomicStart[$m]\:$genomicEnd[$m]"}{count}++; #print " previously found ";
        }

        else
        {
          #If this is the first time that the junction is found, it is defined in the hash
          $junctions{"$genomicStart[$m]\:$genomicEnd[$m]"}{junctionID} = $junctionID[$m]; #print " first time ";
          $junctions{"$genomicStart[$m]\:$genomicEnd[$m]"}{start} = $genomicStart[$m];
          $junctions{"$genomicStart[$m]\:$genomicEnd[$m]"}{end} = $genomicEnd[$m];
          $junctions{"$genomicStart[$m]\:$genomicEnd[$m]"}{donor} = $donor[$m];
          $junctions{"$genomicStart[$m]\:$genomicEnd[$m]"}{acceptor} = $acceptor[$m];
          $junctions{"$genomicStart[$m]\:$genomicEnd[$m]"}{strand} = $strand;
          $junctions{"$genomicStart[$m]\:$genomicEnd[$m]"}{count} = 1;
        }
      }
    }

    #Resetting all arrays
    @junctionID = 0;
    @genomicStart = 0;
    @genomicEnd = 0;
    @donor = 0;
    @acceptor = 0;
    @unit1_size = 0;
    @unit2_size = 0;
    

    $n = 0;
    $totalExons = 0;
    $strand = 0;

  }

  else
  {
    if ($line =~ /^Strand\:\s(\w)/)
    {
      $strand = $1;
    }

    if ($line =~ /^Number\sof\sexons\:\s(\d+)/)
    {
      $totalExons = $1;
    }
   
    #Case for positive (+) strand
    if ($strand eq "p")
    {
      
      #Test for data line
      if ($line =~ /^Exon/)
      {
        
        if ($totalExons == 1) {next();}

        $junctionID[$n] = $n+1;
        
        #Extract info
        unless ($n == $totalExons-1) #Skips the last line when extracting junction start coord, donor canonicalness, and unit1_size
        {
          if ($line =~ /\-(\d+)\s\(gen\)/)
          {
            $genomicStart[$n] = $1;
          }

          if ($line =~ /\(d\s+a\)\:\s(\d)\s+/)
          {
            $donor[$n] = $1;
          }

          if ($line =~ /Exon\s\d+\:\s(\d+)\-(\d+)\s\(gen\)/)
          {
            $unit1_size[$n] = $2-$1;
          }
        }


        unless ($n == 0) #Skips the first line when extracting junction end coord, acceptor canonicalness, and unit2_size
        {
          if ($line =~ /Exon\s\d+\:\s(\d+)\-/)
          {
            $genomicEnd[$n-1] = $1;
          }

          if ($line =~  /\(d\s+a\)\:\s\d\s+(\d)/)
          {
            $acceptor[$n-1] = $1
          }

          if ($line =~ /Exon\s\d+\:\s(\d+)\-(\d+)\s\(gen\)/)
          {
            $unit2_size[$n-1] = $2-$1;
          }

        }

        $n++;
  
      }
    }
    
    #Case for minus (-) strand
    if ($strand eq "m")
    {
      
      #Test for data line
      if ($line =~ /^Exon/)
      {

        if ($totalExons == 1) {next();};

        $junctionID[$n] = $n+1;
      
        #Extract info
        unless ($n == 0) #Skips the first line when extracting junction start coord, acceptor canonicalness, 
        {
          if ($line =~ /\:\s\d+\-(\d+)\s\(gen\)/)
          {
            $genomicStart[$n-1] = $1; #print "$junctionID[$n-1] start: $genomicStart[$n-1]\n";
          }

          if ($line =~ /\(d\s+a\)\:\s\d\s+(\d)/)
          {
            $acceptor[$n-1] = $1;
          }

          if ($line =~ /Exon\s\d+\(\-\)\:\s(\d+)\-(\d+)\s/)
          {
            $unit2_size[$n-1] = $2-$1;
          }

          
        }


        unless ($n == $totalExons-1) #Skips the last line when extracting junction end coord, donor canonicalness, 
        {
          if ($line =~ /Exon\s\d+\(\-\)\:\s(\d+)\-/)
          {
            $genomicEnd[$n] = $1; #print "$junctionID[$n] end: $genomicEnd[$n]\n";
          }

          if ($line =~ /\(d\s+a\)\:\s(\d)\s+/)
          {
            $donor[$n] = $1;
          }

          if ($line =~ /Exon\s\d+\(\-\)\:\s(\d+)\-(\d+)\s/)
          {
            $unit1_size[$n] = $2-$1;
          }


        }

        $n++;
  
      }
    }


  }

}

#print Dumper %junctions;
foreach my $j (sort {$junctions{$a}->{start} <=> $junctions{$b}->{start}} keys %junctions)
{
  print "$j\t$junctions{$j}{start}\t$junctions{$j}{end}\t$junctions{$j}{count}\t$junctions{$j}{donor}\t$junctions{$j}{acceptor}\n";
}

#Output
#open (OUTPUT, ">output.txt");
#print OUTPUT Dumper(%junctions);
#close (OUTPUT);
#print "\nsaved as output.txt\n";

exit();


