#!/usr/bin/perl -w
#
# PROGRAM:  tagRepeats.perl
# PURPOSE:  to put consensus tags on any found ALU or any other 
#           kind of database sequence so that the consed user can
#           easily see whether a region is in such a repeat
#
# HOW IT WORKS:  It extracts the consensus from the ace file and calls 
#           cross_match, and parses the output.  
# INPUTS:   The ace file.  You must also provide a fasta file of the 
#           ALU or any other sequence you want to tag
# REVISION:  Nov 1998 (David Gordon)


$szVersion = "000727";
if ( $#ARGV >= 0 ) {
  if ( $ARGV[0] eq "-V" || $ARGV[0] eq "-v" ) {
    print "$szVersion\n";
    exit( 1 );
  }
}


$szUsage = "Usage: tagRepeats.perl (ace file)";

defined( $szConsedHome = $ENV{'CONSED_HOME'} ) ||
    ( $szConsedHome = "/usr/local/" );


$szCrossMatchExe = $szConsedHome . "/bin/cross_match";
$szRepeatsFile = $szConsedHome . "/lib/screenLibs/repeats.fasta";

die "$szUsage" if ( $#ARGV != 0 );

$szAceFile = $ARGV[0];


die "Can't read $szAceFile" if ( ! -r $szAceFile );
die "Can't execute $szCrossMatchExe" if ( ! -x $szCrossMatchExe );


# strip off the .fasta.screen.*
($szRootName = $szAceFile) =~ s/\.fasta\.screen\..*$//;

$szContigsFile = $szRootName . ".contigs";

($nFormatOfAceFile, %aContigsInAceFile ) = 
  &bDetermineAceFileFormat( $szAceFile );



$szPadsTable = $szRootName . ".pads";
if ( $nFormatOfAceFile == 2) {
  &constructPadsTable( $szAceFile, $szPadsTable );
}
else {
  &constructPadsTableAceFormat1( $szAceFile, $szPadsTable );
}

&buildPaddedFromUnpaddedTable;

if ( $nFormatOfAceFile == 2 ) {
  &constructContigsFile( $szAceFile, $szContigsFile );
}
else {
  &constructContigsFileAceFormat1( $szAceFile, $szContigsFile,
                                   %aContigsInAceFile );
}


$szCrossMatchOutputFilename =
  $szRootName . "_to_alu.cross";

$szCommand = "$szCrossMatchExe $szContigsFile $szRepeatsFile -tags -minmatch 10 >$szCrossMatchOutputFilename";
#$szCommand = "$szCrossMatchExe $szContigsFile $szRepeatsFile -tags >$szCrossMatchOutputFilename";

# not checking output because for some reason !system doesn't work
# on old versions of cross_match or perl
system( $szCommand );

$szNewTagsFile = $szRootName . ".newtags";
open( filNewTags, ">$szNewTagsFile" ) || die "can't open $szNewTagsFile for output";

&readAndParseCrossMatchOutput();

close( filNewTags );

print "about to add tags $szNewTagsFile to $szAceFile...\n";

open( filAce, ">>$szAceFile" ) || die "couldn't open $szAceFile for append";
open( filNewTags, "$szNewTagsFile" ) || die "couldn't open $szNewTagsFile for append";

while( <filNewTags> ) {
    print( filAce  $_ );
}

close( filAce );
close( filNewTags );

print "done tagging repeats\n";
exit( 0 );








sub bDetermineAceFileFormat {
    ( $szAceFile ) = @_;

    my %aContigs = ();

    open( filAce, "$szAceFile" ) || die "couldn't open $szAceFile for reading";
    $szFirstLine = <filAce> || die "0 length file $szAceFile";

    if ( length( $szFirstLine ) > 3 ) {
        if (substr( $szFirstLine, 0, 3 ) eq "AS " ) {
            return( (2, %aContigs ) );
        }
    }
    
    seek( filAce, 0, 0 );  # reposition to beginning of file

    my @aListOfContigsAndReads = ();
    my %aListOfReads = ();

    while( <filAce> ) {
        if ( length($_) >= 4 ) {
            if ( substr( $_, 0, 4 ) eq "DNA " ) {
                @aWords = split;
                push( @aListOfContigsAndReads, $aWords[1] );
            }
            elsif( length( $_ ) >= 16 ) {
                if ( substr( $_, 0, 16 ) eq "Assembled_from* " ) {
                    @aWords = split;
                    $aListOfReads{$aWords[1] } = "";
                }
            }
        }
    }

    foreach $szContigOrRead (@aListOfContigsAndReads ) {
        if ( !exists( $aListOfReads{ $szContigOrRead } ) ) {
            $aContigs{ $szContigOrRead } = "";
        }
    }

    return( (1, %aContigs ) );
}



sub constructContigsFile {
    ( $szAceFile, $szContigsFile ) = @_;

    open( filAce, "$szAceFile" ) || die "couldn't open $szAceFile for reading";
    open( filContigs, ">$szContigsFile" ) || 
        die "couldn't open $szContigsFile for writing";

    while( <filAce>) {
        if ( length( $_ ) > 3 ) {
            if ( substr( $_, 0, 3 ) eq "CO " ) {
                # found a CO line
                @aWords = split;
                $szContigName = $aWords[1];

                print( filContigs ">$szContigName\n" );
            
                while( <filAce>) {
                    last if ( length( $_ ) == 1 );

                    # do an efficient check in case someone
                    # has edited the ace file and put 
                    # some whitespace on the line
                    if ( (substr( $_, 0, 1 ) eq " ") ||
                        (substr( $_, 0, 1 ) eq "\t" )) {
                        last if ( $_ =~ /^\s*$/ );
                    }
                    # remove any pads
                    s/\*//g;
                    print( filContigs $_ );
                }
                
                # end of that contig
                print( filContigs "\n" );
            }
        }
    }
    close( filAce );
    close( filContigs );
}


 


sub constructContigsFileAceFormat1 {
    ( $szAceFile, $szContigsFile, %aContigs ) = @_;

    open( filAce, "$szAceFile" ) || die "couldn't open $szAceFile for reading";
    open( filContigs, ">$szContigsFile" ) || 
        die "couldn't open $szContigsFile for writing";

    
    while( <filAce>) {
        if ( length( $_ ) > 4 ) {
            if ( substr( $_, 0, 4 ) eq "DNA " ) {
                # found a DNA line.  See if it is for a contig
                @aWords = split;
                $szContigName = $aWords[1];
                if ( !exists( $aContigs{ $szContigName } ) ) {
                    next;
                }

                print( filContigs ">$szContigName\n" );
            
                while( <filAce>) {
                    last if ( length( $_ ) == 1 );

                    # do an efficient check in case someone
                    # has edited the ace file and put 
                    # some whitespace on the last line
                    if ( (substr( $_, 0, 1 ) eq " ") ||
                        (substr( $_, 0, 1 ) eq "\t" )) {
                        last if ( $_ =~ /^\s*$/ );
                    }
                    # remove any pads
                    s/\*//g;
                    print( filContigs $_ );
                }
                
                # end of that contig
                print( filContigs "\n" );
            }
        }
    }
    close( filAce );
    close( filContigs );
}



sub readAndParseCrossMatchOutput {
  open( filCrossMatchOutput, $szCrossMatchOutputFilename ) ||
    die "couldn't open cross match output file $szCrossMatchOutputFilename";

  while( <filCrossMatchOutput> ) {
    chomp;
    if ( $_ =~ /^ALIGNMENT/ ) {
      &parseAlignmentLine;
    }
  }

  close( filCrossMatchOutput );
}




sub constructPadsTable {
    ( $szAceFile, $szPadsTable ) = @_;
# INPUT:  consed new_ace file 
# OUTPUT: table of pad positions


    
    open( filAce, "$szAceFile" ) || die "couldn't open $szAceFile for reading";
    open( filPads, ">$szPadsTable" ) || die "couldn't open $szPadsTable for writing";

    while( <filAce>) {
        if ( length( $_ ) > 3 ) {
            if ( substr( $_, 0, 3 ) eq "CO " ) {
                # found a CO line
                @aWords = split;
                $szContigName = $aWords[1];

                print( filPads  ">$szContigName" );
                
                $nPaddedPos = 0;
                $nNumberOfPads = 0;
                while( <filAce>) {
                    last if ( length( $_ ) == 1 );
                    chop;

                    # do an efficient check in case someone
                    # has edited the ace file and put 
                    # some whitespace on the line
                    if ( (substr( $_, 0, 1 ) eq " ") ||
                        (substr( $_, 0, 1 ) eq "\t" )) {
                        last if ( $_ =~ /^\s*$/ );
                    }
                    for( $nPosOnLine = 0; $nPosOnLine < length( $_ ); ++$nPosOnLine ) {
                        ++$nPaddedPos;
                        if ( substr( $_, $nPosOnLine, 1 ) eq "*" ) {
                            ++$nNumberOfPads;
                            if ( ($nNumberOfPads % 10) == 1 ) {
                                print( filPads  "\n$nPaddedPos" );
                            }
                            else {
                                print( filPads  " $nPaddedPos" );
                            }
                        }
                    }
                }

                # end of that contig
                print( filPads  "\n\n" );
            }
        }
    }
    close( filAce );
    close( filPads );
}


    
sub constructPadsTableAceFormat1 {
    ( $szAceFile, $szPadsTable, %aContigs ) = @_;
# INPUT:  consed new_ace file 
# OUTPUT: table of pad positions

    open( filAce, "$szAceFile" ) || die "couldn't open $szAceFile for reading";
    open( filPads, ">$szPadsTable" ) || die "couldn't open $szPadsTable for writing";

    while( <filAce>) {
        if ( length( $_ ) > 4 ) {
            if ( substr( $_, 0, 4 ) eq "DNA " ) {
                # found a DNA line.  Check if it is for a contig.
                @aWords = split;
                $szContigName = $aWords[1];

                if ( !exists( $aContigs{ $szContigName } )) {
                    next;
                }

                print( filPads  ">$szContigName" );
                
                $nPaddedPos = 0;
                $nNumberOfPads = 0;
                while( <filAce>) {
                    last if ( length( $_ ) == 1 );
                    chop;

                    # do an efficient check in case someone
                    # has edited the ace file and put 
                    # some whitespace on the last line
                    if ( (substr( $_, 0, 1 ) eq " ") ||
                        (substr( $_, 0, 1 ) eq "\t" )) {
                        last if ( $_ =~ /^\s*$/ );
                    }
                    for( $nPosOnLine = 0; $nPosOnLine < length( $_ ); ++$nPosOnLine ) {
                        ++$nPaddedPos;
                        if ( substr( $_, $nPosOnLine, 1 ) eq "*" ) {
                            ++$nNumberOfPads;
                            if ( ($nNumberOfPads % 10) == 1 ) {
                                print( filPads  "\n$nPaddedPos" );
                            }
                            else {
                                print( filPads  " $nPaddedPos" );
                            }
                        }
                    }
                }

                # end of that contig
                print( filPads  "\n\n" );
            }
        }
    }
    close( filAce );
    close( filPads );
}






sub parseAlignmentLine {

    
# if reached here, the first token is ALIGNMENT
# Typical examples:
# uncomplemented:
# ALIGNMENT   772  0.49 0.37 0.61  Contig2        5   822 (0)    NewContig2      403 1218 (0)  
# complemented:
# ALIGNMENT   382  0.25 0.76 0.00  Contig1        1   397 (248)  C NewContig2   (818)   400     1  


  @aWords = split;

  $szContigName = $aWords[5];
  $nContigAlignLeft = $aWords[6];
  $nContigAlignRight = $aWords[7];
  $cCompFlag = $aWords[9];
  if ( $cCompFlag eq "C" ) {

    $szThingToTagName = $aWords[10];
    $nThingToTagAlignLeft = $aWords[12];
    $nThingToTagAlignRight = $aWords[13];
  }
  else {
    $cCompFlag = "U";
    $szThingToTagName = $aWords[9];
    $nThingToTagAlignLeft = $aWords[10];
    $nThingToTagAlignRight = $aWords[11];
  }

  $nPaddedLeft = &nPaddedFromUnpadded( $szContigName, $nContigAlignLeft );
  $nPaddedRight = &nPaddedFromUnpadded( $szContigName, $nContigAlignRight );

  
  $szDate = &szGetDateForTag();

  print filNewTags "\nCT{\n";
  $szMainTagLine = "$szContigName repeat tagRepeats.perl $nPaddedLeft $nPaddedRight $szDate NoTrans\n";
  print filNewTags $szMainTagLine;
  print filNewTags "$szThingToTagName\n";
  print filNewTags "}\n\n";
}


sub szGetDateForTag {
  my $szDate;
  ($nSecond, $nMinute, $nHour, $nDayInMonth, $nMonth, $nYear, $wday, $yday, $isdst ) = localtime;

  undef $isdst;
  undef $wday;
  undef $yday;
  
  if ( $nYear >= 100 ) {
    $nYear = $nYear % 100;
  }

  $szDate = sprintf( "%02d%02d%02d:%02d%02d%02d",
           $nYear,
           $nMonth + 1,
           $nDayInMonth,
           $nHour,
           $nMinute,
           $nSecond );

  return( $szDate );
}



sub buildPaddedFromUnpaddedTable {
# now deal with pads in new assembly



    %paddedFromUnpadded = ();


    open( filPadsTable, "$szPadsTable" ) ||
        die "Couldn't open new pads table $szPadsTable";

    while( <filPadsTable> ) {
        chop;
        if ( length( $_ ) > 1 ) {
            if ( substr( $_, 0, 1 ) eq ">" ) {
                # found a start of a contig
                # skip over the ">" and take the contig name                
                $szContig = substr( $_, 1 );
                
                $refArray = $paddedFromUnpadded{ $szContig } = [];
                
                while( <filPadsTable> ) {
                    chop;
                    last if ( length( $_ ) == 0 );
                    @aPadPositions = split;
                    push( @$refArray, @aPadPositions );
                }

                # now we have finished collecting the pads for this contig
                # convert to the unpadded array by subtracting 
                # [0] subtract 1
                # [1] subtract 2
                # etc.
                # Thus we will have the upper bound in unpadded positions 
                # that a particular range is useful for.  Alternatively, this
                # gives the unpadded position base just before each pad.  Hence
                # to find the padded position, find the largest element such 
                # that the unpadded pos is <= the value of the element.  
                # Then add the number of pads before that element.  This will
                # be the padded position.  You can easily find the number of
                # pads before the unpadded pos since it will be the index of the 
                # element.

                for( $nPos = 0; $nPos <= $#$refArray; ++$nPos ) {
                    $$refArray[ $nPos ] -= ($nPos + 1 );
                }
            }
        }
    }

    close( filPadsTable );
    unlink( $szPadsTable );
}




sub nPaddedFromUnpadded {
    die "usage: &nPaddedFromUnpadded( <contigname>, <unpadded position> );" 
        if ( $#_ != 1 );

    $szContig = $_[0];
    $nUnpaddedPos = $_[1];

    if (!exists $paddedFromUnpadded{ $szContig } ) {
        print "Contig $szContig doesn't exist";
        return -10000;
    }
    
    # now convert to padded:
    
    $refArray = $paddedFromUnpadded{ $szContig };
    if ( $#$refArray == -1 ) {
        # case in which there are no pads
        $nNumberOfPadsBeforeBase = 0;
    }
    elsif ( $nUnpaddedPos <= $$refArray[ 0 ] ) {
        # base is before the first pad
        $nNumberOfPadsBeforeBase = 0;
    }
    else {
        $bFoundUnpaddedMax = 0;
        for( $nRange = 0; ($nRange <= $#$refArray) && !$bFoundUnpaddedMax;
            ++$nRange ) {
        
            $nUnpaddedTopOfRange = $$refArray[ $nRange ];
            if ( $nUnpaddedPos <= $nUnpaddedTopOfRange ) {
                $bFoundUnpaddedMax = 1;
                $nNumberOfPadsBeforeBase = $nRange;
                                #  the unpadded pos goes in the
                                #  previous range
            }
        }

        if (!$bFoundUnpaddedMax ) {
            # This is the case in which the unpadded pos is 
            # after the last pad.  So convert to padded just
            # by adding the total number of pads.

            $nNumberOfPadsBeforeBase = $#$refArray + 1;
                                # $#$refArray is the subscript of 
                                # the last pad, and the subscript of 
                                # the first pad is 0, hence the number
                                # of pads is $#$refArray + 1
        }
    }

    
    $nPaddedPos = $nUnpaddedPos + $nNumberOfPadsBeforeBase;

    return $nPaddedPos;
}
            
    
