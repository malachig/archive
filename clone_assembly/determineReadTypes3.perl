#!/usr/bin/perl -w
#  deteremineReadTypes.perl
#  
#  Purpose:  to write into the phd file information about the template 
#            name and forward/reverse information in order to 
#            help phrap and consed/autofinish.   This came about when
#            we found that we couldn't force labs to change their 
#            read naming convention to that expected by phrap.  So we
#            allow you to have any consistent naming convention you want,
#            but you have to do a little perl programming to customize
#            this script to fit your naming convention.
#
#  How to use it:  
#     determineReadTypes.perl -justThisPhdFile <name of phd file with ../phd_dir/>
#        This will just examine and add primer and template WR items
#        to the phd file specified, unless primer and template WR items
#        are already in this phd file
#
#     determineReadTypes.perl (no arguments)
#        This will examine every phd file and add to the ones that do not
#        already have primer and template WR items
#        
#
#  What it writes into the phd file:
#            whole read info items of the following types:
#
#  How to modify it:
#     First of all, you must know a little perl, especially regular 
#     expressions..  If you don't know perl, 
#     don't dispair--it will just take you a couple of hours to learn it.
#     I suggest Randal Schwartz, Learning Perl.
#
#     There are 2 models of interpreting read names:  
#          findTemplateNameAndTemplateTypeStLouisFormat
#          and
#          findTemplateNameAndTemplateTypeSeattleFormat
#     You can start with either one of these (found below), 
#     whichever is closest
#     to your naming convention, and modify it.
#     
#     As you are debugging your script, run it like this:
#     determineReadTypes.perl -justView
#     That will not print anything out--it will just tell you what 
#     it thinks each read is!
#
#  WR{
#  template phredPhrap 990224:045110
#  name: (template name)  
#  type: (bac, cosmid, puc, pbc, pcr, etc)
#  size: (template size, especially useful for PCR products)
#  }
# 
#  each of the above fields would be optional
#  Currently, this program will typically just put the "name:" line
#  into the phd file.  In the case of bac or cosmid sequencing reads,
#  it will write the "type:" line with "bac" or "cos" on the line.
#
#
#  WR{
#  primer phredPhrap 990224:045110
#  type: univ fwd
#  seq: (primer sequence)              (this line is optional)
#  }
# 
#  In the above, "univ fwd" could be replaced by:
#
#  univ rev
#  PCR end
#  walk
#
#  Currently this program does not put in the "seq:" line.
#  If the read is with a custom primer and is not the end of a pcr product,
#  then this program puts "walk" into the "type:" line.
#  
#  If this program determines that the read is a universal primer read,
#  then it tries to determine whether it is a universal primer forward or 
#  universal primer reverse read. 
#
#  Note:  It is ok to run this program over and over again on the same 
#  project since it will not add anything further to a phd file if it finds 
#  it has already added template and primer WR items to that phd file.

###############################################################
#
# main program start
#
###############################################################




$szConsedVersion = "7.45 (beta) (990406)";
$phdDirPath = "../phd_dir";


# the -justView option is used to help debug your modifications of this
# script.  Type:
# determineReadTypes.perl -justView
# and no changes will take place, but this script will tell you
# what it thinks about each read.

$bJustView = 0;
if ( $#ARGV >= 0 ) {
   if ( $ARGV[0] eq "-V" || $ARGV[0] eq "-v" ) {
      print "version: $szConsedVersion\n";
      exit( 1 );
   }
   elsif ( $ARGV[0] eq "-justView" ) {
     $bJustView = 1;
     shift @ARGV;
   }
}


if ( $#ARGV == -1 ) {
  # case of no arguments so examine all phd files
  
  $nNumberOfFilesProcessed = 0;
  $nNumberOfFilesExamined = 0;
  opendir( dirPhdDir, $phdDirPath ) || die "could not even open $phdDirPath";
  
  while( defined( $szPhdFile = readdir( dirPhdDir ) ) ) {
    if ( index( $szPhdFile, ".phd." ) >= 0 ) {

      $bIsPhdFileAlreadyProcessedd = &bIsPhdFileAlreadyProcessed( $szPhdFile );

      if ( ( $bIsPhdFileAlreadyProcessedd == 0) || $bJustView  ) {
        &processOnePhdFile( $szPhdFile );
        ++$nNumberOfFilesProcessed;
      }
      ++$nNumberOfFilesExamined;

      if ( $nNumberOfFilesExamined % 100 == 0 ) {
        
        print "files examined: $nNumberOfFilesExamined, files processed: $nNumberOfFilesProcessed\n";
      }
    }
  }
  
  closedir( dirPhdDir );
}
elsif ( $ARGV[0] eq "-justThisPhdFile" ) {
  die "the phd filename must follow -justThisPhdFile" if ( $#ARGV < 1 );

  $szFullPath = $ARGV[1];

  die "couldn't open $szFullPath for reading" if (!-r $szFullPath );
  die "couldn't open $szFullPath for writing" if (!-w $szFullPath );

  # cut off the path

  use File::Basename;

  my $szPhdFile = basename( $szFullPath );

  $bIsPhdFileAlreadyProcessedd = &bIsPhdFileAlreadyProcessed( $szPhdFile );
  
  if ( ( $bIsPhdFileAlreadyProcessedd == 0 ) || $bJustView ) {
    &processOnePhdFile( $szPhdFile );
    print "$szPhdFile processed\n";
  }
  else {
    print "$szPhdFile already processed\n";
  }
}
else {
  die "unrecognized arguments";
}


exit(0);

## end of program

# S. Zuyderduyn (31 Oct 2000; 1 Dec 2000)
# Y. Butterfield (30 Aug 2001)
# R. Warren (25 Apr 2003)
sub findTemplateNameAndTemplateTypeVancouverFormat {
	die "usage: &findTemplateNameAndTemplateTypeVancouverFormat( <phd file name> );"
	if( $#_ != 0 );

	my $phdFileName = $_[0];

	my $readName;
	( $readName = $phdFileName ) =~ s/\.phd\.[0-9]+$//;

	#$firstDot = index( $readName, "." );
	#if( $firstDot == -1 ) {
	#	print STDERR "Can't figure out template and primer type for $phdFileName since no dot in filename\n";
	#	return(  (1) );
	#}

	#if( $firstDot == length( $readName ) ) {
	#	print STDERR "Can't figure out template and primer type for $phdFileName since no letter after the dot\n";
	#	return( (1) );
	#}

	#my $nameUpToDot = substr( $readName, 0, $firstDot );

	my $templateName;
	my $chemTag;
	my $library;

	my $bNewFormat = 0;

	my $templateType = "pcr";
	my $primerType;

	if( ( $readName =~ /^(\d)EST_(\d+)__/ )|( $readName =~ /^EXTERNAL__(\d?)EST_(\d+)__/ ) ) {
	    $templateName=$2;
	    my $direction=$1;
	    if( $direction == 3 ) { $primerType = "univ rev"; }
	    elsif( $direction == 5 ) { $primerType = "univ fwd"; }
	    else { $primerType = "walk"; }
	    #else { print "Can't figure out primerType. ( $readName )\n"; die "\n"; }
	    $bNewFormat = 1;
	  } elsif( ( $readName =~ /^TRANS__(\w+)_/ )|( $readName =~ /^EXTERNAL__TRANS__(.*)\.?\d?\d?$/ ) ) {
	    $templateName=$1;
	    $primerType="walk";
	    $templateType="bac";
	    $bNewFormat = 1;
	  } elsif( $readName =~ /^\dOLIGO_(\d+)__/ ) {
	    $templateName=$1;
	    $templateType="bac";
	    $primerType = "walk";
	    $bNewFormat = 1;
	  } elsif( $readName =~ /^OLIGO_(\w{5})(\w+)_(\w+)_\d+__\w+_.+/ ) {
	    # file is of format like OLIGO_MGC0115_A11_2__MGC0122G02_CB.1 - used in MGC closure project
	    $templateName="$1$2$3";;
	    $templateType="bac";
	    $primerType = "walk";
	    $bNewFormat = 1;
	  } elsif( $readName =~ /^PolyT_(\d+)__/ ) {
	    $templateName=$1;
	    $primerType = "univ rev";
	    $bNewFormat = 1;
	  } elsif( $phdFileName =~ /^(\w{5})(\w+\w{3})\.(\w+)\.\d?\d?\.?abd\.phd\.\d+$/ ) {
	    # file is of format like TL05A8dF03.EA.abd.phd.1 or TL05A8dF03.EA.1.abd.phd.1
	    $library = $1;
	    $templateName = "$1$2";
	    $chemTag = $3;
	  } elsif( $phdFileName =~ /^(\w{5})(\w+)\.(\w+)\.?\d?\d?_(\w+)___\d+\.ab1\.phd\.\d+$/ ) {
	    # file is of format like TL05A8d.BA_F03___031.ab1.phd.1 or TL05A8d.BA.1_F03___031.ab1.phd.1
	    print "$1\t$2\t$3\t$4\n";
	    $library = $1;
	    $templateName = "$1$2$4";
	    $chemTag = $3;
	  } elsif( $phdFileName =~ /^(\w{5})(\w+)\.(\w+)\.?\d?\d?_(\w+)___x+\.ab1\.phd\.\d+$/ ) {
	    # file is of format like TL05A8d.BA_F03___031.ab1.phd.1 or TL05A8d.BA.1_F03___031.ab1.phd.1
	    $library = $1;
	    $templateName = "$1$2$4";
	    $chemTag = $3;
	  } elsif( $phdFileName =~ /^(\w{5})(\w+)\.(\w+)\.?\d?\d?_(\w+)\.phd\.\d+$/ ) {####standard adopted for trace names without capillary #
            # file is of format like TL05A8d.BA_F03.phd.1 or TL05A8d.BA.1_F03.phd.1
            $library = $1;
            $templateName = "$1$2$4";
            $chemTag = $3;
            print "$library $templateName $chemTag\n";
	  }    elsif( $phdFileName =~ /^([XN]M_)(\d+)(_)(ORF).*\.phd\.\d+$/ ) {#### special case for fake RefSeq ORF phd files #
            # file is of format like Hs_25351_amplicon.phd.1
            $library = $1;
            $templateName = "$1$2$4";
	    $bNewFormat = 1;
	    $primerType = "walk";
	  }elsif( $phdFileName =~ /^HIV_ref_ORF\.phd\.\d+$/ ) {#### special case for fake RefSeq ORF phd files #
            # file is of format like Hs_25351_amplicon.phd.1
            $library = "HIV01";
            $templateName = "HIV01";
	    $bNewFormat = 1;
	    $primerType = "walk";
	  } elsif($phdFileName =~ /^(\w{5})(\w+)\.(\w+)\.?\d?\d?_(\w+)_-_(\d+)\.ab1\.phd\.\d+$/){
	    #Added by MG to deal with reads MAPPED AS 384 WELL FORMAT for MGC Closure project
	    #MGC051.B21.2_A01_-_015.ab1.phd.1
	    #print "\nDEBUG: 1 = $1, 2 = $2, 3 = $3, 4 = $4, 5 = $5\n";
	    $library = $1;
	    $templateName = "$1$2$4";
	    $chemTag = $3;
	  } elsif ($phdFileName =~ /^(\w{5})(\w+)\.(\w+)\.?\d?\d?_-_(\w+)_(\d+)\.ab1\.phd\.\d+$/){
	    #Added by MG to deal with 3100 MACHINE reads MAPPED AS 384 WELL FORMAT! for MGC Closure project
	    #MGC051.BR.2_-_A01_01.ab1.phd.1
	    #print "\nDEBUG 3100 MACHINE: 1 = $1, 2 = $2, 3 = $3, 4 = $4, 5 = $5\n";
	    $library = $1;
	    $templateName = "$1$2$4";
	    $chemTag = $3;

	  } elsif ($phdFileName =~ /^(\w{5})(\w+)\.(\w+)\.?\d?\d?_(\w+)\.ab1\.phd\.\d+$/){
	    #Added by MG to deal with new read names which lack the capillary id value
	    #AVF032.B3.2_A01.ab1.phd.1
	    #print "\nDEBUG NEW GENERIC READ NAME: 1 = $1, 2 = $2, 3 = $3, 4 = $4, 5 = $5\n";
	    $library = $1;
	    $templateName = "$1$2$4";
	    $chemTag = $3;

	  } elsif ($phdFileName =~ /^.*_.*_(\w{5})(\w+)\.(\w+)\.?\d?\d?_(\w+)\.phd\.\d+$/){
	    #Added by MG to deal with AVF read names
	    #HA._01_AVF021.B3_A03.phd.1
	    #print "\nDEBUG NEW GENERIC READ NAME: 1 = $1, 2 = $2, 3 = $3, 4 = $4, 5 = $5\n";
	    $library = $1;
	    $templateName = "$1$2$4";
	    $chemTag = $3;

	  }elsif($phdFileName =~ /^(\w+)\.phd\.\d+/){
            print "THIS APPEARS TO BE A FAKE READ: $phdFileName\n";
            $library = "FAKE";
            $templateName = $1;
            $chemTag = "CB";
          } else {
	    print "Can't figure out file format!  ---> $phdFileName\n";
	    die "\n";
	  }
	
	if( $bNewFormat == 0 ) {
		$chemTag =~ /^(\w)(\w+)/;

		my $sub_chemtype = $1;
		my $priming_site = $2;


		if( $priming_site eq "7" ) {
		  $primerType = "univ fwd";
		} elsif( $priming_site eq "21" ) {
		  $primerType = "univ fwd";
		} elsif( $priming_site =~ /^A/ ) {
		  $primerType = "walk";
		} elsif( $priming_site =~ /^B/ ) {
		  $primerType = "walk";
		} elsif ($priming_site =~ /F/){
		  $primerType = "univ fwd";
		} elsif ($priming_site =~ /R/){
		  $primerType = "univ rev";
		} elsif ($priming_site eq "3" ){ #mg
		  $primerType = "univ rev";   #mg
		}else {
			print "Can't figure out primer information from tag: \"$chemTag\"! PrimerSite = $priming_site\n";
			die "\n";
		}


		# kludge for transposon cDNA
		if( $phdFileName =~ /^TL/ || $phdFileName =~ /^GT/ || $phdFileName =~ /^RT/)#RT added 27Jan03 by Rene Warren
		{
			$primerType = "walk";
			#$templateName = "TLxxx"; #$library;	# "TL05A"
			$templateType = "bac";
		}
	}

print "$phdFileName\t$templateName\t$templateType\t$primerType\n";

	return( (0, $templateName, $templateType, $primerType) );
}


sub findTemplateNameAndTemplateTypeStLouisFormat {
  die "usage: &findTemplateNameAndTemplateTypeStLouisFormat( <phd file name> );"
    if ($#_ != 0 );

  my $phdFileName = $_[0];

  
  my $readName;
  ( $readName = $phdFileName ) =~ s/\.phd\.[0-9]+$//;

  $firstDot = index( $readName, "." );
  if ( $firstDot == -1 ) {
    print STDERR "Can't figure out template and primer type for $phdFileName since no dot in filename\n";
    return( (1) );
  }


  if ( $firstDot == length( $readName ) ) {
    print STDERR "Can't figure out template and primer type for $phdFileName since no letter after the dot\n";
    return( (1) );
  }

  
  if ( $firstDot == 0 ) {
    print STDERR "Can't figure out template for $phdFileName since file starts with a dot\n";
    return( (1) );
  }

  # at St Louis, template names never have underscores in them, except
  # for .c and .a which are ignored by autofinish anyway

  # If there is a walking read, it looks like this:
  # das42f09_10.x1.phd.1
  # das42f09_10.y1.phd.1
  # das42f09_10.b1.phd.1
  # das42f09_10.i1.phd.1
  # das42f09_10.g1.phd.1
  # 
  # Special chemistry custom primer reads are:
  # das42f09_a_10.b1.phd.1
  # das42f09_a_10.g1.phd.1
  # das42f09_a_10.x1.phd.1
  # das42f09_a_10.y1.phd.1

  # Special chemistry universal primer reads are:
  # das42f09_a.b1.phd.1
  # das42f09_a.g1.phd.1
  # das42f09_a.x1.phd.1
  # das42f09_a.y1.phd.1


  my $nameUpToDot = substr( $readName, 0, $firstDot );

  my $firstUnderscore = index( $nameUpToDot, "_" );
  my $templateName;
  my $readType;
  my $universalPrimer = 1;
  my $customPrimer = 2;
  if ( $firstUnderscore == -1 ) {
    # there is no underscore in the name
    # like this: das42f09.x1.phd.1

    $templateName = $nameUpToDot;
    $readType = $universalPrimer;
  }
  else {
    # there is an underscore in the name

    $templateName = $nameUpToDot;
    # cut off from underscore to end
    # das42f09_10.x1.phd.1 has template name das42f09

    $templateName =~ s/_.*//;


    # if there is an underscore
    ( $afterUnderscore = $nameUpToDot ) =~ s/^[^_]*_//;


    if ( $afterUnderscore =~ /^[0-9]+$/ ) {
	# e.g., sl77h06_6.z2.phd.3
      # sl76c01_19.g1.phd.1
      # ot50d12_30.i1.phd.1

      $readType = $customPrimer;
    }
    elsif ( $afterUnderscore =~ /^[a-zA-Z]+$/ ) {
      # e.g., wv81h09_c.g1.phd.1 
      $readType = $universalPrimer;
    }
    elsif ( $afterUnderscore =~ /^[0-9]+[a-zA-Z]+$/ ) {
      # old style oligo walks with special chemistry, like this:
      # e.g., ot54e02_32g.i1.phd.1 oy47e12_31g.i1.phd.1 sl76c01_12g.b1.phd.1
      $readType = $customPrimer;
    }
    elsif ( $afterUnderscore =~ /^[a-zA-Z]+_[0-9]+$/ ) {
      # new style oligo walks with special chemistry
      $readType = $customPrimer;
    }
    else {
       print STDERR "1: Can't figure out read type for $phdFileName since after the underscore is $afterUnderscore\n";
       return( (1) );
    } # if ( $afterUnderscore ...
  } # if ( $firstUnderscore ...

  # if reached here, know whether the read is a universal primer or
  # custom primer
  # In the case of universal primer, we need to know whether it is 
  # universal forward or universal reverse

  $letterAfterFirstDot = substr( $readName, $firstDot + 1, 1 );

  my $templateType = "";
  if ( $letterAfterFirstDot eq 'a' ||
       $letterAfterFirstDot eq 'c' ) {
    $templateType = 'fake';
  }


  my $orientation;
  if ( $letterAfterFirstDot =~ /[ryg]/ ) {
    $orientation = "rev";
  }
  else {
    $orientation = "fwd";
  }


  my $primerType;
  if ( $readType == $universalPrimer ) {
    if ( $orientation eq "fwd" ) {
      $primerType = "univ fwd";
    }
    else {
      $primerType = "univ rev";
    }
  }
  elsif ( $readType == $customPrimer ) {
      if ( $templateType eq "pcr" ) {
         $primerType = "pcr end";
      }
      else {
         $primerType = "walk";
      }
  }
  else {
    die "read $phdFileName should either be a custom primer or universal primer";
  }

  return( (0, $templateName, $templateType, $primerType) );
}





sub findTemplateNameAndTemplateTypeSeattleFormat {

   die "usage: &findTemplateNameAndTemplateTypeSeattleFormat( <phd file name> );"
    if ( $#_ != 0 );

   my $phdFileName = $_[0];


   $firstDot = index( $phdFileName, "." );
   if ( $firstDot == -1 ) {
      print STDERR "Can't figure out template and primer type for $phdFileName since no dot in filename\n";
      return( (1) );
   }

   
   if ( $firstDot == length( $phdFileName )  ) {
      print STDERR "Can't figure out template and primer type for $phdFileName since no letter after the dot\n";
      return( (1) );
   }

   if ( $firstDot == 0 ) {
      print STDERR "Can't figure out template for $phdFileName since filename starts with a dot\n";
      return( (1) );
   }
           

   $templateName = substr( $phdFileName, 0, $firstDot );

   $firstUnderscore = index( $templateName, "_" );
   if ( $firstUnderscore != -1 ) {

      if ( $firstUnderscore == 0 ) {
         print STDERR "Filename must not start with an underscore in phd file $phdFileName\n";
         return( (1) );
      }

      if ( $firstUnderscore == ( length( $templateName ) - 1 ) ) {
         print STDERR "Strange phd file $phdFileName with underscore just before the dot\n";
         return( (1) );
      }

      $lastPartOfTemplateName = substr( $templateName, $firstUnderscore + 1, 
                           length( $templateName ) - $firstUnderscore - 1 );
   

      if ( $lastPartOfTemplateName =~ /^pcr/ ) {
         $templateType = "pcr";
      }
      elsif ($lastPartOfTemplateName eq "bac" ) {
         $templateType = "bac";
      }
      elsif ($lastPartOfTemplateName eq "cos" ) {
         $templateType = "cosmid";
      }
      else {
         # M13 or plasmid, but can't tell which oen
         $templateType = "";
      }
   }
   else {
      $templateType = "";
   }


   $afterDot = substr( $phdFileName, 
                       $firstDot + 1,
                       length( $phdFileName ) - $firstDot - 1 );


   # cut off the .phd.# part of the name

   $afterDot =~ s/\.phd.[0-9]+$//;


   $universalPrimer = 1;
   $customPrimer = 2;


   if ( $afterDot =~ /^[a-z][0-9]+u[0-9]+_/ ) {
      $readType = $universalPrimer;
      print "found $phdFileName\n";
   }
   elsif ( $afterDot =~ /^[a-z][0-9]+u[0-9]+$/ ) {
      # e.g., read name djs233_2255.x1u2
      $readType = $universalPrimer;
   }
   elsif( $afterDot =~ /^[a-z][0-9]+_up/ ) {
      $readType = $universalPrimer;
   }
   elsif( $afterDot =~ /^[a-z][0-9]+r[0-9]+p[0-9]+_/ ) {
      # e.g., read name 
      $readType = $customPrimer;
      print "found $phdFileName\n";
   }
   elsif( $afterDot =~ /^[a-z][0-9]+r[0-9]+p[0-9]+$/ ) {
      # e.g., read name djs233_2857.x1r2p5
      $readType = $customPrimer;
   }
   elsif( $afterDot =~ /^[a-z][0-9]+$/ ) {
      # e.g., read name djs233_1401.s1
      # is this true?  Cindy says "yes".
      $readType = $universalPrimer;
   }
   elsif( $afterDot =~ /^[a-z][0-9][\x7f]$/ ) {
      # e.g., read name djs77_3073.s1^?
      # Cindy says that these are all universal primer reads
      $readType = $universalPrimer;
   }
   elsif( $afterDot =~ /^[a-z][0-9]+_/ ) {
     if( $afterDot =~ /^[a-z][0-9]+_[0-9]+$/ ) {
       # e.g., read name djs233_1357.x1_7
       # is this true?  Cindy says yes.
       $readType = $customPrimer;
     }
     elsif( $afterDot =~ /^[a-z][0-9]+_[0-9]+[abc]$/ ) {
       # e.g., read name djs77_1377.x1_02b
       # is a custom primer, according to Cindy
       $readType = $customPrimer;
     }
     else {
       # e.g., read name djs10_776.x1_excellent
       $readType = $universalPrimer;
     }
   }
   else {
      print STDERR "unrecognized read type $phdFileName with afterDot $afterDot\n";
      print "afterDot length = ", length( $afterDot ), "\n";

      return( (1) );
   }
   

   $letterAfterDot = substr( $afterDot, 0, 1 );

   if ( $letterAfterDot =~ /[sfxzibtped]/ ) {
      $orientation = "fwd";
   }
   elsif( $letterAfterDot =~ /[ryg]/ ) {
      $orientation = "rev";
   }
   elsif( $letterAfterDot =~ /[ac]/ ) {
      # the cases that don't make sense
      $orientation = "fwd";
   }
   else {
      print STDERR "unknown orientation for read $phdFileName\n";
      return( (1) );
   }



   if ( $readType == $universalPrimer ) {
      if ( $orientation eq "fwd" ) {
         $primerType = "univ fwd";
      }
      else {
         $primerType = "univ rev";
      }
   }
   elsif ( $readType == $customPrimer ) {
      if ( $templateType eq "pcr" ) {
         $primerType = "pcr end";
      }
      else {
         $primerType = "walk";
      }
   }
   else {
      die "read $phdFileName should either be a custom primer or universal primer";
   }


  return( (0, $templateName, $templateType, $primerType) );
}
   


sub bIsPhdFileAlreadyProcessed {
   die "usage: &bIsPhdFileAlreadyProcessed( <phd filename> );" 
     if ( $#_ != 0 );
          
   my $szPhdFile = $_[0];
   my $szFullPath = $phdDirPath . "/" . $szPhdFile;
          
   $bSuccess = open( filPhd, $szFullPath );
   if (!$bSuccess ) {
      print STDERR "couldn't open $szFullPath\n";
      # don't try to process this file--something is wrong with it
      return 1;
   }
          
   while( <filPhd> ) {
      if (/^END_DNA/ ) {
         last;
      }
   }

   while( <filPhd> ) {
      if ( /^WR\{/ ) {
         if ( /^WR\{[\s]*$/ ) {
            # found a whole read item
            # let's see if it is a primer or template 
            defined( $_ = <filPhd> ) || die "premature end of phd file $szPhdFile while inside a WR{ item";

            @aWords = split;

            if ( ( $aWords[0] eq "primer" ) ||
                 ( $aWords[0] eq "template" ) ) {
               return( 1 );
            }
         }
      }
   }

   
   # if reached here, there must be no whole read item of type primer or template

   return( 0 );
}



          
          
   
   

   
      
sub szGetDateForTag {
  my $szDate;
  ($nSecond, $nMinute, $nHour, $nDayInMonth, $nMonth, $nYear, $wday, $yday, $isdst ) = localtime;

  undef $isdst;
  undef $yday;
  undef $wday;
  
  if ( $nYear > 100 ) {
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

sub processOnePhdFile {

   die "usage: &processOnePhdFile( <phd filename> )"
     if ( $#_ != 0 );


   my $phdFileName = $_[0];

   my $templateName;
   my $templateType;
   my $primerType;


   # Altered 31 Oct 2000 S. Zuyderduyn
  # ($bErrorFlag, $templateName, $templateType, $primerType ) = &findTemplateNameAndTemplateTypeStLouisFormat( $phdFileName );
  ($bErrorFlag, $templateName, $templateType, $primerType ) = &findTemplateNameAndTemplateTypeVancouverFormat( $phdFileName );


   return if ( $bErrorFlag );

   if ( $bJustView ) {
     print "read $phdFileName is $primerType template: $templateName $templateType\n";
   }
   else {

     $szDate = &szGetDateForTag();

     my $szPhdFileFullPath = $phdDirPath . "/" . $phdFileName;

     open( filPhd, ">>$szPhdFileFullPath" ) || die "couldn't open $szPhdFileFullPath for append";

     print( filPhd "\n" );
     print( filPhd "WR{\n" );
     print( filPhd "template dscript $szDate\n" );
     print( filPhd "name: $templateName\n" );
     if ( $templateType ne "" ) {
       print( filPhd "type: $templateType\n" );
     }
     print( filPhd "}\n" );

     print( filPhd "\n" );
     print( filPhd "WR{\n" );
     print( filPhd "primer dscript $szDate\n" );
     print( filPhd "type: $primerType\n" );
     print( filPhd "}\n" );

     close( filPhd );
   }
}




