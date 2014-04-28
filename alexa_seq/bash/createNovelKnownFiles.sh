#!/bin/bash
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

# This script counts the jobs on the cluster for a user every N seconds
E_BADARGS=65
if [ ! -n "$1" ]
then
  echo "Usage: `basename $0` /projects/malachig/solexa/read_records/ HS04391 53  (i.e. [path] [library] [ensembl version])"
  exit $E_BADARGS
fi  

echo
echo "Using args: $1 $2 $3"

JPATH="$1/$2/Junctions_v$3/Summary/"
JFILE="$2_JunctionExpression_v$3.txt"
KJFILE="$2_KnownJunctionExpression_v$3.txt"
NJFILE="$2_NovelJunctionExpression_v$3.txt"

BPATH="$1/$2/Boundaries_v$3/Summary/"
BFILE="$2_BoundaryExpression_v$3.txt"
KBFILE="$2_KnownBoundaryExpression_v$3.txt"
NBFILE="$2_NovelBoundaryExpression_v$3.txt"

#echo $JPATH
#echo $JFILE
#echo $KJFILE
#echo $NJFILE

#echo $BPATH
#echo $BFILE
#echo $KBFILE
#echo $NBFILE

#Change to working dir for junctions 
cd $JPATH

#Get known for junctions
echo "Creating known junction expression file"
COL=$(head -n 1 $JFILE | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Supporting_EnsEMBL_Count/){print "\n\n$c\n\n"}}')
export COL
cat $JFILE | perl -ne 'chomp($_); $col=$ENV{"COL"}; @line=split("\t", $_); if ($_ =~ /Expressed/){print "$_\n";} if ($line[$col-1]>=1){print "$_\n";}' > $KJFILE

#Get novel for junctions 
echo "Creating novel junction expression file"
cat $JFILE | perl -ne 'chomp($_); $col=$ENV{"COL"}; @line=split("\t", $_); if ($_ =~ /Expressed/){print "$_\n";} if ($line[$col-1] eq "0"){print "$_\n";}' > $NJFILE

#Change to working dir for boundaries
cd $BPATH

#Get known for boundaries
echo "Creating known boundary expression file"
COL=$(head -n 1 $BFILE | perl -ne '@line = split("\t", $_); foreach $data (@line){$c++; if ($data =~/Supporting_EnsEMBL_Count/){print "\n\n$c\n\n"}}')
export COL
cat $BFILE | perl -ne 'chomp($_); $col=$ENV{"COL"}; @line=split("\t", $_); if ($_ =~ /Expressed/){print "$_\n";} if ($line[$col-1]>=1){print "$_\n";}' > $KBFILE

#Get novel for boundaries
echo "Creating novel boundary expression file"
cat $BFILE | perl -ne 'chomp($_); $col=$ENV{"COL"}; @line=split("\t", $_); if ($_ =~ /Expressed/){print "$_\n";} if ($line[$col-1] eq "0"){print "$_\n";}' > $NBFILE

echo
echo
