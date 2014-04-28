#!/bin/bash
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

E_BADARGS=65
if [ ! -n "$1" ]
then
  echo "Usage: `basename $0` summarizeReadAssignments.sh /projects/malachig/solexa/read_records/HS0499"
  exit $E_BADARGS
fi  

zcat $1/*.txt.gz | cut -f 4,5 | perl -ne '$r_count+=2; if ($_ =~/(\S+)\s+(\S+)/){$counts{$1}++; $counts{$2}++;} if (eof){print "\nReads: $r_count\n"; foreach my $c (sort {$a cmp $b} keys %counts){print "$c: $counts{$c}\n"};}';  







