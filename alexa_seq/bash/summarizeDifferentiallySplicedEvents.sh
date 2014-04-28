#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

function getStats() {
  echo
  COUNT=$(wc -l $1)
  echo Line count: $COUNT
}

#Define list of files - Transcripts, Exon Regions, Junctions, Boundaries, Introns, Active Intron Regions, Silent Intron Regions
FILES="
/projects/malachig/solexa/figures_and_stats/SI/Transcripts_v53/Mip101_vs_Mip5FuR_Transcript_SI_Values_Sorted_Cutoff.txt
/projects/malachig/solexa/figures_and_stats/SI/ENST_v53/Mip101_vs_Mip5FuR_ExonRegion_SI_Values_Sorted_Cutoff.txt
/projects/malachig/solexa/figures_and_stats/SI/Junctions_v53/Mip101_vs_Mip5FuR_Junction_SI_Values_Sorted_Cutoff.txt
/projects/malachig/solexa/figures_and_stats/SI/Junctions_v53/Mip101_vs_Mip5FuR_KnownJunction_SI_Values_Sorted_Cutoff.txt
/projects/malachig/solexa/figures_and_stats/SI/Junctions_v53/Mip101_vs_Mip5FuR_NovelJunction_SI_Values_Sorted_Cutoff.txt
/projects/malachig/solexa/figures_and_stats/SI/Boundaries_v53/Mip101_vs_Mip5FuR_Boundary_SI_Values_Sorted_Cutoff.txt
/projects/malachig/solexa/figures_and_stats/SI/Boundaries_v53/Mip101_vs_Mip5FuR_KnownBoundary_SI_Values_Sorted_Cutoff.txt
/projects/malachig/solexa/figures_and_stats/SI/Boundaries_v53/Mip101_vs_Mip5FuR_NovelBoundary_SI_Values_Sorted_Cutoff.txt
/projects/malachig/solexa/figures_and_stats/SI/Introns_v53/Mip101_vs_Mip5FuR_Intron_SI_Values_Sorted_Cutoff.txt
/projects/malachig/solexa/figures_and_stats/SI/Introns_v53/Mip101_vs_Mip5FuR_SilentIntronRegion_SI_Values_Sorted_Cutoff.txt
/projects/malachig/solexa/figures_and_stats/SI/Introns_v53/Mip101_vs_Mip5FuR_ActiveIntronRegion_SI_Values_Sorted_Cutoff.txt
"

for f in $FILES
do
  getStats $f
done




