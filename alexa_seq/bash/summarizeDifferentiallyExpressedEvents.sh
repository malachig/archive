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

#Define list of files - Genes, Transcripts, Exon Regions, Junctions, Boundaries, Introns, Active Intron Regions, Silent Intron Regions, Intergenics, Active Intergenic Regions, Silent Intergenic Regions
FILES="
/projects/malachig/solexa/figures_and_stats/DE/ENST_v53/Mip101_vs_Mip5FuR_Gene_DE_Values_Significant_MTC.txt
/projects/malachig/solexa/figures_and_stats/DE/Transcripts_v53/Mip101_vs_Mip5FuR_Transcript_DE_Values_Significant_MTC.txt
/projects/malachig/solexa/figures_and_stats/DE/ENST_v53/Mip101_vs_Mip5FuR_ExonRegion_DE_Values_Significant_MTC.txt
/projects/malachig/solexa/figures_and_stats/DE/Junctions_v53/Mip101_vs_Mip5FuR_Junction_DE_Values_Significant_MTC.txt
/projects/malachig/solexa/figures_and_stats/DE/Junctions_v53/Mip101_vs_Mip5FuR_KnownJunction_DE_Values_Significant_MTC.txt
/projects/malachig/solexa/figures_and_stats/DE/Junctions_v53/Mip101_vs_Mip5FuR_NovelJunction_DE_Values_Significant_MTC.txt
/projects/malachig/solexa/figures_and_stats/DE/Boundaries_v53/Mip101_vs_Mip5FuR_Boundary_DE_Values_Significant_MTC.txt
/projects/malachig/solexa/figures_and_stats/DE/Boundaries_v53/Mip101_vs_Mip5FuR_KnownBoundary_DE_Values_Significant_MTC.txt
/projects/malachig/solexa/figures_and_stats/DE/Boundaries_v53/Mip101_vs_Mip5FuR_NovelBoundary_DE_Values_Significant_MTC.txt
/projects/malachig/solexa/figures_and_stats/DE/Introns_v53/Mip101_vs_Mip5FuR_Intron_DE_Values_Significant_MTC.txt
/projects/malachig/solexa/figures_and_stats/DE/Introns_v53/Mip101_vs_Mip5FuR_SilentIntronRegion_DE_Values_Significant_MTC.txt
/projects/malachig/solexa/figures_and_stats/DE/Introns_v53/Mip101_vs_Mip5FuR_ActiveIntronRegion_DE_Values_Significant_MTC.txt
/projects/malachig/solexa/figures_and_stats/DE/Intergenics_v53/Mip101_vs_Mip5FuR_Intergenic_DE_Values_Significant_MTC.txt
/projects/malachig/solexa/figures_and_stats/DE/Intergenics_v53/Mip101_vs_Mip5FuR_SilentIntergenicRegion_DE_Values_Significant_MTC.txt
/projects/malachig/solexa/figures_and_stats/DE/Intergenics_v53/Mip101_vs_Mip5FuR_ActiveIntergenicRegion_DE_Values_Significant_MTC.txt
"

for f in $FILES
do
  getStats $f
done




