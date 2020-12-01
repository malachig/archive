#!/usr/bin/perl

##NAME: findSNPRichGenes.pl
##DESCRIPTION: Finds genes between a specified region and the number of SNPS that they contain
##Reports back the single, most SNP rich gene.
##WRITTEN FOR GENOMICS 2007 CBW WORKSHOP
##STUDENT NAME:________________________

##IMPORTS
use lib "/home/obig/lib/bioperl-live";
use lib "/home/obig/lib/ensembl_43_perl_API/ensembl/modules";
use lib "/home/obig/lib/ensembl_43_perl_API/ensembl-variation/modules";
use Bio::EnsEMBL::Registry;
use Getopt::Long;

##GET INPUT PARAMETERS,
##CHROMOSOME, START (ABSOLUTE), END (ABSOLUTE), AND STRAND
my $chromosome;
my $start;
my $end;
my $strand;
my $species;

##LOAD EACH COMMAND-LINE PARAMETER INTO THE APPROPRIATE VARIABLE
GetOptions ('chromosome=s'=>\$chromosome,
	    'start=s'=>\$start,
	    'end=s'=>\$end,
	    'species=s'=>\$species,
	    'strand=s'=>\$strand);

##VERIFY THAT EVERY PARAMETER IS DEFINED
if (!defined($chromosome) || !defined($start) || !defined($species)
        || !defined($end) || !defined($strand)) {

        print "\nSome parameters not defined, ";
        print "usage ./findSNPRichGenes.pl ";
        print "--chromosome 1 --start 1 --end 100000 --strand 1 --species Human\n";
        exit(1);        ##EXIT ERROR CODE, NON-ZERO
} elsif (($end-$start) > 1000000) {
        print "\nRegion too big.  Use less than a MB.  Don't hog the bandwidth!\n";
        exit(1);
}

##CONNECT TO ENSEMBL DATABASE
my $registry = "Bio::EnsEMBL::Registry";
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
    -db_version => '43'
);

##GET A CORE DATABASE ADAPTOR
my $db = Bio::EnsEMBL::Registry->get_DBAdaptor($species, "core");

##GET A VARIATION DATABASE ADAPTOR
my $vdb = Bio::EnsEMBL::Registry->get_DBAdaptor($species, "variation");

##GET A VARIATION FEATURE ADAPTOR
my $variationfeature_adaptor = $vdb->get_VariationFeatureAdaptor();

##GET THE SLICE
my $slice_adaptor = $db->get_SliceAdaptor();
my $slice = $slice_adaptor->fetch_by_region('chromosome', $chromosome, $start, $end, $strand);

#Declare variables to be used below
my $snp_rich_gene_name = undef;
my $non_syn_snp_rich_gene_name = undef;
my $snp_frequency_gene_name = undef;
my $highest_snp_count = 0;
my $highest_non_syn_snp_count = 0;
my $highest_snp_frequency;

##GET ALL GENES IN THE SLICE
my $genes = $slice->get_all_Genes();

##FOR EACH GENE DETERMINE THE NUMBER OF SNPS
foreach my $gene (@$genes) {
        my $snp_count=0;
	my $non_syn_snp_count=0;
	my $snp_frequency=0;
        my $subslice = $slice->sub_Slice($gene->start, $gene->end, $gene->strand);
        if (!defined $subslice) {
                $subslice = $slice_adaptor->fetch_by_region('chromosome', $chromosome, ($slice->start + $gene->start - 1), ($slice->start + $gene->end - 1));
        }

	#Get all SNPs for the gene
        my @snps = @{$variationfeature_adaptor->fetch_all_by_Slice($subslice)};
        foreach my $snp (@snps) {
	  my $start=$snp->{'start'};
	  my $end=$snp->{'end'};
	  my $consequence_type=$snp->{'consequence_type'}[0]; #Consequence type is the type of SNP (UTR, CODING, SYNONYMOUS, NON-SYNONYMOUS, ETC)
	  if ($start == $end) { $snp_count++; }
	  if ($consequence_type eq "NON_SYNONYMOUS_CODING") {$non_syn_snp_count++;}
        }

	#Determing absolute start and end of the gene so that gene length and SNP frequency can be determined
        my $abs_gene_start = $start - 1 + $gene->start;
        my $abs_gene_end = $start - 1 + $gene->end;
	my $gene_length = $abs_gene_end - $abs_gene_start;
	my $gene_strand = $gene->strand;
	my $snp_frequency = $snp_count/$gene_length;
        print "GENE NAME: " . $gene->display_id;
        print " ($abs_gene_start, $abs_gene_end, $strand, $gene_length) ";
        print "SNP COUNT: $snp_count NON-SYNONYMOUS SNP COUNT: $non_syn_snp_count SNP FREQUENCY: $snp_frequency\n";

	#Keep track of Gene with highest SNP count
        if ($snp_count > $highest_snp_count) {
                $snp_rich_gene_name = $gene->display_id;
                $highest_snp_count = $snp_count;
        }

	#Keep track of Gene with highest non-synonymous SNP count
        if ($non_syn_snp_count > $highest_non_syn_snp_count) {
                $non_syn_snp_rich_gene_name = $gene->display_id;
                $highest_non_syn_snp_count = $non_syn_snp_count;
        }

	#Keep track of Gene with highest SNP frequency (SNPs/bp)
        if ($snp_frequency > $highest_snp_frequency) {
                $snp_frequency_gene_name = $gene->display_id;
                $highest_snp_frequency = $snp_frequency;
        }

}

#Print out results
print "SNP RICH GENE IS: " . $snp_rich_gene_name . " WITH: " . $highest_snp_count . " SNPS\n";
print "NON-SYN SNP RICH GENE IS: " . $non_syn_snp_rich_gene_name . " WITH: " . $highest_non_syn_snp_count . " SNPS\n";
print "HIGHEST FREQUENCY SNP RICH GENE IS: " . $snp_frequency_gene_name . " WITH: " . $highest_snp_frequency . " SNPS/BP\n";

