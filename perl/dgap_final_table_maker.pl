#!/perl/bin/perl

use strict;
use Getopt::Long;

my %opts = ();
GetOptions (\%opts,'f=s',
		   'c=s',
		   'd=s',
		   'h=s',
		   't=s',
		   'm=s',
		   'o=s'); 

sub printHelp{
	print "\n";
	print "\tusage:	 dgap_final_table_maker.pl [options]\n";
	print "\t-f  	 TAD analysis file\n";
	print "\t-c  	 TAD position file\n";
	print "\t-d  	 DHS file\n";
	print "\t-h  	 Haplosensitivity file\n";
	print "\t-t  	 Triplosensitivity file\n";
	print "\t-m  	 Phenomatch score file\n";
	print "\t-o 	 Output file\n";
	exit;
}	

####### Variable declaration
my $file = $opts{f};
my $TADfile = $opts{c};
my $DHSfile = $opts{d};
my $haplofile = $opts{h};
my $triplofile = $opts{t};
my $phenoscorefile = $opts{m};
my $outfile = $opts{o};
my ($start,$end,$DosMb);
my @temp_data;
my @id;
my @TADS;
my @name;
my $DosMb=0;
my $TAD=0;
my $DHS=0;
my $haplo=0;
my $triplo=0;
my $pheno=0;
my $maxpheno=0;
my $pheno_perc=0;
my $maxpheno_perc=0;

#File structure	
#0	Breakpoint_ID	
#1	chr	
#2	start	
#3	end	
#4	win_start	
#5	win_end	
#6	Hi_domain	
#7	HiC_chr	
#8	HiC_start	
#9	HiC_end	
#10	HI_gene_chr	
#11	HI_gene_start	
#12	HI_gene_end	
#13	HI_gene_name	
#14	LOD	
#15	HI_prob	
#16	HI_ensembl_gene_ID	
#17	HI_chr	
#18	Gene_start	
#19	gene_end	
#20	HGNC_ID	
#21	HGNC_symbol	
#22	Gene_type	
#23	WikiGene_Description	
#24	Description	
#25	Phenotype_description	
#26	DGAP_ID	
#27	DGAP_karyotype	
#28	DGAP_phenotype

#perl dgap_final_table_maker.pl -f final_list_6mb_run.txt -c hESC_hg37_domains.bed -d DHS_promoter_broken_breakpDGAP.txt -h ClinGen_haploinsufficiency_gene.bed -t ClinGen_triplosensitivity_gene.bed -m percentiles_phenomatch.txt -o TEST_FINAL.txt

###Open file and read positions from chr, start and end

open (FILE,">$outfile") || die "cannot create $outfile\n";

open (VAR,"$file");
while (<VAR>){

	$DosMb=0;
	$TAD=0;
	$DHS=0;
	$haplo=0;
	$triplo=0;
	$pheno=0;
	$maxpheno=0;
	$pheno_perc=0;
	$maxpheno_perc=0;
	$start = 0;
	$end = 0;

	chomp($_);
	print FILE "$_";
	@temp_data = split (/\t/,$_);
	@id = split(/\_/,$temp_data[0]); #$id[0] = DGAP#

	$start = $temp_data[2]-1000000;
	$end = $temp_data[3]+1000000;

	#check the 2Mb window inclusion
	if (($temp_data[11] <= $start && $temp_data[12] >= $start ) || ($temp_data[11] <= $end && $temp_data[12] >= $end ) || ($temp_data[11] >= $start && $temp_data[12] <= $end )){$DosMb=1;}

	#########################################
	#######check gene inclusion in TAD domain
	#########################################
		
	if (($temp_data[8] <= $temp_data[11] && $temp_data[9] >= $temp_data[11]) || ($temp_data[8] <= $temp_data[12] && $temp_data[9] >= $temp_data[12]) || ($temp_data[8] >= $temp_data[11] && $temp_data[9] <= $temp_data[12])){$TAD=1;}


	#####################################################
	#######check DHS enhancer/promoter disrupted contacts
	#####################################################
	#DHS file structure
	#0	1	2	3	4	5	6	7	8	9	10	11	12	13
	#DGAP002_2_B	chr16	75871424	75871425	72871424	78871425	chr16	75529320	75529470	CHST6	chr16	75977780	75977930	0.906709

	open (TAD,"$DHSfile");
	while(<TAD>){
		chomp($_);
		@TADS = split (/\t/,$_);
		@name = split(/\_/,$TADS[0]); #$id[0] = DGAP#
		if (($name[0] eq $id[0]) && ($TADS[9] eq $temp_data[13])){$DHS=1;last;}
	}
	close(TAD);

	################################
	#####check gene haplosensitivity
	################################
	#ClinGen haplosensitive file structure
	#chr1	103342022	103574052	COL11A1	1
	#chr1	110091185	110138452	GNAI3	0
	open (TAD,"$haplofile");
	while(<TAD>){
		chomp($_);
		@TADS = split (/\t/,$_);
		if (($temp_data[1] eq $TADS[0]) && ($TADS[3] eq $temp_data[13])){$haplo=$TADS[4];last;}
	}
	close(TAD);


	################################
	####check gene triplosensitivity
	################################
	#ClinGen triplosensitive file structure
	#chr1	103342022	103574052	COL11A1	0
	#chr1	110091185	110138452	GNAI3	0
	open (TAD,"$triplofile");
	while(<TAD>){
		chomp($_);
		@TADS = split (/\t/,$_);
		if (($temp_data[1] eq $TADS[0]) && ($TADS[3] eq $temp_data[13])){$triplo=$TADS[4];last;}
	}
	close(TAD);

	###################################################
	#########Obtain phenomatch and maxPhenomatch scores
	###################################################
	#NOTE: Please note that the included phenomatch scores file correspond to the non-coding DGAP dataset ONLY. If you want to analyze your own set of subject phenotypes, Jonas Ibn-Salem is currently writing a paper and creating a software that performs this analysis. If you want to urgently analyze your data or have a collaboration with him, please contact Jonas for the phenomatch analysis of your dataset. 

	#Pheno percentile file structure
	#name	gene_symbol	chr	phenoMatchScore	maxPhenoMatchScore	#genes_6Mb_window_DGAP_case	Percentile_Path_Gene	Min_phenoscore	Max_phenoscore	Percentile_MaxScorePath_Gene	Min_Maxscorephenoscore	Max_Maxphenoscore

	open (TAD,"$phenoscorefile");
	while(<TAD>){
		chomp($_);
		@TADS = split (/\t/,$_);
		if (($TADS[0] eq $id[0]) && ($TADS[1] eq $temp_data[13])){
			$pheno=$TADS[3];
			$maxpheno=$TADS[4];
			$pheno_perc=$TADS[6];
			$maxpheno_perc=$TADS[9];
			last;
		}
	}
	close(TAD);
	print FILE "\t1\t$DosMb\t$TAD\t$DHS\t$haplo\t$triplo\t$pheno\t$maxpheno\t$pheno_perc\t$maxpheno_perc\n";

}
close(VAR);
close(FILE);

