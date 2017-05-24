#!/usr/bin/perl

use strict;
use Getopt::Long;

my %opts = ();
GetOptions (\%opts,'f=s',
		   'i=s',
		   'g=s',
		   'c=s',
		   'n=i',
		   'o=s',
		   'help|h'); 

sub printHelp{
	print "\n";
	print "\tusage:	 dgap_features_check.pl [options]\n";
	print "\t-f  	 bed coords input file\n";
	print "\t-i  	 haploinsufficiency file\n";
	print "\t-g  	 ensembl file\n";
	print "\t-c  	 HiC file\n";
	print "\t-n  	 window size\n";
	print "\t-o 	 Output folder\n";
	print "\tExample\n";
	print "\tperl perl dgap_features_check.pl -f DGAP_control_positions_final.bed -i HI_Predictions_Version3.bed -g GRCh37.p13_ensembl_genes.txt -c hESC_hg37_domains.bed -n 3000000 -o /home/czm/Documents/TAD_analysis/DGAP_positivecontrols\n";
	exit;
}	

####### Variable declaration
my $file = $opts{f};
my $haplo = $opts{i};
my $ensemble = $opts{g};
my $HiC = $opts{c};
my $win = $opts{n};
my $outfolder = $opts{o};
my ($chr,$start,$end,$llave,$newchr,$name);
my @temp_data;
my @temp_data2;
my @temp_data3;
my @temp_data4;
my %breakps;
my %genes;

#######################################################################
###Open file and read positions from chr, start and end
#File structure
#0	1		2		3
#chr5	177479707	177479712	DG100_A
#chrX	44900487	44900489	DG100_B

open (VAR,"$file");
while (<VAR>){
	chomp($_);
	@temp_data = split (/\s+/,$_);
		$breakps{$temp_data[3]}[0] = $temp_data[0];		#chr
		$breakps{$temp_data[3]}[1] = $temp_data[1]-$win;	#start	
		$breakps{$temp_data[3]}[2] = $temp_data[2]+$win;	#end
		$breakps{$temp_data[3]}[3] = 0;				#HiC domain inclusion
		$breakps{$temp_data[3]}[4] = 0;				#HiC chr
		$breakps{$temp_data[3]}[5] = 0;				#HiC start
		$breakps{$temp_data[3]}[6] = 0;				#HiC end
}
close(VAR);

#######################################################################
###Open ensembl genes file
#0			1	2		3	4	5		6		7			8		9
#EnsemblGeneID		Chr	GeneStart	GeneEnd	HGNCID	HGNC_symbol	Gene_type	WikiGene_Description	Description	Phenotype_description
#ENSG00000261657	HG991_PATCH	66119285	66465398	20661	SLC25A26	protein_coding	solute carrier family 25 (S-adenosylmethionine carrier), member 26	solute carrier family 25 (S-adenosylmethionine carrier), member 26 [Source:HGNC Symbol;Acc:20661]	
#ENSG00000223116	13	23551994	23552136			miRNA			
#ENSG00000233440	13	23708313	23708703	19121	HMGA1P6	pseudogene		high mobility group AT-hook 1 pseudogene 6 [Source:HGNC Symbol;Acc:19121]	

open (VAR,"$ensemble");
while (<VAR>){
	chomp($_);
	@temp_data = split (/\s+/,$_);
	$genes{$temp_data[5]} = $_;
}
close(VAR);

#######################################################################
#Check HiC inclusion
#0	1	2
#chr1	770137	1250137
#chr1	1250137	1850140
#chr1	1850140	2330140

open (VAR,"$HiC");
open (FILE,">${outfolder}/HiC_list_DGAP.txt") || die "cannot open ${outfolder}/HiC_list_DGAP.txt\n";

while (<VAR>){
	chomp($_);
	@temp_data = split (/\s+/,$_);	

	foreach $llave (sort {$a<=>$b} keys %breakps) {
		if($temp_data[0] eq $breakps{$llave}[0]){

			$start = $breakps{$llave}[1]+$win;
			$end = $breakps{$llave}[2]-$win;

			if(($temp_data[1] <= $start && $temp_data[2] >= $start) || ($temp_data[1] <= $end && $temp_data[2] >= $end) || ($temp_data[1] >= $start && $temp_data[2] <= $end)){
				$breakps{$llave}[3]=1; 		   #inside HiC domain
				$breakps{$llave}[4]=$temp_data[0]; #hic chr
				$breakps{$llave}[5]=$temp_data[1]; #hic domain start
				$breakps{$llave}[6]=$temp_data[2]; #hic domain end
				print FILE "$llave\t$breakps{$llave}[0]\t$start\t$end\t$breakps{$llave}[3]\t$breakps{$llave}[4]\t$breakps{$llave}[5]\t$breakps{$llave}[6]\n";
			}
		}
		else {next;}
	}
}
close(VAR);
close(FILE);

#######################################################################
#Check HI genes overlap
#0	#1		#2		#3				#4
#chrX	37850070	37850569	HYPM|0.000043663|100%		0.000043663	.	37850070	37850569	0,255,0
#chr5	43039335	43043272	ANXA2R|0.000047349|100%		0.000047349	.	43039335	43043272	0,255,0
#chr11	61957688	61961011	SCGB1D1|0.000054551|99.99%	0.000054551	.	61957688	61961011	0,255,0

open (VAR,"$haplo");
open (FILE,">${outfolder}/HI_list_DGAP.txt") || die "cannot open ${outfolder}/HI_list_DGAP.txt\n";

while (<VAR>){

	chomp($_);
	@temp_data = split (/\s+/,$_);
	foreach $llave (sort {$a<=>$b} keys %breakps){

		if($temp_data[0] eq $breakps{$llave}[0]){

			if(($temp_data[1] <= $breakps{$llave}[1] && $temp_data[2] >= $breakps{$llave}[1]) || ($temp_data[1] <= $breakps{$llave}[2] && $temp_data[2] >= $breakps{$llave}[2]) || ($temp_data[1] >= $breakps{$llave}[1] && $temp_data[2] <= $breakps{$llave}[2])){

				@temp_data2 = split (/\|/,$temp_data[3]);
				@temp_data3 = split (/\%/,$temp_data2[2]); 

				#if($temp_data3[0]<=10){ #remove this line's comment and its bracket (below) if you only want to assess genes with stronger HI evidence. See Huang et al., 2010 paper for an explanation on this cut-off value.

				@temp_data4 = split (/\_/,$llave);
				$start = $breakps{$llave}[1]+$win;
				$end = $breakps{$llave}[2]-$win;
				print FILE "$llave\t$breakps{$llave}[0]\t$start\t$end\t$breakps{$llave}[1]\t$breakps{$llave}[2]\t$breakps{$llave}[3]\t$breakps{$llave}[4]\t$breakps{$llave}[5]\t$breakps{$llave}[6]\t$temp_data[0]\t$temp_data[1]\t$temp_data[2]\t$temp_data2[0]\t$temp_data2[1]\t$temp_data2[2]\t$genes{$temp_data2[0]}\n";

				#} #remove this line's comment if you only want to assess genes with stronger HI evidence. See Huang et al., 2010 paper for an explanation on this cut-off value.
			}
		}

		else {next;}
	}
}
close(VAR);
close(FILE);

