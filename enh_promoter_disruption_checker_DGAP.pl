#!/usr/bin/perl

use strict;
use Getopt::Long;

my %opts = ();
GetOptions (\%opts,'f=s',
		   'd=s',
		   'a=i',
		   'o=s',
		   'help|h'); 

sub printHelp{
	print "\n";
	print "\tusage:	 enh_promoter_disruption_checker_DGAP.pl [options]\n";
	print "\t-f  	 analysis regions\n";
	print "\t-d  	 input file with DHS-promoter data\n";
	print "\t-a  	 amount of bp to surround each region\n";
	print "\t-o 	 Output file\n";
	print "\tExample\n";
	print "\tperl enh_promoter_disruption_checker_DGAP.pl -f non_coding_DGAP_positions.bed -d genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8 -a 3000000 -o DHS_promoter_broken_DGAP.txt\n";
	exit;
}	

####### Variable declaration
my $file_original = $opts{f};
my $fileDHS = $opts{d};
my $addition = $opts{a};
my $outfile = $opts{o};
my ($chr,$llave);
my @temp_data;
my %breakps;
my $start;
my $end;

#Get DGAP breakpoint positions
#File structure
#0	1		2		3
#chr5	177479707	177479712	DG100_A
#chrX	44900487	44900489	DG100_B

###Open file and read positions from chr, start and end
open (VAR,"$file_original") || die "Cannot open $file_original";
while (<VAR>){
	chomp($_);
	@temp_data = split (/\s+/,$_);
	$breakps{$temp_data[3]}[0] = $temp_data[0];		#chr
	$breakps{$temp_data[3]}[1] = $temp_data[1]-$addition;	#window start	
	$breakps{$temp_data[3]}[2] = $temp_data[2]+$addition;	#window end
	$breakps{$temp_data[3]}[3] = $temp_data[1];		#start	
	$breakps{$temp_data[3]}[4] = $temp_data[2];		#end
	$breakps{$temp_data[3]}[5] = 0;				#enh/promoter contacts broken count
}
close(VAR);

#DHS enhancer/promoter contacts file structure
#0	1		2		3		4	5		6		7
#chr22	24979420	24979570	GGT1		chr22	25233680	25233830	0.766967
#chr22	24979420	24979570	GGT1		chr22	25335260	25335410	0.794857
#chr22	24979420	24979570	GGT1		chr22	25448040	25448190	0.804424
#chr22	24989040	24989190	C22orf36	chr22	24581040	24581190	0.771504
#chr22	24989040	24989190	C22orf36	chr22	24581660	24581810	0.838512
#chr22	24989040	24989190	C22orf36	chr22	24582040	24582190	0.736335

foreach $llave (sort {$a<=>$b} keys %breakps){

	open (VAR,"$fileDHS") || next;

	while (<VAR>){
	
		chomp($_);
		@temp_data = split (/\s+/,$_);

		if($temp_data[2] <= $temp_data[5]){
			$start=$temp_data[1];
			$end=$temp_data[6];
		}

		else{
			$start=$temp_data[5];
			$end=$temp_data[2];
		}

		#We filter for regions with 0.7 or higher correlation. If you want to change stringency, modify this cut-off value
		if(($temp_data[0] eq $temp_data[4]) && ($breakps{$llave}[0] eq $temp_data[0]) && ($start >= $breakps{$llave}[1]) && ($end <= $breakps{$llave}[2]) && ($temp_data[7] > 0.7)){
		
			print "entro aqui\n";
			if($breakps{$llave}[3] >= $start && $breakps{$llave}[4] <= $end){

				open (FILE,">>$outfile") || die "cannot create $outfile\n";
				print FILE "$llave\t$breakps{$llave}[0]\t$breakps{$llave}[3]\t$breakps{$llave}[4]\t$breakps{$llave}[1]\t$breakps{$llave}[2]\t$_\n";
				close(FILE);
				$breakps{$llave}[5]=$breakps{$llave}[5]+1;

			}

			else{next;}
		}

		else{next;}
	}

	close(VAR);
}



#print summary file
open (FILE,">${outfile}_summary.out") || die "cannot open ${outfile}_summary.out\n";

foreach $llave (sort {$a<=>$b} keys %breakps) {
	print FILE "$llave\t$breakps{$llave}[0]\t$breakps{$llave}[3]\t$breakps{$llave}[4]\t$breakps{$llave}[1]\t$breakps{$llave}[2]\t$breakps{$llave}[5]\n";
}

close(FILE);
