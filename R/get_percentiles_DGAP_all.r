########################################################################################################################
### Calculate percentiles of phenomatch scores for the analyzed gene values per DGAP case 
### Input: phenomatch and max_phenomatch scores file. Remember to delete the _A and _B from the file before processing.
### Output: file with percentiles for the analyzed pheno/maxphenomatch scores per gene.
### how to run: Rscript --vanilla get_percentiles_DGAP_all.r v07_breakpoint_window_with_HPO.bed.6MB_win.bed.annotated.out.overlapped_genes_formatted.txt > percentiles_6Mb_pheno_maxpheno.txt
########################################################################################################################
require(lattice)
require(stringr) # for str_split_fixed()

#arguments: 
args <- commandArgs(trailingOnly = TRUE)

allgenes <- args[1]
rm(args)

#input file format
#1	2	3	4	5		6		7		8
#chr	start	end	name	phenotypes	gene_symbol	phenoMatchScore MaxPhenoMatchScore
#chr18	60566041	66566041	DGAP322	HP:0008736;HP:0011342;HP:0000119;HP:0012856;HP:0001510;HP:0000047	SERPINB7	0.014	0.003
#chr18	60566041	66566041	DGAP322	HP:0008736;HP:0011342;HP:0000119;HP:0012856;HP:0001510;HP:0000047	BCL2	0.821	0.797
#chr10	7161499	13161499	DGAP287	HP:0001264;HP:0001251;HP:0002307;HP:0011343;HP:0001344	GATA3	1.822	1.365
#chr10	7161499	13161499	DGAP287	HP:0001264;HP:0001251;HP:0002307;HP:0011343;HP:0001344	DHTKD1	9.057	2.09

#remember to delete the _A and _B from the file before processing

allgene <- read.table(file = allgenes, header = TRUE, comment.char = "")

# delete the _A and _B from name column
allgene$name <- str_split_fixed(allgene$name, '_', 2)[, 1]

for (i in 1:nrow(allgene)){

	data <- allgene[allgene$name == allgene[i,4],]

	#write.table(data, file="out_matrices_phenomatch_6Mb.txt", append = TRUE, quote = FALSE, sep = "\t", eol = "\n", row.names = FALSE, col.names = TRUE)

	#print(allgene[i,4],max.levels = 0)

	cat(allgene[i,4])
	cat("\t")
	cat(as.character(allgene[i,6]))
	cat("\t")
	cat(as.character(allgene[i,1]))
	cat("\t")
	cat(allgene[i,7])
	cat("\t")
	cat(allgene[i,8])
	cat("\t")
	cat(nrow(data))
	cat("\t")
	cat(mean(data[,7] <= allgene[i,7]))
	cat("\t")
	cat(min(data[,7]))
	cat("\t")
	cat(max(data[,7]))
	cat("\t")
	cat(mean(data[,8] <= allgene[i,8]))
	cat("\t")
	cat(min(data[,8]))
	cat("\t")
	cat(max(data[,8]))
	cat("\n")

}
