
#=======================================================================
# step 0
#=======================================================================
unzip data/data_files.zip -d data/

#=======================================================================
# step 1
#=======================================================================
perl perl/dgap_features_check.pl \
  -f data/non_coding_DGAP_positions.bed \
  -i data/HI_Predictions_Version3.bed \
  -g data/GRCh37.p13_ensembl_genes.txt \
  -c data/hESC_hg37_domains.bed \
  -n 3000000 \
  -o .

#=======================================================================
# step 2
#=======================================================================
perl perl/enh_promoter_disruption_checker_DGAP.pl \
  -f data/non_coding_DGAP_positions.bed \
  -d data/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8 \
  -a 3000000 \
  -o DHS_promoter_broken_DGAP.txt
  

#=======================================================================
# step 3
#=======================================================================

Rscript R/combine_breakpoints_and_phenotypes.R \
  data/HPO_distal_cases.txt \
  data/non_coding_DGAP_positions.bed \
  6000000 \
  breakpoint_window_with_HPO.6MB_win.bed
  
java -jar bin/phenomatch.jar \
  -i breakpoint_window_with_HPO.6MB_win.bed \
  -g data/knownGene.txt.entrez_id.tab.unique \
  -O data/hp.obo \
  -a data/ALL_SOURCES_TYPICAL_FEATURES_genes_to_phenotype.txt  \
  -o breakpoint_window_with_HPO.6MB_win.bed.phenomatch

# java -jar bin/phenomatch.jar \
#   -i breakpoint_window_with_HPO.6MB_win.bed \
#   -g data/knownGene.txt.entrez_id.tab.unique \
#   -O data_raw/hp.obo \
#   -a data_raw/genes_to_phenotype.txt \
#   -o breakpoint_window_with_HPO.6MB_win.bed.phenomatch



#=======================================================================
# step 4
#=======================================================================
Rscript --vanilla R/get_percentiles_DGAP_all.r \
  breakpoint_window_with_HPO.6MB_win.bed.phenomatch.overlapped_genes.txt \
  > percentiles_6Mb_pheno_maxpheno.txt
  
#=======================================================================
# step 5
#=======================================================================
perl perl/dgap_final_table_maker.pl \
  -f HI_list_DGAP.txt \
  -c data/hESC_hg37_domains.bed \
  -d DHS_promoter_broken_breakpDGAP.txt \
  -h data/ClinGen_haploinsufficiency_gene.bed \
  -t data/ClinGen_triplosensitivity_gene.bed \
  -m percentiles_6Mb_pheno_maxpheno.txt \
  -o DGAP_table_summary.txt
