# position_effect
**Prediction of Position Effects of Apparently Balanced Human Chromosome Rearrangements**

The scirpts in this repository document the computational analysis in the followinng manuscript:
Zepeda-Mendoza and Ibn-Salem at al. 2017 "Computational Prediction of Position Effects of
Apparently Balanced Human Chromosome Rearrangements" 

To analyze your own genomic regions for position effects, please follow the steps in order. All included datasets are in hg19 version.

## STEP 1 
Run the dgap_features_check.pl script with the following parameters:

```
perl dgap_features_check.pl \
  -f GENOMIC_REGIONS \
  -i HI_SCORES_FILE \
  -g ENSEMBL_GENE_FILE_SUMMARY \
  -c HIC_DOMAINS \
  -n WINDOW_BP_SIZE \
  -o OUTPUT_FOLDER_PATH
```


**Example:**
```
perl dgap_features_check.pl \
  -f non_coding_DGAP_positions.bed \
  -i HI_Predictions_Version3.bed \
  -g GRCh37.p13_ensembl_genes.txt \
  -c hESC_hg37_domains.bed \
  -n 3000000 \
  -o ~/Desktop/Position_Effects_Analysis
```
The GENOMIC_REGIONS file needs to be in a .bed format with columns chr, start, end, case_id. 
This script generates two files, HiC_list_DGAP.txt (a list of overlapped Hi-C domains by the rearrangement breakpoints) and HI_list_DGAP.txt. 
The latter is a summary of the overlapped analysis windows, Hi-C domains, genes within the analysis windows with HI information, as well as gene type, description, and phenotype associations. 
The output column structure for this file is the following: 

+ 0	Breakpoint_ID	
+ 1	chr	
+ 2	start	
+ 3	end	
+ 4	win_start	
+ 5	win_end	
+ 6	Hi_domain	
+ 7	HiC_chr	
+ 8	HiC_start	
+ 9	HiC_end	
+ 10	HI_gene_chr	
+ 11	HI_gene_start	
+ 12	HI_gene_end	
+ 13	HI_gene_name	
+ 14	LOD	
+ 15	HI_prob	
+ 16	HI_ensembl_gene_ID	
+ 17	HI_chr	
+ 18	Gene_start	
+ 19	gene_end	
+ 20	HGNC_ID	
+ 21	HGNC_symbol	
+ 22	Gene_type	
+ 23	WikiGene_Description	
+ 24	Description	
+ 25	Phenotype_description

You can replace the haploinsufficiency file (HI_Predictions_Version3.bed, Huang et al., 2010), the gene analysis file (GRCh37.p13_ensembl_genes.txt, ENSEMBL GRCh37) and the HiC domains file (hESC_hg37_domains.bed, Dixon et al., 2012) with more recent versions or other custom files as long as the column structure is maintained.

## STEP 2 
Asses the number of DHS enhancer/promoter interactions disrupted by the rearrangement breakpoints

```
perl enh_promoter_disruption_checker_DGAP.pl \
  -f GENOMIC_REGIONS \
  -d DHS_ENH_PROMOTER_FILE \
  -a WINDOW_BP_SIZE \
  -o OUTPUT_FILENAME
```

**Example:**

```
perl enh_promoter_disruption_checker_DGAP.pl \
  -f non_coding_DGAP_positions.bed \
  -d genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8 \
  -a 3000000 \
  -o DHS_promoter_broken_DGAP.txt
```

This script will output the DHS enhancer/promoter contacts affected by the rearrangement breakpoints as well as a summary of the number of the disrupted contacts per analyzed case. The file used in this analysis was obtained from Andresson et al., 2014. The first file will be used for building the final analysis table for the analyzed regions (see Step 4).

## STEP 3 
Here phenomatch scores are calculated for genes located close to rearrangment breakpoints.

### a) Combine breakpoints and phenotype annotation in one file

```
Rscript combine_breakpoints_and_phenotypes.R \
  SUBJECT_TO_PHENOTYPES \
  BREAKPOINTS_BED \
  OUTPUT_FILE
```
This scripts takes two files as input:
`SUBJECT_TO_PHENOTYPES` is a tab delimited file with two columns and column headers:

+ **ID** the subject ID
+ **HPO** a single HPO term (e.g. HP:0000526)

Subjects with more than one HPO annotation have multiple lines in the file with the same ID in the first column

The second input file `BREAKPOINTS_BED` is a BED file wiht breakpoint locations.
In the fourth column should contain the patient ID with additinaol breakpoint identifiers separated by undersocre `_`. For example `DGAP111_A` and `DGAP111_B`.

The output file `OUTPUT_FILE` contains all breakpoints with an additional column with a comma-separated list of phenotypes.

**Example:**

```
Rscript combine_breakpoints_and_phenotypes.R \
  HPO_distal_cases.txt \
  non_coding_DGAP_positions.bed \
  breakpoint_window_with_HPO.bed
```

### b) Compute phenomatch scores

```
java -jar bin/topodombar.commandline-0.0.1-SNAPSHOT-jar-with-dependencies.jar  \
                  -i INPUT_FILE \
                  -d DOMAINS \
                  -g GENES \
                  -O PHENOTYPE_ONTOLOGY
                  -a ANNOTATION_FILE -o OUTPUT_FILE \
                  -e ENHANCERS \
                  -o OUTPUT_FILE

```

Use `-h` option for more usage information and other options including permutation anlaysis.
The calculation of phenomatch score is described in [Ibn-Salem et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0423-1) and the source code of the topodombar tool is available [here](https://github.com/charite/topodombar/tree/java-implementation).

**Example:**

```
java -jar bin/topodombar.commandline-0.0.1-SNAPSHOT-jar-with-dependencies.jar  \
    -i breakpoint_window_with_HPO.bed.6MB_win.bed \
    -d hESC_hg37_domains.bed  \
    -g knownGene.txt.entrez_id.tab.unique \
    -O hp.obo \
    -a ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt \
    -e Fetal_Brain.tab \
    -o breakpoint_window_with_HPO.bed.6MB_win.bed.annotated.out
```

The HPO files `hpo.obo`and `ALL_SOURCES_TYPICAL_FEATURES_genes_to_phenotype.txt` can be downloaded from the following URLs:

+ http://purl.obolibrary.org/obo/hp.obo
+ http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt

The tool will create several output files with the same prefix. The ouput file `breakpoint_window_with_HPO.bed.6MB_win.bed.annotated.out.overlapped_genes.txt` contains phenomatch scores.

## STEP 4:
You can calculate the percentile for the phenomatch and max_phenomatch scores by using the R script get_percentiles_DGAP_all.r

You can run the script with:

```
Rscript --vanilla get_percentiles_DGAP_all.r PHENO_FILE > OUTPUT_FILENAME
```

**Example:**
```
Rscript --vanilla get_percentiles_DGAP_all.r \
  breakpoint_window_with_HPO.bed.6MB_win.bed.annotated.out.overlapped_genes.txt \
  > percentiles_6Mb_pheno_maxpheno.txt
```


## STEP 5 
Summarize results

```
perl dgap_final_table_maker.pl \
  -f TAD_ANALYSIS_FEATURES_FILE \
  -c hESC_hg37_domains.bed \
  -d DHS_ENH_PROMOTER_DISRUPTED_CONTACTS_FILE \
  -h CLINGEN_HAPLO_FILE \
  -t CLINGEN_TRIPLO_FILE \
  -m PHENOMATCH_PERCENTILES_FILE \
  -o OUTPUT_FILENAME
```

**Example:**
```
perl dgap_final_table_maker.pl \
  -f HI_list_DGAP.txt \
  -c hESC_hg37_domains.bed \
  -d DHS_promoter_broken_breakpDGAP.txt \
  -h ClinGen_haploinsufficiency_gene.bed \
  -t ClinGen_triplosensitivity_gene.bed \
  -m percentiles_6Mb_pheno_maxpheno.txt \
  -o DGAP_table_summary.txt
```

The TAD_ANALYSIS_FEATURES_FILE is the HI_list_DGAP.txt file from step 1. The DHS_ENH_PROMOTER_DISRUPTED_CONTACTS_FILE is the output file of step 2 (not the summary!). Running this script will produce an additional section from the features analyzed in step 1, and analyzes the genes' inclusion in the 6Mb (or other sized) windows, the 2Mb window, the gene inclusion within the breakpoint TAD, whether or not predicted enh/promoter contacts (from Anderson et al) were disrupted for the gene, and finally a summary of the genes phenomatch and max_phenomatch scores with their corresponding percentile values within the analyzed dataset (output of step 4). Values of 0 and 1 indicate inclusion within the analyzed regions (6Mb, 2Mb, TADs) or presence of disrupted contacts (enh/promoter).

The output file has the following columns:

+ Breakpoint_ID
+ chr
+ start
+ end
+ win_start
+ win_end
+ Hi_domain
+ HiC_chr
+ HiC_start
+ HiC_end
+ HI_gene_chr
+ HI_gene_start
+ HI_gene_end
+ HI_gene_name
+ LOD
+ HI_prob
+ HI_ensembl_gene_ID
+ HI_chr
+ Gene-start
+ gene_end
+ HGNC_ID
+ HGNC_symbol
+ Gene_type
+ WikiGene_Description
+ Description
+ Phenotype_description
+ 6Mb
+ 2Mb
+ TAD
+ DHS_promoter
+ Haploinsufficiency_score
+ Triplosensitivity_score
+ Phenomatch_score
+ MaxPhenoScore
+ Pheno_percentile
+ MaxPheno_percentile

Because ClinGen haplo/triplo-sensitivity scores are indicated with ranges from 0 (not evidence) to 30 (known), we ask the program users to give weights to the values they want to analyze (i.e. if you want to consider only score values from 2 and above (recommended), give those regions a final table value of 1. You can either do this in the ClinGen file by susbstituting the desires scores with 1 and making everything else a 0, or by adding an extra column to the excel file and making this change with the IF selection formula). The same applies to the phenomatch and max_phenomatch percentiles (i.e. if you want to include only the top quartile, assign a 1 to everything with >=0.75 percentile value). 

Finally, after adding these 0 or 1 selection values to the haplo, triplo, and phenomatch percentiles, you can either use the PERC+DHS+TAD+HAPLO+TRIPLOor the PERC+DHS+2Mb+HAPLO+TRIPLO ranking criteria similar to Table S4A from the manuscript. Please refer to this table if you have any further questions, the fields are self-explanatory and compliment what was described in this section. If still in doubt, please contact Cinthya Zepeda (cinthya.zepeda.m@gmail.com).

You can replace the haploinsufficiency file (ClinGen_haploinsufficiency_gene.bed, ClinGen), the triplosensitivity file (ClinGen_triplosensitivity_gene.bed, ClinGen) and the HiC domains file (hESC_hg37_domains.bed, Dixon et al., 2012) with more recent versions or other custom files as long as the column structure is maintained.


