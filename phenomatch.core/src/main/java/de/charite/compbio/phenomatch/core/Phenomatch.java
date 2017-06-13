/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package de.charite.compbio.phenomatch.core;

import annotation.AnnotateCNVs;
import static annotation.AnnotateGenes.addGeneSymbol;
import genomicregions.CNV;
import genomicregions.Gene;
import genomicregions.GenomicElement;
import genomicregions.GenomicSet;
import io.GeneSymbolParser;
import io.TabFileParser;
import io.TabFileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.StringUtils;
import permutation.PermutedGenePhenotypes;
import phenotypeontology.PhenotypeData;


/**
 * This is the main program class of the phenomatch project.
 * 
 * 
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class Phenomatch {
    
        
    /**
     * The phenotype ontology instance (can be HPO or Uberpheno).
     */
    
    /** access to the phenotype ontology and gene to phenotype associations. */
    private PhenotypeData phenotypeData;

    /** Input copy number variations (CNVs) to annotate with effect mechanism*/ 
    private GenomicSet<CNV> cnvs;
    
    /** Topological domains regions */
    private GenomicSet<GenomicElement> domains;

    /** Genes */
    private GenomicSet<Gene> genes;
    
    /** file paths to input and output files **/
    private String cnvPath;
    
    private String genesPath;
    
    private String domainPath;
    
    private String ontologyPath;
    
    private String annotationPath;

    private Integer genePermutations;    

    private String outputPath;
    
    
    /**
     * Constructor for an instance of the {@link Phenomatch} program with input 
     * data and parameters.
     * 
     * @param argMap a {@link Map} holding the input data and parameters
     * @throws IOException if files cannot be read or write
     */
    public Phenomatch(Map<String, Object> argMap) throws IOException{
        
        // get the individual values
        this.cnvPath = (String) argMap.get("input_file");
        this.genesPath = (String) argMap.get("genes");
        this.ontologyPath = (String) argMap.get("phenotype_ontology");
        this.annotationPath = (String) argMap.get("annotation_file");
        this.outputPath = (String) argMap.get("output_file");
        
        // parse optional arguments:
        this.domainPath = (String) argMap.get("domains");
        this.genePermutations = (Integer) argMap.get("permut_genes");
        
        // read the phenotype ontology
        this.phenotypeData = new PhenotypeData(this.ontologyPath, this.annotationPath);
        System.out.println("[INFO] Ontology and annotation table were parsed.");

        ////////////////////////////////////////////////////////////////////////
        //  CNVs
        ////////////////////////////////////////////////////////////////////////

        // read CNV data from input file:
        TabFileParser cnvParser = new TabFileParser(cnvPath);        
        cnvs = cnvParser.parseCNVwithPhenotypeAnnotation(this.phenotypeData);
        
        ////////////////////////////////////////////////////////////////////////
        //  Domains and Boundaries
        ////////////////////////////////////////////////////////////////////////

        // read topological domain regions and compute boundaries from it.
        
        if(this.domainPath != null){
            
            TabFileParser domainParser = new TabFileParser(this.domainPath);
            // parse topological domains 
            domains = domainParser.parse();
        }
        
        ////////////////////////////////////////////////////////////////////////
        //  Genes
        ////////////////////////////////////////////////////////////////////////
        genes = new TabFileParser(genesPath).parseGeneWithTerms(phenotypeData);        
        
        // add GeneSymbol to genes
        HashMap<String, String> entrezToSymbol = new GeneSymbolParser(annotationPath).parseEntrezToSymbol();
        addGeneSymbol(genes, entrezToSymbol);
        
    }

    /**
     * Run permutation analysis to get significance of actual data
     */
    public void runPermutations() throws IOException{

        
        // save original CNVs and gene to phenotype mapping
        GenomicSet<CNV> orgCNVs = this.cnvs;
                
        if(this.genePermutations > 0){
            
            // get back the original cnvs
            this.cnvs = orgCNVs;

            // permutate gene phenotypes:
            analysePermutedGenePhenotypes(this.genePermutations);

        }
    }
    
    /**
     * Runs the entire analysis.
     */
    public void runAnalysis(){
        
        // annotate CNVs with genes that are completely overlapped by the CNV
        AnnotateCNVs.annotateOverlappedGenes(cnvs, genes);
        
        if (this.domainPath != null){
            // annotate CNVs with genes that are within TAD that have any overlap with the CNV
            AnnotateCNVs.annotateGenesInOverlapTADs(cnvs, domains, genes);
        }
    }
    
    /**
     * Runs permutations of gene phenotypes and report background rates of 
     * effect mechanism classes.
     * @param permutations number of permutations
     */
    private void analysePermutedGenePhenotypes(Integer permutations) throws IOException{
        
        // initialize output lines for overlapped genes with phenoMatchscore:
        ArrayList<String> outLinesOl = new ArrayList<String>();
        ArrayList<String> outLinesOlTAD = new ArrayList<String>();

        // run N times the actual permutations and analyse the datachr2chr2
        for (int i=0; i<permutations; i++){
            
            // permutate the phenotypes of the genes randomly
            PhenotypeData permPhenotypeData = PermutedGenePhenotypes.permuteGenePhenotypes(this.phenotypeData);
            this.phenotypeData = permPhenotypeData;

            try {
                this.genes = new TabFileParser(this.genesPath).parseGeneWithTerms(this.phenotypeData);
            } catch (IOException ex) {
                Logger.getLogger(Phenomatch.class.getName()).log(Level.SEVERE, null, ex);
            }
            
            // rerun the whole analysis
            // TODO run just the needed steps of the analysis
            this.runAnalysis();
            
            // iterate over all CNVs and append gene output lines
            for (CNV cnv : cnvs.values() ){
            
                outLinesOl.addAll(cnv.getOverlappedGenesOutputLine(phenotypeData, cnv.getGenesInOverlap()));
                
                if(this.domainPath != null){
                    outLinesOlTAD.addAll(cnv.getOverlappedGenesOutputLine(phenotypeData, cnv.getGenesInOverlapTADs()));
                }
            }
        }

        // write gene output lines to output file:
        TabFileWriter geneOutWriterOl = new TabFileWriter(outputPath + 
                ".permutGenePT_" + this.genePermutations + ".overlapped_genes.txt");
        geneOutWriterOl.writeLines(outLinesOl);
        System.out.println("[INFO] Wrote all overlapped genes from "
                + "permutated gene to phenotype annotation to output file "
                + "'"+this.outputPath+".permutGenePT_" + this.genePermutations + ".overlapped_genes.txt'.");

        if (this.domainPath != null){
            // write gene output lines to output file:
            TabFileWriter geneOutWriterOlTAD = new TabFileWriter(outputPath + 
                    ".permutGenePT_" + this.genePermutations + ".genes_in_overlapped_TADs.txt");
            geneOutWriterOlTAD.writeLines(outLinesOlTAD);
            System.out.println("[INFO] Topodombar: Wrote all genes in overlapped TADs from "
                    + "permutated gene to phenotype annotation to output file "
                    + "'"+this.outputPath+".permutGenePT_" + this.genePermutations + ".genes_in_overlapped_TADs.txt'.");
        }
    }
    
    
    /**
     * writes the cnvs and each ovelrapped gene per line with pheno match score 
     * to output file.
     * 
     * @throws IOException 
     */
    public void writeGeneOutput() throws IOException{
        
        
        TabFileWriter<CNV> outWriterOl = new TabFileWriter<CNV>(this.outputPath + ".overlapped_genes.txt");

        ArrayList<String> outLines = getOverlappedGenesOutput(this.cnvs, this.phenotypeData);
        outWriterOl.writeLines(outLines);
        System.out.println("[INFO] Wrote all overlapped genes to output file '"+this.outputPath+".overlapped_genes.txt'.");

        if (this.domainPath != null){
            TabFileWriter<CNV> outWriterOlTAD = new TabFileWriter<CNV>(this.outputPath + ".genes_in_overlapped_TADs.txt");

            ArrayList<String> outLinesTADs = getGenesInOverlappedTADs(this.cnvs, this.phenotypeData);
            outWriterOlTAD.writeLines(outLinesTADs);
            System.out.println("[INFO] Wrote all genes in overlapped TADs to output file '"+this.outputPath+".genes_in_overlapped_TADs.txt'.");
        }
    }    

    /**
     * Get lines for each overlapped gene per line for output file.
     * @param cnvs CNVs to be written to the output file
     * @param phenotypeData 
     * @return ArrayList of strings which are tab-separted output liens per gene
     */
    public ArrayList<String> getOverlappedGenesOutput(GenomicSet<CNV> cnvs, PhenotypeData phenotypeData) {
        
        // to get all output lines from the GenomicSet object:
        ArrayList<String> outLines = new ArrayList<>();

        // sort CNVs by there effect mechanism class
        ArrayList<CNV> sortedCNVs = new ArrayList<>(cnvs.values());
        
        for (CNV c : sortedCNVs){
            outLines.addAll(c.getOverlappedGenesOutputLine(phenotypeData, c.getGenesInOverlap()));
        }
        // put togeter all annotation string separated by TAB
        String headerLine = GenomicElement.getOutputHeaderLine()
            + "\t" 
            + StringUtils.join(new String[]{"phenotypes", "gene_symbol", 
                "phenoMatchScore", 
                "maxPhenoMatchScore", 
                "maxPatientMatchTerm", 
                "maxGeneMatchTerm", 
                "maxCommonTerm", 
                "allPatientMatchTerm",
                "allGeneMatchTerm",
                "allCommonTerm",
                "allMatchScore"                
                }, '\t');            

        // add header to beginning of output lines
        outLines.add(0, headerLine);
        
        // write all lines to the output file.
        return(outLines);
    }

    /**
     * Retunrs CNVs with each overlapped gene as line for output file.
     * @param cnvs CNVs to be written to the output file
     * @param phenotypeData 
     * @return ArrayList of Strings which are tab-separated output liens  
     */
    public ArrayList<String> getGenesInOverlappedTADs(GenomicSet<CNV> cnvs, PhenotypeData phenotypeData) {
        
        // to get all output lines from the GenomicSet object:
        ArrayList<String> outLines = new ArrayList<String>();

        // sort CNVs by there effect mechanism class
        ArrayList<CNV> sortedCNVs = new ArrayList<CNV>(cnvs.values());

        for (CNV c : sortedCNVs){
            outLines.addAll(c.getOverlappedGenesOutputLine(phenotypeData, c.getGenesInOverlapTADs()));
        }

        // put togeter all annotation string separated by TAB
        String headerLine = GenomicElement.getOutputHeaderLine()
            + "\t" 
            + StringUtils.join(new String[]{"phenotypes", "gene_symbol", 
                "phenoMatchScore", 
                "maxPhenoMatchScore", 
                "maxPatientMatchTerm", 
                "maxGeneMatchTerm", 
                "maxCommonTerm", 
                "allPatientMatchTerm",
                "allGeneMatchTerm",
                "allCommonTerm",
                "allMatchScore"
                    }, '\t');            

        // add header to beginning of output lines
        outLines.add(0, headerLine);
        
        // return all lines.
        return(outLines);
    }

}
