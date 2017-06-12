/*
 * Copyright (c) 2014, Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

package genomicregions;
   
import annotation.AnnotateCNVs;
import io.Utils;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import ontologizer.go.Term;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils; // provides a join(iterable, char) function
import phenotypeontology.PhenotypeData;
import phenotypeontology.TermPair;
import static phenotypeontology.TermPair.toScore;

/**
 * This class implements a copy number variation (CNV) object. The CNV has a 
 * defined type and location in a reference genome. Furthermore it is associated to 
 * phenotypes observed in the patient that carries the CNV. Several members of 
 * this object can be used to annotate the CNV with respect to other genomic
 * elements and potential effect mechanisms. 
 * 
 * @author jonas
 */
public class CNV extends GenomicElement {
    
    /**
     * Type of CNV ("loss", "gain", "inversion"). This field can be later used to indicate 
     * more complex structural variations.
     */
    private final String type;
    
    /**
     * Phenotype terms of the CNV carrier as {@link HashSet} of {@link Term} objects.
     */
    private final HashSet<Term> phenotypes;

    
    /** Target term or phenotype category as single general HPO term ID. */    
    private final Term targetTerm;
    
    
    /** List of overlapping boundaries. */
    private GenomicSet<GenomicElement> boundaryOverlap;    
   
    /** List of overlapping genes (any overlap) */
    private GenomicSet<Gene> genesInOverlap;
    
    /** List of overlapping genes in overlapping TADs */
    private GenomicSet<Gene> genesInOverlapTADs;
    
    /** Adjacent genomic region on the left (5') site of the CNV. */
    private GenomicElement leftAdjacentRegion;
    
    /** Adjacent genomic region on the right (3') site of the CNV. */
    private GenomicElement rightAdjacentRegion;
    
    /** {@link GenomicSet} of {@link Genes} in the left adjacent region of this CNV. */
    private GenomicSet<Gene> genesInLeftRegion;
    
    /** {@link GenomicSet} of {@link Genes} in the right adjacent region of this CNV. */
    private GenomicSet<Gene> genesInRightRegion;
    
    /** {@link GenomicSet} of enhancers in the left adjacent region of this CVN. */
    private GenomicSet<GenomicElement> enhancersInLeftRegion;
    
    /** {@link GenomicSet} of enhancers in the right adjacent region of this CVN. */
    private GenomicSet<GenomicElement> enhancersInRightRegion;

    /** Overlapped genomic region in the domain overlapping the 3' end of the CNV */
    private GenomicElement leftOverlappedDomainRegion;
    
    /** Overlapped genomic region in the domain overlapping the 5' end of the CNV */
    private GenomicElement rightOverlappedDomainRegion;
    
    /** Phenogram score of all genes overlapped by the CNV. */
    private Double overlapPhenogramScore;
    
    /** Phenogram score of genes in the left adjacent region. */
    private Double leftAdjacentPhenogramScore;

    /** Phenogram score of genes in the right adjacent region. */
    private Double rightAdjacentPhenogramScore;
    
    /**
     * Possible CNV effect mechanism classes maped to the possible effect annotations.
     */
    private final  static HashMap<String, String []> effectMechansim2effects;
    static{
        effectMechansim2effects = new HashMap();
        effectMechansim2effects.put("TDBD", new String [] {"TDBD", "Mixed", "GDE", "NoData", "NA"});    
        effectMechansim2effects.put("newTDBD", new String [] {"TDBD", "Mixed", "GDE", "NoData", "NA"});    
        effectMechansim2effects.put("EA", new String [] {"EA", "Mixed", "GDE", "NoData", "NA"});    
        effectMechansim2effects.put("EAlowG", new String [] {"EAlowG", "Mixed", "GDE", "NoData", "NA"});    
        effectMechansim2effects.put("InvEA", new String [] {"EandGInvEA", "EnhancerInvEA", "GeneInvEA", "NoInvEA", "NA"});    
        /** 
         * Inversion enhancer removal (InvER).
         * Here the inversion removes an enhancer that is normally associated to phenotypically relevant gene
         * Or alternatively the gene gets inverted (moves to an other domain).
         */
        effectMechansim2effects.put("InvER", new String [] {"InvertedEnhancer", "InvertedGene", "NoInvER", "NA"});    
        effectMechansim2effects.put("TanDupEA", new String [] {"TanDupEA", "onlyGDE", "NoData", "NA"});    
    }
   
    /**
     * Constructor for CNV object.
     * Construct a {@link GenomicElement} and sets all {@link CNV} specific 
     * annotations to default values.
     * 
     * @param chr   Chromosome identifier
     * @param start Start coordinate (0-based)
     * @param end   End coordinate (0-based, half-open)
     * @param name  Name or ID of the CNV
     * @param type  CNV type (e.g. loss or gain)
     * 
     * @throws IllegalArgumentException 
     */
    public CNV(String chr, int start, int end, String name, String type) throws IllegalArgumentException {
        super(chr, start, end, name);
        
        this.overlapPhenogramScore = -1.0;
        this.genesInOverlap = new GenomicSet();
        this.boundaryOverlap = new GenomicSet();
        
        // set default values for annotations
        this.type = type;
        this.phenotypes = new HashSet<Term>();
        this.targetTerm = null;
        this.boundaryOverlap = new GenomicSet<GenomicElement>();
        this.genesInOverlap = new GenomicSet<Gene>();
        this.overlapPhenogramScore = -1.0;
        this.genesInLeftRegion = new GenomicSet<Gene>();
        this.leftAdjacentPhenogramScore = -1.0;
        this.genesInRightRegion = new GenomicSet<Gene>();
        this.rightAdjacentPhenogramScore = -1.0;
        this.enhancersInLeftRegion = new GenomicSet<GenomicElement>();
        this.enhancersInRightRegion = new GenomicSet<GenomicElement>();
        
    }
    
    /**
     * Constructor for {@link CNV} object.
     * Construct an CNV object by taking all annotations as arguments.
     * 
     * @param chr   Chromosome identifier
     * @param start Start coordinate (0-based)
     * @param end   End coordinate (0-based, half-open)
     * @param name  Name or ID of the CNV
     * @param type  CNV type (loss or gain)
     * @param phenotypes    List of HPO term IDs that represent the phenotypes used to annotate the patient carrying the CNV.
     * @param targetTerm A unspecific target term as a HPO term ID that is used to group patients in cohort or associate patients to tissues for which enhancer data his available.
     */
    public CNV(String chr, int start, int end, String name, 
                String type, HashSet<Term> phenotypes, Term targetTerm){

        // consturct an CVN object using the constructor of the {@link GenomicElement} super calss
        super(chr, start, end, name);
        
        this.overlapPhenogramScore = -1.0;
        this.genesInOverlap = new GenomicSet();
        this.boundaryOverlap = new GenomicSet();
        
        // add annotations
        this.type = type;
        this.phenotypes = phenotypes;
        this.targetTerm = targetTerm;        
        
        // set default annotaton for other fealds
        this.boundaryOverlap = new GenomicSet<GenomicElement>();
        this.genesInOverlap = new GenomicSet<Gene>();
        this.genesInOverlapTADs = new GenomicSet<Gene>();
        this.overlapPhenogramScore = -1.0;
        this.genesInLeftRegion = new GenomicSet<Gene>();
        this.leftAdjacentPhenogramScore = -1.0;
        this.genesInRightRegion = new GenomicSet<Gene>();
        this.rightAdjacentPhenogramScore = -1.0;
        this.enhancersInLeftRegion = new GenomicSet<GenomicElement>();
        this.enhancersInRightRegion = new GenomicSet<GenomicElement>();        

    }
    

    /**
     * Create a header line as TAB-separated string with all column identifiers 
     * for the simple output file format.
     * 
     * @return simple output file header as TAB separated column descriptions
     */
    public static String getSimpleOutputHeader(){
        
        String [] cnvColumns = new String[]{
                    "#chr",
                    "start",
                    "end",
                    "name",
                    "mostLikelyEffect", 
                    "gene", 
                    "enhancer",
                    "affectedGenes",
                };
        // add the effect mechanism columns and return a TAB-separated string.
        return StringUtils.join(cnvColumns, '\t');
    }

    /**
     * This function constructs  a ArrayList of {@link String} that represents 
     * output lines for each gene overlapped by this CNV.
     * 
     * @param phenotypeData a {@link PhenotypeData} object to calculate phenoMatch scores
     * @param genes set of genes to consider
     * @return a TAB-separated output line to write BED like files.
     */
    public ArrayList<String> getOverlappedGenesOutputLine(PhenotypeData phenotypeData, GenomicSet<Gene> genes){
        
        ArrayList<String> outLines = new ArrayList<>();
        
        //convert phenotpye terms to Strings
        HashSet<String> phenotypesIDs = new HashSet<>();
        for (Term t : this.phenotypes){
            phenotypesIDs.add(t.getIDAsString()); 
        }
        
        // For columns with multiple elements, separate them by semiclon ';'
        String phenotypeCol = StringUtils.join(phenotypesIDs, ';');

        // initialize maximal score per CNV
        Double maxScore = 0.0;
        
        // if CNV does overlap any gene add output line for each gene with score > 0
        if (!genes.isEmpty() ){
                        
            for (Gene g : genes.values()){


                String geneSymbol = g.getSymbol();
                //double geneScore = phenotypeData.phenoMatchScore(this.phenotypes, g);
                
                ArrayList<TermPair> termMatching = phenotypeData.phenoMatchScoreWithMatching(this.phenotypes, g);

                double maxGeneScore = 0.0;
                double sumGeneScore = 0.0;
                String maxPatientMatchTerms = "";
                String maxGeneMatchTerms = "";
                String maxLca = "";
                
                if(termMatching.size() > 0){
                    
                    TermPair maxPair = Collections.max(termMatching, TermPair.TERM_PAIR_SCORE_ORDER);
                    maxGeneScore = maxPair.getS();
                    for (TermPair tp: termMatching){
                        sumGeneScore += tp.getS();
                    }
                    maxPatientMatchTerms = maxPair.getPp().getIDAsString();
                    maxGeneMatchTerms = maxPair.getGp().getIDAsString();
                    maxLca = maxPair.getLca().getIDAsString();
                }

                String sumScoreStr = Utils.roundToString(sumGeneScore);
                String maxScoreStr = Utils.roundToString(maxGeneScore);

                // only if there is a score larger than zero output the gene
                if (maxGeneScore > 0){
                    
                    String allPatientMatchTerms = "";
                    String allGeneMatchTerms = "";
                    String allLca = "";
                    String allMatchScores = "";

                    String sep = "";
                    
                    for (TermPair tp: termMatching){
                        
                        if (allPatientMatchTerms.length() > 0){
                            sep=";";                            
                        }
                        allPatientMatchTerms += sep + tp.getPp().getIDAsString();
                        allGeneMatchTerms += sep + tp.getGp().getIDAsString();
                        allLca += sep + tp.getLca().getIDAsString();
                        allMatchScores += sep + Utils.roundToString(tp.getS());
                        
                    }

                    String [] cnvAnnotations = new String[]{
                            phenotypeCol, 
                            geneSymbol,
                            sumScoreStr,
                            maxScoreStr,
                            maxPatientMatchTerms,
                            maxGeneMatchTerms,
                            maxLca,
                            allPatientMatchTerms,
                            allGeneMatchTerms,
                            allLca,
                            allMatchScores
                        };

                    // put togeter all annotation string separated by TAB
                    String outLineGene = super.toOutputLine()
                        + "\t" 
                        + StringUtils.join(cnvAnnotations, '\t');            

                    // append to output lines
                    outLines.add(outLineGene);
                    
                    // upate maxScore
                    if (maxGeneScore > maxScore){
                        maxScore = maxGeneScore;
                    }
                }
            }
        }
        
        // if CNV does not overlap any gene or overlapped genes have score 0
        if( maxScore == 0.0){
            String outLineCNV = super.toOutputLine()
                    + "\t" 
                    + StringUtils.join(new String[]{
                        phenotypeCol, 
                        ".",
                        Utils.roundToString(0.0),
                        Utils.roundToString(0.0),
                        ".",
                        ".",
                        ".",
                        ".",
                        ".",
                        ".",
                        Utils.roundToString(0.0)
                    }, '\t');
            // append to output lines
            outLines.add(outLineCNV);
        }
        
        return outLines;
        
    }

    /**
     * Type of CNV ("loss", "gain" or "inversion"). This field can be later used to indicate 
     * more complex structural variations.
     * 
     * @return the type
     */
    public String getType() {
        return type;
    }

    /**
     * Phenotype terms of the {@link CNV} carrier as {@link HashSet} of {@link Term} objects.
     * @return 
     */
    public HashSet<Term> getPhenotypes() {
        return this.phenotypes;
    }

    /**
     * Target term or phenotype category as single general HPO term ID
     * @return the targetTerm
     */
    public Term getTargetTerm() {
        return this.targetTerm;
    }

    /**
     * {@link GenomicSet} of overlapping boundaries.
     * @return the boundaryOverlap
     */
    public GenomicSet<GenomicElement> getBoundaryOverlap() {
        return boundaryOverlap;
    }

    /**
     * {@link GenomicSet} of overlapping boundaries.
     * @param boundaryOverlap the boundaryOverlap to set
     */
    public void setBoundaryOverlap(GenomicSet<GenomicElement> boundaryOverlap) {
        this.boundaryOverlap = boundaryOverlap;
    }

    /**
     * Should be {@code true} if CNV overlaps a boundary element.
     * @return the hasBoundaryOverlap
     */
    public boolean hasBoundaryOverlap() {
        return ! boundaryOverlap.isEmpty();
    }

    /**
     * {@link GenomicSet} of overlapping {@link Gene}s (any overlap)
     * @return the genesInOverlap
     */
    public GenomicSet<Gene> getGenesInOverlap() {
        return genesInOverlap;
    }

    /**
     * {@link GenomicSet} of overlapping {@link Gene}s (any overlap)
     * @param genesInOverlap the genesInOverlap to set
     */
    public void setGenesInOverlap(GenomicSet<Gene> genesInOverlap) {
        this.genesInOverlap = genesInOverlap;
    }
    
    /**
     * {@link GenomicSet} of {@link Gene}s overlapping TADs that overlap the CNV.
     * @return 
     */
    public GenomicSet<Gene> getGenesInOverlapTADs() {
        return genesInOverlapTADs;
    }

    /**
     * {@link GenomicSet} of {@link Gene}s overlapping TADs that overlap the CNV.
     * @param genesInOverlapTADs 
     */
    public void setGenesInOverlapTADs(GenomicSet<Gene> genesInOverlapTADs) {
        this.genesInOverlapTADs = genesInOverlapTADs;
    }

    /**
     * Phenogram score of all genes overlapped by the CNV.
     * @return the overlapPhenogramScore
     */
    public Double getOverlapPhenogramScore() {
        return overlapPhenogramScore;
    }

    /**
     * Phenogram score of all genes overlapped by the CNV.
     * @param overlapPhenogramScore the overlapPhenogramScore to set
     */
    public void setOverlapPhenogramScore(Double overlapPhenogramScore) {
        this.overlapPhenogramScore = overlapPhenogramScore;
    }

    /**
     * Add a phenotype {@link Term} to the set of {@Term}s.
     * This method is just for creating artificial CNV objects for testing.
     * @param t 
     */
    public void addPhenotypeTerm(Term t) {
        this.phenotypes.add(t);
    }

    /**
     * Adjacent genomic region on the left (5') site of the CNV.
     * @return the leftAdjacentRegion
     */
    public GenomicElement getLeftAdjacentRegion() {
        return leftAdjacentRegion;
    }

    /**
     * Adjacent genomic region on the left (5') site of the CNV.
     * @param leftAdjacentRegion the leftAdjacentRegion to set
     */
    public void setLeftAdjacentRegion(GenomicElement leftAdjacentRegion) {
        this.leftAdjacentRegion = leftAdjacentRegion;
    }

    /**
     * Adjacent genomic region on the right (3') site of the CNV.
     * @return the rightAdjacentRegion
     */
    public GenomicElement getRightAdjacentRegion() {
        return rightAdjacentRegion;
    }

    /**
     * Adjacent genomic region on the right (3') site of the CNV.
     * @param rightAdjacentRegion the rightAdjacentRegion to set
     */
    public void setRightAdjacentRegion(GenomicElement rightAdjacentRegion) {
        this.rightAdjacentRegion = rightAdjacentRegion;
    }

    /**
     * Phenogram score of genes in the left adjacent region.
     * @return the leftAdjacentPhenogramScore
     */
    public Double getLeftAdjacentPhenogramScore() {
        return leftAdjacentPhenogramScore;
    }

    /**
     * Phenogram score of genes in the left adjacent region.
     * @param leftAdjacentPhenogramScore the leftAdjacentPhenogramScore to set
     */
    public void setLeftAdjacentPhenogramScore(Double leftAdjacentPhenogramScore) {
        this.leftAdjacentPhenogramScore = leftAdjacentPhenogramScore;
    }

    /**
     * Phenogram score of genes in the right adjacent region.
     * @return the rightAdjacentPhenogramScore
     */
    public Double getRightAdjacentPhenogramScore() {
        return rightAdjacentPhenogramScore;
    }

    /**
     * Phenogram score of genes in the right adjacent region.
     * @param rightAdjacentPhenogramScore the rightAdjacentPhenogramScore to set
     */
    public void setRightAdjacentPhenogramScore(Double rightAdjacentPhenogramScore) {
        this.rightAdjacentPhenogramScore = rightAdjacentPhenogramScore;
    }

    /**
     * {@link GenomicSet} of {@link Genes} in the left adjacent region of this CNV.
     * @return the genesInLeftRegion
     */
    public GenomicSet<Gene> getGenesInLeftRegion() {
        return genesInLeftRegion;
    }

    /**
     * {@link GenomicSet} of {@link Genes} in the left adjacent region of this CNV.
     * @param genesInLeftRegion the genesInLeftRegion to set
     */
    public void setGenesInLeftRegion(GenomicSet<Gene> genesInLeftRegion) {
        this.genesInLeftRegion = genesInLeftRegion;
    }

    /**
     * {@link GenomicSet} of {@link Genes} in the right adjacent region of this CNV.
     * @return the genesInRightRegion
     */
    public GenomicSet<Gene> getGenesInRightRegion() {
        return genesInRightRegion;
    }

    /**
     * {@link GenomicSet} of {@link Genes} in the right adjacent region of this CNV.
     * @param genesInRightRegion the genesInRightRegion to set
     */
    public void setGenesInRightRegion(GenomicSet<Gene> genesInRightRegion) {
        this.genesInRightRegion = genesInRightRegion;
    }

    /**
     * {@link GenomicSet} of enhancers in the left adjacent region of this CVN.
     * @return the enhancersInLeftRegion
     */
    public GenomicSet<GenomicElement> getEnhancersInLeftRegion() {
        return enhancersInLeftRegion;
    }

    /**
     * {@link GenomicSet} of enhancers in the left adjacent region of this CVN.
     * @param enhancersInLeftRegion the enhancersInLeftRegion to set
     */
    public void setEnhancersInLeftRegion(GenomicSet<GenomicElement> enhancersInLeftRegion) {
        this.enhancersInLeftRegion = enhancersInLeftRegion;
    }

    /**
     * {@link GenomicSet} of enhancers in the right adjacent region of this CVN.
     * @return the enhancersInRightRegion
     */
    public GenomicSet<GenomicElement> getEnhancersInRightRegion() {
        return enhancersInRightRegion;
    }

    /**
     * {@link GenomicSet} of enhancers in the right adjacent region of this CVN.
     * @param enhancersInRightRegion the enhancersInRightRegion to set
     */
    public void setEnhancersInRightRegion(GenomicSet<GenomicElement> enhancersInRightRegion) {
        this.enhancersInRightRegion = enhancersInRightRegion;
    }
    
    /**
     * Overlapped genomic region in the domain overlapping the 3' end of the CNV
     * @return the leftOverlappedDomainRegion
     * @see AnnotateCNVs.defineOverlappedDomainRegions
     */
    public GenomicElement getLeftOverlappedDomainRegion() {
        return leftOverlappedDomainRegion;
    }

    /**
     * Overlapped genomic region in the domain overlapping the 3' end of the CNV
     * @param leftOverlappedDomainRegion the leftOverlappedDomainRegion to set
     * @see AnnotateCNVs.defineOverlappedDomainRegions
     */
    public void setLeftOverlappedDomainRegion(GenomicElement leftOverlappedDomainRegion) {
        this.leftOverlappedDomainRegion = leftOverlappedDomainRegion;
    }

    /**
     * Overlapped genomic region in the domain overlapping the 5' end of the CNV
     * @return the rightOverlappedDomainRegion
     * @see AnnotateCNVs.defineOverlappedDomainRegions
     */
    public GenomicElement getRightOverlappedDomainRegion() {
        return rightOverlappedDomainRegion;
    }

    /**
     * Overlapped genomic region in the domain overlapping the 5' end of the CNV
     * @param rightOverlappedDomainRegion the rightOverlappedDomainRegion to set
     * @see AnnotateCNVs.defineOverlappedDomainRegions
     */
    public void setRightOverlappedDomainRegion(GenomicElement rightOverlappedDomainRegion) {
        this.rightOverlappedDomainRegion = rightOverlappedDomainRegion;
    }

    public void debugPrint(){
        System.out.println("DEBUG CNV: " + this.toString());
        System.out.println("DEBUG CNV type " + this.type);
        
        System.out.println("DEBUG CNV phenotype: " + this.phenotypes);
        System.out.println("DEBUG CNV targetTerm: " + this.targetTerm);
        
        System.out.println("DEBUG CNV boundaryOverlap: " + this.boundaryOverlap);
        System.out.println("DEBUG CNV genesInOverlap: " + this.genesInOverlap);
        System.out.println("DEBUG CNV overlapPhenogramScore: " + this.overlapPhenogramScore);
        
        System.out.println("DEBUG CNV leftAdjacentRegion: " + this.leftAdjacentRegion);
        System.out.println("DEBUG CNV enhancersInLeftRegion: " + this.enhancersInLeftRegion);
        System.out.println("DEBUG CNV genesInLeftRegion: " + this.genesInLeftRegion);
        System.out.println("DEBUG CNV leftAdjacentPhenogramScore: " + this.leftAdjacentPhenogramScore);
        
        System.out.println("DEBUG CNV rightAdjacentRegion: " + this.rightAdjacentRegion);
        System.out.println("DEBUG CNV enhancersInRightRegion: " + this.enhancersInRightRegion);
        System.out.println("DEBUG CNV genesInRightRegion: " + this.genesInRightRegion);
        System.out.println("DEBUG CNV rightAdjacentPhenogramScore: " + this.rightAdjacentPhenogramScore);
        
        System.out.println("DEBUG CNV leftOverlappedDomainRegion: " + this.leftOverlappedDomainRegion);
        System.out.println("DEBUG CNV rightOverlappedDomainRegion: " + this.rightOverlappedDomainRegion);
        
    }
}
