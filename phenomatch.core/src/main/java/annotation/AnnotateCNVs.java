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

package annotation;

import genomicregions.CNV;
import genomicregions.Gene;
import genomicregions.GenomicElement;
import genomicregions.GenomicSet;
import phenotypeontology.PhenotypeData;

/**
 * This class implements functionality to annotate CNVs and predict their most
 * likely effect mechanisms.
 *
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class AnnotateCNVs {
    
    /**
     * This function annotates the input CNVs with overlapped topological domain
     * boundaries and genes together with the corresponding phenogram score of 
     * these genes overlapped by the CNV.
     * 
     * @param cnvs
     * @param boundaries
     * @param genes
     * @param enhancers
     * @param phenotypeData 
     */
    

    /**
     * Annotates CNVs with overlapping boundaries.
     * 
     * @param cnvs  copy number variations (CVNs)
     * @param boundaries Topological domain boundaries
     */
    public static void boundaryOverlap(GenomicSet<CNV> cnvs, GenomicSet<GenomicElement> boundaries){
        
        // iterate over all CNVs:
        for (CNV cnv : cnvs.values()){
            
            GenomicSet<GenomicElement> overlap = boundaries.completeOverlap(cnv);
            cnv.setBoundaryOverlap( overlap );
        }
    }
    
    /**
     * Annotates all input CNVs with all genes that have any overlap with the CNV.
     * For each {@link CNV} object the variable {@link CNV.geneOverlap} is filled 
     * with a {@link GenomicSet} of {@link Gene} objects.
     * 
     * @param cnvs CNVs that should be annotated
     * @param genes Set of genes
     */
    public static void annotateOverlappedGenes(GenomicSet<CNV> cnvs, GenomicSet<Gene> genes){
        
        // iterate over all CNVs:
        for (CNV cnv : cnvs.values()){
            
            GenomicSet<Gene> overlap = genes.anyOverlap(cnv);
            cnv.setGenesInOverlap( overlap );
            
        }
    }
    
    /**
     * Annotates all input CNVs with genes that are within TADs that are overlapped with the CNV.
     * For each {@link CNV} object the variable {@link CNV.genesInOverlapTADs} is filled 
     * with a {@link GenomicSet} of {@link Gene} objects.
     * @param cnvs CNVs that should be annotated
     * @param domains set of TADs
     * @param genes set of genes
     */
    public static void annotateGenesInOverlapTADs(GenomicSet<CNV> cnvs, GenomicSet<GenomicElement> domains, GenomicSet<Gene> genes){

        // iterate over all CNVs:
        for (CNV cnv : cnvs.values()){
            
            // get subset of TADs overlapping with the CNV
            GenomicSet<GenomicElement> overlapTADs = domains.anyOverlap(cnv);
            
            // initialize set of genes
            GenomicSet<Gene> overlapGenes = new GenomicSet<Gene>();
            
            // iterate over all overlpping TADs and update the set of genes
            for (GenomicElement tad: overlapTADs.values()){
                overlapGenes.putAll( genes.anyOverlap(tad) );                
            }

            // set the genes in overlapping TAD annotation
            cnv.setGenesInOverlapTADs( overlapGenes );
            
        }

    }
    
   
    /**
     * Compute the phenogram score for genes overlapped by the input {@link CNV}s.
     * It writes the memeber variables {@link CNV.overlapPhenogramScore} in each {@link CNV} object
     * Note, this function assumes that the CNVs are annotated with overlped and
     * adjacent genes by the functions {@link annotateOverlappedGenes}.
     * 
     * @param cnvs CNVs for which the phenogram score should be calculated.
     * @param phenotypeData the phenotype ontology
     */
    public static void overlapPhenogramScore(GenomicSet<CNV> cnvs, PhenotypeData phenotypeData){
        
        for (CNV cnv: cnvs.values()){
            // overlap PhenogramScore
            cnv.setOverlapPhenogramScore( 
                phenotypeData.phenoGramScore( cnv.getPhenotypes(), cnv.getGenesInOverlap()) 
            );
        }
    }
    
    /**
     * Compute the phenogram score for genes in the left and right adjacent regions.
     * It writes the memeber variables {@link CNVleftAdjacentPhenogramScore} and 
     * {@link CNVrightAdjacentPhenogramScore} in each {@link CNV} object
     * Note, this function assumes that the CNVs are annotated adjacent genes by 
     * the function {@link annotateAdjacentGenes}.
     * 
     * @param cnvs CNVs for which the phenogram score should be calculated.
     * @param phenotypeData the phenotype ontology
     */
    public static void adjacentPhenogramScore(GenomicSet<CNV> cnvs, PhenotypeData phenotypeData){
        
        for (CNV cnv: cnvs.values()){
        
            // left adjacent PhenogramScore:
            cnv.setLeftAdjacentPhenogramScore( 
                phenotypeData.phenoGramScore( cnv.getPhenotypes(), cnv.getGenesInLeftRegion()) 
            );

            // right adjacent PhenogramScore:
            cnv.setRightAdjacentPhenogramScore(
               phenotypeData.phenoGramScore( cnv.getPhenotypes(), cnv.getGenesInRightRegion())
            );
        }
    }
    
    /**
     * Annotates all input CNVs with enhancers in the left and right adjacent regions of the CNV.
     * For each {@link CNV} object the variable {@link CNV.enhancersInLeftRegion} 
     * and {@link CNV.enhancersInRightRegion} is filled with a 
     * {@link GenomicSet} of {@link GenomicElemnt} objects representing the enhancers.
     * Note, his function assumes pre-defined the the adjacent regions 
     * (e.g by the function {@link defineAdjacentRegionsByDomains}).  
     * Note, this function requires that an enhancer lies completely in the adjacent regions,
     * that is the adjacent regions overlaps completely the entire enhancer.
     * 
     * @param cnvs CNVs that should be annotated
     * @param enhancers Set of enhancers
     */    
    public static void annotateAdjacentEnhancers(GenomicSet<CNV> cnvs, GenomicSet<GenomicElement> enhancers){

        // iterate over all input CNVs:
        for (CNV cnv : cnvs.values()){

            // get all enhancers in the left adjacent regions
            GenomicSet<GenomicElement> leftEnhancers = enhancers.completeOverlap(cnv.getLeftAdjacentRegion());
            cnv.setEnhancersInLeftRegion( leftEnhancers );

            // get all enhancers in the right adjacent regions
            GenomicSet<GenomicElement> rightEnhancers = enhancers.completeOverlap(cnv.getRightAdjacentRegion());
            cnv.setEnhancersInRightRegion( rightEnhancers );
            
        }
        
    }
    

    /**
     * This function defines the overlapped inner-domain regions by the CNV.
     * There is a region for each CNV breakpoint that is defined from the break-point
     * to the next boundary that is overlapped by the cnv.
     * <pre>
Domains:               /'''''''''\   /'''\  /''''''''''''\
Boundary:                         ---     --       
CNV:                         ======================
overlappedDomainRegions:     *****          *******
     </pre>
     * These regions are needed to interpret tandem duplications and inversions for 
     * enhancer adoption mechanisms (e.g relevant enhancer in left overlap region
     * might come close and interact with a duplicated gene in the right overlapped region).
     * 
     * @param cnvs
     * @param domains 
     */
    public static void defineOverlappedDomainRegions(GenomicSet<CNV> cnvs, GenomicSet<GenomicElement> domains){
        
        for (CNV cnv: cnvs.values()){
            
            int start = cnv.getStart();
            int end = cnv.getEnd();
            String chr = cnv.getChr();

            // check if CNV ovelaps at leaste one boundary element
            if (cnv.hasBoundaryOverlap()){
                
                // get domain regions underling the left (start) and right (end) borders of the CNV
                GenomicSet<GenomicElement> leftDomains = domains.anyOverlap(new GenomicElement(chr, start-1, start, "cnvStart"));
                GenomicSet<GenomicElement> rightDomains = domains.anyOverlap(new GenomicElement(chr, end, end+1, "cnvEnd"));

                // if left CNV border lies not in a domain region (but in boundary or unorganized chromatin),
                // an zero length region will be defined.
                if (leftDomains.isEmpty()){
                    cnv.setLeftOverlappedDomainRegion(new GenomicElement(chr, start, start, "leftOverlapped"));
                }else{

                    // get end of overlapped region as end of the domain overlapped by the CNV start
                    GenomicElement leftDomain = leftDomains.values().iterator().next();
                    int leftOverlapRegionEnd = leftDomain.getEnd();

                    // construct element for the left adjacent regions
                    cnv.setLeftOverlappedDomainRegion(new GenomicElement(chr, cnv.getStart(), leftOverlapRegionEnd, "leftOverlapped"));
                }

                // same for the right site:
                if(rightDomains.isEmpty()){
                    // zero-lenght region
                    cnv.setRightOverlappedDomainRegion(new GenomicElement(chr, end, end, "rightOverlapped"));
                }else{

                    // get start of the query region as start of right adjacten domain
                    GenomicElement rightDomain = rightDomains.values().iterator().next();
                    int rightOverlapRegionStart = rightDomain.getStart();

                    // construct query element for the right adjacent regions
                    cnv.setRightOverlappedDomainRegion( new GenomicElement(chr, rightOverlapRegionStart, cnv.getEnd(), "rightOverlapped"));

                }
            
            // in case of no boundary overlap set default regions of zero size
            }else{
                cnv.setLeftOverlappedDomainRegion(new GenomicElement(chr, start, start, "leftOverlapped"));
                cnv.setRightOverlappedDomainRegion(new GenomicElement(chr, end, end, "rightOverlapped"));
            }
        }
    }
    
 
    
    
}
