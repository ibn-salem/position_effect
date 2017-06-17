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
 
}
