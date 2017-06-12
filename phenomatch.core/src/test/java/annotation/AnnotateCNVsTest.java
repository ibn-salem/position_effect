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
import io.TabFileParser;
import io.TabFileParserTest;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import ontologizer.go.Term;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;
import phenotypeontology.PhenotypeData;
import toyexampledata.ExampleData;

/**
 *
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class AnnotateCNVsTest {
    
    // declare some variables to be knwon in this test
    private static GenomicSet<CNV> cnvs; 
    private static GenomicSet<CNV> exampleCNVs; 
    private static GenomicSet<GenomicElement> domains;
    private static GenomicSet<GenomicElement> boundaries;
    //private static GenomicElementSet<GenGenomicSet  public AnnotateCNVsTest() {
    private static GenomicSet<Gene> genes;
    
    private static ExampleData exampleData;
    
    @BeforeClass
    public static void setUpClass() throws IOException {

        // read sample CNVs:
        String cnvPath = TabFileParserTest.class.getResource("/sample_CNV_chr22.tab").getPath();
        TabFileParser cnvParser = new TabFileParser(cnvPath);
        cnvs = cnvParser.parseCNV();

        // read sample boundary elements
        String boundaryPath = TabFileParserTest.class.getResource("/sample_boundary_chr22.tab.addOne").getPath();
        TabFileParser boundaryParser = new TabFileParser(boundaryPath);
        boundaries = boundaryParser.parse();
                
        // read sample genes elements
        String genePath = TabFileParserTest.class.getResource("/sample_genes_chr22.tab").getPath();
        TabFileParser geneParser = new TabFileParser(genePath);
        genes = geneParser.parseGene();

        // parse toy example data set
        exampleData = new ExampleData();
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of boundaryOverlap method, of class AnnotateCNVs.
     */
    @Test
    public void testBoundaryOverlap() {
        System.out.println("boundaryOverlap");
        
        //System.out.println(inversions);
        //System.out.println(boundaries);
        
        AnnotateCNVs.boundaryOverlap(cnvs, boundaries);
        
        // check if all boundaries are read.
        // The sample file containes 30 elements
        assertEquals(boundaries.size(), 30);
        
        AnnotateCNVs.boundaryOverlap(cnvs, boundaries);
        
        // count number of inversions that overlap any boundary
        Integer cnt = 0;
        for (CNV c : cnvs.values()){
            if (! c.getBoundaryOverlap().isEmpty()){
                cnt++;
            }
        }
        // 34 of the sample CNVs should overlap completely at leaste one boundary, 
        // as calulated with the python script 'filter_overlap.py' from the barrier project.
        Integer expCnt = 34;
        assertEquals(expCnt, cnt);


    }

    /**
     * Test of annotateOverlappedGenes method, of class AnnotateCNVs.
     */
    @Test
    public void testGeneOverlap() {
        System.out.println("geneOverlap");
        AnnotateCNVs.annotateOverlappedGenes(cnvs, genes);
        
    }

    /**
     * Test of defineOverlappedDomainRegions method, of class AnnotateCNVs.
     */
    @Test
    public void testDefineOverlappedDomainRegions() {
        System.out.println("defineOverlappedDomainRegions");
        GenomicSet<CNV> cnvs = exampleData.getCnvs();
        GenomicSet<GenomicElement> domains = exampleData.getDomains();
        GenomicSet<GenomicElement> boundaries = exampleData.getBoundaries();
        
        AnnotateCNVs.boundaryOverlap(cnvs, boundaries);        
        AnnotateCNVs.defineOverlappedDomainRegions(cnvs, domains);
        
        GenomicElement leftOverlappedCnv1 = new GenomicElement("chr1", 9, 12, "leftOverlapped");
        GenomicElement rightOverlappedCnv1 = new GenomicElement("chr1", 15, 19, "rightOverlapped");
        
        assertEquals(leftOverlappedCnv1, cnvs.get("cnv1").getLeftOverlappedDomainRegion());
        assertEquals(rightOverlappedCnv1, cnvs.get("cnv1").getRightOverlappedDomainRegion());
        
        // cnv4 left border overlaps domain
        GenomicElement leftOverlappedCnv4 = new GenomicElement("chr1", 14, 14, "leftOverlapped");
        assertEquals(leftOverlappedCnv4, cnvs.get("cnv4").getLeftOverlappedDomainRegion());
    }

    /**
     * Test of annotateOverlappedGenes method, of class AnnotateCNVs.
     */
    @Test
    public void testAnnotateOverlappedGenes() {
        System.out.println("annotateOverlappedGenes");
        GenomicSet<CNV> cnvs = exampleData.getCnvs();
        GenomicSet<Gene> genes = exampleData.getGenes();

        AnnotateCNVs.annotateOverlappedGenes(cnvs, genes);

        // ovelaped genes of cnv1 is only geneB and geneD (becasue any overlap is relevant here)
        GenomicSet<Gene> cnv1genes = new GenomicSet<Gene>();
        cnv1genes.put("geneB", genes.get("geneB"));
        cnv1genes.put("geneD", genes.get("geneD"));
        assertEquals(cnv1genes, cnvs.get("cnv1").getGenesInOverlap());
        
        // cnv4 does not ovelrap any gene
        assertEquals(new GenomicSet<Gene>(), cnvs.get("cnv4").getGenesInOverlap());
    }


    /**
     * Test of overlapPhenogramScore method, of class AnnotateCNVs.
     */
    @Test
    public void testOverlapPhenogramScore() {
        System.out.println("overlapPhenogramScore");
        GenomicSet<CNV> cnvs = exampleData.getCnvs();
        GenomicSet<Gene> genes = exampleData.getGenes();
        PhenotypeData phenotypeData = exampleData.getPhenotypeData();
        AnnotateCNVs.overlapPhenogramScore(cnvs, phenotypeData);

        // cnv1 overlaps only genes B and D with phenomatch scores each of 0.29
        assertEquals(0.29, cnvs.get("cnv1").getOverlapPhenogramScore(), 0.01);
    
        // cnv4 does not overlap any gene
        assertEquals(0.0, cnvs.get("cnv4").getOverlapPhenogramScore(), 0.01);

    }

    /**
     * Test of annotateGenesInOverlapTADs method, of class AnnotateCNVs.
     */
    @Test
    public void testAnnotateGenesInOverlapTADs() {
        System.out.println("annotateGenesInOverlapTADs");
        GenomicSet<CNV> cnvs = exampleData.getCnvs();
        GenomicSet<GenomicElement> domains = exampleData.getDomains();
        GenomicSet<Gene> genes = exampleData.getGenes();

        AnnotateCNVs.annotateGenesInOverlapTADs(cnvs, domains, genes);
        
        // cnv4 should have two gens, D and A, in overlapping TADs
        assertEquals(2, cnvs.get("cnv4").getGenesInOverlapTADs().size());
        assert(cnvs.get("cnv4").getGenesInOverlapTADs().containsKey("geneA"));
        assert(cnvs.get("cnv4").getGenesInOverlapTADs().containsKey("geneD"));
        
        // cnv3 should have two gens, D and A, in overlapping TADs
        assertEquals(2, cnvs.get("cnv3").getGenesInOverlapTADs().size());
        assert(cnvs.get("cnv3").getGenesInOverlapTADs().containsKey("geneA"));
        assert(cnvs.get("cnv3").getGenesInOverlapTADs().containsKey("geneD"));
 
        assertFalse( cnvs.get("cnv3").getGenesInOverlapTADs().containsKey("geneB"));

        // cnv1 should have all four gens, A,B,C,D in overlapping TADs
        assertEquals(4, cnvs.get("cnv1").getGenesInOverlapTADs().size());
        assert(cnvs.get("cnv1").getGenesInOverlapTADs().containsKey("geneA"));
        assert(cnvs.get("cnv1").getGenesInOverlapTADs().containsKey("geneB"));
        assert(cnvs.get("cnv1").getGenesInOverlapTADs().containsKey("geneC"));
        
    }

}
