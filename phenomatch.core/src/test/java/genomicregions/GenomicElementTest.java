/*
 * Copyright (C) 2014 Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package genomicregions;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import java.util.ArrayList;
import de.charite.compbio.jannovar.impl.intervals.Interval;
import de.charite.compbio.jannovar.impl.intervals.MutableInterval;

/**
 *
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class GenomicElementTest {
    
    public GenomicElementTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
        GenomicElement gi5_16 = new GenomicElement("chr1", 5, 16, "instance");
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testGenmicInterval() {
        
        GenomicElement gi = new GenomicElement("chr1", 1, 10, "aName");
        assertTrue("Start smaller or equal end?", gi.getStart() <= gi.getEnd());
                

    }

    @Test(expected=IllegalArgumentException.class)
    public void testGenmicIntervalNegativCoordinate() {
            GenomicElement wrongGI = new GenomicElement("chr1", -1, 1, "wrongGI");
    }
    
    @Test(expected=IllegalArgumentException.class)
    public void testGenmicIntervalInvertedInterval() {
            GenomicElement wrongGI = new GenomicElement("chr1", 10, 1, "wrongGI");
    }

    /**
     * Test of hasOverlap method, of class GenomicElement.
     */
    @Test
    public void testHasOverlap() {
        System.out.println("hasOverlap");
        GenomicElement instance = new GenomicElement("chr1", 5, 16, "instance");
        GenomicElement giPositive = new GenomicElement("chr1", 10, 20, "pos");
        GenomicElement giNegative = new GenomicElement("chr1", 20, 30, "neg");
        
        boolean posResult = instance.hasOverlap(giPositive);
        boolean negResult = instance.hasOverlap(giNegative);
        assertTrue("has overlap", posResult);
        assertFalse("no overlap", negResult);
    }

    /**
     * Test of getChr method, of class GenomicElement.
     */
    @Test
    public void testGetChr() {
        System.out.println("getChr");
        GenomicElement instance = new GenomicElement("chr1", 0, 1, "name");
        String expResult = "chr1";
        String result = instance.getChr();
        assertEquals(expResult, result);
    }

    /**
     * Test of getStart method, of class GenomicElement.
     */
    @Test
    public void testGetStart() {
        System.out.println("getStart");
        GenomicElement instance = new GenomicElement("chr1", 0, 1, "name");
        int expResult = 0;
        int result = instance.getStart();
        assertEquals(expResult, result);
    }

    /**
     * Test of getEnd method, of class GenomicElement.
     */
    @Test
    public void testGetEnd() {
        System.out.println("getEnd");
        GenomicElement instance = new GenomicElement("chr1", 0, 1, "name");
        int expResult = 1;
        int result = instance.getEnd();
        assertEquals(expResult, result);
    }

    /**
     * Test of getName method, of class GenomicElement.
     */
    @Test
    public void testGetName() {
        System.out.println("getName");
        GenomicElement instance = new GenomicElement("chr1", 0, 1, "name");
        String expResult = "name";
        String result = instance.getName();
        assertEquals(expResult, result);
    }

    /**
     * Test of toString method, of class GenomicElement.
     */
    @Test
    public void testToString() {
        System.out.println("toString");
        GenomicElement instance = new GenomicElement("chr1", 0, 1, "name");
        //String expResult = "name|chr1:0-0";
        String expResult = "name:chr1:[0,1)";
        String result = instance.toString();
        assertEquals(expResult, result);
    }

    /**
     * Test of equals method, of class GenomicElement.
     */
    @Test
    public void testEquals() {
        System.out.println("equals");
        GenomicElement instance = new GenomicElement("chr1", 0, 1, "name");
        GenomicElement other = new GenomicElement("chr1", 0, 1, "name");
        boolean expResult = true;
        boolean result = instance.equals(other);
        assertEquals(expResult, result);
    }
    
   

    /**
     * Test of toOutputLine method, of class GenomicElement.
     */
    @Test
    public void testToOutputLine() {
        System.out.println("toOutputLine");
        GenomicElement instance = new GenomicElement("chr1", 0, 1, "name");
        String expResult = "chr1\t0\t1\tname";
        String result = instance.toOutputLine();
        assertEquals(expResult, result);
    }

    /**
     * Test of getOutputHeaderLine method, of class GenomicElement.
     */
    @Test
    public void testGetOutputHeaderLine() {
        System.out.println("getOutputHeaderLine");
        GenomicElement instance =  new GenomicElement("chr1", 0, 1, "name");
        String expResult = "#chr\tstart\tend\tname";
        String result = instance.getOutputHeaderLine();
        assertEquals(expResult, result);
    }

    /**
     * Test of toInterval method, of class GenomicElement.
     * @throws java.lang.Exception
     */
    @Test
    public void testToInterval() throws Exception {
        System.out.println("toInterval");
        GenomicElement instance =  new GenomicElement("chr1", 0, 10, "name");
        MutableInterval expResult = new MutableInterval(0, 9, instance);
        MutableInterval result = instance.toInterval();
        assertEquals(expResult.toString(), result.toString());
    }

    /**
     * Test of toInterval method, of class GenomicElement with zero-length element.
     * @throws java.lang.Exception
     */
    @Test(expected = Exception.class)
    public void testToIntervalZeroLenght() throws Exception {
        System.out.println("toInterval with lengtz zero");
        GenomicElement zeroLength = new GenomicElement("chr1", 10, 10, "name");
        
        // toInterval() method should now throw an Excpeiont
        MutableInterval zeroIV = zeroLength.toInterval();

    }

    /**
     * Test of completeOverlapped method, of class GenomicElement.
     */
    @Test
    public void testCompleteOverlapped() {
        System.out.println("completeOverlapped");
        GenomicElement other1 = new GenomicElement("chr1", 0, 100, "name");
        GenomicElement other2 = new GenomicElement("chr1", 15, 100, "name");
        GenomicElement instance = new GenomicElement("chr1", 10, 20, "name");
        boolean result1 = instance.completeOverlapped(other1);
        boolean result2 = instance.completeOverlapped(other2);
        assertEquals(true, result1);
        assertEquals(false, result2);
    }

    /**
     * Test of compareTo method, of class GenomicElement.
     */
    @Test
    public void testCompareTo() {
        System.out.println("compareTo");
        GenomicElement instance = new GenomicElement("chr1", 0, 100, "name");
        GenomicElement same = new GenomicElement("chr1", 0, 100, "name");
        GenomicElement other1 = new GenomicElement("chr1", 0, 200, "name");
        GenomicElement other2 = new GenomicElement("chr1", 0, 100, "name2");
        GenomicElement other3 = new GenomicElement("chr2", 0, 100, "name");
        GenomicElement other4 = new GenomicElement("chr0", 0, 100, "name");
        assertEquals(0, instance.compareTo(same));
        assertEquals(-1, instance.compareTo(other1));
        assertTrue(instance.compareTo(other2) > 0);
        assertTrue(instance.compareTo(other3) < 0 );
        assertTrue(instance.compareTo(other4) > 0);
    }

    /**
     * Test of hashCode method, of class GenomicElement.
     */
    @Test
    public void testHashCode() {
        System.out.println("hashCode");
        GenomicElement instance = new GenomicElement("chr1", 0, 100, "name");
        GenomicElement same = new GenomicElement("chr1", 0, 100, "name");
        GenomicElement other1 = new GenomicElement("chr1", 0, 200, "name");
        GenomicElement other2 = new GenomicElement("chr1", 0, 100, "name2");
        int expResult = 0;
        int result = instance.hashCode();
        assertEquals(instance.hashCode(), same.hashCode());
        assertNotSame(instance.hashCode(), other1.hashCode());
        assertNotSame(instance.hashCode(), other2.hashCode());
    }

    /**
     * Test of reciprocalOverlap method, of class GenomicElement.
     */
    @Test
    public void testReciprocalOverlap() {
        System.out.println("reciprocalOverlap");
        GenomicElement instance = new GenomicElement("chr1", 0, 100, "name");
        GenomicElement other1 = new GenomicElement("chr1", 10, 110, "name1");
        GenomicElement other2 = new GenomicElement("chr1", 0, 300, "name2");
        GenomicElement other3 = new GenomicElement("chr1", 500, 600, "name3");
        GenomicElement other4 = new GenomicElement("chr1", 99, 100, "name4");
        double fraction = 0.5;

        assertTrue(instance.reciprocalOverlap(other1, fraction));
        assertFalse(instance.reciprocalOverlap(other2, fraction));
        assertFalse(instance.reciprocalOverlap(other3, fraction));
        assertTrue(instance.reciprocalOverlap(other4, 0.01));
    }

    /**
     * Test of length method, of class GenomicElement.
     */
    @Test
    public void testLength() {
        System.out.println("length");
        GenomicElement instance = new GenomicElement("chr1", 0, 100, "name");
        GenomicElement other1 = new GenomicElement("chr1", 50, 200, "name1");
        GenomicElement other2 = new GenomicElement("chr1", 100, 100, "name2");
        
        assertEquals(100, instance.length());
        assertEquals(150, other1.length());
        assertEquals(0, other2.length());
    }

}
