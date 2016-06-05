package com.georgikolishovski.bioalgorithms.salib;

import static org.junit.Assert.*;
import org.junit.Test;

import com.georgikolishovski.bioalgorithms.salib.StringHelper;

import edu.emory.mathcs.backport.java.util.Arrays;

public class StringHelperTest {
	@Test
	public void testClumpFinding() {
		assertTrue(true);
	}
	
	@Test 
	public void testComputingFrequencies() {
		StringHelper sht = new StringHelper();
		String[] expectedOutput = new String[]{"2", "1", "0", "0", "0", "0", "2", "2", "1", "2", "1", "0", "0", "1", "1", "0"};
		String[] observedOutput = sht.computingFrequencies("ACGCGGCTCTGAAA", 2).split(" ");
		assertArrayEquals(expectedOutput, observedOutput);
	}
	
	/*
	@Test
	public void testFasterFrequentWords() {
		StringHelper sht = new StringHelper();
		String[] expectedOutput = new String[]{"CATG", "GCAT"};
		String[] observedOutput = sht.fasterFrequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4).split(" ");
		Arrays.sort(observedOutput);
		
		assertArrayEquals(expectedOutput, observedOutput);
	} */
	
	@Test
	public void testNeighbors() {
		StringHelper sht = new StringHelper();
		String[] expectedOutput = new String[]{"AAG", "ACA", "ACC", "ACG", "ACT", "AGG", "ATG", "CCG", "GCG", "TCG"};
		String[] observedOutput = sht.neighbors("ACG", 1).split(" ");
		Arrays.sort(observedOutput);
		assertArrayEquals(expectedOutput, observedOutput);
	}
	
	@Test
	public void testNumberToPattern() {
		StringHelper sht = new StringHelper();
		String expectedOutput = "CCCATTC";
		assertEquals(expectedOutput, sht.numberToPattern(5437, 7));
	}
	
	@Test
	public void testPatternToNumber() {
		StringHelper sht = new StringHelper();
		int expectedOutput = 912;
		assertEquals(expectedOutput, sht.patternToNumber("ATGCAA"));
	}
	
	@Test
	public void testReverseComplement() {
		StringHelper sht = new StringHelper();
		String expectedOutput = "ACCGGGTTTT";
		assertEquals(expectedOutput, sht.reverseComplement("AAAACCCGGT"));
	}
}