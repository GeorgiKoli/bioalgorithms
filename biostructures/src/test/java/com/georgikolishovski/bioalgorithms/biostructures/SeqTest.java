package com.georgikolishovski.bioalgorithms.biostructures;

import static org.junit.Assert.*;
import org.junit.Test;

import com.georgikolishovski.bioalgorithms.biostructures.Seq;

import edu.emory.mathcs.backport.java.util.Arrays;

public class SeqTest {
	
	@Test
	public void testPeptideEncoding() {
		Seq s = new Seq("ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA", "DNA");
		String[] expectedOutput = new String[] {"ATGGCC", "ATGGCC", "GGCCAT"};
		String[] observedOutput = s.peptideEncoding("MA");
		Arrays.sort(observedOutput);
		assertArrayEquals(expectedOutput, observedOutput);
	}
}
