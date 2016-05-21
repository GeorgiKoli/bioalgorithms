package com.georgikolishovski.bioalgorithms;

import java.util.Arrays;
import java.util.Collections;

import com.georgikolishovski.bioalgorithms.salib.*;

public class DnaReplication {
	protected StringHelper sh = null;
	
	public DnaReplication() {
		sh = new StringHelper();
	}
	
	/**
	 * 
	 * @param text
	 * @param k
	 * @param d
	 * @return
	 */
	public String FrequentWordsWithMismatchesRC(String text, int k, int d) {
		StringBuilder freqPatterns = new StringBuilder();
		String strArr[] = sh.computingFrequencies(text, k).split(" ");
		Integer freqArr[] = new Integer[strArr.length];
		
		for(int i = 0; i < freqArr.length; i++) {
			freqArr[i] = 0;
		}
		
		for(int i = 0; i < freqArr.length; i++) {
			int nsum = 0;
			String pattern = sh.numberToPattern(i, k);
			String neighbors[] = sh.neighbors(pattern, d).split(" ");
			String neighborsRC[] = sh.neighbors(sh.reverseComplement(pattern), d).split(" ");
			
			for(int j = 0; j < neighbors.length; j++) {
				nsum = nsum + Integer.parseInt(strArr[(int) sh.patternToNumber(neighbors[j])]);
				nsum = nsum + Integer.parseInt(strArr[(int) sh.patternToNumber(neighborsRC[j])]);
			}
			freqArr[i] = nsum;
		}
		int maxCount = (int) Collections.max(Arrays.asList(freqArr));
		
		for(int i = 0; i < freqArr.length; i++) {
			if(freqArr[i] >= maxCount) {
				freqPatterns.append(sh.numberToPattern(i, k) + " ");
			}
		}
		return freqPatterns.toString().trim();
	}
}