package com.georgikolishovski.bioalgorithms.salib;

public class StringApproxMatching extends StringMatching {
	
	/**
	 * NAIVE IMPLEMENTATION
	 * @param pattern string
	 * @param text string
	 * @param d 
	 * @return all starting positions in 'text' where 'pattern' appears as a substring with at most 'd' mismatches
	 * @author Georgi Kolishovski
	 */
	public String approximateStringMatching(String pattern, String text, int d) {
		StringBuilder indices = new StringBuilder();
		int numwindows = text.length() - pattern.length() + 1;
		
		for(int i = 0; i < numwindows; i++) {
			int hammingDist = hammingDistance(pattern, text.substring(i, i+pattern.length()));
			
			if(hammingDist <= d) {
				indices.append(String.valueOf(i) + " ");
			}
		}
		return indices.toString().trim();
	}
	
	/**
	 * 
	 * @param dna1
	 * @param dna2
	 * @return an integer value representing the Hamming distance between dna1 & dna2
	 * @author Georgi Kolishovski
	 */
	public int hammingDistance(String dna1, String dna2) {
		int distance = 0;
		
		// 'dan1' and 'dna2' must be of equal length; if not - throw exception 
		if(dna1.length() != dna2.length()) {
			throw new IllegalArgumentException(" the function arguments must be of equal length");
		}
		
		for(int i = 0; i < dna1.length(); i++) {
			if(dna1.charAt(i) != dna2.charAt(i)) {
				distance += 1;
			}
		}
		
		return distance;
	}
}
