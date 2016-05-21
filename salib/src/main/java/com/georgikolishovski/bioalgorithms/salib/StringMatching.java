package com.georgikolishovski.bioalgorithms.salib;

public class StringMatching extends StringHelper {
	
	public StringMatching() { };
	
	/**
	 * NAIVE IMPLEMENTATION
	 * @param pattern
	 * @param text
	 * @return all starting positions in 'text' where 'pattern' appears as a substring
	 * @author Georgi Kolishovski
	 */
	public String stringMatching(String pattern, String text) {
		StringBuilder indices = new StringBuilder();
		int numwindows = text.length() - pattern.length() + 1;
		
		for(int i = 0; i < numwindows; i++) {
			if(text.substring(i).equals(pattern)) {
				indices.append(String.valueOf(i) + " ");
			}
		}
		return indices.toString().trim();
	}
}
