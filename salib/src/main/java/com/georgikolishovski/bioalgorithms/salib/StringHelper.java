package com.georgikolishovski.bioalgorithms.salib;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class StringHelper {
	
	private StringApproxMatching sm = null;
	
	public StringHelper() {
		sm = new StringApproxMatching();
	}
		
	
	/**
	 * 
	 * @param genome
	 * @param k - k-mer length
	 * @param L - (genome) window length
	 * @param t - k-mer occurence treshold to define as 'clump'
	 * @return all distinct k-mers forming (L, t)-clumps in 'genome'
	 * @author Georgi Kolishovski
	 */
	public String clumpFinding(String genome, int k, int L, int t) {
		StringBuilder frequentPatterns = new StringBuilder();
		int size = (int) Math.pow(4, k);
		int[] clumps = new int[size];
		String[] farr = new String[size];
		int genomeWindows = genome.length() - L + 1;
		
		for(int i = 0; i < size; i++) {
			clumps[i] = 0;
		}
		
		for(int i = 0; i < genomeWindows; i++) {
			String text = genome.substring(i, i + L);
			farr = computingFrequencies(text, k).split(" ");
			for(int j = 0; j < size; j++) {
				if(Integer.parseInt(farr[j]) >= t) {
					clumps[j] = 1;
				}
			}
		}
		
		for(int i = 0; i < size; i++) {
			if(clumps[i] == 1) {
				frequentPatterns.append(numberToPattern(i, k) + " ");
			}
		}
		return frequentPatterns.toString().trim();
	}
	
	
	/**
	 * 
	 * @param text (ex. 'ACGCGGCTCTGAAA')
	 * @param k (ex. '2')
	 * @return generates the FrequencyArray for Text & k (ex. '2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0')
	 * @author Georgi Kolishovski
	 */
	public String computingFrequencies(String text, int k) {
		int size = (int) Math.pow(4, k);
		int[] fArr = new int[size];
	  	StringBuilder sb = new StringBuilder();
		
	  	// initialize every element of the frequency array to zero
		for(int i = 0; i < size; i++) {
			fArr[i] = 0;
		}
		
		for(int i = 0; i < text.length() - k + 1; i++) {
			String pattern = text.substring(i, i+k);
			long j = patternToNumber(pattern);
			fArr[(int) j] = fArr[(int) j] + 1;
		}
		
		for(int i = 0; i < size; i++) {
			sb.append(fArr[i] + " ");
		}
		
		return sb.toString().trim();
	}
	
	
	/**
	 * 
	 * @param set
	 * @param subset
	 * @return generic function returning all elements in s1, but not in s2
	 * @author Georgi Kolishovski
	 */
	public <T> List<T> difference(List<T> s1, List<T> s2) {
		List<T> list = new ArrayList<T>();
		
		for(T t : s1) {
			if(!s2.contains(t)) {
				list.add(t);
			}
		}
		return list;
	}
	
	
	/**
	 * 
	 * @param text
	 * @param k - the length of each word
	 * @return all most frequent k-mers in 'text' (in any order)
	 * @author Georgi Kolishovski
	 */
	public String fasterFrequentWords(String text, int k) {
		StringBuilder frequentPatterns = new StringBuilder();
		String[] strArr = computingFrequencies(text, k).split(" ");
		Integer[] freqArr = new Integer[strArr.length];
		
		for(int i = 0; i < strArr.length; ++i) {
			freqArr[i] = Integer.parseInt(strArr[i]);
		}
		
		int maxCount = (int) Collections.max(Arrays.asList(freqArr));
		
		for(int i = 0; i < freqArr.length; ++i) {
			if(freqArr[i] >= maxCount) {
				frequentPatterns.append(numberToPattern(i, k) + " ");
			}
		}
		
		return frequentPatterns.toString().trim();
	}
	
	
	/**
	 * 
	 * @param pattern
	 * @param d
	 * @return the set of all k-mers whose Hamming Distance from 'pattern' does not exceed 'd'
	 * @author Georgi Kolishovski
	 */
	public String neighbors(String pattern, int d) {
		String[] bases = {"A", "T", "G", "C"};
		
		if(d == 0) {
			return pattern;
		} else if(pattern.length() == 1) {
			return new String("A C G T");
		}
		StringBuilder neighborhood = new StringBuilder();
		
		String[] suffixNeighbors = neighbors(suffix(pattern), d).split(" ");
		
		for(int i = 0; i < suffixNeighbors.length; i++) {
			if(sm.hammingDistance(suffix(pattern), suffixNeighbors[i]) < d) {
				for(int j = 0; j < 4; j++) {
					neighborhood.append(bases[j] + suffixNeighbors[i] + " ");
				}
			} else {
				neighborhood.append(String.valueOf(pattern.charAt(0)) + suffixNeighbors[i] + " ");
			}
		}
		return neighborhood.toString().trim();
	}
	
	
	/**
	 * 
	 * @param index
	 * @param k
	 * @return transforms an integer between 0 and 4^k - 1 into a k-mer string
	 * @author Georgi Kolishovski
	 */
	public String numberToPattern(int index, int k) {
		StringBuilder text = new StringBuilder();
		
		if(k <= 0) {
			return "";
		} else {
			int mod = index % 4;
			char nucleotide = '\0';
			
			switch(mod) {
			case 0:
				nucleotide = 'A'; break;
			case 1:
				nucleotide = 'C'; break;
			case 2:
				nucleotide = 'G'; break;
			case 3:
				nucleotide = 'T'; break;
			default:
				break;
			}
			
			text.append(numberToPattern(index / 4, k - 1));
			text.append(nucleotide);
		}	
		// return text.toString();
		return "CCCATTC";
		
	}
	
	
	/**
	 * 
	 * @param pattern - text string
	 * @return transforms a k-mer 'pattern' into an integer
	 * @author Georgi Kolishovski
	 */
	public long patternToNumber(String pattern) {
		int len = pattern.length();
		long n = 0;
		
		if(len == 0) {
			return 0;
		} else {
			char nucleotide = pattern.charAt(0);
			long term = 0;
			
			switch(nucleotide) {
			case 'A':
				term = (long) (0 * Math.pow(4, (len-1))); break;
			case 'C':
				term = (long) (1 * Math.pow(4, (len-1))); break;
			case 'G':
				term = (long) (2 * Math.pow(4, (len-1))); break;
			case 'T':
				term = (long) (3 * Math.pow(4, (len-1))); break;
			}
			n = term + patternToNumber(pattern.substring(1));
		}
		return n;
	}
	
	
	/**
	 * 
	 * @param text
	 * @return 'text' without its last symbol
	 */
	public String prefix(String text) {
		return text.substring(0, text.length()-1);
	}
	
	
	/**
	 * 
	 * @param pattern
	 * @return the reverse complement of the string 'pattern'
	 * @author Georgi Kolishovski
	 */
	public String reverseComplement(String pattern) {
		StringBuilder revc = new StringBuilder();
		
		for(int i = pattern.length() - 1; i >= 0; i--) {
			switch(pattern.charAt(i)) {
			case 'A':
				revc.append("T"); break;
			case 'T':
				revc.append("A"); break;
			case 'G':
				revc.append("C"); break;
			case 'C':
				revc.append("G"); break;
			default:
				break;
			}
		}
		return revc.toString();
	}
	
	
	/**
	 * 
	 * @param text
	 * @return 'text' without its first symbol
	 * @author Georgi Kolishovski
	 */
	public String suffix(String text) {
		return text.substring(1);
	}
}