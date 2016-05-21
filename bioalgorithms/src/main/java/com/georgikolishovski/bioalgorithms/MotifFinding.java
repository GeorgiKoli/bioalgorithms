package com.georgikolishovski.bioalgorithms;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import com.georgikolishovski.bioalgorithms.salib.*;
import com.georgikolishovski.bioalgorithms.includes.*;

public class MotifFinding {
	protected StringHelper sh = null;
	
	protected StringApproxMatching sam = null;
	
	public MotifFinding() {
		sh = new StringHelper();
		
		sam = new StringApproxMatching();
	}
	
	/**
	 * 
	 * @param dna
	 * @param k
	 * @param t
	 * @return a collection of k-mers BestMotifs resulting from running GreedyMotifSearch(Dna, k, t)
	 * @author Georgi Kolishovski
	 */
	public String[] GreedyMotifSearch(String[] dna, int k, int t) {
		List<String> bestMotifs = new ArrayList<String>();
		
		for(int i = 0; i < dna.length; i++) {
			bestMotifs.add(dna[i].substring(0, k));
		}
		
		for(int i = 0; i < dna[0].length() - k + 1; i++) {
			ArrayList<String> currentMotifs = new ArrayList<String>();
			currentMotifs.add(dna[0].substring(i, i+k));
			
			for(int j = 1; j < t; j++) {
				double[][] profile = generateProfile(currentMotifs.toArray(new String[currentMotifs.size()]), true);
				currentMotifs.add(profileMostProbableKmer(dna[j], k, profile));
			}
			
			if(score(currentMotifs.toArray(new String[0])) < score(bestMotifs.toArray(new String[0]))) {
				bestMotifs.removeAll(bestMotifs);
				
				for(String s : currentMotifs) {
					bestMotifs.add(s);
				}
			}
		}
		return bestMotifs.toArray(new String[0]);
	}
	
	public String Consensus(String[] motifs) {
		return null;
	}
	
	/**
	 * 
	 * @param dna
	 * @param k
	 * @param t
	 * @return a collection BestMotifs resulting from running the function n times
	 * @author Georgi Kolishovski
	 */
	public String[] randomizedMotifSearch(String[] dna, int k, int t) {
		List<String> bestMotifs = new ArrayList<String>();
		List<String> currentMotifs = new ArrayList<String>();
		int len = dna[0].length() - k + 1;
		
		Random randNum = new Random();
		
		// collection of randomly chosen k-mers to begin with
		for(int i = 0; i < t; i++) {
			int start = randNum.nextInt(len);
			String motif = dna[i].substring(start, start + k);
			bestMotifs.add(motif);
			currentMotifs.add(motif);
		}
		
		int bestScore = score(bestMotifs.toArray(new String[0]));
		
		while(true) {
			double[][] profile = generateProfile(currentMotifs.toArray(new String[0]), true);
			currentMotifs.clear();
			
			for(int j = 0; j < t; j++) {
				currentMotifs.add(this.profileMostProbableKmer(dna[j], k, profile));
			}
			
			int currentScore = score(currentMotifs.toArray(new String[0]));
			
			if(currentScore < bestScore) {
				bestMotifs.clear();
				bestScore = currentScore;
				
				for(String s : currentMotifs) {
					bestMotifs.add(s);
				}
			} 
			else {
				break;
			}
		}
		
		return bestMotifs.toArray(new String[0]);
	}
	
	/**
	 * 
	 * @param dna
	 * @param k
	 * @param t
	 * @param N
	 * @return 
	 */
	public String[] GibbsSampler(String[] dna, int k, int t, int N) {
		List<String> bestMotifs = new ArrayList<String>();
		List<String> currentMotifs = new ArrayList<String>();
		List<String> reducedMotifs = new ArrayList<String>();
		int len = dna[0].length() - k + 1;
		
		Random randNum = new Random();
		// set random start positions
		for(int i = 0; i < t; i++) {
			int start = randNum.nextInt(len);
			String motif = dna[i].substring(start, start + k);
			bestMotifs.add(motif);
			currentMotifs.add(motif);
		}
		
		int bestScore = k * t;
		
		for(int i = 0; i < N; i++) {
			int r = randNum.nextInt(t); // r = the strand that will be replaced
			// generate a profile from the remaining strands and pick the weighted-random k-mer from
			// the strand that was left out
			for(int j = 0; j < t; j++) {
				if(j == r) {
					continue;
				}
				reducedMotifs.add(currentMotifs.get(j));
			}
			double[][] profile = generateProfile(reducedMotifs.toArray(new String[0]), true);
			
			reducedMotifs.clear();
			
			String motif = profileRandomKmer(dna[r], k, profile);
			// important to insert the newly picked motif back into its old position
			currentMotifs.set(r, motif);
			int currentScore = score(currentMotifs.toArray(new String[0]));
			
			if(currentScore < bestScore) {
				bestMotifs.clear();
				bestScore = currentScore;
				
				for(String s : currentMotifs) {
					bestMotifs.add(s);
				}
			}
		}
		
		return bestMotifs.toArray(new String[0]);
	}
	
	/**
	 * 
	 * @param k
	 * @param d
	 * @param dna
	 * @return all (k, d)-motifs in 'dna'. A k-mer is a (k, d)-motif if it appears in every string from
	 * 'dna' with at most d mismatches
	 */
	public String MotifEnumeration(int k, int d, String[] dna) {
		StringBuilder patterns = new StringBuilder();
		// create an empty set for each DNA strand and fill it with all k-mers with at most d mismatches
		List<Set<String>> strandKmers = new ArrayList<Set<String>>();
		for(int i = 0; i < dna.length; i++) {
			strandKmers.add(new HashSet<String>());
		}
		
		for(int i = 0; i < dna.length; i++) {
			for(int j = 0; j < dna[i].length() - k + 1; j++) {
				String kmer = dna[i].substring(j, j+k);
				strandKmers.get(i).addAll(Arrays.asList(sh.neighbors(kmer, d).split(" ")));
			}
		}
		
		// intersect all sets to end up with only the sequences that occur in all strands
		Set<String> result = strandKmers.get(0);
		
		for(int i = 0; i < strandKmers.size(); i++) {
			result.retainAll(strandKmers.get(i));
		}
		
		Iterator<String> iter = result.iterator();
		while(iter.hasNext()) {
			patterns.append(iter.next() + " ");
		}
		
		return patterns.toString().trim();
	}
	
	/**
	 * 
	 * @param k
	 * @param dna - a collection of strings
	 * @return a k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers  
	 * @author Georgi Kolishovski
	 */
	public String MedianString(int k, String[] dna) {
		String medstr = null;
		int minDist = k * dna.length;
		int numKmers = (int) Math.pow(4, 6);
		
		for(int i = 0; i < numKmers; i++) {
			String pattern = sh.numberToPattern(i, k);
			int dist = DistancePatternStrings(pattern, dna);
			
			if(dist < minDist) {
					minDist = dist;
					medstr = pattern;
			}
		}
		return medstr;
	}
	
	/**
	 * 
	 * @param motifs - list of DNA strings
	 * @param pseudocounts - boolean to determine if pseudocounts should be used
	 * @return profile matrix for which P(i,j) is the frequency of the i-th nucleotide in the j-th 
	 * column of the motif matrix
	 * @author Georgi Kolishovski
	 */
	private double[][] generateProfile(String[] motifs, boolean pseudocounts) {
		double[][] matrix = new double[4][motifs[0].length()];
		int count = 0;
		Map<String, Integer> trans = new HashMap<String, Integer>();
		trans.put("A", 0); 
		trans.put("C", 1); 
		trans.put("G", 2); 
		trans.put("T", 3);
		
		if(pseudocounts == true) {
			count = motifs.length + 4;
			for(int i = 0; i < 4; i++) {
				for(int j = 0; j < motifs[0].length(); j++) {
					matrix[i][j] = 1.0 / count;
				}
			}
		}
		else {
			count = motifs.length;
			for(int i = 0; i < 4; i++) {
				for(int j = 0; j < motifs[0].length(); j++) {
					matrix[i][j] = 0.0;
				}
			}
		}
		
		double baseval = 1.0 / count;
		for(int i = 0; i < motifs.length; i++) {
			for(int j = 0; j < motifs[0].length(); j++) {
				matrix[trans.get(motifs[i].substring(j, j+1))][j] += baseval;
			}
		}
		
		return matrix;
	}
	
	/**
	 * 
	 * @param text
	 * @param k
	 * @param profile - 4 x k matrix
	 * @return a Profile-most probable k-mer in Text ( if multiple answers exist, any one can be returned)
	 * @author Georgi Kolishovski
	 */
	public String profileMostProbableKmer(String text, int k, double[][] profile) {
		double maxProb = -1.0;
		String mostProbKmer = null;
		
		for(int i = 0; i < text.length() - k + 1; i++) {
			double prob = 1.0;
			for(int j = 0; j < k; j++) {
				switch(text.charAt(i+j)) {
				case 'A':
					prob = prob * profile[0][j];
					break;
				case 'C':
					prob = prob * profile[1][j];
					break;
				case 'G':
					prob = prob * profile[2][j];
					break;
				case 'T':
					prob = prob * profile[3][j];
					break;
				default:
					break;
				}
				
				if(prob < maxProb) {
					break;
				}
			}
			
			if(prob > maxProb) {
				maxProb = prob;
				mostProbKmer = text.substring(i, i+k);
			}
		}
		
		return mostProbKmer;
	}
	
	/**
	 * 
	 * @param text
	 * @param k
	 * @param profile
	 * @return
	 */
	public String profileRandomKmer(String text, int k, double[][] profile) {
		int len = text.length() - k + 1;
		double[] probabilities = new double[len]; 
		
		for(int i = 0; i < len; i++) {
			double prob = 1.0;
			for(int j = 0; j < k; j++) {
				switch(text.charAt(i+j)) {
				case 'A':
					prob = prob * profile[0][j];
					break;
				case 'C':
					prob = prob * profile[1][j];
					break;
				case 'G':
					prob = prob * profile[2][j];
					break;
				case 'T':
					prob = prob * profile[3][j];
					break;
				default:
					break;
				}
			}
			probabilities[i] = prob;
		}
		
		int randomN = RandomNumGenerator.weightedRandom(probabilities);
		
		return text.substring(randomN, randomN+k);
	}
	
	/**
	 * 
	 * @param pattern
	 * @param dna
	 * @return the sum of distances between 'pattern' and each string in 'dna' = {Dna1, ..., Dnat}
	 */
	public int DistancePatternStrings(String pattern, String[] dna) 
	{	
		int k = pattern.length();
		int totdist = 0;
		
		for(int i = 0; i < dna.length; i++) {
			int stringMinHammDist = k;
			for(int j = 0; j < dna[i].length() - k + 1; j++) {
				int hd = sam.hammingDistance(pattern, dna[i].substring(j, j+k));
				
				if(hd == 0) {
					stringMinHammDist = 0; break;
				} else if(stringMinHammDist > hd) {
					stringMinHammDist = hd;
				}
			}
			totdist = totdist + stringMinHammDist;
		}
		return totdist;
	}
	
	/**
	 * 
	 * @param motifs - collection of DNA strings
	 * @return counts the total number of unpopular (lower case) symbols in the motifs matrix
	 */
	public int score(String[] motifs) {
		Map<String, Integer> counts = new HashMap<String, Integer>();
		int score = 0; 
		
		counts.put("A", 0); counts.put("C", 1); counts.put("G", 2); counts.put("T", 3);
		
		for(int i = 0; i < motifs[0].length(); i++) {
			int[] l = new int[4];
			for(int j = 0; j < motifs.length; j++) {
				l[counts.get(motifs[j].substring(i,i+1))] += 1;
			}
			int max = 0;
			int sum = 0;
			for(int a : l) {
				sum += a;
				if(a > max) {
					max = a;
				}
			}
			score += sum - max;
		}
		
		return score;
	}
}