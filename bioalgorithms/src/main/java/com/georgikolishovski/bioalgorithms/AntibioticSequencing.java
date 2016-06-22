package com.georgikolishovski.bioalgorithms;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import com.georgikolishovski.bioalgorithms.biostructures.includes.IntegerMassTable;

import edu.emory.mathcs.backport.java.util.Collections;

public class AntibioticSequencing {
	public AntibioticSequencing() { };
	
	/**
	 * 
	 * @param peptide - an amino-acid string representing a cyclic peptide
	 * @return generates and returns the theoretical spectrum of 'peptide'
	 * @author Georgi Kolishovski
	 */
	public List<Integer> cyclospectrum(String peptide) {
		List<Integer> spectrum = new ArrayList<Integer>();
		// list of the amino acid masses in peptide (ex. NQEL = 114,128,129,113)
		List<Integer> masses = new ArrayList<Integer>();
		int peptideTotalMass = 0; // total sum of the masses of all amino acids in peptide
		
		for(int i = 0; i < peptide.length(); ++i) {
			int m = IntegerMassTable.getIntegerMass(peptide.substring(i, i+1));
			masses.add(m);
			peptideTotalMass += m;
		}
		
		// adding first (always 0) and last (sum of all individual amino acid) masses to spectrum
		spectrum.add(0); spectrum.add(peptideTotalMass);
		
		masses.addAll(masses); // for example: 'NQEL' + 'NQEL' = 'NQELNQEL'
		
		for(int i = 0; i < peptide.length(); ++i) {
			int m = masses.get(i);
			spectrum.add(m);
			
			for(int j = 1; j < peptide.length()-1; ++j) {
				m = m + masses.get(i + j);	
				spectrum.add(m);
			}
		}
		
		Collections.sort(spectrum);
		return spectrum;
	}
	
	/**
	 * 
	 * @param peptide string representing an amino acid sequence
	 * @param experimental list of numbers representing the generated experimental spectrum
	 * @return the score of the generated experimental spectrum versus the theoretical spectrum
	 * mass spectrometers generate spectra that is not ideal - they are characterized by having 
	 * both false masses and missing masses.
	 * @author Georgi Kolishovski
	 */
	public int cycloScore(String peptide, ArrayList<Integer> experimental) {
		int score = 0; int t = 0; int e = 0;
		List<Integer> theoretical = cyclospectrum(peptide); // peptide's theoretical spectrum
		
		boolean end = (theoretical.size() > 0) ? false : true;
		
		// the loop moves 2 indices along the theoretical & experimental spectrums/lists respectively
		while(end == false) {
			if(theoretical.get(t).intValue() == experimental.get(e).intValue()) {
				// matching masses - increase score
				score += 1; t++; e++;
			} else if(theoretical.get(t).intValue() > experimental.get(e).intValue()) {
				e++; // false mass in experimental
			} else if(theoretical.get(t).intValue() < experimental.get(e).intValue()) {
				t++; // missing mass in experimental
			}
			
			if(theoretical.get(t-1).intValue() == theoretical.get(theoretical.size()-1)) {
				end = true;
			}
		}	
		return score;
	}
	
	
	/**
	 * 
	 * @param peptide - an amino acid string representing a linear peptide
	 * @return generates and returns the theoretical spectrum of 'peptide'
	 * @author Georgi Kolishovski
	 */
	public List<Integer> linearspectrum(String peptide) {
		List<Integer> spectrum = new ArrayList<Integer>();
		// list of the amino acid masses in peptide (ex. NQEL = [114,128,129,113])
		List<Integer> masses = new ArrayList<Integer>();
		int peptideTotalMass = 0; // total sum of the masses of all amino acids in peptide
		
		for(int i = 0; i < peptide.length(); ++i) {
			int m = IntegerMassTable.getIntegerMass(peptide.substring(i, i+1));
			masses.add(m);
			peptideTotalMass += m;
		}
		// adding first (always 0) and last (sum of all individual amino acid) masses to spectrum
		spectrum.add(0); spectrum.add(peptideTotalMass);
		
		for(int i = 0; i < peptide.length(); ++i) {
			int m = masses.get(i);
			spectrum.add(m);
			
			for(int j = i+1; j < peptide.length(); ++j) {
				m = m + masses.get(j);
				spectrum.add(m);
			}
		}
		
		Collections.sort(spectrum);
		return spectrum;
	}
	
	/**
	 * 
	 * @param spectrum - an array of (possibly repeated) integers corresponding 
	 * to an ideal (no mistakes) experimental spectrum
	 * @return the integer mass strings of a cyclic peptide whose theoretical spectrum 
	 * matches the given experimental spectrum   
	 * @author Georgi Kolishovski
	 */
	public List<ArrayList<Integer>> cyclopeptideSequencing(ArrayList<Integer> spectrum) {
		int parentMass = spectrum.get(spectrum.size()-1);
		List<ArrayList<Integer>> results = new ArrayList<ArrayList<Integer>>();
		
		Set<ArrayList<Integer>> peptides = new HashSet<ArrayList<Integer>>();
		// add an initial element containing just an empty set
		peptides.add(new ArrayList<Integer>());
		
		while(peptides.size() > 0) {
			peptides = expand(peptides);
			
			Iterator<ArrayList<Integer>> it = peptides.iterator();
			
			while(it.hasNext()) {
				StringBuilder peptide = new StringBuilder();
				ArrayList<Integer> peptideMasses = it.next();
				int peptideMass = 0;
				
				for(int i = 0; i < peptideMasses.size(); ++i) {
					String k = IntegerMassTable.getKeyByValue(peptideMasses.get(i));
					
					if(k != null) 
						peptide.append(k);
					
					peptideMass += peptideMasses.get(i);
				}
				// get this peptide's theoretical spectrum
				ArrayList<Integer> cycloSpectrum = (ArrayList<Integer>) cyclospectrum(peptide.toString());
				
				if(peptideMass == parentMass) {
					if(spectrum.containsAll(cycloSpectrum)) {
						results.add(peptideMasses);
					}
					it.remove();
				}
				else if(!spectrum.containsAll(linearspectrum(peptide.toString()))) {
					it.remove();
				}
			}
		}
		return results;
	}
	
	/**
	 * 
	 * @param kmers - peptides list - each peptide is of length K
	 * @return new collection containing all possible extensions of the peptides in 
	 * 'kmers' by a single amino acid
	 * @author Georgi Kolishovski
	 */
	public Set<ArrayList<Integer>> expand(Set<ArrayList<Integer>> kmers) {
		List<Integer> masses = IntegerMassTable.getIntegerMassValues();
		masses.remove(masses.indexOf(113)); // 113 repeats twice in the integer mass table
		masses.remove(masses.indexOf(128)); // 128 repeats twice in the integer mass table
		
		Set<ArrayList<Integer>> expandedKmers = new HashSet<ArrayList<Integer>>();
		Iterator<ArrayList<Integer>> it = kmers.iterator();
		
		while(it.hasNext()) { // for each peptide in 'kmers'
			ArrayList<Integer> kmer = it.next();
			// add one additional amino acid (for all possible amino acids)
			for(int i = 0; i < masses.size(); ++i) {
				ArrayList<Integer> p = new ArrayList<Integer>();
				for(int k : kmer) {
					p.add(k);
				}
				p.add(masses.get(i));
				expandedKmers.add(p);
			}
		}
		return expandedKmers; // size of expandedKmers is: length(kmers) * 18
	}
}
