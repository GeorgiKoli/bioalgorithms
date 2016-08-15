package com.georgikolishovski.bioalgorithms.util;

import java.util.ArrayList;
import java.util.List;

import com.georgikolishovski.bioalgorithms.biostructures.includes.IntegerMassTable;

import edu.emory.mathcs.backport.java.util.Collections;

public class CyclicSpectrum implements Spectrum {
	
	/**
	 * 
	 */
	public List<Integer> getSpectrum(String peptide) {
		List<Integer> spectrum = new ArrayList<Integer>();
		// list of the amino acid masses of the peptide (ex. NQEL = 114,128,129,113)
		List<Integer> intmasses = new ArrayList<Integer>();
		// total sum of the masses of all amino acids in peptide
		int peptideTotalMass = 0;
		
		for(int i = 0; i < peptide.length(); ++i) {
			int m = IntegerMassTable.getIntegerMass(peptide.substring(i, i+1)); 
			intmasses.add(m);
			peptideTotalMass += m;
		}
		
		// adding also the first (always 0) and the last (sum of all individual amino acid) 
		// masses to the spectrum
		spectrum.add(0); spectrum.add(peptideTotalMass);
		
		intmasses.addAll(intmasses); //for example: 'NQEL' + 'NQEL' = 'NQELNQEL'
		
		for(int i = 0; i < peptide.length(); i++) {
			int m = intmasses.get(i);
			spectrum.add(m);
			
			for(int j = 1; j < peptide.length()-1; j++) {
				m = m + intmasses.get(i + j);	
				spectrum.add(m);
			}
		}
		
		Collections.sort(spectrum);
		
		return spectrum;
	}
	
	/**
	 * 
	 * @param peptide - string representing an amino acid sequence
	 * @param experimental - list of numbers representing the generated experimental spectrum
	 * @return the score of the generated experimental spectrum versus the theoretical spectrum
	 * mass spectrometers generate spectra that is not ideal - they are characterized by having 
	 * both false masses and missing masses.
	 * @author Georgi Iliev
	 */
	public int getPeptideScore(String peptide, ArrayList<Integer> experimental) {
		int score = 0; int t = 0; int e = 0;
		
		List<Integer> theoretical = getSpectrum(peptide); // peptide's theoretical spectrum
		
		boolean end = (theoretical.size() > 0) ? false : true;
		
		// moving 2 indices along the theoretical & experimental spectrums/lists respectively
		while(end == false) {
			if(theoretical.get(t).intValue() == experimental.get(e).intValue()) {
				// matching masses - increase score
				score += 1; t++; e++;
			} else if(theoretical.get(t).intValue() > experimental.get(e).intValue()) {
				e++; // false mass in experimental
			} else if(theoretical.get(t).intValue() < experimental.get(e).intValue()) {
				t++; // missing mass in experimental
			}
			
			if(theoretical.get(t-1).intValue() == theoretical.get(theoretical.size()-1)
					|| experimental.get(e-1).intValue() == experimental.get(experimental.size()-1)) {
				end = true;
			}
		}
		
		return score;
	}
}