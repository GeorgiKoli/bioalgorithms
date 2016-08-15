package com.georgikolishovski.bioalgorithms.util;

import java.util.ArrayList;
import java.util.List;

import com.georgikolishovski.bioalgorithms.biostructures.includes.IntegerMassTable;

import edu.emory.mathcs.backport.java.util.Collections;

public class LinearSpectrum implements Spectrum {
	/**
	 * 
	 * @param peptide - an amino acid string representing a linear peptide
	 * @return generates and returns the theoretical spectrum of 'peptide'
	 * @author Georgi Kolishovski
	 */
	public List<Integer> getSpectrum(String peptide) {
		List<Integer> spectrum = new ArrayList<Integer>();
		spectrum.add(0);
		List<Integer> prefixMass = new ArrayList<Integer>();
		prefixMass.add(0);
		
		List<Integer> masses = IntegerMassTable.getIntegerMassValues();
		masses.remove(masses.indexOf(113)); // 113 repeats twice in the integer mass table
		masses.remove(masses.indexOf(128)); // 128 repeats twice in the integer mass table
		
		for(int i = 0; i < peptide.length(); i++) {
			for(int j = 0; j < masses.size(); j++) {
				if(masses.get(j).equals(IntegerMassTable.getIntegerMass(peptide.substring(i, i+1)))) {
					prefixMass.add(prefixMass.get(prefixMass.size()-1) + masses.get(j));
				}
			}
		}
		
		for(int i = 0; i < peptide.length(); i++) {
			for(int j = i + 1; j < peptide.length()+1; j++) {
				spectrum.add(prefixMass.get(j) - prefixMass.get(i));
			}
		}
		
		Collections.sort(spectrum);
		return spectrum;
	}

	public int getPeptideScore(String peptide, ArrayList<Integer> experimental) {
		int score = 0; int t = 0; int e = 0;
		
		List<Integer> theoretical = getSpectrum(peptide); // peptide's theoretical spectrum
		System.out.println("Theoretical: " + theoretical);
		
		System.out.println("Experimental: " + experimental);
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