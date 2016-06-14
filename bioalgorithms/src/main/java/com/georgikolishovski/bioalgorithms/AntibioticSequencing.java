package com.georgikolishovski.bioalgorithms;

import java.util.ArrayList;
import java.util.List;

import com.georgikolishovski.bioalgorithms.biostructures.includes.IntegerMassTable;

import edu.emory.mathcs.backport.java.util.Collections;

public class AntibioticSequencing {
	
	public AntibioticSequencing() { };
	
	/**
	 * 
	 * @param peptide - an amino-acid string representing a cyclic peptide
	 * @return the theoretical spectrum of 'peptide'
	 * @author Georgi Kolishovski
	 */
	public String cyclospectrum(String peptide) {
		List<Integer> masses = new ArrayList<Integer>(); 
		List<Integer> spectrum = new ArrayList<Integer>();
		
		int peptideTotalMass = 0; // sum of the masses of all elements
		
		for(int i = 0; i < peptide.length(); ++i) {
			int m = IntegerMassTable.getIntegerMass(peptide.substring(i, i+1));
			masses.add(m);
			
			peptideTotalMass += m;
		}
		
		spectrum.add(0); spectrum.add(peptideTotalMass);
		
		masses.addAll(masses);
		
		for(int i = 0; i < peptide.length(); ++i) {
			int m = masses.get(i);
			spectrum.add(m);
			
			for(int j = 1; j < peptide.length()-1; ++j) {
				m = m + masses.get(i + j);
				spectrum.add(m);
			}
		}
		
		Collections.sort(spectrum);
		
		StringBuilder sbspectrum = new StringBuilder();
		
		for(int i : spectrum) {
			sbspectrum.append(Integer.toString(i) + " ");
		}
		
		return sbspectrum.toString().trim();
	}
}
