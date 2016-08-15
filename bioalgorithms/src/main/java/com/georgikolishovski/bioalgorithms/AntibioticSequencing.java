package com.georgikolishovski.bioalgorithms;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.georgikolishovski.bioalgorithms.biostructures.includes.IntegerMassTable;
import com.georgikolishovski.bioalgorithms.includes.*;
import com.georgikolishovski.bioalgorithms.util.*;

public class AntibioticSequencing {
	Spectrum cs = new CyclicSpectrum();
	Spectrum ls = new LinearSpectrum();
	
	public AntibioticSequencing() { };
	
	/**
	 * 
	 * @param spectrum - an array of (possibly repeated) integers corresponding 
	 * to an ideal (no mistakes) experimental spectrum
	 * @return the integer mass strings of a cyclic peptide whose theoretical spectrum 
	 * matches the given experimental spectrum   
	 * @author Georgi Iliev
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
				ArrayList<Integer> cycloSpectrum = (ArrayList<Integer>) cs.getSpectrum(peptide.toString());
				
				if(peptideMass == parentMass) {
					if(spectrum.containsAll(cycloSpectrum)) {
						results.add(peptideMasses);
					}
					it.remove();
				}
				else if(!spectrum.containsAll(ls.getSpectrum(peptide.toString()))) {
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
	 * @author Georgi Iliev
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
	
	/**
	 * 
	 * @param lb
	 * @param N
	 * @param spectrum
	 * @return
	 */
	public ArrayList<String> leaderboardCut(ArrayList<String> lb, int N, ArrayList<Integer> spectrum) {
		// if the size is less than ...
		if(lb.size() <= N) {
			return lb;
		}
		
		List<String> g = new ArrayList<String>();
		
		Map<String, Integer> mmm = new HashMap<String, Integer>();
		
		for(int i = 0; i < lb.size(); i++) {
			mmm.put(lb.get(i), cs.getPeptideScore(lb.get(i), spectrum));
		}
		
		Map<String, Integer> res = new HashMap<String, Integer>();
		
		res = LeaderBoard.sortByValue(mmm, N);
		
		Iterator<?> it = res.entrySet().iterator();
		
		while(it.hasNext()) {
			Map.Entry pair = (Map.Entry)it.next();
			g.add((String) pair.getKey());
			//System.out.println(pair.getKey() + " " + pair.getValue());
		}
		
		return (ArrayList<String>) g;
	}
	
	/**
	 * 
	 * @param spectrum
	 * @return
	 */
	public String leaderboardCyclopeptideSequencing(ArrayList<Integer> spectrum, int N) {
		int parentMass = spectrum.get(spectrum.size() - 1);
		// 
		Set<ArrayList<Integer>> leaderBoard = new HashSet<ArrayList<Integer>>();
		leaderBoard.add(new ArrayList<Integer>());
		
		String leaderPeptide = "";
		
		while(leaderBoard.size() > 0) {
			
			leaderBoard = expand(leaderBoard);
			
			Iterator<ArrayList<Integer>> it = leaderBoard.iterator();
			
			while(it.hasNext()) {
				StringBuilder peptide = new StringBuilder();
				ArrayList<Integer> peptideMasses = it.next();
				int peptideMass = 0;
				
				for(int i = 0; i < peptideMasses.size(); ++i) {
					// String k = IntegerMassTable.getKeyByValue(peptideMasses.get(i));
					
					// if(k != null) {
						// peptide.append(k);
					// }
					peptide.append(peptideMasses.get(i) + " ");
					
					peptideMass += peptideMasses.get(i);
				}
				
				// ArrayList<Integer> cycloSpectrum = (ArrayList<Integer>) cyclospectrum(peptide.toString());
				
				if(peptideMass == parentMass) {
					if(cs.getPeptideScore(IntegerMassTable.getPeptideFromMasses(peptide.toString().trim()), spectrum) > cs.getPeptideScore(leaderPeptide, spectrum)) {
						leaderPeptide = IntegerMassTable.getPeptideFromMasses(peptide.toString().trim());
					}
				}
				else if(peptideMass > parentMass) {
					it.remove();
				}
			}
			
			System.out.println("BANG: " + leaderBoard.size());
			// leaderBoard = leaderboardCut(leaderBoard, N, spectrum);
			
			
		}
		
		return IntegerMassTable.strMassFromPeptide(leaderPeptide);
	}
}
