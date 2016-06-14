package com.georgikolishovski.bioalgorithms.biostructures;

import java.util.ArrayList;
import java.util.List;

import com.georgikolishovski.bioalgorithms.biostructures.includes.CodonTable;
import com.georgikolishovski.bioalgorithms.salib.*;

public class Seq {
	// valid alphabets
	private enum Alphabet {DNA, RNA, PROTEIN};
	int flag = 0;
	private final String alphabet;
	private final String sequence;
	
	public Seq(String s, String a) throws IllegalArgumentException { 
		this.sequence = s; // assign the sequence
		
		for(Alphabet al : Alphabet.values()) {
			if(a.equalsIgnoreCase(al.name())) {
				flag = 1; break;
			}
		}
		if(flag == 0) {
			throw new IllegalArgumentException("Invalid alphabet... accepted options are: 'DNA', 'RNA', 'PROTEIN'");
		}
		
		this.alphabet = Alphabet.valueOf(a).name(); // assign the alphabet
	}
	
	/**
	 * GETTER METHOD
	 * @return sequence
	 */
	public String getSequence() {
		return sequence;
	}
	
	/**
	 * 
	 * @return
	 * @throws UnsupportedOperationException
	 */
	public String revc() throws UnsupportedOperationException {
		if(!alphabet.equals(Alphabet.DNA.name())) {
			String msg = alphabet + " sequences cannot have reverse complement.";
			throw new UnsupportedOperationException(msg);
		}
		StringHelper sp = new StringHelper();
		
		return sp.reverseComplement(sequence);
	}
	
	/**
	 * 
	 * @param peptide - an amino acid string
	 * @return all kmers (substrings of length k) from the DNA sequence that encode 'peptide'
	 * @author Georgi Kolishovski
	 */
	public String[] peptideEncoding(String peptide) throws UnsupportedOperationException {
		int offset = peptide.length() * 3;
		List<String> kmers = new ArrayList<String>();
		
		if(!alphabet.equals(Alphabet.DNA.name())) {
			String msg = "";
			
			throw new UnsupportedOperationException(msg);
		}
		
		for(int j = 0; j < 3; ++j) {
			for(int i = j; i < (sequence.length() - offset + 1); i = i+3) {
				Seq s1 = new Seq(sequence.substring(i, i + offset), "DNA");
				// System.out.println(s1.getSequence());
				String str = s1.translate();
				if(str.equals(peptide)) {
					kmers.add(s1.getSequence());
				}
				
				Seq s2 = new Seq(s1.revc(), alphabet);
				if(s2.translate().equals(peptide)) {
					kmers.add(s1.getSequence());
				}
			}
		}
		return kmers.toArray(new String[0]);
	}
	
	/**
	 * 
	 * @return the transcribed RNA string
	 * @author Georgi Kolishovski
	 */
	public String transcribe() throws UnsupportedOperationException {
		
		if(!alphabet.equals(Alphabet.DNA.name())) {	
			String msg = alphabet + " sequences cannot be transcribed to RNA."; 
			
			throw new UnsupportedOperationException(msg);
		}
		return sequence.replaceAll("T", "U");
	}
	
	/**
	 * 
	 * @return the translation of the RNA sequence into the amino acid string 'peptide'
	 * @author Georgi Kolishovski
	 */
	public String translate() throws UnsupportedOperationException {
		String rnaseq = "";
		String protein = "";
		
		if(alphabet.equals(Alphabet.DNA.name())) {
			rnaseq = transcribe();
		} else if(alphabet.equals(Alphabet.RNA.name())) {
			rnaseq = sequence;
		} else {
			String msg = alphabet + " sequences cannot be translated to proteins.";
			throw new UnsupportedOperationException(msg);
		}
		
		for(int i = 0; i < rnaseq.length(); i = i + 3) {
			protein = protein + CodonTable.convertRnaToAminoAcid(rnaseq.substring(i, i+3));
		}
		
		return protein;
	}
}
