package com.georgikolishovski.bioalgorithms.biostructures.includes;

import java.util.Hashtable;

public class CodonTable {
	
	private static Hashtable<String, String> codonTable = new Hashtable<String, String>();
	
	private static void setupCodonTable() {
		codonTable.put("AUU","I"); codonTable.put("AUC","I"); codonTable.put("AUA","I");
		
		codonTable.put("CUU","L"); codonTable.put("CUC","L"); codonTable.put("CUA","L"); 
		codonTable.put("CUG","L"); codonTable.put("UUA","L"); codonTable.put("UUG","L");
		
		codonTable.put("GUU","V"); codonTable.put("GUC","V");
		codonTable.put("GUA","V"); codonTable.put("GUG","V");
		
		codonTable.put("UUU","F"); codonTable.put("UUC","F");
		
		codonTable.put("AUG","M");
		
		codonTable.put("UGU","C"); codonTable.put("UGC","C");
		
		codonTable.put("GCU","A"); codonTable.put("GCC","A");
		codonTable.put("GCA","A"); codonTable.put("GCG","A");
		
		codonTable.put("GGU","G"); codonTable.put("GGC","G");
		codonTable.put("GGA","G"); codonTable.put("GGG","G");
		
		codonTable.put("CCU","P"); codonTable.put("CCC","P");
		codonTable.put("CCA","P"); codonTable.put("CCG","P");
		
		codonTable.put("ACU","T"); codonTable.put("ACC","T");
		codonTable.put("ACA","T"); codonTable.put("ACG","T");
		
		codonTable.put("UCU","S"); codonTable.put("UCC","S"); codonTable.put("UCA","S");
		codonTable.put("UCG","S"); codonTable.put("AGU","S"); codonTable.put("AGC","S");
		
		codonTable.put("UAU","Y"); codonTable.put("UAC","Y");
		
		codonTable.put("UGG","W");
		
		codonTable.put("CAA","Q"); codonTable.put("CAG","Q");
		
		codonTable.put("AAU","N"); codonTable.put("AAC","N");
		
		codonTable.put("CAU","H"); codonTable.put("CAC","H");
		
		codonTable.put("GAA","E"); codonTable.put("GAG","E");
		
		codonTable.put("GAU","D"); codonTable.put("GAC","D");
		
		codonTable.put("AAA","K"); codonTable.put("AAG","K");
		
		codonTable.put("CGU","R"); codonTable.put("CGC","R"); codonTable.put("CGA","R");
		codonTable.put("CGG","R"); codonTable.put("AGA","R"); codonTable.put("AGG","R");
		
		codonTable.put("UAA","");
		codonTable.put("UAG","");
		codonTable.put("UGA","");
	}
	
	/**
	 * 
	 * @param codon
	 * @return the amino acid encoded by the RNA codon(3 nucleobases)
	 */
	public static String convertRnaToAminoAcid(String codon) {
		if(codonTable.isEmpty()) {
			setupCodonTable();
		}
		return codonTable.get(codon);
	}
	
	/**
	 * 
	 * @param codon
	 * @return checks if the RNA codon is a STOP codon for translation 
	 */
	public static boolean isStopCodon(String codon) {
		if(codon.equals("UAA") || codon.equals("UAG") || codon.equals("UGA")) {
			return true;
		}
		
		return false;
	}
}
