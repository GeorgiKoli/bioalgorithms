package com.georgikolishovski.bioalgorithms.biostructures;

public class DnaMolecule implements Molecule {
	static String[] bases = {"A", "T", "G", "C"};
	
	/**
	 * 
	 * @return "A", "T", "G", and "C"
	 * @author Georgi Kolishovski
	 */
	public static String[] getBases() {
		return bases;
	}
}
