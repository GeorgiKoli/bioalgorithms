package com.georgikolishovski.bioalgorithms.biostructures.includes;

import java.util.Hashtable;

public class IntegerMassTable {
	
	private static Hashtable<String, Integer> intMassTable = new Hashtable<String, Integer>();
	
	// 'static constructor'
	static {
		intMassTable.put("G", 57);
		intMassTable.put("A", 71);
		intMassTable.put("S", 87);
		intMassTable.put("P", 97);
		intMassTable.put("V", 99);
		intMassTable.put("T", 101);
		intMassTable.put("C", 103);
		intMassTable.put("I", 113);
		intMassTable.put("L", 113);
		intMassTable.put("N", 114);
		intMassTable.put("D", 115);
		intMassTable.put("K", 128);
		intMassTable.put("Q", 128);
		intMassTable.put("E", 129);
		intMassTable.put("M", 131);
		intMassTable.put("H", 137);
		intMassTable.put("F", 147);
		intMassTable.put("R", 156);
		intMassTable.put("Y", 163);
		intMassTable.put("W", 186);
	};
	
	/**
	 * GETTER METHOD
	 * @param aa - an amino acid symbol
	 * @return the integer mass of 'aa'
	 */
	public static int getIntegerMass(String aa) {
		return intMassTable.get(aa);
	}
}
