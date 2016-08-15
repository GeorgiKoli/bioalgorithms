package com.georgikolishovski.bioalgorithms.util;

import java.util.ArrayList;
import java.util.List;

public interface Spectrum {
	
	public List<Integer> getSpectrum(String peptide);
	
	public int getPeptideScore(String peptide, ArrayList<Integer> experimental);
}