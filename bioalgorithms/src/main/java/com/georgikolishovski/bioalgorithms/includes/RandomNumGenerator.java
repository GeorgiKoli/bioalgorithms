package com.georgikolishovski.bioalgorithms.includes;

public class RandomNumGenerator {
	public RandomNumGenerator() { }
	
	/**
	 * 
	 * @param distributions
	 * @return generates a random number based on the probabilities distribution in 'probabilities' 
	 */
	static public int weightedRandom(double[] distributions) {
		double totalWeight = 0.0;
		
		for(int i = 0; i < distributions.length; i++) {
			totalWeight += distributions[i];
		}
		int randomIndex = -1;
		double randomNum = Math.random() * totalWeight;
		
		for(int i =0; i < distributions.length; i++) {
			randomNum -= distributions[i];
			if(randomNum <= 0.0) {
				randomIndex = i;
				break;
			}
		};
		return randomIndex;
	}
}
