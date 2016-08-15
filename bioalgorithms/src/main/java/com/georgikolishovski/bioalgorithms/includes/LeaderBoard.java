package com.georgikolishovski.bioalgorithms.includes;

import java.util.*;

public class LeaderBoard {
	
	public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> map, int N) {
		
		List<Map.Entry<K, V>> ls = new LinkedList<Map.Entry<K, V>>(map.entrySet());
		
		Collections.sort(ls, new Comparator<Map.Entry<K, V>>()
		{
			public int compare(Map.Entry<K, V> o1, Map.Entry<K, V> o2) {
				return (o2.getValue().compareTo(o1.getValue()));
			}
		});
		
		V v = ls.get(N-1).getValue(); // GK
		
		Map<K, V> result = new LinkedHashMap<K, V>();
		
		for(Map.Entry<K, V> e : ls) {
			if(e.getValue().compareTo(v) >= 0) {
				result.put(e.getKey(), e.getValue());
			}
		}
		
		// System.out.println("LB: " + v);
		
		return result;
	}
}