package com.georgikolishovski.bioalgorithms;

import java.text.Collator;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import com.georgikolishovski.bioalgorithms.salib.*;

public class GenomeAssembling {
	
	protected StringHelper sh = null;
	
	public GenomeAssembling() {
		sh = new StringHelper();
	}
	
	/**
	 * 
	 * @param k
	 * @param text
	 * @return constructs the de Bruijn graph of a string in the form of an adjacency list
	 */
	public Map<String, List<String>> deBruijnGraph(int k, String text) {
		int len = text.length() - k + 2; // number of (k-1)-mers
		Map<String, List<String>> graph = new HashMap<String, List<String>>(); 
		
		for(int i = 0; i < len - 1; ++i) {
			String key = text.substring(i, i + k - 1);
			if(!graph.containsKey(key)) {
				graph.put(key, new ArrayList<String>());
			}
			graph.get(key).add(text.substring(i+1, i+k));
		}
		
		return graph;
	}
	
	/**
	 * 
	 * @param patterns
	 * @return
	 */
	public Map<String, List<String>> deBruijnGraphList(String[] patterns) {
		Map<String, List<String>> graph = new HashMap<String, List<String>>();
		
		for(int i = 0; i < patterns.length; ++i) {
			String prefix = sh.prefix(patterns[i]);
			if(!graph.containsKey(prefix)) {
				graph.put(prefix, new ArrayList<String>());
			}
			graph.get(prefix).add(sh.suffix(patterns[i]));
		}
		
		return graph;
	}
	
	/**
	 * 
	 * @param adjlist
	 * @return 
	 */
	public List<Integer> eulerianCycle(String[] adjlist) {
		int n = adjlist.length;
		@SuppressWarnings("unchecked")
		List<Integer>[] g = (List<Integer>[]) new List[n];
		for(int i = 0; i < n; ++i) {
			g[i] = new ArrayList<Integer>();
		}
		
		for(int i = 0; i < n; ++i) {
			String[] st = adjlist[i].split("->");
			String[] v = st[1].split(",");
			
			for(int j = 0; j < v.length; ++j) {
				g[Integer.parseInt(st[0].trim())].add(Integer.parseInt(v[j].trim()));
			}
		}
		
		int[] curEdge = new int[n];
		List<Integer> res = new ArrayList<Integer>();
		Stack<Integer> stack = new Stack<Integer>();
		int v = 0;
		stack.add(v);
		
		while (!stack.isEmpty()) {
			 v = stack.pop();
			 while (curEdge[v] < g[v].size()) {
				 stack.push(v);
				 v = g[v].get(curEdge[v]++);
			 }
			 res.add(v);
		}
		
		Collections.reverse(res);
		
		return res;
	}
	
	/**
	 * 
	 * @param text
	 * @param k
	 * @return the collection of all k-mer substrings of 'text' (including repeated k-mers).
	 * The k-mers are listed in lexicographic order
	 */
	public String[] stringComposition(String text, int k) {
		int kmerNumber = text.length() - k+1;
		List<String> composition = new ArrayList<String>();
		
		for(int i = 0; i < kmerNumber; ++i) {
			composition.add(text.substring(i, i + k));
		}
		
		return lexicographicOrder(composition).toArray(new String[0]);
	}
	
	/**
	 * 
	 * @param patterns
	 */
	public List<String[]> overlapGraph(String[] patterns) {
		Map<String, List<String>> targetnodes = new HashMap<String, List<String>>();
		List<String[]> graph = new ArrayList<String[]>();
		
		// group target nodes by prefix
		for(int i = 0; i < patterns.length; ++i) {
			String prefix = sh.prefix(patterns[i]);
			
			if(!targetnodes.containsKey(prefix)) {
				targetnodes.put(prefix, new ArrayList<String>());
			}
			targetnodes.get(prefix).add(patterns[i]);
		}
		// map target nodes to suffix
		for(int i = 0; i < patterns.length; i++) {
			String suffix = sh.suffix(patterns[i]);
			
			if(targetnodes.containsKey(suffix)) {
				for(String p : targetnodes.get(suffix)) {
					String[] a = new String[2];
					a[0] = patterns[i];
					a[1] = p;
					graph.add(a);
				}
			}
		}
		return graph;
	}
	
	/**
	 * 
	 * @param str
	 * @return
	 */
	public List<String> lexicographicOrder(List<String> str) {
		Collections.sort(str, Collator.getInstance());
		
		return str;
	}
}