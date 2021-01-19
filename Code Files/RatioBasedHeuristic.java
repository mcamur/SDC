import java.io.IOException;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Set;
import java.util.Map;

import static java.util.stream.Collectors.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;


public class RatioBasedHeuristic {
	public static final String graphType = "..."; 
	
	public static void main(String[] args) throws IOException {

		int firstNElements = 0; //specify how many star that you would like to use to warm start
		int numNodes =0; //number of nodes in the network
		int parameter =0; // parameter used to generate each network
		double p =0; //probability value (only used in small world networks)

		String DATADIR =null;	// Specify the file location
				
		long start = System.currentTimeMillis();
		Graph<Integer, DefaultEdge> myGraph =  new SimpleGraph<>(DefaultEdge.class);
		HashMap<Integer, HashSet<Integer>> adjacencyList = new HashMap<Integer, HashSet<Integer>>();
		
		Path filePath = Paths.get(DATADIR).resolve(Paths.get("adjacencyList.txt"));         
		try {
            Files.lines(filePath).forEach(line -> {
                String[] columnValues = line.split("\\s+");
                adjacencyList.computeIfAbsent(Integer.valueOf(columnValues[0]),
                		k -> new HashSet<>()).add(Integer.valueOf(columnValues[1]));
            });
        } catch (IOException e) {
            e.printStackTrace();
        }

		
		for(int i = 0 ; i <numNodes; i++) //add nodes
			myGraph.addVertex(i);
				
		
		for(int i = 0 ; i< numNodes; i++) {	//add edges
			if(adjacencyList.get(i) != null) {
				for(Integer j : adjacencyList.get(i)) {
					myGraph.addEdge(i, j);
				}
			}
		}
		
		
		HashMap<Integer, Set<Integer>> secDegreeNodes = new HashMap<Integer, Set<Integer>>(); 
		// store the second degree nodes for the center node
		for(int i=0;i<numNodes; i++) {
			Set<Integer> myList = new HashSet<>();
			if(adjacencyList.get(i) != null) { 
				for(Integer j : adjacencyList.get(i)) {
						for(Integer k : adjacencyList.get(j)) {
							if(!adjacencyList.get(i).contains(k) && i !=k)
							  myList.add(k);
						}
				}
			}
			secDegreeNodes.put(i, myList);
		}
		
		HashMap<Integer, Set<Integer>> stars = new HashMap<Integer, Set<Integer>>();  // store all the stars	
		HashMap<Integer,Integer> objectives = new HashMap<Integer, Integer>();  	  // store objective of each star

		for(Integer center : myGraph.vertexSet()) {				
			Integer leaf = -1;	
			Set<Integer> objective = new HashSet<Integer>(); //store each neighbor node
			Set<Integer> leafNodes = new HashSet<Integer>(); 
			Set<Integer> openNeighborhood = Graphs.neighborSetOf(myGraph, center); 
			
			//Graph<Integer, DefaultEdge> subGraph =new AsSubgraph<Integer, DefaultEdge>(myGraph, openNeighborhood); //create the induced sub-graph
			
			HashMap<Integer, Integer> f2 = calculateF2(secDegreeNodes.get(center), openNeighborhood , adjacencyList);
			Set<Integer> candidate1 = new HashSet<Integer>();
			Set<Integer> candidate2 = new HashSet<Integer>();
			
			for (Map.Entry<Integer, Integer> entry : f2.entrySet()) {
			   if(entry.getValue() == 0) {
				   candidate1.add(entry.getKey());
				}else {
					candidate2.add(entry.getKey());
				}
			}
			
			HashMap<Integer, Integer> weights = new HashMap<Integer, Integer>();
			
			
			while(!candidate1.isEmpty() ||  !candidate2.isEmpty()) {
				
				if(!candidate1.isEmpty()) {
					
					weights.clear();
					weights = calculateWeights(secDegreeNodes.get(center), candidate1 , adjacencyList);
					leaf = Collections.max(weights.entrySet(), Map.Entry.comparingByValue()).getKey(); //get the key producing the max value
					
					if(weights.get(leaf) > 1) {
						leafNodes.add(leaf);
						candidate1.remove(leaf);
						secDegreeNodes.put(center,  updateSecDegreeNodes(secDegreeNodes.get(center), leaf , adjacencyList));
					}else {
						candidate1.clear();
					}						
				}else {			
					
					f2.clear();
					weights.clear();
					
					
					weights = calculateWeights(secDegreeNodes.get(center), candidate2 , adjacencyList);
					f2 = calculateF2(secDegreeNodes.get(center), candidate2 , adjacencyList);
					float val =0;						
					boolean flag= false;
					for(Integer e : candidate2) {
						if(weights.get(e) == 0) {
							continue;
						}else {
							if(weights.get(e)/ f2.get(e) > val) {
								val= weights.get(e)/ f2.get(e);
								leaf = e;
								flag =true;
							}
						}
					}
					
					if(!flag) {
						break;
					}
					leafNodes.add(leaf);
					candidate2.remove(leaf);
					secDegreeNodes.put(center,  updateSecDegreeNodes(secDegreeNodes.get(center), leaf , adjacencyList));
					
				
				   Iterator<Integer> iterator = candidate2.iterator();
				    while(iterator.hasNext()) {
				        Integer nodesToRemove = iterator.next();
						if(adjacencyList.get(leaf).contains(nodesToRemove))
							iterator.remove();	
				    }	
				    
				}
				if(candidate2.isEmpty()) {
					continue;
				}else {
					f2.clear();
					f2 = calculateF2(secDegreeNodes.get(center), candidate2 , adjacencyList);
					for (Map.Entry<Integer, Integer> entry : f2.entrySet()) {
						   if(entry.getValue() == 0) {
							   candidate1.add(entry.getKey());
							   candidate2.remove(entry.getKey());
							}
						}
				}

			}
			stars.put(center, leafNodes);
			objective.addAll(Graphs.neighborSetOf(myGraph, center)); // first, put all the nodes around the center
			if(!leafNodes.isEmpty()) {
				for(Integer ele: leafNodes) {
					objective.addAll(Graphs.neighborSetOf(myGraph, ele)); //then add all the nodes reachable from the leaves
				}
				objective.removeAll(leafNodes); // finally, remove leaves
				objectives.put(center, objective.size()-1); //remove the center
			}else {
				objectives.put(center, objective.size()); 
			}
		}
		
		
		 Map<Integer, Integer> sorted = objectives.entrySet() .stream() .sorted(Collections.reverseOrder(Map.Entry.comparingByValue())) .collect(toMap(Map.Entry::getKey, Map.Entry::getValue, 
				 (e1, e2) -> e2, LinkedHashMap::new)); 
		
		long timeTaken =System.currentTimeMillis() - start;
		FileWriter writer = new FileWriter(DATADIR+"/stars.txt", true);
		
		
		int temp =1;
		
		for (Map.Entry<Integer, Integer> entry : sorted.entrySet()) {		
			if(temp <= firstNElements) {
				temp ++;
				writer.write(entry.getKey()  + " "); //add center
				if(stars.get(entry.getKey()).isEmpty()) {
					writer.write( "\n");
				}else {
					for(Integer val : stars.get(entry.getKey())) {
						writer.write(val + " ");
					}
					writer.write( "\n");
				}
				
			}
    	} 
		

		writer.close();
		System.out.println("Time taken -> "+ (double) timeTaken /  1000.0 );
		FileWriter timeWriter = new FileWriter(DATADIR+"/starTime.txt", true);
		timeWriter.write((double) timeTaken /  1000.0 + "\n");
		timeWriter.close();
	}
	
	private static HashMap<Integer, Integer> calculateF2(Set<Integer> secDegreeNodes, Set<Integer> nei,
		HashMap<Integer, HashSet<Integer>> adjacencyList) {
		HashMap<Integer, Integer> f2 = new HashMap<Integer, Integer>();
		for(Integer j : nei) { 
			Set<Integer> count =  new HashSet<Integer>();
			for(Integer adjacent : nei) {
				if(j != adjacent) {
					if(adjacencyList.get(j).contains(adjacent)) {
						for(Integer k : secDegreeNodes) {
							if(adjacencyList.get(adjacent).contains(k)) {
								count.add(k); //duplicates are not allowed. You're safe.
							}
						}
					}
				}
			}
			f2.put(j, count.size());
		}		
		return f2;
	}

	private static Set<Integer> updateSecDegreeNodes(Set<Integer> previousSecDegreeNodes, Integer leaf, HashMap<Integer, HashSet<Integer>> adjacencyList) {
		Set<Integer> newSecDegreeNodes = new HashSet<Integer>();
		for(Integer k : previousSecDegreeNodes) {
			if(!adjacencyList.get(leaf).contains(k)) {
				newSecDegreeNodes.add(k);
			}
		}
		return newSecDegreeNodes;
	}
	private static HashMap<Integer, Integer> calculateWeights(Set<Integer> secDegreeNodes, Set<Integer> nei, HashMap<Integer, HashSet<Integer>> adjacencyList) {
		
		HashMap<Integer, Integer> weights = new HashMap<Integer, Integer>();
		for(Integer j : nei) { 
			int weight =0;
			for(Integer k : secDegreeNodes) {
				if(adjacencyList.get(j).contains(k)) {
					weight++;
				}
			}
			weights.put(j, weight);
		}
		return weights;
	}
}
