
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/*This code produces an upper bound on the objective function for every node.
* You should only update the path and parameters. 
*/

public class UpperBound {
	public static final String graphType = "..."; 

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub		
		int numNodes =0; //number of nodes in the network
		int parameter =0; // parameter used to generate each network
		double p =0; //probability value (only used in small world networks)

		String DATADIR =null;	// Specify the file location
		
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
		
		HashMap<Integer, Set<Integer>> newSec = new HashMap<Integer, Set<Integer>>();
		for(int i=0;i<numNodes; i++) {
			Set<Integer> myList = new HashSet<>();
			if(adjacencyList.get(i) != null) { //if there is no adjacent node, then UB must be zero.
				for(Integer j : adjacencyList.get(i)) {
					if(adjacencyList.get(j) != null) { 
						for(Integer k : adjacencyList.get(j)) {
							if(!adjacencyList.get(i).contains(k) && i !=k)
							  myList.add(k);
						}
					}
				}
			}
			newSec.put(i, myList);
		}
		

		
		HashMap<Integer, Integer> UB = new HashMap<Integer, Integer>();
		
		
		long start = System.currentTimeMillis();
		//Heuristic starts here
		for(int i =0; i< numNodes; i++) {
			Set<Integer> secondNodes = new HashSet<>();
			Set<Integer> adjacentNodes = new HashSet<Integer>();
			
			if(adjacencyList.get(i) != null) { //if there is no adjacent node, then UB must be zero.
				adjacentNodes = adjacencyList.get(i);					
			}else {
				UB.put(i, 0);
				continue;
			}
			
			if(newSec.get(i) != null) { // if there is no second degree node, then UB is the neighborhood.
				secondNodes = newSec.get(i);				
			}
			else {
				UB.put(i, adjacentNodes.size());
				continue;
			}
			
			
			HashMap<Integer, Integer> pred = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> found = new HashMap<Integer, Integer>();
			
			for(Integer sec : secondNodes) {
				pred.put(sec, -1);
				found.put(sec, 0);
			}			
			
			HashMap<Integer, Integer> remaining = new HashMap<Integer, Integer>();
			int rem=0;
				for(Integer adjacent: adjacentNodes) {
					for( Integer sec : secondNodes) {
						if(adjacencyList.get(adjacent).contains(sec)) {
							rem +=1;
						}	    				
	    			}
					remaining.put(adjacent,rem);
					rem=0;
	    		}		

				for(Integer adjacent: adjacentNodes) {
					for( Integer sec : secondNodes) {
						if(adjacencyList.get(adjacent).contains(sec)) {
							 if(found.get(sec) < 0.5) {
								 pred.put(sec, adjacent);
								 found.put(sec,1);
							 }else if( 0.5 < found.get(sec) &&  found.get(sec) < 1.5) {
								 remaining.put(adjacent,remaining.get(adjacent)-1);
								 remaining.put(pred.get(sec), remaining.get(pred.get(sec)) -1);
								 found.put(sec,2);
							 }else {
								 found.put(sec,found.get(sec)+1);
								 remaining.put(adjacent,remaining.get(adjacent)-1);
							 }
						}	    				
	    			}
					
	    		}		
			
			  int counter=0;
			  for (Integer k : remaining.values()) {			        	
	             if(k > 0.5) {
	            	 counter +=1;
	             }
		      }
			  if(counter > 0.5) {
				  UB.put(i, adjacentNodes.size() -counter + secondNodes.size() ); 
			  }else {
				  UB.put(i, adjacentNodes.size() + secondNodes.size() - 1); 
			  }
			  
		}
		long timeTaken =System.currentTimeMillis() - start;
		System.out.println();
		
		//Create a file called upper bound. This will be used in the Benders implementation.
		FileWriter writer = new FileWriter(DATADIR+"upperBound.txt", true);
		UB.entrySet().forEach(entry -> {
			try {
				writer.write(entry.getKey() + " " + entry.getValue() + "\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		});
		writer.close();

		FileWriter timeWriter = new FileWriter(DATADIR+"heurTime.txt", true);
		timeWriter.write((double) timeTaken /  1000.0 + "\n");
		timeWriter.close();
			
 }		       
}


