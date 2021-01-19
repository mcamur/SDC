import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;

import ilog.cplex.*;
import ilog.concert.*;

public class VCIP {
	public static final boolean warmStart = false;
	public static final boolean upperBoundImplement = true;
	public static final String graphType = "..."; //Specify the graph type for DATADIR
	public static final double MSTOSEC = 0.001;

	public static void main(String[] args) throws IOException {
		
		FileWriter writer = new FileWriter("...", true); //Specify the location where you want to save the results
		
		if(graphType.equals("Barabasi") || graphType.equals("ErdosRenyi")  ) { //Specify what information you would like to save for each network instance
			writer.write("..." );
		}else if(graphType.equals("SmallWorld")){
			writer.write("..");
		}else {
			writer.write("...");
		}

		int numNodes =0; //number of nodes in the network
		int parameter =0; // parameter used to generate each network
		double p =0; //probability value (only used in small world networks)

		
		String DATADIR =null;	// Specify the file location
	

		
		HashMap<Integer, HashSet<Integer>> adjacencyList = new HashMap<Integer, HashSet<Integer>>();				
		Scanner scanner = new Scanner(new FileReader(DATADIR+"adjacencyList.txt"));
			//create the adjacency list
	        while (scanner.hasNextLine()) {
	            final String[] line = scanner.nextLine().split("\\s+"); 
	            adjacencyList.computeIfAbsent(Integer.valueOf(line[0]),
	            		k -> new HashSet<>()).add(Integer.valueOf(line[1]));
	        }
	        scanner.close();
	
		        
			try {
			IloCplex cplex = new IloCplex();
			double heurTime=0; 
			
	        HashMap<Integer, Integer> upperBounds = new HashMap<>();
	    	Path myFilePath = Paths.get(DATADIR).resolve(Paths.get("upperBound.txt"));    //upper bound for the open neighborhood   	
	        try {
	            Files.lines(myFilePath).forEach(line -> {
	                String[] columnValues = line.split("\\s+");
	                upperBounds.put(Integer.parseInt(columnValues[0]), Integer.parseInt(columnValues[1]));
	            });
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
	
			
			if(upperBoundImplement) {
		        scanner = new Scanner(new FileReader(DATADIR+"heurTime.txt"));				               
		        while (scanner.hasNextLine()) {
		            final String[] line = scanner.nextLine().split("\\s+"); 
		            heurTime = Double.valueOf(line[0]);
		        }
		        scanner.close();
			}
			cplex.setParam(	IloCplex.Param.TimeLimit,  2*1800-heurTime);
			
	        IloIntVar[] x = new IloIntVar[numNodes]; //center node
	        IloIntVar[] y = new IloIntVar[numNodes]; //leaf node
	        IloNumVar[] z = new IloNumVar[numNodes]; //neighbor node
	
	        
	        IloLinearNumExpr expr = cplex.linearNumExpr();
	        for(int i=0; i< numNodes; i++) { //Define the variables and generate the first constraint
	        	x[i] = cplex.boolVar("x_"+ i);
	        	y[i] = cplex.boolVar("y_"+ i);
	        	z[i] = cplex.numVar(0,1,"z_"+ i);         	 
	        	cplex.add(x[i]);	  			        	
	        	cplex.add(y[i]);	        				        	
	        	cplex.add(z[i]);
	        	
	        	cplex.addLe(cplex.sum(y[i],z[i]), 1, "first_"+ i);
	        	expr.addTerm(z[i], 1);
	        }
	        cplex.addMaximize(expr);
	        
	        if(upperBoundImplement) {
		        for (int i = 0; i < numNodes; i++) { //upper bound for the size of open neighborhood
		        	expr.addTerm(x[i], upperBounds.get(i));		                	
		        	}
		        cplex.addLe(cplex.sum(z) , expr, "upper_Bound");
		        expr.clear();
	        }
	        
	        for (int i = 0; i < numNodes; i++) { //to be neighbor, you need to be adjacent to a node in star
	        	expr.clear();	
	        		adjacencyList.get(i).forEach((val) -> { 
	        			try { 
	        				expr.addTerm(y[val], 1);  
	        			}catch (IloException e) {
	    					e.printStackTrace();
	    				} 
	        		});
	        	cplex.addLe(z[i],expr, "second_"+i);		                	
	        }
	           
	        for (int i = 0; i < numNodes; i++) { //constraint to be in star meaning that be adjacent to the center node
	        	expr.clear();	
	        		adjacencyList.get(i).forEach((val) -> { 
	        			try { 
	        				expr.addTerm(x[val], 1);
	        			}catch (IloException e) {
	    					e.printStackTrace();
	    				} 
	        		});
	        		expr.addTerm(x[i], 1);
	        	cplex.addLe(y[i],expr, "third_"+i);		                	
	        }
	        
	        for (int i = 0; i < numNodes; i++) { //no leaf is adjacent to another node
	        	int node = i;
	        		adjacencyList.get(i).forEach((val) -> { 
	        			if(val > node) {
	        			try { 
	        				expr.clear();
	        				expr.addTerm(y[node], 1);
	        				expr.addTerm(y[val], 1);
	        				expr.addTerm(x[node], -1);
	        				expr.addTerm(x[val], -1);
	        				cplex.addLe(expr,1, "fourth_"+node);	
	        			}catch (IloException e) {
	    					e.printStackTrace();
	    				} 
	        		 }	
	        		});
	        		                	
	        }
	        
	        expr.clear();
	      for(int i=0; i< numNodes; i++)  //only one node can be center
	    	expr.addTerm(x[i], 1);        
	      
	         cplex.addEq(expr, 1, "Fifth");
	       
	         expr.clear();
	         for(int i=0; i< numNodes; i++)  
	        	 cplex.addLe(x[i], y[i]);
	         
	         if(warmStart) {		        	 
	     		TreeMap<Integer, Set<Integer>> stars = new TreeMap<Integer, Set<Integer>>();  
	     		myFilePath = Paths.get(DATADIR).resolve(Paths.get("stars.txt"));    	
	            try {
	                Files.lines(myFilePath).forEach(line -> {
	                    String[] columnValues = line.split("\\s+");
	                    Set<Integer> leaves = new HashSet<Integer>();
	                    for(int i =1 ; i < columnValues.length; i++) {
	                    	leaves.add(Integer.valueOf(columnValues[i]));
	                    }
	                      stars.put(Integer.valueOf(columnValues[0]),leaves );
	                  
	                });
	            } catch (IOException e) {
	                e.printStackTrace();
	            }
	         
		         for (Map.Entry<Integer,Set<Integer>> e : stars.entrySet()) {
		        	  IloNumVar[] vars = new IloNumVar[2 + e.getValue().size()];
		        	  double[] vals = new double[vars.length];
		        	  int idx = 0;
		        	  vars[idx] = x[e.getKey()];
		        	  vals[idx++] = 1;
		        	  vars[idx] = y[e.getKey()];
		        	  vals[idx++] = 1;
		        	  for (Integer yIdx : e.getValue()) {
		        	    vars[idx] = y[yIdx];
		        	    vals[idx++] = 1;
		        	  } 
		        	  cplex.addMIPStart(vars, vals);
		        	}
		          
	         }
	         
	
	        cplex.setOut(null);
	        long start = System.currentTimeMillis();
	        if (cplex.solve()) {
	        	//save any information you would like to analyze
	     		writer.write((System.currentTimeMillis() - start) * MSTOSEC + ",");	
	            writer.write("\n");
	            writer.flush();
	        } else {
	        	writer.write("Solution status: " +   cplex.getStatus() + "\n");
	        }
		
		        cplex.end();
		        expr.clear();
	 
			}catch (IloException exc) {
		         System.err.println("Concert exception '" + exc + "' caught");
		      }   
			
			writer.close();
			System.out.println("Experiments are over!");
	}
}


