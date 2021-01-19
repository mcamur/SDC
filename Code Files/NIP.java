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

public class NIP {
	public static final double MSTOSEC = 0.001;
	public static final boolean warmStart = false;
	public static final boolean validIneqaulity = true;
	public static final boolean upperBoundImplement = true;
	public static final String graphType = "..."; //Specify the graph type for DATADIR
	
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
		HashMap<Integer, Integer> degree = new HashMap<>(); // degree of each node
		HashMap<Integer, Integer> bounds = new HashMap<>();	 // RHS of the constraint, which is the upper bound for the max independent set
		HashMap<Integer, Integer> upperBounds = new HashMap<>(); //upper bounds explained in Section 5.3

	   
		Scanner scanner = new Scanner(new FileReader(DATADIR+"adjacencyList.txt"));
		//create the adjacency list
        while (scanner.hasNextLine()) {
            final String[] line = scanner.nextLine().split("\\s+"); 
            adjacencyList.computeIfAbsent(Integer.valueOf(line[0]),
            		k -> new HashSet<>()).add(Integer.valueOf(line[1]));
        }	    	
        scanner.close();

        
    	Path myFilePath = Paths.get(DATADIR).resolve(Paths.get("upperBound.txt"));    //upper bound for the open neighborhood   	
        try {
            Files.lines(myFilePath).forEach(line -> {
                String[] columnValues = line.split("\\s+");
                upperBounds.put(Integer.parseInt(columnValues[0]), Integer.parseInt(columnValues[1]));
            });
        } catch (IOException e) {
            e.printStackTrace();
        }

			
    	
        if(validIneqaulity) {       	
        	myFilePath = Paths.get(DATADIR).resolve(Paths.get("bound.txt"));       	
            try {
                Files.lines(myFilePath).forEach(line -> {
                    String[] columnValues = line.split("\\s+");
                    bounds.put(Integer.parseInt(columnValues[0]), Integer.parseInt(columnValues[1]));
                });
            } catch (IOException e) {
                e.printStackTrace();
            }	           	     	
        }else { 
            for (int i = 0; i < numNodes; i++) {	
            	degree.put(i, adjacencyList.get(i).size());         	          
            }          
        }	   
        
		try {
			IloCplex cplex = new IloCplex();
		
			double heurTime=0; 
			if(upperBoundImplement) {
		        scanner = new Scanner(new FileReader(DATADIR+"heurTime.txt"));				               
		        while (scanner.hasNextLine()) {
		            final String[] line = scanner.nextLine().split("\\s+"); 
		            heurTime = Double.valueOf(line[0]);
		        }
		        scanner.close();
			}
			cplex.setParam(IloCplex.Param.TimeLimit, 2*1800-heurTime);
			cplex.setOut(null);
			
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
	        	
	            cplex.addLe(cplex.sum(x[i],y[i],z[i]), 1, "first_"+ i); // a node can be either center, leaf, neighbor or has nothing to do with the star
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
        
        for (int i = 0; i < numNodes; i++) { //constraint for neighbors meaning that you need to be adjacent to the star
        	expr.clear();	
        		if(adjacencyList.get(i) != null) {
	        		adjacencyList.get(i).forEach((val) -> { 
	        			try { 
	        				expr.addTerm(y[val], 1);  
	        				expr.addTerm(x[val], 1);
	        			}catch (IloException e) {
	    					e.printStackTrace();
	    				} 
	        		});
	        		
        		}	
        		cplex.addLe(z[i],expr, "second_"+i);	
        }
       
        for (int i = 0; i < numNodes; i++) { //constraint to be a leaf node you need to be adjacent to center
        	expr.clear();	
        	if(adjacencyList.get(i) != null) {
        		adjacencyList.get(i).forEach((val) -> { 
        			try { 
        				expr.addTerm(x[val], 1);
        			}catch (IloException e) {
    					e.printStackTrace();
    				} 
        		});
        	}
        	cplex.addLe(y[i],expr, "third_"+i);		                	
        }
        
        if(validIneqaulity) {
            for (int i = 0; i < numNodes; i++) { //constraint to make sure that no leaf pair share an edge
            	expr.clear();	
            	if(adjacencyList.get(i) != null) {
            		adjacencyList.get(i).forEach((val) -> { 
            			try { 
            				expr.addTerm(y[val], 1);
            			}catch (IloException e) {
        					e.printStackTrace();
        				} 
            		});
            	}
            	cplex.addLe(expr, cplex.prod(bounds.get(i) , cplex.diff(1, y[i])), "fourth_"+i);		                	
            }
        }else { 
            for (int i = 0; i < numNodes; i++) { //constraint to make sure that no leaf pair share an edge
            	expr.clear();	
            	if(adjacencyList.get(i) != null) {
            		adjacencyList.get(i).forEach((val) -> { 
            			try { 
            				expr.addTerm(y[val], 1);
            			}catch (IloException e) {
        					e.printStackTrace();
        				} 
            		});
            	}
            	cplex.addLe(expr, cplex.prod(degree.get(i) , cplex.diff(1, y[i])), "fourth_"+i);		                	
            }
        }
        expr.clear();  
        
        for(int i=0; i< numNodes; i++)  //only one node can be center
    		expr.addTerm(x[i], 1);        
      
         cplex.addEq(expr, 1, "Fifth");	         
        
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
	        	  IloNumVar[] vars = new IloNumVar[1 + e.getValue().size()];
	        	  double[] vals = new double[vars.length];
	        	  int idx = 0;
	        	  vars[idx] = x[e.getKey()];
	        	  vals[idx++] = 1;
	        	  for (Integer yIdx : e.getValue()) {
	        	    vars[idx] = y[yIdx];
	        	    vals[idx++] = 1;
	        	  }
	        	  cplex.addMIPStart(vars, vals);
	        	}
	          
         }

         
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
