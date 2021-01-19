import ilog.cplex.*;
import ilog.concert.*;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

public class BendersNIP {
	public static final double MSTOSEC = 0.001;
	public static final double FUZZ = 1e-6;
	public static final boolean validIneqaulity = true; 
	public static final boolean seperateFractional = true; 
	public static final boolean warmStart = false; 
	public static final String graphType = "..."; //Specify the graph type for DATADIR
	public static final double GAMMA = 1; // Gamma 1 
	private static final class Neighbor {
	      private final int numNodes;
	      private final HashMap<Integer, HashSet<Integer>>  adjacencyList;

	      
		public Neighbor(int numNodes, HashMap<Integer, HashSet<Integer>>  adjacencyList) throws IloException {
	         this.numNodes = numNodes;
	         this.adjacencyList = adjacencyList;			 
     	          
		}

		public IloRange[] separate(IloNumVar[] x, IloNumVar[] y, IloNumVar[] estObj, double[] xSol, double[] ySol, double[] zMaster) throws IloException {
			IloRange[] cut = new IloRange[numNodes];
		    double[] alpha = new double[numNodes];
     		double[] beta =new double[numNodes];
     		double[] sum = new double[1];    
     		double[] dualSol = new double[numNodes]; //objective of dual problem
     		IloCplex cplex = new IloCplex();	
		      for(int i =0; i< numNodes; i++) {
		    	  sum[0] =0;
		    	  if(adjacencyList.get(i) != null) {
	        		adjacencyList.get(i).forEach((val) -> { 	        			
	        			 sum[0] += (xSol[val] + ySol[val]);	        			 
	        		});
		    	  }

		    	  if( (1- xSol[i] - ySol[i]) > FUZZ ) {
		        		if( (1- xSol[i] - ySol[i]) >   sum[0] +FUZZ ) {
		        			alpha[i] = 0;
		        			beta[i] = 1;
		        		}else if((1- xSol[i] - ySol[i]) <   sum[0] -FUZZ ) {
		        			alpha[i] = 1;
		        			beta[i]  =  0;
		        		}else {
		        			alpha[i] =  GAMMA;
		        			beta[i]  =  1-GAMMA;
		        		}
		    	  }else {
		        		if(sum[0] < 0.5) {        	
		        			alpha[i] = GAMMA;
		        			beta[i]  =  1-GAMMA;	        
		        		}else {
		        			alpha[i] =  1;
		        			beta[i]  =  0;
		        		}
		    	  }
		      }
		      for(int i =0; i< numNodes; i++) {
		    	  dualSol[i]=0;		    	  
		      }
		      
		      for(int i =0; i< numNodes; i++) { //dual objective for each cut
		    	  sum[0] =0;
		    	  dualSol[i]  += (1- xSol[i] - ySol[i])* alpha[i];
		    	  if(adjacencyList.get(i) != null) {
	        		adjacencyList.get(i).forEach((val) -> { 
	        			 sum[0] += xSol[val] + ySol[val];
	        		});
		    	  }
	        		dualSol[i] +=  (sum[0] * beta[i]);
		      }
		  
		      

    	   for (int i = 0; i < numNodes; i++) {
	        	 if (zMaster[i] > dualSol[i]  + FUZZ) {	
	        		IloLinearNumExpr expr = cplex.linearNumExpr();
        		 	expr.setConstant(alpha[i]);	
     	        	expr.addTerm(-alpha[i],x[i]);
     	        	expr.addTerm(-alpha[i],y[i]);
     	        	int node = i;
     	        	if(adjacencyList.get(i) != null) {
  	     	        	adjacencyList.get(i).forEach((val) -> { 
      	        			try {
      	        				expr.addTerm(beta[node], x[val] );
      							expr.addTerm(beta[node], y[val]);
      						} catch (IloException e) {		
      							e.printStackTrace();
      						}
      	        		});	 	
     	        	}
	     	       cut[i] = (IloRange) cplex.le(estObj[i],expr);
	     	       expr.clear();
	  		        	 }
	        	   }
	        cplex.end();
			return cut;
		}
	}
	
	
	 private static final class BendersPPINCallback implements IloCplex.Callback.Function {
	      private final IloIntVar[] x;
	      private final IloIntVar[] y;
	      private final IloNumVar[] estObj;
	      private final Neighbor[] neighbors;
	      private final HashMap<Integer, HashSet<Integer>> adjacencyList;
		  ;
	      public BendersPPINCallback(IloIntVar x[], IloIntVar[] y,  IloNumVar[] estObj, HashMap<Integer, HashSet<Integer>> adjacencyList, int numThreads) {
	    	  this.x = x;
	    	  this.y = y;
	    	  this.estObj = estObj;
	    	  this.neighbors = new Neighbor[numThreads];
	    	  this.adjacencyList = adjacencyList;
	      }
	      
		@Override
		public void invoke(IloCplex.Callback.Context context) throws IloException {
			try {
			int threadNo = context.getIntInfo(IloCplex.Callback.Context.Info.ThreadId);
			int numNodes = x.length;
			 // setup	
	         if (context.inThreadUp()) {
	        	 neighbors[threadNo] = new Neighbor(numNodes, adjacencyList);
	             return;
	          }
	       
	      // teardown
	         if (context.inThreadDown()) {	        	 
	        	 neighbors[threadNo] = null;	        	 
	             return;
	          }
	         
	         double[] xSol = new double[numNodes];
	         double[] ySol = new double[numNodes];
	         double[] zMaster =new double[numNodes];
	         
	         if ( context.inCandidate() ) {
	             if ( !context.isCandidatePoint() ) // The model is always bounded
	                throw new IloException("Unbounded solution");
	             for (int i = 0; i < numNodes; ++i) {	            	 
	                xSol[i] = context.getCandidatePoint(x[i]);
	                ySol[i] = context.getCandidatePoint(y[i]);	 
	                zMaster[i] = context.getCandidatePoint(estObj[i]);	                
	             }
	             
	               
	          }else if( context.inRelaxation() ) {
	                  for (int i = 0; i < numNodes; ++i) {
	                	    xSol[i] = context.getRelaxationPoint(x[i]);
	    	                ySol[i] = context.getRelaxationPoint(y[i]);	 
	    	                zMaster[i] = context.getRelaxationPoint(estObj[i]);	
	                  }
	               
	          } else {
	              throw new IloException("Unexpected contextID");
	          }
	         
	         Neighbor neighbor = neighbors[threadNo];
	         IloRange[] violated = new IloRange[numNodes];
	         violated = neighbor.separate(x, y,estObj, xSol, ySol, zMaster);
	         for(int i =0; i< numNodes ; i++) {
	         if (violated[i] != null) {
	             // Add the cut	        	
	             if ( context.inCandidate() ) {
	            	// System.out.println(violated);
	            	 context.rejectCandidate(violated[i]);
	             }else if( context.inRelaxation() ){ 
	                 context.addUserCut(violated[i],
                             IloCplex.CutManagement.UseCutPurge,
                             false);
	             }else
	                throw new IloException("Unexpected contextID");
	             }
	         }
	         
		   } catch (Exception e) {
			      System.out.println("Exception from callback: " + e);
			      throw e;
			   }  
	         
		}
		 
	 }
	
	public static void main(String[] args) throws IloException, IOException {		

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
	    scanner = new Scanner(new FileReader(DATADIR+"heurTime.txt"));
	    double heurTime=0; // There is no need for this	        
	    while (scanner.hasNextLine()) {
	        final String[] line = scanner.nextLine().split("\\s+"); 
	        heurTime = Double.valueOf(line[0]);
	    }
	    scanner.close();
		
	    IloCplex masterCplex = new IloCplex();		
		masterCplex.setParam(IloCplex.Param.TimeLimit, 2*1800-heurTime);
		masterCplex.setParam(IloCplex.Param.MIP.Strategy.VariableSelect, 3); //Strong branching
		masterCplex.setParam(IloCplex.Param.Emphasis.MIP, 2); // Optimality over Feasibility
		masterCplex.setParam(IloCplex.Param.MIP.Strategy.RINSHeur, 1000);
		masterCplex.setOut(null);
		
		try {								  		   
	        
			HashMap<Integer, Integer> degree = new HashMap<>();
			HashMap<Integer, Integer> bounds = new HashMap<>();	
			HashMap<Integer, Integer> upperBounds = new HashMap<>();	
	
	    	Path filePath = Paths.get(DATADIR).resolve(Paths.get("bound.txt")); //right hand side of hard constraint     	
	        try {
	            Files.lines(filePath).forEach(line -> {
	                String[] columnValues = line.split("\\s+");
	                bounds.put(Integer.parseInt(columnValues[0]), Integer.parseInt(columnValues[1]));
	            });
	        } catch (IOException e) {
	            e.printStackTrace();
	        }  
	        
	    	Path myFilePath = Paths.get(DATADIR).resolve(Paths.get("upperBound.txt"));    //upper bound for the open neighborhood   	
	        try {
	            Files.lines(myFilePath).forEach(line -> {
	                String[] columnValues = line.split("\\s+");
	                upperBounds.put(Integer.parseInt(columnValues[0]), Integer.parseInt(columnValues[1]));
	            });
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
	        
			
	         IloIntVar[] x = new IloIntVar[numNodes]; //center node
	         IloIntVar[] y = new IloIntVar[numNodes]; //leaf node	        
	         IloNumVar[] estObj = new IloNumVar[numNodes]; // contribution of each node 	        
	         
	        for (int i = 0; i < numNodes; i++) { //degree of each node   		
	        	degree.put(i, adjacencyList.get(i).size());
	        }
	
			for(int i=0; i< numNodes; i++) { //Define the variables
	        	x[i] = masterCplex.boolVar("x_"+ i);
	        	masterCplex.add(x[i]);
	        	y[i] = masterCplex.boolVar("y_"+ i);
	        	masterCplex.add(y[i]);
	        	estObj[i]= masterCplex.numVar(0.0, degree.get(i),"t_"+ i );  
	        	masterCplex.add(estObj[i]);
	        }
			
			masterCplex.addMaximize(masterCplex.sum(estObj), "obj");
	  
	        
	        IloLinearNumExpr expr = masterCplex.linearNumExpr();
	        
	        for (int i = 0; i < numNodes; i++) { //constraint to be a leaf node
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
	        	masterCplex.addLe(y[i],expr, "leaf_"+i);		                	
	        }
	        expr.clear();
	        
	        for (int i = 0; i < numNodes; i++) { //upper bound for the size of open neighborhood
	        	expr.addTerm(x[i], upperBounds.get(i));		                	
	        }
	        masterCplex.addLe(masterCplex.sum(estObj) , expr, "upper_Bound");         
	        expr.clear();
	        
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
	            	masterCplex.addLe(expr, masterCplex.prod(bounds.get(i) , masterCplex.diff(1, y[i])), "noEdgeShared_"+i);		                	
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
	            	masterCplex.addLe(expr, masterCplex.prod(degree.get(i) , masterCplex.diff(1, y[i])), "fourth_"+i);		                	
	            }
	        }
	        
	       expr.clear();	
	       
	      for(int i=0; i< numNodes; i++)  //only one node can be center
	    	expr.addTerm(x[i], 1);  
	     
	      masterCplex.addEq(expr, 1, "Fifth");
	      expr.clear();
	        
	         
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
	                 TreeSet<Integer> treeSet = new TreeSet<Integer>(leaves); 
	                   stars.put(Integer.valueOf(columnValues[0]),treeSet );
	               
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
	        	  masterCplex.addMIPStart(vars, vals);
	        	}
	      }
	      
	        int numThreads = masterCplex.getNumCores();
	        
	         final BendersPPINCallback cb = new BendersPPINCallback(x,y,estObj, adjacencyList, numThreads);
	         long contextmask = IloCplex.Callback.Context.Id.Candidate
	            | IloCplex.Callback.Context.Id.ThreadUp
	            | IloCplex.Callback.Context.Id.ThreadDown;
	         
	         if(seperateFractional)
	                 contextmask |= IloCplex.Callback.Context.Id.Relaxation;
	         masterCplex.use(cb, contextmask);
	
		        long start = System.currentTimeMillis();
		        if (masterCplex.solve()) {
		        	//save any information you would like to analyze
		     		writer.write((System.currentTimeMillis() - start) * MSTOSEC + ",");	
		            writer.write("\n");
		            writer.flush();
	                masterCplex.end();
	                
		        }else {
		        	writer.write("Solution status: " +   masterCplex.getStatus() + "\n");
		        	masterCplex.end();
		        }
		}finally {
	         masterCplex.end();
	      }
	
		writer.close();
		System.out.println("Experiments are over!");
	}
	
	
}
