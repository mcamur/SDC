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


public class BendersVCIP {
	public static final double MSTOSEC = 0.001;
	public static final double FUZZ = 1e-6;
	public static final boolean seperateFractional = true; 
	public static final boolean warmStart = false; 
	public static final String graphType = "..."; //Specify the graph type for DATADIR
	public static final double GAMMA = 1; 
	private static final class Neighbor {
	      private final int numNodes;
	      private final HashMap<Integer, HashSet<Integer>>  adjacencyList;
	      
		public Neighbor(int numNodes, HashMap<Integer, HashSet<Integer>>  adjacencyList) throws IloException {
	         this.numNodes = numNodes;
	         this.adjacencyList = adjacencyList;			 
     	          
		}

		public IloRange[] separate(IloNumVar[] y, IloNumVar[] estObj,  double[] ySol, double[] zMaster) throws IloException {
			IloRange[] cut = new IloRange[numNodes];
		    double[] alpha = new double[numNodes];
     		double[] beta =new double[numNodes];
     		double[] sum = new double[1];    //to solve the problem with an algorithm rather than using CPLEX
     		double[] dualSol = new double[numNodes]; //objective of dual problem
     		
     		IloCplex cplex = new IloCplex();	
		      for(int i =0; i< numNodes; i++) {
		    	  sum[0] =0;
	        		adjacencyList.get(i).forEach((val) -> { 	        			
	        			 sum[0] += (ySol[val]);	        			 
	        		});
		    	  if( (1 - ySol[i]) > 0.5 ) {
		        		if( (1-  ySol[i]) >   sum[0] +FUZZ ) {
		        			alpha[i] = 0;
		        			beta[i] = 1;
		        		}else if((1- ySol[i]) <   sum[0] -FUZZ ) {
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
		    	  dualSol[i]  += (1- ySol[i])* alpha[i];
	        		adjacencyList.get(i).forEach((val) -> { 
	        			 sum[0] +=  ySol[val];
	        		});
	        		dualSol[i] +=  (sum[0] * beta[i]);
		      }
		      

        	   for (int i = 0; i < numNodes; i++) {
  		        	 if (zMaster[i] > dualSol[i]  + FUZZ) {	
  		        		IloLinearNumExpr expr = cplex.linearNumExpr();
  	        		 	expr.setConstant(alpha[i]);	
  	     	        	expr.addTerm(-alpha[i],y[i]);
  	     	        	int node = i;
  	     	        	adjacencyList.get(i).forEach((val) -> { 
      	        			try {
      							expr.addTerm(beta[node], y[val]);
      						} catch (IloException e) {		
      							e.printStackTrace();
      						}
      	        		});	 	    	        		  
  		     	       cut[i] = (IloRange) cplex.le(estObj[i],expr);
  		     	       expr.clear();
  		        	 }
        	   }
	        cplex.end();
			return cut;
		}
	}
	
	
	 private static final class BendersPPINCallback implements IloCplex.Callback.Function {
	      private final IloIntVar[] y;
	      private final IloNumVar[] estObj;
	      private final Neighbor[] neighbors;
	      private final HashMap<Integer, HashSet<Integer>> adjacencyList;
		  ;
	      public BendersPPINCallback(IloIntVar[] y,  IloNumVar[] estObj, HashMap<Integer, HashSet<Integer>> adjacencyList, int numThreads) {
	    	  this.y = y;
	    	  this.estObj = estObj;
	    	  this.neighbors = new Neighbor[numThreads];
	    	  this.adjacencyList = adjacencyList;
	      }
	      
		@Override
		public void invoke(IloCplex.Callback.Context context) throws IloException {
			try {
			int threadNo = context.getIntInfo(IloCplex.Callback.Context.Info.ThreadId);
			int numNodes = y.length;
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
	         
	         double[] ySol = new double[numNodes];
	         double[] zMaster =new double[numNodes];
	         
	         if ( context.inCandidate() ) {
	             if ( !context.isCandidatePoint() ) // The model is always bounded
	                throw new IloException("Unbounded solution");
	             for (int i = 0; i < numNodes; ++i) {	            	 
	                ySol[i] = context.getCandidatePoint(y[i]);	 
	                zMaster[i] = context.getCandidatePoint(estObj[i]);	                
	             }	                
	          }else if( context.inRelaxation() ) {
                  for (int i = 0; i < numNodes; ++i) {
  	                ySol[i] = context.getRelaxationPoint(y[i]);	 
  	                zMaster[i] = context.getRelaxationPoint(estObj[i]);	
                }
             
              }else {
	              throw new IloException("Unexpected contextID");
	          }
	         
	         Neighbor neighbor = neighbors[threadNo];
	         IloRange[] violated = new IloRange[numNodes];
	         violated = neighbor.separate(y,estObj,ySol, zMaster);
	         for(int i =0; i< numNodes ; i++) {
	         if (violated[i] != null) {
	             // Add the cut	        	
	             if ( context.inCandidate() ) {
	            	 context.rejectCandidate(violated[i]);
	             }else if( context.inRelaxation() ){ 
	                 context.addUserCut(violated[i],
                             IloCplex.CutManagement.UseCutPurge,
                             false);
	             }    
	             else
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

		final IloCplex masterCplex = new IloCplex();		
		masterCplex.setParam(IloCplex.Param.TimeLimit, 2*1800 - heurTime);
		masterCplex.setParam(IloCplex.Param.MIP.Strategy.VariableSelect, 3); //Strong branching
		masterCplex.setParam(IloCplex.Param.Emphasis.MIP, 2); // Optimality over Feasibility
		masterCplex.setParam(IloCplex.Param.MIP.Strategy.RINSHeur, 1000);
		masterCplex.setOut(null);
			
			try {							
				HashMap<Integer, Integer> degree = new HashMap<>();
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
		      		        
		        final IloIntVar[] x = new IloIntVar[numNodes]; //center node
		        final IloIntVar[] y = new IloIntVar[numNodes]; //leaf node	        
		        final IloNumVar[] estObj = new IloNumVar[numNodes]; // contribution of each node 
		        for (int i = 0; i < numNodes; i++) //degree of each node  		
		            	degree.put(i, adjacencyList.get(i).size());

				for(int i=0; i< numNodes; i++) { //Define the variables
		        	x[i] = masterCplex.boolVar("x_"+ i);
		        	masterCplex.add(x[i]);
		        	y[i] = masterCplex.boolVar("y_"+ i);
		        	masterCplex.add(y[i]);
		        	estObj[i]= masterCplex.numVar(0.0, degree.get(i),"t_"+ i );  // a node can contribute to objective with at most 1 but leaving this loose gives better results.
		        	masterCplex.add(estObj[i]);
		        }
				
				masterCplex.addMaximize(masterCplex.sum(estObj), "obj");

		        IloLinearNumExpr expr = masterCplex.linearNumExpr();
		        
		        for (int i = 0; i < numNodes; i++) { // you should be adjacent to the center node in order to be in star! 
		    		expr.clear();	
		    		adjacencyList.get(i).forEach((val) -> { 
		    			try { 
		    				expr.addTerm(x[val], 1);
		    			}catch (IloException e) {
							e.printStackTrace();
						} 
		    		});
		    		expr.addTerm(x[i], 1); // Closed neighborhood
		    		masterCplex.addLe(y[i],expr, "third_"+i);		                	
		        }
		        expr.clear();
		        
		        for(int i=0; i< numNodes; i++)  //Missing constraint in the original formulation. The center must be in star.
		        	masterCplex.addLe(x[i], y[i]);
		        
		        for (int i = 0; i < numNodes; i++) { //upper bound for the size of open neighborhood
		        	expr.addTerm(x[i], upperBounds.get(i));		                	
		        }
		        masterCplex.addLe(masterCplex.sum(estObj) , expr, "upper_Bound");
		       expr.clear();
		       
		       for (int i = 0; i < numNodes; i++) { //no leaf is adjacent to another leaf. 
		       	int node = i;
		       		adjacencyList.get(i).forEach((val) -> { 
		       			if(val > node) {
		       			try { 
		       				expr.clear();
		       				expr.addTerm(y[node], 1);
		       				expr.addTerm(y[val], 1);
		       				expr.addTerm(x[node], -1);
		       				expr.addTerm(x[val], -1);
		       				masterCplex.addLe(expr,1, "fourth_"+node);	
		       			}catch (IloException e) {
		   					e.printStackTrace();
		   				} 
		       		 }	
		       		});      		                	
		       }
		       expr.clear();
		       
		      for(int i=0; i< numNodes; i++)  //only one node can be center
		    	expr.addTerm(x[i], 1);  
		     
		      masterCplex.addEq(expr, 1, "Fifth");
		      
		      expr.clear();
		        
      
		        int numThreads = masterCplex.getNumCores();
		        
		         final BendersPPINCallback cb = new BendersPPINCallback(y,estObj, adjacencyList, numThreads);
		         long contextmask = IloCplex.Callback.Context.Id.Candidate
		            | IloCplex.Callback.Context.Id.ThreadUp
		            | IloCplex.Callback.Context.Id.ThreadDown;
		         if(seperateFractional)
	                 contextmask |= IloCplex.Callback.Context.Id.Relaxation;
		         
		        if(warmStart) {		        	 
		      		TreeMap<Integer, Set<Integer>> stars = new TreeMap<Integer, Set<Integer>>();  
		      		 myFilePath = Paths.get(DATADIR).resolve(Paths.get("stars.txt"));      	
		             try {
		                 Files.lines(myFilePath).forEach(line -> {
		                     String[] columnValues = line.split("\\s+");
		                     Set<Integer> yVar = new HashSet<Integer>();
		                     for(int i =0 ; i < columnValues.length; i++) {
		                    	 yVar.add(Integer.valueOf(columnValues[i]));
		                     }
		                     TreeSet<Integer> treeSet = new TreeSet<Integer>(yVar); 
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
		          
		         masterCplex.use(cb, contextmask);
		        //masterCplex.exportModel("multiCutMasterProblem.lp");
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
