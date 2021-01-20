# SDC
This repository contains the code for the manuscript entitled "The Star Degree Centrality Problem: A  Decomposition Approach", by Mustafa C. Camur,  Thomas C. Sharkey, and Chrysafis Vogiatzis. The paper has been accepted at INFORMS Journal on Computing and a preprint may be downloaded at:

http://www.optimization-online.org/DB_HTML/2020/03/7694.html

In each data file, you can reach the adjacency list as a txt file. In addition, we share "bound.txt" file which can be used as a constraint tightening in NIP formulation. 

NIP.java : Solves the problem with the new IP formulation proposed by using CPLEX solver. 
VCIP.java : Solves the problem with the existing IP formulation in the literature by using CPLEX solver.
BendersNIP.java : Solves the problem by using branch-and-cut approach motivated by Benders Decomposition on NIP formulation. 
BendersVCIP.java : Solves the problem by using branch-and-cut approach motivated by Benders Decomposition on VCIP formulation. 
UpperBound.java : Calculates an upper bound on the objective value for each node in the network. 
RatioBasedHeuristic.java : Determines a feasible star for each node through the ratio-based heuristic.
