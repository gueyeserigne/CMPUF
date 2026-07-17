/**
 * Author: Serigne Gueye
 * Date: 17 July 2026
 */

/****************************************************************/
/* Serigne Gueye						*/
/* Projet IEF FLOUE MIT-ZLC					*/
/****************************************************************/


//////////////////////////////////////////////////////////////////////////////////
// But : Classe pmedian								//
//////////////////////////////////////////////////////////////////////////////////
class Pmedian
{
	FloatVarMatrix zVar;			// variables z
	FloatMatrix z;				// valeurs des variables z
	IloNumVarArray XVar;			// variables x x
	IloNumArray X;				// solution x
	IloNum VOPT;				// Optimal value
	IloInt n; 				// Number of clients
	IloInt m; 				// Number of facilities
	IloObjective obj;                       // Objective function
	IloModel model;                         // ILOG Model
	IloCplex cplex;                         // Cplex object
	IloNumArray h;                          // Client demands vector
	FloatMatrix profit;                     // Profit matrix (distances)
	FloatMatrix xpmedian;                   // Optimal solution

	public :
	
	Pmedian(const IloEnv  env, char* fichier, std::string type, IloInt m_arg);        // Conqstructor 
	IloNum Opt();						        // Get optimal value
	void Init_z(const IloEnv  env);                                 // z constructor
	void Init_X(const IloEnv  env);                                 // x constructor
	void Init_contraintes(const IloEnv  env);                       // Build model constraints 
	void Init_obj(const IloEnv  env);                               // Build objective function
	void Resolution();                                              // Solve the problem 
	FloatMatrix XtoXpmedianSolution(const IloEnv env);              // Convert location vector X to the matrix form xpmedian
	void Conversion(const IloEnv  env);                             // Convert X to binary variables
	IloNumArray CapaciteMedian(const IloEnv  env);                  // Return an array of the sum of the demands assigned to each facility in the optimal 
	                                                                // m-median solution 
};
