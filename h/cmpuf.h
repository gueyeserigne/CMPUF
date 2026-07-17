/**
 * Author: Serigne Gueye
 * Date: 17 July 2026
 */

/****************************************************************/
/* Projet IEF FLOUE MIT-ZLC					*/
/****************************************************************/
/* MPUF : Unreliable Capacitated Facility Location Problem
/****************************************************************/
/* Abtract : 
In this code we aim at solving the problem of locating facilities knowing that their
respective capacities may be disrupted and, according to some probability distribution, 
capacity may be reduced from its theoretical (planned) value. Service facilities can and usually, 
at some point in time, do fail. Failures may happen with high frequency,
*/
/****************************************************************/
/*  We have:                                                    */
/****************************************************************
- n  : numer of clients. Each client has a demand w[i]. w is the vector of the client demands.
****************************************************************/
/****************************************************************
- H is the total demand. H = w[1]+w[2]+...+w[n]
****************************************************************/
/****************************************************************
- m facilities : each facility i may have "ncap" different capacity ncap[i][j], j=1,2..,ncap. These capacity are proportional
to Cm[i] (a value computed by solving a m-median problem).
-----------------------------------------------
Details :  
We call "median capacity" of a facility the capacity corresponding to the sum of the demand assigned to the facility in a classical m-median facility 
location problem. We first solve a m-median problem where facility has unlimited capacity. The sum of the demand assigned to the facility
i  is tored in Cm[i]. Then in our unreliable facility location problem the possible capacity of each facility i depend on Cm[i] and on the parameters in  
wbeta.
****************************************************************/
/****************************************************************
- S : For each facility i, S[i] = gamma*Cm[i]
****************************************************************/
/****************************************************************
- ncap : number of possible capacities per facility
****************************************************************/
/****************************************************************
- p is the vector of the capacity probabilities. For each facility the probability that it has the capacity cap[i][j] ( j = 1,2,..,ncap) is p[j].
****************************************************************/
/****************************************************************
- cap : matrix of possible capacities. For each facility index i = 1,2,...,n and each capacity index j = 1,2,...,ncap-1, cap[i][j] = wbeta[j]*S[i]
and cap[i][ncap-1] = S[i]
****************************************************************/
/****************************************************************
- K  = ncap^m is the total number of capacity disruption scenarios. 
The capacity disruption scenarios is defined as the cartesian product of all possible capacity of each facility. 
****************************************************************/
/****************************************************************
- C : Scenario capacity matrix. C has K lines. Each index scenario k and each facility j, C[k][j] is the capacity of j in the scenario  k.
The probability of realization of the scenatio k is the product of the probability of realization of each facility capacity of the scenario.
****************************************************************/


//////////////////////////////////////////////////////////////////////////////////
// But : Classe representant le modele lineaire de resolution du probleme MPUF  //
//       (voir rapport technique pour le modèle       				//
//										//
// Attribut :                                        				//
// Attribut :                                        				//
//////////////////////////////////////////////////////////////////////////////////

class CMpuf : public IloModel
{
	IloInt K;                   // number of scenarios = ncap^m
	IloInt m;                   // number of facilities
	IloInt n;                   // number of clients
	IloNum H;                   // sum of client demands
	IloInt ncap;                // number of capacities
	FloatMatrix C;              // Scenario capacity matrix
	IloNumArray S;              // Vector des capacités par scénarios
	IloNumArray Cm;             // Cm[i]  = sum of the demands of the clients assigned to the facility i in the optimal solution of the uncapacittated
	                            // m-median problem
	IloNumArray p;              // Probability vector. p[j] = probability of any facility i to have the capacity cap[i][j]  
	FloatMatrix cap;            // Capacity Matrix
	IloNumArray probaScenarios; // Vecteur scenario probabilities
	IloNumArray w;              // vecteur of client demands
	IloNum VOPT;                // Optimal value of the the model coded in model.cpp
	IloNum VOPTTRANSPORT;       //
	IloNum VGREEDY;             // Optimal value of the greedu heuristic solution
	IloNum VLS;                 //  Optimal value of a local search procedure	
	IloNum alpha;             // Parameter used to compute p
	IloNumArray walpha;       // Weights used in p computation
	IloNum gamma;             // Parameter used to compute S
	IloNum beta;              // 
	IloNumArray wbeta;        // Parameter array used to compute cap
	FloatMatrix T;            // Matrix of capacity indices of each scenario. For each scenario k = 0,1...,K-1  , T[k][j] the index (i = 0,1,..,ncap-1) of the capacity 
	                          // of the facility j = 0,1,..,m-1 in the scenario k.
	FloatMatrix MatriceCoutTransport;
	FloatMatrix SolScenarios; // Matrix containing in each line k = 0,1...,K the optimal facility location solution of the scenario k.
	FloatMatrix SolSupply;    // Matrix containing in each line k = 0,1...,K the quantity S[i] that the facility i must supply to its associated clients.
	FloatMatrix V;
	  FloatVarMatrixMatrix zVar;  // Variables z of the model. zVar[l][j][k] = qauntity supply to the client demand in j by a facility located in l in the scenario k
	                            // and 0 otherwise
	FloatMatrixMatrix z;        // Values of z
        FloatVarMatrix xVar;        // Variables x. x[i[l] = 1 if the facility is located on the node l, 0 otherwise.
	FloatMatrix x;              // Values x
	FloatMatrix profit;         // Matrice des profits. The "profit" red is the distance between nodes
	IloRangeArray dcon;         // Liste of supply constraints. For each cllent j and scenario k dcon contain a supply constraint.
	IloRangeArray ct;           // Local Search constraint
	IloNumArray d;              // Slack values for the constraints in dcon

	public :
	
	CMpuf(const IloEnv & env , char* fichier, std::string type, char* param, IloNum alpha_arg, IloNum beta_arg, IloNum gamma_arg, IloInt m_arg); 	// Constructeur de la classe.
	void Init_cap(const IloEnv  & env, IloNumArray Cm);                 // cap initialization
	void Init_z(const IloEnv  & env);                                   // z initialization
	void Init_x(const IloEnv  & env);                                   // x initialization
	void Init_xpmedian(const IloEnv & env, FloatMatrix & xpmedian);   // Init xpmedian to the values computed in x
	void Init_x(IloNumArray L);                                         // x initialization with location in L
	void Init_contraintes(const IloEnv  & env);                         // Model constraints
	void Init_contraintes_x(const IloEnv & env);                        // Model constraints on location variables x
	void Init_contraintes_x(const IloEnv & env, int i);                 // Model constraints on location variables x of the facility index i
	void Init_contraintes_x_x(const IloEnv & env);                      // Model constraints forbidding co-location
	void Init_obj(const IloEnv & env);                                  // Objective function initialization
	void Init_C(const IloEnv & env);                                    // Matrix C initialization
        IloNumArray Init_L(const IloEnv & env);                             // Initialization of L with the content of x 
	IloNumArray Init_L(const IloEnv & env,  FloatMatrix & xpmedian);    // Initialization of L with the content of xpmedian
	void Init_T(const IloEnv  & env, int i, IloNumArray numeros);       // Recursive method to initialized T
	///////////////////////////////////////////////
	void AffichageDonnees();                                            // Attributes printing                                       
	void AffichageSolutions();                                          // Solution printing
	void AffichageSolutions(IloNumArray L, char* filename);             // Writing solution  Lin the file filename
	void AffichageSolutionsScenarios(IloNumArray L,char* filename);     // Writing solution  Lin the file filename
	void Solution_z(IloInt K_prec, IloCplex cplex);                     // Checking purpose of what variable is extracted
	void Conversion(const IloEnv & env);                                // Conversion of variables x to integer
	IloNum Resolution(const IloEnv & env, IloCplex cplex);              // Running cplex optimize
	IloNum Resolution(const IloEnv & env, IloCplex cplex, FloatMatrix sol); // Fixing the variables x to the location given in sol. And running cplex
	void AffichageSeparateurs(IloInt taille);                               // Print some separating characters
	///////////////////////////////////////////////
	// List of getter. Method to access to attribute
	///////////////////////////////////////////////
	FloatMatrix Profit();                   // Attribute profit
	IloNumArray Capacites(IloInt k);        // return C[k]  
	FloatMatrix LireT();                    // return T
	IloNumArray W();                        // return w
	IloNum Alpha();                         // return alpha
	IloNum Beta();                          // return beta
	IloNum Gamma();                         // return gamm
	IloNum Opt();                           // return VOPT
	IloNum Greedy();                        // return VGREEDY
	IloNum LS();                            // return VLS
	IloNum Pmedian();                       // return VOPTTRANSPORT
	IloInt M();                             // return m
	IloInt LireK();                         // return K
	FloatMatrix Lirex();                    // return x
	IloInt N();                             // retunr n

	IloNum Centralization(IloNumArray L, FloatMatrix profit);         // Compute centralization value of L
	IloNum Colocation(const IloEnv & env, IloNumArray L);             // Compute colocation value of L
	

	void Fixation_x(const IloEnv & env, IloNumArray Fixed);           // Fixed contains facility location. The variables x are fixed with Fixed
	void Fixation_x(IloNumArray Fixed, int i);                        // Fix facility i  to the location Fixed[i]
	void Defixation_x(IloNumArray Fixed, int i);                      // Defix facility i  to the location Fixed[i]
	void Defixation_x(IloNumArray Fixed);                             // Defixed facility location.
	void Fixation_z(FloatMatrixMatrix z_prec, IloInt K_prec);         // Fix variables z to the values in z_prec for the first K_prec scenarios
	IloInt Localisation(IloInt i);                                    // return the location of the facility i
	FloatMatrixMatrix Z(const IloEnv & env);                          // Build and return a copy of the matrix z
	void Destructeur(const IloEnv & env);                             // Destructor
	IloNumArray Localsearch(const IloEnv & env, IloCplex cplex, IloNumArray L, IloNum VL, IloNum rayon);  // Local Search algorithm
	void Annotations(IloCplex cplex);                                 // Annotations of the variables used for the benders decompistion scheme
	void Zero_x(int i, int l);                                        // Fix to 0 the variable xVar[i][l];
	void Zero_x(int i);                                               // Fix to 0 the line  xVar[i][l];
	void Zero_x();                                                    // Fix to 0 xVar
	void Un_x(int i, int l);                                          // Fix to 1 the variable xVar[i][l];
	void Zero_Un_x(int i, int l);                                     // Set to 0 the lower bound of xVar[i][l] and to 1 its upper bound      
	void Zero_Un_x(int i);                                            // Run Zero_Un_x(int i, int l) in the line xVar[i]
	void Zero_Un_x();                                                 // Run Zero_Un_x(int i) in the xVar
	IloNumArray GreedyHeuristic(const IloEnv & env, IloCplex cplex);  // Greedy Heuristic
	IloNumArray MyGreedyHeuristic(const IloEnv & env, IloCplex cplex);  // Greedy Heuristic using benders decompoosition annotations
	void Solution(FloatMatrix sol);                                     // Set xVar with the solution in sol
	void MIPStart(const IloEnv env, IloCplex cplex, IloNumArray L);     // addMIPStart with the location L   
	void AddLocalConstraint(const IloEnv & env, IloNumArray L, IloInt rayon); // Add constraint for the local search
	void RemoveLocalConstraint();                                             // Remove local constraint
	IloNum LocalBranchAndBound(const IloEnv & env, IloCplex cplex, IloNumArray L, IloNum rayon); // Branch-and-bound with local constraint
	IloNum CoutTransport(IloInt k);                                           // Compute the transportation cost with the assignment z of the scenatio k
	void AffichageCoutTransport(const IloEnv & env, IloNum cent, char* filename); // Write int the file filename beta, gamma, VOPT and centraliation
	void AffichageCoutTransport();                // Print probaScenarios, transportation cost
	void ResolutionScenarios(const IloEnv & env); // Solve a m-median problem for each scenario k
	void GenerationSolutionsScenarios(const IloEnv & env, IloNumArray L,char* filename); // Wirte in the filename ..... (to finsh) 
};

