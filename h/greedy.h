/**
 * Author: Serigne Gueye
 * Date: 17 July 2026
 */

/****************************************************************/
/* Serigne Gueye						*/
/****************************************************************/

//////////////////////////////////////////////////////////////////////////////////
// But : Greey Heuristic class                                                  //
//////////////////////////////////////////////////////////////////////////////////
class Greedy
{
	IloNumArray L; // Location of each facility
	char filename[100];
	IloInt ncap;  // Number of capacites per facility
	// Parameters
	IloNum alpha;
        IloNum beta;
	IloNum gamma;
	IloInt m;       // Number of facilities
	IloNum VOPT;    // Optimal value
	
	public :
	// Constructor
	Greedy(const IloEnv &  env , char* fichier_arg, IloInt ncap_arg, IloNum alpha_arg, IloNum beta_arg, IloNum gamma_arg, IloInt m_arg); 	
	void GreedyHeuristic(IloNumArray Cm); // Greedy heuristic 
	void Localsearch(IloNumArray Cm); 
	IloNum Opt();
	IloNum Centralization(FloatMatrix profit);
	IloNum Colocation(const IloEnv & env);
	void AffichageSolutions(char* filename);
};

