/**
 * Author: Serigne Gueye
 * Date: 17 July 2026
 */

/****************************************************************/
/* Serigne Gueye						*/
/* Projet IEF FLOUE MIT-ZLC					*/
/****************************************************************/


//////////////////////////////////////////////////////////////////////////////////
// But : Classe representant le modele lineaire de resolution du probleme MPUF  //
//       (voir rapport technique pour le modèle       				//
//										//
// Attribut :                                        				//
// Attribut :                                        				//
//////////////////////////////////////////////////////////////////////////////////

class Greedy : public IloModel
{
	IloInt K; // Nombre de scenarios
	IloInt m; // Nombre de services
	IloInt n; // Nombre de clients
	IloInt ncap; // 
	FloatMatrix C; // Matrice des capacités par scénarios
	IloNumArray p; // Vecteur des probabilités
	IloNumArray cap; // Vecteur des capacites
	IloNumArray probaScenarios; // Vecteur des probabilités
	IloNumArray w; // vecteur des demandes
	IloNum VOPT;
	IloNum alpha;
	IloNum gamma;
	IloNum beta;
	FloatMatrix T;

	IloNumArray L;
	
	FloatVarMatrixMatrix zVar; // Variables distance
	FloatMatrixMatrix z; // Variables distance

	FloatVarMatrix xVar; // Variables distance
	FloatMatrix x; // Variables distance

	FloatMatrix profit; // Matrice des profits

	public :
	
	Greedy(const IloEnv  env , char* fichier, IloInt ncap_arg, IloNum alpha_arg, IloNum beta_arg, IloNum gamma_arg); 	// Constructeur de la classe.

	void Init_z(const IloEnv  env);

	void Init_x(const IloEnv  env);

	void Init_contraintes(const IloEnv  env);

	void Init_obj(const IloEnv env);

	void Init_C(const IloEnv env);

	void Conversion(const IloEnv env);

	void AffichageDonnees();

	void AffichageSolutions();

	void Resolution(const IloEnv env, IloCplex cplex);

	IloNum Alpha();

	IloNum Beta();

	IloNum Gamma();

	IloNum Opt();

	IloInt M();

	void Init_L(const IloEnv env);

	void Init_T(const IloEnv & env, int i, int * numeros);

	void GreedyHeuristic(const IloEnv env, IloCplex cplex);

	void GreedyHeuristicNew(const IloEnv env, IloCplex cplex);

	void AffichageSeparateurs(IloInt taille);

	IloNum Centralization();

	IloNum Colocation(const IloEnv env);
};

