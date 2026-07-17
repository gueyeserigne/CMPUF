/**
 * Author: Serigne Gueye
 * Date: 17 July 2026
 */

/****************************************************************/
/* 21 Avril 2024						*/
/* Problème de Transport
/* Serigne Gueye						*/
/* Projet IEF FLOUE MIT-ZLC					*/
/****************************************************************/


//////////////////////////////////////////////////////////////////////////////////
// But : Classe Transport							//
//////////////////////////////////////////////////////////////////////////////////

class Transport: public IloModel
{
	FloatVarMatrix zVar;			// variables z
	FloatVarMatrix xVar;			// variables z
	FloatMatrix z;				// valeurs des variables z
	FloatMatrix x;				// valeurs des variables x
	IloNum VOPT;				// Valeur Optimale

	IloInt n; 				// Nombre de lignes de x
	IloInt m; 				// Nombre de colonnes de x
	IloObjective obj;
	IloNumArray h;
	FloatMatrix p;
	IloNumArray c; // Capacité des services

	public :
	
	Transport(const IloEnv  & env, IloInt n_arg, IloInt m_arg, FloatMatrix p_arg, IloNumArray h_arg, IloNumArray c_arg); 
	IloNum Opt();						
	void Init_z(const IloEnv  & env);
	void Init_x(const IloEnv  & env);
	void Init_contraintes(const IloEnv  & env);
	void Init_obj(const IloEnv  & env);
	void Resolution(const IloEnv & env, IloCplex cplex);
	IloNumArray Init_L(const IloEnv & env);
	IloNumArray Init_Supply(const IloEnv & env, IloNumArray L, IloNum prob);
};
