/**
 * Author: Serigne Gueye
 * Date: 17 July 2026
 */

/****************************************************************/
/* Serigne Gueye						*/
/* Projet IEF FLOUE MIT-ZLC					*/
/****************************************************************/

//#include <ilcplex/ilocplex.h>
//#include <ilsolver/ilosolverint.h>

#include <ilcplex/ilocplex.h>
//#include <ilcp/cp.h>
#include <iostream>
#include <string>

typedef IloArray<IloModel> IloModelArray;
typedef IloArray<IloCplex> IloCplexArray;
typedef IloArray<IloObjective> IloObjectiveArray;
typedef IloArray<IloRangeArray> RangeMatrix;
typedef IloArray<IloNumArray> FloatMatrix;
typedef IloArray<FloatMatrix> FloatMatrixMatrix;
typedef IloArray<FloatMatrixMatrix> FloatMatrixMatrixMatrix;
typedef IloArray<FloatMatrixMatrixMatrix> FloatMatrixMatrixMatrixMatrix;
typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<IloBoolArray> BoolMatrix;
typedef IloArray<IloNumVarArray> FloatVarMatrix;
typedef IloArray<FloatVarMatrix> FloatVarMatrixMatrix;
typedef IloArray<FloatVarMatrixMatrix> FloatVarMatrixMatrixMatrix;
typedef IloArray<FloatVarMatrixMatrixMatrix> FloatVarMatrixMatrixMatrixMatrix;
typedef IloArray <FloatMatrix> FloatMatrixMatrix;
typedef IloArray<IloInt> IntVector;
typedef IloArray<IntVector> IntMatrix;
typedef IloArray <IntMatrix> IntMatrixMatrix;
typedef IloArray<IloRangeArray> IloRangeArray1;
typedef IloArray<IloRangeArray1> IloRangeArray2;
typedef IloArray<IloRangeArray2> IloRangeArray3;

#include "../h/greedy.h"
#include "../h/cmpuf.h"
#include "../h/tri.h"
#include <time.h>

#define INFINI 1e+30
#define EPS 1e-10

ILOSTLBEGIN

//////////////////////////////////////////////////////////////////
// Constructor								//
//////////////////////////////////////////////////////////////////
Greedy::Greedy(const IloEnv & env , 
		char* fichier_arg, 
		IloInt ncap_arg, 
		IloNum alpha_arg, 
		IloNum beta_arg, 
		IloNum gamma_arg, 
		IloInt m_arg)
{

	//cout << "debut Greedy" << endl;
	//cin.get();  

	strcpy(filename,fichier_arg);
	alpha = alpha_arg;
	gamma = gamma_arg;
	beta = beta_arg;
	ncap = ncap_arg;
	m = m_arg;

	L = IloNumArray(env);

	//cout << "Fin Greedy" << endl;
	//cin.get();  
}
//////////////////////////////////////////////////////////////////
// Return VOPT
//////////////////////////////////////////////////////////////////
IloNum Greedy::Opt()
{
	return(VOPT);
}
//////////////////////////////////////////////////////////////////
// Greedy Heuristic
//////////////////////////////////////////////////////////////////
void Greedy::GreedyHeuristic(IloNumArray Cm)
{
	IloInt i;
	FloatMatrixMatrix z_prec;
	IloInt K_prec = 0;

	//cout << "Begin Greedy Heuristic..." << endl;
	//cin.get();

	//IloEnv env0;

	for(i = 0 ; i < m ; i++) // For each facility
	{
		IloEnv env;		
		//cout << "Placement service = " << i << endl;
		//cin.get();
		CMpuf cmpuf(env,filename,2,alpha,beta,gamma,i+1);
		//cout << "cmpuf" << endl;
		//cin.get();
		//CMpuf cmpuf(env,filename,3,alpha,beta,gamma,i+1);
		cmpuf.Init_cap(env,Cm);
		//cout << "Init_cap" << endl;
		//cin.get();

		//cout << cmpuf.M()+1 << endl;

		IloNumArray numeros(env,cmpuf.M());

		//cout << "Begin init_T" << endl;
		cmpuf.Init_T(env,i,numeros);
		//cout << "End init_T" << endl;

		//cout << "Begin init_C" << endl;
		cmpuf.Init_C(env);
		//cout << "End init_C" << endl;

		//cmpuf.AffichageDonnees();
		cmpuf.Init_x(env);
		cmpuf.Fixation_x(env,L);
		cmpuf.Init_z(env);

		//cmpuf.Fixation_z(z_prec,K_prec);

		//env0.end();
		//z_prec.end();

		cmpuf.Init_contraintes(env);
		cmpuf.Init_obj(env);
		
		IloCplex cplex(env);
		cplex.extract(cmpuf); 

		cmpuf.Conversion(env);

		cplex.setParam(IloCplex::HeurFreq,-1);
		cplex.setParam(IloCplex::FlowPaths,-1);
		cplex.setParam(IloCplex::FlowCovers,-1);
		cplex.setParam(IloCplex::CutPass,-1);
		cplex.setParam(IloCplex::MIPDisplay,0);
		//cplex.exportModel("cmpuf_greedy.lp");
		//cout << "export du model cmpuf_greedy.lp" << endl;
		//cin.get();

		cmpuf.Resolution(env,cplex);
		//cmpuf.AffichageSolutions();
		//cin.get();

		VOPT = cmpuf.Opt();

		//cmpuf.Solution_z(K_prec,cplex);

		L.add(cmpuf.Localisation(i));

		//z_prec = cmpuf.Z(env0);

		//K_prec = cmpuf.LireK();

		/*
		cout << "i = " << i << endl; 
		cout << "L = " << L << endl;
		cout << "VOPT = " << VOPT << endl;
		cout << "------------------------"  << endl;
		cin.get();*/
				
		//cout << "VOPT heuristique= " << VOPT << endl;
		//cin.get();

		//free(numeros);

		//cmpuf.Destructeur(env);
		//cmpuf.end();
		env.end();

	}


	//env0.end();



	//cout << "End Greedy Heuristic..." << endl;
}
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
void Greedy::Localsearch(IloNumArray Cm)
{
	IloBool stop = false;

	while(!stop){	

	IloEnv env;
	stop = true;

	//cout << "ICI1" << endl;
	//cin.get();

	CMpuf cmpuf(env,filename,2,alpha,beta,gamma,m);
	//CMpuf cmpuf(env,filename,3,alpha,beta,gamma,i+1);
	cmpuf.Init_cap(env,Cm);

	//cout << cmpuf.M()+1 << endl;

	//cout << "ICI2" << endl;
	//cin.get();

	IloNumArray numeros(env,cmpuf.M());

	//cout << "Begin init_T" << endl;
	//cin.get();
	cmpuf.Init_T(env,m-1,numeros);
	//cout << "End init_T" << endl;
	//cin.get();

	//cout << "Begin init_C" << endl;
	//cin.get();
	cmpuf.Init_C(env);
	//cout << "End init_C" << endl;
	//cin.get();

	//cmpuf.AffichageDonnees();
	cmpuf.Init_x(env);

	//cout << "Init_x" << endl;
	//cin.get();

	cmpuf.Fixation_x(env,L);

	//cout << "Fixation_x" << endl;
	//cin.get();

	cmpuf.Init_z(env);

	//cout << "Init_z" << endl;
	//cin.get();

	cmpuf.Init_contraintes(env);

	//cout << "Init_cont" << endl;
	//cin.get();

	cmpuf.Init_obj(env);

	//cout << "Init_obj" << endl;
	//cin.get();

		
	IloCplex cplex(env);
	cplex.extract(cmpuf); 

	//cout << "Init_cplex" << endl;
	//cin.get();


	cmpuf.Conversion(env);

	//cout << "conversion" << endl;
	//cin.get();


	cplex.setParam(IloCplex::HeurFreq,-1);
	cplex.setParam(IloCplex::FlowPaths,-1);
	cplex.setParam(IloCplex::FlowCovers,-1);
	cplex.setParam(IloCplex::CutPass,-1);
	cplex.setParam(IloCplex::MIPDisplay,0);
	//cplex.exportModel("cmpuf_local_search1.lp");
	//cout << "export du model local_search1.lp" << endl;
	//cin.get();
	IloNum max = VOPT;
	IloNum tmp;

	IloInt besti;
	IloInt bestli;
	IloInt bestj;
	IloInt bestlj;

	//cout << "VOPT = " << VOPT << endl;
	//cin.get();

		for(int i = 0 ; i < m ; i++){
			//cout << "i = " << i << endl;
			cmpuf.Defixation_x(L,i);
			for(int j = i+1 ; j < m ; j++){
				cmpuf.Defixation_x(L,j);
				//cout << "Defixation " << endl;
				//cplex.exportModel("cmpuf_local_search2.lp");
				//cout << "export du model local_search2.lp" << endl;
				//cin.get();
				cmpuf.Resolution(env,cplex);
				//cmpuf.AffichageSolutions();
				//cin.get();

				tmp = cmpuf.Opt();
				//cout << "tmp = " << tmp << endl;

				if(tmp > max){
					besti = i;
					bestli = cmpuf.Localisation(i);
					bestj = j;
					bestlj = cmpuf.Localisation(j);
					max = tmp;

				}
				cmpuf.Fixation_x(L,j);
			}

			cmpuf.Fixation_x(L,i);
			//cout << "Fixation " << endl;

		}

		if(max > VOPT){
			VOPT = max;
			//cout << "besti " << besti << endl;
			//cout << "bestli " << bestli << endl;
			L[besti] = bestli;
			L[bestj] = bestlj;
			stop  = false;
		}

	//cout << "Fin Localsearch" << endl;
	//cin.get();

	//free(numeros);

//	cmpuf.Destructeur(env);
	//env.end();
	env.end();

	}




}
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
IloNum Greedy::Centralization(FloatMatrix profit)
{
	IloInt i,j,l;
	IloNum somme = 0;

	//cout << "Debut centralization" << endl;
	//cin.get();

	//cout << "L Centralization h = " << L << endl;
	//cin.get();


	for(i = 0 ; i < m ; i++)
		for(j = i+1 ; j < m ; j++){

			//cout << "L[ " << i << "] = " << L[i]  << " - " <<  "L[ " << j << "] = " << L[j] << endl;
			//cout << "profit = " << profit[L[i]][L[j]] << endl;
			somme += profit[L[i]][L[j]];
		}

	//cout << "Fin centralization" << endl;
	//cin.get();

	//cout << "Centralization h = " << somme << endl;

	//cout << "Centralization h = " << somme/(m*(m-1)/2) << endl;
	//cin.get();

	return(somme/(m*(m-1)/2));
}
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
IloNum Greedy::Colocation(const IloEnv & env)
{
	IloInt i;
	IloNum somme = 1;
	IntVector ordre(env,m);

	for(i = 0 ; i < m ; i++)
		ordre[i] = i;

	//cout << "Debut Colocation" << endl;
	//cin.get();

	//cout << "L = " << L << endl;
	//cin.get();

	ROPTri(L,ordre,0,m-1,1);

	//cout << "L = " << L << endl;
	//cout << "ordre = " << ordre << endl;
	//cin.get();

	for(i = 1 ; i < m ; i++)
		if( (IloInt) L[ordre[i]] != (IloInt) L[ordre[i-1]])
			somme ++;

	//cout << "Fin colocation" << endl;
	//cin.get();

	return( ( ((m - somme)/m)*100 ) );
}


void Greedy::AffichageSolutions(char* filename)
{
	ofstream fichier(filename, ios::out | ios::app);

	for(int i = 0 ; i < m-1 ; i++)
		fichier << L[i] << " ";

	fichier << L[m-1] << endl;


	fichier.close();

}

