/**
 * Author: Serigne Gueye
 * Date: 17 July 2026
 */

/****************************************************************/
/* Serigne Gueye						*/
/* Projet IEF FLOUE MIT-ZLC					*/
/****************************************************************/

#include <ilcplex/ilocplex.h>
//#include <ilcp/cp.h>
#include <bitset>


typedef IloArray<IloNumArray> FloatMatrix;
typedef IloArray<FloatMatrix> FloatMatrix3d;
typedef IloArray<FloatMatrix3d> FloatMatrixMatrixMatrix;
typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<IloBoolArray> BoolMatrix;
typedef IloArray<IloNumVarArray> FloatVarMatrix;
typedef IloArray<FloatVarMatrix> FloatVarMatrixMatrix;
typedef IloArray<FloatVarMatrixMatrix> FloatVarMatrixMatrixMatrix;
typedef IloArray <FloatMatrix> FloatMatrixMatrix;
typedef IloArray<IloInt> IntVector;
typedef IloArray<IntVector> IntMatrix;
typedef IloArray <IntMatrix> IntMatrixMatrix;
typedef IloArray<IloRangeArray> IloRangeArray1;
typedef IloArray<IloRangeArray1> IloRangeArray2;
typedef IloArray<IloRangeArray2> IloRangeArray3;

#include "../h/cmpuf.h"
#include "../h/tri.h"
#include "../h/transport.h"
#include <time.h>

ILOSTLBEGIN

//////////////////////////////////////////////////////////////////
// But : Creation et Dimensionnement des variables z	        //
// Variables : z_{ljk}						//
//////////////////////////////////////////////////////////////////
void CMpuf::Init_z(const IloEnv & env)
{
	int j,l,k;
	char nom[100];
	char tmp[100];

	//cout << "debut Init z" << endl;
	//cin.get();

	zVar = FloatVarMatrixMatrix(env,n);
	z = FloatMatrixMatrix(env,n);

	for(l = 0 ; l < n ; l++)
	{
		zVar[l] = FloatVarMatrix(env,n);
		z[l] = FloatMatrix(env,n);

		for(j = 0 ; j < n ; j ++)
		{
			zVar[l][j] = IloNumVarArray(env);
			z[l][j] = IloNumArray(env,K);

			for(k = 0 ; k < K ; k++)
			{
				strcpy(nom,"z");
				sprintf(tmp,"%d.%d.%d",l,j,k);
				strcat(nom,tmp);
				zVar[l][j].add(IloNumVar(env,0,IloInfinity,nom));
			}
		}  
	}

	//cout << "Fin Init z" << endl;
	//cin.get();
}
//////////////////////////////////////////////////////////////////
// Fix the variables zVar to values in z_prec
//////////////////////////////////////////////////////////////////
void CMpuf::Fixation_z(FloatMatrixMatrix z_prec, IloInt K_prec)
{
	IloInt l,j,k;

	// Pour tous les scénarios
	for(k = 0 ; k < K_prec ; k++)
		for(l = 0 ; l < n ; l++)
			for(j = 0 ; j < n ; j++)
			{
				// La variable zVar[l][j][k] est fixée à la valeur z_prec[l][j][k]
				zVar[l][j][k].setUB(z_prec[l][j][k]);
				zVar[l][j][k].setLB(z_prec[l][j][k]);
			}

}
//////////////////////////////////////////////////////////////////
// Check extracted zVar variables
//////////////////////////////////////////////////////////////////
void CMpuf::Solution_z(IloInt K_prec, IloCplex cplex)
{
	IloInt l,j,k;

	//cout << "Begin solution_z" << endl;

	for(k = 0 ; k < K ; k ++) 
		for(l = 0 ; l < n ; l ++) 
			for(j = 0 ; j < n; j ++)  
				if(!cplex.isExtracted(zVar[l][j][k]))
					cout << "Variable non extraite !" << endl;

	//cout << "End solution_z" << endl;
}
//////////////////////////////////////////////////////////////////
// Creation of x variables
//////////////////////////////////////////////////////////////////
void CMpuf::Init_x(const IloEnv & env)
{
	int i,l;
	char nom[100];
	char tmp[100];

	//cout << "debut Init x" << endl;
	//cin.get();

	//cout << "n = " << n;
	//cout << "m = " << m;
	//cin.get();


	xVar = FloatVarMatrix(env,m);
	x = FloatMatrix(env,m);


	for(i = 0 ; i < m ; i++)
	{
		xVar[i] = IloNumVarArray(env);
		x[i] = IloNumArray(env,n);

		for(l = 0 ; l < n ; l++)
		{
			strcpy(nom,"x");
			sprintf(tmp,"%d.%d",i,l);
			strcat(nom,tmp);
			xVar[i].add(IloNumVar(env,0,1,nom));
			//xVar[i].add(IloNumVar(env,nom));
		}
	}  

	//cout << "Fin Init x" << endl;
	//cin.get();
}
//////////////////////////////////////////////////////////////////
// Variables annotation for Benders decomposition
//////////////////////////////////////////////////////////////////
void CMpuf::Annotations(IloCplex cplex)
{
	IloCplex::LongAnnotation decomp = cplex.newLongAnnotation(IloCplex::BendersAnnotation,CPX_BENDERS_MASTERVALUE + 1);

	for (IloInt i = 0; i < m; i++){
		for(IloInt l = 0 ; l < n ; l++)
			cplex.setAnnotation(decomp, xVar[i][l], CPX_BENDERS_MASTERVALUE);
	}

	for(IloInt k = 0 ; k < K ; k ++) 
		for(IloInt l = 0 ; l < n ; l ++) 
			for(IloInt j = 0 ; j < n; j ++)
				cplex.setAnnotation(decomp, zVar[l][j][k], CPX_BENDERS_MASTERVALUE+k+1);
}
//////////////////////////////////////////////////////////////////
// Fix xVar to location in the array fixed
//////////////////////////////////////////////////////////////////
void CMpuf::Fixation_x(const IloEnv & env, IloNumArray Fixed)
{
	IloInt i;

	IloInt k = Fixed.getSize();

	for(i = 0 ; i < k ; i++)
	{
		xVar[i][Fixed[i]].setLB(1.0);
		xVar[i][Fixed[i]].setUB(1.0);
	}

	//cout << "Fin Fixation_x" << endl;
	//cin.get();
}
//////////////////////////////////////////////////////////////////
// Fix xVar[i][Fixed[i]].setLB(1.0) to 1
//////////////////////////////////////////////////////////////////
void CMpuf::Fixation_x(IloNumArray Fixed, int i)
{
	IloInt k = Fixed.getSize();

	if(i < k){
		xVar[i][Fixed[i]].setLB(1.0);
		xVar[i][Fixed[i]].setUB(1.0);
	}
}
//////////////////////////////////////////////////////////////////
// Fix xVar[i][l] to 0
//////////////////////////////////////////////////////////////////
void CMpuf::Zero_x(int i, int l)
{
	xVar[i][l].setLB(0.0);
	xVar[i][l].setUB(0.0);
}
//////////////////////////////////////////////////////////////////
// Fix the line xVar[i] to 0
//////////////////////////////////////////////////////////////////
void CMpuf::Zero_x(int i)
{
	for (int l = 0; l  < n; l++)
		Zero_x(i,l);
}
//////////////////////////////////////////////////////////////////
// Fix xVar to 0
//////////////////////////////////////////////////////////////////
void CMpuf::Zero_x()
{
	for (int i = 0; i  < m; i++){
		for (int l = 0; l  < n; l++)
			Zero_x(i,l);
	}
}
//////////////////////////////////////////////////////////////////
// Fix xVar[i][l] to 1
//////////////////////////////////////////////////////////////////
void CMpuf::Un_x(int i, int l)
{
	xVar[i][l].setLB(1.0);
	xVar[i][l].setUB(1.0);
}
//////////////////////////////////////////////////////////////////
// Defix (free) xVar[i][l]
//////////////////////////////////////////////////////////////////
void CMpuf::Zero_Un_x(int i, int l)
{
	xVar[i][l].setLB(0.0);
	xVar[i][l].setUB(1.0);
}
//////////////////////////////////////////////////////////////////
// Defix xVar[i]
//////////////////////////////////////////////////////////////////
void CMpuf::Zero_Un_x(int i)
{
	for (int l = 0; l  < n; l++)
		Zero_Un_x(i,l);
}
//////////////////////////////////////////////////////////////////
// Defix xVar
//////////////////////////////////////////////////////////////////
void CMpuf::Zero_Un_x()
{
	for (int i = 0; i  < m; i++){
		for (int l = 0; l  < n; l++)
			Zero_Un_x(i,l);
	}
}
//////////////////////////////////////////////////////////////////
// Fix xVar to the solution in sol
//////////////////////////////////////////////////////////////////
void CMpuf::Solution(FloatMatrix sol)
{
	IloNum eps = 1e-10;

	for(int i = 0 ; i < m ; i++)
	{

		for(int l = 0 ; l < n ; l++){
			if(sol[i][l] > 1-eps)
				Un_x(i,l);
			else
				Zero_x(i,l);

		}
	}

}
//////////////////////////////////////////////////////////////////
// Defix xVar[i][Fixed[i]]
//////////////////////////////////////////////////////////////////
void CMpuf::Defixation_x(IloNumArray Fixed, int i)
{
	IloInt k = Fixed.getSize();

	if(i < k){
		xVar[i][Fixed[i]].setLB(0.0);
		xVar[i][Fixed[i]].setUB(1.0);
	}
}
//////////////////////////////////////////////////////////////////
// Defix xVar to the values previously fixed in Fixed
//////////////////////////////////////////////////////////////////
void CMpuf::Defixation_x(IloNumArray Fixed)
{
	IloInt k = Fixed.getSize();


	for(int i = 0; i < k; i++)
		Defixation_x(Fixed,i);
}
//////////////////////////////////////////////////////////////////////////
// But : Creation et Dimensionnement des variables xpemdian	        //
// Fixation de la solution du problème pmedian
//////////////////////////////////////////////////////////////////////////
void CMpuf::Init_xpmedian(const IloEnv & env, FloatMatrix & xpmedian)
{
	int i,l;

	//cout << "debut Init xpmedian" << endl;
	//cin.get();

	//cout << "n = " << n;
	//cout << "m = " << m;
	//cin.get();

	xpmedian = FloatMatrix(env,m);

	for(i = 0 ; i < m ; i++)
		xpmedian[i] = IloNumArray(env,n);


	for(i = 0 ; i < m ; i ++) 
		for(l = 0 ; l < n; l++)  
			xpmedian[i][l] =  x[i][l];

	//cout << "Fin Init xpmedian" << endl;
	//cin.get();
}
//////////////////////////////////////////////////////////////////
// Add the location constraint saying that i must be located in exactly one of the location 
// l = 0,1,...n-1
//////////////////////////////////////////////////////////////////
void CMpuf::Init_contraintes_x(const IloEnv & env, int i)
{
	IloExpr expr(env);

	for(int l = 0 ; l < n ; l++)
		expr += xVar[i][l];

	add(IloRange(env,1,expr,1));
}
//////////////////////////////////////////////////////////////////
// Add for each facility i a location constraint
//////////////////////////////////////////////////////////////////
void CMpuf::Init_contraintes_x(const IloEnv & env)
{
	for(int i = 0 ; i < m ; i++)
		Init_contraintes_x(env,i);
}
//////////////////////////////////////////////////////////////////
// Add the constraint saying that for each location l at max one facility i may be allocated in l 
// l = 0,1,...n-1
//////////////////////////////////////////////////////////////////
void CMpuf::Init_contraintes_x_x(const IloEnv & env)
{

	for(int l = 0 ; l < n ; l++){
		IloExpr expr(env);
		for(int i = 0 ; i < m ; i++)
			expr += xVar[i][l];


		add(IloRange(env,expr,1));
	}

}
//////////////////////////////////////////////////////////////////
// Add the constraints on z and coupling z and x
//////////////////////////////////////////////////////////////////
void CMpuf::Init_contraintes(const IloEnv & env)
{
	int i,j,k,l;
	IloNum eps = 1e-10;
	char nom[100];
	char tmp[100];
        ////////////////////////////////////////////
	//cout << "Debut Init contraintes" << endl;
	//cin.get();
        ////////////////////////////////////////////
	int cpt = 0;

	for(k = 0 ; k < K ; k++)
	{
		for(l = 0 ; l < n ; l++)
		{
			IloExpr expr(env);

			for(j = 0 ; j < n ; j++)
				expr +=(zVar[l][j][k]);

			for(i = 0 ; i < m ; i++)
				expr -= C[k][i]* xVar[i][l];

			add(IloRange(env,expr,0));
			cpt++;
		}

		for(j = 0 ; j < n ; j++)
		{
			IloExpr expr(env);

			for(l = 0 ; l < n ; l++)
				expr +=(zVar[l][j][k]);

			strcpy(nom,"Demand.");
			sprintf(tmp,"%d.%d",j,k);
			strcat(nom,tmp);

			dcon.add(IloRange(env,expr,w[j],nom));
			d.add(0);

			//add(IloRange(env,expr,w[j]));
			add(dcon[dcon.getSize()-1]);
			cpt++;
			//add(IloRange(env,w[j],expr,w[j]));
		}		
	}

	//cout << "Fin Init contraintes = " << cpt << endl;
	//cin.get();
}
////////////////////////////////////////////
// Creates the objective function
////////////////////////////////////////////
void CMpuf::Init_obj(const IloEnv & env)
{
	IloInt i,j,k,l;
	IloObjective obj;
	IloExpr expr(env);
	IloNum max = 0;

	for(l = 0 ; l < n ; l++)
		for(j = 0 ; j < n ; j++)
			if(profit[l][j] > max)
				max = profit[l][j];
      
	for(k = 0 ; k < K ; k++){
		for(l = 0 ; l < n ; l++)
			for(j = 0 ; j < n ; j++)
				expr +=  probaScenarios[k] * (max -profit[l][j]) * zVar[l][j][k];
	}

	//cout << "obj = " << expr << endl;
	//cin.get();

	obj = IloObjective(env,expr,IloObjective::Maximize);  

	add(obj);
}
//////////////////////////////////////////////////////////////////
// Computation of the transportation cost associated to the solution store 
// in z
//////////////////////////////////////////////////////////////////
IloNum CMpuf::CoutTransport(IloInt k)
{
	IloInt i,j,l;
	IloNum max = 0;
	IloNum S = 0;

	for(l = 0 ; l < n ; l++)
		for(j = 0 ; j < n ; j++)
			if(profit[l][j] > max)
				max = profit[l][j];      

	for(l = 0 ; l < n ; l++)
		for(j = 0 ; j < n ; j++){
			S +=  (max -profit[l][j]) * z[l][j][k];
			//cout << "S = " << (int) S << endl;
			//cin.get();
		}

	return(S);
}
//////////////////////////////////////////////////////////////////
// Print transportation probability and transportation cost for each scenario k
//////////////////////////////////////////////////////////////////
void CMpuf::AffichageCoutTransport()
{
	IloInt k;

	cout << "Scenario \t" << "Proba \t" << "Cout Transport" << endl;
	for(k = 0 ; k < K ; k++)
		cout << k << " \t" << probaScenarios[k] << " \t" << CoutTransport(k) << endl;
}
//////////////////////////////////////////////////////////////////
// Save VOPT and parameters alpha, beta,gamma in filename
//////////////////////////////////////////////////////////////////
void CMpuf::AffichageCoutTransport(const IloEnv & env, IloNum cent, char* filename)
{
	ofstream fichier(filename, ios::out | ios::app);
	
      
	//for(int k = 0 ; k < K ; k++)
	//	fichier << (int) CoutTransport(k) << " ";

	fichier << fixed; // forces output streams to write floating-point numbers in fixed-point notation instead of scientific notation

	fichier << alpha << "\t" << beta << "\t" << gamma << "\t" << VOPT << "\t" << cent << endl;

	//fichier << cent << endl;

	fichier.close();
}
//////////////////////////////////////////////////////////////////
// convert variables xVar to be binaries
//////////////////////////////////////////////////////////////////
void CMpuf::Conversion(const IloEnv & env)
{
	IloInt i,l;


	for(i = 0 ; i < m ; i++)
		for(l = 0 ; l < n ; l++)
			add(IloConversion(env,xVar[i][l],ILOINT));
}
//////////////////////////////////////////////////////////////////
// solve the model
// Retrieve : optimal value (VOPT), optimal solution (x,z), 
// slack variable values of the constraints in dcon
//////////////////////////////////////////////////////////////////
IloNum CMpuf::Resolution(const IloEnv & env, IloCplex cplex)
{
	IloInt i,j,k,l;
	IloNum res;

	//cout << "Resolution " ;
	//cin.get();
        //cplexk.setParam(IloCplex::TiLim,100);


	cplex.solve();

	//cout << "Begin Resolution " << endl;
	//cout << "status = " << cplex.getCplexStatus();
	//cin.get();

	if((cplex.getCplexStatus() == IloCplex::Optimal) || (cplex.getCplexStatus() == IloCplex::OptimalTol)|| (cplex.getCplexStatus() == IloCplex::SolLim) 
			|| (cplex.getCplexStatus() == IloCplex::AbortTimeLim)  )
	{
		cplex.getSlacks(d,dcon);
		res = cplex.getObjValue();
		//cout << "Res = " << res << endl;
		//cin.get();

		// lecture solution x

		for(i = 0 ; i < m ; i ++) 
			for(l = 0 ; l < n; l++)  
				if(cplex.isExtracted(xVar[i][l]))
						x[i][l] =  cplex.getValue(xVar[i][l]);
					else
						x[i][l] =  0;
		// lecture solution z

		for(k = 0 ; k < K ; k ++) 
			for(l = 0 ; l < n ; l ++) 
				for(j = 0 ; j < n; j ++)  
					if(cplex.isExtracted(zVar[l][j][k]))
							z[l][j][k] =  cplex.getValue(zVar[l][j][k]);
						else
							z[l][j][k] =  0;
	}

	VOPT = res;
	return(res);
	//cout << "End Resolution " << endl;
}


//////////////////////////////////////////////////////////////////
// Greedy Heuristic to solve the model
//////////////////////////////////////////////////////////////////
IloNumArray CMpuf::GreedyHeuristic(const IloEnv & env, IloCplex cplex)
{
	IloNum tmp;
	IloNumArray L(env);
	
	Conversion(env);
	//Annotations(cplex);
	cplex.setParam(IloCplex::Param::Benders::Strategy,-1);	
	//cplex.setParam(IloCplex::Param::RootAlgorithm,2);
	cplex.setParam(IloCplex::Param::NodeAlgorithm,2);
	//cplex.setParam(IloCplex::Param::Benders::WorkerAlgorithm,2);
	cplex.setParam(IloCplex::MIPDisplay,2);
	//cplex.setParam(IloCplex::Param::Preprocessing::Presolve,0);
	/*cplex.setParam(IloCplex::Param::MIP::Strategy::Branch,1);
	cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq,-1);
	cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicEffort,0);
	cplex.setParam(IloCplex::Param::MIP::Cuts::Implied,-1);*/
	//IloNumArray Lstart(env);


	Zero_x(); // Fix all the variables to tes les variables x sont fixées à 0

	for(int i = 0; i < m; i++){ // For each facility i
		//Lstart.add(i);
		//MIPStart(env,cplex,Lstart);
		//cplex.setParam(IloCplex::Param::Advance,0);
		Zero_Un_x(i); // Defix all variables x[i][l] deviennent libre
		Init_contraintes_x(env,i); // Add the location constraint for the facility i
		VGREEDY = Resolution(env,cplex); // Solve the optimization problem and get the optimal value
		//cout << "VGREDDY = " << VGREEDY << endl;
		//AffichageSolutions();
		L.add(Localisation(i));         // Get the optimal location of i

		Un_x(i,L[i]);                   // Fix x[i][L[i]] to 1
	}

	Zero_Un_x(); // Deix all variables

	return(L);
}

//////////////////////////////////////////////////////////////////
// Greedy heuristic using benders decomposition to solve the mathematical 
// program
//////////////////////////////////////////////////////////////////
IloNumArray CMpuf::MyGreedyHeuristic(const IloEnv & env, IloCplex cplex)
{
	IloNum tmp;
	IloNumArray L(env);
	
	Conversion(env);
	Annotations(cplex);
	cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq,-1);
	//cplex.setParam(IloCplex::Param::Benders::Strategy,-1);	
	//cplex.setParam(IloCplex::Param::RootAlgorithm,2);
	//cplex.setParam(IloCplex::Param::NodeAlgorithm,2);
	//cplex.setParam(IloCplex::Param::Benders::WorkerAlgorithm,3);
	//cplex.setParam(IloCplex::MIPDisplay,2);
	//cplex.setParam(IloCplex::Param::Preprocessing::Presolve,0);
	/*cplex.setParam(IloCplex::Param::MIP::Strategy::Branch,1);
	cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq,-1);
	cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicEffort,0);
	cplex.setParam(IloCplex::Param::MIP::Cuts::Implied,-1);*/
	//IloNumArray Lstart(env);

	for(int i = 0; i < m; i++){ // For each facility 

		IloRange ct;
		IloExpr expr(env);
                // Fix to zeor all variables xVar[j][l] for which j > i 
		for(int j = i+1 ; j < m ; j++)
			for(int l = 0 ; l < n ; l++)
				expr +=xVar[j][l];

		ct = IloRange(env,expr,0); 
		add(ct);

		//Lstart.add(i);
		//MIPStart(env,cplex,Lstart);
		//cplex.setParam(IloCplex::Param::Advance,0);
		Init_contraintes_x(env,i);          // Add a location constraint for i
		VGREEDY = Resolution(env,cplex);    // Solve the optimization problem and get the optimal value
		//cout << "VGREDDY = " << VGREEDY << endl;
		//AffichageSolutions();
		L.add(Localisation(i));           // Add the new optimal location

		Un_x(i,L[i]);                     // Fix x[i][L[i]] to 1

		remove(ct);
		ct.end();
		expr.end();
	}

	Zero_Un_x(); // Toutes les variables x redeviennent libre

	return(L);
}
//////////////////////////////////////////////////////////////////
// Objective : computing the optimal value of the model corresponding to the solution "sol"
// Input :
// - env : ILOG environmental variable
// - cplex : cplex object
// - sol : location solution. sol[i][l] = 1 if the facility i is located in the location l
// Output : 
// - optimal value of the problem where the variables x are fixed as in sol
//////////////////////////////////////////////////////////////////
IloNum CMpuf::Resolution(const IloEnv & env, IloCplex cplex, FloatMatrix sol)
{

	/*cplex.setParam(IloCplex::HeurFreq,-1);
	cplex.setParam(IloCplex::FlowPaths,-1);
	cplex.setParam(IloCplex::FlowCovers,-1);
	cplex.setParam(IloCplex::CutPass,-1);
	cplex.setParam(IloCplex::MIPDisplay,0);*/

	//Fixation_x(env,L);
	//cmpuf.Init_L(env,xpmedian);

	Solution(sol);
	//cplex.setParam(IloCplex::Param::RootAlgorithm,3);
	IloNum res = Resolution(env,cplex);

	Zero_Un_x();

	return(res);
}
//////////////////////////////////////////////////////////////////////////////////////////
// Solution printing
//////////////////////////////////////////////////////////////////////////////////////////
void CMpuf::AffichageSolutions()
{
	IloInt i,j,k,l;

	cout << "VOPT =" << VOPT << endl;
	//cin.get();

	//AffichageSeparateurs(n);

	cout << "x" << endl;

	for(i = 0 ; i < m ; i++)
		cout << x[i] << endl;

	for(i = 0 ; i < m ; i++)
	{
		cout << i << ":" << "\t";

		for(l = 0 ; l < n ; l++)
			if(x[i][l] > 0)
				cout << "(" << i << "," << l << ")" << "\t" ;

		cout << endl;
	}

/*
	cout << "z = " << endl;
	for(k = 0 ; k < K ; k++){
		cout << "Scenario " << k << endl;
		for(l = 0 ; l < n ; l ++){
			cout << "l:" << l << "->";
			for(j = 0 ; j < n ; j ++)  
				cout << " " << z[l][j][k];
			cout << endl;
		}			
	}
	cin.get();*/



	k = d.getSize();

/*
	int nonsatisfait = 0;

	for(i = 0 ; i < k ; i++){
		if(d[i] > 0){
			cout << dcon[i].getName() << " = " << d[i] << " " << dcon[i].getUB() << endl;		
			nonsatisfait++;
		}
	}

	cout << "Demande non satisfaite = " << nonsatisfait << endl;*/

	cin.get();
}
//////////////////////////////////////////////////////////////////////////////////////////
// Save the location solution L in filename
//////////////////////////////////////////////////////////////////////////////////////////
void CMpuf::AffichageSolutions(IloNumArray L,char* filename)
{
	/*ofstream fichier(filename, ios::out | ios::app);

	for(int i = 0 ; i < m ; i++)
	{
		for(int l = 0 ; l < n ; l++)
			if(x[i][l] > 0)
				fichier << i << " " << l ;

		fichier << endl;
	}

	fichier.close();*/	


	ofstream fichier(filename, ios::out | ios::app);

	for(int i = 0 ; i < m-1 ; i++)
		fichier << L[i] << " ";

	fichier << L[m-1] << endl;


	fichier.close();

}
//////////////////////////////////////////////////////////////////////////////////////////
// Save the location solution L in filename
//////////////////////////////////////////////////////////////////////////////////////////
void CMpuf::AffichageSolutionsScenarios(IloNumArray L,char* filename)
{
	/*ofstream fichier(filename, ios::out | ios::app);

	for(int i = 0 ; i < m ; i++)
	{
		for(int l = 0 ; l < n ; l++)
			if(x[i][l] > 0)
				fichier << i << " " << l ;

		fichier << endl;
	}

	fichier.close();*/	

	ofstream fichier(filename, ios::out | ios::app);

	for(int i = 0 ; i < m-1 ; i++)
		fichier << L[i] << " ";

	fichier << L[m-1] << endl;

	fichier.close();
}
//////////////////////////////////////////////////////////////////////////////////////////
// Save the optimal location solution of each scenario in filename
//////////////////////////////////////////////////////////////////////////////////////////
void CMpuf::GenerationSolutionsScenarios(const IloEnv & env, IloNumArray L,char* filename)
{
	IntVector ordre(env,K);

	cout << "Solutions Scenarios" << endl;

	for(int k = 0; k < K; k++){
		cout <<  fixed << setprecision (0) << k << "\t";

		for(int i = 0; i < m; i++)
			cout << SolScenarios[k][i] << "\t";
		cout  <<  SolSupply[k] <<  "\t-" << C[k] << endl;
			
	}

 
	cout << "Opt" << "\t";
	for(int i = 0; i < m; i++)
		cout << L[i] << "\t";

	cin.get();

	V =  FloatMatrix(env,m);
	IloNumArray Vtmp(env,K);
	int cpt;

	for(int i = 0; i < m; i++){
		V[i] = IloNumArray(env);

		for(int k = 0; k < K; k++)
			Vtmp[k] = SolScenarios[k][i];

		//cout << "Vtmp = " << Vtmp << endl;
		//cin.get();

		for(int j = 0 ; j < K ; j++)
			ordre[j] = j;

		ROPTri(Vtmp,ordre,0,K-1,1);

		//cout << "ordre = " << ordre << endl;
		//cin.get();


		V[i].add(Vtmp[ordre[0]]);

		for(int k = 1 ; k < K ; k++)
			if( (IloInt) Vtmp[ordre[k]] != (IloInt) Vtmp[ordre[k-1]])
				V[i].add(Vtmp[ordre[k]]);					
	}



	cout << "V" << "\t";
	for(int i = 0; i < m; i++)
		cout << V[i] << endl;

	cin.get();






}
//////////////////////////////////////////////////////////////////////////////////////////
// Return the location index of the faciltity corresponding to the solution x
//////////////////////////////////////////////////////////////////////////////////////////
IloInt CMpuf::Localisation(IloInt i)
{
	IloInt l;

	if(i < m)
	{
		for(l = 0 ; l < n ; l++)
			if(x[i][l] > 1 - 1e-3)
				return(l);
	}

	cout <<"Probleme dans Localisation" << endl;

	return(-1);
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
FloatMatrixMatrix CMpuf::Z(const IloEnv & env)
{
	IloInt l,j,k;

	//cout << "Begin Z..." << endl;
	//cout << "Begin Z..." << endl;
	FloatMatrixMatrix z_prec(env,n);

	//cout << "End Z 0..." << endl;

	for(l = 0 ; l < n ; l++)
	{
		//cout << "l = " << l << endl;
		z_prec[l] = FloatMatrix(env,n);

		for(j = 0 ; j < n ; j ++)
		{
			//cout << "j = " << j << endl;
			z_prec[l][j] = IloNumArray(env,K);
		}
	}

	//cout << "End Z 1..." << endl;

	for(l = 0 ; l < n ; l++)
		for(j = 0 ; j < n ; j++)
			for(k = 0 ; k < K ; k++)
				z_prec[l][j][k] = z[l][j][k];

	//cout << "End Z 2..." << endl;

	return(z_prec);
}
//////////////////////////////////////////////////////////////////
// Local search algorithm
//////////////////////////////////////////////////////////////////
IloNumArray CMpuf::Localsearch(const IloEnv & env, IloCplex cplex, IloNumArray L, IloNum VL, IloNum rayon)
{
	IloBool stop = false;
	/*Conversion(env);
	cplex.setParam(IloCplex::HeurFreq,-1);
	cplex.setParam(IloCplex::FlowPaths,-1);
	cplex.setParam(IloCplex::FlowCovers,-1);
	cplex.setParam(IloCplex::CutPass,-1);
	cplex.setParam(IloCplex::MIPDisplay,0);*/

	IloInt k = L.getSize();
	IloNumArray res(env,k);
        /////////////////////////////////////
        // Fix res to the value in L 
	for(int i = 0; i < k; i++)
		res[i] = L[i];
        /////////////////////////////////////
	Fixation_x(env,res); // Fix the variables xVar to the location in res
        /////////////////////////////////////
        // Add local constraints
        /////////////////////////////////////
	if(rayon < n)
		AddLocalConstraint(env,L,rayon);
        /////////////////////////////////////
	while(!stop){	
		stop = true;
		IloNum max = VL;
		IloNum tmp;
		//cout << "Res = " << res << endl;
		//cout << "max = " << max << endl;
		//cin.get();
                /////////////////////////////////////////////////////////
		IloInt besti;
		IloInt bestli;
                /////////////////////////////////////////////////////////
                // For each facility i 
                // Defix the location of i
                // Solve the optimization problem
                // Then update max
                // Fix the locatio of at its initial value and iterates
                /////////////////////////////////////////////////////////
		for(int i = 0 ; i < m ; i++){
			Defixation_x(res,i);
			tmp = Resolution(env,cplex);
			//cout << "tmp = " << tmp;
			//cin.get();

			if(tmp > max){
				besti = i;
				bestli = Localisation(i);
				max = tmp;

			}				
			Fixation_x(res,i);
		}
		////////////////////////////////////////////////////////
		// VL is the optimal value of the loca search procedure
		// If max > val then 
                ////////////////////////////////////////////////////////
		if(max > VL){
			VL = max;
			//cout << "besti " << besti << endl;
			//cout << "bestli " << bestli << endl;
			Defixation_x(res,besti);          // Defix the location of the facility besti previously fixed
			res[besti] = bestli;              // Fix the best location of besti to bestli
			Fixation_x(res, besti);           // Fix the variable x[besti][bestli] accordingly
			stop  = false;
		}
                ////////////////////////////////////////////////////////
		//free(numeros);
		//cmpuf.Destructeur(env);
		//env.end();
	}

	Defixation_x(res);
	//cout << "Fin Localsearch" << endl;
	//cin.get();
	VLS = VL;

	if(rayon < n)
		RemoveLocalConstraint();

	return(res);
}
//////////////////////////////////////////////////
// Construction a local search constraint around the location L 
//////////////////////////////////////////////////
void CMpuf::AddLocalConstraint(const IloEnv & env, IloNumArray L, IloInt rayon)
{
	IntVector ordre(env,n);
	//cout << "Local Constraint" << endl;
	//cout << "L = " << L << endl;

	ct = IloRangeArray(env);

	for(int i = 0 ; i < m ; i++)
	{

		for(int j = 0 ; j < n ; j++)
			ordre[j] = j;

		//cout << "Profit[" << L[i] << " = " << profit[L[i]] << endl;
                                                              // L[i] is the location of the facility i
                                                              // profit[L[i]] contain the distance  between L(i] and all other nodes
		ROPTri(profit[L[i]],ordre,0,n-1,1);           // The array ordre is sorted in increasing nrder of the values in profit[L[i]]
		//cout << "Ordre = " << ordre << endl;
		//cin.get();


                //////////////////////////////////////////////////////////////
                // The constraint imposed the facility i must be located in one of the k-th (with k = rayon) closest node from i (included i)
                //////////////////////////////////////////////////////////////
		IloExpr expr(env);
		for(int l = 0 ; l < rayon ; l++)
				expr +=xVar[i][ordre[l]];
                //////////////////////////////////////////////////////////////
                ct.add(IloRange(env,1,expr,1));			
		add(ct[ct.getSize()-1]);
	}
}

//////////////////////////////////////////////////
// Remove all local search constraint
//////////////////////////////////////////////////
void CMpuf::RemoveLocalConstraint()
{
	int s = ct.getSize();
	for(int i = 0 ; i < s ; i++)
		remove(ct[i]);

	ct.end();
}
//////////////////////////////////////////////////////////////////
// Branch-and-bound in which a local search constraint has been added
//////////////////////////////////////////////////////////////////
IloNum CMpuf::LocalBranchAndBound(const IloEnv & env, IloCplex cplex, IloNumArray L, IloNum rayon)
{
	IloBool stop = false;
	IloNum res;

	/*Conversion(env);
	cplex.setParam(IloCplex::HeurFreq,-1);
	cplex.setParam(IloCplex::FlowPaths,-1);
	cplex.setParam(IloCplex::FlowCovers,-1);
	cplex.setParam(IloCplex::CutPass,-1);
	cplex.setParam(IloCplex::MIPDisplay,0);*/

	IloInt k = L.getSize();

	if(rayon < n)
		AddLocalConstraint(env,L,rayon);

	res = Resolution(env,cplex);

	if(rayon < n)
		RemoveLocalConstraint();

	return(res);
}
//////////////////////////////////////////////////////////////////
// Fix the MIP start cplex parameter to the values in L
//////////////////////////////////////////////////////////////////
void CMpuf::MIPStart(const IloEnv  env, IloCplex cplex, IloNumArray L)
{
	//cout << "L = " << L << endl;
	//cin.get();
	Init_x(L);
	int k = L.getSize();
	//cplex.setParam(IloCplex::MIPDisplay,2);
	IloNumVarArray startVar(env);
	IloNumArray startVal(env);
	for (int i = 0; i < k; ++i)
		for (int l = 0; l < n; ++l) {
			startVar.add(xVar[i][l]);
			startVal.add(x[i][l]);
         }
	
	cplex.addMIPStart(startVar, startVal);
	startVal.end();
	startVar.end();
}
//////////////////////////////////////////////////////////////////
// For each scenario k = 0,1,...,K-1. The capacity of each facility is known
// Then a transportation problem (Class Transport) is solved
//////////////////////////////////////////////////////////////////
void CMpuf::ResolutionScenarios(const IloEnv & env)
{
	SolScenarios =  FloatMatrix(env);
	SolSupply =  FloatMatrix(env);
/*
	IloNum K = cmpuf.LireK();
	IloInt n = cmpuf.N();
	IloInt m = cmpuf.M();
	IloNumArray w = cmpuf.W();

	FloatMatrix profit = cmpuf.Profit();
*/
	for(int k = 0; k < K; k++){
		IloEnv  envsc;
		Transport scenario(envsc,n,m,profit,w,C[k]);

		scenario.Init_x(envsc);
		scenario.Init_z(envsc);
		scenario.Init_contraintes(envsc);
		scenario.Init_obj(envsc);
		IloCplex cplex(envsc);
		cplex.extract(scenario);
		cplex.setOut(envsc.getNullStream()); 
		cplex.setParam(IloCplex::MIPDisplay,0);
		scenario.Resolution(envsc,cplex);
		IloNumArray L = scenario.Init_L(env);
		IloNumArray S = scenario.Init_Supply(env,L,probaScenarios[k]);
		SolScenarios.add(L);
		SolSupply.add(S);

		envsc.end();
	}
}



