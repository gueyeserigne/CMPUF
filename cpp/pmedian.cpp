/**
 * Author: Serigne Gueye
 * Date: 17 July 2026
 */

/****************************************************************/
/* Author : Serigne Gueye						*/
/****************************************************************/

//#include <ilcplex/ilocplex.h>
//#include <ilsolver/ilosolverint.h>

#include <ilcplex/ilocplex.h>
//#include <ilcp/cp.h>
#include <iostream>
#include <string>
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_matrix.h>

typedef IloArray<IloEnv> IloEnvArray;
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


#include "../h/pmedian.h"
#include "../h/tri.h"
#include <time.h>

#define INFINI 1e+30
#define EPS 1e-10


ILOSTLBEGIN

//////////////////////////////////////////////////////////////////
// Constructor  						//
//////////////////////////////////////////////////////////////////

Pmedian::Pmedian(const IloEnv  env, char* fichier,std::string type, IloInt m_arg)
{
	IloInt i,j;
	IloNumArray ptmp;

	//cout << "n = " << n << endl;
	//cout << "m = " << m << endl;

	// Lecture des donnees du probleme
	
	ifstream file(fichier);

	//cout << "fichier =" << fichier << endl;
	//cin.get();
	
  	if(! file)
    		cerr << "Erreur a l'ouverture  du fichier " << fichier << endl;
  	else
  	{
    		file >> n;
    		h = IloNumArray(env,n);
    		if(type.compare("geometric") == 0){
    		    int mtmp;
    		    IloNumArray ptmp;
    		    file >> mtmp;
    		    ptmp = IloNumArray(env,mtmp);
    		    file >> ptmp;
    		    file >> h;
		}
		else{		
        		for(i = 0 ; i < n ; i++)
	        		file >> h[i];
	      	}

		profit = FloatMatrix(env,n);
		for(i = 0 ; i < n ; i++)
			profit[i]  = IloNumArray(env,n);	
		
 		if(type.compare("geometric") == 0)
		    file >> profit;
		else{
        	    for(i = 0 ; i < n ; i++){
			    for(j = 0 ; j < n ; j++)
				    file >> profit[i][j];
		    }
		}

		file.close();
	}

	m = m_arg;

	model = IloModel(env);
	cplex = IloCplex(model);
/*
	cout << "h = " << h << endl;
	cout << "profit = " << endl;
	for(i = 0 ; i < n ; i ++) 
		cout << profit[i] << endl;
	cin.get();
*/
	//cout << "Demande totale = " << IloSum(h) << endl;
	//cin.get();
}
//////////////////////////////////////////////////////////////////
// Return VOPT
//////////////////////////////////////////////////////////////////
IloNum Pmedian::Opt()
{
	return(VOPT);
}
//////////////////////////////////////////////////////////////////
// Build X							//
//////////////////////////////////////////////////////////////////
void Pmedian::Init_X(const IloEnv  env)
{
	int l;
	char nom[100];
	char tmp[100];

	//cout << "init_X beginning ... " << endl;
	//cin.get();

	// Allocation x

	X = IloNumArray(env,n);

	// Allocation xVar

	XVar = IloNumVarArray(env);

	for(l = 0 ; l < n ; l++)
	{
		strcpy(nom,"X");
		sprintf(tmp,"%d",l);
		strcat(nom,tmp);
		XVar.add(IloNumVar(env,0,m,nom));
	}

	//cout << "init_X end ... " << endl;
	//cin.get();
}
//////////////////////////////////////////////////////////////////
// Build z								//
//////////////////////////////////////////////////////////////////
void Pmedian::Init_z(const IloEnv  env)
{
	int l,j;
	char nom[100];
	char tmp[100];

	//cout << "init_z beginning ... " << endl;
	//cin.get();

	// Allocation z

	z = FloatMatrix(env,n);

	for(l = 0 ; l < n ; l++)
		z[l] = IloNumArray(env,n);

	// Allocation zVar

	zVar = FloatVarMatrix(env,n);

	for(l = 0 ; l < n ; l++)
	{
		zVar[l] = IloNumVarArray(env);
		for(j = 0 ; j < n ; j++)
		{
			strcpy(nom,"z");
			sprintf(tmp,"%d.%d",l,j);
			strcat(nom,tmp);
			zVar[l].add(IloNumVar(env,0,IloInfinity,nom));
		}
	}

	//cout << "init_z end ... " << endl;
	//cin.get();
}
//////////////////////////////////////////////////////////////////
// Add m-median constraints					//
//////////////////////////////////////////////////////////////////
void Pmedian::Init_contraintes(const IloEnv  env)
{
	IloInt j,l;

	//cout << "contraintes begin ... " << endl;
	//cin.get();
	
	IloExpr expr(env);

	for(l = 0 ; l < n ; l++)
		expr += XVar[l];

	model.add(IloRange(env,m,expr,m));

	//cout << "contraintes de capacités... " << endl;
	//cin.get();
	
	for(l = 0 ; l < n ; l++)
		for(j = 0 ; j < n ; j++){
			IloExpr expr(env);
			expr = zVar[l][j] - h[j]*XVar[l];
			model.add(IloRange(env,expr,0)); // Add the following : il there is no facility in l then l cannot serve no client
		}

	//cout << "contraintes de demande... " << endl;
	//cin.get();
	
	for(j = 0 ; j < n ; j++){
		IloExpr expr(env);
		for(l = 0 ; l < n ; l++)
			expr += zVar[l][j];
		model.add(IloRange(env,h[j],expr,h[j])); // Each client j must receive h[j] unity of goods for facilities located in the nodes l
	}

	//cout << "contraintes end ... " << endl
	//cin.get();

}
//////////////////////////////////////////////////////////////////
// Build objective fucntion       				//
//////////////////////////////////////////////////////////////////
void Pmedian::Init_obj(const IloEnv  env)
{
	IloInt i,j,l;
	IloExpr expr(env);
	IloNum max = 0;

	//cout << "init_obj beginning ... " << endl;
	//cin.get();	

	for(l = 0 ; l < n ; l++)
		for(j = 0 ; j < n ; j++)
			if(profit[l][j] > max)
				max = profit[l][j];

	for(l = 0 ; l < n ; l++)
		for(j = 0 ; j < n ; j++)
			expr += zVar[l][j]*(max - profit[l][j]);

	obj = IloObjective(env,expr,IloObjective::Maximize);  

	model.add(obj);

	//cout << "init_obj end ... " << endl;
	//cin.get();
}
/////////////////////////////////////////////////////////////////
// Convert xVar variables to binary
//////////////////////////////////////////////////////////////////
void Pmedian::Conversion(const IloEnv env)
{
	IloInt i;

	for(i = 0 ; i < n ; i++)
		model.add(IloConversion(env,XVar[i],ILOINT));
}
//////////////////////////////////////////////////////////////////
// Solve the problem, get optimal solution and value
//////////////////////////////////////////////////////////////////
void Pmedian::Resolution()
{
	IloInt j,l;

	//cout << "Begin resolution p-median" << endl;

	cplex.setParam(IloCplex::MIPDisplay,0);
	//cplex.exportModel("cmpuf_pmedian.lp");

	cplex.solve();

	//cout << "status = " << cplex.getCplexStatus();
	//cin.get();

	if((cplex.getCplexStatus() == IloCplex::Optimal) || 
		(cplex.getCplexStatus() == IloCplex::OptimalTol)|| 
			(cplex.getCplexStatus() == IloCplex::SolLim) ){
			
		VOPT = cplex.getObjValue();
		// lecture solution x
		for(l = 0 ; l < n ; l ++) 
			if(cplex.isExtracted(XVar[l]))
						X[l] =  cplex.getValue(XVar[l]);
					else
						X[l] =  0;
		// lecture solution z
		for(l = 0 ; l < n ; l ++) 
			for(j = 0 ; j < n; j ++)  
				if(cplex.isExtracted(zVar[l][j]))
					z[l][j] =  cplex.getValue(zVar[l][j]);
				else
					z[l][j] =  0;
	}
/*
	cout << "x = " << X << endl;
	cout << "z = " << endl;
	for(l = 0 ; l < n ; l ++) 
		cout << z[l] << endl;
	cin.get();
*/
	//cout << "End resolution p-median" << endl;
}
//////////////////////////////////////////////////////////////////
// Convert location vector X to the matrix form xpmedian
//////////////////////////////////////////////////////////////////
FloatMatrix Pmedian::XtoXpmedianSolution(const IloEnv env)
{
	int i,l,k,d,r;
	double tmp;

	//cout << "xtopmedian 1" << endl;
	//cin.get();

	xpmedian = FloatMatrix(env,m);

	for(i = 0 ; i < m ; i++)
		xpmedian[i] = IloNumArray(env,n);

	for(i = 0 ; i < m ; i++)
		for(l = 0 ; l < n ; l ++) 
			xpmedian[i][l] = 0;

	//cout << "xtopmedian 2" << endl;
	//cin.get();

	d = 0;

	for(l = 0 ; l < n ; l ++) {
		tmp = (X[l]);
		r = 0;
		//cout << "X[" << l << "] = " << X[l] << endl;
		//cout << "tmp = " << tmp << endl;
		for(k = 0 ; k < (tmp - 1e-3) ; k ++) {
			//cout << "k = " << k << endl;
			//cout << "k+d = " << k+d << endl;
			xpmedian[k+d][l] = 1;
			r++;
		}
		
		/*
		for(i = 0 ; i < m ; i ++) 
			cout << xpmedian[i] << endl;

		cin.get();*/
		
		d += r;
	}
/*
	cout << "xpmedian = " << endl;
	for(i = 0 ; i < m ; i ++) 
		cout << xpmedian[i] << endl;
*/
	return(xpmedian);
}
//////////////////////////////////////////////////////////////////
// Return an array of the sum of the demands assigned to each facility 
// in the optimal m-median solution
//////////////////////////////////////////////////////////////////
IloNumArray Pmedian::CapaciteMedian(const IloEnv env)
{
	IloInt i,j,l;
	IloBool found;

	IloNumArray Cm  = IloNumArray(env,m);

	for(i = 0 ; i < m ; i++)
	{
		Cm[i] = 0;
		l = 0;
		found = IloFalse;

		while((l < n) && (! found)) 
		{
			if(xpmedian[i][l] > 0)
			{
				found = IloTrue;
				

				for(j = 0 ; j < n ; j++)
					Cm[i] += z[l][j];
			}

			l++;
		}
	}

	//cout << "Cm = " << Cm << endl;
	//cin.get();

	//cout << "Capacité totale = " << IloSum(Cm) << endl;
	//cin.get();


	return(Cm);
}

