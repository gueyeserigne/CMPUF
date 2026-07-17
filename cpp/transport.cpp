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
#include <iostream>
#include <string>

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


#include "../h/transport.h"
#include "../h/tri.h"
#include <time.h>

#define INFINI 1e+30
#define EPS 1e-10


ILOSTLBEGIN

//////////////////////////////////////////////////////////////////
//	Constructor						//
//////////////////////////////////////////////////////////////////
Transport::Transport(const IloEnv  & env, IloInt n_arg, IloInt m_arg, FloatMatrix p_arg, IloNumArray h_arg, IloNumArray c_arg):IloModel(env)
{
	//cout << "Transport " << endl;
	//cin.get();
	n = n_arg;
	m = m_arg;
	p = FloatMatrix(env);
	c = IloNumArray(env);
	h = IloNumArray(env);

	h = h_arg;
	p = p_arg;
	c = c_arg;

/*
	cout << "n =" << n << endl;
	cout << "m =" << m << endl;
	cin.get();
	cout << "h =" << h << endl;
	cin.get();

	cout << "profit = " << endl;

	for(int i = 0 ; i < n ; i++)
		cout << p[i] << endl;

	cin.get();

	cout << "C = " << c << endl;

	cin.get();
*/
}



/////////////////////////////////////////////////////////////////// 
// Return optimal value VOPT 
//////////////////////////////////////////////////////////////////
IloNum Transport::Opt()
{
	return(VOPT);
}
//////////////////////////////////////////////////////////////////
// But : Creation et Dimensionnement des variables x	        //
//////////////////////////////////////////////////////////////////
void Transport::Init_x(const IloEnv & env)
{
	int i,l;
	char nom[100];
	char tmp[100];

	/*cout << "debut Init x" << endl;
	cin.get();

	cout << "n = " << n;
	cout << "m = " << m;
	cin.get();
*/

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
/*
	cout << "Fin Init x" << endl;
	cin.get();
*/
}
//////////////////////////////////////////////////////////////////
// But : Creation et Dimensionnement des variables z	        //
//////////////////////////////////////////////////////////////////
void Transport::Init_z(const IloEnv  & env)
{
	int l,j;
	char nom[100];
	char tmp[100];

/*	cout << "init_z beginning ... " << endl;
	cin.get();
*/
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
/*
	cout << "init_z end ... " << endl;
	cin.get();
*/
}
//////////////////////////////////////////////////////////////////
// But : Buil and add constraints                               //
//////////////////////////////////////////////////////////////////
void Transport::Init_contraintes(const IloEnv  & env)
{
	IloInt i,j,l;

/*	cout << "Init contraintes " << endl;
	cin.get();
*/
	
	for(i = 0 ; i < m ; i++){
		IloExpr expr(env);

		for(int l = 0 ; l < n ; l++)
			expr += xVar[i][l];

//		cout << "expr = " << expr << endl;

		add(IloRange(env,1,expr,1)); // Each facility i is located in exactly one location l
	}

/*	cout << "contraintes x... " << endl;
	cin.get();
*/
        // Capacity constraint
	for(l = 0 ; l < n ; l++){
		IloExpr expr(env);

		for(j = 0 ; j < n ; j++)
			expr +=zVar[l][j];

		for(i = 0 ; i < m ; i++)
			expr -= c[i]* xVar[i][l];

		add(IloRange(env,expr,0));
	}

/*	cout << "contraintes de demande... " << endl;
	cin.get();
*/	
        // Demand constraints 
	for(j = 0 ; j < n ; j++)
	{
		IloExpr expr(env);

		for(l = 0 ; l < n ; l++)
			expr += zVar[l][j];

		add(IloRange(env,expr,h[j]));
	}

/*	cout << "contraintes end ... " << endl;
	cin.get();
*/
}

//////////////////////////////////////////////////////////////////
// Objective function   					//
//////////////////////////////////////////////////////////////////
void Transport::Init_obj(const IloEnv  & env)
{
	IloInt i,j,l;
	IloExpr expr(env);
	IloNum max = 0;
/*
	cout << "init_obj beginning ... " << endl;
	cin.get();	
*/
	for(l = 0 ; l < n ; l++)
		for(j = 0 ; j < n ; j++)
			if(p[l][j] > max)
				max = p[l][j];


	for(l = 0 ; l < n ; l++)
		for(j = 0 ; j < n ; j++)
			expr += zVar[l][j]*(max - p[l][j]);

	obj = IloObjective(env,expr,IloObjective::Maximize);  

	add(obj);
/*
	cout << "init_obj end ... " << endl;
	cin.get();
*/
}
//////////////////////////////////////////////////////////////////
// Solve the problem, and get optimal solution and value
//////////////////////////////////////////////////////////////////
void Transport::Resolution(const IloEnv & env, IloCplex cplex)
{
	IloInt i,j,l;

	//cout << "Begin resolution du problème de transport" << endl;
	for(int i = 0 ; i < m ; i++)
		for(int l = 0 ; l < n ; l++)
			add(IloConversion(env,xVar[i][l],ILOINT));


	cplex.solve();

	//cout << "status = " << cplex.getCplexStatus();
	//cin.get();

	if((cplex.getCplexStatus() == IloCplex::Optimal) || 
		(cplex.getCplexStatus() == IloCplex::OptimalTol)|| 
			(cplex.getCplexStatus() == IloCplex::SolLim) )
	{
		VOPT = cplex.getObjValue();

		// lecture solution x

		for(i = 0 ; i < m ; i ++) 
			for(l = 0 ; l < n; l++)  
				if(cplex.isExtracted(xVar[i][l]))
						x[i][l] =  cplex.getValue(xVar[i][l]);
					else
						x[i][l] =  0;


		// lecture solution z

		for(l = 0 ; l < n ; l ++) 
			for(j = 0 ; j < n; j ++)  
				if(cplex.isExtracted(zVar[l][j]))
					z[l][j] =  cplex.getValue(zVar[l][j]);
				else
					z[l][j] =  0;
	}
/*
	cout << "z = " << endl;
	for(l = 0 ; l < n ; l ++) 
		cout << z[l] << endl;
	cin.get();
*/
	//cout << "End resolution p-median" << endl;
}
//////////////////////////////////////////////////////////////////
// Initialized an array L according to the values in the variables xVar
// L[i] is the location of i
//////////////////////////////////////////////////////////////////
IloNumArray Transport::Init_L(const IloEnv  & env)
{
	IloInt i,l;
	IloNum eps = 1e-2;

	//L = IloNumArray(env,m);
	IloNumArray L;

	L = IloNumArray(env,m);

	for(i = 0 ; i < m ; i++)
	{
		L[i] = 0;

		for(l = 0 ; l < n ; l++){
			if(x[i][l] > 1-eps)
			//L[i] += l*x[i][l];
				L[i] = l;
		}
	}

	return(L);

	//cout << "L = " << L << endl;
}
//////////////////////////////////////////////////////////////////
// ???
//////////////////////////////////////////////////////////////////
IloNumArray Transport::Init_Supply(const IloEnv  & env, IloNumArray L, IloNum prob)
{
	IloInt i,j,l;
	IloNum eps = 1e-2;
	IloNum max = 0;


	//L = IloNumArray(env,m);
	IloNumArray S;

	for(l = 0 ; l < n ; l++)
		for(j = 0 ; j < n ; j++)
			if(p[l][j] > max)
				max = p[l][j];


	S = IloNumArray(env,m);

	for(i = 0 ; i < m ; i++)
	{
		S[i] = 0;
		l = L[i];
		for(j = 0 ; j < n; j ++)  
			S[i] += z[l][j]*(max-p[l][j]);

		S[i] = (int) (S[i]*prob);

		
	}

	return(S);

	//cout << "L = " << L << endl;
}

