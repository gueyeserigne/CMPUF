/**
 * Author: Serigne Gueye
 * Date: 17 July 2026
 */

/****************************************************************/
/* Author  : Serigne Gueye						*/
/****************************************************************/

//#include <ilcplex/ilocplex.h>
//#include <ilsolver/ilosolverint.h>

#include <ilcplex/ilocplex.h>
//#include <ilcp/cp.h>
#include <iostream>
#include <string>
#include <cstring>

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


#include "../h/cmpuf.h"
#include "../h/pmedian.h"
#include "../h/tri.h"
#include "../h/transport.h"
#include <ctime>

#define INFINI 1e+30
#define EPS 1e-10

ILOSTLBEGIN


//////////////////////////////////////////////////////////////////
// Constructeur d'un objet de type CMPUF 							//
//////////////////////////////////////////////////////////////////
CMpuf::CMpuf(const IloEnv & env , char* fichier, std::string type, char* param, IloNum alpha_arg, IloNum beta_arg, IloNum gamma_arg, IloInt m_arg):IloModel(env)
{
	IloInt i,j,k,l;
	IloNumArray ptmp;

	//cout << "debut Mpuf" << endl;
	//cin.get();  
        // Read the parameters
	alpha = alpha_arg;
	gamma = gamma_arg;
	beta = beta_arg;
	m = m_arg;
	 
	probaScenarios = IloNumArray(env);
	T = FloatMatrix(env);
	dcon = IloRangeArray(env);
	d = IloNumArray(env);
	///////////////////////////////////////////
	// Read the file containing the some problem parameters : ncap,walpha,wbeta
	///////////////////////////////////////////
	ifstream fileparam(param);
        if(! fileparam)
    		cerr << "Erreur a l'ouverture  du fichier " << param << endl;
  	else
  	{
    		fileparam >> ncap;
    		walpha = IloNumArray(env,ncap-1);
    		wbeta = IloNumArray(env,ncap-1);
    		
    		for(i = 0 ; i < ncap-1 ; i++)
	    		fileparam >> walpha[i];

    		for(i = 0 ; i < ncap-1 ; i++)
	    		fileparam >> wbeta[i];
          }
          K = IloPower(ncap,m); // Nombre total de scenarios d'indisponibilité de capacités
          fileparam.close();
	///////////////////////////////////////////
	// Read the file containing the some problem parameters : n,profit
	///////////////////////////////////////////
	ifstream file(fichier);	
	/*
	cout << "fichier =" << fichier << endl;
	cout << "ncap = " << ncap << endl;
	cout << "walpha = " << walpha << endl;
	cout << "wbeta = " << wbeta << endl;
	cin.get();
	*/
	
  	if(! file)
    		cerr << "Erreur a l'ouverture  du fichier " << fichier << endl;
  	else
  	{
    		file >> n;
    		w = IloNumArray(env,n);
    		if(type.compare("geometric") == 0){
    		    int mtmp;
    		    IloNumArray ptmp;
    		    file >> mtmp;
    		    ptmp = IloNumArray(env,mtmp);
    		    file >> ptmp;
    		    file >> w;
    		}
    		else{		
        		for(i = 0 ; i < n ; i++)
	        		file >> w[i];
	      	}
	      	
		H = IloSum(w);
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

        ////////////////////////////////////////////
        // Compute probability vector
        ////////////////////////////////////////////
	p = IloNumArray(env,ncap); // Vecteur des probabilités de         
        for(i = 0 ; i < ncap-1 ; i++)
          p[i] = walpha[i]*alpha;
        p[i] = (1-alpha);          
        
	//cout << "S = " << S << endl;
	//cout << "m = " << m << endl;
	//cout << "beta = " << beta << endl;
	//cin.get();

	//cout << "Fin Mpuf" << endl;
	//cin.get();  
}



///////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////

void CMpuf::Destructeur(const IloEnv & env)
{
	dcon.end();
	d.end();
	T.end();
	probaScenarios.end();
	p.end();
	walpha.end();
	wbeta.end();
}


//////////////////////////////////////////////////////////////////
// Initialization of cap
//////////////////////////////////////////////////////////////////
void CMpuf::Init_cap(const IloEnv & env, IloNumArray Cm)
{
	IloInt i,j;
        /*
	cout << "Begin init_cap" << endl;
        cout << "ncap  = " << ncap << endl;
        cout << "m = " << m << endl;
        cin.get();
        */
	S = IloNumArray(env,m);
	cap = FloatMatrix(env,m);

	for(i = 0 ; i < m ; i++)
		S[i] = gamma*Cm[i];
      /*
        cout << "S = " << S << endl;
        cin.get();
        */
        for(i = 0 ; i < m ; i++)
		cap[i] = IloNumArray(env,ncap);
        /*
        cout << "wbeta = " << wbeta << endl;
        cin.get();
        */
        for(i = 0 ; i < m ; i++)
	{			
		for(j = 0 ; j < ncap-1 ; j++)
		    cap[i][j] = wbeta[j]*S[i];
		cap[i][j] = S[i];
	}
        /*
	cout << "cap = " << cap << endl;
	cin.get();
	*/
	//cout << "End init_cap" << endl;
}
//////////////////////////////////////////////////////////////////
// Initialization of T
//////////////////////////////////////////////////////////////////
void CMpuf::Init_T(const IloEnv  & env, int i, IloNumArray numeros)
{
	IloInt j,k;
	//cout << "i- = " << i << endl;
	//if(i <= m)
	if(i >= 0){
		for(j = 0; j < ncap ; j++){
			numeros[i] = j;
			Init_T(env,i-1,numeros);
		}		
	}
	else{
		IloNumArray temp(env,m);
		IloNum proba = 1;

		for(j = 0 ; j < m ;j++)
			temp[j] = numeros[j];

		T.add(temp);

		for(j = 0 ; j < m ; j++)
			proba *= p[numeros[j]];

		probaScenarios.add(proba);
	}
}
//////////////////////////////////////////////////////////////////
// Initialization of C
//////////////////////////////////////////////////////////////////
void CMpuf::Init_C(const IloEnv & env)
{
	IloInt i,k;

	C = FloatMatrix(env);

	for(k = 0 ; k < K ; k++){
		IloNumArray temp(env,m);

		for(i = 0 ; i < m ; i++)
			temp[i] = cap[i][T[k][i]] ;

		C.add(temp);
	}
}
//////////////////////////////////////////////////////////////////
// But : Affichage des caracteres de separation		        //
// 	permettant de rendre la visualisation plus conviviale	//
//////////////////////////////////////////////////////////////////
void CMpuf::AffichageSeparateurs(IloInt taille)
{
	IloInt i;

	for(i = 0 ; i < taille ; i++)
		cout << "-";
	cout << endl;
}
//////////////////////////////////////////////////////////////////
// But : Print data                                             //
//////////////////////////////////////////////////////////////////
void CMpuf::AffichageDonnees()
{
	IloInt i,j,k;

	AffichageSeparateurs(10);
	cout << "n =" << n << endl;
	cout << "m =" << m << endl;
	cout << "K =" << K << endl;
	AffichageSeparateurs(10);
	cin.get();
	cout << "p =" << p << endl;
	AffichageSeparateurs(10);
	cin.get();
	cout << "w =" << w << endl;
	AffichageSeparateurs(10);
	cin.get();

	cout << "profit = " << endl;

	for(i = 0 ; i < n ; i++)
		cout << profit[i] << endl;

	AffichageSeparateurs(n);
	cin.get();

	cout << K << " scenarios" << endl;
	cout << "C = " << endl;


	for(k = 0 ; k < K ; k++)
		cout << C[k] << endl;

	cout << "Probas = " << endl;
	cout << probaScenarios << endl;

	cin.get();
}
//////////////////////////////////////////////////////////////////
// Return profit attribute
//////////////////////////////////////////////////////////////////
FloatMatrix CMpuf::Profit()
{
	return(profit);
}
//////////////////////////////////////////////////////////////////
// Return profit facility capacities in the scenario k
//////////////////////////////////////////////////////////////////
IloNumArray CMpuf::Capacites(IloInt k)
{
	return(C[k]);
}
//////////////////////////////////////////////////////////////////
// Return demands w
//////////////////////////////////////////////////////////////////
IloNumArray CMpuf::W()
{
	return(w);
}
//////////////////////////////////////////////////////////////////
// Return demands T
//////////////////////////////////////////////////////////////////
FloatMatrix CMpuf::LireT()
{
	return(T);
}
//////////////////////////////////////////////////////////////////
// Return demands alpha
//////////////////////////////////////////////////////////////////
IloNum CMpuf::Alpha()
{
	return(alpha);
}
//////////////////////////////////////////////////////////////////
// Return demands beta
//////////////////////////////////////////////////////////////////
IloNum CMpuf::Beta()
{
	return(beta);
}
//////////////////////////////////////////////////////////////////
// Return demands gamma
//////////////////////////////////////////////////////////////////
IloNum CMpuf::Gamma()
{
	return(gamma);
}
//////////////////////////////////////////////////////////////////
// Return VOPT
//////////////////////////////////////////////////////////////////
IloNum CMpuf::Opt()
{
	return(VOPT);
}
//////////////////////////////////////////////////////////////////
// Return VGREDDY
//////////////////////////////////////////////////////////////////
IloNum CMpuf::Greedy()
{
	return(VGREEDY);
}
//////////////////////////////////////////////////////////////////
// Return VLS
//////////////////////////////////////////////////////////////////
IloNum CMpuf::LS()
{
	return(VLS);
}
//////////////////////////////////////////////////////////////////
// Return VOPTTRANSPORT
//////////////////////////////////////////////////////////////////
IloNum CMpuf::Pmedian()
{
	return(VOPTTRANSPORT);
}
//////////////////////////////////////////////////////////////////
// Return m
//////////////////////////////////////////////////////////////////
IloInt CMpuf::M()
{
	return(m);
}
//////////////////////////////////////////////////////////////////
// Return demands K
//////////////////////////////////////////////////////////////////
IloInt CMpuf::LireK()
{
	return(K);
}
//////////////////////////////////////////////////////////////////
// Return demands n
//////////////////////////////////////////////////////////////////
IloInt CMpuf::N()
{
	return(n);
}
//////////////////////////////////////////////////////////////////
//  Use x to compute an return an array L
//  L[i] is the location of the facility i
//////////////////////////////////////////////////////////////////
IloNumArray CMpuf::Init_L(const IloEnv & env)
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
				L[i] = l;
		}
	}

	return(L);

	//cout << "L = " << L << endl;
}
//////////////////////////////////////////////////////////////////
//  Use xpmedian to compute an return an array L
//  L[i] is the location of the facility i
//////////////////////////////////////////////////////////////////
IloNumArray CMpuf::Init_L(const IloEnv & env,  FloatMatrix & xpmedian)
{
	IloInt i,l;
	IloNum eps = 1e-2;

	IloNumArray L;

	L = IloNumArray(env,m);

	for(i = 0 ; i < m ; i++)
	{
		L[i] = 0;

		for(l = 0 ; l < n ; l++){
			if(xpmedian[i][l] > 1-eps)
			//L[i] += l*x[i][l];
				L[i] = l;
		}
	}

	return(L);
	//cout << "L = " << L << endl;
}
//////////////////////////////////////////////////////////////////
//  Compute the matrix x corresponding to the locations in L
//////////////////////////////////////////////////////////////////
void CMpuf::Init_x(IloNumArray L)
{
	int k = L.getSize();
	for(int i = 0 ; i < k ; i++)
		for(int l = 0 ; l < n ; l++)
			x[i][l] = 0;

	for(int i = 0 ; i < k ; i++)
		x[i][L[i]] = 1;
}
//////////////////////////////////////////////////////////////////
//  Compute the centralization value of the location L
//////////////////////////////////////////////////////////////////
IloNum CMpuf::Centralization(IloNumArray L,FloatMatrix profit)
{
	IloInt i,j;
	IloNum somme = 0;
	IloNum max = 0;
	
	//double a,b;
	//int c,d;
	//cout << "Debut centralization" << endl;
	//cin.get();
	//cout << "L Centralization e = " << L << endl;
	//cin.get();
	//cout << "profit = " << endl;

	for(int l = 0 ; l < n ; l++)
		for(int j = 0 ; j < n ; j++)
			if(profit[l][j] > max)
				max = profit[l][j];

	for(int i = 0 ; i < m ; i++)
		for(int j = i+1 ; j < m ; j++){
			//a = L[i];
			//b = L[j];
			//cout << round(a) << "-" << (int) b << endl;
			//cout << "L = " << L << endl;
			//cout << "L[ " << i << "] = " << a  << " - " <<  "L[ " << j << "] = " << b << endl;
			//cout << "profit = " << profit[a][b] << endl;
			somme += profit[L[i]][L[j]];
			//somme += (max-profit[L[i]][L[j]]);
		}

	//cout << "Fin centralization" << endl;
	//cin.get();

	//cout << "Centralization e = " << somme << endl;
	//cout << "Centralization e = " << somme/(m*(m-1)/2) << endl;
	//cin.get();

	return(somme/(m*(m-1)/2));
}
//////////////////////////////////////////////////////////////////
//  Compute the colocation value of the location L
//////////////////////////////////////////////////////////////////
IloNum CMpuf::Colocation(const IloEnv & env,IloNumArray L)
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


const char* NomApresDernierSeparateur(const char* chemin)
{
    const char* dernierSeparateur = std::strrchr(chemin, '/');
    return dernierSeparateur == nullptr ? chemin : dernierSeparateur + 1;
}
///////////////////////////////////////////
// Concatene in "nom" dir with argv[..]
///////////////////////////////////////////
void concat(char* nom, char* dir, char** argv)
{
	strcpy(nom,dir);
	strcat(nom,NomApresDernierSeparateur(argv[1]));
        strcat(nom,"_");
	strcat(nom,argv[3]);
        strcat(nom,"_");
	strcat(nom,argv[4]);
        strcat(nom,"_");
	strcat(nom,argv[5]);
        strcat(nom,"_");
	strcat(nom,argv[6]);
        strcat(nom,"_");
	strcat(nom,argv[7]);
}
/****************************************************************/
/*								*/
/*			PROGRAMME PRINCIPAL			*/
/*								*/
/****************************************************************/

int main(int argc, char ** argv)
{
	if( argc == 10 )    
	{
		IloInt m,k;
		char filename[200];   // Filename of the input data file
		char nom[300];        // Name of the result ouputs file
		char solgreedy[300];  // Name of greddy solution ouputs file
		char solopt[300];     // Name of the optimal solution outputs file
		char solcent[300];    //  Name of the ouput files containing centralization and colocation values
		char type[] = "e";    // Optimization type
                //////////////////////////////////////        
		IloNum alpha = 0;
		IloNum beta = 0;
		IloNum gamma = 0;
		IloNum pas_alpha = 0.1;
		IloNum pas_gamma = 0.25;
		//IloNum pas_gamma = 0.1;
		IloNum pas_beta = 0.33;
		IloInt pas_m = 2;
		IloInt m_max;
		//IloNum pas_beta = 0.25;
		double topt = -1;
		double tgreedy = -1;
		double tgreedyls = -1;
		double tpmedian1 = -1;
		double tpmedian2 = -1;
		double gap;
		IloNum rayon;
		IloNum rayonbb;
		IloNum duree;
		IloNumArray bestL;
		IloNum bestV;
	        
		char * dir = "/home/serigne-gueye/RECHERCHE/CMPUF/CODES/V10-15-07-2026/results/hospital/results2026_v1_";
		char * dirgreedy = "/home/serigne-gueye/RECHERCHE/CMPUF/CODES/V10-15-07-2026/results/hospital/solgreedy_v1_";
		char * diropt = "/home/serigne-gueye/RECHERCHE/CMPUF/CODES/V10-15-07-2026/results/hospital/solopt_v1_";
		char * dircent = "/home/serigne-gueye/RECHERCHE/CMPUF/CODES/V10-15-07-2026/results/hospital/solcent_v1_";
		char * param = argv[8];
		std::string typeinput(argv[9]);
	
		FILE * file1;
		FILE * file2;

     	        char * format = "%s & %d & %d & %3.2f &  %3.2f & %3.2f & %10.2f & %10.2f & %10.2f & %10.2f &%10.2f & %3.3f & %3.3f & %3.3f & %3.3f & %3.3f & %3.3f & %3.3f & %3.3f & %3.3f & %3.3f & %3.3f & %3.3f & %3.3f & %3.3f & %3.3f \\\\ \n";
     	        char * format1 = "%d&%d&%3.2f&%3.2f&%3.2f&%10.2f&%10.2f&%10.2f&%10.2f&%10.2f&%3.3f&%3.3f&%3.3f&%3.3f&%3.3f&%3.3f&%3.3f&%3.3f&%3.3f&%3.3f&%3.3f&%3.3f&%3.3f&%3.3f&%3.3f\n";
		
                int ncap = 4;
		if(file1 = fopen(argv[1],"r"))
		{
			concat(nom,dir,argv);
                        concat(solgreedy,dirgreedy,argv);
                        concat(solopt,diropt,argv);
                        concat(solcent,dircent,argv);
			//cout << "Nom = " << nom << endl;
			//cin.get();			
                        // Nom du fichier où seront stockés les résultats des tests
			file2 = fopen(nom,"a");
			
                        fprintf(file2, "%s&%s&%s&%s&%s&%s&%s&%s&%s&%s&%s&%s&%s&%s&%s&%s&%s&%s&%s&%s&%s&%s&%s&%s&%s",
				"n",
				"m",
				"Alpha",
			        "Beta",
				"Gamma",
				"VOPT",
			        "VGREEDY",
				"VGREEDYLS",
				"VMEDIAN1",
				"VMEDIAN2",
				"cent_opt",
				"cent_greedy",
				"cent_greedyls",
				"cent_median1",
				"cent_median2",
				"coloc_opt",
				"coloc_greedy",
				"coloc_greedyls",
				"coloc_median1",
				"coloc_median2",
				"topt",
				"tgreedy",
				"tgreedyls",
				"tpmedian1",
				"tpmedian2");			

			while(! feof(file1) )
			{
				fscanf(file1,"%s\n",filename);

				m = atoi(argv[3]);        // Number of facilities
				gap = atof(argv[4]);      // Minimum Gap parameter for CPLEX. Under this value the CPLEX Branch-and-bound stop
				rayon = atof(argv[5]);    // Ray used in the local search constraint
				rayonbb = atof(argv[6]);  // Ray used in the local search constraint
				duree = atof(argv[7]);    // Maximum Branch-and-Bound processing time
				
				IloEnv env1;
				FloatMatrix xpmedian;      // Solution optimale p-median : alpha = 0, beta  = 0, gamma = 0.
				IloNumArray Cm;            // Capacité théorique. Cm[i]  = sum of the client demands assigned to i in the m-median optimal solution

				clock_t c_start = clock();
				/////////////////////////////////////
				// Solving m-median problem
                                // resolution du problème p-median.
				//cout << "Pmedian" << endl;
				//cin.get();
                                /////////////////////////////////////
				Pmedian cmpuf_pmedian(env1,filename,typeinput,m);
				cmpuf_pmedian.Init_X(env1);
				cmpuf_pmedian.Init_z(env1);
				cmpuf_pmedian.Init_contraintes(env1);
				cmpuf_pmedian.Init_obj(env1);
				cmpuf_pmedian.Conversion(env1);
				cmpuf_pmedian.Resolution();
				xpmedian = cmpuf_pmedian.XtoXpmedianSolution(env1); // Get p-median optimal solution
				Cm = cmpuf_pmedian.CapaciteMedian(env1);            // Get capacities of the facility
                                /////////////////////////////////////
				clock_t c_end = clock();
				tpmedian1 = (c_end - c_start) / CLOCKS_PER_SEC;
				//////////////////////////////////////////////
                                // For each couple of parameters gamma and alpha a 2-stage stochastic program is solved
                                //////////////////////////////////////////////
				for(gamma = 1 ; gamma <= 1.75 ; gamma = gamma + pas_gamma){
					beta = 0; 
					for(alpha = 0 ; alpha <= 1 ; alpha = alpha + pas_alpha){ // Version originale
						IloEnv env;
					        IloEnv env_heur;
					        
						IloNum centralizationopt = -1;
						IloNum centralizationgreedy = -1;
						IloNum centralizationgreedyls = -1;
						IloNum centralizationmedian1 = -1;
						IloNum centralizationmedian2 = -1;
						IloNum colocationopt  = -1;
						IloNum colocationgreedy = -1;
						IloNum colocationgreedyls = -1;
						IloNum colocationmedian1 = -1;
						IloNum colocationmedian2 = -1;
						IloNum VMEDIAN1 = -1;
						IloNum VMEDIAN2 = -1;
						IloNum VGREEDY = -1;
						IloNum VGREEDYLS = -1;
						IloNum VOPT = -1;
						bestL = IloNumArray(env);
						///////////////////////////////////
                                                // Object cmpuf construction
                                                ///////////////////////////////////
						CMpuf cmpuf(env,filename,typeinput,param,alpha,beta,gamma,m);						
						CMpuf cmpuf_heur(env_heur,filename,typeinput,param,alpha,beta,gamma,m);
						//cout << "Cmpuf heur : valeur objectif du CMPUF de la solution optimale pmedian" << endl;
						//cin.get();
                                                ///////////////////////////////////
                                                // Attributes initialization
                                                //////////////////////////////////
						IloNumArray numeros_heur(env_heur,cmpuf_heur.M());
						/********************************/
						cmpuf_heur.Init_cap(env_heur,Cm);
						cmpuf_heur.Init_T(env_heur,m-1,numeros_heur);
						cmpuf_heur.Init_C(env_heur);
						//cmpuf_heur.AffichageDonnees();
						 ///////////////////////////////////
                                                // Model construction
                                                //////////////////////////////////
						cmpuf_heur.Init_x(env_heur);
						cmpuf_heur.Init_z(env_heur);
						cmpuf_heur.Init_contraintes(env_heur);
						cmpuf_heur.Init_obj(env_heur);
						IloCplex cplex_heur(env_heur);
						cplex_heur.setOut(env_heur.getNullStream()); 
						cplex_heur.extract(cmpuf_heur);
						cplex_heur.setParam(IloCplex::SimDisplay,0);
						/////////////////////////////////////////
						// Solving m-median problem
                                                /////////////////////////////////////////
						VMEDIAN1 = cmpuf_heur.Resolution(env,cplex_heur,xpmedian); // The model in cmpuf by fixing the location to the values in xpmedian
						                                                           // Get the optimal value 
						IloNumArray Lm = cmpuf_heur.Init_L(env,xpmedian);          // Get optimal m-median location

						bestL = Lm;
						bestV = VMEDIAN1;
						/*
						cout << "VMEDIAN1 =\t"  << VMEDIAN1 << endl;
						cout << " L median =\t" << Lm << endl;
						*/
						//cmpuf_heur.AffichageCoutTransport();
						//cin.get();
						/////////////////////////////////////////
						// Init object attributes
                                                /////////////////////////////////////////
						IloNumArray numeros(env,cmpuf.M());
						/********************************/
						cmpuf.Init_cap(env,Cm);
						cmpuf.Init_T(env,m-1,numeros);
						cmpuf.Init_C(env);
						//cmpuf_heur.AffichageDonnees();
						cmpuf.Init_x(env);
						cmpuf.Init_z(env);
						cmpuf.Init_contraintes(env);
						cmpuf.Init_obj(env);
						/////////////////////////////////////////
						// Cplex object associated to cmpuf
                                                /////////////////////////////////////////
						IloCplex cplex(env);
						cplex.extract(cmpuf);
						cplex.setOut(env.getNullStream());
						cplex.setParam(IloCplex::MIPDisplay,0);

						double tcplex;
						tcplex = cplex.getCplexTime();
						//cout << "Greedy Heuristic" << endl;
						//cin.get();

						IloNumArray Lg = cmpuf.MyGreedyHeuristic(env,cplex);  //Solve the model with the greedy heuristic. Get the optimal location
						VGREEDY = cmpuf.Greedy();                             // Get the optimal value
						/*
						cout << " ! VGREDDY ! = " << VGREEDY << endl;
						cout << " ! L GREDDY ! = " << Lg << endl;
						cin.get();
						*/
						if(VGREEDY > bestV){
							bestV = VGREEDY;
							bestL = Lg;
						}
						tgreedy = cplex.getCplexTime() - tcplex;
						/*
						cout << "Local Search à partir de la solution Greedy" << endl;
						cin.get();
						*/
						tcplex = cplex.getCplexTime();
						
						IloNumArray Lls1 = cmpuf.Localsearch(env,cplex,Lg,cmpuf.Greedy(),rayon); // Local Search Algorithm aroung the greedy solution
						tgreedyls = cplex.getCplexTime() - tcplex;
						VGREEDYLS = cmpuf.LS(); // Get optimal value of the local search algorithm 
						
						if(VGREEDYLS > bestV){
						      bestV = VGREEDYLS;
						      bestL = Lls1;
						  }
						  /*
						  cout << " ! VLS Greedy ! = " << VGREEDYLS << endl;
						  cout << " ! L LS Greedy ! = " << Lls1 << endl;
						  cin.get();
						  */
						  //cout << "Local Search à partir de la solution pmedian" << endl;
						  //cin.get();
						  tcplex = cplex.getCplexTime();
						  IloNumArray Lls2 = cmpuf.Localsearch(env,cplex,Lm,VMEDIAN1,rayon);  // Local Search Algorithm aroung the median solution
						  tpmedian2 = cplex.getCplexTime() - tcplex;
						  VMEDIAN2 = cmpuf.LS();

						  if(VMEDIAN2 > bestV){
							bestV = VMEDIAN2;
							bestL = Lls2;
						  }

						  //cout << " ! VLS Pmedian! = " << VMEDIAN2 << endl;
						  //cout << " ! L LS Pmedian ! = " << Lls2 << endl;
						  //cin.get();
						  /////////////////////////////////////////////////////////
						  // Computation of centralization and colocation values
						  /////////////////////////////////////////////////////////
						  centralizationgreedy = cmpuf.Centralization(Lg,cmpuf.Profit());
						  centralizationgreedyls = cmpuf.Centralization(Lls1,cmpuf.Profit());
						  centralizationmedian1 = cmpuf.Centralization(Lm,cmpuf.Profit());
						  centralizationmedian2 = cmpuf.Centralization(Lls2,cmpuf.Profit());
  						  colocationgreedy = cmpuf.Colocation(env,Lg);
						  colocationgreedyls = cmpuf.Colocation(env,Lls1);
						  colocationmedian1 = cmpuf.Colocation(env,Lm);
						  colocationmedian2 = cmpuf.Colocation(env,Lls2);
						  /////////////////////////////////////////////////
                                                  // Save greedy solution in the file "solgreedy
                                                  /////////////////////////////////////////////////
						  cmpuf.AffichageSolutions(Lls1,solgreedy);
						  //cplex_heur.setParam(IloCplex::SimDisplay,0);

						  if(strcmp(type,argv[2]) == 0) // Optimal solving
						  {
							//cout << endl << " ! MIP ! " << endl;
							//cplex.setParam(IloCplex::MIPDisplay,2);
							//cin.get();
							cplex.setParam(IloCplex::Param::Benders::Strategy,-1);
							//cplex.setParam(IloCplex::Param::Benders::Strategy,2);
							//cplex.setParam(IloCplex::Param::Emphasis::MIP,5);
							cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq,-1);
							//cplex.setParam(IloCplex::Param::Preprocessing::Symmetry,1);
							//cplex.setParam(IloCplex::Param::MIP::Cuts::Covers,3);
							//cplex.setParam(IloCplex::Param::MIP::Cuts::FlowCovers,2);
							//cplex.setParam(IloCplex::Param::Advance,0);		
							cplex.setParam(IloCplex::Param::TimeLimit,duree);
							cplex.setParam(IloCplex::EpGap,gap);
							//cplex_opt.setParam(IloCplex::CutLo,cmpuf_greedy.Opt());
							//cplex_opt.setParam(IloCplex::CutLo,VOPTTRANSPORT);
							//cplex.exportModel("cmpuf.lp");
							//cout << "export du model" << endl;
							//cin.get();
								
							//cout << "BestL = " << bestL << endl;
							//cout << "BestV = " << bestV << endl;								      
							//cmpuf.MIPStart(env,cplex,bestL);

							double tcplex;
							tcplex = cplex.getCplexTime();

							//cmpuf.Init_contraintes_x_x(env);
							VOPT = cmpuf.Resolution(env,cplex);
							//VOPT = cmpuf.LocalBranchAndBound(env,cplex,Lm,rayonbb);
							//VOPT = cmpuf.LocalBranchAndBound(env,cplex,Lls1,rayonbb);
							//cmpuf_opt.AffichageSolutions();
							//cin.get();
									
							topt = cplex.getCplexTime() - tcplex;
							IloNumArray L  = cmpuf.Init_L(env);

							//cout << "VOPT =\t"  << VOPT << endl;
							//cout << " L opt =\t" << L << endl;
							//cin.get();
		
							centralizationopt = cmpuf.Centralization(L,cmpuf.Profit());
							cmpuf.AffichageCoutTransport(env,centralizationopt,solcent);
							colocationopt = cmpuf.Colocation(env,L);
							cmpuf.ResolutionScenarios(env);
							cmpuf.AffichageSolutionsScenarios(L,solopt);
							//cmpuf.GenerationSolutionsScenarios(env,L,solopt);
							//cout << "ici 6" << endl;
							//cin.get();

					      }

					      printf( format,
						      NomApresDernierSeparateur(filename),
						      cmpuf.N(),
    						      cmpuf.M(),
						      cmpuf.Alpha(),
						      cmpuf.Beta(),
						      cmpuf.Gamma(),
						      VOPT,
						      VGREEDY,
						      VGREEDYLS,
						      VMEDIAN1,
						      VMEDIAN2,
						      centralizationopt,
						      centralizationgreedy,
						      centralizationgreedyls,
						      centralizationmedian1,
						      centralizationmedian2,
						      colocationopt,
						      colocationgreedy,
						      colocationgreedyls,
						      colocationmedian1,
						      colocationmedian2,
						      topt,
						      tgreedy,
						      tgreedyls,
						      tpmedian1,
						      tpmedian2);

									                                          
					      fprintf(file2, format1,
							cmpuf.N(),
							cmpuf.M(),
							cmpuf.Alpha(),
							cmpuf.Beta(),
							cmpuf.Gamma(),
							VOPT,
						        VGREEDY,
							VGREEDYLS,
							VMEDIAN1,
							VMEDIAN2,
							centralizationopt,
							centralizationgreedy,
							centralizationgreedyls,
							centralizationmedian1,
							centralizationmedian2,
							colocationopt,
							colocationgreedy,
							colocationgreedyls,
							colocationmedian1,
							colocationmedian2,
							topt,
							tgreedy,
							tgreedyls,
							tpmedian1,
							tpmedian2);

						fflush(file2);
								
						env.end();
						env_heur.end();
				    }
				}
				env1.end();
			}
		}

		fclose(file1);
		fclose(file2);      			
	}
	else
		cout << "Erreur sur le nombre d'arguments" << endl;
}
