/**
 * Author: Serigne Gueye
 * Date: 17 July 2026
 */

/****************************************************************/
/* Serigne Gueye						*/
/* Projet IEF FLOUE MIT-ZLC					*/
/****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <cstring>
#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace boost;
using namespace std;

typedef property<edge_weight_t, double> weight;
typedef adjacency_list<listS, vecS, undirectedS, no_property, weight > Graph;
typedef std::pair<int, int> Edge;

int n;
const int infinity = 1e20;
//const int ns = 3; // Nombre de classes
const int ns = 3; // Nombre de classes
const int np = 3; // Nombre de probabilites
const int nm = 2; // Nombre de services
const int npsmax = 500;
int nps = 0;
int *P[npsmax];
int S[ns] = {1,2,3};
//int S[ns] = {1};
//double p[np] ={0.05,0.1,0.3,0.5};
double p[np] ={0.05,0.1,0.2};
//int m[nm] = {3,5,10};
int m[nm] = {5,10};
int ** d;
int ** distances;
double * h;
int seed = 13345678;
bool partie[np];


void demandes()
{
	//cout << "Debut demandes....." << endl;
	//cin.get();

	int i;
	double somme = 0;

	h = new double[n];

	for(i = 0 ; i < n ; i++)
	{
		h[i] = rand() % 100;		
		somme += h[i];
	}

	//for(i = 0 ; i < n ; i++)
	//	h[i] = (h[i] / somme);		

	//cout << "fin demandes....." << endl;
	//cin.get();

}

void combinaison(int i, int s, int  prof)
{
	//cout << "Debut combinaison....." << endl;
	//cin.get();

	int j,k;

	if(prof < s)
	{
		for(j = i+1 ; j < np ; j++)
		{
			partie[j] = 1;
			combinaison(j,s,prof+1);
			partie[j] = 0;
		}
	}
	else
	{
		P[nps] = new int[s];
		k = 0;

		for(j = 0 ; j < np ; j++)
			if(partie[j] == 1)
			{
				P[nps][k] = j;
				k++;
			}

		nps++;
	}

	//cout << "fin combinaison....." << endl;
	//cin.get();

}


void lecture(char * filename)
{
	//cout << "Debut lecture....." << endl;
	//cin.get();

	int i,j;
	double w;
	char* format = "%d %d %lf\n";

	FILE * file = fopen(filename,"r");
	
  	if(! file)
    		std::cerr << "Erreur a l'ouverture du fichier " << filename << std::endl;
  	else
  	{
		fscanf(file,"%d\n",&n);

		//cout << "n = " << n << endl;
		//cin.get();

		d = new int*[n];

		for(i = 0 ; i < n ; i++)
			d[i] = new int[n];

		for(i = 0 ; i < n ; i++)
			for(j = 0 ; j < n ; j++)
				d[i][j] = infinity;

		while(! feof(file) )
		{
			fscanf(file,"%d %d %lf\n",&i,&j,&w);


			// Geometric

			if((100*w - floor(100*w)) < (ceil(100*w) - 100*w))
				d[i][j] = d[j][i] = floor(100*w);
			else
				d[i][j] = d[j][i] = ceil(100*w);

			// Grid

			//d[i][j] = d[j][i] = w;
			
			// Tree

			//d[i][j] = d[j][i] = floor(100*w);
		}

		fclose(file);
	}

	//cout << "fin lecture....." << endl;
	//cin.get();
}


void ecrire(char * filename, int m, double * p)
{
	//cout << "Debut ecrire....." << endl;
	//cin.get();

	int i,j;

	FILE * file = fopen(filename,"w");
	
  	if(! file)
    		std::cerr << "Erreur a l'ouverture du fichier " << filename << std::endl;
  	else
  	{
		fprintf(file,"%d\n",n);
		fprintf(file,"%d\n",m);

		fprintf(file,"[%lf",h[0]);

		for(i = 1 ; i < (n-1) ; i++)
			fprintf(file,",%lf",h[i]);

		fprintf(file,",%lf]\n",h[n-1]);

		fprintf(file,"[%lf",p[0]);

		for(i = 1 ; i < (m-1) ; i++)
			fprintf(file,",%lf",p[i]);

		fprintf(file,",%lf]\n",p[m-1]);

		fprintf(file,"[");

		for(i = 0 ; i < (n-1) ; i++)
		{
			fprintf(file,"[");

			for(j = 0 ; j < (n-1) ; j++)
				fprintf(file,"%d,",distances[i][j]);

			fprintf(file,"%d],\n",distances[i][n-1]);
		}

		fprintf(file,"[");

		for(j = 0 ; j < (n-1) ; j++)
			fprintf(file,"%d,",distances[n-1][j]);

		fprintf(file,"%d]",distances[n-1][n-1]);


		fprintf(file,"]");
	}

	//cout << "fin ecrire....." << endl;
	//cin.get();


}


void shortestpaths()
{
	//cout << "Debut shortest paths....." << endl;
	//cin.get();

	int i,j;

	Graph clique(n);

	for(i = 0 ; i < n ; i ++)
		for(j = i+1 ; j < n ; j ++)
		   add_edge(i, j, weight(d[i][j]), clique);

	srand( time(NULL) );

	property_map<Graph, edge_weight_t>::type poids = get(edge_weight, clique);

	// vector for storing distance property
	std::vector<double> dist(n);
	std::vector<double> pred(n);

	//cout << "i = " << i << " - j = " << j << endl;

	distances = new int*[n];

	for(i = 0 ; i < n ; i++)
		distances[i] = new int[n];
	/*
	for(i = 0 ; i < n ; i ++)
		for(j = i+1 ; j < n ; j ++)
		   cout << "d[" << i << "]" << "[" << j << "] = " << d[i][j] << endl;

	cin.get();*/

        graph_traits<Graph>::vertex_iterator vi;

	for(i = 0 ; i < n ; i ++)
	{
		dijkstra_shortest_paths(clique, i, distance_map(&dist[0]));
	        //dijkstra_shortest_paths(clique, i, predecessor_map(&pred[0]).distance_map(&dist[0]));

		for(vi = vertices(clique).first; vi != vertices(clique).second; ++vi)
		{
			distances[i][*vi] = dist[*vi];

			//cout << "distance(" << *vi << ") = " << dist[*vi] << endl;
		}

		//std::cin.get();
	}

	//cout << "fin shortest paths....." << endl;
	//cin.get();
}



/****************************************************************/
/*								*/
/*			PROGRAMME PRINCIPAL			*/
/*								*/
/****************************************************************/
int main(int argc, char ** argv)
{
	int i,j,pos,s,cpt,k,l,r,buf,buf1;
	double * pbis;

	for(i = 0 ; i < np ; i++)
		partie[i] = 0;

	if( argc == 3 )    
	{
		char filename[500];
		char newfilename[500];
		char tmp[10];
		FILE * file1;
		FILE * file2;
		char * dir = "/home/gueye/RECHERCHE/FLP/MPUF/CODES/BENCHMARK/GEOMETRIC_V5/";
		//char * dir = "/home/gueye/RECHERCHE/FLP/MPUF/CODES/BENCHMARK/GEOMETRIC_V1/";
		//char * dir = "/home/gueye/RECHERCHE/FLP/MPUF/CODES/BENCHMARK/GRID/";
		//char * dir = "/home/gueye/RECHERCHE/FLP/MPUF/CODES/BENCHMARK/TREE/";
		char * dir1 = "/home/gueye/RECHERCHE/GRAPHS/GRAPHGENERATOR/Randgraph/TREE/";

		if( ( file1 = fopen(argv[1],"r")) && ( file2 = fopen(argv[2],"a")) )
		{
			while(! feof(file1) )
			{
				int res = fscanf(file1,"%s\n",filename);
				
				if(res != EOF)
				{
					lecture(filename);
					shortestpaths();

					for(j = 0 ;  j < nm ; j++)
					{
						if(m[j] <= (n/2))
						{
							pbis = new double[m[j]];

							for(s = 0 ;  s < ns ; s++)
							{
								//cout << "S[" << s << "] = " << S[s] << endl;
								//cin.get();

								nps = 0;

								combinaison(-1,S[s],0);

								for(pos = 0 ;  pos < nps ; pos++)
								{
									l = 0;
									buf = m[j];
									buf1 = S[s];

									//cout << "S[" << s << "] = " << S[s] << endl;
									//cout << "buf " << buf << endl;

									for(k = 0 ; k < (S[s]-1) ; k++)
									{				
										r = (rand() % (buf - buf1+1));
										r++;
										buf -= r;
										buf1--;
			
										//for(i = 0 ;  i < m[j] ; i++)
										//	pbis[i] = p[P[pos][rand() % S[s]]];	
										for(i = 0 ;  i < r ; i++)
											pbis[l+i] = p[P[pos][k]];

										l += r;

										//cout << "l- = " << l << endl;
										//cout << "r- = " << r << endl;
										//cin.get();
									}

									r = m[j] - l;
									//cout << "r = " << r << endl;
									//cin.get();

									for(i = 0 ;  i < r ; i++)
										pbis[l+i] = p[P[pos][S[s]-1]];

									demandes();

									strcpy(newfilename,dir);
									//strcat(newfilename,filename+59);
									strcat(newfilename,filename+67);
									//cout << "newfilename = " << newfilename << endl;
									//cin.get();

									sprintf(tmp,"_%d_%d_%d",m[j],S[s],pos);
									strcat(newfilename,tmp);

									//cout << "newfilename = " << newfilename << endl;
									//cin.get();

									ecrire(newfilename,m[j],pbis);

									fprintf(file2,"%s\n",newfilename);
									fflush(file2);
								}

								for(pos = 0 ;  pos < nps ; pos++)
									free(P[pos]);
							}

							free(pbis);
						}
					}

					free(h);

					for(i = 0 ; i < n ; i++)
					{
						free(distances[i]);
						free(d[i]);
					}

					free(distances);
					free(d);
				}
			}
			
			fclose(file1);
			fclose(file2);
		}
		else
		{
			std::cout << " ERREUR a l'ouverture de FILE 1 ou FILE 2 " << std::endl;
			std::cin.get();
		}
	}
	else
	{
		std::cout << "Nombre d'argument insuffisant !" << std::endl;
		std::cin.get();
	}
}
