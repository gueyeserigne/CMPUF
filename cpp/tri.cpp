/**
 * Author: Serigne Gueye
 * Date: 17 July 2026
 */

#include <ilcplex/ilocplex.h>

ILOSTLBEGIN
typedef IloArray<IloNumArray> FloatMatrix;
typedef IloArray<IloInt> IntVector;
typedef IloArray<IntVector> IntMatrix;

#include "../h/tri.h"

#define INFINI 100000



//////////////////////////////////////////////////////////
// But : Echange entre les elements vect[i] et vect[j]  //
//       et entre ordre[i] et ordre[j]                  //
//                                                      //
// Parametres :                                         //
//      ordre : Vecteur des indices                     //
//       vect : Vecteur des valeurs                     //
//          i : Indice de la valeur a echanger          //
//          j : Indice de la valeur a echanger          //
// Retours :                                            //
//      ordre est argument de retour ou ordre[i] a ete  //
//      echanger avec ordre[j].                         //
//      v est argument de retour ou ordre[i] a ete      //
//      echanger avec v[j].                             //
//////////////////////////////////////////////////////////

void ROPEchanger(IntVector & ordre, IloNumArray & vect, IloInt i, IloInt j)
{
	IloInt tamp;
	IloNum tampon;
	
	tampon = vect[i];
	vect[i] = vect[j];
	vect[j] = tampon;
				
	tamp = ordre[i];
	ordre[i] = ordre[j];
	ordre[j] = tamp;

	
}
 

//////////////////////////////////////////////////////////
// But : Tri des elements de v en gardant               //
//       l'ordre des indices des elements tries dans    //
//       le vecteur ordre. Ceci grace a la methode      //
//       QuickSort                                      // 
//                                                      //
// Parametres :                                         //
//       vect : Vecteur des valeurs                     //
//      ordre : Vecteur des indices                     //
//          i : Indice de la valeur a echanger          //
//          j : Indice de la valeur a echanger          //
// Retours :                                            //
//      ordre est argument de retour ou ordre[i] a ete  //
//      echanger avec ordre[j].                         //
//      v est argument de retour ou ordre[i] a ete      //
//      echanger avec v[j].                             //
//////////////////////////////////////////////////////////

void ROPTriRapide(IloNumArray & vect, IntVector & ordre, 
                  IloInt gauche, IloInt droite,
		  IloInt flag)
{
	IloInt i,dernier;
	
	
	if(gauche >= droite)
		return;
		
	ROPEchanger(ordre,vect,gauche,(gauche+droite)/2);
	
	dernier = gauche;
	
	if(flag)
	{
	    for(i = gauche+1; i<= droite;i++)
		if(vect[i] < vect[gauche])
		    ROPEchanger(ordre,vect,++dernier,i);
	}
	else
	{
	    for(i = gauche+1; i<= droite;i++)
		if(vect[i] >= vect[gauche])
		   ROPEchanger(ordre,vect,++dernier,i);
	}

	
	ROPEchanger(ordre,vect,gauche,dernier);
	
	ROPTriRapide(vect,ordre,gauche,dernier-1,flag);
	ROPTriRapide(vect,ordre,dernier+1,droite,flag);
}



//////////////////////////////////////////////////////////
// But : Tri des elements de v en gardant               //
//       l'ordre des indices des elements tries dans    //
//       le vecteur ordre. On copie les elements de v   //
//       dans un vecteur tampon pour eviter sa          //
//       modification et on initiale ordre              //
//                                                      //
// Parametres :                                         //
//          v : Vecteur des valeurs                     //
//      ordre : Vecteur des indices                     //
//          i : Indice de la valeur a echanger          //
//          j : Indice de la valeur a echanger          //
// Retours :                                            //
//      ordre est argument de retour ou ordre[i] a ete  //
//      echanger avec ordre[j].                         //
//      v est argument de retour ou ordre[i] a ete      //
//      echanger avec v[j].                             //
//////////////////////////////////////////////////////////



void ROPTri(IloNumArray & vect, IntVector & ordre,
	    IloInt gauche, IloInt droite,
	    IloInt flag)
{
        IloEnv env;
	
	IloNumArray tvect(env,(droite+1));
	IloInt i;
	

	for(i=0;i<(droite+1);i++)
	  tvect[i]  = vect[i];
	
	//for(i=0;i<(droite+1);i++)
	  //ordre[i] = i;

	ROPTriRapide(tvect,ordre,gauche,droite,flag);
	
 	env.end();
	
}



/////////////////////////////////////////////////////
// But : Determiner l'element maximum d'un vecteur //
//                                                 //
// Parametres :                                    //
//          V : Vecteur                            //
//          n : Taille de V                        //
//                                                 //
// Retour :                                        //
//   ROPMaxVectFloat est la valeur maximum des     //
//   composantes de V                              //
/////////////////////////////////////////////////////

IloNum ROPMaxVectFloat(IloNumArray & V , IloInt n)
{
  IloNum max = - INFINI;
  IloInt i;

  for(i = 0 ; i < n ; i++)
    if(V[i] > max)
      max = V[i];

  return(max);

}
