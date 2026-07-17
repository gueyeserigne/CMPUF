/**
 * Author: Serigne Gueye
 * Date: 17 July 2026
 */

// Fonction d'echange entre les elements vect[i] et vect[j]
// et entre ordre[i] et ordre[j]

void ROPEchanger(IntVector & ordre, IloNumArray & vect, IloInt i, IloInt j);

// Fonction de tri des elements de vect en gardant
// l'ordre des indices des elements tries dans 
// le vecteur ordre. Ceci grace a la methode
// QuickSort 

void ROPTriRapide(IloNumArray & vect, IntVector & ordre, 
                  IloInt gauche, IloInt droite,
		  IloInt flag);


// Fonction de tri des elements de vect en gardant
// l'ordre des indices des elements tries dans 
// le vecteur ordre.  On copie les elements de v 
// dans un vecteur tampon pour eviter sa          
// modification et on initiale ordre       

void ROPTri(IloNumArray & vect, IntVector & ordre,
	    IloInt gauche, IloInt droite,
	    IloInt flag);

// Fonction qui retourne le maximum des composantes 
// d'un vecteur V

IloNum ROPMaxVectFloat(IloNumArray & V , IloInt n);
