# Compilation et exécution de CMPUF

Ce document décrit la compilation et l'exécution de la cible active
`cmpuf`.

## 1. Prérequis

Le projet est prévu pour Linux x86-64 et nécessite :

- un compilateur C++ compatible C++11 (GCC/g++ recommandé) ;
- GNU Make 4 ou une version récente de CMake ;
- IBM ILOG CPLEX Studio avec CPLEX et Concert ;
- une licence CPLEX valide au moment de l'exécution.

Par défaut, les fichiers de compilation recherchent CPLEX Studio dans :

```text
/opt/ibm/ILOG/CPLEX_Studio201
```

L'installation doit notamment contenir :

```text
cplex/include
concert/include
cplex/lib/x86-64_linux/static_pic
concert/lib/x86-64_linux/static_pic
```

Le Makefile référence également `cpoptimizer` pour les anciennes cibles.

## 2. Compilation avec CMake (recommandée)

Depuis la racine du projet :

```bash
cmake --preset release
cmake --build --preset cmpuf --parallel
```

Le premier appel configure un build Release dans :

```text
build/cmake-release
```

Le second compile la cible active et crée :

```text
bin/cmpuf
```

Pour utiliser une installation CPLEX située ailleurs :

```bash
cmake --preset release \
  -DCPLEX_STUDIO_DIR=/chemin/vers/CPLEX_Studio
cmake --build --preset cmpuf --parallel
```

Pour repartir d'une configuration CMake vierge, supprimer le répertoire
`build/cmake-release`, puis relancer la configuration.

## 3. Compilation avec GNU Make

Le Makefile peut être appelé depuis la racine du projet :

```bash
make -j2 cmpuf
```

La cible par défaut étant `cmpuf`, ceci est équivalent à :

```bash
make -j2
```

Les objets sont écrits dans `build/obj` et l'exécutable dans
`bin/cmpuf`.

Pour indiquer une autre installation CPLEX :

```bash
make -j2 cmpuf \
  CPLEX_STUDIO_DIR=/chemin/vers/CPLEX_Studio
```

Pour nettoyer les objets et les répertoires de build :

```bash
make clean
```

Cette commande ne supprime pas les exécutables déjà présents dans `bin`.

## 4. Ligne de commande

L'exécutable attend exactement neuf arguments :

```bash
./bin/cmpuf \
  LISTE_INSTANCES \
  MODE \
  M \
  GAP \
  RAYON_LOCAL \
  RAYON_BB \
  LIMITE_TEMPS \
  FICHIER_PARAMETRES \
  TYPE_ENTREE
```

Signification des arguments :

1. `LISTE_INSTANCES` : fichier texte contenant un chemin d'instance par
   ligne ;
2. `MODE` : utiliser `e` pour activer la résolution exacte CPLEX, en plus
   des heuristiques ;
3. `M` : nombre d'établissements à localiser ;
4. `GAP` : gap relatif d'arrêt du MIP CPLEX, par exemple `0.00` ;
5. `RAYON_LOCAL` : nombre de nœuds proches autorisés dans la recherche
   locale ;
6. `RAYON_BB` : rayon prévu pour le branch-and-bound local ; ce chemin est
   actuellement désactivé dans le programme principal, mais l'argument reste
   obligatoire ;
7. `LIMITE_TEMPS` : limite de temps CPLEX, en secondes ;
8. `FICHIER_PARAMETRES` : fichier définissant les états probabilistes et
   les multiplicateurs de capacité ;
9. `TYPE_ENTREE` : utiliser `geometric` pour le format géométrique ; utiliser
   `hospital` pour le format standard avec un vecteur de demandes suivi de la
   matrice de distances.

Exemple depuis la racine du projet :

```bash
./bin/cmpuf \
  ../BENCHMARK/calgaryhospital15 \
  e \
  4 \
  0.00 \
  5 \
  40 \
  30000 \
  parameters/param \
  hospital
```

L'ancien exemple du README est lancé depuis `bin` :

```bash
cd bin
./cmpuf \
  ../../BENCHMARK/calgaryhospital15 \
  e 4 0.00 5 40 30000 \
  ../parameters/param \
  hospital
```

Les chemins relatifs sont interprétés depuis le répertoire courant au moment
de l'exécution.

## 5. Format de la liste d'instances

Le premier argument n'est pas directement une instance. Il s'agit d'un
fichier contenant un chemin d'instance par ligne :

```text
/chemin/vers/instance1
/chemin/vers/instance2
```

Éviter les espaces dans les chemins : le programme lit chaque entrée avec
`fscanf(..., "%s", ...)`.

## 6. Format d'une instance

Chaque fichier d'instance doit contenir, dans cet ordre :

1. le nombre de nœuds/clients `n` ;
2. `n` demandes ;
3. une matrice carrée de `n × n` coûts ou distances.

Schématiquement :

```text
n
w[0] w[1] ... w[n-1]
d[0][0] ... d[0][n-1]
...
d[n-1][0] ... d[n-1][n-1]
```

Les séparateurs sont des espaces ou des retours à la ligne.

## 7. Format du fichier de paramètres

Le fichier fourni avec le projet, `parameters/param`, contient :

```text
4
0.10 0.40 0.50
0.00 0.33 0.66
```

Son format est :

1. `ncap`, nombre d'états de capacité par établissement ;
2. `ncap - 1` poids probabilistes `walpha` ;
3. `ncap - 1` multiplicateurs de capacité `wbeta`.

Le dernier état de capacité correspond automatiquement à la pleine capacité,
avec une probabilité `1 - alpha`.

Le nombre total de scénarios générés est :

```text
K = ncap^M
```

La mémoire et le temps de calcul peuvent donc augmenter très rapidement avec
`M`.

## 8. Calculs exécutés

Pour chaque instance, le programme :

1. résout un problème p-médiane ;
2. calcule les capacités nominales des établissements ;
3. parcourt plusieurs valeurs de `gamma` et `alpha` ;
4. construit tous les scénarios de capacité ;
5. évalue la solution p-médiane ;
6. lance l'heuristique gloutonne intégrée à `CMpuf` ;
7. lance deux recherches locales ;
8. en mode `e`, résout le modèle exact avec CPLEX ;
9. écrit les valeurs et les solutions dans `results/hospital`.

Attention : les chemins de sortie sont actuellement codés en dur dans
`cpp/cmpuf.cpp`. Le projet doit donc rester à son emplacement
actuel, ou ces chemins doivent être adaptés dans le code.

## 9. Vérifications rapides

Vérifier que l'exécutable existe :

```bash
file bin/cmpuf
```

Afficher l'erreur d'utilisation sans lancer de calcul :

```bash
./bin/cmpuf
```

La réponse attendue est :

```text
Erreur sur le nombre d'arguments
```

Cette vérification confirme seulement le démarrage de l'exécutable. Une
exécution complète nécessite une instance valide et peut lancer des calculs
CPLEX longs.

## 10. Problèmes fréquents

### En-têtes CPLEX introuvables

Erreur typique :

```text
fatal error: ilcplex/ilocplex.h: No such file or directory
```

Renseigner le bon chemin avec `CPLEX_STUDIO_DIR`.

### Bibliothèques CPLEX introuvables à l'édition de liens

Vérifier les sous-répertoires :

```text
cplex/lib/x86-64_linux/static_pic
concert/lib/x86-64_linux/static_pic
```

Si l'installation emploie une autre plateforme ou un autre format, adapter
`CPLEX_SYSTEM` et `CPLEX_LIBFORMAT` dans CMake ou fournir les variables
correspondantes à Make.

### Erreur de licence

Une compilation réussie ne garantit pas l'accès au solveur. Vérifier la
configuration de licence IBM ILOG CPLEX sur la machine d'exécution.

### Calcul très volumineux

Avec quatre états de capacité, le nombre de scénarios vaut `4^M`. Réduire
`M`, le nombre d'états ou la limite de temps pour les premiers essais.
