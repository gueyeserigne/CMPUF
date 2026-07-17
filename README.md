"In the memory of Mozart BC Menezes"

Versions gérées par Git

Compilation avec Make depuis la racine :

make -j2 cmpuf

Compilation avec CMake depuis la racine :

cmake --preset release
cmake --build --preset cmpuf --parallel

Exécution depuis la racine :

./bin/cmpuf ../BENCHMARK/calgaryhospital15 e 4 0.00 5 40 30000 parameters/param hospital

