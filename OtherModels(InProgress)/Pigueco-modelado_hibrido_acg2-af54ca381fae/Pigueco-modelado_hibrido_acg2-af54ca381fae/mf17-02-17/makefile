EXE : Main.o Intracellular.o Read.o Utilities.o PDESolver.o
	gcc -o EXE Main.o Intracellular.o Read.o Utilities.o PDESolver.o 

Main.o : Main.c Header.h
	gcc -c Main.c

Intracellular.o : Intracellular.c Header.h
	gcc -c Intracellular.c

Read.o : Read.c Header.h
	gcc -c Read.c

Utilities.o : Utilities.c Header.h
	gcc -c Utilities.c

PDESolver.o : PDESolver.c Header.h
	gcc -c PDESolver.c
