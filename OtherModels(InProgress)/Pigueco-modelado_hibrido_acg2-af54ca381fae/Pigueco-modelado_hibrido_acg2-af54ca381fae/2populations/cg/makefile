EXE : MainV0.o Intracellular.o Read.o Utilities.o PDESolver.o
	gcc -o EXE MainV0.o Intracellular.o Read.o Utilities.o PDESolver.o 

Main.o : MainV0.c Header.h
	gcc -c MainV0.c

Intracellular.o : Intracellular.c Header.h
	gcc -c Intracellular.c

Read.o : Read.c Header.h
	gcc -c Read.c

Utilities.o : Utilities.c Header.h
	gcc -c Utilities.c

PDESolver.o : PDESolver.c Header.h
	gcc -c PDESolver.c
