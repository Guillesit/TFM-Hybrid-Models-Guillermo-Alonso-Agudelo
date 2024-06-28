gcc Main.c Read.c PDESolver.c Utilities.c Gillespie.c -o EXE -lm
#gcc -g -O3 -Wall -o tau.exe main_space.c birth.c sh 
# -lm -lgsl -lgslcblas -lrt
# ./tau.exe edadesJuan.dat Init_oxy Param_stochastic Parsim
./EXE Probando/pop12 Probando/oxy12 Probando/pars_pop Probando/pars_sim a
