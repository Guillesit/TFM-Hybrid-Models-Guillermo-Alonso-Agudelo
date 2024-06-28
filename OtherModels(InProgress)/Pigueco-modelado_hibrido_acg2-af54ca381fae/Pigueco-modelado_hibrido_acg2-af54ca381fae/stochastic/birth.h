#include <math.h>							
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

/*int ode(double t, const double y[], double f[], void *params);
int jac();
int birth(double tini, double tau, double y[],  double datos[]);
int odeOx (double t,const double y[], double f[], void *params);
int Oxygen(double tini, double tfin,double y[],double datos[] );*/



int jac();
double birth(double ox);
//double birth2(double ox);
int odeOx (double t,const double y[], double f[], void *params);
int Oxygen(double tini, double tfin,double y[],double datos[] );


double reactions(int end, double **gilAux, int *J, double **age, double **b, double nu, double R, double L);
void gillespie(int i, int j, double **age, int *J, double **ag, double **b);







