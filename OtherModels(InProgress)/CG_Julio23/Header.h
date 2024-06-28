

#define positive_part(A) ((A)>0.0 ? (A) : 0.0)
#define min(A,B) ( (A) < (B) ? (A) : (B) )
#define max(A,B) ( (A) < (B) ? (B) : (A) )
#define INFINITE_TIME 1000000000000000000.0 //A time instant we will hopefully never reach
#define N_OUTPUTS 100 
#define TOL 1e-10 //Para controlar Newton


//Some system constants
#define c_0 1.1
#define r_cr 1.004
#define aminus 30600 //8.5 HORAS EN SEGUNDOS //8250 era el valor primitivo (uds. desconocidas)
#define typical_oxy 0.00001 //Valor tipico de la concentracion del oxigeno, INVENTADO
                    //Medido en moles por metro
#define beta 0.2
#define severe_hypoxia 0.01 //1% del valor tipico


//Some functions determining division rates
#define c_cr(x) (1.0-0.4*log((0.51+(1.0-1.0/(1.0-0.45185/(x*x)))/2.0)/0.0085))







#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>



////////////////****************************////////////////////////////////

//STRUCTS

struct sim_parameters{
    double death_rate_hat; //dimensionless 
    double tau_p; //dimensional (recoprocal of birth rate), SI units
    double p6_over_p3; //ratio of p's
    double diff_coef_pop; //dimensional, SI units, MAYBE NOT NEEDED
    double survival_rate;
    double aplus; //DEPRECATED
    double aminus_hat; //dimensionless
    double k_decay_hat;
    double k_consumption_hat; //dimensionless
    double source_oxygen_hat; ///dimensionless
    double diff_coef_oxygen; //dimensional, SI units, MAYBE NOT NEEDED
    double diff_coefs_ratio; //dimensionless, for the adimensional oxygen evolution
    double diff_coef_eff; //dimensionless, for the adimensional population evolution
    long n_xslots;
    double Delta_x_hat; //dimensionless
    double Delta_t_hat; //dimensionless
    double critical_oxy_hat; //c_crit del modelo de ciclo celular, adimensional
    double typical_lengthscale; //dimensional, SI units
    double CFL_number; //dimensionless time (check aagin)
};




//SUBROUTINES


///////////////////
//READ FILE

long Read_Init_Space(
                     FILE *DATA_FILE,
                     double **pinitial_population
                     );


long Read_Init_Space_Oxygen(
                            FILE *DATA_FILE,
                            double **pinitial_oxygen
                            );


void Read_params_population(struct sim_parameters *params,
                            FILE *PARAMETERS
                            );


void Read_params_sim(struct sim_parameters *params,
                     FILE *SIM_DATA,
                     double *ptstop,
                     double *pDelta_x,
                     double *pDelta_t,
                     long *pn_files
                     );



///////////////////////
//PDE SOLVER FILE


//WE always include a survival fraction (take =1 for no therapy regimes)
void EvalPopulation(struct sim_parameters *params,
                    double *kpop, //rhs readout
                    double *density_previous,
                    double *threshold
                    );

void EvalOxygen(struct sim_parameters *params,
                    double *koxy, //rhs readout
                    double *oxygen_concentration, //old status oxygen
                    double *density_of_individuals
                    );


void  Eulerpop_oxy(struct sim_parameters *params,
                double delta_t,
                double *density_of_individuals,
                double *division_threshold,
                double *oxygen_concentration
                );


void GlobalPDE_Handler( 
                struct sim_parameters *params,
                double *density_of_individuals,
                double *division_threshold,
                double *oxygen_concentration
                );




//////////////////////
//UTILITIES FILE

void reverse(char s[]);

void itoa(int n,
          char s[]
          );

void Print_Vector(
				FILE *OUTPUT_DATA,
				char *output_path, 
				int iteration, //Si ponemos char no podemos pasar de 127 ficheros
				long vector_size, 
				double *vector
				);

void Print_VectorLocation(
        FILE *OUTPUT_DATA,
        char *output_path, 
        char iteration, 
        long vector_size, 
        double *vector
        );


double Get_Eigenvalue(
                    double *division_threshold, 
                    long xslot, 
                    double death_rate,
                    double tau_p,
                    double survival_rate
                    );

double GetMax(
			double *number_of_cells,
			long n_xslot
			);


void Compute_division_threshold(struct sim_parameters *params,
                                double *division_threshold,     //array encoding that info vs spatial location
                                double *oxygen_concentration,   //current oxygen concentration vs space
                                long right_boundary 
                                );

void Compute_equilibria(struct sim_parameters *params,
                        double *pneq,
                        double *pceq
    );


