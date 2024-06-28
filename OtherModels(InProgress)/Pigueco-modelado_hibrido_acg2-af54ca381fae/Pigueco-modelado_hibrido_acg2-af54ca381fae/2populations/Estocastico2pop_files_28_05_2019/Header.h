



#define positive_part(A) ((A)>0.0 ? (A) : 0.0)
#define min(A,B) ( (A) < (B) ? (A) : (B) )
#define max(A,B) ( (A) < (B) ? (B) : (A) )
#define NEGLIGIBILITY_THRESHOLD 1e-20 //Populations below this value are capped to zero
#define INFINITE_TIME 1000000000000000000.0 //A time instant we will hopefully never reach
#define N_OUTPUTS  100 //number of output files that will be produced
//Of lately we read this from the input...
#define MEMORY_BATCH 50000 //number of extra memory locations that will be requested 
  //in a memory extension round (aka TAM in earlier versions)



//Some system constants
#define c_0 1.1
#define r_cr 1.004
#define aminus 8250
#define beta 0.2
#define INIT_MASS 5.0
#define DIFFUSION_COEF_CELLS 0.0000006 //mm^2/(6seg), taken from Anderson
#define DIFFUSION_COEF_O2 0.006 //mm^2/(6seg), taken from Anderson (initially in cm^2/s)
#define NUMBER_OF_REACTIONS 4 //birth, death, left and right diffusion

//Some functions determining division rates
#define c_cr(x) (1.0-0.4*log(1.0/0.017*(2.02-1.0/(1.0-0.45185/(x*x)))))
//(see formula (3) in the JCopmPhys paper)

#define First_division_threshold(a,c)    a*exp(-c/c_0)
//here a is aplus (which we must compute explicitly at the beginning of the simulation) 
//and c is the current oxygen level (see formula (2) in the JCopmPhys paper)

#define Second_division_threshold(x,c)  aminus*pow((c/c_cr(x))-1,-beta)
//here x is p6_over_p3 and c is the current oxygen level
//(see formula (2) in the JCopmPhys paper)




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>



////////////////****************************////////////////////////////////

//SUBROUTINES



//READ FILE

long Read_Init_Population(
                          FILE *INITIAL_POPULATION,
                          long ***pStochasticNcells,
                          double ***pStochasticAge,
                          long **pnspaces
                          );

long Read_Init_Space_Oxygen(
                            FILE *DATA_FILE,
                            double **pinitial_oxygen
                            );


void Read_params_population(
                            FILE *PARAMETERS,
                            double *pdeath_rate,
                            double *ptau_p,
                            double *pp6_over_p3
                            );


void Read_params_sim(
                     FILE *SIM_DATA,
                     double *psource_oxygen,
                     double *pconsumption_oxygen,
                     double *pdelta_t,
                     double *pdelta_x,
                     double *ptstop,
                     int *pn_files,
                     double *psurvival_rate
                     );



//PDE SOLVER FILE


void OxygenSolverDamp(
                    double *oxygen_level,
                    long xsize, 
                    double *number_of_cells, 
                    double delta_t, 
                    double k_oxygen, 
                    double source_oxygen,
                    // double peak_number_of_cells,
                    double delta_x
                    );


//BIRTH FILE



double Reactions(int n_spatial_slots,
                 double **gilAux,
                 long *J,
                 long **StochasticNcells,
                 double **b,
                 double death_rate,
                 double diffusion_R, 
                 double diffusion_L
                 );


//Includes therapy (effective modification of birth-death rates,
//instead of Alarcon's a posteriori rule)
double Reactions_two_pop(int n_spatial_slots,
                         double **gilAux,
                         long *J_host,
                         long **StochasticNcells_host,
                         double **b_host,
                         long *J_inv,
                         long **StochasticNcells_inv,
                         double **b_inv, 
                         double death_rate_host,
                         double death_rate_invader,
                         double diffusion_R, 
                         double diffusion_L,
                         double survival_rate
                         );


void Gillespie(int i,
               int j,
               long **StochasticNcells,
               long *J,
               double **StochasticAge,
               double **b
               );




double birth(double ox);






//INTRACELLULAR FILE

void Book_Vector(
        double **pn, 
        long S
        );

double Logistic(
        double time, 
        double init,
        double alpha, 
        double m_0
        );

void Eval_F(
      double *outcome,
      double *spatial_arguments,
      double *constants,
      double time
      );

void RK4(
    double *previous_status,
    double *updated_status,
    double *constants,
    double t,
    double delta_t, 
    int n_eqs
    );

double Compute_aplus(double p6_over_p3);




//UTILITIES FILE

double Compute_density(
                     long n_spatial_slots,
                     long *J,
                     long **StochasticNcells,
                     double *densidad
                     );

double Get_division_threshold(double oxygen_level,
                              double p6_over_p3,
                              double aplus
                              );


double Cummulative_Cells(
						double *number_of_cells, 
						long age_slots, 
						double **Actual_population, 
						long n_xslot
						);


void reverse(char s[]);

void itoa(int n,
          char s[]
          );

void Print_Vector(
				FILE *OUTPUT_DATA,
				char *output_path, 
				char iteration, 
				long vector_size, 
				double *vector
				);

double Estimate_Limit_O2(
						double death_rate_inv,
						double tau_p,
						double p6_over_p3
						);

double Get_Eigenvalue(
					double *threshold, 
					long xslot, 
					double death_rate_inv,
					double tau_p
					);

double GetMax(
			double *number_of_cells,
			long n_xslot
			);


void Compute_division_threshold(
                                double *division_threshold,     //array encoding that info vs spatial location
                                double *oxygen_level,   //current oxygen concentration vs space
                                double p6_over_p3,
                                double aplus,
                                long right_boundary  //rightmost slot entering the computation (see usage below)
);
