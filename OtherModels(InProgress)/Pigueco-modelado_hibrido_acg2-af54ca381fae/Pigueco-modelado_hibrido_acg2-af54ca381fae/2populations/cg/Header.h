



#define positive_part(A) ((A)>0.0 ? (A) : 0.0)
#define min(A,B) ( (A) < (B) ? (A) : (B) )
#define max(A,B) ( (A) < (B) ? (B) : (A) )
#define NEGLIGIBILITY_THRESHOLD 1e-20 //Populations below this value are capped to zero
#define INFINITE_TIME 1000000000000000000.0 //A time instant we will hopefully never reach
#define N_OUTPUTS  100 //number of output files that will be produced



//Some system constants
#define c_0 1.1
#define r_cr 1.004
#define aminus 8250
#define beta 0.2
#define INIT_MASS 5.0
//#define DIFFUSION_COEF_CELLS 0.0000006 //mm^2/(6seg), taken from Anderson
//#define DIFFUSION_COEF_O2 0.006 //mm^2/(6seg), taken from Anderson (initially in cm^2/s)


//Some functions determining division rates
#define c_cr(x) (1.0-0.4*log(1.0/0.017*(2.02-1.0/(1.0-0.45185/(x*x)))))

#define First_division_threshold(a,c)    a*exp(-c/c_0)
//here a is aplus (which we must compute explicitly at the beginning of the simulation) and c is the current oxygen level

#define Second_division_threshold(x,c)  aminus*pow((c/c_cr(x))-1,-beta)
//here x is p6_over_p3 and c is the current oxygen level





#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>



////////////////****************************////////////////////////////////

//SUBROUTINES



//READ FILE

long Read_Init_Space(
                     FILE *DATA_FILE,
                     double **pinitial_population//,
                     //double **pbackup_initial_population
                     );


long Read_Init_Space_Oxygen(
                            FILE *DATA_FILE,
                            double **pinitial_oxygen
                            );


void Read_params_population(
                            FILE *PARAMETERS,
                            double *pdeath_rate_inv,
                            double *ptau_p,
                            double *pp6_over_p3,
                            double *pdiff_coef
                            );


void Read_params_sim(
                     FILE *SIM_DATA,
                     double *psource_oxygen,
                     double *pk_oxygen,
                     double *pdiff_coef,
                     double *pdelta_t,
                     double *pdelta_x,
                     long *psampling_period
                     );


//PDE SOLVER FILE


void EvalPopulation(
                    double *kpop, //rhs readout
                    double *N_previous, //old status population
                    double death_rate_inv,
                    double *threshold,
                    double tau_p,
                    double diff_coef,
                    long x_size
                    );

void EvalOxygen(
                    double *koxy, //rhs readout
                    double *oxygen_level, //old status oxygen
                    long xsize, 
                    double *number_of_cells, //old status population (sum of both ones)
                    double k_oxygen, 
                    double source_oxygen,
                    double diff_coef
                    );

double Refine_time_step(
          double tau,   
          double Diff_Coef_Pop_host, 
          double Diff_Coef_Pop_invader,
          double Diff_Coef_O2, 
          double k_oxygen, 
          double peak_number_of_cells
          );

void  RK4pop_oxy(
                double delta_t,
                long n_xslots,
                double *number_of_cells_host,
                double *number_of_cells_invader,
                double *total_number_of_cells,
                double death_rate_inv_host,
                double death_rate_inv_invader,
                double *division_threshold_host,
                double *division_threshold_invader,
                double tau_p_host,
                double tau_p_invader,
                double aplus_host,
                double aplus_invader,
                double p6_over_p3_host,
                double p6_over_p3_invader,
                double diff_coef_host,
                double diff_coef_invader,
                double diff_coef_oxygen,
                double *oxygen_level,
                double k_oxygen,
                double source_oxygen
                );


void RK4Handler( 
                double delta_t,
                long n_xslots,
                //double delta_x,
                double *number_of_cells_host,
                double *number_of_cells_invader,
                double death_rate_inv_host,
                double death_rate_inv_invader,
                double *division_threshold_host,
                double *division_threshold_invader,
                double tau_p_host,
                double tau_p_invader,
                double aplus_host,
                double aplus_invader,
                double p6_over_p3_host,
                double p6_over_p3_invader,
                double diff_coef_host,
                double diff_coef_invader,
                double diff_coef_oxygen,
                double *oxygen_level,
                double k_oxygen,
                double source_oxygen
                );
/*
void Finite_Difference_Solver(
                              double delta_t,
                              double *N_previous,
                              double death_rate_inv,
                              double *threshold,
                              double tau_p,
                              long x_size,
                              double delta_x,
                              double diffcoef
                              );


void OxygenSolverDamp(double *oxygen_level,
                      long xsize,
                      double *number_of_cells,
                      double delta_t,
                      double k_oxygen,
                      double source_oxygen,
                      double peak_number_of_cells,
                      double delta_x,
                      double diffcoef
                      );
*/


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