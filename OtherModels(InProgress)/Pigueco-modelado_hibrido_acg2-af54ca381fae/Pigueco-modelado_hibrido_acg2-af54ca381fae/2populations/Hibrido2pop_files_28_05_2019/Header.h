



#define positive_part(A) ((A)>0.0 ? (A) : 0.0)
#define min(A,B) ( (A) < (B) ? (A) : (B) )
#define max(A,B) ( (A) < (B) ? (B) : (A) )
#define POPULATION_THRESHOLD 500 //To tell between mean field and stochastic regimes
#define NEGLIGIBILITY_THRESHOLD 1e-20 //Populations below this value are capped to zero
#define INFINITE_TIME 1000000000000000000.0 //A time instant we will hopefully never reach
#define N_OUTPUTS  100 //number of output files that will be produced
//NOTE THAT...
//Of lately we read this from the input...
#define MEMORY_BATCH 100000 //number of extra memory locations that will be requested 
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

// void EvalPopulation(
//                     double *kpop, //rhs readout
//                     double *N_previous, //old status population
//                     double death_rate,
//                     double *threshold,
//                     double tau_p,
//                     double diff_coef,
//                     long init_index,
//                     long end_index
//                     );

void EvalPopulation_targeted(
                    double *kpop, //rhs readout
                    double *N_previous, //old status population
                    double death_rate,
                    double *threshold,
                    double tau_p,
                    double diff_coef,
                    long init_index,
                    long end_index,
                    double survival_rate
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


double Refine_time_step_pop(
                    double tau,     
                    double Diff_Coef
                    );


double Refine_time_step_oxy(
          double tau,   
          double Diff_Coef_O2, 
          double k_oxygen, 
          double peak_number_of_cells
          );


void  RK4pop(
                double delta_t,
                long n_xslots,
                long init_index,
                long end_index,
                double *density_of_cells,
                double death_rate,
                double *division_threshold,
                double tau_p,
                double diff_coef,
                double survival_rate
                );


void  RK4oxy(
                double delta_t,
                long n_xslots,
                double *total_number_of_cells,
                double diff_coef_oxygen,
                double *oxygen_level,
                double k_oxygen,
                double source_oxygen
                );


void RK4HandlerPop( 
                double delta_t,
                long n_xslots,
                long init_index,
                long end_index,
                double *density_of_cells,
                double death_rate,
                double *division_threshold,
                double tau_p,
                double diff_coef,
                double survival_rate
                );


void RK4HandlerOxy(    
                double delta_t,
                long n_xslots,
                double *number_of_cells_host, //DE HECHO TENEMOS EL TOTAL NUMBER OF CELLS EN LA RUTINA PRINCIPAL
                double *number_of_cells_invader,
                double diff_coef_oxygen,
                double *oxygen_level,
                double k_oxygen,
                double source_oxygen
                );




//BIRTH FILE

double Reactions_invader( //sum of propensities for the population at the left
                        int n_spatial_slots,
                        double **Propensities_inv,
                        long *J_inv,
                        long **StochasticNcells_inv,
                        double **b_inv,
                        long interface_invader,
                        double death_rate_invader,
                        double diffusion_invader,
                        double survival_rate
                        );



double Reactions_host( //sum of propensities for the population at the right
                        int n_spatial_slots,
                        double **Propensities_host,
                        long *J_host,
                        long **StochasticNcells_host,
                        double **b_host,
                        long interface_host,
                        double death_rate_host,
                        double diffusion_host,
                        double survival_rate
                        );


double Reactions_two_pop(
        int n_spatial_slots,
        double **gilAux,  
        long *J_host,
        long **StochasticNcells_host,
        double **b_host,
        long *J_inv,
        long **StochasticNcells_inv,
        double **b_inv, 
        long interface_host,
        long interface_invader,
        double death_rate_host,
        double death_rate_invader,
        double diffusion_host, 
        double diffusion_invader,
        double survival_rate //associated with a therapy, affects both populations
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

double Compute_ncells( //Returns total number of cells, although we will not always use that value in the main body
                     long n_spatial_slots,
                     long *J,
                     long **StochasticNcells,
                     double *ncells
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
				int iteration, 
				long vector_size, 
				double *vector
				);

double Estimate_Limit_O2(
						double death_rate_inv,
						double tau_p,
						double p6_over_p3
						);

// double Get_Eigenvalue(
// 					double *threshold, 
// 					long xslot, 
// 					double death_rate_inv,
// 					double tau_p
// 					);

double Get_Eigenvalue_survival_rate(
                    double *threshold, 
                    long xslot, 
                    double death_rate,
                    double tau_p,
                    double survival_rate
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





//HYBRID FILE

double Sample_age_distribution(
                                double ag1s, //scalar, nor vector
                                double eigenvalue,
                                double nu,
                                double taup,
                                double diffusion_coefficient,
                                char *pflag_above_threshold
                                );

void Transfer_particle(
                        long Interface_location,
                        double *ag1s,
                        double death_rate, 
                        double tau_p, 
                        long **StochasticNcells, 
                        long *J, 
                        double **b, 
                        double **StochasticAge, 
                        double survival_rate,
                        double diffusion_coef
                        );

void Eliminate_particle(
                        long location,
                        long **StochasticNcells, 
                        long *J, 
                        double **b, 
                        double **StochasticAge 
                        );




void Renormalize_center(
                        double *number_of_cells,  
                        long Interface_location, 
                        long init_index,
                        long end_index,
                        double *ag1s,
                        double death_rate, 
                        double tau_p, 
                        long **StochasticNcells, 
                        long *J, 
                        double **b, 
                        double **StochasticAge,
                        double survival_rate,
                        double diffusion_coef,
                        double cells_at_interface
                        );


void Renormalize_inwards(
                        long n_spatial_slots,
                        long Interface_location,
                        long aux_interface,
                        double *number_of_cells,
                        double *division_threshold,
                        double death_rate, 
                        double tau_p, 
                        long **StochasticNcells, 
                        long *J, 
                        double **b, 
                        double **StochasticAge,
                        double survival_rate,
                        double diffusion_coef
                        );

