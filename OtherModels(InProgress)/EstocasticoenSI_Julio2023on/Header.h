

#define positive_part(A) ((A)>0.0 ? (A) : 0.0)
#define min(A,B) ( (A) < (B) ? (A) : (B) )
#define max(A,B) ( (A) < (B) ? (B) : (A) )
#define INFINITE_TIME 1000000000000000000.0 //A time instant we will hopefully never reach
#define N_OUTPUTS 100 //We may check this from the input files instead
#define TOL 1e-10 //Para controlar Newton
#define MEMORY_BATCH 50000 //100 //memory margin
#define MEMORY_GAP 3
#define NUMBER_OF_REACTIONS 4 //1=birth, 2=death, 3= left_diff, 4=right_diff (substract 1 for vector indices)


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
        //(1.0-0.4*log(1.0/0.017*(2.02-1.0/(1.0-0.45185/(x*x)))))

//Gillespie tags
#define BIRTH 0
#define DEATH 1
#define LDIFF 2
#define RDIFF 3



//LIBRARIES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>



//STRUCTS

struct age_structure{ //one for each spatial slot
    double *age_distribution;
    long ncells;
    long size_allocated_memory;
    long ncells_ready_to_divide;//ready_to_divide_index;
    double propensity_vector[NUMBER_OF_REACTIONS];
    		/*0=birth, 1=death, 2=left diffusion, 3=right diffusion*/
};


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
    double CFL_number; //dimensionless time (check again)
};

////////////////****************************////////////////////////////////

//SUBROUTINES


///////////////////
//READ FILE

long Read_Init_Population(
                          FILE *INITIAL_POPULATION,
                          struct age_structure ***pStochasticAges,
                          long **pncells_per_slot,  //lo natural seria long y no double
                          struct sim_parameters *params,
                          double *oxygen_concentration
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
                     long *pn_files,
                     double *pDelta_x,
                    double *pDelta_t
                     );



///////////////////////
//PDE SOLVER FILE



void EvalOxygen(    struct sim_parameters *params,
                    double *koxy, //rhs readout
                    double *oxygen_concentration, //old status oxygen
                    long *number_of_cells //old status population (sum of both ones)
                    );


void  Euler_oxy(struct sim_parameters *params,
                double delta_t,
                long *number_of_cells,
                double *oxygen_concentration
                );


void GlobalPDE_Handler(struct sim_parameters *params,
                double tau,
                double *oxygen_concentration,
                long *number_of_cells
                );




//////////////////////
//UTILITIES FILE

void reverse(char s[]);

void itoa(int n,
          char s[]
          );

int cmpfunc (const void * a, const void * b);

void Print_Vector_Double(
				FILE *OUTPUT_DATA,
				char *output_path, 
				int iteration, //Si ponemos char no podemos pasar de 127 ficheros
				long vector_size, 
				double *vector
				);

void Print_Vector_Long(
                FILE *OUTPUT_DATA,
                char *output_path, 
                int iteration, //Si ponemos char no podemos pasar de 127 ficheros
                long vector_size, 
                long *vector
                );

void Print_VectorLocation(
        FILE *OUTPUT_DATA,
        char *output_path, 
        char iteration, 
        long vector_size, 
        double *vector
        );


double Get_Eigenvalue(
                    double threshold,  
                    double death_rate,
                    //double tau_p,
                    double survival_rate
                    );

double GetMax(
			double *number_of_cells,
			long n_xslot
			);

double GetMax_Long(long *number_of_cells,long n_xslot);


void Compute_division_threshold(struct sim_parameters *params,
                                double *division_threshold,     //array encoding that info vs spatial location
                                double *oxygen_concentration,   //current oxygen concentration vs space
                                long right_boundary  //rightmost slot entering the computation (aka nslots here)
);


void Check_for_extra_memory(struct age_structure **StochasticAges,
                            long n_xslots
                            );

void Compute_equilibria(struct sim_parameters *params,
                        double *pneq,
                        double *pceq
                        );

void Initialize_adimensional_ages(struct sim_parameters *params,
                                struct age_structure ***pStochasticAges
                                    );

///////////////////
//GILLESPIE FILE

void Initialize_propensities(struct sim_parameters *params,
                            struct age_structure **StochasticAges,
                            double *oxygen_concentration
                            );


double Compute_total_propensity(struct age_structure **StochasticAges,
                                long n_xslots);


void Determine_event(struct sim_parameters *params,
                        long *firing_slot,
                        int *firing_reaction,
                        long *firing_cell,
                        struct age_structure **StochasticAges,
                        double fraction_of_total_propensity
                        );


void Advance_ages(struct age_structure **StochasticAges,
                    long n_xslots,
                    double tau);


void Remove_cell(double *vector,
                long location,
                long n_elements);


void Insert_cell(double *vector,
                double age,
                long n_elements);


void Handle_event(struct sim_parameters *params,
                struct age_structure **StochasticAges,
                  long firing_slot,    
                  int firing_reaction,
                  long firing_cell,
                  long *ncells_per_slot,
                  double r3);


void Update_birth_issues(struct sim_parameters *params,
                        struct age_structure **StochasticAges,
                        double *oxygen_concentration
                        );


