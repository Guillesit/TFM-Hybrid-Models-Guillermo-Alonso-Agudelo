
//LIBRARIES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <random>
#include <ctime>
#include <chrono>


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

struct age_structure{ //one for each spatial slot
    double *age_distribution;
    long ncells;
    long size_allocated_memory;
    long ncells_ready_to_divide;//ready_to_divide_index;
    double propensity_vector[NUMBER_OF_REACTIONS];
    		/*0=birth, 1=death, 2=left diffusion, 3=right diffusion*/
};


struct all_parameters{
    double delta_x; 
    double delta_t; 
    double tstop; 
    double D;
    double r; 
    double K; 
    double n_xslots; 
    double death_rate_hat; 

}; 


void Read_sim_params(struct all_parameters *params,
                            FILE *PARAMETERS
                            ); 

 void Read_eq_params(struct all_parameters *params,
                            FILE *PARAMETERS
                            ); 
long Read_Init_Space_Oxygen(
                            FILE *DATA_FILE,
                            double **pinitial_oxygen
                            );
long Read_Init_Population(
                          FILE *INITIAL_POPULATION,
                          struct age_structure ***pStochasticAges,
                          double **pncells_per_slot,  //lo natural seria long y no double
                          struct all_parameters *params,
                          double *oxygen_concentration
                          );

void reverse(char s[]);

void itoa(int n,
          char s[]
          );
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
                
void euler(struct all_parameters *params, 
                double *concentration
                    );
void eval_c(struct all_parameters *params, double *concentration, double *rd);

void PDE_Handler( struct all_parameters *params,
                double *concentration); 


//Gillespie

void Initialize_propensities(struct all_parameters *params,
                            struct age_structure **StochasticAges,
                            double *oxygen_concentration
                            );


double Compute_total_propensity(struct age_structure **StochasticAges,
                                long n_xslots);


void Determine_event(struct all_parameters *params,
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


void Handle_event(struct all_parameters *params,
                struct age_structure **StochasticAges,
                  long firing_slot,    
                  int firing_reaction,
                  long firing_cell,
                  double *ncells_per_slot,
                  double r3);


void Update_birth_issues(struct all_parameters *params,
                        struct age_structure **StochasticAges,
                        double *oxygen_concentration
                        );