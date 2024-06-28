
//LIBRARIES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <chrono>
#include <iostream>


#define N_OUTPUTS 100 



struct all_parameters{
    double delta_x; 
    double delta_t; 
    double tstop; 
    double D;
    double r; 
    double K; 
    double n_xslots; 


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
                          FILE *DATA_FILE,
                        
                          double **cells_concentration  //lo natural seria long y no double
                         
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



                    