

//Handy macros

#define positive_part(A) ((A)>0.0 ? (A) : 0.0)
#define min(A,B) ( (A) < (B) ? (A) : (B) )


//General-purpose parameters

#define INFINITE_TIME 1000000000000000000.0 //A time instant we will hopefully never reach
#define POPULATION_THRESHOLD 2000 //To tell between mean field and stochastic regimes
#define TAM 100000 //memory issues
//NOTES: 1) It is probably more useful to get the threshold from the user or from some script (06/09/2016), 
//2) Dynamic memory extension DOES NOT WORK so far (06/09/2016). It is advisable then to use a high initial memory size
// 3) INFINITE_TIME is used specifically in connection with age division thresholds. TO COMPLETE (06/09/2016)
#define N_OUTPUTS  100 //number of output files that will be produced
 




//Some model constants

#define c_0 1.1
#define r_cr 1.004
#define aminus 8250
#define beta 0.2
#define INIT_MASS 5.0
#define DIFFUSION_COEF_CELLS 0.0000006 //mm^2/(6seg), taken from Anderson
#define DIFFUSION_COEF_O2 0.006 //mm^2/(6seg), taken from Anderson (initially in cm^2/s)
//NOTES: 1) check notations at WHICH? reference document. TO COMPLETE (06/09/2016)
//2) probably diffusion coefficients should be read by means of some script, for the sake of comparisions (06/09/2016)




//Some functions determining division rates (used for the intracellular part)

#define c_cr(x) (1.0-0.4*log(1.0/0.017*(2.02-1.0/(1.0-0.45185/(x*x)))))

#define First_division_threshold(a,c)    a*exp(-c/c_0)
//here a is aplus (which we must compute explicitly at the beginning of the simulation) and c is the current oxygen level,

#define Second_division_threshold(x,c)  aminus*pow((c/c_cr(x))-1,-beta)
//here x is p6_over_p3 and c is the current oxygen level


//Constants for the random number generator
//What ranum.h used to have (we moved the function ran2 to the utilities file and we link no more to ranum.h)

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
//NOTES: this is done in that way to avoid compilation errors. Roberto used to have a different 
// random number generator, but so far we do not know how to have both of them operating at different 
//parts of the code in a compatible way (06/09/2016)



/**********************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>




////////////////****************************////////////////////////////////

//SUBROUTINES



//READ FILE

//master routine to read the input
long Read_Init_Space(				//returns the number of spatial slots
	FILE *DATA_FILE,      			//Initial population and oxygen values (user-provided)
	FILE *PARAM_FILE, 				//Values of some relevant parameters (user-provided) --those below.
	double *pdeath_rate_inv,
	double *ptau_p,
	double *psource_oxygen,			//source term in oxygen's evolution equation
	double *pk_oxygen,				//consumption term in oxygen's evolution equation
	double **pinitial_oxygen,		//initial oxygen array
	double *pp6_over_p3, 			//p6/p3 trait (the same for the whole population, not changing during evolution)
	double **pinitial_population  //initial population array
	);
//NOTES: 1) the use of "death_rate_inv" is extremely confusing at some parts of the code,
//we'd better change names and get a consistent notation (06/09/2016) YET TO BE DONE




//SOLVER FILE
//All these are routines to simulate the mean field module

double Refine_time_step(		//returns adapted time step fulfilling stability criteria
	double tau, 
	double Diff_Coef_Pop,		//normalized diffusion coefficient (population)
	 double Diff_Coef_O2, 		//normalized diffusion coefficient (oxygen)
	 double k_oxygen, 			 //consumption term in oxygen's evolution equation
	 double peak_number_of_cells //#cells at the most populated slot
	 );

void PopulationSolver(
	double delta_t, 
	double *N_previous, 
	double *threshold, 
	double death_rate_inv, 
	double tau_p, 
	long x_size, 
	double Diff_Coef
	);

void OxygenSolver(
	double delta_t,
	double *oxygen_level, 
	double *number_of_cells,  
	double k_oxygen, 			//consumption term in oxygen's evolution equation
	double source_oxygen, 		//source term in oxygen's evolution equation
	double Diff_Coef, 
	long xsize
	);

void OxygenSolverDamp(
	double delta_t,
	double *oxygen_level, 
	double *number_of_cells,  
	double k_oxygen, 			//consumption term in oxygen's evolution equation
	double source_oxygen, 		//source term in oxygen's evolution equation
	//double damping_rate,  //CURRENTLY (11/11/16) we use a very specific damping rate for which we do not need to ask
	double Diff_coef, 
	long xsize);

void Global_Mean_Field_Handler(
	double tau,  
	double *number_of_cells, 
	double *density_of_cells, 
	double *division_threshold,  
	double *oxygen_level, 
	double Diff_Coef_Pop, 
	double Diff_Coef_O2, 
	double k_oxygen,			//consumption term in oxygen's evolution equation
	double source_oxygen,		//source term in oxygen's evolution equation
	double delta_x, 
	double tau_p, 
	double death_rate_inv,
	double p6_over_p3,
	long Interface_location, 
	long n_xslots
	);




//INTRACELLULAR FILE
//homemade way of implementing the generalization of Bedessem-Stephanou's module

void Book_Vector(
	double **pn, 	//brand new dynamic array
	long S 			//array size
	);

double Logistic(		//returns evaluation of a specific logistic function
	double time, 
	double init,
	double alpha, 
	double m_0
	);

void Eval_F(
	double *outcome,
	double *spatial_arguments,	//current unknown values
	double *constants, 		//ODE system constant's
	double time 		//current time instant
	);

void RK4(
	double *previous_status,
	double *updated_status,
	double *constants,	//ODE system constant's
	double t,			//current time instant
	double delta_t, 	//desired time step
	int n_eqs			//number of equations
	);

double Compute_aplus(double p6_over_p3);
//NOTES: 1) Roberto has a different implementation for this part, using a god-given RK4 method
//2) Notations for this part can be obtained in WHICH? reference document. TO COMPLETE (06/09/2016)




//UTILITIES FILE
//homemade routines to perform miscellaneous operations

double Cummulative_Cells(	//returns total number of cells in the system? YET TO SCAN PROPERLY (06/09/2016)
	double *number_of_cells, 
	long age_slots, 
	double **Actual_population, 
	long n_xslot
	);


void Handle_memory_extension( //IT DOES NOT WORK (06/09/2016)
	double **Actual_population,
	long *page_slots, 
	long space_slots
	);


//Reverses a string
//This is an example in Kernigan-Ritchie's book
void reverse(char s[]);


//Converts an integer n into a string of characters s
//This is an example in Kernigan-Ritchie's book
void itoa(int n,
          char s[]
          );


void Print_Vector(
    FILE *OUTPUT_DATA,  //where to print the stuff
    char *output_path,  //where is that file located
    int iteration,      //current iteration of the hybrid method
    long vector_size,   
    double *vector
    );



double Estimate_Limit_O2( //returns mean field prediction for the limit oxygen value as time grows big
	double death_rate_inv,
	double tau_p,
	double p6_over_p3
	);


double Get_Eigenvalue( //returns the dominant eigenvalue for the proliferation term after coarsening
	double lower_limit, 
	double death_rate_inv,
	double tau_p
	);


double GetMax( //returns the maximum of an array (crapy implementation)
	double *number_of_cells,	//array under scrutiny
	long n_xslot 	//rightmost slot
	);


void Compute_division_threshold(
    double *division_threshold,     //array encoding that info vs spatial location
    double *oxygen_level,   //current oxygen concentration vs space
    double p6_over_p3,
    long right_boundary  //rightmost slot entering the computation
    );


float ran2(long *idum); //random number generator at ranum.h





//ROBERTO FILE
//Clumsy adaptation of Gillespie's module as coded by Roberto to the general hybrid structure

double reactions(
	long first, 
	long end, 
	double **gilAux, 
	long *J, 
	double **age, 
	double **b, 
	double nu, 
	double delta_x
	);

void gillespie(double tau, 
                long first,    
                long end, 
                double **gilAux, 
                double ** ageAux,
                double **age, 
                long *J, 
                double **b, 
                double **ag, 
                double suma, 
                long *pseed, 
                double *number_of_cells, 
                double *division_threshold, 
                double *oxygen_level, 
                double p6_over_p3, 
                double tau_p, 
                long *ptam
                );
//NOTES: 
//2) It would be HIGHLY desirable to have some self-explanatory names for the variables here





//HYBRID HANDLING FILE

void Transfer_particle(
	long location,
	long *pseed,
	double ag1s,
	double death_rate_inv, 
	double tau_p, 
	double **age, 
	long *J, 
	double **b, 
	double **ag, 
	double **gilAux, 
	double **ageAux, 
	long *ptam
	);

void Transfer_particle2(
	long location,
	long *pseed,
	double ag1s,
	double death_rate_inv,
	double *number_of_cells,
	double **age,
	long *J, 
	double **b, 
	double **ag, 
	double **gilAux, 
	double **ageAux, 
	long *ptam
	);


long Realocate_interface(
	double *number_of_cells, 
	long n_xslots
	);

void Renormalize_center(
	double *number_of_cells,  
	long Interface_location, 
	long *pseed,
	double ag1s,
	double death_rate_inv, 
	double tau_p, 
	double **age, 
	long *J, 
	double **b, 
	double **ag,
	double **gilAux, 
	double **ageAux, 
	long *ptam,
	double newAux
	);

void R_l(
	long position,
	double *number_of_cells,
	long *pseed,
	double ag1s,
	double death_rate_inv, 
	double tau_p, 
	double **age, 
	long *J, 
	double **b, 
	double **ag,
	double **gilAux, 
	double **ageAux, 
	long *ptam,
	double newAux
	);

void Renormalize_left(
	long Interface_location,
	long aux_interface,
	double *number_of_cells,
	long *pseed,
	double *division_threshold,
	double death_rate_inv, 
	double tau_p, 
	double **age, 
	long *J,
	double **b, 
	double **ag,
	double **gilAux, 
	double **ageAux, 
	long *ptam,
	double newAux
	);



