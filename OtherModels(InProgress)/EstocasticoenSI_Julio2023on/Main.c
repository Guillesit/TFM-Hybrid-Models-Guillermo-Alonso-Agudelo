

//One spatial dimension integrating through age structure

#include "Header.h"




int main(int argc, char **argv)
{
    FILE *INITIAL_POPULATION;
    FILE *INITIAL_OXYGEN;
    FILE *SIM_DATA;
    FILE *POPULATION;
    FILE *PARAMETERS;
    FILE *OUTPUT_DATA;
    FILE *MESSAGES;
    
    
    long i, x_slot, n_xslots=0;  //counters and tracers
    
    long n_files;
    double sampling_time_window = 0.0;
    long counter=1; 
    
    //SI units: meters and seconds

    double t=0.0,  tau=0.0; //Internal units
    double tstop=0.0; //SI units

    struct age_structure **StochasticAges; //SI units at input, then converted to internal units
    
    //depending on time and space
    //double *oxygen_level; //SI units, just for input-output
    double *oxygen_concentration; //Internal units
    long *ncells_per_slot; //number of cells at each spatial slot, regardless of their age
    //At the only place where ew need densities of cells
    //We will rescale to get an effective rate constant
    long total_number_of_cells; //overall number of cells in the system

    struct sim_parameters parameters; //SI units at input, then converted to internal units
    double Delta_x; //grid size (SI units)
    double Delta_t; //Suggested time increment (SI units)
    double n_eq=-1.0;
    double c_eq=-1.0;

	//Stuff to handle output files
    int output_iter=0;
    
    char output_path1[200]="OutputValuesPopulation";
    char output_path2[200]="OutputValuesOxygen";
    char usertag[100]="";
    char label1[200]="OutputValuesPopulation";
    char label2[200]="OutputValuesOxygen";
    char label3[200]="messages";
    char label4[200]="total_population_vs_time";
    char create_folder[300]="mkdir ";

    int hour; //to seed the random number generator
    double r1, r2, r3; //random numbers

    double total_propensity; //Internal units
    long firing_slot=-1;
    long firing_cell=-1;
    int firing_reaction=0;//0=birth, 1=death, 2= left_diff, 3=right_diff

    
    /************************/
    
    /*****************/
    //OPENING FILES
    
    if(argc!=6){
		printf("Error, incorrect argument number\n");
		printf("Use as: \n exec_name  \n initial_population_file \n initial_oxygen_file \n population_parameters_file \n simulation_parameters_file \n output_files_tag \n");
		return(1);
	};
    
    if((INITIAL_POPULATION=fopen(argv[1],"rt"))==NULL){
		printf("Error: initial_population_file could not be opened \n");
		return(1) ;
	};
    
    if((INITIAL_OXYGEN=fopen(argv[2],"rt"))==NULL){
		printf("Error: initial_oxygen_file could not be opened \n");
		return(1) ;
	};
    
    if((PARAMETERS=fopen(argv[3],"rt"))==NULL){
		printf("Error: population_parameters_file could not be opened \n");
		return(1) ;
	};
    
    if((SIM_DATA=fopen(argv[4],"rt"))==NULL){
		printf("Error: simulation_parameters_file could not be opened \n");
		return(1) ;
	};
    
    strcpy(usertag,argv[5]);
    strcat(label4,usertag);
    
    if((POPULATION=fopen(label4,"w"))==NULL){
		printf("Error: output_population_file could not be opened \n");
		return(1) ;
	};
  
   
    /*********************/
    //GATHERING ADDITIONAL INITIAL INFORMATION AND INITIALIZING

    //Read population parameters 
    Read_params_population(&parameters,
                            PARAMETERS);

    //Read parameters (simulation infos)
    Read_params_sim(&parameters,
                    SIM_DATA,
                    &tstop,
                    &n_files,
                    &Delta_x,
                    &Delta_t
                    ); 
    
    /*******************************/
    //PRINTING SIMULATION INFOS

    //Check again: 
    sampling_time_window=tstop/n_files;
  
    //create "messages" (plus user tag)
    strcat(label3,usertag);
    if((MESSAGES=fopen(label3,"w"))==NULL){ //Master output file
		printf("Error: output messages file could not be created \n");
		return(1) ;
	};
    
    fprintf(MESSAGES,"#Debug info: \n \n");
    fprintf(MESSAGES,"#Data generated with the following .exe file:%s\n",argv[0]);
    fprintf(MESSAGES,"#Initial condition file (population ages):%s\n",argv[1]);
    fprintf(MESSAGES,"#Initial condition file (oxygen):%s\n",argv[2]);
    fprintf(MESSAGES,"#Population parameters file:%s\n",argv[3]);
    fprintf(MESSAGES,"#Parameters file (simulation):%s\n",argv[4]);

    fprintf(MESSAGES,"#Default unit system as given by the S.I. (seconds, meters, kilograms),\n#let X for space and T for time\n");
    fprintf(MESSAGES, "#Grid step (aka delta_x): %lf X \n",Delta_x);
    fprintf(MESSAGES, "#Suggested time increment (aka delta_t): %lf T \n",Delta_t);
    fprintf(MESSAGES,"#Sampling time window: %lf \n",sampling_time_window);
    fprintf(MESSAGES,"#User-provided final time: %lf \n",tstop);
    fprintf(MESSAGES,"#Number of output files: %ld \n",n_files);
    
    fprintf(MESSAGES,"#survival rate: %lf \n",(parameters.survival_rate));
    fprintf(MESSAGES,"#Death rate: %.10lf T^-1 \n",parameters.death_rate_hat);
    fprintf(MESSAGES,"#tau_p: %.10lf T\n",parameters.tau_p);
    fprintf(MESSAGES,"#aplus: %.10lf T\n",parameters.aplus);
    fprintf(MESSAGES,"#p6/p3: %.10lf \n",parameters.p6_over_p3);
    fprintf(MESSAGES,"#Diffusion coefficients (pop, oxy): %.10lf %lf X^2/T\n",parameters.diff_coef_pop,parameters.diff_coef_oxygen);
    fprintf(MESSAGES,"#k_decay: %.10lf T^-1\n",parameters.k_decay_hat);
    fprintf(MESSAGES,"#k_consumption: %.10lf X T^-1\n",parameters.k_consumption_hat);
    fprintf(MESSAGES,"#k_source: %.10lf X^-1 T^-1\n",parameters.source_oxygen_hat);
 
     
    Compute_equilibria(&parameters,&n_eq,&c_eq); //En el sistema de unidades primitivo
    //cuidado que son valores de densidades
    //El nº individuos no lo vamos a normalizar

    fprintf(MESSAGES,"#neq=%.10lf X^-1, ceq=%.10lf X^-1\n",n_eq,c_eq);
    fprintf(MESSAGES,"#Note: spatial outputs are dimensionless");

    //////////////////
    //ADIMENSIONALIZE 
    //Esto de aqui hasta el final de la simulacion, incluso para las salidas

    ////////////////
    //Rescale constants to the internal unit system (ADIMENSIONAL FROM NOW ON)
    parameters.death_rate_hat=(parameters.death_rate_hat)*(parameters.tau_p);
    parameters.aminus_hat=aminus/(parameters.tau_p);
    parameters.diff_coefs_ratio=(parameters.diff_coef_oxygen)/(parameters.diff_coef_pop);
    parameters.k_decay_hat=(parameters.k_decay_hat)*(parameters.tau_p);

    parameters.k_consumption_hat=(parameters.k_consumption_hat)*(parameters.tau_p)*Delta_x;//*n_eq;
    //Please note that we should multiply by n_eq and we are not to involve Delta_x.
    // However, since we shall use the number of individuals to mimic the population concentration, 
    //the effective rate constant would be obtained by multiplication with Delta_ x and dividing by n_eq afterwards 
    //Hence we would have a cancellation and therefore we do nothing. 

    parameters.source_oxygen_hat=(parameters.source_oxygen_hat)*(parameters.tau_p)/c_eq;
    parameters.typical_lengthscale=sqrt((parameters.diff_coef_pop)*(parameters.tau_p));

    parameters.critical_oxy_hat=(typical_oxy*severe_hypoxia)*(c_cr(parameters.p6_over_p3)/c_cr(1))/c_eq;

    parameters.Delta_x_hat=Delta_x/(parameters.typical_lengthscale);
    parameters.Delta_t_hat=Delta_t/(parameters.tau_p);

    //Effective diffusion coefficients to take into account second order finite differences:
    parameters.diff_coef_eff=1.0/((parameters.Delta_x_hat)*(parameters.Delta_x_hat)); //diff pop
    parameters.diff_coefs_ratio=parameters.diff_coefs_ratio/((parameters.Delta_x_hat)*(parameters.Delta_x_hat));
  
    //////Easier to determine CFL condition here (at least the part that does not depend on the iteration)
    parameters.CFL_number=min(1.0/(6*(parameters.diff_coefs_ratio)),
                                min(1.0/(1.0+(parameters.k_decay_hat)),1.0/(1.0+(parameters.death_rate_hat)))     
                            );

    /////////////////////
    //Reading and initializing the stochastic age distribution+oxygen

    //Read oxygen
    n_xslots=Read_Init_Space_Oxygen(INITIAL_OXYGEN,&oxygen_concentration);

    for (i = 0; i < n_xslots;i++){oxygen_concentration[i]=oxygen_concentration[i]/c_eq;};
  
    //Read initial population 
    //Generate adimensional age equilibrium distribution while we are at it

    //POR PROBAR DE FORMA EXHAUSTIVA
    if (n_xslots!=Read_Init_Population(INITIAL_POPULATION,
                                &StochasticAges,
                                &ncells_per_slot,
                                &parameters,
                                oxygen_concentration)
        ){
        fprintf(stderr,"Mismatched number of spatial cells, aborting execution\n");
        exit(1);
    };
   
    parameters.n_xslots=n_xslots;

    //Ordering from bigger to smaller ages within each spatial cell
    for (i = 0; i < n_xslots; i++){
        //printf("%ld\n",ncells_per_slot[i]); //DEBUG BUSSINESS
        qsort((StochasticAges[i]->age_distribution), ncells_per_slot[i], sizeof(double), cmpfunc);
    };
    

    /*******************************/
    //INITIALIZING and printing the rest of the stuff

     //create "OutputValuesPopulation" and "OutputValuesOxygen" (plus user tag)
     strcat(label1,usertag);
     strcat(label2,usertag);
     
     strcat(create_folder,label1);
     system(create_folder);
     strcpy(create_folder,"mkdir ");
     strcat(create_folder,label2);
     system(create_folder);
     
     
     //Printing spatial number of cells regardless of their age
     strcat(label1,"/Out"); //Final label for population folder
     strcpy(output_path1,label1);
     Print_Vector_Long(OUTPUT_DATA,output_path1,output_iter,n_xslots,ncells_per_slot);
   
     //Printing oxygen spatial distribution
     strcat(label2,"/Out"); //Final label for oxygen folder
     strcpy(output_path2,label2);
     Print_Vector_Double(OUTPUT_DATA,output_path2,output_iter,n_xslots,oxygen_concentration);

    
    //computing and printing total number of cells 
    total_number_of_cells=0.0;
    for(x_slot=0;x_slot<n_xslots;x_slot++){
        total_number_of_cells=total_number_of_cells+ncells_per_slot[x_slot];
    };
    fprintf(POPULATION,"%lf %ld \n",t,total_number_of_cells);

    //master handler   
    Initialize_propensities(&parameters,
                            StochasticAges,
                            oxygen_concentration
                            );

    //seeding the random number generator
    hour = time(NULL);
    srand48(hour);

    /*****************/



    /******************************/
    
    
    /******************************/
    //MAIN LOOP
    //here t is measured in internal time units
   
    while(t<tstop){

        //Hemos movido esta sentencia del final al principio del bucle
        //Pero sigue dando segmentation fault    
        Check_for_extra_memory(StochasticAges,n_xslots);

        r1=drand48(); 
        r2=drand48();

        total_propensity=Compute_total_propensity(StochasticAges,n_xslots);
        //Here we check for zero total population
        tau=log(1/r1)/total_propensity;
       
        //tau=1.0;
        //print the last state prior surpassing the time mark
        if(t+tau>=counter*sampling_time_window){
            //CUIDADO: hay que ver cuantos saltos de counter pegamos
            //En simulaciones con poca poblacion damos varios saltos con un solo tau
            //y la salida no es fiable
            //POR MEJORAR (y modularizar)
            
            //Printing spatial number of cells
            strcpy(output_path1,label1);
            Print_Vector_Long(OUTPUT_DATA,output_path1,counter,n_xslots,ncells_per_slot);
            
            //Printing oxygen spatial distribution
            strcpy(output_path2,label2);
            Print_Vector_Double(OUTPUT_DATA,output_path2,counter,n_xslots,oxygen_concentration);

            total_number_of_cells=0.0; 
            for(x_slot=0;x_slot<n_xslots;x_slot++){
                total_number_of_cells=total_number_of_cells+ncells_per_slot[x_slot];
            };
            fprintf(POPULATION,"%lf %ld \n",t+tau,total_number_of_cells);

            fprintf(stdout,"Time mark t=%lf\n",counter*sampling_time_window);
            
            counter++; //Not accurate for small simulations

        };//end if() for printing stuff,
        
        
        GlobalPDE_Handler(&parameters,
                        tau,
                        oxygen_concentration,
                        ncells_per_slot);
    
        Determine_event(&parameters,
        				&firing_slot,
                        &firing_reaction,
                        &firing_cell,
                        StochasticAges,
                        total_propensity*r2
                        );
    
        Advance_ages(StochasticAges,
                    n_xslots,
                    tau); 
   
        //////////////////////////////
        r3=drand48();
        Handle_event(&parameters,
        			StochasticAges,
                    firing_slot,    
                    firing_reaction,
                    firing_cell,
                    ncells_per_slot,
                    r3); 
        //Birth propensities are not updated here
    
        Update_birth_issues(&parameters,
        					StochasticAges,
                            oxygen_concentration
                            );                                    

        //safety checks and memory extensions
        ////Check_for_extra_memory(StochasticAges,n_xslots);
        //I do not see the need for anything else

        t=t+tau;

    }; //end main while

    fprintf(MESSAGES,"Tiempo final=%lf \n",t);

    
    /***************************************/
    

     
    /***************************************/
    
   
    //CLOSING FILES AND RETRIEVING MEMORY
    fclose(MESSAGES);
    fclose(POPULATION);

    for(i=0;i<n_xslots;i++){
        free(StochasticAges[i]->age_distribution);
        free(StochasticAges[i]);
    };
    free(StochasticAges);
    free(ncells_per_slot);
    free(oxygen_concentration);

    
	return(0);  
	
}      //end main




