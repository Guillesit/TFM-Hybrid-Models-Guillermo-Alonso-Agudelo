

//One spatial dimension integrating through age structure

#include "Header.h"




int main(int argc, char **argv)
{
    FILE *INITIAL_DATUM;
    FILE *INITIAL_OXYGEN;
    FILE *SIM_DATA;
    FILE *POPULATION;
    FILE *PARAMETERS;
    FILE *OUTPUT_DATA;
    FILE *MESSAGES;
    
    
    long i, x_slot, n_xslots=0;  //counters and tracers

    long n_files=N_OUTPUTS;
    double sampling_time_window = 0.0;
    long counter=1; 
    
    double *density_of_individuals; //Internal units
    double *oxygen_concentration; //Internal units
    double *number_of_individuals_SI; //SI units
    double *oxygen_concentration_SI; //SI units
    struct sim_parameters parameters; //SI units at input, then converted to internal units
  
    double t=0.0; //SI units
    double tstop=0.0; //SI units
    int output_iter=0;

    double Delta_x, Delta_t; //Dimensional quantities (slot size, suggested time increment)
    double n_eq=-1.0;
    double c_eq=-1.0;
    
    double *division_threshold; //threshold parameter at birth rates
    double total_number_of_cells; //overall number of cells in the system
    
    /**/
    char output_path1[200]="OutputValuesPopulation";
    char output_path2[200]="OutputValuesOxygen";
    char usertag[100]="";
    char label1[200]="OutputValuesPopulation";
    char label2[200]="OutputValuesOxygen";
    char label3[200]="messages";
    char label4[200]="total_population_vs_time";
    char create_folder[300]="mkdir ";

   	
    /************************/
    
    /*****************/
    //OPENING FILES
    
    if(argc!=6){
		printf("Error, incorrect argument number\n");
		printf("Use as: \n exec_name  \n initial_population_file \n initial_oxygen_file \n population_parameters_file \n simulation_parameters_file \n output_files_tag \n");
		return(1);
	};
    
    if((INITIAL_DATUM=fopen(argv[1],"rt"))==NULL){
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
    //Read population
    n_xslots=Read_Init_Space(INITIAL_DATUM,&density_of_individuals);
    
    //Read oxygen
    if (n_xslots!=Read_Init_Space_Oxygen(INITIAL_OXYGEN,&oxygen_concentration)){
        fprintf(stderr,"Mismatched number of spatial cells, aborting execution\n");
        exit(1);
    };


    //define auxiliar vectors
    if((number_of_individuals_SI=(double *) malloc(n_xslots*sizeof(double)))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    }; 
    if((oxygen_concentration_SI=(double *) malloc(n_xslots*sizeof(double)))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    }; 
    
    //Read population parameters 
    Read_params_population(&parameters,
                            PARAMETERS
                            );
    
    //Read parameters (simulation infos)
    Read_params_sim(&parameters,
                    SIM_DATA,
                    &tstop,  
                    &Delta_x,
                    &Delta_t,
                    &n_files
                    );
    parameters.n_xslots=n_xslots;

    sampling_time_window=tstop/n_files;

    
    
    /*******************/
    //PRINTING SIMULATION INFOS

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
    fprintf(MESSAGES,"#Parameters file (simulation):%s\n \n",argv[4]);
    
    fprintf(MESSAGES,"#Default unit system as given by the S.I. (seconds, meters, moles),\n#let X for space, T for time\n#and M for quantity of substance\n");
    fprintf(MESSAGES,"#Grid step (aka delta_x): %lf X \n",Delta_x);
    fprintf(MESSAGES,"#Suggested time increment (aka delta_t): %lf T \n",Delta_t);
    fprintf(MESSAGES,"#Sampling time window: %lf T \n",sampling_time_window);
    fprintf(MESSAGES,"#User-provided final time: %lf T \n",tstop);
    fprintf(MESSAGES,"#Number of output files: %ld \n \n",n_files);
    
    fprintf(MESSAGES,"#survival rate: %lf \n",(parameters.survival_rate));
    fprintf(MESSAGES,"#Death rate: %.10lf T^-1 \n",parameters.death_rate_hat);
    fprintf(MESSAGES,"#tau_p: %.10lf T\n",parameters.tau_p);
    fprintf(MESSAGES,"#aplus: %.10lf T\n",parameters.aplus);
    fprintf(MESSAGES,"#p6/p3: %.10lf \n",parameters.p6_over_p3);
    fprintf(MESSAGES,"#Diffusion coefficients (pop, oxy): %.10lf %lf X^2/T\n",parameters.diff_coef_pop,parameters.diff_coef_oxygen);
    fprintf(MESSAGES,"#k_decay: %.10lf T^-1\n",parameters.k_decay_hat);
    fprintf(MESSAGES,"#k_consumption: %.10lf X T^-1 M^-1\n",parameters.k_consumption_hat);
    fprintf(MESSAGES,"#k_source: %.10lf M X^-1 T^-1\n \n",parameters.source_oxygen_hat);

    Compute_equilibria(&parameters,&n_eq,&c_eq); //En el sistema de unidades primitivo

   	fprintf(MESSAGES,"#neq=%lf M X^-1, ceq=%.10lf M X^-1\n",n_eq,c_eq);
    //fprintf(MESSAGES,"Note: spatial outputs are dimensionless");

    //ADIMENSIONALIZE POPULATION & OXYGEN
    //Esto de aqui hasta el final de la simulacion, incluso para las salidas

    ////////////////
    //Rescale constants to the internal unit system (ADIMENSIONAL FROM NOW ON)
    parameters.death_rate_hat=(parameters.death_rate_hat)*(parameters.tau_p);
    parameters.aminus_hat=aminus/(parameters.tau_p);
    parameters.diff_coefs_ratio=(parameters.diff_coef_oxygen)/(parameters.diff_coef_pop);
    parameters.k_decay_hat=(parameters.k_decay_hat)*(parameters.tau_p);
    parameters.k_consumption_hat=(parameters.k_consumption_hat)*(parameters.tau_p)*n_eq; 
    parameters.source_oxygen_hat=(parameters.source_oxygen_hat)*(parameters.tau_p)/c_eq;
    parameters.typical_lengthscale=sqrt((parameters.diff_coef_pop)*(parameters.tau_p));

    parameters.critical_oxy_hat=(typical_oxy*severe_hypoxia)*(c_cr(parameters.p6_over_p3)/c_cr(1))/c_eq;

    parameters.Delta_x_hat=Delta_x/(parameters.typical_lengthscale);
    parameters.Delta_t_hat=Delta_t/(parameters.tau_p);

    //Effective diffusion coefficients to take into account second order finite differences:
    parameters.diff_coef_eff=1.0/((parameters.Delta_x_hat)*(parameters.Delta_x_hat));
    parameters.diff_coefs_ratio=parameters.diff_coefs_ratio/((parameters.Delta_x_hat)*(parameters.Delta_x_hat));
  
    //////Easier to determine CFL condition here (at least the part that does not depend on the iteration)
    parameters.CFL_number=min(1.0/(6*(parameters.diff_coefs_ratio)),
                    			min(1.0/(6*(parameters.diff_coef_eff)),
                            		min(1.0/(1.0+(parameters.k_decay_hat)),1.0/(1.0+(parameters.death_rate_hat)))    
                            	) 
                        	);
                  
               
     
    /*****************/
    //ADIMENSIONALIZING AND PRINTING INITIAL CONDITIONS
  
     //create "OutputValuesPopulation" and "OutputValuesOxygen" (plus user tag)
     strcat(label1,usertag);
     strcat(label2,usertag);
     
     strcat(create_folder,label1);
     system(create_folder);
     strcpy(create_folder,"mkdir ");
     strcat(create_folder,label2);
     system(create_folder);
     
     //Printing cell concentration regardless of their age
     strcat(label1,"/Out"); //Final label for population folder
     strcpy(output_path1,label1);
     //Adimensionalizing the population:
     // Save the number of individuals
     //for (i = 0; i < n_xslots;i++){
     //           number_of_individuals_SI[i]=density_of_individuals[i];
     //   };
     //Print_Vector(OUTPUT_DATA,output_path1,output_iter,n_xslots,number_of_individuals_SI);

     Print_Vector(OUTPUT_DATA,output_path1,output_iter,n_xslots,density_of_individuals);
     for (i = 0; i < n_xslots;i++){density_of_individuals[i]=density_of_individuals[i]/n_eq;};
    
     //Printing oxygen spatial distribution
     strcat(label2,"/Out"); //Final label for oxygen folder
     strcpy(output_path2,label2);
     //Adimensionalizing the oxygen:
     Print_Vector(OUTPUT_DATA,output_path2,output_iter,n_xslots,oxygen_concentration);
     for (i = 0; i < n_xslots;i++){oxygen_concentration[i]=oxygen_concentration[i]/c_eq;};
    
    
    //computing and printing total number of cells 
    //DOUBLE CHECK THIS
    total_number_of_cells=0.0;
    for(x_slot=0;x_slot<n_xslots;x_slot++){
        total_number_of_cells=total_number_of_cells+density_of_individuals[x_slot]*n_eq*Delta_x;
    };
    fprintf(POPULATION,"%lf %lf \n",t,total_number_of_cells);
   
    /*******************************/
    //INITIALIZING DIVISION RATES

    //book division_threshold 
    if((division_threshold= (double *) malloc(sizeof(double)*n_xslots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    Compute_division_threshold(&parameters,
                                division_threshold,
                                oxygen_concentration,
                                n_xslots-1);
    //Already adimensional times
    
    
    /******************************/
    
    
    /******************************/
    //MAIN LOOP
    
    while(t<tstop){   
        //Control time is dimensional, simulation time is not 
        
        //Explicit Euler time marching.
        
        //Main submodule that:
        //ADVANCEs THE POPULATION
        //ADVANCEs THE OXYGEN
        //RECOMPUTEs SETTINGS FOR NEXT ITERATION


        GlobalPDE_Handler(&parameters,
                        density_of_individuals,
                        division_threshold,
                        oxygen_concentration
                        );
    
        //recompute division_threshold
        Compute_division_threshold(&parameters,
                                division_threshold,
                                oxygen_concentration,
                                n_xslots-1);
        
        
       //PRINTING DATA ON SELECTED ITERATIONS
        if(t+parameters.Delta_t_hat>=counter*sampling_time_window){ //control in dimensional time (seconds)

            //Printing spatial distribution of cells 
            strcpy(output_path1,label1);
            for (i = 0; i < n_xslots;i++){
                number_of_individuals_SI[i]=density_of_individuals[i]*n_eq*Delta_x;
                //fprintf(stdout, "%lf %lf\n", number_of_individuals_SI[i],density_of_individuals[i] );
            };
            Print_Vector(OUTPUT_DATA,output_path1,counter,n_xslots,number_of_individuals_SI);

            //Printing oxygen spatial distribution
            strcpy(output_path2,label2);
            for (i = 0; i < n_xslots;i++){oxygen_concentration_SI[i]=oxygen_concentration[i]*c_eq;};
            Print_Vector(OUTPUT_DATA,output_path2,counter,n_xslots,oxygen_concentration_SI);


            //printing total number of cells 
            total_number_of_cells=0.0;
            for(x_slot=0;x_slot<n_xslots;x_slot++){
                total_number_of_cells=total_number_of_cells+density_of_individuals[x_slot]*n_eq*Delta_x;
            };
            fprintf(POPULATION,"%lf %lf \n",t+parameters.Delta_t_hat,total_number_of_cells);

            fprintf(stdout,"Time mark t=%lf\n",counter*sampling_time_window);
            //This output time is dimensional (measured in seconds, right?)
        
            counter++; //Not accurate for small simulations???

        }; //end if printing outputs
        
        t=t+parameters.Delta_t_hat;
       

    }; //end main loop

    
    /***************************************/
    

     
    /***************************************/
    
    
    //CLOSING FILES AND RETRIEVING MEMORY
    
    free(density_of_individuals);
    free(oxygen_concentration);
    free(number_of_individuals_SI);
    free(oxygen_concentration_SI);
   	free(division_threshold);
    fclose(POPULATION);
    
	return(0);  
	
}      //end main




