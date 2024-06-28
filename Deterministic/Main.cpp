#include "Header.h"

int main(int argc, char **argv)
{
    // Get the starting point for the timer
    auto start = std::chrono::high_resolution_clock::now();

    FILE *INITIAL_POPULATION;
    FILE *INITIAL_OXYGEN;
    FILE *SIM_PARAMETERS;
    //FILE *POPULATION;
    FILE *EQ_PARAMETERS;
    FILE *OUTPUT_DATA;
    FILE *MESSAGES;

    long n_files=N_OUTPUTS;
    long n_xslots, i = 0;  
    double t=0.0,  tau=0.0; //Internal units
    double tstop=0.0; //SI units
    double sampling_time_window = 0.0;
    long counter=1; 
    double delta_t; 
    double delta_x; 

    double *oxygen_concentration; //concentration of oxygen at each spatial slot
    double *cells_population; //concentration of cells at each spatial slot
 
    struct all_parameters parameters; //SI units at input, then converted to internal units
     
    int output_iter = 0; 

    char output_path1[200]="OutputValuesPopulation";
    char output_path2[200]="OutputValuesOxygen";
    char usertag[100]="";
    char label1[200]="OutputValuesPopulation";
    char label2[200]="OutputValuesOxygen";
    char label3[200]="messages";

    char create_folder[300]="mkdir ";

    if(argc!=6){
		printf("Error, incorrect argument number\n");
		printf("Use as: \n exec_name  \n initial_population_file \n initial_oxygen_file \n sim_parameters_file \n eq_parameters \n output_files_tag \n");
		return(1);
	};
    
    printf("Reading initial data"); 
    if((INITIAL_POPULATION=fopen(argv[1],"rt"))==NULL){
		printf("Error: initial_population_file could not be opened \n");
		return(1) ;
	};
    
    if((INITIAL_OXYGEN=fopen(argv[2],"rt"))==NULL){
		printf("Error: initial_oxygen_file could not be opened \n");
		return(1) ;
	};
    printf("Hola"); 
    if((SIM_PARAMETERS=fopen(argv[3],"rt"))==NULL){
		printf("Error: equation_parameters_file could not be opened \n");
		return(1) ;
	};
    
    if((EQ_PARAMETERS=fopen(argv[4],"rt"))==NULL){
		printf("Error: simulation_parameters_file could not be opened \n");
		return(1) ;
	};
    
    strcpy(usertag,argv[5]);

    // READING INITIAL DATA 

    Read_sim_params(&parameters, SIM_PARAMETERS);
    Read_eq_params(&parameters, EQ_PARAMETERS);

    delta_t = parameters.delta_t; 
    tstop = parameters.tstop; 

    sampling_time_window=tstop/n_files;
    printf("Sampling time window is %lf\n", sampling_time_window); 

    strcat(label3,usertag);
    if((MESSAGES=fopen(label3,"w"))==NULL){ 
		printf("Error: output messages file could not be created \n");
		return(1) ;
	};

    fprintf(MESSAGES,"#Debug info: \n \n");
    fprintf(MESSAGES,"#Data generated with the following .exe file:%s\n",argv[0]);
    fprintf(MESSAGES,"#Initial condition file (population ages):%s\n",argv[1]);
    fprintf(MESSAGES,"#Initial condition file (oxygen):%s\n",argv[2]);
    fprintf(MESSAGES,"#Equation parameters file:%s\n",argv[3]);
    fprintf(MESSAGES,"#Simulation parameters file:%s\n",argv[4]);

    fprintf(MESSAGES, "#Grid step (aka delta_x): %lf X \n",parameters.delta_x);
    fprintf(MESSAGES, "#Suggested time increment (aka delta_t): %lf T \n",parameters.delta_t);
    fprintf(MESSAGES, "#Final time of the simulation: %lf T \n",parameters.tstop);
    fprintf(MESSAGES, "#Growth rate (r): %lf T \n",parameters.r);
    fprintf(MESSAGES, "#Carrying capacity (K): %lf T \n",parameters.K);
    fprintf(MESSAGES, "#Diffusion coefficient: %lf T \n",parameters.D);

    // Rescale constants 

    n_xslots=Read_Init_Space_Oxygen(INITIAL_OXYGEN,&oxygen_concentration);
    printf("slots de oxygen %ld\n ", n_xslots); 

    // Rescale oxygen concentration 
    //for (i = 0; i < n_xslots;i++){oxygen_concentration[i]=oxygen_concentration[i]/c_eq;};

    if (n_xslots!=Read_Init_Population(INITIAL_POPULATION,
                                &cells_population)){
        fprintf(stderr,"Mismatched number of spatial cells, aborting execution\n");
        exit(1);
    };

    parameters.n_xslots=n_xslots;

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
    Print_Vector_Double(OUTPUT_DATA,output_path1,output_iter,n_xslots,cells_population);

    //Printing oxygen spatial distribution
    strcat(label2,"/Out"); //Final label for oxygen folder
    strcpy(output_path2,label2);
    Print_Vector_Double(OUTPUT_DATA,output_path2,output_iter,n_xslots,oxygen_concentration);

    printf("Time loop starting \n"); 
    printf("t is %lf \n", t); 
    printf("delta_t is %lf \n", delta_t); 
    printf("tstop is %lf \n", tstop); 


    // MAIN LOOP 
    while (t<tstop)
    {

      PDE_Handler(&parameters, cells_population); 

      if(t+delta_t>=counter*sampling_time_window){ //control in dimensional time (seconds)
    
            //Printing spatial distribution of cells 
            strcpy(output_path1,label1);
            Print_Vector_Double(OUTPUT_DATA,output_path1,counter,n_xslots,cells_population);
            
            //Printing oxygen spatial distribution
            // strcpy(output_path2,label2);
            // Print_Vector(OUTPUT_DATA,output_path2,counter,n_xslots,oxygen_concentration);

            //printing total number of cells 
            // total_number_of_cells=0.0;
            // for(x_slot=0;x_slot<n_xslots;x_slot++){
            //     total_number_of_cells=total_number_of_cells+density_of_individuals[x_slot]*n_eq*Delta_x;
            // };
            // fprintf(POPULATION,"%lf %lf \n",t+Delta_t,total_number_of_cells);

            fprintf(stdout,"Time mark t=%lf\n",counter*sampling_time_window);
            //This output time is dimensional (measured in seconds, right?)
        
            counter++; //Not accurate for small simulations???

      }; //end if printing outputs
      
      t = t + delta_t; 
    }
    
    free(cells_population); 
    free(oxygen_concentration); 
    //free(concentration); 

    // Get the ending point
  auto end = std::chrono::high_resolution_clock::now();
  // Calculate the duration
  std::chrono::duration<double> elapsed = end - start;
  // Print the elapsed time
  std::cout << "Elapsed time: " << elapsed.count() << " seconds\n" << std::endl;

} 