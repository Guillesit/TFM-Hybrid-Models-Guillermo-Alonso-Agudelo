#include "Header.h"
#include <iostream>
#include <cmath> // For std::round

int main(int argc, char **argv)
{   
    // Get the starting point for the timer
    auto start = std::chrono::high_resolution_clock::now();
    
    FILE *INITIAL_POPULATION;
    FILE *INITIAL_OXYGEN;
    FILE *SIM_PARAMETERS;
    FILE *POPULATION;
    FILE *EQ_PARAMETERS;
    FILE *OUTPUT_DATA;
    FILE *MESSAGES;

    long n_files=N_OUTPUTS;
    //long n_xslots, i = 0;  
    double t=0.0,  tau=0.0; //Internal units
    double tstop=0.0; //SI units
    double sampling_time_window = 0.0;
    long counter=1; 
    double c=0.0;
    double delta_t; 
    double delta_x; 

    long i, x_slot, n_xslots=0;  //counters and tracers
    
    //long n_files;
    //double sampling_time_window = 0.0;
    //long counter=1; 
    



    //SI units: meters and seconds

    
    struct age_structure **StochasticAges; //SI units at input, then converted to internal units
    
    //depending on time and space
    //double *oxygen_level; //SI units, just for input-output
    //double *oxygen_concentration; //Internal units
    double *ncells_per_slot; //number of cells at each spatial slot, regardless of their age
    

    //At the only place where ew need densities of cells
    //We will rescale to get an effective rate constant
    double total_number_of_cells; //overall number of cells in the system

    //struct sim_parameters parameters; //SI units at input, then converted to internal units
    //double Delta_x; //grid size (SI units)
    //double Delta_t; //Suggested time increment (SI units)
    double n_eq=-1.0;
    double c_eq=-1.0;
    int hour; //to seed the random number generator
    double r1, r2, r3, r4; //random numbers

    double total_propensity; //Internal units
    long firing_slot=-1;
    long firing_cell=-1;
    int firing_reaction=0;//0=birth, 1=death, 2= left_diff, 3=right_diff




    double *oxygen_concentration; //concentration of oxygen at each spatial slot
    //double *ncells_per_slot; //concentration of cells at each spatial slot
 
    struct all_parameters parameters; //SI units at input, then converted to internal units
     
    int output_iter = 0; 

    char output_path1[200]="OutputValuesPopulation";
    char output_path2[200]="OutputValuesOxygen";
    char usertag[100]="";
    char label1[200]="OutputValuesPopulation";
    char label2[200]="OutputValuesOxygen";
    char label3[200]="messages";
    char label4[200]="total_population_vs_time";

    char create_folder[300]="mkdir ";

    if(argc!=6){
		printf("Error, incorrect argument number\n");
		printf("Use as: \n exec_name  \n initial_population_file \n initial_oxygen_file \n sim_parameters_file \n eq_parameters \n output_files_tag \n");
		return(1);
	};
    
    printf("Reading initial data\n"); 
    if((INITIAL_POPULATION=fopen(argv[1],"rt"))==NULL){
		printf("Error: initial_population_file could not be opened \n");
		return(1) ;
	};
    
    if((INITIAL_OXYGEN=fopen(argv[2],"rt"))==NULL){
		printf("Error: initial_oxygen_file could not be opened \n");
		return(1) ;
	};
    if((SIM_PARAMETERS=fopen(argv[3],"rt"))==NULL){
		printf("Error: equation_parameters_file could not be opened \n");
		return(1) ;
	};
    
    if((EQ_PARAMETERS=fopen(argv[4],"rt"))==NULL){
		printf("Error: simulation_parameters_file could not be opened \n");
		return(1) ;
	};
    
    strcpy(usertag,argv[5]);
    strcat(label4,usertag);

    if((POPULATION=fopen(label4,"w"))==NULL){
		printf("Error: output_population_file could not be opened \n");
		return(1) ;
	};
  
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
    delta_x=parameters.delta_x;
    printf("ESTE ES EL DELTA DE X %lf\n",delta_x);
    fprintf(MESSAGES, "#Suggested time increment (aka delta_t): %lf T \n",parameters.delta_t);
    fprintf(MESSAGES, "#Final time of the simulation: %lf T \n",parameters.tstop);
    fprintf(MESSAGES, "#Growth rate (r): %lf T \n",parameters.r);
    fprintf(MESSAGES, "#Carrying capacity (K): %lf T \n",parameters.K);
    fprintf(MESSAGES, "#Diffusion coefficient: %lf T \n",parameters.D);

    // Rescale constants 
    //Read oxygen

    n_xslots=Read_Init_Space_Oxygen(INITIAL_OXYGEN,&oxygen_concentration);
    printf("slots de oxygen %ld\n ", n_xslots); 

    //for (i = 0; i < n_xslots;i++){oxygen_concentration[i]=oxygen_concentration[i]/c_eq;};
  
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

    /*
    //Ordering from bigger to smaller ages within each spatial cell
    for (i = 0; i < n_xslots; i++){
        //printf("%ld\n",ncells_per_slot[i]); //DEBUG BUSSINESS
        qsort((StochasticAges[i]->age_distribution), ncells_per_slot[i], sizeof(double), cmpfunc);
    };
    */
    /*
    // Rescale oxygen concentration 
    //for (i = 0; i < n_xslots;i++){oxygen_concentration[i]=oxygen_concentration[i]/c_eq;};

    if (n_xslots!=Read_Init_Population(INITIAL_POPULATION,
                                &ncells_per_slot)){
        fprintf(stderr,"Mismatched number of spatial cells, aborting execution\n");
        exit(1);
    };

    parameters.n_xslots=n_xslots;
    */
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
    //printf(" este %ld %f\n",n_xslots,ncells_per_slot[0]);

    Print_Vector_Double(OUTPUT_DATA,output_path1,output_iter,n_xslots,ncells_per_slot);//error aqui(print cada uno)


    //Printing oxygen spatial distribution
    strcat(label2,"/Out"); //Final label for oxygen folder
    strcpy(output_path2,label2);
    Print_Vector_Double(OUTPUT_DATA,output_path2,output_iter,n_xslots,oxygen_concentration);
  

    printf("Time loop starting \n"); 
    printf("t is %lf \n", t); 
    printf("delta_t is %lf \n", delta_t); 
    printf("tstop is %lf \n", tstop); 
    double CFL=delta_x*delta_x/(2*parameters.D);
    printf("CFL condition is %lf \n", CFL); 

    if (delta_t>CFL){
      printf("\nWARNING: CFL condition not satisfied, PDE solver may be unstable!!\n\n");
    }


    //computing and printing total number of cells 
    
    total_number_of_cells=0.0;
    for(x_slot=0;x_slot<n_xslots;x_slot++){
        total_number_of_cells=total_number_of_cells+ncells_per_slot[x_slot];
    };

    fprintf(POPULATION,"%lf %lf \n",t,total_number_of_cells);
    

    //master handler   
    Initialize_propensities(&parameters,
                            StochasticAges,
                            oxygen_concentration
                            );


    // Create a random number generator
    std::mt19937 rng;

    // Seed the generator
    rng.seed(static_cast<unsigned long>(std::time(nullptr)));

    // Create a distribution from 0.0 to 1.0
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    //double *ncells_stoch;
    //double *ncells_deter;
    //double density[n_xslots]={};
    long interface=0;
    double resto;
    delta_t=parameters.delta_t;
    double time_advance;
    double prepare_next_reaction=1;
    c=0;
    Update_interface(&parameters,StochasticAges,ncells_per_slot,interface);

  // MAIN LOOP 
  while (t<tstop)
  {
      
    //Check if the interface is in the correct spot
    while((ncells_per_slot[interface]>parameters.K*0.5)&&(interface<=n_xslots-1)&&(prepare_next_reaction)){
      interface++;
      printf("Nueva frontera en %ld\n",interface);
    }

    if ((interface<=n_xslots-1)&&(ncells_per_slot[interface]>=1)){ //Stochastic part(we must make sure that there is a stoch part with any cell inside)
      
      if (prepare_next_reaction){//If no reaction is waiting its turn, we set the tau for the next reaction

        // Convert to integer with higher prob of floor or ceil depending on if the closest integer number is up or down

        r4=dist(rng);
        resto=std::modf(ncells_per_slot[interface],&ncells_per_slot[interface]);
      
        if (r4<resto){ 
          ncells_per_slot[interface]++;
          for (int i = 0; i < interface; ++i) {
            ncells_per_slot[i]-=(1-resto)/interface;
          }
        } else {
          for (int i = 0; i < interface; ++i) {
            ncells_per_slot[i]+=(resto)/interface;
          }
        }
    
        //if (ncells_per_slot[interface]!=std::floor(ncells_per_slot[interface])){
          //printf("Non integer number of cells in the interface\n");
        //}
        
        //printf("N of cells in the interface %lf\n",ncells_per_slot[interface]);
        prepare_next_reaction=0;
        // Update propensities in the interface to make up for the changes due to the deterministic flow inbetween reactions
        Update_interface(&parameters,StochasticAges,ncells_per_slot,interface);
        c=0;
        r1=dist(rng); 
        total_propensity=Compute_total_propensity(StochasticAges,n_xslots,interface); 
        //printf("total prop: %lf\n",total_propensity);
        tau=log(1/r1)/total_propensity;
        //printf("tau: %lf\n",tau);

      }

      if (c+delta_t<tau){//The time for the reaction hasnt come yet
        c+=delta_t;
        time_advance=delta_t;
      
      } else {//Do the reaction


        r2=dist(rng);
        //printf("r2: %lf, frac de prop: %lf\n",r2,r2*total_propensity);
        

        Determine_event(&parameters,
        				&firing_slot,
                        &firing_reaction,
                        &firing_cell,
                        StochasticAges,
                        total_propensity*r2,
                        interface
                        );

        r3=dist(rng);
        Handle_event(&parameters,
        			StochasticAges,
                    firing_slot,    
                    firing_reaction,
                    firing_cell,
                    ncells_per_slot,
                    r3); 


        time_advance=tau-c; //The remaining time to reach tau even if there have been deterministic advances in the meantime
        prepare_next_reaction=1;


      }
      
    } else {

      time_advance=delta_t;
    }
    
    if (interface>0){ //Deterministic part

      PDE_Handler(&parameters, ncells_per_slot,interface,time_advance);


    }

    t=t+time_advance;
    
          if(t>=counter*sampling_time_window){ //control in dimensional time (seconds)
    
            //Printing spatial distribution of cells 
            strcpy(output_path1,label1);
            Print_Vector_Double(OUTPUT_DATA,output_path1,counter,n_xslots,ncells_per_slot);
            
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




  

  };
    
  if (delta_t>CFL){
      printf("\nWARNING: CFL condition not satisfied, PDE solver may be unstable!!\n\n");
  } 
  free(ncells_per_slot); 
  free(oxygen_concentration); 
  //free(concentration); 

  // Get the ending point
  auto end = std::chrono::high_resolution_clock::now();
  // Calculate the duration
  std::chrono::duration<double> elapsed = end - start;
  // Print the elapsed time
  std::cout << "Elapsed time: " << elapsed.count() << " seconds\n" << std::endl;
} 