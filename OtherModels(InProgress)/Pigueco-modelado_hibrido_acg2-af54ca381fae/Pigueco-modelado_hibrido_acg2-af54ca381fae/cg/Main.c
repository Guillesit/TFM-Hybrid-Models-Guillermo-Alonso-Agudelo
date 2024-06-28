//We don't include here a damping term at oxigen's equation so far


//One spatial dimension integrating through age structure
//Lots of generalizations still to come...
//Part of the process could be paralelized if we knew how


#include "Header.h"




int main(int argc, char **argv)
{
    FILE *INITIAL_DATUM, *POPULATION, *PARAMETERS, *OUTPUT_DATA, *MESSAGES;
    
    long age_slots=2; //Number of age slots (we will handle this with malloc-realloc)
    long i_iteration, i_slot,i,x_slot, n_xslots=0;  //counters and tracers
    
    double t=0.0;
    double delta_t=1.0;
    double delta_x=1.0; //user's value for this quantity, in mm.
    

    //QUANTITIES TO MONITOR AND SYSTEM PARAMETERS
    
    double death_rate_inv; //death rates^-1
    double tau_p; //parameter related to stepwise birth rates
    double *division_threshold; //threshold parameter at birth rates
    double source_oxygen, k_oxygen; //parameters in the equation for the oxigen
    
    
    double *oxygen_level; //depending on time and space
    double p6_over_p3; //ratio of p's
    double aplus; //switch at zero oxigen, if applicable
    double *number_of_cells,*backup_number_of_cells; //number of cells at each spatial slot, regardless of their age
    double total_number_of_cells; //overall number of cells in the system

    char case_flag; //true if we need to compute aplus, false otherwise
    double peak_number_of_cells;
    
    
	//Stuff to handle output files
    long sampling_period=1; //sampling period
    long number_of_iterations; //100*sampling_period as designed
    char output_path[70]="OutputValuesPopulation/Out";
    int output_iter=0;

   
   
    /************************/
    
    /*****************/
    //OPENING FILES
    
    if(argc!=3){
		printf("Error, incorrect argument number\n");
		printf("Use as: \n exec_name  \n initial_data_file \n parameters_file \n");
		return(1);
	};
    
    if((INITIAL_DATUM=fopen(argv[1],"rt"))==NULL){
		printf("Error: initial_data_file could not be opened \n");
		return(1) ;
	};
    
    if((PARAMETERS=fopen(argv[2],"rt"))==NULL){
		printf("Error: parameters_file could not be opened \n");
		return(1) ;
	};
    
    if((POPULATION=fopen("total_population_vs_time","w"))==NULL){
		printf("Error: output_population_file could not be opened \n");
		return(1) ;
	};

   
    
    /*********************/
    //GATHERING ADDITIONAL INITIAL INFORMATION AND INITIALIZING
    
    n_xslots=Read_Init_Space(INITIAL_DATUM,PARAMETERS,&death_rate_inv,&tau_p,&source_oxygen,&k_oxygen,&oxygen_level,&p6_over_p3,&number_of_cells,&backup_number_of_cells);
    
    //ask for sampling period

    printf("Enter the sampling period\n");
    scanf("%ld",&sampling_period);
    //sampling_period=1;// 30;//10000;//5000;
    //number_of_iterations=100*sampling_period;
    number_of_iterations=((int)N_OUTPUTS)*sampling_period; //N_OUTPUTS is set in the header file
    
    delta_x=1; //Future work: ask user or get this from a data file
    
    
    /*******************/
    //PRINTING INITIAL CONDITION
    
    strcpy(output_path,"OutputValuesPopulation/Out");
    Print_Vector(OUTPUT_DATA,output_path,output_iter,n_xslots,number_of_cells);
   
    //computing and printing total number of cells
    total_number_of_cells=0.0;
    for(x_slot=0;x_slot<n_xslots;x_slot++){
        total_number_of_cells=total_number_of_cells+number_of_cells[x_slot];
    };
    fprintf(POPULATION,"%lf %lf \n",t,total_number_of_cells);

    //Printing oxygen spatial distribution
    strcpy(output_path,"OutputValuesOxygen/Out");
    Print_Vector(OUTPUT_DATA,output_path,output_iter,n_xslots,oxygen_level);
 
    
    /*******************************/
    //PRINTING SIMULATION INFOS
    
    if((MESSAGES=fopen("messages","w"))==NULL){ //Master output file
		printf("Error: output messages file could not be created \n");
		return(1) ;
	};
    
    fprintf(MESSAGES,"#Debug info: \n \n");
    fprintf(MESSAGES,"#Data generated with the following .exe file:%s\n",argv[0]);
    fprintf(MESSAGES,"#Initial condition file:%s\n",argv[1]);
    fprintf(MESSAGES,"#Parameters file:%s\n",argv[2]);
    fprintf(MESSAGES,"#Event number: %ld \n \n",number_of_iterations);
    fprintf(MESSAGES,"#Specific infos: oxygen damping, eigenvalue shortcut, a single cellular line\n");
    
    fprintf(MESSAGES,"#Death rate: %lf \n",1.0/death_rate_inv);
    fprintf(MESSAGES,"#tau_p: %lf\n",tau_p);
    fprintf(MESSAGES,"#source_oxygen: %lf\n",source_oxygen);
    fprintf(MESSAGES,"#k_oxygen: %lf\n",k_oxygen);
    fprintf(MESSAGES,"#p6/p3: %lf\n",p6_over_p3);
    fprintf(MESSAGES,"#slot size: %lf mm\n \n \n",delta_x);
     
    fprintf(MESSAGES,"#List of issues:\n");
    //So far "issues" are printed in stderr, to avoid passing *MESSAGES to subroutines
    //Think about this...
    
    
    /*******************************/
    //INITIALIZING DIVISION RATES
   
    //We need to compute the division at zero oxigen for the sake of formulas
    //Just call the ODE solver once
    
    //book division_threshold 
    if((division_threshold= (double *) malloc(sizeof(double)*n_xslots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    
     //using interpolation formula by de la Cruz et al.
    if(p6_over_p3>r_cr){
        case_flag=1;
        aplus=Compute_aplus(p6_over_p3);
        for(x_slot=0;x_slot<n_xslots;x_slot++){
            division_threshold[x_slot]=First_division_threshold(aplus,oxygen_level[x_slot]);
        };
    }//end if p6/p3 above threshold
    else{
        case_flag=0;
        for(x_slot=0;x_slot<n_xslots;x_slot++){
            if(oxygen_level[x_slot]<=c_cr(p6_over_p3)){division_threshold[x_slot]=INFINITE_TIME;}
            //note that this WILL NOT freeze cells when they enter quiescence
            else{
                division_threshold[x_slot]=Second_division_threshold(p6_over_p3,oxygen_level[x_slot]);
            };
        }; //end spatial for
    }; //end case p6/p3 below threshold
     
   
     
    
    /******************************/
    
    
    /******************************/
    //MAIN LOOP
    
    
    for(i_iteration=1;i_iteration<=number_of_iterations;i_iteration++){
        
        
        //Explicit Euler time marching.
        //First we advance the population, then the oxygen
        //We DO NOT use the updated population to advance the oxygen

        delta_t=1;
        peak_number_of_cells=GetMax(number_of_cells,n_xslots); //Useful at some point to assess stability


        
        //ADVANCE THE POPULATION
        
        //boundary conditions are determined in this module
        Finite_Difference_Solver(delta_t,number_of_cells,death_rate_inv,division_threshold,tau_p,n_xslots,delta_x);
    
        /************/
        
        
        //ADVANCE THE OXIGEN
    
        OxygenSolverDamp(oxygen_level,n_xslots,backup_number_of_cells,delta_t,k_oxygen,source_oxygen,peak_number_of_cells,delta_x);
        
        
        /********************/
        //RECOMPUTE SETTINGS FOR NEXT ITERATION

        for(x_slot=0;x_slot<n_xslots;x_slot++){
            backup_number_of_cells[x_slot]=number_of_cells[x_slot];
        };
    
        t=t+delta_t;
        

        //recompute division_threshold
        if(case_flag){
            for(x_slot=0;x_slot<n_xslots;x_slot++){
                division_threshold[x_slot]=First_division_threshold(aplus,oxygen_level[x_slot]);
            };
        }
        else{
            for(x_slot=0;x_slot<n_xslots;x_slot++){
                if(oxygen_level[x_slot]<=c_cr(p6_over_p3)){
                    division_threshold[x_slot]=INFINITE_TIME;
                }
                else{
                    division_threshold[x_slot]
                    =Second_division_threshold(p6_over_p3,oxygen_level[x_slot]);
                };
            };// end for
        };
       
        
        
        /******************/
       //PRINTING DATA ON SELECTED ITERATIONS
        
        if((i_iteration%sampling_period)==0){
            
            output_iter=i_iteration/sampling_period;
            
            //printing total number of cells
            total_number_of_cells=0.0;
            for(x_slot=0;x_slot<n_xslots;x_slot++){
                total_number_of_cells=total_number_of_cells+number_of_cells[x_slot];
            };
            fprintf(POPULATION,"%lf %lf \n",t,total_number_of_cells);
            
            //Printing spatial distribution of cells
            strcpy(output_path,"OutputValuesPopulation/Out");
            Print_Vector(OUTPUT_DATA,output_path,output_iter,n_xslots,number_of_cells);
            
            //Printing oxygen spatial distribution
            strcpy(output_path,"OutputValuesOxygen/Out");
            Print_Vector(OUTPUT_DATA,output_path,output_iter,n_xslots,oxygen_level);
            
        
        }; //end if printing outputs
        
    }; //end main loop

    
    /***************************************/
    

     
    /***************************************/
    
    
    //CLOSING FILES AND RETRIEVING MEMORY
    
    free(number_of_cells);
    free(backup_number_of_cells);
    free(oxygen_level);
    free(division_threshold);
    
    fclose(POPULATION);
    fclose(MESSAGES);
    
	return(0);  
	
}      //end main




