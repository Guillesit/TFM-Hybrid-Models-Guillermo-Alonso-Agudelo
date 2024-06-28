

//One spatial dimension integrating through age structure

#include "Header.h"




int main(int argc, char **argv)
{
    FILE *INITIAL_DATUM_host,*INITIAL_DATUM_invader;
    FILE *INITIAL_OXYGEN,*SIM_DATA;
    FILE *POPULATION_host, *POPULATION_invader;
    FILE *PARAMETERS_host, *PARAMETERS_invader;
    FILE *OUTPUT_DATA, *MESSAGES;
    
    
    long i_iteration, i_slot,i,x_slot, n_xslots=0;  //counters and tracers
    
    double *number_of_cells_host; //number of cells at each spatial slot, regardless of their age
    double *number_of_cells_invader;
    double *total_number_of_cells;
    double *oxygen_level; //depending on time and space
    
    
    
    
    double t=0.0;
    
    double delta_t=1.0;
    double delta_x=1.0; //user's value for this quantity, in mm.
    

    //QUANTITIES TO MONITOR AND SYSTEM PARAMETERS
    
    double death_rate_inv_host; //death rates^-1
    double tau_p_host; //parameter related to stepwise birth rates
    double p6_over_p3_host; //ratio of p's
    double diff_coef_host;
    double *division_threshold_host; //threshold parameter at birth rates
    
    double death_rate_inv_invader; //death rates^-1
    double tau_p_invader; //parameter related to stepwise birth rates
    double p6_over_p3_invader; //ratio of p's
    double diff_coef_invader;
    double *division_threshold_invader; //threshold parameter at birth rates
    
    double total_number_of_cells_host, total_number_of_cells_invader; //overall number of cells in the system
    
    
    
    double source_oxygen, k_oxygen, diff_coef_oxygen; //parameters in the equation for the oxygen
    
    double aplus_host=-1, aplus_invader=-1; //switch at zero oxygen, if applicable
    
    //char case_flag_host, case_flag_invader; //true if we need to compute aplus, false otherwise
    double peak_number_of_cells;
    
    
    
	//Stuff to handle output files
    long sampling_period=1; //sampling period
    long number_of_iterations; //typically 100*sampling_period
    char output_path[80]="OutputValuesPopulationHost/Out";
    int output_iter=0;
    
    /**/
    char output_path1h[200]="OutputValuesPopulationHost";
    char output_path1i[200]="OutputValuesPopulationInvader";
    char output_path2[200]="OutputValuesOxygen";
    char usertag[100]="";
    char label1h[200]="OutputValuesPopulationHost";
    char label1i[200]="OutputValuesPopulationInvader";
    char label2[200]="OutputValuesOxygen";
    char label3[200]="messages";
    char label4h[200]="total_population_vs_time_host";
    char label4i[200]="total_population_vs_time_invader";
    char create_folder[300]="mkdir ";

   
   
    /************************/
    
    /*****************/
    //OPENING FILES
    
    if(argc!=8){
		printf("Error, incorrect argument number\n");
		printf("Use as: \n exec_name  \n initial_data_file_host \n initial_data_file_invader  \n initial_oxygen \n parameters_file_host \n parameters_file_invader \n simulation_info \n output_files_tag \n");
		return(1);
	};
    
    if((INITIAL_DATUM_host=fopen(argv[1],"rt"))==NULL){
		printf("Error: initial_data_file_host could not be opened \n");
		return(1) ;
	};
    
    if((INITIAL_DATUM_invader=fopen(argv[2],"rt"))==NULL){
		printf("Error: initial_data_file_invader could not be opened \n");
		return(1) ;
	};
    
    if((INITIAL_OXYGEN=fopen(argv[3],"rt"))==NULL){
		printf("Error: initial_oxygen could not be opened \n");
		return(1) ;
	};
    
    if((PARAMETERS_host=fopen(argv[4],"rt"))==NULL){
		printf("Error: parameters_file_host could not be opened \n");
		return(1) ;
	};
    
    if((PARAMETERS_invader=fopen(argv[5],"rt"))==NULL){
		printf("Error: parameters_file_invader could not be opened \n");
		return(1) ;
	};
    
    if((SIM_DATA=fopen(argv[6],"rt"))==NULL){
		printf("Error: simulation_info could not be opened \n");
		return(1) ;
	};
    
    
    strcpy(usertag,argv[7]);
    strcat(label4h,usertag);
    strcat(label4i,usertag);
    
    if((POPULATION_host=fopen(label4h,"w"))==NULL){
		printf("Error: output_population_file could not be opened \n");
		return(1) ;
	};
    
    if((POPULATION_invader=fopen(label4i,"w"))==NULL){
		printf("Error: output_population_file could not be opened \n");
		return(1) ;
	};

   
   
    /*********************/
    //GATHERING ADDITIONAL INITIAL INFORMATION AND INITIALIZING
    
    
    //Read first population
    n_xslots=Read_Init_Space(INITIAL_DATUM_host,&number_of_cells_host);
    
    
    //Read second population
    if (n_xslots!=Read_Init_Space(INITIAL_DATUM_invader,&number_of_cells_invader)
        ){
        fprintf(stderr,"Mismatched number of spatial cells, aborting execution\n");
        exit(1);
    };
    
    //Read oxygen
    if (n_xslots!=Read_Init_Space_Oxygen(INITIAL_OXYGEN,&oxygen_level)){
        fprintf(stderr,"Mismatched number of spatial cells, aborting execution\n");
        exit(1);
    };
    
    //Read parameters (first population)
    Read_params_population(PARAMETERS_host,&death_rate_inv_host,&tau_p_host,&p6_over_p3_host,&diff_coef_host);
    
    //Read parameters (second population)
    Read_params_population(PARAMETERS_invader,&death_rate_inv_invader,&tau_p_invader,&p6_over_p3_invader,&diff_coef_invader);
    
    
    //Read parameters (simulation infos)
    Read_params_sim(SIM_DATA,&source_oxygen,
                         &k_oxygen,
                         &diff_coef_oxygen,
                         &delta_t,
                         &delta_x,
                         &sampling_period
                    );
    
    //Book total number of cells
    if((total_number_of_cells= (double *) malloc(sizeof(double)*n_xslots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    
    
    //sampling_period=1;// 30;//10000;//5000;
    //number_of_iterations=100*sampling_period;
    number_of_iterations=((int)N_OUTPUTS)*sampling_period; //N_OUTPUTS is set in the header file
    
    
    
    /*******************/
    //PRINTING INITIAL CONDITION
    
    
     //create "OutputValuesPopulation" and "OutputValuesOxygen" (plus user tag)
     strcat(label1h,usertag);
     strcat(label1i,usertag);
     strcat(label2,usertag);
     
     strcat(create_folder,label1h);
     system(create_folder);
     strcpy(create_folder,"mkdir ");
     strcat(create_folder,label1i);
     system(create_folder);
     strcpy(create_folder,"mkdir ");
     strcat(create_folder,label2);
     system(create_folder);
     
     
     //Printing spatial DENSITY of cells regardless of their age
     strcat(label1h,"/Out"); //Final label for population folder
     strcpy(output_path1h,label1h);
     Print_Vector(OUTPUT_DATA,output_path1h,output_iter,n_xslots,number_of_cells_host);
     
     strcat(label1i,"/Out"); //Final label for population folder
     strcpy(output_path1i,label1i);
     Print_Vector(OUTPUT_DATA,output_path1i,output_iter,n_xslots,number_of_cells_invader);
     
     //Printing oxygen spatial distribution
     strcat(label2,"/Out"); //Final label for oxygen folder
     strcpy(output_path2,label2);
     Print_Vector(OUTPUT_DATA,output_path2,output_iter,n_xslots,oxygen_level);
    
    
    //computing and printing total number of cells (host)
    total_number_of_cells_host=0.0;
    for(x_slot=0;x_slot<n_xslots;x_slot++){
        total_number_of_cells_host=total_number_of_cells_host+number_of_cells_host[x_slot];
    };
    fprintf(POPULATION_host,"%lf %lf \n",t,total_number_of_cells_host);
    
    //computing and printing total number of cells (invader)
    total_number_of_cells_invader=0.0;
    for(x_slot=0;x_slot<n_xslots;x_slot++){
        total_number_of_cells_invader=total_number_of_cells_invader+number_of_cells_invader[x_slot];
    };
    fprintf(POPULATION_invader,"%lf %lf \n",t,total_number_of_cells_invader);
    

 
    
    /*******************************/
    //PRINTING SIMULATION INFOS
    
    //create "messages" (plus user tag)
    strcat(label3,usertag);
    if((MESSAGES=fopen(label3,"w"))==NULL){ //Master output file
		printf("Error: output messages file could not be created \n");
		return(1) ;
	};
    
    
    fprintf(MESSAGES,"#Debug info: \n \n");
    fprintf(MESSAGES,"#Data generated with the following .exe file:%s\n",argv[0]);
    fprintf(MESSAGES,"#Initial condition file (host):%s\n",argv[1]);
    fprintf(MESSAGES,"#Initial condition file (invader):%s\n",argv[2]);
    fprintf(MESSAGES,"#Initial condition file (oxygen):%s\n",argv[3]);
    fprintf(MESSAGES,"#Parameters file (host):%s\n",argv[4]);
    fprintf(MESSAGES,"#Parameters file (invader):%s\n",argv[5]);
    fprintf(MESSAGES,"#Parameters file (simulation):%s\n",argv[6]);
    
    fprintf(MESSAGES,"#Event number: %ld \n \n",number_of_iterations);
    fprintf(MESSAGES,"#Specific infos: oxygen damping, eigenvalue shortcut\n");
    
    fprintf(MESSAGES,"#Death rates: %lf %lf \n",1.0/death_rate_inv_host,1.0/death_rate_inv_invader);
    fprintf(MESSAGES,"#tau_p's: %lf %lf\n",tau_p_host,tau_p_invader);
    fprintf(MESSAGES,"#p6/p3: %lf %lf\n",p6_over_p3_host,p6_over_p3_invader);
    fprintf(MESSAGES,"#Diffusion coefficients: %lf %lf %lf\n",diff_coef_host,diff_coef_invader,diff_coef_oxygen);
    fprintf(MESSAGES,"#source_oxygen: %lf\n",source_oxygen);
    fprintf(MESSAGES,"#k_oxygen: %lf\n",k_oxygen);
    fprintf(MESSAGES,"#slot size: %lf mm\n \n \n",delta_x);
     
    fprintf(MESSAGES,"#List of issues:\n");
    //So far "issues" are printed in stderr, to avoid passing *MESSAGES to subroutines
    //Think about this...

    /**/
    //Rescale diffusion coefficients to be fed to the spatial solver
    diff_coef_host=diff_coef_host/(delta_x*delta_x);
    diff_coef_invader=diff_coef_invader/(delta_x*delta_x);
    diff_coef_oxygen=diff_coef_oxygen/(delta_x*delta_x);

    
    
    /*******************************/
    //INITIALIZING DIVISION RATES
    
    //We need to compute the division at zero oxigen for the sake of formulas
    //Just call the ODE solver once
    
    //book division_threshold 
    if((division_threshold_host= (double *) malloc(sizeof(double)*n_xslots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    if((division_threshold_invader= (double *) malloc(sizeof(double)*n_xslots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    //Compute aplus and keep it through the simulation
    if(p6_over_p3_host>r_cr){
        aplus_host=Compute_aplus(p6_over_p3_host);
    };
    
    if(p6_over_p3_invader>r_cr){
        aplus_invader=Compute_aplus(p6_over_p3_invader);
    };
    
    //Initialize division rates
    Compute_division_threshold(division_threshold_host,oxygen_level,p6_over_p3_host,aplus_host,n_xslots-1);
    Compute_division_threshold(division_threshold_invader,oxygen_level,p6_over_p3_invader,aplus_invader,n_xslots-1);
    //This may be a bit less efficient than using the flags. However the code looks cleaner
    
    
    /******************************/
    
    
    /******************************/
    //MAIN LOOP
    
    
    for(i_iteration=1;i_iteration<=number_of_iterations;i_iteration++){
        
        
        //Explicit RK4 time marching.
        
        //Main submodule that:
        //ADVANCEs THE POPULATION
        //ADVANCEs THE OXYGEN
        //RECOMPUTEs SETTINGS FOR NEXT ITERATION

        RK4Handler( 
                delta_t,
                n_xslots,
                number_of_cells_host,
                number_of_cells_invader,
                death_rate_inv_host,
                death_rate_inv_invader,
                division_threshold_host,
                division_threshold_invader,
                tau_p_host,
                tau_p_invader,
                aplus_host,
                aplus_invader,
                p6_over_p3_host,
                p6_over_p3_invader,
                diff_coef_host,
                diff_coef_invader,
                diff_coef_oxygen,
                oxygen_level,
                k_oxygen,
                source_oxygen
                );
        
        
        t=t+delta_t;
        

        //recompute division_threshold
        Compute_division_threshold(division_threshold_host,oxygen_level,p6_over_p3_host,aplus_host,n_xslots-1);
        Compute_division_threshold(division_threshold_invader,oxygen_level,p6_over_p3_invader,aplus_invader,n_xslots-1);
 
        
        /******************/
       //PRINTING DATA ON SELECTED ITERATIONS
       
        if((i_iteration%sampling_period)==0){
            
            output_iter=i_iteration/sampling_period;
            
            //printing total number of cells (both species)
            total_number_of_cells_host=0.0;
            for(x_slot=0;x_slot<n_xslots;x_slot++){
                total_number_of_cells_host=total_number_of_cells_host+number_of_cells_host[x_slot];
            };
            fprintf(POPULATION_host,"%lf %lf \n",t,total_number_of_cells_host);
            
            total_number_of_cells_invader=0.0;
            for(x_slot=0;x_slot<n_xslots;x_slot++){
                total_number_of_cells_invader=total_number_of_cells_invader+number_of_cells_invader[x_slot];
            };
            fprintf(POPULATION_invader,"%lf %lf \n",t,total_number_of_cells_invader);
            
    
            //Printing spatial distribution of cells (both populations)
            strcpy(output_path1h,label1h);
            Print_Vector(OUTPUT_DATA,output_path1h,output_iter,n_xslots,number_of_cells_host);
            
            strcpy(output_path1i,label1i);
            Print_Vector(OUTPUT_DATA,output_path1i,output_iter,n_xslots,number_of_cells_invader);
            
            //Printing oxygen spatial distribution
            strcpy(output_path2,label2);
            Print_Vector(OUTPUT_DATA,output_path2,output_iter,n_xslots,oxygen_level);
        
        }; //end if printing outputs
        
    }; //end main loop

    
    /***************************************/
    

     
    /***************************************/
    
    
    //CLOSING FILES AND RETRIEVING MEMORY
    
    free(number_of_cells_host);
    free(number_of_cells_invader);
    free(total_number_of_cells);
    free(oxygen_level);
    free(division_threshold_host);
    free(division_threshold_invader);
    
     
    fclose(POPULATION_host);
    fclose(POPULATION_invader);
    fclose(MESSAGES);
     
    
	return(0);  
	
}      //end main




