#include <math.h>                           
#include <stdio.h>
#include <stdlib.h>     
#include <time.h>                   

#include "Header.h"

//De momento no tenemos extensiones de memoria programadas


int main (int argc, char **argv){
    
    FILE *INITIAL_POPULATION_HOST, *INITIAL_OXYGEN, *PARAMETERS_HOST, *SIM_DATA;
    FILE *INITIAL_POPULATION_INVADER, *PARAMETERS_INVADER;
    FILE *MESSAGES, *OUTPUT_DATA;

    long i,j,k;    //counters
    
    //various times
    double t = 0.0, tstop;
    double tau=0.0;
    double delta_t;

    double delta_x;

    int n_files;
    double sampling_time_window = 0.0;
    int counter=1; //(int) MAY cause trouble if we want a large number of output files
    
    int hour; //to use the random number generator
    double r1, r2, rtot; //random numbers

    double time_spent=-1; //simulation time spent
    
    long **StochasticNcells_host;
    double **StochasticAge_host;
    long **StochasticNcells_invader;
    double **StochasticAge_invader;
    
    long n_spatial_slots=0;

    double *densidad_n_slots_host,*densidad_n_slots_invader,*densidad_n_slots_total;
    //ACTUALLY those are not densities buth rather cell numbers per slot
    long *nspaces,*nspacesbis;
    
    //population parameters
    double diffusion_R,diffusion_L;
    double tau_p_host,tau_p_invader;
    double death_rate_host,death_rate_invader; //death rate
    double p6_over_p3_host,p6_over_p3_invader;
    
	//parameters in oxygen's equation
    double *oxygen_level;
    double source_oxygen;
    double consumption_oxygen; //k_oxygen in other files
    double suma=0.0;
    
    //Age data structures
    double **b_host,**b_invader;
    long *J_host,*J_invader;
    double **gilAux; //Cambiar por Propensities
    //Global, metemos las de las reacciones de las dos poblaciones
    double *ag1s_host,*ag1s_invader;

    double aplus_host=-1.0,aplus_invader=-1.0;//controls switch at zero oxygen
        //see formula (2) in JCompPhys paper

    char output_path1[200]="OutputValuesHostPopulation";
    char output_path2[200]="OutputValuesInvaderPopulation";
    char output_path3[200]="OutputValuesOxygen";
    char usertag[100]="";
    char label1[200]="OutputValuesHostPopulation";
    char label2[200]="OutputValuesInvaderPopulation";
    char label3[200]="OutputValuesOxygen";
    char label4[200]="messages";
    char create_folder[300]="mkdir ";
    
    //Therapy: survival fraction
    double survival_rate=1.0;
   
    
      /*****************/
    //OPENING FILES and reading inputs
    
    if(argc!=8){
		printf("Error, incorrect argument number\n");
		printf("Use as: \n exec_name  \n initial_host_file \n initial_invader_file \n initial_oxygen \n host_parameters_file \n invader_parameters_file \n simulation_info \n output_files_tag \n");
		return(1);
	};
    
    //Open age distribution file 
    if((INITIAL_POPULATION_HOST=fopen(argv[1],"rt"))==NULL){
        printf("Error: initial_host_file could not be opened \n");
        return(1) ;
    };
    
    if((INITIAL_POPULATION_INVADER=fopen(argv[2],"rt"))==NULL){
        printf("Error: initial_invader_file could not be opened \n");
        return(1) ;
    };
    
    if((INITIAL_OXYGEN=fopen(argv[3],"rt"))==NULL){
        printf("Error: STOCHASTIC_file could not be opened \n");
        return(1) ;
    };

    if((PARAMETERS_HOST=fopen(argv[4],"rt"))==NULL){
		printf("Error: host_parameters_file could not be opened \n");
		return(1) ;
	};
    
    if((PARAMETERS_INVADER=fopen(argv[5],"rt"))==NULL){
		printf("Error: invader_parameters_file could not be opened \n");
		return(1) ;
	};

	if((SIM_DATA=fopen(argv[6],"rt"))==NULL){
		printf("Error: parameters_file could not be opened \n");
		return(1) ;
	};

    strcpy(usertag,argv[7]);
    /*****************************/
    //Read initial age distribution
    
    n_spatial_slots=Read_Init_Population(
                                         INITIAL_POPULATION_HOST,
                                         &StochasticNcells_host,
                                         &StochasticAge_host,
                                         &nspaces
                                         );
    
    
    //Doing the stuff for the INVADER POPULATION
    
    if(Read_Init_Population(
                            INITIAL_POPULATION_INVADER,
                            &StochasticNcells_invader,
                            &StochasticAge_invader,
                            &nspacesbis
                            )
       !=n_spatial_slots){
        fprintf(stderr,"Mismatched number of lines, aborting execution\n");
        exit(1);
    };
    
    
    /*****************************/
    //Read the rest of the stuff
    //Read parameters 
    Read_params_population(PARAMETERS_HOST,
                            &death_rate_host,
                            &tau_p_host,
                            &p6_over_p3_host);
    
    Read_params_population(PARAMETERS_INVADER,
                           &death_rate_invader,
                           &tau_p_invader,
                           &p6_over_p3_invader);
    
    //Read parameters (simulation infos)
    Read_params_sim(SIM_DATA,
                        &source_oxygen,
                         &consumption_oxygen,
                         &delta_t,
                         &delta_x,
                         &tstop,
                         &n_files,
                         &survival_rate);
    
    //printf("Survival rate: %lf\n",survival_rate);
    
    sampling_time_window=tstop/n_files;

    //Read initial oxygen
    if (n_spatial_slots!= Read_Init_Space_Oxygen(INITIAL_OXYGEN, &oxygen_level)){
        fprintf(stderr, "Mismatched number of spactial cells, aborting execution \n");
        exit(1);
    };

	/********/
    //Setting age data structures
    
    J_host=(long *)calloc(n_spatial_slots,sizeof(long));
    for (i=0;i<n_spatial_slots;i++){
        if((nspaces[2*i]==1)&&(StochasticAge_host[i][0]==0.0)){J_host[i]=0;}
        //empty slot situation, according to input format convention
        else{J_host[i]=nspaces[2*i];}
    };
    
    J_invader=(long *)calloc(n_spatial_slots,sizeof(long));
    for (i=0;i<n_spatial_slots;i++){
        if((nspacesbis[2*i]==1)&&(StochasticAge_invader[i][0]==0.0)){J_invader[i]=0;}
        //empty slot situation, according to input format convention
        else{J_invader[i]=nspacesbis[2*i];}
    };
    
    b_host=(double **)malloc(n_spatial_slots*sizeof(double *));  //birth rate
    for(i=0;i<n_spatial_slots;i++){b_host[i]= (double*)calloc(MEMORY_BATCH,sizeof(double));}
    
    b_invader=(double **)malloc(n_spatial_slots*sizeof(double *));  //birth rate
    for(i=0;i<n_spatial_slots;i++){b_invader[i]= (double*)calloc(MEMORY_BATCH,sizeof(double));}
    
    gilAux=(double **)malloc(n_spatial_slots*sizeof(double *));
    for(i=0;i<n_spatial_slots;i++){
        gilAux[i]= (double*)calloc(2*NUMBER_OF_REACTIONS*MEMORY_BATCH,sizeof(double));
        //We set twice NUMBER_OF_REACTIONS since we want to handle two populations
    }
    
    ag1s_host=(double *)calloc(n_spatial_slots,sizeof(double));
    ag1s_invader=(double *)calloc(n_spatial_slots,sizeof(double));
    
    //Remaining stuff to set before we start
    densidad_n_slots_host=(double *)calloc(n_spatial_slots,sizeof(double));
    densidad_n_slots_invader=(double *)calloc(n_spatial_slots,sizeof(double));
    densidad_n_slots_total=(double *)calloc(n_spatial_slots,sizeof(double));
    diffusion_R=DIFFUSION_COEF_CELLS/(delta_x*delta_x); 
    diffusion_L=DIFFUSION_COEF_CELLS/(delta_x*delta_x);

    
    
    Compute_density(n_spatial_slots,
                    J_host,
                    StochasticNcells_host,
                    densidad_n_slots_host);
    
    Compute_density(n_spatial_slots,
                    J_invader,
                    StochasticNcells_invader,
                    densidad_n_slots_invader);

    if(p6_over_p3_host>r_cr){aplus_host=Compute_aplus(p6_over_p3_host);};
    if(p6_over_p3_invader>r_cr){aplus_invader=Compute_aplus(p6_over_p3_invader);};
    
    for (i=0; i<n_spatial_slots;i++){
        densidad_n_slots_total[i]=densidad_n_slots_host[i]+densidad_n_slots_invader[i];
    };

    //seeding the random number generator
    hour = time(NULL);
    srand48(hour);

    
    //print initial conditions
    //create "OutputValuesOxygen" (plus user tag), etc
    strcat(label1,usertag);
    strcat(label2,usertag);
    strcat(label3,usertag);
    
    strcat(create_folder,label1);
    system(create_folder);
    strcpy(create_folder,"mkdir ");
    strcat(create_folder,label2);
    system(create_folder);
    strcpy(create_folder,"mkdir ");
    strcat(create_folder,label3);
    system(create_folder);
    
    
    //Printing spatial NUMBER of cells regardless of their age
    strcat(label1,"/Out"); //Final label for population folder 1
    strcpy(output_path1,label1);
    Print_Vector(OUTPUT_DATA,output_path1,0,n_spatial_slots,densidad_n_slots_host);
    strcat(label2,"/Out"); //Final label for oxygen folder
    strcpy(output_path2,label2);
    Print_Vector(OUTPUT_DATA,output_path2,0,n_spatial_slots,densidad_n_slots_invader);
    
    //Printing oxygen spatial distribution
    strcat(label3,"/Out"); //Final label for population folder 1
    strcpy(output_path3,label3);
    Print_Vector(OUTPUT_DATA,output_path3,0,n_spatial_slots,oxygen_level);

    
    

    /*************/
    //Print simulation infos
    strcat(label4,usertag);
    if((MESSAGES=fopen(label4,"w"))==NULL){ //Master output file
		fprintf(stderr,"Error: output messages file could not be created \n");
		return(1) ;
	};

    fprintf(MESSAGES,"#Debug info: \n \n");
    fprintf(MESSAGES,"#Data generated with the following .exe file:%s\n",argv[0]);
    fprintf(MESSAGES,"#Initial condition file (host, age):%s\n",argv[1]);
    fprintf(MESSAGES,"#Initial condition file (invader, age):%s\n",argv[2]);
    fprintf(MESSAGES,"#Initial condition file (oxygen):%s\n",argv[3]);
    fprintf(MESSAGES,"#Host population parameters file:%s\n",argv[4]);
    fprintf(MESSAGES,"#Invader population parameters file:%s\n",argv[5]);
    fprintf(MESSAGES,"#Parameters file (simulation):%s\n",argv[6]);
    
    fprintf(MESSAGES,"#Death rates (host,invader): %.10lf %.10lf\n",death_rate_host,death_rate_invader);
    fprintf(MESSAGES,"#tau_p's (host,invader): %.10lf %.10lf \n",tau_p_host,tau_p_invader);
    fprintf(MESSAGES,"#p6/p3's (host,invader): %.10lf %.10lf\n",p6_over_p3_host,p6_over_p3_invader);
    fprintf(MESSAGES,"#Diffusion coefficient cells: %.10lf \n",DIFFUSION_COEF_CELLS/(delta_x*delta_x));
    fprintf(MESSAGES,"#Diffusion coefficient oxygen: %.10lf\n",DIFFUSION_COEF_O2/(delta_x*delta_x));
    fprintf(MESSAGES,"#source_oxygen: %.10lf\n",source_oxygen);
    fprintf(MESSAGES,"#consumption_oxygen: %.10lf\n",consumption_oxygen);
    fprintf(MESSAGES,"#slot size: %.10lf mm\n",delta_x);
    fprintf(MESSAGES,"#survival fraction: %lf \n \n \n",survival_rate);

    
    
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////Start time step/////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
    
    clock_t begin = clock();
    while (t<tstop){ //main loop
        suma=Reactions_two_pop(n_spatial_slots,
                               gilAux,
                               J_host,
                               StochasticNcells_host,
                               b_host,
                               J_invader,
                               StochasticNcells_invader,
                               b_invader,
                               death_rate_host,
                               death_rate_invader,
                               diffusion_R,
                               diffusion_L,
                               survival_rate
                               );
        
        r1=drand48();
        r2=drand48();
        tau=log(1/r1)/suma;
        
        //print the last state prior surpassing the time mark
        if(t+tau>=counter*sampling_time_window){
    
            //Printing spatial number of cells
            strcpy(output_path1,label1);
            Print_Vector(OUTPUT_DATA,output_path1,counter,n_spatial_slots,densidad_n_slots_host);
            strcpy(output_path2,label2);
            Print_Vector(OUTPUT_DATA,output_path2,counter,n_spatial_slots,densidad_n_slots_invader);
            
            //Printing oxygen spatial distribution
            strcpy(output_path3,label3);
            Print_Vector(OUTPUT_DATA,output_path3,counter,n_spatial_slots,oxygen_level);
            
            //printf("t=%lf\n",counter*sampling_time_window);
        	
        	counter++;
        };//if() for printing stuff,
        //update system state right after
        
        
        //ESto deberia ir antes que la actualizacion de la poblacion???
        OxygenSolverDamp(oxygen_level,
                         n_spatial_slots,
                         densidad_n_slots_total, //////SUM of the two populations
                         tau,
                         consumption_oxygen,
                         source_oxygen,
                         delta_x);
        
        rtot=suma*r2;
        suma=0.0;
        i=-1;
        while(suma<=rtot)
        {
            i=i+1;
            for(j=0;j<NUMBER_OF_REACTIONS*(J_host[i]+J_invader[i]);j++)
            {
                suma=suma+gilAux[i][j];
                if(suma>rtot){break;};
            };
        };
        
        
        if(j<NUMBER_OF_REACTIONS*J_host[i]){
            Gillespie(i,
                       j,
                       StochasticNcells_host,
                       J_host,
                       StochasticAge_host,
                       b_host);
        }
        else{
            Gillespie(i,
                      j-NUMBER_OF_REACTIONS*J_host[i],
                      StochasticNcells_invader,
                      J_invader,
                      StochasticAge_invader,
                      b_invader);
        };
        
        
        for (k=0;k<n_spatial_slots;k++){
            ag1s_host[k]=Get_division_threshold(oxygen_level[k],
                                           p6_over_p3_host,
                                           aplus_host);

            for (i=0;i<J_host[k];i=i+1){
                if(oxygen_level[k]>c_cr(p6_over_p3_host))
                    //Esta condicion vienen a ser meter las celulas en quiescencia (17/05/2018)
                    //Corroborado con Tomas
                {
                    StochasticAge_host[k][i]=StochasticAge_host[k][i]+tau;
                };
                //birth condition (X<threshlod): 
                if(StochasticAge_host[k][i]>ag1s_host[k]){
                    b_host[k][i]=1./tau_p_host;
                }
                else{
                    b_host[k][i]=0;
                };
            }; //end for i
            
            
            ag1s_invader[k]=Get_division_threshold(oxygen_level[k],
                                                   p6_over_p3_invader,
                                                   aplus_invader);
            
            for (i=0;i<J_invader[k];i=i+1){
                if(oxygen_level[k]>c_cr(p6_over_p3_invader)){
                    StochasticAge_invader[k][i]=StochasticAge_invader[k][i]+tau;
                };
                //birth condition (X<threshlod):
                if(StochasticAge_invader[k][i]>ag1s_invader[k]){
                    b_invader[k][i]=1./tau_p_invader;
                }
                else{
                    b_invader[k][i]=0;
                };
            }; //end for i
        }; //end for k
        
        

        if(Compute_density(n_spatial_slots,
                           J_host,
                           StochasticNcells_host,
                           densidad_n_slots_host)
                    <=0.0){break;};
        
        if(Compute_density(n_spatial_slots,
                           J_invader,
                           StochasticNcells_invader,
                           densidad_n_slots_invader)
                    <=0.0){break;};
        
        
        
        for(i=0;i<n_spatial_slots;i++){
            densidad_n_slots_total[i]=densidad_n_slots_invader[i]+densidad_n_slots_host[i];
        };
        
        t+=tau;
        
    }; //end while main loop

    clock_t end = clock();
    time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
    fprintf(stdout,"%lf \n",time_spent);
    //printf("#Tiempo final=%lf \n",t);

    fprintf(MESSAGES,"#end time: %lf \n \n \n",t);
    fprintf(MESSAGES,"#simulation time spent (sec): %lf \n \n \n",time_spent);


    fclose(MESSAGES);    /////////////////////////////////////////
    //free memory all matrix and vectors//////  
    //////////////////////////////////////////
    
    free(nspaces);
    free(nspacesbis);
    for(i=0;i<n_spatial_slots;i++){
        free(StochasticNcells_host[i]);
        free(StochasticAge_host[i]);
        free(StochasticNcells_invader[i]);
        free(StochasticAge_invader[i]);
    };
    free(StochasticNcells_host);
    free(StochasticAge_host);
    free(StochasticNcells_invader);
    free(StochasticAge_invader);
    
    for(k=0;k<n_spatial_slots;k++){
        free(b_host[k]);
        free(gilAux[k]);
        free(b_invader[k]);
    }

    free(J_host);
    free(J_invader);
    free(ag1s_host);
    free(b_host);
    free(ag1s_invader);
    free(b_invader);

    free(gilAux);
    free(densidad_n_slots_host);
    free(densidad_n_slots_invader);
    free(densidad_n_slots_total);
    free(oxygen_level);
    
    return 0;
} 