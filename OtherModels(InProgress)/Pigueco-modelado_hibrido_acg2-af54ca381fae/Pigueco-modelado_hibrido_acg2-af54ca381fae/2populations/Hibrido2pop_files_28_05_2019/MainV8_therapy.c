

//CODE FOR TWO COMPETING WAVES
//no other situation will work fine (interfaces are managed in a specific way 
//that is only valid for this specific case)


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
    double delta_t; //RECORDAR en que uds se mide esto
    
    int hour; //to use the random number generator
    double r1, r2, rtot; //random numbers

    double time_spent=-1; //simulation time spent
    
    long **StochasticNcells_host; //matrix where row encodes spatial location 
                        //and column encodes number of cells with a given age (to be read from the matrix 
                        //StochasticAge_host, whose structure coresponds with this)
    double **StochasticAge_host; //matrix where row encodes spatial location 
                        //and column encodes a given age that shows up at that slot (the number of cells 
                        //with that age is to be read from StochasticNcells_host)
    long **StochasticNcells_invader;
    double **StochasticAge_invader;

    double *number_of_cells_host; //number of cells at each spatial slot, regardless of their age
    double *number_of_births_host;
    double *density_of_cells_host;
    double *number_of_cells_invader;
    double *number_of_births_invader;
    double *density_of_cells_invader;
    double *total_number_of_cells;
    
    long n_spatial_slots=0;
    double delta_x; //APARENTEMENTE medido en milimetros


    long *nspaces,*nspacesbis;  //WHAT FOR?

    //population parameters
    double diffusion_host,diffusion_invader;
    double tau_p_host,tau_p_invader;
    double death_rate_host,death_rate_invader; //death rate
    double p6_over_p3_host,p6_over_p3_invader;
    
	//parameters in oxygen's equation
    double *oxygen_level;
    double diffusion_oxygen;
    double source_oxygen;
    double consumption_oxygen; //k_oxygen in other files


    double total_propensity=0.0;
    double suma_aux=0.0;
    double birth_aux;

    //Therapy: survival fraction
    double survival_rate=1.0;
    
    //Age data structures
    double **b_host,**b_invader;
    long *J_host,*J_invader;
    double **Propensities_host,**Propensities_invader;//**gilAux; //Cambiar por Propensities
                //Global, metemos las de las reacciones de las dos poblaciones
    double *ag1s_host,*ag1s_invader;

    double aplus_host=-1.0,aplus_invader=-1.0;//controls switch at zero oxygen
        //see formula (2) in JCompPhys paper

    //parametros necesarios para manejar las interfaces del hibrido
    long Interface_location_host, aux_interface_host;
    long Interface_location_invader, aux_interface_invader;
    double cells_at_interface_host;
    double cells_at_interface_invader;

    //Stuff to handle output files
    int n_files;
    double sampling_time_window = 0.0;
    int counter=1; //(int) MAY cause trouble if we want a large number of output files
    char output_path1[200]="OutputValuesHostPopulation";
    char output_path2[200]="OutputValuesInvaderPopulation";
    char output_path3[200]="OutputValuesOxygen";
    char output_path4[200]="OutputBirthsInvaderPopulation";
    char output_path5[200]="OutputBirthsHostPopulation";
    char usertag[100]="";
    char label4[200]="messages";
    char create_folder[300]="mkdir ";
    
    
   
    /******************/
    /**********INITIALIZATION AND PROCESSING***********/

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
    //READ initial age distribution
    
    n_spatial_slots=Read_Init_Population(
                                         INITIAL_POPULATION_HOST,
                                         &StochasticNcells_host,
                                         &StochasticAge_host,
                                         &nspaces
                                         );    
    //Doing the same stuff for the INVADER POPULATION
    if(n_spatial_slots!=Read_Init_Population(
                                            INITIAL_POPULATION_INVADER,
                                            &StochasticNcells_invader,
                                            &StochasticAge_invader,
                                            &nspacesbis
                                            )
       ){
        fprintf(stderr,"Mismatched number of lines, aborting execution\n");
        exit(1);
    };
    
    
    /*****************************/
    //READ the rest of the stuff

    //Read host parameters 
    Read_params_population(PARAMETERS_HOST,
                            &death_rate_host,
                            &tau_p_host,
                            &p6_over_p3_host);
    
    Read_params_population(PARAMETERS_INVADER,
                           &death_rate_invader,
                           &tau_p_invader,
                           &p6_over_p3_invader);
    
    //Read simulation parameters 
    Read_params_sim(SIM_DATA,
                        &source_oxygen,
                         &consumption_oxygen,
                         &delta_t,
                         &delta_x,
                         &tstop,
                         &n_files,
                         &survival_rate);
    
    //Read initial oxygen
    if (n_spatial_slots!= Read_Init_Space_Oxygen(INITIAL_OXYGEN, &oxygen_level)){
        fprintf(stderr, "Mismatched number of spactial cells, aborting execution \n");
        exit(1);
    };


	/********/
    //SETTING AGE data structures
    
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
    
    b_host=(double **)calloc(n_spatial_slots,sizeof(double *));  //birth rate
    for(i=0;i<n_spatial_slots;i++){b_host[i]= (double*)calloc(MEMORY_BATCH,sizeof(double));}
    
    b_invader=(double **)calloc(n_spatial_slots,sizeof(double *));  //birth rate
    for(i=0;i<n_spatial_slots;i++){b_invader[i]= (double*)calloc(MEMORY_BATCH,sizeof(double));}
    
    Propensities_host=(double **)calloc(n_spatial_slots,sizeof(double *));
    for(i=0;i<n_spatial_slots;i++){
        Propensities_host[i]= (double*)calloc(2*NUMBER_OF_REACTIONS*MEMORY_BATCH,sizeof(double));
        //We set twice NUMBER_OF_REACTIONS since we want to handle two populations
    };

    Propensities_invader=(double **)calloc(n_spatial_slots,sizeof(double *));
    for(i=0;i<n_spatial_slots;i++){
        Propensities_invader[i]= (double*)calloc(2*NUMBER_OF_REACTIONS*MEMORY_BATCH,sizeof(double));
        //We set twice NUMBER_OF_REACTIONS since we want to handle two populations
    };
    
    ag1s_host=(double *)calloc(n_spatial_slots,sizeof(double)); //division_threshold_host en otras versiones
    ag1s_invader=(double *)calloc(n_spatial_slots,sizeof(double));
//NO ENTIENDO donde se inicializan estos vectores (o quiza es que usemos banderas)

    if(p6_over_p3_host>r_cr){aplus_host=Compute_aplus(p6_over_p3_host);};
    if(p6_over_p3_invader>r_cr){aplus_invader=Compute_aplus(p6_over_p3_invader);};



    /****************/
    //SETTING additional structures and simulation parameters
    
    number_of_cells_host=(double *)calloc(n_spatial_slots,sizeof(double));
    number_of_births_host=(double *)calloc(n_spatial_slots,sizeof(double));
    density_of_cells_host=(double *)calloc(n_spatial_slots,sizeof(double));
    number_of_cells_invader=(double *)calloc(n_spatial_slots,sizeof(double));
    number_of_births_invader=(double *)calloc(n_spatial_slots,sizeof(double));
    density_of_cells_invader=(double *)calloc(n_spatial_slots,sizeof(double));
    total_number_of_cells=(double *)calloc(n_spatial_slots,sizeof(double));

    diffusion_host=DIFFUSION_COEF_CELLS/(delta_x*delta_x); 
    diffusion_invader=DIFFUSION_COEF_CELLS/(delta_x*delta_x);
    diffusion_oxygen=DIFFUSION_COEF_O2/(delta_x*delta_x);

    
    //we extract here the number of cells per spatial slot
    Compute_ncells(n_spatial_slots,
                    J_host,
                    StochasticNcells_host,
                    number_of_cells_host);
    
    Compute_ncells(n_spatial_slots,
                    J_invader,
                    StochasticNcells_invader,
                    number_of_cells_invader);
    
    for (i=0; i<n_spatial_slots;i++){
        total_number_of_cells[i]=number_of_cells_host[i]+number_of_cells_invader[i];
    };

    //now computing the associated densities
    for (i=0; i<n_spatial_slots;i++){
        density_of_cells_host[i]=number_of_cells_host[i]/delta_x;
        density_of_cells_invader[i]=number_of_cells_invader[i]/delta_x;
    };


    sampling_time_window=tstop/n_files;

    //seeding the random number generator
    hour = time(NULL);
    srand48(hour);

    /***********/

    //COMPUTE SYSTEMS INITIAL INTERFACE(s)
    //Say that we have two competing waves. Each population has a single interface, and the meaning
    //is different (one decreases to the left, the other to the right) 

    //for the population that decreases to the right (say the invader)
    Interface_location_invader=0;
    while((number_of_cells_invader[Interface_location_invader]>=POPULATION_THRESHOLD)&&(Interface_location_invader<n_spatial_slots-1)){
        Interface_location_invader++;
    }; //Vector index, nor the natural one
    
    Interface_location_invader=Interface_location_invader-1; //Sets the leftmost interface slot at the last slot above threshold
    // printf("Invader interface: %ld\n",Interface_location_invader);

    //for the population that decreases to the left (say the host)
    Interface_location_host=n_spatial_slots-1;
    while((number_of_cells_host[Interface_location_host]>=POPULATION_THRESHOLD)&&(Interface_location_host>0)){
        Interface_location_host--;
    }; //Vector index, nor the natural one
    
    Interface_location_host=Interface_location_host+1; //Sets the leftmost interface slot at the last slot above threshold 
     //(last is actually first if we look from left to right)
    // printf("Host interface: %ld\n",Interface_location_host);

    //Cases badly covered by the previous procedure:
    //i) a fully deterministic population
    //then it is n_spatial_slots-2 (next to the rightmost slot) for the invader
    //and 1 for the host (next to the leftmost slot)

    //ii) a fully stochastic population (THIS ONE IS DANGEROUS)
    //this is -1 (out of range) for the invader
    //and n_spatial_slots for the host (out of range as well)

    

    /***********/
    //CREATE output folders and paths
    strcat(output_path1,usertag);
    strcat(output_path2,usertag);
    strcat(output_path3,usertag);
    strcat(output_path4,usertag);
    strcat(output_path5,usertag);
    
    strcat(create_folder,output_path1);
    system(create_folder);
    strcpy(create_folder,"mkdir ");
    strcat(create_folder,output_path2);
    system(create_folder);
    strcpy(create_folder,"mkdir ");
    strcat(create_folder,output_path3);
    system(create_folder);
    strcpy(create_folder,"mkdir ");
    strcat(create_folder,output_path4);
    system(create_folder);
    strcpy(create_folder,"mkdir ");
    strcat(create_folder,output_path5);
    system(create_folder);

    strcat(label4,usertag);
    if((MESSAGES=fopen(label4,"w"))==NULL){ //Master output file
        fprintf(stderr,"Error: output messages file could not be created \n");
        return(1) ;
    };

    
    /***********/
    //PRINT initial infos
    //printf("Survival rate: %lf\n",survival_rate);

    //print initial conditions
    //Printing spatial NUMBER of cells regardless of their age
    strcat(output_path1,"/Out"); //Final label for population folder 1
    Print_Vector(OUTPUT_DATA,output_path1,0,n_spatial_slots,number_of_cells_host);
    strcat(output_path2,"/Out"); //Final label for population folder 2
    Print_Vector(OUTPUT_DATA,output_path2,0,n_spatial_slots,number_of_cells_invader);
    strcat(output_path4,"/Out"); //Final label for population folder 4
    Print_Vector(OUTPUT_DATA,output_path4,0,n_spatial_slots,number_of_births_invader);
    strcat(output_path5,"/Out"); //Final label for population folder 5
    Print_Vector(OUTPUT_DATA,output_path5,0,n_spatial_slots,number_of_births_host);
    

    //Printing oxygen spatial distribution
    strcat(output_path3,"/Out"); //Final label for oxygen folder
    Print_Vector(OUTPUT_DATA,output_path3,0,n_spatial_slots,oxygen_level);



    /*************/
    //PRINT simulation infos

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
    fprintf(MESSAGES,"#Diffusion coefficient cells (host, invader): %.10lf %.10lf \n",diffusion_host,diffusion_invader);
    fprintf(MESSAGES,"#Diffusion coefficient oxygen: %.10lf\n",diffusion_oxygen);
    fprintf(MESSAGES,"#source_oxygen: %.10lf\n",source_oxygen);
    fprintf(MESSAGES,"#consumption_oxygen: %.10lf\n",consumption_oxygen);
    fprintf(MESSAGES,"#slot size: %.10lf mm\n",delta_x);
    fprintf(MESSAGES,"#survival fraction: %lf \n \n \n",survival_rate);

    ///ag1s inicialisation 
    for (k=0;k<n_spatial_slots;k++){   
            
            ag1s_invader[k]=Get_division_threshold(oxygen_level[k],
                                                   p6_over_p3_invader,
                                                   aplus_invader);
            ag1s_host[k]=Get_division_threshold(oxygen_level[k],
                                                   p6_over_p3_host,
                                                   aplus_host);
        } // end for k
    
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////Start time step/////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
    
    clock_t begin = clock();
    while (t<tstop){ //main loop

        /*****************/
        //GILLESPIE 1: compute propensities, firing channel and waiting time
//Es necesario pasar la localizacion de las interfaces para 
//(aunque conozcamos el numero de celulas en todos los voxels)
//solo meter en el computo las secciones de baja densidad
    
        total_propensity=Reactions_invader(n_spatial_slots,
                               Propensities_invader,
                               J_invader,
                               StochasticNcells_invader,
                               b_invader,
                               Interface_location_invader,
                               death_rate_invader,
                               diffusion_invader,
                               survival_rate
                               );

        //sum of propensities for the population at the right and total sum
        total_propensity=total_propensity+Reactions_host(n_spatial_slots,
                               Propensities_host,
                               J_host,
                               StochasticNcells_host,
                               b_host,
                               Interface_location_host,
                               death_rate_host,
                               diffusion_host,
                               survival_rate);
        // printf("suma total=%lf\n",total_propensity);
        r1=drand48();
        r2=drand48();
        tau=log(1/r1)/total_propensity;
   
        /**************/
        //SHOWING OUTPUTS on selected iterations

        //print the last state prior surpassing the time mark
        if(t+tau>=counter*sampling_time_window){

            //printing birth events (only within the stochastic domains and invader pop)
            for (k=Interface_location_invader;k<n_spatial_slots;k++){   
                birth_aux=0;
                for (i=0;i<J_invader[k];i=i+1){
                    //birth condition (X<threshlod):
                    if(StochasticAge_invader[k][i]>ag1s_invader[k]){
                        birth_aux = birth_aux + StochasticNcells_invader[k][i];
                    }
                    
                }; //end for i
                number_of_births_invader[k]=birth_aux;
            }; //end for k

            //printing birth events (only within the stochastic domains and host pop)
            for (k=0;k<=Interface_location_host;k++){   
                birth_aux=0;
                for (i=0;i<J_host[k];i=i+1){
                    //birth condition (X<threshlod):
                    if(StochasticAge_host[k][i]>ag1s_host[k]){
                        birth_aux = birth_aux + StochasticNcells_host[k][i];
                    }
                    
                }; //end for i
                number_of_births_host[k]=birth_aux;
            }; //end for k
    
            //Printing spatial number of cells
            Print_Vector(OUTPUT_DATA,output_path1,counter,n_spatial_slots,number_of_cells_host);
            Print_Vector(OUTPUT_DATA,output_path2,counter,n_spatial_slots,number_of_cells_invader);
            
            //Printing oxygen spatial distribution
            Print_Vector(OUTPUT_DATA,output_path3,counter,n_spatial_slots,oxygen_level);
            
            //Printing number of births at the last step 
            Print_Vector(OUTPUT_DATA,output_path4,counter,n_spatial_slots,number_of_births_invader);
            Print_Vector(OUTPUT_DATA,output_path5,counter,n_spatial_slots,number_of_births_host);


            //printf("t=%lf, counter= %ld\n",counter*sampling_time_window,counter);
        	
        	counter++;
        };//if() for printing stuff,



        /******************/
        //PERFORMING TEMPORAL EVOLUTION

        RK4HandlerOxy(tau,
                n_spatial_slots,
                number_of_cells_host,
                number_of_cells_invader,
                diffusion_oxygen,
                oxygen_level,
                consumption_oxygen,
                source_oxygen
                );
        //////TRY TO SUPPLY the SUM of the two populations instead //total_number_of_cells CORRECTO!!!!


//PUNTO IMPORTANTE: segun el articulo de 2017, en la ec para el oxigeno se usa el numero de celulas 
//para el termino de consumo. Tal como hacemos aqui. 
//Lo que no entiendo entonces es que hicimos en la version del campo medio... (salvaguarda via delta_x=1)
                         
//Tras avanzar el oxigeno hacemos lo propio con la poblacion.
//El esquema ideal seria:
//1)haz el update de la parte estocastica (incluyendo estructuras de edad y umbrales de division)

        //ESTO QUE SIGUE se supone que hace el update de la poblacion en el codigo estocastico puro
        //PASO 1

        rtot=total_propensity*r2;
        //printf("random number %lf, rtot %lf \n",r2,rtot);
        suma_aux=0.0;

        //Check whether reaction fired in the invader stochastic domain
        //Given the case, update the stuff
        i=Interface_location_invader-1;
        while((suma_aux<=rtot)&&(i<n_spatial_slots-1))
        {
            i=i+1;
            for(j=0;j<NUMBER_OF_REACTIONS*J_invader[i];j++)
            {
                //printf("suma=%lf\n",suma);
                suma_aux=suma_aux+Propensities_invader[i][j];
                if(suma_aux>rtot){
                        // printf("invader, slot %ld, channel %ld\n",i,j);
                    Gillespie(i,
                      j,
                      StochasticNcells_invader,
                      J_invader,
                      StochasticAge_invader,
                      b_invader);
                        //update everything
                        //and then
                    break;};
            };
        };
        //printf("suma_aux tras salir del bucle invader %lf\n",suma_aux);

        //If we entrer this part, the reaction did fire in the host stochastic domain
        //We locate it and update the stuff
        i=Interface_location_host+1; //to check
        while(suma_aux<=rtot)//&&(control space) //We drop this because it should be authomatically satisfied
        {
            i=i-1;
            for(j=0;j<NUMBER_OF_REACTIONS*J_host[i];j++)
            {
                suma_aux=suma_aux+Propensities_host[i][j];
                if(suma_aux>rtot){
                    // printf("host, slot %ld, channel %ld\n",i,j);
                    Gillespie(i,
                       j,
                       StochasticNcells_host,
                       J_host,
                       StochasticAge_host,
                       b_host);
                        //update everything
                        //and then
                    break;};
            };
        };

        
        /****************************/
        
        
        //SAFETY TEST for negative population values
        
        if(Compute_ncells(n_spatial_slots,
                           J_host,
                           StochasticNcells_host,
                           number_of_cells_host
                           )
                            <=0.0){printf( "Host pop 0\n" ); break;};
        
        if(Compute_ncells(n_spatial_slots,
                           J_invader,
                           StochasticNcells_invader,
                           number_of_cells_invader
                           )
                            <=0.0){printf("Invader pop 0\n" ); break;};
        
        
        for(i=0;i<n_spatial_slots;i++){
            total_number_of_cells[i]=number_of_cells_invader[i]+number_of_cells_host[i];
        };
        
        t+=tau;
        

        /****************/
        //2)Calcula flujo de estocastico a determinista (¡facil!)
        //Esto quiere decir, si el firing channel cayo sobre alguna de las interfaces, se recalcula la densidad en esa casill
        
        /**/
        //PASO 2, flujo estocastico a determinista

        //Mas burro pero igual de efectivo que la cascada de condicionales
        density_of_cells_host[Interface_location_host]=number_of_cells_host[Interface_location_host]/delta_x;
        density_of_cells_invader[Interface_location_invader]=number_of_cells_invader[Interface_location_invader]/delta_x;
        cells_at_interface_host=number_of_cells_host[Interface_location_host];
        cells_at_interface_invader=number_of_cells_invader[Interface_location_invader];
        // printf("%lf %ld\n",lrint(0.3),lrint(0.3) );


        /******/
        //PASO 3, deterministic evolution
        //here we feed density_of_cells_host to the subroutines,
        //then convert it again to cell number globally
        //3)Se calcula la evolucion temporal en las mean field regions
        //Como las poblaciones solo estan acopladas a traves de la ec. para el oxigeno, que ya se trato antes,
        //basta aplicar los PDE solvers a cada mean field region (idealmente solo una por cada poblacion) por separado
        //RECUERDA: hay que aplicar el survival rate en el calculo del eigenvalue

        
        RK4HandlerPop(tau,
                n_spatial_slots,
                Interface_location_host,
                n_spatial_slots-1,
                density_of_cells_host,
                death_rate_host,
                ag1s_host,
                tau_p_host,
                diffusion_host,
                survival_rate
                );
        RK4HandlerPop(tau,
                n_spatial_slots,
                0,
                Interface_location_invader,
                density_of_cells_invader,
                death_rate_invader,
                ag1s_invader,
                tau_p_invader,
                diffusion_invader,
                survival_rate
                );




        //4)Calcula flujo de determinista a estocastico
        //primero aprovechamos para convertir los vecorres de densidad a particle numbers


        //CONVERT density_of_cells_invader TO Number_of_cells_invader BEFORE INVOKING Renormalize_center
        //Same with host
        for(i=Interface_location_host;i<n_spatial_slots;i++){
            number_of_cells_host[i]=density_of_cells_host[i]*delta_x;
        };

        for(i=0;i<=Interface_location_invader;i++){
            number_of_cells_invader[i]=density_of_cells_invader[i]*delta_x;
        };



        //segundo miramos los particle numbers que han resultado sobre cada interfaz
        //si no son numeros enteros se aplica el procedimiento de renormalizacion, sobre los vectores de numero de particulas
        //y tras concluir lo convertimos todo a densidades otra vez
        //YA que estamos gestionar las extensines de memoria ANTES de entrar en Renormalize_center

        Renormalize_center(number_of_cells_host,
                            Interface_location_host,
                            Interface_location_host,
                            n_spatial_slots-1,
                            ag1s_host, 
                            death_rate_host, 
                            tau_p_host,
                            StochasticNcells_host,
                            J_host,
                            b_host,
                            StochasticAge_host,
                            survival_rate,
                            diffusion_host,
                            cells_at_interface_host
                            );

         Renormalize_center(number_of_cells_invader,
                            Interface_location_invader,
                            0,
                            Interface_location_invader,
                            ag1s_invader, 
                            death_rate_invader, 
                            tau_p_invader,
                            StochasticNcells_invader,
                            J_invader,
                            b_invader,
                            StochasticAge_invader,
                            survival_rate,
                            diffusion_invader,
                            cells_at_interface_invader
                            );
        

        
        //Now we recompute again cell densities (as masses may have been renormalized)
        
        for(i=0;i<n_spatial_slots;i++){ //Notice that we are doing this globally, check whether this is the sensible thing to do
            density_of_cells_host[i]=number_of_cells_host[i]/delta_x;
            density_of_cells_invader[i]=number_of_cells_invader[i]/delta_x;
        };


        //5)Recalcular la posicion de las interfaces (la parte mas delicada)
        //for the population that decreases to the right (say the invader)
        aux_interface_invader=0;
        
        while((number_of_cells_invader[aux_interface_invader]>=POPULATION_THRESHOLD)&&(aux_interface_invader<n_spatial_slots-1)){
            aux_interface_invader++;
        }; //Vector index, nor the natural one
    
        aux_interface_invader=aux_interface_invader-1; //Sets the leftmost interface slot at the last slot above threshold
        //printf("Invader interface: %ld %ld\n",Interface_location_invader, aux_interface_invader);
        //Check the potential new interface has a valid value!!! (include it fixs the preious segmentation fault error)
        if (aux_interface_invader==-1){//printf("Invader interface: %ld %ld\n",Interface_location_invader, aux_interface_invader);
        aux_interface_invader=aux_interface_invader+1;}

        if(aux_interface_invader<Interface_location_invader){ //case the interface moves inwards
                Renormalize_inwards(
                        n_spatial_slots,
                        Interface_location_invader,
                        aux_interface_invader,
                        number_of_cells_invader, 
                        ag1s_invader,
                        death_rate_invader, 
                        tau_p_invader, 
                        StochasticNcells_invader, 
                        J_invader, 
                        b_invader, 
                        StochasticAge_invader,
                        survival_rate,
                        diffusion_invader
                        );
        };
        //we do not need to take into account the case of having the interface moving outwards
        //(we already recomputed densities globally)

        Interface_location_invader=aux_interface_invader;

        //Recalcular umbrales de division


        //SAME THING for the host
        aux_interface_host=n_spatial_slots-1;
        
        while((number_of_cells_host[aux_interface_host]>=POPULATION_THRESHOLD)&&(aux_interface_host>0)){
            aux_interface_host--;
        }; //Vector index, nor the natural one
    
        aux_interface_host=aux_interface_host+1; //Sets the leftmost interface slot at the last slot above threshold
        //printf("host interface: %ld %ld\n",Interface_location_host, aux_interface_host);
        //Check the potential new interface has a valid value!!! (include it fixs the preious segmentation fault error)
        if (aux_interface_host==n_spatial_slots){//printf("host interface: %ld %ld\n",Interface_location_host, aux_interface_host);
        aux_interface_host=aux_interface_host-1;}

        if(aux_interface_host>Interface_location_host){ //case the interface moves inwards
                Renormalize_inwards(
                        n_spatial_slots,
                        Interface_location_host,
                        aux_interface_host,
                        number_of_cells_host, 
                        ag1s_host,
                        death_rate_host, 
                        tau_p_host, 
                        StochasticNcells_host, 
                        J_host, 
                        b_host, 
                        StochasticAge_host,
                        survival_rate,
                        diffusion_host
                        );
        };
        //we do not need to take into account the case of having the interface moving outwards
        //(we already recomputed densities globally)

        Interface_location_host=aux_interface_host;
        
        ///Update ag1s at time t+tau  
        for (k=0;k<n_spatial_slots;k++){   
                
                ag1s_invader[k]=Get_division_threshold(oxygen_level[k],
                                                       p6_over_p3_invader,
                                                       aplus_invader);
                ag1s_host[k]=Get_division_threshold(oxygen_level[k],
                                                       p6_over_p3_host,
                                                       aplus_host);
            }
        //Updating ages (only within the stochastic domains)
        for (k=Interface_location_invader;k<n_spatial_slots;k++){   
            for (i=0;i<J_invader[k];i=i+1){
                if(oxygen_level[k]>c_cr(p6_over_p3_invader)){
                    StochasticAge_invader[k][i]=StochasticAge_invader[k][i]+tau;
                };
                //birth condition (X<threshlod):
                if(StochasticAge_invader[k][i]>ag1s_invader[k]){
                    b_invader[k][i]=1./tau_p_invader;
                }
                else{
                    b_invader[k][i]=0.0;
                };
            }; //end for i
        }; //end for k
        
        for (k=0;k<=Interface_location_host;k++){

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
                    b_host[k][i]=0.0;
                };
            }; //end for i
         };
      
        //ES POSIBLE que parte del update que se dice de hacer en 1) convenga posponerlo hasta este punto
        
   }; //end while main loop
   clock_t end = clock();
   time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
   fprintf(stdout,"%lf \n",time_spent);
    //printf("#Tiempo final=%lf \n",t);

    fprintf(MESSAGES,"#end time: %lf \n \n \n",t);
    fprintf(MESSAGES,"#simulation time spent (sec): %lf \n \n \n",time_spent);


    fclose(MESSAGES);

    /////////////////////////////////////////
    //free memory, all matrix and vectors//////  
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
        free(Propensities_host[k]);
        free(b_invader[k]);
        free(Propensities_invader[k]);
    }

    free(J_host);
    free(J_invader);
    free(ag1s_host);
    free(b_host);
    free(ag1s_invader);
    free(b_invader);

    free(Propensities_host);
    free(Propensities_invader);
    free(number_of_cells_host);
    free(number_of_births_host);
    free(density_of_cells_host);
    free(number_of_cells_invader);
    free(number_of_births_invader);
    free(density_of_cells_invader);
    free(total_number_of_cells);
    free(oxygen_level);

    
    
    return 0;
} 