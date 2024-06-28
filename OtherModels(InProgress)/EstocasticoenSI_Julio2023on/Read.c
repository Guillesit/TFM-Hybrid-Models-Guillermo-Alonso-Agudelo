

#include "Header.h"





//SUBROUTINES

//with struct
long Read_Init_Population(
                          FILE *INITIAL_POPULATION,
                          struct age_structure ***pStochasticAges,
                          long **pncells_per_slot,  //lo natural seria long y no double
                          struct sim_parameters *params,
                          double *oxygen_concentration
                          ){
    
    long n_spatial_slots=0;
    int c;
    long i,j;

    struct sim_parameters *p;

    double death_rate;
    //double p6_over_p3=1.0;
    double aminus_hat;
    double c_critic;
    double ag1s;
    double survival_rate;

    double partition_function=1.0;
    double threshold=1.0;
    double sampled_age=0.0;
    double eigenvalue=0.0;

    p= (struct sim_parameters *) params;

    death_rate=(p->death_rate_hat);
    //p6_over_p3=(p->p6_over_p3);
    aminus_hat=(p->aminus_hat);
    c_critic=(p->critical_oxy_hat);
    survival_rate=(p->survival_rate);


    int hour; //to use the random number generator
    double r1; //random numbers
   
    hour = time(NULL);
    srand48(hour);


    /*****************/
    
    //Open file & compute number of spatial slots
    while((c=getc(INITIAL_POPULATION))!=EOF){if(c=='\n'){n_spatial_slots++;};};
    n_spatial_slots++; //We DO NOT jump from last line
    rewind(INITIAL_POPULATION);
    
    //book dynamic vectors
    if((*pncells_per_slot=(long *) calloc(n_spatial_slots,sizeof(long)))==NULL){
      fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
   
    //Read initial number of individuals per slot
    for(i=0;i<n_spatial_slots;i++){
      fscanf(INITIAL_POPULATION,"%ld",*pncells_per_slot+i);
      //
        printf("%ld, %ld\n",*(*pncells_per_slot+i), i); //debug 04/09/23
    };
    
    
    //declara la estructura donde guardas la info   
    if((*pStochasticAges=(struct age_structure**) 
            calloc(n_spatial_slots,sizeof(struct age_structure*)))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    

    //Main spatial loop
    for(i=0;i<n_spatial_slots;i++){

        //Hay que inicializar cada estructura
        if(((*pStochasticAges)[i]=(struct age_structure *) malloc(sizeof(struct age_structure)))==NULL)
        {
            fprintf(stderr,"Error, memory could not be assigned \n");
            exit(1);
        };

        (((*pStochasticAges)[i])->ncells)=(*pncells_per_slot)[i]; 

        (((*pStochasticAges)[i])->size_allocated_memory)=(*pncells_per_slot)[i]+MEMORY_BATCH;

        if((((*pStochasticAges)[i]->age_distribution)=(double *)
                calloc((*pncells_per_slot)[i]+MEMORY_BATCH,sizeof(double)))==NULL)
        {
            fprintf(stderr,"Error, memory could not be assigned \n");
            exit(1);
        };

        /////////////////////////////////////////
        /*Determining age distribution*/

        //Determine division threshold
        //This is for p6/p3 below threshold (specially when it equals one)
        if(oxygen_concentration[i]<=c_critic){ag1s=INFINITE_TIME;}
        else{ag1s=aminus_hat*pow(oxygen_concentration[i]/c_critic-1,-beta);};

        eigenvalue=Get_Eigenvalue(ag1s,
                            death_rate,
                            survival_rate);

        if(eigenvalue>TOL-death_rate){//to avoid dividing by zero or close to when computing the normalization
            partition_function= (1-exp(-ag1s*(eigenvalue+death_rate))/(1+eigenvalue+death_rate))
                                /(eigenvalue+death_rate);
            threshold=(1-exp(-ag1s*(eigenvalue+death_rate))/(eigenvalue+death_rate))/partition_function;

            //This loop to sample the equilibrium age distribution
            for(j=0;j<(*pncells_per_slot)[i];j++){

                r1=drand48();
                if(r1<=threshold){
                   sampled_age=-log(1.0-r1*(eigenvalue+death_rate)*partition_function)
                                    /(eigenvalue+death_rate);
                }
                else{
                   sampled_age=-(log(partition_function*(1-r1)*(1+eigenvalue+death_rate))-ag1s)
                                /(1+eigenvalue+death_rate);
                };
           
                //then write it
                ((*pStochasticAges)[i]->age_distribution)[j]=sampled_age;

            };//end age sampler j-loop


        } //end if nice case
        else{ //close to zero division, different formulae then
            partition_function=1+ag1s;
            threshold=ag1s/partition_function;

            //This loop to sample the equilibrium age distribution
            for(j=0;j<(*pncells_per_slot)[i];j++){

                r1=drand48();
                if(r1<=threshold){ //Double-check all this!!!!!!
                   sampled_age=r1*partition_function;
                }
                else{
                   sampled_age=ag1s-log((1-r1)*partition_function);
                };
           
                //then write it
                ((*pStochasticAges)[i]->age_distribution)[j]=sampled_age;

            };//end age sampler j-loop

        }; //end else strange case

    };//end spatial i-loop
    
    fclose(INITIAL_POPULATION);
    
    return(n_spatial_slots);
    
};




long Read_Init_Space_Oxygen(
                     FILE *DATA_FILE,
                     double **pinitial_oxygen
                     ){
    
    int c;
    long i, nlines=0;
    
    //Compute number of spatial slots
    while((c=getc(DATA_FILE))!=EOF){if(c=='\n'){nlines++;};};
    nlines++; //We do not jump from last line
    rewind(DATA_FILE);
    
    //book dynamic vectors
    if((*pinitial_oxygen=(double *) malloc(nlines*sizeof(double)))==NULL){
      fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
   
    
    //Read the info
    //store initial_oxygen in oxygen_level
    for(i=0;i<nlines;i++){
        fscanf(DATA_FILE,"%lf",*pinitial_oxygen+i);
    };
    
    fclose(DATA_FILE);
    
    return(nlines);
    
}; //End Read_Init_Space_Oxygen



void Read_params_population(struct sim_parameters *params,
                            FILE *PARAMETERS
                            ){

  struct sim_parameters *p;

  p= (struct sim_parameters *) params;

    fscanf(PARAMETERS,"%lf %lf %lf %lf %lf %lf",
       &(p->death_rate_hat),&(p->tau_p),&(p->aplus),&(p->p6_over_p3),
          &(p->diff_coef_pop),&(p->survival_rate));
    fclose(PARAMETERS);


}; //End Read_params_population



void Read_params_sim(struct sim_parameters *params,
                     FILE *SIM_DATA,
                     double *ptstop,
                     long *pn_files,
                     double *pDelta_x,
                    double *pDelta_t
                     ){

  struct sim_parameters *p;

  p= (struct sim_parameters *) params;
    
    fscanf(SIM_DATA,"%lf %lf %lf %lf %lf %ld %lf %lf",
           &(p->k_decay_hat),&(p->k_consumption_hat),&(p->diff_coef_oxygen),
           &(p->source_oxygen_hat),ptstop,pn_files,pDelta_x,pDelta_t);
    fclose(SIM_DATA);
    
}; //End Read_params_sim






