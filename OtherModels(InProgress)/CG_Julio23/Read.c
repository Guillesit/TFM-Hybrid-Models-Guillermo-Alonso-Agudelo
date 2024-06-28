

#include "Header.h"





//SUBROUTINES


long Read_Init_Space(
                    FILE *DATA_FILE, 
                    double **pinitial_population
                    ){

    int c;
    long i, nlines=0;
    
    //Compute number of spatial slots
    while((c=getc(DATA_FILE))!=EOF){if(c=='\n'){nlines++;};};
    nlines++; //We do not jump from last line
    rewind(DATA_FILE);
    
    
    if((*pinitial_population=(double *) malloc(nlines*sizeof(double)))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    //Read the info
    for(i=0;i<nlines;i++){
        fscanf(DATA_FILE,"%lf ",*pinitial_population+i);
    };
    
    fclose(DATA_FILE);
    
    return(nlines);
    
}; //End Read_Init_Space




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
                     double *pDelta_x,
                     double *pDelta_t,
                     long *pn_files
                     ){
    

    struct sim_parameters *p;

    p= (struct sim_parameters *) params;

    fscanf(SIM_DATA,"%lf %lf %lf %lf %lf %ld %lf %lf",
           &(p->k_decay_hat),&(p->k_consumption_hat),&(p->diff_coef_oxygen),
           &(p->source_oxygen_hat),ptstop,pn_files,pDelta_x,pDelta_t
           );
    fclose(SIM_DATA);
    
}; //End Read_params_sim






