
#include "Header.h"


void Read_sim_params(struct all_parameters *params,
                            FILE *PARAMETERS
                            ){

  struct all_parameters *p;

  p= (struct all_parameters *) params;

    fscanf(PARAMETERS,"%lf %lf %lf",
       &(p->delta_x),&(p->delta_t),&(p->tstop));
    fclose(PARAMETERS);


}; 

void Read_eq_params(struct all_parameters *params,
                            FILE *PARAMETERS
                            ){

  struct all_parameters *p;

  p= (struct all_parameters *) params;

    fscanf(PARAMETERS,"%lf %lf %lf",
       &(p->D),&(p->r),&(p->K));
    fclose(PARAMETERS);


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

long Read_Init_Population(FILE *DATA_FILE, double **cells_population){


    
    int c;
    long i, nlines=0;
    
    //Compute number of spatial slots
    while((c=getc(DATA_FILE))!=EOF){if(c=='\n'){nlines++;};};
    nlines++; //We do not jump from last line
    rewind(DATA_FILE);
    
    
    if((*cells_population=(double *) malloc(nlines*sizeof(double)))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    //Read the info
    for(i=0;i<nlines;i++){
        fscanf(DATA_FILE,"%lf ",*cells_population+i);
    };
    
    fclose(DATA_FILE);
    
    return(nlines);

}