

#include "Header.h"





//SUBROUTINES


long Read_Init_Space(
                    FILE *DATA_FILE, 
                    //FILE *PARAM_FILE,
                    //double *pdeath_rate_inv,
                    //double *ptau_p,
                    //double *psource_oxygen,
                    //double *pk_oxygen,
                    //double **pinitial_oxygen,
                    //double *pp6_over_p3,
                    double **pinitial_population//,
                    //double **pbackup_initial_population
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
    //if((*pbackup_initial_population=(double *) malloc(nlines*sizeof(double)))==NULL){
    //    fprintf(stderr,"Error, memory could not be assigned \n");
    //    exit(1);
    //};
    
    //Read the info
    for(i=0;i<nlines;i++){
        fscanf(DATA_FILE,"%lf ",*pinitial_population+i);
        //*(*pbackup_initial_population+i)=*(*pinitial_population+i);
    };
    
    fclose(DATA_FILE);
    //fscanf(PARAM_FILE,"%lf %lf %lf %lf %lf",
      //     pdeath_rate_inv,ptau_p,psource_oxygen,pk_oxygen,pp6_over_p3);
    //fclose(PARAM_FILE);
    
    return(nlines);
    
}; //End Read_Init_Space




long Read_Init_Space_Oxygen(
                     FILE *DATA_FILE,
                     double **pinitial_oxygen
                     ){
    //store initial_oxygen in oxygen_level
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
    for(i=0;i<nlines;i++){
        fscanf(DATA_FILE,"%lf",*pinitial_oxygen+i);
    };
    
    fclose(DATA_FILE);
    
    return(nlines);
    
}; //End Read_Init_Space_Oxygen



void Read_params_population(
                            FILE *PARAMETERS,
                            double *pdeath_rate_inv,
                            double *ptau_p,
                            double *pp6_over_p3,
                            double *pdiff_coef
                            ){

    fscanf(PARAMETERS,"%lf %lf %lf %lf",
         pdeath_rate_inv,ptau_p,pp6_over_p3,pdiff_coef);
    fclose(PARAMETERS);


}; //End Read_params_population



void Read_params_sim(
                     FILE *SIM_DATA,
                     double *psource_oxygen,
                     double *pk_oxygen,
                     double *pdiff_coef,
                     double *pdelta_t,
                     double *pdelta_x,
                     long *psampling_period
                     ){
    
    fscanf(SIM_DATA,"%lf %lf %lf %lf %lf %ld",
           psource_oxygen,pk_oxygen,pdiff_coef,pdelta_t,pdelta_x,psampling_period);
    fclose(SIM_DATA);
    
}; //End Read_params_sim






