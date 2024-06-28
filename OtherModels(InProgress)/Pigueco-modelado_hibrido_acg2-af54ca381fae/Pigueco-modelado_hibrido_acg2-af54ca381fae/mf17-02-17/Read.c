

#include "Header.h"





//SUBROUTINES


long Read_Init_Space(FILE *DATA_FILE, FILE *PARAM_FILE, double *pdeath_rate_inv,double *ptau_p,double *psource_oxygen,double *pk_oxygen,double **pinitial_oxygen,double *pp6_over_p3, double **pinitial_population){
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
    if((*pinitial_population=(double *) malloc(nlines*sizeof(double)))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    //Read the info
    for(i=0;i<nlines;i++){
        fscanf(DATA_FILE,"%lf %lf",*pinitial_oxygen+i,*pinitial_population+i);
    };
    
    fclose(DATA_FILE);
    fscanf(PARAM_FILE,"%lf %lf %lf %lf %lf",
           pdeath_rate_inv,ptau_p,psource_oxygen,pk_oxygen,pp6_over_p3);
    fclose(PARAM_FILE);
    
    return(nlines);
    
}; //End Read_Init_Space














