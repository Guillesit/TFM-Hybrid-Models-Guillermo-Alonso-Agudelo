

#include "Header.h"





//SUBROUTINES

long Read_Init_Space(               //returns the number of spatial slots
    FILE *DATA_FILE,                //Initial population and oxygen values (user-provided)
    FILE *PARAM_FILE,               //Values of some relevant parameters (user-provided) --those below.
    double *pdeath_rate_inv,
    double *ptau_p,
    double *psource_oxygen,         //source term in oxygen's evolution equation
    double *pk_oxygen,              //consumption term in oxygen's evolution equation
    double **pinitial_oxygen,       //initial oxygen array
    double *pp6_over_p3,            //p6/p3 trait (the same for the whole population, not changing during evolution)
    double **pinitial_population  //initial population array
    ){


    //store initial_oxigen in oxigen_level
    int c;
    long i, nlines=0;
    
    //Compute number of spatial slots
    while((c=getc(DATA_FILE))!=EOF){if(c=='\n'){nlines++;};};
    nlines++; 
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
//NOTES: 1) the use of "death_rate_inv" is extremely confusing at some parts of the code,
//we'd better change names and get a consistent notation (06/09/2016) YET TO BE DONE
//2) The structure of the input file is ASSUMED to be such that the last line has no '\n' before 'EOF'












