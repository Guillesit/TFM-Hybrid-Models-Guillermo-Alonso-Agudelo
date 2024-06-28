

#include "Header.h"





//SUBROUTINES



long Read_Init_Population(
                          FILE *INITIAL_POPULATION,
                          long ***pStochasticNcells,
                          double ***pStochasticAge,
                          long **pnspaces
                          ){
    
    long nlines=0;
    long n_spatial_slots=0;
    int c;
    long i,j;
    
    
    while((c=getc(INITIAL_POPULATION))!=EOF){if(c=='\n'){nlines++;};};
    //nlines++; //We WILL jump from last line
    if(nlines%2!=0){
        fprintf(stderr,"Odd number of lines, aborting execution\n");
        exit(1);
    };
    n_spatial_slots=nlines/2; //plays the role of the former Stochastic_lines
    
    rewind(INITIAL_POPULATION);
    
    //reserva vector para numero de espacios por linea del fichero de entrada
    if((*pnspaces=(long *) calloc(nlines,sizeof(long)))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    nlines=0;
    while((c=getc(INITIAL_POPULATION))!=EOF){
        if(c==' '){
            (*pnspaces)[nlines]++;
        };
        if(c=='\n'){
            //nspaces[nlines]++;//CAREFUL HERE, depends on the input format
            nlines++;
        };
    };
    rewind(INITIAL_POPULATION);
    
    //declara la estructura donde guardas la info
    if((*pStochasticNcells=(long **) calloc(n_spatial_slots,sizeof(long*)))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    if((*pStochasticAge=(double **) calloc(n_spatial_slots,sizeof(double*)))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    for(i=0;i<n_spatial_slots;i++){//*pStochasticNcells+i does not compute
        //nor does (*pStochasticNcells)+i
        //vs (*pStochasticNcells)[i] and *pStochasticNcells[i]
        if(((*pStochasticNcells)[i]=(long *)
            calloc(max((*pnspaces)[2*i],MEMORY_BATCH),sizeof(long)))==NULL){
            fprintf(stderr,"Error, memory could not be assigned \n");
            exit(1);
        };
        if(((*pStochasticAge)[i]=(double *)
            calloc(max((*pnspaces)[1+2*i],MEMORY_BATCH),sizeof(double)))==NULL){
            fprintf(stderr,"Error, memory could not be assigned \n");
            exit(1);
        };
        if((*pnspaces)[2*i]!=(*pnspaces)[1+2*i]){
            fprintf(stderr,"Mismatched number of cells and ages, aborting execution\n");
            exit(1);
        };
    };
    
    //lee tercera vez y guarda la informacion
    //(int-long-double issues??)
    for(i=0;i<n_spatial_slots;i++){
        
        for(j=0;j<(*pnspaces)[2*i];j++){
            fscanf(INITIAL_POPULATION,"%ld",*(*pStochasticNcells+i)+j);
        };
        
        for(j=0;j<(*pnspaces)[1+2*i];j++){
            fscanf(INITIAL_POPULATION,"%lf",*(*pStochasticAge+i)+j);
        };
        
    };
    fclose(INITIAL_POPULATION);
    
    return(n_spatial_slots);
    
};



long Read_Init_Space_Oxygen(
                     FILE *DATA_FILE,
                     double **pinitial_oxygen
                     ){
    //store initial_oxygen in oxygen_level
    int c;
    long i, nlines=0;
    
    //Compute number of spatial slots
    while((c=getc(DATA_FILE))!=EOF){if(c=='\n'){nlines++;};};
    nlines++; //IMPORTANT: We do not jump from last line
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
                            double *pdeath_rate,
                            double *ptau_p,
                            double *pp6_over_p3
                            ){

    fscanf(PARAMETERS,"%lf %lf %lf",
         pdeath_rate,ptau_p,pp6_over_p3);
    fclose(PARAMETERS);


}; //End Read_params_population



void Read_params_sim(
                     FILE *SIM_DATA,
                     double *psource_oxygen,
                     double *pconsumption_oxygen,
                     double *pdelta_t,
                     double *pdelta_x,
                     double *ptstop,
                     int *pn_files,
                     double *psurvival_rate
                     ){
    
    fscanf(SIM_DATA,"%lf %lf %lf %lf %lf %i %lf",
           psource_oxygen,pconsumption_oxygen,pdelta_t,pdelta_x,ptstop,pn_files,psurvival_rate);
    fclose(SIM_DATA);
    
}; //End Read_params_sim






