#include "Header.h"


double Refine_time_step(double tau, double Diff_Coef_Pop, double Diff_Coef_O2, double k_oxygen, double peak_number_of_cells){
    
    double delta_t;
    

    delta_t=min(min(tau,min(1.0/(6*Diff_Coef_Pop),min(1.0/(6*Diff_Coef_O2),1.0/(2*k_oxygen*peak_number_of_cells)))),1.0);
    //printf("Candidates: %lf %lf %lf %lf\n",tau,1.0/(6*Diff_Coef_Pop),1.0/(6*Diff_Coef_O2),1.0/(2*k_oxygen*peak_number_of_cells));
    //printf("Refined tstep: %lf\n \n",delta_t);
    
    //Piensa si en vez de cambiar el umbral justo para el Euler del oxigeno --que es 2/(kO2*peakncells)-- por lo que hemos puesto aqui (pensando en tau muy largo y que el peakncells dinamico cambie bastante en este tiempo) no sera mas preciso hacer una estimacion lineal usando tau y el eigenvalue
    
    return (delta_t);
    
};// end Refine_time_step




//Case in which we have already integrated out the age structure
//There is no troubleshooting in case that the deterministic domain is too small, be careful with this
void PopulationSolver(double delta_t, double *N_previous, double *threshold, double death_rate_inv, double tau_p, long x_size, double Diff_Coef){
    
    //Here x_size is now the interface
    
    double *aux;
    long x_slot;
    double  eigenvalue;
    //Note that the eigenvalue depends on space actually
    
    
    //book aux
    if((aux= (double *) malloc(sizeof(double)*(x_size+1)))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    
    //Do explicit finite differences
    for(x_slot=1;x_slot<x_size;x_slot++){
        eigenvalue=Get_Eigenvalue(threshold[x_slot],death_rate_inv,tau_p);
        aux[x_slot]=N_previous[x_slot]*(1+delta_t*eigenvalue)
                +delta_t*Diff_Coef*(N_previous[x_slot-1]+N_previous[x_slot+1]-2*N_previous[x_slot]);
    };
    //left end
    N_previous[0]=aux[1];
    
    //right end minus one (leftmost interface compartment, we DO NOT include proliferation here --apparently this is already taken into account in Gillespie's part)
    aux[x_slot]=N_previous[x_slot]+delta_t*Diff_Coef
            *(N_previous[x_slot-1]-N_previous[x_slot]);
    
    //right end, recall that x_size-1=Interface_location
   // N_previous[x_size]=N_previous[x_size]
    //    +delta_t*Diff_Coef*(N_previous[x_size-1]-N_previous[x_size]);
            
    //swap stuff
    for(x_slot=1;x_slot<=x_size;x_slot++){
        N_previous[x_slot]=aux[x_slot];
    };
    
    free(aux);
    
}; //end PopulationSolver





//This handles both temporal dynamics and diffusion using explicit finite differences to account for oxigen's evolution
void OxygenSolver(double delta_t, double *oxygen_level, double *number_of_cells,  double k_oxygen, double source_oxygen, double Diff_Coef, long xsize){

    
    double diffusion, reaction;
    long x_slot;
   
    
    for(x_slot=1;x_slot<xsize-1;x_slot++){
        reaction=source_oxygen-k_oxygen*oxygen_level[x_slot]*number_of_cells[x_slot];
            
        diffusion= oxygen_level[x_slot-1]+oxygen_level[x_slot+1]-2*oxygen_level[x_slot];
            
        oxygen_level[x_slot]=oxygen_level[x_slot]
                +delta_t*reaction
                +delta_t*Diff_Coef*diffusion;
    };//end loop over space
    //include spatial boundary conditions at the end (crapy zero Neumann)
    //left end:
    oxygen_level[0]=oxygen_level[1];
    //right end:
    oxygen_level[xsize-1]=oxygen_level[xsize-2];
    
}; //End OxygenSolver





void OxygenSolverDamp(
                      double delta_t,
                      double *oxygen_level,
                      double *number_of_cells,
                      double k_oxygen,
                      double source_oxygen,
                      double Diff_coef,
                      long xsize
                      ){
    
    
    double diffusion, reaction;
    long x_slot,i;
    double *temp;
    
    
    //book temp
    if((temp=(double *)malloc(xsize*sizeof(double)))==NULL){
        fprintf(stderr,"Error, memory could not be allocated");
        exit(1);
    };
    
    for(x_slot=1;x_slot<xsize-1;x_slot++){
        
        reaction=source_oxygen-k_oxygen*oxygen_level[x_slot]*number_of_cells[x_slot]
        -oxygen_level[x_slot]*source_oxygen;  //NOTE the very specific damping rate here (11/11/16), after the meeting at CRM on 7/11/16
        
        diffusion= oxygen_level[x_slot-1]+oxygen_level[x_slot+1]-2*oxygen_level[x_slot];
        
        temp[x_slot]=oxygen_level[x_slot]
        +delta_t*reaction
        +delta_t*Diff_coef*diffusion;
    };
    
    //boundary values
    oxygen_level[0]=temp[1];
    oxygen_level[xsize-1]=temp[xsize-2];
    
    //swap stuff
    for(i=1;i<xsize-1;i++){oxygen_level[i]=temp[i];};
    
    free(temp);
    
}; //End OxygenSolverDamp




//This encodes all the mean field handling: time step refinement after we get it from Gillespie's part plus splitting iteration between population and oxigen, and finite differencing
//We advance the oxigen GLOBALLY, both for the mean field and stochastic parts
//Whereas we only update the population in the mean field region
void Global_Mean_Field_Handler(double tau,  double *number_of_cells, double *density_of_cells, double *division_threshold,  double *oxygen_level, double Diff_Coef_Pop, double Diff_Coef_O2, double k_oxygen,double source_oxygen,double delta_x, double tau_p, double death_rate_inv,double p6_over_p3, long Interface_location, long n_xslots){
    
    double delta_t, peak_number_of_cells;
    long x_slot,j=1;
    
    //Determine a refined time step and initialize auxiliar parameters:
    for(x_slot=0;x_slot<Interface_location+1;x_slot++){
        number_of_cells[x_slot]=density_of_cells[x_slot]*delta_x;
    };
    
    peak_number_of_cells=GetMax(number_of_cells,n_xslots); //PODRIA calcular el peak de la densidad y luego multiplicar por delta_x...(05/11/16)
    delta_t=Refine_time_step(tau,Diff_Coef_Pop,Diff_Coef_O2,k_oxygen,peak_number_of_cells);
    
    /***/
    
    //Ahora Tomas sugiere que integremos el oxigeno de forma exacta en vez de con un esquema en diferencias...
    
    //loop over our two solvers advancing succesive time steps until we are close to the full value of tau:
    while(j*delta_t<tau){ //I think we address properly the case delta_t=tau
        //No split in delta_t/2 anymore (07/11/16)
        
        //Population:
        PopulationSolver(delta_t
                         ,density_of_cells,division_threshold,death_rate_inv,tau_p,Interface_location,Diff_Coef_Pop); //eliminacion +1
        
        //Oxygen:
        OxygenSolverDamp(delta_t,oxygen_level,number_of_cells,k_oxygen,source_oxygen,Diff_Coef_O2,n_xslots);
        
        j++;
    };//end while intermediate time steps
    
    //last thrust until we meet tau exactly:
    delta_t=tau-(j-1)*delta_t;
    PopulationSolver(delta_t
                     ,density_of_cells,division_threshold,death_rate_inv,tau_p,Interface_location,Diff_Coef_Pop);//eliminacion +1
    OxygenSolverDamp(delta_t,oxygen_level,number_of_cells,k_oxygen,source_oxygen,Diff_Coef_O2,n_xslots);
    for(x_slot=0;x_slot<=Interface_location;x_slot++){
        number_of_cells[x_slot]=density_of_cells[x_slot]*delta_x;
    };
    
}; //end Global_Mean_Field_Handler





