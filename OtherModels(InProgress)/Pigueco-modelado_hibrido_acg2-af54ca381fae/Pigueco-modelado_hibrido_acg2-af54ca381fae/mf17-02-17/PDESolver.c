#include "Header.h"


//Spatial solver for transport equation (age-structured)
//We follow the straight characteristics and integrate the ODE using the same time step and age structure
void FakeTransportSolver (double delta_t,
                          double *N_previous,
                          long number_age_slots,
                          double death_rate_inv,
                          double division_threshold,
                          double tau_p){
    
    double diffusion, reaction,transport;
    long a_slot;
    double *N_updated;
    double sum;

    Book_Vector(&N_updated, number_age_slots);
    sum=0.0;
    for (a_slot=1; a_slot<number_age_slots; a_slot++) {
        //(comentario antiguo) BEWARE! This works only for delta_t=1. If not the movement is not exactly from one slot to the next at its right. In that case: 1) we are to refine the age grid or 2) we are to rescale system constants to keep delta_t=1
        
        //NO, puede trabajar con delta_t distinto de 1 (como de hecho es el caso ya que al hacer el splitting le estamos pasando un delta_t de 1/2), tiene de hecho en cuenta el que los umbrales de division se miden en unidades de "one time unit" y no en las de delta_t. ESTO PUEDE SER CONFUSO AL LEER Y HAY QUE ANOTAR MEJOR EN SUCESIVAS VERSIONES DEL CODIGO
        
        //Hasta ahora no ha dado ningun problema este lio de la doble malla porque la condicion inicial tenia todas las edades cero. Ahora al leer de Roberto hay que encajar sus edades en un mallado de "0.5 time unit" (ya que fuera tomamos delta_t=1 pero luego al splitting le pasamos delta_t=0.5).
        
        N_updated[a_slot]= N_previous[a_slot-1]*exp(-delta_t*death_rate_inv);
       
        
        if(a_slot*delta_t>=division_threshold){
            N_updated[a_slot]=N_updated[a_slot]*exp(-delta_t/tau_p);
            if(N_updated[a_slot]<NEGLIGIBILITY_THRESHOLD){N_updated[a_slot]=0.0;};
            sum=sum+2*(1-exp(-delta_t/tau_p))*N_previous[a_slot-1];
        }
        else if(N_updated[a_slot]<NEGLIGIBILITY_THRESHOLD){N_updated[a_slot]=0.0;};
        
    };
    
    N_updated[0]=sum; //boundary condition
    
    
    //swap step
    for (a_slot=0; a_slot<number_age_slots; a_slot++) {
        N_previous[a_slot]=N_updated[a_slot];
    };
    
    free(N_updated);
    
    
}; //End FakeTransportSolver




//Simulates the diffusion term in the equation for the population
void DiffusionSolver(double delta_t,double **Actual_population,long a_slot, long x_size, double delta_x){
    
    double *aux;
    long x_slot,j;
    double Diff_Coef, CFLnumber, DeltaT;
    
    
    //book aux
    if((aux= (double *) malloc(sizeof(double)*x_size))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    
    //rescale diffusion coefficient
    Diff_Coef=DIFFUSION_COEF_CELLS/(delta_x*delta_x);
    CFLnumber=delta_t*Diff_Coef;
    
    /****/
    //Case in which CFL condition is violated
    if(CFLnumber>=0.5){  //refine time step and proceed in several steps until delta_t is reached
        //fprintf(stderr,"CFL condition is not verified for this particular setting\n");
        //fprintf(stderr,"Refining internal time step\n");
        
        DeltaT= delta_t/(6*CFLnumber); //redefine time step in order to fulfill CFL condition with some margin
        CFLnumber=Diff_Coef*DeltaT;// /(DeltaX*DeltaX);
        j=1;
        
        while(j*DeltaT<delta_t){
            //Do explicit finite differences
            for(x_slot=1;x_slot<x_size-1;x_slot++){
                aux[x_slot]=Actual_population[x_slot][a_slot]+CFLnumber*(Actual_population[x_slot-1][a_slot]+Actual_population[x_slot+1][a_slot]-2*Actual_population[x_slot][a_slot]);
            };
            
            //Zero Neumann boundary conditions, crapy implementation
            Actual_population[0][a_slot]=aux[1];
            Actual_population[x_size-1][a_slot]=aux[x_size-2];
            //swap stuff
            for(x_slot=1;x_slot<x_size-1;x_slot++){
                Actual_population[x_slot][a_slot]=aux[x_slot];
            };
            //move on to next intermediate step
            j++;
        }; //end while intermediate time steps
        
        //Last time step matching with a total advance of delta_t
        DeltaT=delta_t-(j-1)*DeltaT;
        CFLnumber=Diff_Coef*DeltaT;
        
        for(x_slot=1;x_slot<x_size-1;x_slot++){
            aux[x_slot]=Actual_population[x_slot][a_slot]+CFLnumber*(Actual_population[x_slot-1][a_slot]+Actual_population[x_slot+1][a_slot]-2*Actual_population[x_slot][a_slot]);
        };
        
        Actual_population[0][a_slot]=aux[1];
        Actual_population[x_size-1][a_slot]=aux[x_size-2];
    
        for(x_slot=1;x_slot<x_size-1;x_slot++){
            Actual_population[x_slot][a_slot]=aux[x_slot];
        };

    }// end if CFL violada
    /*******/
    //Case in which CFL condition is fulfilled
    else{
    
        //Do explicit finite differences
        for(x_slot=1;x_slot<x_size-1;x_slot++){
            aux[x_slot]=Actual_population[x_slot][a_slot]+delta_t*Diff_Coef*(Actual_population[x_slot-1][a_slot]+Actual_population[x_slot+1][a_slot]-2*Actual_population[x_slot][a_slot]);
        };
    
        //Zero Neumann boundary conditions, crapy implementation
        Actual_population[0][a_slot]=aux[1];
        Actual_population[x_size-1][a_slot]=aux[x_size-2];
        //swap stuff
        for(x_slot=1;x_slot<x_size-1;x_slot++){
            Actual_population[x_slot][a_slot]=aux[x_slot];
        };
    
    }; //end else CFL fulfilled
    
    free(aux);

}; //End DiffusionSolver








void OxygenSolverDamp(
                      double *oxygen_level,
                      long xsize,
                      double *number_of_cells,
                      double delta_t,
                      double k_oxygen,
                      double source_oxygen,
                      double peak_number_of_cells,
                      double delta_x
                      ){
    
    double refined_tstep, aux_tstep;
    double diffusion, reaction;
    long x_slot,i;
    double Diff_Coef, CFLnumber;
    double aux1,aux2;
    double *temp;
    
    
    //book temp
    if((temp=(double *)malloc(xsize*sizeof(double)))==NULL){
        fprintf(stderr,"Error, memory could not be allocated");
        exit(1);
    };
    
    Diff_Coef=DIFFUSION_COEF_O2/(delta_x*delta_x);
    CFLnumber=delta_t*Diff_Coef;
    
    refined_tstep=delta_t;
    aux1=k_oxygen*peak_number_of_cells;
    aux2=2*CFLnumber;
    
    //Checking some sort of stability condition
    //This is NOT TO BE TRUSTED (far too light stability condition)
    if((aux2>=1)||(aux1>=1)){
        refined_tstep=min(1.0/(lrint(aux1)+1),1.0/(lrint(aux2)+1));
    };
    aux_tstep=0;
    
    
    while(aux_tstep<delta_t){
        for(x_slot=1;x_slot<xsize-1;x_slot++){
            
            reaction= source_oxygen-k_oxygen*oxygen_level[x_slot]*number_of_cells[x_slot]
            -oxygen_level[x_slot]*source_oxygen; //change discussed on (11/11/16)
            
            diffusion= oxygen_level[x_slot-1]+oxygen_level[x_slot+1]-2*oxygen_level[x_slot];
            
            temp[x_slot]=oxygen_level[x_slot]
            +refined_tstep*reaction
            +refined_tstep*Diff_Coef*diffusion;
        }; //end for
        
        oxygen_level[0]=temp[1];
        oxygen_level[xsize-1]=temp[xsize-2];
        aux_tstep=aux_tstep+refined_tstep;
        
        //swap stuff
        for(i=1;i<xsize-1;i++){oxygen_level[i]=temp[i];};
        
    }; //end while
    
    //Last time step matching with a total final advance of delta_t
    refined_tstep=delta_t-aux_tstep+refined_tstep; //time defect when we exit the while loop
    
    for(x_slot=1;x_slot<xsize-1;x_slot++){
        
        reaction= source_oxygen-k_oxygen*oxygen_level[x_slot]*number_of_cells[x_slot]
        -oxygen_level[x_slot]*source_oxygen; //change discussed on (11/11/16)
        
        diffusion= oxygen_level[x_slot-1]+oxygen_level[x_slot+1]-2*oxygen_level[x_slot];
        
        temp[x_slot]=oxygen_level[x_slot]
        +refined_tstep*reaction
        +refined_tstep*Diff_Coef*diffusion;
    }; //end for
    
    oxygen_level[0]=temp[1];
    oxygen_level[xsize-1]=temp[xsize-2];
    
    for(i=1;i<xsize-1;i++){oxygen_level[i]=temp[i];};
    
    
    free(temp);
    
}; //End OxygenSolverDamp