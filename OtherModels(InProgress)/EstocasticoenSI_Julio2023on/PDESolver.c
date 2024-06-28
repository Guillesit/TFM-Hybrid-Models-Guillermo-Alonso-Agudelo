#include "Header.h"






void EvalOxygen(struct sim_parameters *params,
                    double *koxy, //rhs readout
                    double *oxygen_concentration, //old status oxygen
                    long *number_of_cells //old status population 
                    ){
    
    double diffusion, reaction;
    long x_slot;
    //long i;

    struct sim_parameters *p;

    double k_consumption;
    double k_decay;
    double diff_coef;
    double source_oxygen;
    //double delta_x;
    long n_xslots;

    p= (struct sim_parameters *) params;
  
    k_consumption=(p->k_consumption_hat);
    k_decay=(p->k_decay_hat);
    diff_coef=(p->diff_coefs_ratio);
    source_oxygen=(p->source_oxygen_hat);
    n_xslots=(p->n_xslots);
    //delta_x=(p->Delta_x);
    
    for(x_slot=1;x_slot<n_xslots-1;x_slot++){
            
            reaction= source_oxygen 
              -k_consumption*oxygen_concentration[x_slot]*number_of_cells[x_slot]
                -k_decay*oxygen_concentration[x_slot];
            
            diffusion= oxygen_concentration[x_slot-1]+oxygen_concentration[x_slot+1]-2*oxygen_concentration[x_slot];
            
            koxy[x_slot]=reaction+diff_coef*diffusion;
    }; //end for

}; //End EvalOxygen


void  Euler_oxy(struct sim_parameters *params, 
                double delta_t,
				long *number_of_cells,
                double *oxygen_concentration
                ){

    
    double *stat1_oxy;
    long i=0;
    struct sim_parameters *p;
    long n_xslots;

    p= (struct sim_parameters *) params;
    n_xslots=(p->n_xslots);
    
    //book memory  
    if((stat1_oxy= (double *) malloc(sizeof(double)*n_xslots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned (PDE evolution) \n");
        exit(1);
    };
    
    //Step 1  //Eval_F(k1,previous_status,constants,t)
    
    EvalOxygen(p,
                stat1_oxy,
                oxygen_concentration,
                number_of_cells//total_number_of_cells,
                );
    
    
    //Step2: Euler's formula
    for(i=1;i<n_xslots-1;i++){
        oxygen_concentration[i]=oxygen_concentration[i]+delta_t*stat1_oxy[i];
    };
    
    //Boundary conditions treated after bulk advancement is made
    //Zero Neumann at both sides
    oxygen_concentration[0]=(4*oxygen_concentration[1]-oxygen_concentration[2])/3.0;
    oxygen_concentration[n_xslots-1]=(4*oxygen_concentration[n_xslots-2]-oxygen_concentration[n_xslots-3])/3.0;
    //In case we want inflow conditions, use something like this: (BEWARE OF THE SPATIAL STEP,to be included)
    //oxygen_level[0]=(4*oxygen_level[1]-oxygen_level[2]+2*incoming_oxygen_flux)/3.0;
    //oxygen_level[n_xslots-1]=(4*oxygen_level[n_xslots-2]-oxygen_level[n_xslots-3]+2*incoming_oxygen_flux)/3.0;
    //The incoming flux needs to be rescaled before that (...)

    free(stat1_oxy);
    
}; //end Euler_oxy




void GlobalPDE_Handler(struct sim_parameters *params,
                double tau,
                double *oxygen_concentration,
                long *number_of_cells
                ){
    
    double DeltaT; //refined time step
    double peak_number_of_cells;
    double stiff_consumption;
    long j=1;
    double TIME_INCREMENT=1.0;

    struct sim_parameters *p;

    double k_consumption;
    long n_xslots;
    double CFL;

    p= (struct sim_parameters *) params;

    TIME_INCREMENT=(p->Delta_t_hat); //Time interval to advance in this call (dimensionless)
    k_consumption=(p->k_consumption_hat); //recall that this is an effective rate constant (density conversion included)
    n_xslots=(p->n_xslots);
    CFL=(p->CFL_number);
    

    /////////////
    //Refine time step according to stability considerations
    peak_number_of_cells=GetMax_Long(number_of_cells,n_xslots);
    stiff_consumption=0.5/(1.0+k_consumption*peak_number_of_cells); 
    //we hope that during our time window "peak_number_of_cells" will not double, despite how fast it changes through
    
    DeltaT=min(TIME_INCREMENT,
                min(tau,
                    min(CFL,stiff_consumption)));
  
    //////////////////////////////////
    //Perform temporal iteration
    //Create the sub-iteration structure (on each such step we perform Euler marching)
    while(j*DeltaT<tau){
        
        Euler_oxy(p,
                   DeltaT,
                   number_of_cells,
                   oxygen_concentration
                   );
        
        j++;
        
    };//end while intermediate time steps
    
    //last thrust until we meet tau exactly:
    DeltaT=tau-(j-1)*DeltaT;
    
    //do the stuff again for this special time step
    
    Euler_oxy(p,
            DeltaT,
            number_of_cells,
            oxygen_concentration
            );
	
}; //end GlobalPDE_Handler
