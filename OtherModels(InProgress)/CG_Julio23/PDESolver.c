#include "Header.h"



void EvalPopulation(struct sim_parameters *params,
                    double *kpop, //rhs readout
                    double *density_previous,
                    double *division_threshold
                    ){
 
    long x_slot;
    double eigenvalue;
    //Note that the eigenvalue depends on space actually
    long n_xslots;
    double death_rate;
    double tau_p;
    double survival_rate;
    double diff_coef_pop;
    struct sim_parameters *p;

    p= (struct sim_parameters *) params;

    n_xslots=(p->n_xslots);
    death_rate=(p->death_rate_hat);
    tau_p=(p->tau_p);
    survival_rate=(p->survival_rate);
    diff_coef_pop=(p->diff_coef_eff);
    //delta_x=(p->Delta_x_hat);

    
    for(x_slot=1;x_slot<n_xslots-1;x_slot++){
        
            eigenvalue=Get_Eigenvalue(division_threshold,
                            x_slot,
                            death_rate,
                            tau_p,
                            survival_rate);
          
            kpop[x_slot]=density_previous[x_slot]*eigenvalue
                +(density_previous[x_slot-1]+density_previous[x_slot+1]-2*density_previous[x_slot])
                    *diff_coef_pop;
    };
   
}; //end EvalPopulation(_targeted)




void EvalOxygen(struct sim_parameters *params,
                    double *koxy, //rhs readout
                    double *oxygen_concentration, //old status oxygen
                    double *density_of_individuals
                    ){
    
    double diffusion, reaction;
    long x_slot;

    struct sim_parameters *p;

    double k_consumption;
    double k_decay;
    double diff_coef;
    double source_oxygen;
    long n_xslots;

    p= (struct sim_parameters *) params;
  
    k_consumption=(p->k_consumption_hat);
    k_decay=(p->k_decay_hat);
    diff_coef=(p->diff_coefs_ratio);
    source_oxygen=(p->source_oxygen_hat);
    n_xslots=(p->n_xslots);
    //delta_x=(p->Delta_x_hat);
    
    for(x_slot=1;x_slot<n_xslots-1;x_slot++){
            
            reaction= source_oxygen 
              -k_consumption*oxygen_concentration[x_slot]*density_of_individuals[x_slot]
                -k_decay*oxygen_concentration[x_slot];
            
            diffusion= oxygen_concentration[x_slot-1]+oxygen_concentration[x_slot+1]-2*oxygen_concentration[x_slot];
            
            koxy[x_slot]=reaction+diff_coef*diffusion;//(delta_x*delta_x);

    }; //end for

}; //End EvalOxygen



void  Eulerpop_oxy( struct sim_parameters *params,
                double delta_t,
                double *density_of_individuals,
                double *division_threshold,
                double *oxygen_concentration
                ){

    
    double *stat1_pop;
    double *stat1_oxy;
    long i=0;

    struct sim_parameters *p;
    long n_xslots;

    p= (struct sim_parameters *) params;
    n_xslots=(p->n_xslots);

    //book memory  
    if((stat1_pop= (double *) malloc(sizeof(double)*n_xslots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned (PDE evolution) \n");
        exit(1);
    };
    if((stat1_oxy= (double *) malloc(sizeof(double)*n_xslots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned (PDE evolution) \n");
        exit(1);
    };
    
    //Step 1  //Eval_F(k1,previous_status,constants,t)

    EvalPopulation(p,
                    stat1_pop, 
                    density_of_individuals,
                    division_threshold
                    );

    EvalOxygen(p,
                stat1_oxy,
                oxygen_concentration,
                density_of_individuals
                );

    
    //Step2: Euler's formula
    for(i=1;i<n_xslots-1;i++){
        density_of_individuals[i]=density_of_individuals[i]+delta_t*stat1_pop[i];
        oxygen_concentration[i]=oxygen_concentration[i]+delta_t*stat1_oxy[i];
    };
    
    //Boundary conditions treated after bulk advancement is made
    //Zero Neumann at both sides for both species
    density_of_individuals[0]=(4*density_of_individuals[1]-density_of_individuals[2])/3.0;
    density_of_individuals[n_xslots-1]=(4*density_of_individuals[n_xslots-2]-density_of_individuals[n_xslots-3])/3.0;
    oxygen_concentration[0]=(4*oxygen_concentration[1]-oxygen_concentration[2])/3.0;
    oxygen_concentration[n_xslots-1]=(4*oxygen_concentration[n_xslots-2]-oxygen_concentration[n_xslots-3])/3.0;
    //In case we want inflow conditions, use something like this: (BEWARE OF THE SPATIAL STEP,to be included)
    //oxygen_level[0]=(4*oxygen_level[1]-oxygen_level[2]+2*incoming_oxygen_flux)/3.0;
    //oxygen_level[n_xslots-1]=(4*oxygen_level[n_xslots-2]-oxygen_level[n_xslots-3]+2*incoming_oxygen_flux)/3.0;
    //The incoming flux needs to be rescaled before that (...)
    
   
    //Update division threshold and total number of cells 
    
    Compute_division_threshold(p,
                                division_threshold,
                                oxygen_concentration,
                                n_xslots-1);
    
    //free stuff
    free(stat1_pop);
    free(stat1_oxy);
    
}; //end Eulerpop_oxy




void GlobalPDE_Handler(	
                struct sim_parameters *params,
                double *density_of_individuals,
                double *division_threshold,
                double *oxygen_concentration
                ){
    
    double DeltaT; //refined time step
    double peak_density_of_cells;
    double stiff_consumption;
    long j=1;
    double TIME_INCREMENT=1.0;

    struct sim_parameters *p;

    double k_consumption;
    long n_xslots;
    double CFL;

    p= (struct sim_parameters *) params;
  
    TIME_INCREMENT=(p->Delta_t_hat); //Time interval to advance in this call (dimensionless)
    k_consumption=(p->k_consumption_hat);
    n_xslots=(p->n_xslots);
    CFL=(p->CFL_number);
    

    /////////////
    //Refine time step according to stability considerations
    peak_density_of_cells=GetMax(density_of_individuals,n_xslots);
    stiff_consumption=0.5/(1.0+k_consumption*peak_density_of_cells); 
    //we hope that during our time window "peak_number_of_cells" will not double, despite how fast it changes through
    
    DeltaT=min(TIME_INCREMENT,min(CFL,stiff_consumption));
    
    //////////////////////////////////
    //Perform temporal iteration
    //Create the sub-iteration structure (on each such step we perform Euler marching)

    while(j*DeltaT<TIME_INCREMENT){ 

        Eulerpop_oxy(p,
                DeltaT,
                density_of_individuals,
                division_threshold,
                oxygen_concentration
                );
        
        j++;
        
    };//end while intermediate time steps
    
    //last thrust until we meet TIME_INCREMENT =1 (internal uds) exactly:
    DeltaT=TIME_INCREMENT-(j-1)*DeltaT;

    
        Eulerpop_oxy(p,
                DeltaT,
                density_of_individuals,
                division_threshold,
                oxygen_concentration
                );

    
}; //end GlobalPDE_Handler
