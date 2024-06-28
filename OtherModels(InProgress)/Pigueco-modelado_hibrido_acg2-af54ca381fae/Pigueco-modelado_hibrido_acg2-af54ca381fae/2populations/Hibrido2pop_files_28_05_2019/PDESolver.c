#include "Header.h"



 
                    


// //Case in which we have already integrated out the age structure
// void EvalPopulation(
//                     double *kpop, //rhs readout
//                     double *N_previous, //old status population
//                     double death_rate,
//                     double *threshold,
//                     double tau_p,
//                     double diff_coef,
//                     long init_index,
//                     long end_index
//                     ){
 
//     long x_slot;
//     double eigenvalue;
//     //Note that the eigenvalue depends on space actually
    
//     for(x_slot=init_index+1;x_slot<end_index;x_slot++){
//             eigenvalue=Get_Eigenvalue(threshold,x_slot,death_rate,tau_p);
//             //if(x_slot==1){printf("eigenvalue1=%lf\n",eigenvalue);};
//             kpop[x_slot]=N_previous[x_slot]*eigenvalue
//                 +diff_coef*(N_previous[x_slot-1]+N_previous[x_slot+1]-2*N_previous[x_slot]);
//     };
    
// }; //end EvalPopulation


//For those populations undergoing theraphy described by means of a survival fraction
void EvalPopulation_targeted(
                    double *kpop, //rhs readout
                    double *N_previous, //old status population
                    double death_rate,
                    double *threshold,
                    double tau_p,
                    double diff_coef,
                    long init_index,
                    long end_index,
                    double survival_rate
                    ){
 
    long x_slot;
    double eigenvalue;
    //Note that the eigenvalue depends on space actually
    
    for(x_slot=1+init_index;x_slot<end_index;x_slot++){
            eigenvalue=Get_Eigenvalue_survival_rate(threshold,x_slot,death_rate,tau_p,survival_rate);
            //printf("eigenvalue1=%lf\n",eigenvalue);
            kpop[x_slot]=N_previous[x_slot]*eigenvalue
                +diff_coef*(N_previous[x_slot-1]+N_previous[x_slot+1]-2*N_previous[x_slot]);
    };
    
    
}; //end EvalPopulation_targeted




void EvalOxygen(
                    double *koxy, //rhs readout
                    double *oxygen_level, //old status oxygen
                    long xsize, 
                    double *number_of_cells, //old status population (sum of both ones)
                    double k_oxygen, 
                    double source_oxygen,
                    double diff_coef
                    ){
    
    double diffusion, reaction;
    long x_slot,i;
    
    
    
    for(x_slot=1;x_slot<xsize-1;x_slot++){
            
            reaction= source_oxygen-k_oxygen*oxygen_level[x_slot]*number_of_cells[x_slot]
                -oxygen_level[x_slot]*source_oxygen; //change discussed on (11/11/16)
            
            diffusion= oxygen_level[x_slot-1]+oxygen_level[x_slot+1]-2*oxygen_level[x_slot];
            
            koxy[x_slot]=reaction+diff_coef*diffusion;
    }; //end for
    

}; //End EvalOxygen



double Refine_time_step_pop(
                    double tau,     
                    double Diff_Coef
                    ){
    
    double delta_t;
    
    //Here we check if tau needs to be refined in order to fulfill:
    //- CFL conditions for both population and oxygen
    //- Euler stability condition (sort of...) for the time solver at oxygen's equation
    //There was not a great advantage in the coarse-grained coded when A-stability criteria
    //were taken into account (roughly same outputs, longer computation times)
    
    delta_t=min(1.0,
                min(tau,1.0/(6*Diff_Coef)
                    )
                );
    
    //Posibles alternativas:
    
    //(a) 2.0*max(k_oxygen*peak_number_of_cells,source_oxygen); //in case you want to give it a try
    
    //(b) Piensa si en vez de cambiar el umbral justo para el Euler del oxigeno --que es 2/(kO2*peakncells)-- por lo que hemos puesto aqui (pensando en tau muy largo y que el peakncells dinamico cambie bastante en este tiempo) no sera mas preciso hacer una estimacion lineal usando tau y el eigenvalue
    
    
    return (delta_t);
    
};// end Refine_time_step_pop



double Refine_time_step_oxy(
					double tau, 	
					double Diff_Coef_O2, 
					double k_oxygen, 
					double peak_number_of_cells
					){
    
    double delta_t;
    
    //Here we check if tau needs to be refined in order to fulfill:
    //- CFL conditions for both population and oxygen
    //- Euler stability condition (sort of...) for the time solver at oxygen's equation
    //There was not a great advantage in the coarse-grained coded when A-stability criteria
    //were taken into account (roughly same outputs, longer computation times)
    
    delta_t=min(1.0,
    			min(tau,
	    				min(1.0/(6*Diff_Coef_O2),1.0/(2*k_oxygen*peak_number_of_cells)
    					   )
    				)
    			);
    
    //Posibles alternativas:
    
    //(a) 2.0*max(k_oxygen*peak_number_of_cells,source_oxygen); //in case you want to give it a try
    
    //(b) Piensa si en vez de cambiar el umbral justo para el Euler del oxigeno --que es 2/(kO2*peakncells)-- por lo que hemos puesto aqui (pensando en tau muy largo y que el peakncells dinamico cambie bastante en este tiempo) no sera mas preciso hacer una estimacion lineal usando tau y el eigenvalue
    
    
    return (delta_t);
    
};// end Refine_time_step





void  RK4pop(
                double delta_t,
                long n_xslots,
				long init_index,
                long end_index,
				double *density_of_cells,
                double death_rate,
                double *division_threshold,
                double tau_p,
                double diff_coef,
                double survival_rate
                ){


    double *k1pop,*k2pop,*k3pop,*k4pop,*auxpop;
 
    long i=0;
    
    //book memory
    
    Book_Vector(&k1pop,n_xslots);
    Book_Vector(&k2pop,n_xslots);
    Book_Vector(&k3pop,n_xslots);
    Book_Vector(&k4pop,n_xslots);
    Book_Vector(&auxpop,n_xslots);
    
   
    
    
    //Intermediate steps
    
    //Step 1  //Eval_F(k1,previous_status,constants,t)
    EvalPopulation_targeted(k1pop, 
                    density_of_cells, 
                    death_rate,
                    division_threshold,
                    tau_p,
                    diff_coef,
                    init_index,
                    end_index,
                    survival_rate
                    );
    
    for(i=1+init_index;i<end_index;i++){
        auxpop[i]=density_of_cells[i]+delta_t*k1pop[i]/2.0;
    };

    if(init_index==0){ //case mean field region is at the left of the simulation domain
        auxpop[init_index]=auxpop[1+init_index]; 
        k1pop[end_index]=diff_coef*(density_of_cells[end_index-1]-density_of_cells[end_index]);
        auxpop[end_index]=density_of_cells[end_index]+delta_t*k1pop[end_index]/2.0;
    }
    else{ //case mean field region is at the right of the simulation domain
        auxpop[end_index]=auxpop[end_index-1];
        k1pop[init_index]=diff_coef*(density_of_cells[init_index+1]-density_of_cells[init_index]);
        auxpop[init_index]=density_of_cells[init_index]+delta_t*k1pop[init_index]/2.0;
    }; 
    
    //Step 2  //Eval_F(k2,aux,constants,t+delta_t/2.0);
    EvalPopulation_targeted(k2pop, 
                     auxpop, 
                    death_rate,
                    division_threshold,
                    tau_p,
                    diff_coef,
                    init_index,
                    end_index,
                    survival_rate
                    );
   
    
    for(i=1+init_index;i<end_index;i++){
        auxpop[i]=density_of_cells[i]+delta_t*k2pop[i]/2.0;
    };

    if(init_index==0){ //case mean field region is at the left of the simulation domain
        auxpop[init_index]=auxpop[1+init_index]; 
        k2pop[end_index]=diff_coef*(auxpop[end_index-1]-auxpop[end_index]);
        auxpop[end_index]=density_of_cells[end_index]+delta_t*k2pop[end_index]/2.0;
    }
    else{ //case mean field region is at the right of the simulation domain
        auxpop[end_index]=auxpop[end_index-1];
        k2pop[init_index]=diff_coef*(auxpop[init_index+1]-auxpop[init_index]);
        auxpop[init_index]=density_of_cells[init_index]+delta_t*k2pop[init_index]/2.0;
    }; 
   


    //Step 3   //Eval_F(k3,aux,constants,t+delta_t/2.0);
    EvalPopulation_targeted(k3pop, 
                     auxpop, 
                    death_rate,
                    division_threshold,
                    tau_p,
                    diff_coef,
                    init_index,
                    end_index,
                    survival_rate
                    );
    
    for(i=1+init_index;i<end_index;i++){
      	auxpop[i]=density_of_cells[i]+delta_t*k3pop[i];
    };
    
    if(init_index==0){ //case mean field region is at the left of the simulation domain
        auxpop[init_index]=auxpop[1+init_index]; 
        k3pop[end_index]=diff_coef*(auxpop[end_index-1]-auxpop[end_index]);
        auxpop[end_index]=density_of_cells[end_index]+delta_t*k3pop[end_index];
    }
    else{ //case mean field region is at the right of the simulation domain
        auxpop[end_index]=auxpop[end_index-1];
        k3pop[init_index]=diff_coef*(auxpop[init_index+1]-auxpop[init_index]);
        auxpop[init_index]=density_of_cells[init_index]+delta_t*k3pop[init_index]/2.0;
    }; 
    


    //Step 4
    //Eval_F(k4,aux,constants,t+delta_t);
    EvalPopulation_targeted(k4pop, 
                     auxpop, 
                    death_rate,
                    division_threshold,
                    tau_p,
                    diff_coef,
                    init_index,
                    end_index,
                    survival_rate
                    );
    
    
    //RK4 formula
    for(i=1+init_index;i<end_index;i++){
        density_of_cells[i]=density_of_cells[i]+delta_t*(k1pop[i]+2*k2pop[i]+2*k3pop[i]+k4pop[i])/6.0;
    };
    
     
    //Boundary conditions treated after bulk advancement is made
    if(init_index==0){ //case mean field region is at the left of the simulation domain
        density_of_cells[init_index]=density_of_cells[1+init_index]; 
        k4pop[end_index]=diff_coef*(auxpop[end_index-1]-auxpop[end_index]);
        density_of_cells[end_index]=density_of_cells[end_index]
                    +delta_t*(k1pop[end_index]+2*k2pop[end_index]+2*k3pop[end_index]+k4pop[end_index])/6.0;
    }
    else{ //case mean field region is at the right of the simulation domain
        density_of_cells[end_index]=density_of_cells[end_index-1];
        k4pop[init_index]=diff_coef*(auxpop[init_index+1]-auxpop[init_index]);
        density_of_cells[init_index]=density_of_cells[init_index]
                    +delta_t*(k1pop[init_index]+2*k2pop[init_index]+2*k3pop[init_index]+k4pop[init_index])/6.0;
    }; 

    
    //
    
    

    
    //free stuff
    free(k1pop);
    free(k2pop);
    free(k3pop);
    free(k4pop);
    free(auxpop);
   
    
    
}; //end RK4pop



void  RK4oxy(
                double delta_t,
                long n_xslots,
                double *total_number_of_cells,
                double diff_coef_oxygen,
                double *oxygen_level,
                double k_oxygen,
                double source_oxygen
                ){


    
    double *k1oxy,*k2oxy,*k3oxy,*k4oxy,*auxoxy;
    long i=0;
    
    //book memory
    
    Book_Vector(&k1oxy,n_xslots);
    Book_Vector(&k2oxy,n_xslots);
    Book_Vector(&k3oxy,n_xslots);
    Book_Vector(&k4oxy,n_xslots);
    Book_Vector(&auxoxy,n_xslots);
    
    //Intermediate steps
    
    //Step 1  //Eval_F(k1,previous_status,constants,t)
    
    EvalOxygen(k1oxy,oxygen_level,n_xslots,total_number_of_cells,k_oxygen,source_oxygen,diff_coef_oxygen);
    
    for(i=1;i<n_xslots-1;i++){
        auxoxy[i]=oxygen_level[i]+delta_t*k1oxy[i]/2.0;
    };
    auxoxy[0]=auxoxy[1];
    auxoxy[n_xslots-1]=auxoxy[n_xslots-2];

    
    //Step 2  //Eval_F(k2,aux,constants,t+delta_t/2.0);
   
    EvalOxygen(k2oxy,auxoxy,n_xslots,total_number_of_cells,k_oxygen,source_oxygen,diff_coef_oxygen);
    
    for(i=1;i<n_xslots-1;i++){
        auxoxy[i]=oxygen_level[i]+delta_t*k2oxy[i]/2.0;
    };
    auxoxy[0]=auxoxy[1];
    auxoxy[n_xslots-1]=auxoxy[n_xslots-2];
    
   


    //Step 3   //Eval_F(k3,aux,constants,t+delta_t/2.0);
 
    EvalOxygen(k3oxy,auxoxy,n_xslots,total_number_of_cells,k_oxygen,source_oxygen,diff_coef_oxygen);
    
    for(i=1;i<n_xslots-1;i++){
        auxoxy[i]=oxygen_level[i]+delta_t*k3oxy[i];
    };
    
    auxoxy[0]=auxoxy[1];
    auxoxy[n_xslots-1]=auxoxy[n_xslots-2];
    
    


    //Step 4
    //Eval_F(k4,aux,constants,t+delta_t);
    
    EvalOxygen(k4oxy,auxoxy,n_xslots,total_number_of_cells,k_oxygen,source_oxygen,diff_coef_oxygen);
    
    
    //RK4 formula
    for(i=1;i<n_xslots-1;i++){
        
        oxygen_level[i]=oxygen_level[i]+delta_t*(k1oxy[i]+2*k2oxy[i]+2*k3oxy[i]+k4oxy[i])/6.0;
    };
    
     
    //Boundary conditions treated after bulk advancement is made
    //Zero Neumann at both sides for both species
    
    oxygen_level[0]=oxygen_level[1];
    oxygen_level[n_xslots-1]=oxygen_level[n_xslots-2];
    
   
    //Update division threshold and total number of cells for the last time
   
    
    
    //free stuff
    
    free(k1oxy);
    free(k2oxy);
    free(k3oxy);
    free(k4oxy);
    free(auxoxy);
    
}; //end RK4oxy



void RK4HandlerPop(	
                double delta_t,
                long n_xslots,
                long init_index,
                long end_index,
                double *density_of_cells,
                double death_rate,
                double *division_threshold,
                double tau_p,
                double diff_coef,
                double survival_rate
                ){
    
    double DeltaT; //refined time step
    long i,j=1;
    
    
    //Determine refined time step
    DeltaT=Refine_time_step_pop( 
                            delta_t,
                            diff_coef
                            );
    
    
    //Create the sub-iteration structure (on each such step we perform a RK4)
    while(j*DeltaT<delta_t){
        
        RK4pop(
                   DeltaT,
                   n_xslots,
                   init_index,
                   end_index,
                   density_of_cells,
                   death_rate,
                   division_threshold,
                   tau_p,
                   diff_coef,
                   survival_rate
                   );
        
        j++;
        
    };//end while intermediate time steps
    
    //last thrust until we meet tau exactly:
    DeltaT=delta_t-(j-1)*DeltaT;
    
    //do the stuff again for this special time step
    RK4pop(
                   DeltaT,
                   n_xslots,
                   init_index,
                   end_index,
                   density_of_cells,
                   death_rate,
                   division_threshold,
                   tau_p,
                   diff_coef,
                   survival_rate
                   );
	
    
}; //end RK4HandlerPop


void RK4HandlerOxy(    
                double delta_t,
                long n_xslots,
                double *number_of_cells_host, //DE HECHO TENEMOS EL TOTAL NUMBER OF CELLS EN LA RUTINA PRINCIPAL
                double *number_of_cells_invader,
                double diff_coef_oxygen,
                double *oxygen_level,
                double k_oxygen,
                double source_oxygen
                ){
    
    double DeltaT; //refined time step
    double peak_number_of_cells_host,peak_number_of_cells_invader,peak_number_of_cells;
    long i,j=1;
    double *total_number_of_cells;
    
    


    //Determine refined time step
    //AQUI LO QUE HAY QUE HACER ES EL PEAK DE LA SUMA (que es lo que realmente entra ahora en la ecuacion del oxigeno)
    Book_Vector(&total_number_of_cells,n_xslots);
    for(i=0;i<n_xslots;i++){
            total_number_of_cells[i]=number_of_cells_host[i]+number_of_cells_invader[i];
        };

    peak_number_of_cells=GetMax(total_number_of_cells,n_xslots);

    
    DeltaT=Refine_time_step_oxy( 
                            delta_t,
                            diff_coef_oxygen,
                            k_oxygen,
                            peak_number_of_cells
                            );
    
    
    
    //Create the sub-iteration structure (on each such step we perform a RK4)
    while(j*DeltaT<delta_t){
        
        RK4oxy(
                   DeltaT,
                   n_xslots,
                   total_number_of_cells,
                   diff_coef_oxygen,
                   oxygen_level,
                   k_oxygen,
                   source_oxygen
                   );
        
        j++;
        
    };//end while intermediate time steps
    
    //last thrust until we meet tau exactly:
    DeltaT=delta_t-(j-1)*DeltaT;
    
    //do the stuff again for this special time step
    RK4oxy(
                   DeltaT,
                   n_xslots,
                   total_number_of_cells,
                   diff_coef_oxygen,
                   oxygen_level,
                   k_oxygen,
                   source_oxygen
                   );
    

    free(total_number_of_cells);
    
}; //end RK4HandlerOxy

