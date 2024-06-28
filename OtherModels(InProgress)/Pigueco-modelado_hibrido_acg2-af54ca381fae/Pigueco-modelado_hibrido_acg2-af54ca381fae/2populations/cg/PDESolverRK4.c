#include "Header.h"



//Case in which we have already integrated out the age structure
void EvalPopulation(
                    double *kpop, //rhs readout
                    double *N_previous, //old status population
                    double death_rate_inv,
                    double *threshold,
                    double tau_p,
                    double diff_coef,
                    long x_size
                    ){
 
    long x_slot;
    double eigenvalue;
    //Note that the eigenvalue depends on space actually
    
    for(x_slot=1;x_slot<x_size-1;x_slot++){
            eigenvalue=Get_Eigenvalue(threshold,x_slot,death_rate_inv,tau_p);
            kpop[x_slot]=N_previous[x_slot]*eigenvalue
                +diff_coef*(N_previous[x_slot-1]+N_previous[x_slot+1]-2*N_previous[x_slot]);
    };
    
}; //end EvalPopulation




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



double Refine_time_step(
					double tau, 	
					double Diff_Coef_Pop_host, 
					double Diff_Coef_Pop_invader,
					double Diff_Coef_O2, 
					double k_oxygen, 
					double peak_number_of_cells
					){
    
    double delta_t;
    
    //Here we check if tau needs to be refined in order to fulfill:
    //- CFL conditions for both population and oxygen
    //- Euler stability condition (sort of...) for the time solver at oxygen's equation
    //There was not a great advantage in the coarse-grained coded when A-stability criteria
    //where taken into account (roughly same outputs, longer computation times)
    
    delta_t=min(1.0,
    			min(tau,
    				min(1.0/(6*Diff_Coef_Pop_host),
    					min(1.0/(6*Diff_Coef_Pop_invader),
	    					min(1.0/(6*Diff_Coef_O2),1.0/(2*k_oxygen*peak_number_of_cells)
	    						)
    						)
    					)
    				)
    			);
    
    //Posibles alternativas:
    
    //(a) 2.0*max(k_oxygen*peak_number_of_cells,source_oxygen); //in case you want to give it a try
    
    //(b) Piensa si en vez de cambiar el umbral justo para el Euler del oxigeno --que es 2/(kO2*peakncells)-- por lo que hemos puesto aqui (pensando en tau muy largo y que el peakncells dinamico cambie bastante en este tiempo) no sera mas preciso hacer una estimacion lineal usando tau y el eigenvalue
    
    
    return (delta_t);
    
};// end Refine_time_step



void  RK4pop_oxy(
                double delta_t,
				long n_xslots,
				double *number_of_cells_host,
                double *number_of_cells_invader,
                double *total_number_of_cells,
                double death_rate_inv_host,
                double death_rate_inv_invader,
                double *division_threshold_host,
                double *division_threshold_invader,
                double tau_p_host,
                double tau_p_invader,
                double aplus_host,
                double aplus_invader,
                double p6_over_p3_host,
                double p6_over_p3_invader,
                double diff_coef_host,
                double diff_coef_invader,
                double diff_coef_oxygen,
                double *oxygen_level,
                double k_oxygen,
                double source_oxygen
                ){


    
    
    double *k1poph,*k2poph,*k3poph,*k4poph,*auxpoph;
    double *k1popi,*k2popi,*k3popi,*k4popi,*auxpopi;
    double *k1oxy,*k2oxy,*k3oxy,*k4oxy,*auxoxy;
    long i=0;
    
    //book memory
    Book_Vector(&k1poph,n_xslots);
    Book_Vector(&k2poph,n_xslots);
    Book_Vector(&k3poph,n_xslots);
    Book_Vector(&k4poph,n_xslots);
    Book_Vector(&auxpoph,n_xslots);
    
    Book_Vector(&k1popi,n_xslots);
    Book_Vector(&k2popi,n_xslots);
    Book_Vector(&k3popi,n_xslots);
    Book_Vector(&k4popi,n_xslots);
    Book_Vector(&auxpopi,n_xslots);
    
    Book_Vector(&k1oxy,n_xslots);
    Book_Vector(&k2oxy,n_xslots);
    Book_Vector(&k3oxy,n_xslots);
    Book_Vector(&k4oxy,n_xslots);
    Book_Vector(&auxoxy,n_xslots);
    
    //Intermediate steps
    
    //Step 1  //Eval_F(k1,previous_status,constants,t)
    EvalPopulation(k1poph, number_of_cells_host, death_rate_inv_host,division_threshold_host,tau_p_host,diff_coef_host,n_xslots);
    EvalPopulation(k1popi, number_of_cells_invader, death_rate_inv_invader,division_threshold_invader,tau_p_invader,diff_coef_invader,n_xslots);
	EvalOxygen(k1oxy,oxygen_level,n_xslots,total_number_of_cells,k_oxygen,source_oxygen,diff_coef_oxygen);
    
    for(i=1;i<n_xslots-1;i++){
        auxpoph[i]=number_of_cells_host[i]+delta_t*k1poph[i]/2.0;
        auxpopi[i]=number_of_cells_invader[i]+delta_t*k1popi[i]/2.0;
        auxoxy[i]=oxygen_level[i]+delta_t*k1oxy[i]/2.0;
    };
    auxpoph[0]=auxpoph[1];
    auxpoph[n_xslots-1]=auxpoph[n_xslots-2];
    auxpopi[0]=auxpopi[1];
    auxpopi[n_xslots-1]=auxpopi[n_xslots-2];
    auxoxy[0]=auxoxy[1];
    auxoxy[n_xslots-1]=auxoxy[n_xslots-2];
     //actualiza el umbral de division y la poblacion total ya que estas
    Compute_division_threshold(division_threshold_host,auxoxy,p6_over_p3_host,aplus_host,n_xslots-1);
    Compute_division_threshold(division_threshold_invader,auxoxy,p6_over_p3_invader,aplus_invader,n_xslots-1);

    for(i=0;i<n_xslots;i++){ //for the sake of brevity we do not declare an auxiliary array to record this
            total_number_of_cells[i]=auxpoph[i]+auxpopi[i];
    };
     
    
    //Step 2  //Eval_F(k2,aux,constants,t+delta_t/2.0);
    EvalPopulation(k2poph, auxpoph, death_rate_inv_host,division_threshold_host,tau_p_host,diff_coef_host,n_xslots);
    EvalPopulation(k2popi, auxpopi, death_rate_inv_invader,division_threshold_invader,tau_p_invader,diff_coef_invader,n_xslots);
    EvalOxygen(k2oxy,auxoxy,n_xslots,total_number_of_cells,k_oxygen,source_oxygen,diff_coef_oxygen);
    
    for(i=1;i<n_xslots-1;i++){
        auxpoph[i]=number_of_cells_host[i]+delta_t*k2poph[i]/2.0;
        auxpopi[i]=number_of_cells_invader[i]+delta_t*k2popi[i]/2.0;
        auxoxy[i]=oxygen_level[i]+delta_t*k2oxy[i]/2.0;
    };
    auxpoph[0]=auxpoph[1];
    auxpoph[n_xslots-1]=auxpoph[n_xslots-2];
    auxpopi[0]=auxpopi[1];
    auxpopi[n_xslots-1]=auxpopi[n_xslots-2];
    auxoxy[0]=auxoxy[1];
    auxoxy[n_xslots-1]=auxoxy[n_xslots-2];
    
    Compute_division_threshold(division_threshold_host,auxoxy,p6_over_p3_host,aplus_host,n_xslots-1);
    Compute_division_threshold(division_threshold_invader,auxoxy,p6_over_p3_invader,aplus_invader,n_xslots-1);
    for(i=0;i<n_xslots;i++){ 
        total_number_of_cells[i]=auxpoph[i]+auxpopi[i];
    };


    //Step 3   //Eval_F(k3,aux,constants,t+delta_t/2.0);
    EvalPopulation(k3poph, auxpoph, death_rate_inv_host,division_threshold_host,tau_p_host,diff_coef_host,n_xslots);
    EvalPopulation(k3popi, auxpopi, death_rate_inv_invader,division_threshold_invader,tau_p_invader,diff_coef_host,n_xslots);
    EvalOxygen(k3oxy,auxoxy,n_xslots,total_number_of_cells,k_oxygen,source_oxygen,diff_coef_oxygen);
    
    for(i=1;i<n_xslots-1;i++){
      	auxpoph[i]=number_of_cells_host[i]+delta_t*k3poph[i];
      	auxpopi[i]=number_of_cells_invader[i]+delta_t*k3popi[i];
      	auxoxy[i]=oxygen_level[i]+delta_t*k3oxy[i];
    };
    auxpoph[0]=auxpoph[1];
    auxpoph[n_xslots-1]=auxpoph[n_xslots-2];
    auxpopi[0]=auxpopi[1];
    auxpopi[n_xslots-1]=auxpopi[n_xslots-2];
    auxoxy[0]=auxoxy[1];
    auxoxy[n_xslots-1]=auxoxy[n_xslots-2];
    
    Compute_division_threshold(division_threshold_host,auxoxy,p6_over_p3_host,aplus_host,n_xslots-1);
    Compute_division_threshold(division_threshold_invader,auxoxy,p6_over_p3_invader,aplus_invader,n_xslots-1);
    for(i=0;i<n_xslots;i++){ 
            total_number_of_cells[i]=auxpoph[i]+auxpopi[i];
    };


    //Step 4
    //Eval_F(k4,aux,constants,t+delta_t);
    EvalPopulation(k4poph, auxpoph, death_rate_inv_host,division_threshold_host,tau_p_host,diff_coef_host,n_xslots);
    EvalPopulation(k4popi, auxpopi, death_rate_inv_invader,division_threshold_invader,tau_p_invader,diff_coef_invader,n_xslots);
    EvalOxygen(k4oxy,auxoxy,n_xslots,total_number_of_cells,k_oxygen,source_oxygen,diff_coef_oxygen);
    
    
    //RK4 formula
    for(i=1;i<n_xslots-1;i++){
        number_of_cells_host[i]=number_of_cells_host[i]+delta_t*(k1poph[i]+2*k2poph[i]+2*k3poph[i]+k4poph[i])/6.0;
        number_of_cells_invader[i]=number_of_cells_invader[i]+delta_t*(k1popi[i]+2*k2popi[i]+2*k3popi[i]+k4popi[i])/6.0;
        oxygen_level[i]=oxygen_level[i]+delta_t*(k1oxy[i]+2*k2oxy[i]+2*k3oxy[i]+k4oxy[i])/6.0;
    };
    
     
    //Boundary conditions treated after bulk advancement is made
    //Zero Neumann at both sides for both species
    number_of_cells_host[0]=number_of_cells_host[1];
    number_of_cells_host[n_xslots-1]=number_of_cells_host[n_xslots-2];
    number_of_cells_invader[0]=number_of_cells_invader[1];
    number_of_cells_invader[n_xslots-1]=number_of_cells_invader[n_xslots-2];
    oxygen_level[0]=oxygen_level[1];
    oxygen_level[n_xslots-1]=oxygen_level[n_xslots-2];
    
   
    //Update division threshold and total number of cells for the last time
    Compute_division_threshold(division_threshold_host,oxygen_level,p6_over_p3_host,aplus_host,n_xslots-1);
    Compute_division_threshold(division_threshold_invader,oxygen_level,p6_over_p3_invader,aplus_invader,n_xslots-1);
    for(i=0;i<n_xslots;i++){ 
            total_number_of_cells[i]=number_of_cells_host[i]+number_of_cells_invader[i];
    };
    
    
    //free stuff
    free(k1poph);
    free(k2poph);
    free(k3poph);
    free(k4poph);
    free(auxpoph);
    free(k1popi);
    free(k2popi);
    free(k3popi);
    free(k4popi);
    free(auxpopi);
    free(k1oxy);
    free(k2oxy);
    free(k3oxy);
    free(k4oxy);
    free(auxoxy);
    
}; //end RK4pop_oxy




void RK4Handler(	
                double delta_t,
                long n_xslots,
                //double delta_x,
                double *number_of_cells_host,
                double *number_of_cells_invader,
                double death_rate_inv_host,
                double death_rate_inv_invader,
                double *division_threshold_host,
                double *division_threshold_invader,
                double tau_p_host,
                double tau_p_invader,
                double aplus_host,
                double aplus_invader,
                double p6_over_p3_host,
                double p6_over_p3_invader,
                double diff_coef_host,
                double diff_coef_invader,
                double diff_coef_oxygen,
                double *oxygen_level,
                double k_oxygen,
                double source_oxygen
                ){
    
    double DeltaT; //refined time step
    //double Diff_Coef_Pop_host, Diff_Coef_Pop_invader, Diff_Coef_O2=2; //diffusion coefficients
    double peak_number_of_cells_host,peak_number_of_cells_invader,peak_number_of_cells;
    long i,j=1;
    double *total_number_of_cells;
    
    
    /*
    Diff_Coef_O2=diff_coef_oxygen/(delta_x*delta_x);  //It would save some time to do this in Main.c
    Diff_Coef_Pop_host=diff_coef_host/(delta_x*delta_x);
    Diff_Coef_Pop_invader=diff_coef_invader/(delta_x*delta_x);
    //De hecho ninguna subrutina usa delta_x si montamos bien esta parte (solo sirve para reescalar los coeficientes de difusiion)
    */


    //Determine refined time step
    //AQUI LO QUE HAY QUE HACER ES EL PEAK DE LA SUMA (que es lo que realmente entra ahora en la ecuacion del oxigeno)
    Book_Vector(&total_number_of_cells,n_xslots);
    for(i=0;i<n_xslots;i++){
            total_number_of_cells[i]=number_of_cells_host[i]+number_of_cells_invader[i];
        };

    peak_number_of_cells=GetMax(total_number_of_cells,n_xslots);

    
    DeltaT=Refine_time_step( //CHANGE PROTOTYPE
                            delta_t,
                            diff_coef_host,
                            diff_coef_invader,
                            diff_coef_oxygen,
                            k_oxygen,
                            peak_number_of_cells
                            );
    
    
    
    //Create the sub-iteration structure (on each such step we perform a RK4)
    while(j*DeltaT<delta_t){
        
        RK4pop_oxy(
                   DeltaT,
                   n_xslots,
                   //delta_x,
                   number_of_cells_host,
                   number_of_cells_invader,
                   total_number_of_cells,
                   death_rate_inv_host,
                   death_rate_inv_invader,
                   division_threshold_host,
                   division_threshold_invader,
                   tau_p_host,
                   tau_p_invader,
                   aplus_host,
                   aplus_invader,
                   p6_over_p3_host,
                   p6_over_p3_invader,
                   diff_coef_host,
                   diff_coef_invader,
                   diff_coef_oxygen,
                   oxygen_level,
                   k_oxygen,
                   source_oxygen
                   //case_flag,
                   );
        
        j++;
        
    };//end while intermediate time steps
    
    //last thrust until we meet tau exactly:
    DeltaT=delta_t-(j-1)*DeltaT;
    
    //do the stuff again for this special time step
    RK4pop_oxy(
                   DeltaT,
                   n_xslots,
                   number_of_cells_host,
                   number_of_cells_invader,
                   total_number_of_cells,
                   death_rate_inv_host,
                   death_rate_inv_invader,
                   division_threshold_host,
                   division_threshold_invader,
                   tau_p_host,
                   tau_p_invader,
                   aplus_host,
                   aplus_invader,
                   p6_over_p3_host,
                   p6_over_p3_invader,
                   diff_coef_host,
                   diff_coef_invader,
                   diff_coef_oxygen,
                   oxygen_level,
                   k_oxygen,
                   source_oxygen
                   );
	

    free(total_number_of_cells);
    
}; //end RK4Handler
