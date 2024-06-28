#include "Header.h"



//Case in which we have already integrated out the age structure
void Finite_Difference_Solver(
                            double delta_t, 
                            double *N_previous, 
                            double death_rate_inv,
                            double *threshold, 
                            double tau_p, 
                            long x_size, 
                            double delta_x
                            ){
    
    double *aux;
    long x_slot,j;
    double Diff_Coef, CFLnumber, DeltaT, eigenvalue;
    //Note that the eigenvalue depends on space actually
    
    
    //book aux
    if((aux= (double *) malloc(sizeof(double)*x_size))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    
    //rescale diffusion coefficient
    Diff_Coef=DIFFUSION_COEF_CELLS/(delta_x*delta_x);
    CFLnumber=delta_t*Diff_Coef;
    
    //This is NOT TO BE TRUSTED (far too light stability condition, not taking into account reaction terms)
    if(CFLnumber>=0.5){  //refine time step and proceed in several steps until delta_t is reached
        fprintf(stderr,"CFL condition is not verified for this particular setting\n");
        fprintf(stderr,"Refining internal time step\n");
        
        DeltaT= delta_t/(6*CFLnumber); //redefine time step in order to fulfill CFL condition with some margin
        CFLnumber=Diff_Coef*DeltaT;// /(DeltaX*DeltaX);
        j=1;
        
        while(j*DeltaT<delta_t){
            //Do explicit finite differences
            for(x_slot=1;x_slot<x_size-1;x_slot++){
                eigenvalue=Get_Eigenvalue(threshold,x_slot,death_rate_inv,tau_p); //to do this here is time consuming
                aux[x_slot]=N_previous[x_slot]*(1+DeltaT*eigenvalue)
                    +CFLnumber*(N_previous[x_slot-1]+N_previous[x_slot+1]-2*N_previous[x_slot]);
            };
            
            //Zero Neumann boundary conditions, crapy implementation
            N_previous[0]=aux[1];
            N_previous[x_size-1]=aux[x_size-2];
            //swap stuff
            for(x_slot=1;x_slot<x_size-1;x_slot++){
                N_previous[x_slot]=aux[x_slot];
            };
            //move on to next intermediate step
            j++;
        }; //end while intermediate time steps
        
        //Last time step matching with a total advance of delta_t
        DeltaT=delta_t-j*DeltaT;
        CFLnumber=Diff_Coef*DeltaT;
        
        for(x_slot=1;x_slot<x_size-1;x_slot++){
            eigenvalue=Get_Eigenvalue(threshold,x_slot,death_rate_inv,tau_p);
            aux[x_slot]=N_previous[x_slot]*(1+DeltaT*eigenvalue)
                +CFLnumber*(N_previous[x_slot-1]+N_previous[x_slot+1]-2*N_previous[x_slot]);
        };
        
        N_previous[0]=aux[1];
        N_previous[x_size-1]=aux[x_size-2];
        
        for(x_slot=1;x_slot<x_size-1;x_slot++){
            N_previous[x_slot]=aux[x_slot];
        };
        
    }// end if CFL violada
    /*******/
    //Case in which CFL condition is fulfilled
    else{
   
        for(x_slot=1;x_slot<x_size-1;x_slot++){
            eigenvalue=Get_Eigenvalue(threshold,x_slot,death_rate_inv,tau_p);
            aux[x_slot]=N_previous[x_slot]*(1+delta_t*eigenvalue)
                +delta_t*Diff_Coef*(N_previous[x_slot-1]+N_previous[x_slot+1]-2*N_previous[x_slot]);
        };
        
        //Zero Neumann boundary conditions, crapy implementation
        N_previous[0]=aux[1];
        N_previous[x_size-1]=aux[x_size-2];
        //swap stuff
        for(x_slot=1;x_slot<x_size-1;x_slot++){
            N_previous[x_slot]=aux[x_slot];
        };
        
    }; //end else CFL fulfilled
    
    free(aux);
    
}; //end Finite_Difference_Solver





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


/*
//ARGUABLY REFINED VERSION:


double refined_tstep, aux_tstep;
    double diffusion, reaction;
    long x_slot;
    double Diff_Coef;

    double CFLconstraint=1.0;
    double coefs_size=1.0;
    double stability_constant=1.0;
    long i,j,niter=1;
    double *temp;
    
    Diff_Coef=DIFFUSION_COEF_O2/(delta_x*delta_x);
    //CFLnumber=delta_t*Diff_Coef;

    //book temp
    if((temp=(double *)malloc(xsize*sizeof(double)))==NULL){
        fprintf(stderr,"Error, memory could not be allocated");
        exit(1);
    };

    CFLconstraint=0.5*(delta_x*delta_x)/DIFFUSION_COEF_O2;
    coefs_size=2.0*max(k_oxygen*peak_number_of_cells,source_oxygen);
    stability_constant=1.0/max(coefs_size,CFLconstraint);
        
    if(delta_t>=stability_constant){

        niter=lrint(1.5*delta_t/stability_constant)+1;
        refined_tstep=delta_t/niter;

        for(j=0;j<niter;j++){
            //sweep space
            for(x_slot=1;x_slot<xsize-1;x_slot++){
            
            reaction=source_oxygen-k_oxygen*oxygen_level[x_slot]*number_of_cells[x_slot]
                -oxygen_level[x_slot]*source_oxygen; //change discussed on (11/11/16)
            
            diffusion= oxygen_level[x_slot-1]+oxygen_level[x_slot+1]-2*oxygen_level[x_slot];
            
            temp[x_slot]=oxygen_level[x_slot]
                +refined_tstep*reaction
                +refined_tstep*Diff_Coef*diffusion;
            };
        
            //left boundary
            oxygen_level[0]=temp[1];
            //right boundary
            oxygen_level[xsize-1]=temp[xsize-2];
            //swap stuff
            for(i=1;i<xsize-1;i++){oxygen_level[i]=temp[i];};


        };//end for j

    } //end if CFL and/or time stability not fulfilled
    else{ //CFL stuff and stability is OK

        for(x_slot=1;x_slot<xsize-1;x_slot++){
            
            reaction=source_oxygen-k_oxygen*oxygen_level[x_slot]*number_of_cells[x_slot]
                -oxygen_level[x_slot]*source_oxygen; //change discussed on (11/11/16)
            
            diffusion= oxygen_level[x_slot-1]+oxygen_level[x_slot+1]-2*oxygen_level[x_slot];
            
            temp[x_slot]=oxygen_level[x_slot]
                +delta_t*reaction
                +delta_t*Diff_Coef*diffusion;
        };
        
        oxygen_level[0]=temp[1];
        oxygen_level[xsize-1]=temp[xsize-2];
        for(i=1;i<xsize-1;i++){oxygen_level[i]=temp[i];};

    }; //end else

    free(temp);

}; //End OxygenSolverDamp
*/
