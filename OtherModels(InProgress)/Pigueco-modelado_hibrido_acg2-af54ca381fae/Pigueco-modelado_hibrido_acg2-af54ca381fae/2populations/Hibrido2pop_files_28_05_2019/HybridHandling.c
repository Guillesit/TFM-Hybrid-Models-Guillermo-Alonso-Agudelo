#include "Header.h"
//#include "ranum.h" //ran2() generator...



double Sample_age_distribution(
                                double ag1s, //scalar, nor vector
                                double eigenvalue,
                                double nu,
                                double taup,
                                double diffusion_coefficient,
                                char *pflag_above_threshold
                                ){

    double partition_function,eps,threshold,r1;
    double sampled_age, taupnu, lambda;

    taupnu=nu*taup;
    eps=nu/diffusion_coefficient; 
    lambda=fabs(eigenvalue);


    partition_function= 1.0/(lambda+eps)
                + exp(-(lambda+eps)*ag1s)*
                (1.0/(lambda+eps*(1.0+1.0/taupnu))-1.0/(lambda+eps));

    threshold=(1.0-exp(-ag1s*(lambda+eps)))/(partition_function*(lambda+eps));

    r1=drand48();
    if(r1<=threshold){
        sampled_age=log(1.0/(1-r1*(lambda+eps)*partition_function))/(lambda+eps);        
    }
    else{
        sampled_age=(1.0/(lambda+eps*(1+1.0/taupnu)))*
                        (ag1s*eps/taupnu
                         +log(1.0/(partition_function*(1-r1)*(lambda+eps*(1+1.0/taupnu))))
                         );
        *pflag_above_threshold=1;
    };

    return(sampled_age*eps/nu);

} //end Sample_age_distribution
           
           

//Insertion and age computation when we add an extra particle to the discrete region coming from roundoffs at the continuous region
void Transfer_particle(
                        long Interface_location,
                        double *ag1s,
                        double death_rate, 
                        double tau_p, 
                        long **StochasticNcells, 
                        long *J, 
                        double **b, 
                        double **StochasticAge, 
                        double survival_rate,
                        double diffusion_coef
                        ){

    double eigenvalue,r;
    double cells_above_division_threshold,cells_below_division_threshold;
    double sampled_age;
    char flag_above_threshold=0;

    
    //now do the stuff

    eigenvalue=Get_Eigenvalue_survival_rate(ag1s,Interface_location,death_rate,tau_p,survival_rate);
    
    r= drand48();

    sampled_age=Sample_age_distribution(
                                ag1s[Interface_location], //scalar, nor vector
                                eigenvalue,
                                death_rate,
                                tau_p,
                                diffusion_coef,
                                &flag_above_threshold
                                );

    if(flag_above_threshold){
        b[Interface_location][J[Interface_location]]=1/tau_p;
    }
    else{
        b[Interface_location][J[Interface_location]]=0; 
    };
    
    
    //insert one cell with sampled_age into Roberto's structure...
    StochasticNcells[Interface_location][J[Interface_location]]=1;  //APPEARS TO BE OK
    StochasticAge[Interface_location][J[Interface_location]]=sampled_age;
    J[Interface_location]=J[Interface_location]+1;
    
    

};  //end Transfer_particle



//Eliminacion particula 
void Eliminate_particle(
                        long location,
                        long **StochasticNcells, 
                        long *J, 
                        double **b, 
                        double **StochasticAge
                        ){

	double r,rCel;
	double auxSuma;
    long suma;
	long j,k;

    
   
    r= drand48();
    suma=0;
    for(k=0;k<J[location];k++){
        suma=suma+StochasticNcells[location][k];
    };


    rCel=suma*r;
    auxSuma=0;
    j=-1;
    while(auxSuma<rCel){
    	j=j+1;
    	auxSuma=auxSuma+StochasticNcells[location][j];
    }
    StochasticNcells[location][j]=StochasticNcells[location][j]-1;
    if(StochasticNcells[location][j]==0)
        {
            for(k=j;k<J[location];k++)
            {
                StochasticNcells[location][k]=StochasticNcells[location][k+1];
                StochasticAge[location][k]=StochasticAge[location][k+1];
                b[location][k]=b[location][k+1];
            }
            J[location]=J[location]-1;
        }
    

};//end Transfer2



//
//Deterministic-to-stochastic flux
void Renormalize_center(
                        double *number_of_cells,  
                        long Interface_location, 
                        long init_index,
                        long end_index,
                        double *ag1s,
                        double death_rate, 
                        double tau_p, 
                        long **StochasticNcells, 
                        long *J, 
                        double **b, 
                        double **StochasticAge,
                        double survival_rate,
                        double diffusion_coef,
                        double cells_at_interface //number of cells just after evolving the stochastic part
                        ){
  
    
    
    double frac_part, r;
    double mean_field_mass=0.0, rescaled_mass=0.0;//
    double temp;
    long i;
    
    temp=number_of_cells[Interface_location];
    frac_part=fmod(temp,1);
    r= drand48();
    
    //compute total mass inside (purely) mean field compartment

    if(Interface_location==end_index){ //interface is at the right

        for(i=init_index;i<end_index;i++){ //compute mass in the purely deterministic region
            mean_field_mass=mean_field_mass+number_of_cells[i];
        };

        //decide how to round the fractional particle number
        if(r<frac_part){//We round up
            number_of_cells[Interface_location]=lrint(temp+1-frac_part); //this is a double with no decimal part
            rescaled_mass=1-(1-frac_part)/mean_field_mass;

            for(i=init_index;i<end_index;i++){ //THIS MAY INTRODUCE SMALL ARTIFACTS, discuss with Tomas
                number_of_cells[i]=rescaled_mass*number_of_cells[i];
            };

            if(fabs(number_of_cells[Interface_location]-cells_at_interface)>0.00001){ 
                    //practically speaking we test that those quantities do not coincide 
                 //(note that both quantities represent integer numbers but are managed as doubles internally)
                Transfer_particle(Interface_location,
                                    ag1s,
                                    death_rate,
                                    tau_p,
                                    StochasticNcells,
                                    J,
                                    b,
                                    StochasticAge,
                                    survival_rate,
                                    diffusion_coef
                                    );
            }; //else we do nothing as we fall back to the state we had BEFORE evolving the mean field part   
            //which seems to be the case 99% of the times 


        }//end if we round up

        else{//We round down
            number_of_cells[Interface_location]=lrint(temp-frac_part); //this is a double with no decimal part
            rescaled_mass=1+frac_part/mean_field_mass;
            for(i=init_index;i<end_index;i++){  //THIS MAY INTRODUCE SMALL ARTIFACTS, discuss with Tomas
                number_of_cells[i]=rescaled_mass*number_of_cells[i];
            };
        
            if(fabs(number_of_cells[Interface_location]-cells_at_interface)>0.00001){ 
                Eliminate_particle(
                        Interface_location,
                        StochasticNcells, 
                        J, 
                        b, 
                        StochasticAge
                        );    
            }; //else we stay put
        
        };//end else we round down


    }//end if interface is at the right


    else{ //interface is at the left
        for(i=init_index+1;i<=end_index;i++){
            mean_field_mass=mean_field_mass+number_of_cells[i];
        };

         //decide how to round the fractional particle number
        if(r<frac_part){//We round up
            number_of_cells[Interface_location]=lrint(temp+1-frac_part);
            rescaled_mass=1-(1-frac_part)/mean_field_mass;

            for(i=init_index+1;i<=end_index;i++){
                number_of_cells[i]=rescaled_mass*number_of_cells[i];
            };
        
            if(fabs(number_of_cells[Interface_location]-cells_at_interface)>0.00001){     
                Transfer_particle(Interface_location,
                                    ag1s,
                                    death_rate,
                                    tau_p,
                                    StochasticNcells,
                                    J,
                                    b,
                                    StochasticAge,
                                    survival_rate,
                                    diffusion_coef
                                    );
            };    

        }//end if we round up

        else{//We round down
            number_of_cells[Interface_location]=lrint(temp-frac_part);
            rescaled_mass=1+frac_part/mean_field_mass;
            for(i=init_index+1;i<=end_index;i++){
                number_of_cells[i]=rescaled_mass*number_of_cells[i];
            };
        
            if(fabs(number_of_cells[Interface_location]-cells_at_interface)>0.00001){ 
                Eliminate_particle(
                        Interface_location,
                        StochasticNcells, 
                        J, 
                        b, 
                        StochasticAge
                        );
            };    
        
        };//end else we round down
        
    
    }; //end else interface is at the left

    
    
}; //End Renormalize_center





//Main subroutine
void Renormalize_inwards(
                        long n_spatial_slots,
                        long Interface_location,
                        long aux_interface,
                        double *number_of_cells,
                        double *division_threshold,
                        double death_rate, 
                        double tau_p, 
                        long **StochasticNcells, 
                        long *J, 
                        double **b, 
                        double **StochasticAge,
                        double survival_rate,
                        double diffusion_coef
                        ){

    long position;

    if(aux_interface<Interface_location){ //moving to the left
        for(position=Interface_location-1;position>=aux_interface;position--){
           
            Renormalize_center(
                        number_of_cells,  
                        position, 
                        0,
                        position,
                        division_threshold,//////////
                        death_rate, 
                        tau_p, 
                        StochasticNcells, 
                        J, 
                        b, 
                        StochasticAge,
                        survival_rate,
                        diffusion_coef,
                        -1  //now there is no previous stochastic evolution to compare with
                        );
        };
    }
    else{ //moving to the right
        for(position=Interface_location+1;position<=aux_interface;position++){

            Renormalize_center(
                        number_of_cells,  
                        position, 
                        position,
                        n_spatial_slots-1,
                        division_threshold,//////////
                        death_rate, 
                        tau_p, 
                        StochasticNcells, 
                        J, 
                        b, 
                        StochasticAge,
                        survival_rate,
                        diffusion_coef,
                        -1
                        );
        };
    };
    
    
    
}; //end Renormalize_inwards



