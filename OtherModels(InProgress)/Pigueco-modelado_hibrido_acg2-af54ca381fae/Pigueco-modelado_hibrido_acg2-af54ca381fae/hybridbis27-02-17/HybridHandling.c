#include "Header.h"
//#include "ranum.h" //ran2() generator...




//Insertion and age computation when we add an extra particle to the discrete region coming from roundoffs at the continuous region
void Transfer_particle(long location,long *pseed,double ag1s,double death_rate_inv, double tau_p, double **age, long *J, double **b, double **ag, double **gilAux, double **ageAux, long *ptam){

    double eigenvalue,r;
    double cells_above_division_threshold,cells_below_division_threshold;
    double sampled_age;

    //deal with possible memory extensions:
    //creemos que no funciona bien
    
    if(J[location]>=*ptam){
    
        *ptam=*ptam+TAM;
        
        if((b[location]= (double *) realloc(b[location],sizeof(double)*(*ptam)))==NULL){
            fprintf(stderr,"Error, memory could not be assigned \n");
            exit(1);
        };
        
        if((age[location]= (double *) realloc(age[location],sizeof(double)*(*ptam)))==NULL){
            fprintf(stderr,"Error, memory could not be assigned \n");
            exit(1);
        };
        
        if((ag[location]= (double *) realloc(ag[location],sizeof(double)*(*ptam)))==NULL){
            fprintf(stderr,"Error, memory could not be assigned \n");
            exit(1);
        };
        
        if((gilAux[location]= (double *) realloc(gilAux[location],sizeof(double)*(*ptam)*4))==NULL){
            fprintf(stderr,"Error, memory could not be assigned \n");
            exit(1);
        };
        
        if((ageAux[location]= (double *) realloc(ageAux[location],sizeof(double)*(*ptam)))==NULL){
            fprintf(stderr,"Error, memory could not be assigned \n");
            exit(1);
        };
    
    
    };

    //now do the stuff

    eigenvalue=Get_Eigenvalue(ag1s,death_rate_inv,tau_p);
    
    //THESE TWO FORMULAS ARE WRONG (check draft) (11/11/16)
    cells_above_division_threshold=exp(-(death_rate_inv+eigenvalue)*ag1s)/(death_rate_inv+eigenvalue+1/tau_p); //B
    cells_below_division_threshold=(1-exp(-(death_rate_inv+eigenvalue)*ag1s))/(death_rate_inv+eigenvalue); //A
    
    r= ran2(pseed);
    
    if(r<cells_below_division_threshold/
        (cells_below_division_threshold+cells_above_division_threshold))
    { //case below division threshold
        r= ran2(pseed);
        sampled_age=-log(r)/(death_rate_inv+eigenvalue);
        b[location][J[location]]=0; //Roberto duda //SEEMS OK
    }
    else{ //case above division threshold
        r= ran2(pseed);
        sampled_age=-log(r)/(death_rate_inv+eigenvalue+1/tau_p);
        b[location][J[location]]=1/tau_p;  //APPEARS TO BE OK
    };
    //insert one cell with sampled_age into Roberto's structure...
    age[location][J[location]]=1;  //APPEARS TO BE OK
    ag[location][J[location]]=sampled_age;
    J[location]=J[location]+1;
    
    
    

};  //end Transfer_particle

//Eliminacion particula (Tranfer_particle2)
void Transfer_particle2(long location,long *pseed,double ag1s,double death_rate_inv, double *number_of_cells, double **age, long *J, double **b, double **ag, double **gilAux, double **ageAux, long *ptam){

	double r,rCel;
	double auxSuma;
	long j,k;
    //deal with possible memory extensions:
    //creemos que no funciona bien
    
    if(J[location]>=*ptam){
    
        *ptam=*ptam+TAM;
        
        if((b[location]= (double *) realloc(b[location],sizeof(double)*(*ptam)))==NULL){
            fprintf(stderr,"Error, memory could not be assigned \n");
            exit(1);
        };
        
        if((age[location]= (double *) realloc(age[location],sizeof(double)*(*ptam)))==NULL){
            fprintf(stderr,"Error, memory could not be assigned \n");
            exit(1);
        };
        
        if((ag[location]= (double *) realloc(ag[location],sizeof(double)*(*ptam)))==NULL){
            fprintf(stderr,"Error, memory could not be assigned \n");
            exit(1);
        };
        
        if((gilAux[location]= (double *) realloc(gilAux[location],sizeof(double)*(*ptam)*4))==NULL){
            fprintf(stderr,"Error, memory could not be assigned \n");
            exit(1);
        };
        
        if((ageAux[location]= (double *) realloc(ageAux[location],sizeof(double)*(*ptam)))==NULL){
            fprintf(stderr,"Error, memory could not be assigned \n");
            exit(1);
        };
    
    
    };

   
    r= ran2(pseed);
    rCel=number_of_cells[location]*r;
    auxSuma=0;
    j=-1;
    while(auxSuma<rCel){
    	j=j+1;
    	auxSuma=auxSuma+age[location][j];
    }
    age[location][j]=age[location][j]-1;
    if(age[location][j]==0)
        {
            for(k=j;k<J[location];k++)
            {
                age[location][k]=age[location][k+1];
                ag[location][k]=ag[location][k+1];
                b[location][k]=b[location][k+1];
            }
            J[location]=J[location]-1;
        }
    
    //age[location][J[location]]=1;  //APPEARS TO BE OK
    //ag[location][J[location]]=sampled_age;
    //J[location]=J[location]+1;
    
    
    

};

//
//Deterministic-to-stochastic flux
void Renormalize_center(double *number_of_cells,  long Interface_location, long *pseed,
                        double ag1s,double death_rate_inv, double tau_p, double **age, long *J, double **b, double **ag,double **gilAux, double **ageAux, long *ptam, double newAux){
    //ag1s above is invoked as division_threshold[Interface_location]
    //pseed is invoked as &seed from main body
    
    double frac_part, r, mean_field_mass=0.0;//
    double temp;
    long i;
    
    temp=number_of_cells[Interface_location];
    frac_part=fmod(temp,1);
    r= ran2(pseed);
    
    //compute total mass inside (purely) mean field compartment
    for(i=0;i<Interface_location;i++){mean_field_mass=mean_field_mass+number_of_cells[i];};
    
   /*if(r<frac_part){//We transfer a whole particle
        number_of_cells[Interface_location]=lrint(temp+1-frac_part);
        mean_field_mass=1-(1-frac_part)/mean_field_mass;
        for(i=0;i<Interface_location;i++){number_of_cells[i]=mean_field_mass*number_of_cells[i];};
        
        //Decide in which age bin does the particle fall and transfer it
        Transfer_particle(Interface_location,pseed,ag1s,death_rate_inv,tau_p,age,J,b,ag,gilAux,ageAux,ptam);
    }
    else{//No particle is transferred
        number_of_cells[Interface_location]=lrint(temp-frac_part);
        mean_field_mass=1+frac_part/mean_field_mass;
        for(i=0;i<Interface_location;i++){number_of_cells[i]=mean_field_mass*number_of_cells[i];};
    };*/
    
    //fprintf(stderr,"%e \n",newAux);
    //fprintf(stderr,"vamos \n");
   if(r<frac_part){//We round up
		number_of_cells[Interface_location]=lrint(temp+1-frac_part);
        mean_field_mass=1-(1-frac_part)/mean_field_mass;
        for(i=0;i<Interface_location;i++){number_of_cells[i]=mean_field_mass*number_of_cells[i];};
        
        if(newAux<number_of_cells[Interface_location]){
        //new cell 
        //fprintf(stderr,"1 \n");
        Transfer_particle(Interface_location,pseed,ag1s,death_rate_inv,tau_p,age,J,b,ag,gilAux,ageAux,ptam);}
        else{//según Tomás quitar una célular para poner otra (con distinta edad) --->NO 
			//fprintf(stderr,"2 \n");
			}
	}
	else{//We round down
		number_of_cells[Interface_location]=lrint(temp-frac_part);
        mean_field_mass=1+frac_part/mean_field_mass;
        for(i=0;i<Interface_location;i++){number_of_cells[i]=mean_field_mass*number_of_cells[i];};
        
        if(newAux<number_of_cells[Interface_location]){
        //quitar una célula 
       	Transfer_particle2(Interface_location,pseed,ag1s,death_rate_inv,number_of_cells,age,J,b,ag,gilAux,ageAux,ptam);
		}
        else{//fprintf(stderr,"4 \n");
        }//no hacer nada?
        
	}
    
    
    
} //End Renormalize_center




//Routine to check where should we position the new interface
//Assumes a single interface, mean field at the left, stochastic at the right
long Realocate_interface(double *number_of_cells, long n_xslots){
    
   long aux=0;
    
    while((number_of_cells[aux]>=POPULATION_THRESHOLD)&&(aux<n_xslots-2)){
        aux++;
    };
    //Troubleshooting for aux=n_xslots-1 is in the main body. Is it?????
    
  return(aux-1); 
    
} //end Realocate_interface




//Converting deterministic compartments into stochastic ones
//Most of the job done by R_l which handles unit displacements
void R_l(long position,double *number_of_cells,long *pseed,double ag1s,double death_rate_inv, double tau_p, double **age, long *J, double **b, double **ag, double **gilAux, double **ageAux, long *ptam, double newAux){
    
    double frac_part, r,temp;
    long i;
    double mean_field_mass=0.0;
    
    temp=number_of_cells[position];// *DeltaX;
    frac_part=fmod(temp,1);
    r= ran2(pseed);
 
    //Compute current mass in mean field compartment
    for(i=0;i<position;i++){mean_field_mass=mean_field_mass+number_of_cells[i];};
    
    //If a whole particle is transferred
   /* if(r<frac_part){
        number_of_cells[position]=lrint(temp+1-frac_part);
        mean_field_mass=1-(1-frac_part)/mean_field_mass;
        for(i=0;i<position;i++){number_of_cells[i]=mean_field_mass*number_of_cells[i];};
        
        //Decide in which age bin does the particle fall and transfer it
        Transfer_particle(position,pseed,ag1s,death_rate_inv,tau_p,age,J,b,ag,gilAux,ageAux,ptam);
    }
    else{ //If no particle is transferred
        number_of_cells[position]=lrint(temp-frac_part);
        mean_field_mass=1+frac_part/mean_field_mass;
        for(i=0;i<position;i++){number_of_cells[i]=mean_field_mass*number_of_cells[i];};
    };*/

    if(r<frac_part){//We round up
        number_of_cells[position]=lrint(temp+1-frac_part);
        mean_field_mass=1-(1-frac_part)/mean_field_mass;
        for(i=0;i<position;i++){number_of_cells[i]=mean_field_mass*number_of_cells[i];};
        
        if(newAux<number_of_cells[position]){
        //new cell 
        Transfer_particle(position,pseed,ag1s,death_rate_inv,tau_p,age,J,b,ag,gilAux,ageAux,ptam);}
        else{//según Tomás quitar una célular para poner otra (con distinta edad)  
            //fprintf(stderr,"2 \n");
            }
    }
    else{//We round down
        number_of_cells[position]=lrint(temp-frac_part);
        mean_field_mass=1+frac_part/mean_field_mass;
        for(i=0;i<position;i++){number_of_cells[i]=mean_field_mass*number_of_cells[i];};
        if(newAux<number_of_cells[position]){
        //quitar una célula 
        Transfer_particle2(position,pseed,ag1s,death_rate_inv,number_of_cells,age,J,b,ag,gilAux,ageAux,ptam);
        }
        else{//fprintf(stderr,"4 \n");
        }//no hacer nada?
    }
    //Recompute parameters for the Gillepie procedure accordingly
 
    
}; //end R_l



//Main subroutine
void Renormalize_left(long Interface_location,long aux_interface,double *number_of_cells,long *pseed, double *division_threshold,double death_rate_inv, double tau_p, double **age, long *J, double **b, double **ag,double **gilAux, double **ageAux, long *ptam,double newAux){

    long position;
    
    for(position=Interface_location-1;position>=aux_interface;position--){
        R_l(position,number_of_cells,pseed,division_threshold[position],death_rate_inv, tau_p, age, J, b, ag,gilAux,ageAux,ptam,newAux); //conversion in the position-th slot
    };
    
} //end Renormalize_left


/*
//Converting stochastic compartments into deterministic ones
 //WE WONT USE SUCH A THING, THE EQUIVALENT OF THIS IS DONE ON THE MAIN CODE
void Renormalize_right(long aux,long Interface_location,double *n_anterior,long *N
            ){
    
    long i;
    
    for(i=Interface_location+1;i<=aux;i++){ 
        n_anterior[i]=N[i];
    };
    
};//end Renormalize_right*/
/****************************************/

