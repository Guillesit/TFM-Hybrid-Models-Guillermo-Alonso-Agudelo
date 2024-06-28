

#include "Header.h"



//SUBROUTINES



//Reverses a string
//This is an example in Kernigan-Ritchie's book
void reverse(char s[]){
    
    int c,i,j;
    
    for(i=0,j=strlen(s)-1;i<j;i++,j--){
        c=s[i];
        s[i]=s[j];
        s[j]=c;
    };
}; //end reverse



//Converts an integer n into a string of characters s
//This is an example in Kernigan-Ritchie's book
void itoa(int n,
          char s[]
          ){
    
    int i, sign;
    
    if((sign=n)<0){n=-n;};
    i=0;
    
    do{
        s[i++]=n%10+'0';
    } while ((n/=10)>0);
    
    if(sign<0){s[i++]='-';};
    s[i]='\0';
    reverse(s);

}; //end itoa



//Data output handling
void Print_Vector(
                FILE *OUTPUT_DATA, //where to print the stuff
                char *output_path, //where is that file located
                int iteration,  //current iteration of the method
                long vector_size, 
                double *vector
                ){

        
    char extension[5]=".txt";
    char tag[7]="";
    long slot;
    

    itoa(iteration,tag);
    strcat(output_path,tag);
    strcat(output_path,extension);
    
    if((OUTPUT_DATA=fopen(output_path,"w"))==NULL){
        fprintf(stderr,"Error: output file could not be opened\n");
        exit(1);
    };
    //print the stuff:
    for(slot=0;slot<vector_size-1;slot++){
        fprintf(OUTPUT_DATA,"%.10lf\n",vector[slot]);
    };
    fprintf(OUTPUT_DATA,"%.10lf",vector[slot]);
    
    //close the file
    fclose(OUTPUT_DATA);
    
}; //end Print_Vector



//Same as before but oncluding abcissae
void Print_VectorLocation(
                FILE *OUTPUT_DATA, //where to print the stuff
                char *output_path, //where is that file located
                char iteration,  //current iteration of the method
                long vector_size, 
                double *vector
                ){

        
    char extension[5]=".txt";
    char tag[7]="";
    long slot;
    

    itoa(iteration,tag);
    strcat(output_path,tag);
    strcat(output_path,extension);
    
    if((OUTPUT_DATA=fopen(output_path,"w"))==NULL){
        fprintf(stderr,"Error: output file could not be opened\n");
        exit(1);
    };
    //print the stuff:
    for(slot=0;slot<vector_size;slot++){
        fprintf(OUTPUT_DATA,"%ld %lf\n",slot,vector[slot]);
    };
    //close the file
    fclose(OUTPUT_DATA);
    
}; //end Print_VectorLocation





//Computing the coarse-grained proliferation rate lambda
//Including a survival rate (see eq. 48 de la Cruez et al 2016)
//We removed tau_p since we perform with adimensional units (see implementation notes)
double Get_Eigenvalue(
                    double *threshold, 
                    long xslot, 
                    double death_rate,
                    double tau_p,//deprecated
                    double survival_rate
                    ){
    
    double aux=0.0;
    double eigenvalue=0.0;
    int i=0;
    

    //Newton's procedure
    
    if(threshold[xslot]>=INFINITE_TIME){eigenvalue=-death_rate;}
    else{
        do{
            eigenvalue=aux;

            aux=eigenvalue
               -(eigenvalue + death_rate+1)
                 *(threshold[xslot]*(eigenvalue + death_rate)
                    +log(1+eigenvalue + death_rate)-log(2*survival_rate)
                 )
                /(1 + threshold[xslot]*(eigenvalue + death_rate+1));

            i++;
        }while((fabs(aux-eigenvalue)>TOL)&&(i<7));
        eigenvalue=aux;
        
    };//end if-else

    return(eigenvalue);
    
}; //end Get_Eigenvalue


//Specific routine to get the highest value in a given array
//Useful for those cases in which we do not loop through age structure
double GetMax(double *number_of_cells,long n_xslot){
    
    double temp_store=0;
    long x_slot;
    
    for(x_slot=0;x_slot<n_xslot;x_slot++){
        if(number_of_cells[x_slot]>temp_store){
            temp_store=number_of_cells[x_slot];
        };
    };
    
    return(temp_store);
    
}; //end GetMax


//Stuff to compute division thresholds depending on the spatial location
void Compute_division_threshold(
            struct sim_parameters *params,
            double *division_threshold,     //array encoding that info vs spatial location
            double *oxygen_concentration,   //current oxygen concentration vs space
            long right_boundary  //rightmost slot entering the computation (aka nslots here)
){
    
    long x_slot;
    struct sim_parameters *p;
    double p6_over_p3;
    double c_crit;
    double aminus_hat;

    p= (struct sim_parameters *) params;
    p6_over_p3=(p->p6_over_p3);
    c_crit=(p->critical_oxy_hat);
    aminus_hat=(p->aminus_hat);

    
    if(p6_over_p3>r_cr){
        for(x_slot=0;x_slot<=right_boundary;x_slot++){
            //division_threshold[x_slot]=aplus*exp(-oxygen_concentration[x_slot]/cecero)/Delta_t;
            fprintf(stderr, "Vuelva usted manana\n");
        };
    }
    else{
        for(x_slot=0;x_slot<=right_boundary;x_slot++){
            if(oxygen_concentration[x_slot]<=c_crit){division_threshold[x_slot]=INFINITE_TIME;}
            else{
                division_threshold[x_slot]=aminus_hat*pow(oxygen_concentration[x_slot]/c_crit-1,-beta);
            };
        }; 
    };
    
};// end Compute_division_threshold

//NOTES: 1) as it stands, cells do not freeze when they enter quiescence
//PERO ESTO YA ESTA CONFIRMADO CON TOMAS

void Compute_equilibria(struct sim_parameters *params,
                        double *pneq,
                        double *pceq
                        ){
        struct sim_parameters *p;
        double k_decay;
        double k_consumption;
        double source_oxygen;
        double survival_rate;
        double death_rate; //death rates 
        double tau_p;
        double p6_over_p3;
        double ccritico;

        p= (struct sim_parameters *) params;
        k_decay=(p->k_decay_hat);
        k_consumption=(p->k_consumption_hat);
        source_oxygen=(p->source_oxygen_hat);
        survival_rate=(p->survival_rate);
        death_rate=(p->death_rate_hat);
        tau_p=(p->tau_p);
        p6_over_p3=(p->p6_over_p3);

        ccritico=(typical_oxy*severe_hypoxia)*(c_cr(p6_over_p3)/c_cr(1));

        *pceq=ccritico*(1+pow(
                        log(2*survival_rate/(1+death_rate*tau_p))/(death_rate*aminus),
                        -1/beta));

        *pneq=(source_oxygen-k_decay*(*pceq))/(k_consumption*(*pceq));


};//end compute equilibria                    


