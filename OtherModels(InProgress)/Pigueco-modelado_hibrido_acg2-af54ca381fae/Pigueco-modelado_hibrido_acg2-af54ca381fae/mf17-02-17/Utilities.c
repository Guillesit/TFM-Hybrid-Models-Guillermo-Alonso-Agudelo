//Utilities for the mean field code


#include "Header.h"





//SUBROUTINES


//compute cummulative number of cells for each spatial slot in a given age-space array
//returns the peak value of the number of cells while we are at it
double Cummulative_Cells(
                         double *number_of_cells,
                         long age_slots,
                         double **Actual_population,
                         long n_xslot
                         ){
    
    long x_slot,i_slot;
    double temp_store;
    
    temp_store=0;
    for(x_slot=0;x_slot<n_xslot;x_slot++){
        number_of_cells[x_slot]=0;
        for(i_slot=0;i_slot<age_slots;i_slot++){
            number_of_cells[x_slot]=number_of_cells[x_slot]+Actual_population[x_slot][i_slot];
        };
        if(number_of_cells[x_slot]>temp_store){
            temp_store=number_of_cells[x_slot];
        };
    };
    
    return(temp_store);
    
};//end Cummulative_Cells




//Books additional age slots as the system ages
void Handle_memory_extension(double **Actual_population,long *page_slots, long space_slots){
    
    long aux,x_slot=0;
    
    aux=*page_slots;
    
    while((x_slot<space_slots)&&(!Actual_population[x_slot][aux-1])){x_slot++;};

    if(x_slot<space_slots){//enlarge all slots
        
        *page_slots=*page_slots+MEMORY_BATCH;
    
        for(x_slot=0;x_slot<space_slots;x_slot++){
            if((Actual_population[x_slot]= (double *) realloc(Actual_population[x_slot],sizeof          (double)*(*page_slots)))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            Actual_population[x_slot][aux]=0.0;
        };
    };//end if additional memory

}; //end Handle_memory_extension



//Data output handling
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
    
}; //end Print_Vector



//From formula (88)??? of the draft we can compute the mean field prediction for the limit oxygen value as time grows big
double Estimate_Limit_O2(
                         double death_rate_inv,
                         double tau_p,
                         double p6_over_p3
                         ){
    
    double rhs;
    
    rhs=(log(2.0)-log(tau_p*(death_rate_inv+1.0/tau_p)))/death_rate_inv;
    
    if(p6_over_p3>r_cr){return(-c_0*log(rhs/Compute_aplus(p6_over_p3)));}
    else{return(1+c_cr(p6_over_p3)*pow(rhs/aminus,beta));};
    
}; //end Estimate_Limit_O2


//Stuff to compute division thresholds depending on the spatial location
void Compute_division_threshold(
                                double *division_threshold,     //array encoding that info vs spatial location
                                double *oxygen_level,   //current oxygen concentration vs space
                                double p6_over_p3,
                                long right_boundary  //rightmost slot entering the computation (see usage below)
){
    
    //right_boundary= Interface_location (to sweep the mean field zone, including the interface (set +1 in the two-compartment case))
    //or
    //right_boundary=n_xslots-1 (to sweep all the system)
    
    double aplus;
    long x_slot;
    
    //using interpolation formula by de la Cruz et al. PLEASE REFER (06/09/2016)
    if(p6_over_p3>r_cr){
        aplus=Compute_aplus(p6_over_p3);
        for(x_slot=0;x_slot<=right_boundary;x_slot++){
            division_threshold[x_slot]=First_division_threshold(aplus,oxygen_level[x_slot]);
        };
    }//end if p6/p3 above threshold
    else{
        for(x_slot=0;x_slot<=right_boundary;x_slot++){
            if(oxygen_level[x_slot]<=c_cr(p6_over_p3)){division_threshold[x_slot]=INFINITE_TIME;}
            else{
                division_threshold[x_slot]=Second_division_threshold(p6_over_p3,oxygen_level[x_slot]);
            };
        }; //end spatial for
    }; //end case p6/p3 below threshold
    
};// end Compute_division_threshold
//NOTES: 1) as it stands, cells do not freeze when they enter quiescence
