

#include "Header.h"





//SUBROUTINES




//compute cummulative number of cells for each spatial slot in a given age-space array
double Cummulative_Cells(   //returns the peak value of the number of cells 
    double *number_of_cells,  //number of cells in each spatial slot
    long age_slots,             //largest number of occupied age slots at current system's state??
    double **Actual_population, //age-space-distribution of cells
    long n_xslot  //rightmost slot entering the computation
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
void Handle_memory_extension( //IT DOES NOT WORK (06/09/2016)
    double **Actual_population,
    long *page_slots, 
    long space_slots
    ){
    
    long x_slot=0;
    
    while((x_slot<space_slots)&&(!Actual_population[x_slot][*page_slots-1])){x_slot++;};

    if(x_slot<space_slots){//enlarge all slots
        
        (*page_slots)++;
    
        for(x_slot=0;x_slot<space_slots;x_slot++){
            if((Actual_population[x_slot]= (double *) realloc(Actual_population[x_slot],sizeof          (double)*(*page_slots)))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            Actual_population[x_slot][(*page_slots)-1]=0.0;
        };
    };//end if additional memory

}; //end Handle_memory_extension


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
    FILE *OUTPUT_DATA,  //where to print the stuff
    char *output_path,  //where is that file located
    int iteration,      //current iteration of the hybrid method
    long vector_size,   
    double *vector
    ){

    char extension[5]=".txt";
    char tag[7]="";
    long slot;
    
    //strcat(output_path,tag[iteration]);
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




//From formula (88) of the draft -numbering has probably changed (27/02/17)- we can compute the mean field prediction for the limit oxygen value as time grows big
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
//NOTES: Do not trust the numbering (88) above, please update (06/09/2016)




//Use equilibrium formula in delaCruz et al to get the eigenvalue after we integrate out the age structure
//Uses Newton-Rhapson solver...
double Get_Eigenvalue(
    double lower_limit,     
    double death_rate_inv,
    double tau_p
    ){
   
    double aux,eigenvalue=0.0;
    int i;
    
    //Newton's procedure to solve WHAT? equation. Please update (06/09/2016)
    //The function seems well conditioned, 20 iterations should be more than enough to converge
    if(lower_limit>=INFINITE_TIME){eigenvalue=-death_rate_inv;}
    else{
        for(i=0;i<20;i++){
            
            aux=eigenvalue
            -(eigenvalue + death_rate_inv+1/tau_p)
            *(lower_limit*(eigenvalue + death_rate_inv)
              +log(1+tau_p*(eigenvalue + death_rate_inv))-log(2)
              )
            /(1 + lower_limit*(eigenvalue + death_rate_inv+1/tau_p));
            
            eigenvalue=aux;
        };// end for
    };//end if-else
    
    return(eigenvalue);
    
}; //end Get_Eigenvalue





//Specific routine to get the higher value in a given array
//Useful for those cases in which we do not loop through age structure
double GetMax( 
    double *number_of_cells,    //array under scrutiny (quite biased...)
    long n_xslot    //rightmost slot
    ){
    
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
    


//The famous ran_2 routine in ranum.h ...
float ran2(long *idum)
{
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    
    if (*idum <= 0) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
};


