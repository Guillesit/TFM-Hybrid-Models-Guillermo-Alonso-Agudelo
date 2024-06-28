

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




//Use equilibrium formula in delaCruz et al to get the eigenvalue after we integrate out the age structure
//Uses Newton-Rhapson solver...
double Get_Eigenvalue(
                    double *threshold, 
                    long xslot, 
                    double death_rate_inv,
                    double tau_p
                    ){
    
    double aux,eigenvalue=0.0;
    int i;
    

    //Newton's procedure
    
    if(threshold[xslot]>=INFINITE_TIME){eigenvalue=-death_rate_inv;}
    else{
        for(i=0;i<20;i++){
            
            aux=eigenvalue
               -(eigenvalue + death_rate_inv+1/tau_p)
                 *(threshold[xslot]*(eigenvalue + death_rate_inv)
                    +log(1+tau_p*(eigenvalue + death_rate_inv))-log(2)
                 )
                /(1 + threshold[xslot]*(eigenvalue + death_rate_inv+1/tau_p));
            
            eigenvalue=aux;
        };// end for
    };//end if-else
    
    //The function seems well conditioned, 20 iterations should be more than enough to converge
    return(eigenvalue);
    
}; //end Get_Eigenvalue





//Specific routine to get the higher value in a given array
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



