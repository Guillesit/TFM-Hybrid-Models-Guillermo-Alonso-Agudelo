


#include <stdio.h>


int main(int argc, char **argv)
{
    FILE  *OUTPUT_DATA, *OUTPUT_PARAMETERS;
    
    
    
    double mu_inverse; //death rate^-1
    double tau_p; //parameter related to stepwise birth rates
    
    double source_oxigen, k_oxigen; //parameters in the equation for the oxigen
    double initial_oxigen;
    double p6_over_p3=1; //ratio of p's
    double init_pop;
   
    long i;

   
   
    /************************/
    
    /*****************/
    //OPENING FILES
    
    if(argc!=3){
		printf("Error, incorrect argument number\n");
		printf("Use as: \n exec_name  \n output_data_file \n output_parameters_file \n");
		return(1);
	};

    if((OUTPUT_DATA=fopen(argv[1],"w"))==NULL){
		printf("Error: output_data_file could not be opened \n");
		return(1) ;
	};	

    if((OUTPUT_PARAMETERS=fopen(argv[2],"w"))==NULL){
		printf("Error: output_parameters_file could not be opened \n");
		return(1) ;
	};
   
    
    /*********************/
    //Parameters:
    
     //Death rate^-1:
    //tau_p:
    //Source term for the oxigen:
    //Consumption rate for the oxigen:
    //Ratio of p6 over p3:
    fprintf(OUTPUT_PARAMETERS,"%.10lf %.10lf %.10lf %.10lf %.10lf",
            0.0000416667,480.0000000000,0.0157000000,0.0001570000,1.0000000000);
    
    
    
    //Data:
    
    //Initial oxigen level:
     //Initial population:
    
    /*************/
    fprintf(OUTPUT_DATA,"%lf %lf",0.0233,4259.0);//4259.0); //In order not to print "\n" after the last datum in the file...
    
    for(i=1;i<10;i++){
        fprintf(OUTPUT_DATA,"\n%lf %lf",0.0233,4259.0);//4259.0);
    };
    
    for(i=10;i<100;i++){
        fprintf(OUTPUT_DATA,"\n%lf %lf",0.0233,0.0);
    };
    
    
    
    
    
    
    /**********/
    fclose(OUTPUT_DATA);
    fclose(OUTPUT_PARAMETERS);
    
    
	return(0);  
	
}      /*end main*/




