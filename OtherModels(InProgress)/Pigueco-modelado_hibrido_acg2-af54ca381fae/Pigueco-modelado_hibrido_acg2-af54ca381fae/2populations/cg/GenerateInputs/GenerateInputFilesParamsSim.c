


#include <stdio.h>


int main(int argc, char **argv)
{
    FILE  *OUTPUT_PARAMETERS;

   
   
    long i;

   
   
    /************************/
    
    /*****************/
    //OPENING FILES
    
    if(argc!=2){
		printf("Error, incorrect argument number\n");
		printf("Use as: \n exec_name  \n  output_parameters_file \n");
		return(1);
	};


    if((OUTPUT_PARAMETERS=fopen(argv[1],"w"))==NULL){
		printf("Error: output_parameters_file could not be opened \n");
		return(1) ;
	};
   
    
    /*********************/
    //Parameters:
    
    //Source term for the oxigen:
    //Consumption rate for the oxygen:
    //Diff coef oxygen (//mm^2/(6seg), taken from Anderson (initially in cm^2/s))
    //Time step
    //Spatial slot size
    //Sampling period
    
    fprintf(OUTPUT_PARAMETERS,"%.10lf %.10lf %.10lf %.10lf %.10lf %ld",
            0.0157000000,0.0001570000,0.006,1.0,1.0,1);
    
    
    fclose(OUTPUT_PARAMETERS);
    
    
	return(0);  
	
}      /*end main*/




