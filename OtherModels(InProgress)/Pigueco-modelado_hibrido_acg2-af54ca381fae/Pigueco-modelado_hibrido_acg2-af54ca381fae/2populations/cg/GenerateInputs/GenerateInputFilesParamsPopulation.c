


#include <stdio.h>


int main(int argc, char **argv)
{
    FILE   *OUTPUT_PARAMETERS;
    
  
   
    long i;

   
   
    /************************/
    
    /*****************/
    //OPENING FILES
    
    if(argc!=2){
		printf("Error, incorrect argument number\n");
		printf("Use as: \n exec_name  \n output_parameters_file \n");
		return(1);
	};

    if((OUTPUT_PARAMETERS=fopen(argv[1],"w"))==NULL){
		printf("Error: output_parameters_file could not be opened \n");
		return(1) ;
	};
   
    
    /*********************/
    //Parameters:
    
     //Death rate^-1:
    //tau_p:
    //Ratio of p6 over p3:
    //Diffusion coefficient (mm^2/(6seg), taken from Anderson)
    fprintf(OUTPUT_PARAMETERS,"%.10lf %.10lf %.10lf %.10lf",
            0.0000416667,480.0000000000,1.0000000000,0.0000006);
    
    
    
    fclose(OUTPUT_PARAMETERS);
    
    
	return(0);  
	
}      /*end main*/




