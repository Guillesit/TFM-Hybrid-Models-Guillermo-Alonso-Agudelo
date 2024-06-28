


#include <stdio.h>


int main(int argc, char **argv)
{
    FILE  *OUTPUT_DATA, *OUTPUT_PARAMETERS;
    
  
   
    long i;

   
   
    /************************/
    
    /*****************/
    //OPENING FILES
    
    if(argc!=2){
		printf("Error, incorrect argument number\n");
		printf("Use as: \n exec_name  \n output_data_file \n");
		return(1);
	};

    if((OUTPUT_DATA=fopen(argv[1],"w"))==NULL){
		printf("Error: output_data_file could not be opened \n");
		return(1) ;
	};	

    
    
    /*************/
    
    
    fprintf(OUTPUT_DATA,"%lf",0.0233); //In order not to print "\n" after the last datum in the file...
    for(i=1;i<25;i++){
        fprintf(OUTPUT_DATA,"\n%lf",0.0233);
    };
    
    
    
    
    
    
    /**********/
    fclose(OUTPUT_DATA);
    
    
    
	return(0);  
	
}      /*end main*/




