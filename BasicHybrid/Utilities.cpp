#include "Header.h"


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


void Print_Vector_Double(
                    FILE *OUTPUT_DATA, 
                    char *output_path, 
                    int iteration, 
                    long vector_size, 
                    double *vector){


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
    
    fclose(OUTPUT_DATA); 

}
void Print_Vector_Long(
                FILE *OUTPUT_DATA, //where to print the stuff
                char *output_path, //where is that file located
                int iteration,  //current iteration of the method
                long vector_size, 
                long *vector
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
        fprintf(OUTPUT_DATA,"%ld\n",vector[slot]);
    };
    fprintf(OUTPUT_DATA,"%ld",vector[slot]);
    
    //close the file
    fclose(OUTPUT_DATA);
    
}; //end Print_Vector

