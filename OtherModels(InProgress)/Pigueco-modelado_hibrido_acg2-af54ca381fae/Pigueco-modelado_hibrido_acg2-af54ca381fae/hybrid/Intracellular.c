
#include "Header.h"

//RK4 to simulate the modification of Bedessem-Stephanou's model

//NOTES: 1) Roberto has a different implementation for this part, using a god-given RK4 method
//2) Notations for this part can be obtained in WHICH? reference document. TO COMPLETE (06/09/2016)
//3) Apparently Roberto has changed twice the constants and/or specifc equations in the cell cycle model.
// Nothing ensures that this is the last stable version... (06/09/2016)

/************************/
//SUBROUTINES


//Book_vector
void Book_Vector(
    double **pn,    //brand new dynamic array
    long S          //array size
    ){
    
    if((*pn=(double *) malloc(S*sizeof(double)))==NULL){
        printf("Error, memory could not be assigned \n");
        exit(1);
    };
    
} //End Book_Vector





//Specific logistic function
double Logistic(        //returns evaluation of an specific logistic function
    double time, 
    double init,
    double alpha, 
    double m_0
    ){

    return(m_0*init*exp(alpha*time)/(m_0+init*(exp(alpha*time)-1)));
    
} //End Logistic




//Function representing the ODE system
void Eval_F(
    double *outcome,
    double *spatial_arguments,  //current unknown values
    double *constants,      //ODE system constant's
    double time         //current time instant
    ){
    //Indices correspond to generalized coordinates as follows:
    //0=q1, 1=q8, 2=q5, 3=q9, 4=q10
    
    //eq for q1
    outcome[0]=constants[0]-constants[1]
    -constants[2]*spatial_arguments[0];
    
    //eq for q8
    outcome[1]= -constants[3]*spatial_arguments[1]
        -constants[4]*constants[5]*spatial_arguments[1]*spatial_arguments[2]
    +Logistic(time,INIT_MASS,constants[23],constants[24])*constants[6]*spatial_arguments[4]*(1-spatial_arguments[3]/constants[7]);
    
    //eq for q5
    outcome[2]=constants[8]*constants[9]*constants[10]*(constants[11]-spatial_arguments[2])/(constants[12]+constants[11]-spatial_arguments[2])
        -constants[13]*constants[14]*constants[15]*spatial_arguments[1]*spatial_arguments[2]/(constants[16]+spatial_arguments[2]);
 
    
    //eq for q9
    outcome[3]=constants[17]-(constants[18]+constants[19]*constants[20]*spatial_arguments[0])*spatial_arguments[3];
    
    //eq for q10
    outcome[4]=constants[21]-constants[22]*spatial_arguments[4];
  
    
};//End Eval_F






//RK4 abstract method
 void RK4(
    double *previous_status,
    double *updated_status,
    double *constants,  //ODE system constant's
    double t,           //current time instant
    double delta_t,     //desired time step
    int n_eqs           //number of equations
    ){

    double *k1,*k2,*k3,*k4,*aux;
    long i=0;
    
    //book memory
    Book_Vector(&k1,n_eqs);
    Book_Vector(&k2,n_eqs);
    Book_Vector(&k3,n_eqs);
    Book_Vector(&k4,n_eqs);
    Book_Vector(&aux,n_eqs);
    
    //Intermediate steps
    Eval_F(k1,previous_status,constants,t);
 
    for(i=0;i<n_eqs;i++){
        aux[i]=previous_status[i]+delta_t*k1[i]/2.0;
    };
    Eval_F(k2,aux,constants,t+delta_t/2.0);
 
    for(i=0;i<n_eqs;i++){
        aux[i]=previous_status[i]+delta_t*k2[i]/2.0;
    };
    Eval_F(k3,aux,constants,t+delta_t/2.0);
 
    for(i=0;i<n_eqs;i++){
        aux[i]=previous_status[i]+delta_t*k3[i];
    };
    Eval_F(k4,aux,constants,t+delta_t);
 
    //RK4 formula
    for(i=0;i<n_eqs;i++){
        updated_status[i]=previous_status[i]+delta_t*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
    };
     
 
    //free stuff
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(aux);

 }; //End RK4




//Function involved in the computation of the division threshold as a function of oxigen and p6/p3
double Compute_aplus(double p6_over_p3){
    
    double constants[25];
    int n_eqs=5;
    double previous_status[5], updated_status[5];
    double E,S,k7; //scaling parameters (fast and slow variables)
    
    double delta_t=1;
    double t=0.0;
    long i;
    
    //Give initial conditions and constants for the ODE system
    //Taken from Roberto somewhere in late 2014
    n_eqs=5;
    
    S=10.0; //scaling parameter (slow variables)
    E=1.0; //scaling parameter (fast variables)
    k7=1.0;
    
    //System constants
    constants[0]=(0.51*S)/(k7*E*S*S); //kap1
    constants[1]=(0.0085*S*exp(2.5))/(k7*E*S*S); //kap2 for zero oxigen
    constants[2]=1.0/(k7*E*S); //kap3
    
    constants[3]=0.5/(k7*E*S); //kap11
    constants[4]=(1.0/S)/(k7*E); //kap12
    constants[5]=1; //p5 (que valga 1 p. ej)
    constants[6]=0.018/(k7*E*S); //kap10
    constants[7]=1.0; //e2ft
    
    constants[8]=(1.0*S/E)/(k7*S*S); //kap6
    constants[9]=1.0; //p3 (then p6 takes the value of their ratio)
    constants[10]=constants[9]; //pe1
    constants[11]=1.0; //pc
    constants[12]=0.04; //J2
    constants[13]=(14.0/E)/(k7*S); //kap9
    constants[14]=p6_over_p3; //p6
    constants[15]=constants[14]; //pe2
    constants[16]=0.04; //J1
    
    constants[17]=(0.1*S)/(k7*E*S*S); //kap13
    constants[18]=0.1/(k7*E*S); //kap14
    constants[19]=(0.2/S)/(k7*E); //kap15
    constants[20]=1.0; //p1 (que valga 1 p. ej)
    
    constants[21]=0.016*1.0*S/(k7*E*S*S); //kap16
    constants[22]=0.016/(k7*E*S); //kap17
    
    constants[23]=0.005/(k7*E*S); //alpha
    constants[24]=10.0; //m_
    
    //Initial conditions
    //0=q1, 1=q8, 2=q5, 3=q9, 4=q10
    previous_status[0]=0.1;
    previous_status[1]=0.1;
    previous_status[2]=0.9;
    previous_status[3]=1.0;
    previous_status[4]=0.1;
    
    delta_t=0.001; //This should be more than enough to ensure solver's proficiency
    
    //Solve until switch takes place
    while(previous_status[2]>0.1){
        
        RK4(previous_status,updated_status,constants,t,delta_t,n_eqs);
        
        //swap
        for(i=0;i<n_eqs;i++){
            previous_status[i]=updated_status[i];
        };
        
        t=t+delta_t;
    }; //end while
    
    return(t);
    
}; //end Compute_aplus




