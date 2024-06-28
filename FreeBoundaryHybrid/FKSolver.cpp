
#include "Header.h"


void eval_c(struct all_parameters *params, double *concentration, double *rd,long interface) //interface isnt used for the moment
                    {
    double reaction, diffusion;
    struct all_parameters *p;
    p= (struct all_parameters *) params;

    double K; 
    double r; 
    double D;
    long n_xslots; 
    double delta_x; 

    K=(p->K);
    r=(p->r);
    D=(p->D);
    n_xslots=(p->n_xslots);
    delta_x=(p->delta_x); 

    long x_slot; 

    long limit;

    n_xslots=(p->n_xslots);
    if (interface<=n_xslots-1){
        limit=interface;
        rd[interface]=D*(( concentration[interface-1] - concentration[interface])/(delta_x*delta_x));

    }else{
        limit=n_xslots-1;
        rd[n_xslots-1]=r * concentration[n_xslots-1] - (r/K)*concentration[n_xslots-1]*concentration[n_xslots-1]+D*(( concentration[n_xslots-2] - concentration[n_xslots-1])/(delta_x*delta_x));
    }

    rd[0]=r * concentration[0] - (r/K)*concentration[0]*concentration[0]+D*(( concentration[1] - concentration[0])/(delta_x*delta_x));
    
    for(x_slot=1;x_slot<limit;x_slot++){
        reaction = r * concentration[x_slot] - (r/K)*concentration[x_slot]*concentration[x_slot]; 
        // be careful when rescaling we need to rescale also delta_x and the dif coef 
        diffusion = (concentration[x_slot-1] + concentration[x_slot+1] - 2*concentration[x_slot])/(delta_x*delta_x); 
        rd[x_slot] = reaction + D*diffusion; 
        //printf("%lf\n",rd[x_slot]);
    }
    
    //rd[n_xslots-1]=r * concentration[n_xslots-1] - (r/K)*concentration[n_xslots-1]*concentration[n_xslots-1]+D*(( concentration[n_xslots-2] - concentration[n_xslots-1])/(delta_x*delta_x));
    //rd[interface]=D*(( concentration[interface-1] - concentration[interface])/(delta_x*delta_x));
    } //TO DO: CONSIDER WHEN ALL DETER
    
    
   

void euler(struct all_parameters *params, 
                double *concentration,
                long interface,
                double tau
                    ){
    long i = 0; 
    double *rd; 
    long n_xslots;  
    
    struct all_parameters *p;
    p= (struct all_parameters *) params;

    long limit;
    n_xslots=(p->n_xslots);
    if (interface<=n_xslots-1){
        limit=interface;
    }else{
        limit=n_xslots-1;
    }
     
      // Allocate memory for rd
    rd = (double *)malloc(n_xslots * sizeof(double));
    if (rd == NULL) {
        fprintf(stderr, "Error: memory allocation failed\n");
        exit(1);
    }
    eval_c(p, concentration,rd,interface) ; 
    //printf("%lf\n",tau); //Errorazo si tau grande
    // Rewrite concentration with the concentration of the new time step 
    for(i=0;i<=limit;i++){ //TO DO: CONSIDER WHEN ALL DETER
        // Advance in space except for the boundaries 
        // calculate the function using the previous step
        //printf("%lf\n",concentration[i]);
        
        concentration[i] = concentration[i] + tau * rd[i];    
        //printf("%lf\n",concentration[i]);
        }
        // Neumann Boundary Conditions set to 0 
        //concentration[0]=(4*concentration[1]-concentration[2])/3.0;
        //concentration[n_xslots-1]=(4*concentration[n_xslots-2]-concentration[n_xslots-3])/3.0;
        // Save data of each time iteration

    // Update number of cells for the next iteration
   
    free(rd); 


      }

void PDE_Handler( struct all_parameters *params,
                double *concentration,
                long interface,
                double tau){
                   

     
    
    euler(params,concentration,interface,tau);
}