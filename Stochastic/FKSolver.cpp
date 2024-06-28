
#include "Header.h"


void eval_c(struct all_parameters *params, double *concentration, double *rd)
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

    for(x_slot=1;x_slot<n_xslots-1;x_slot++){
        reaction = r * concentration[x_slot] - (r/K)*concentration[x_slot]*concentration[x_slot]; 
        // be careful when rescaling we need to rescale also delta_x and the dif coef 
        diffusion = (concentration[x_slot-1] + concentration[x_slot+1] - 2*concentration[x_slot])/(delta_x*delta_x); 
        rd[x_slot] = reaction + D*diffusion; 
    }
    
    }
   

void euler(struct all_parameters *params, 
                double *concentration
                    ){
    long i = 0; 
    double *rd; 
    long n_xslots;  
    double delta_t; 
    struct all_parameters *p;
    p= (struct all_parameters *) params;
    delta_t = (p->delta_t); 
    n_xslots=(p->n_xslots);
    
     
      // Allocate memory for rd
    rd = (double *)malloc(n_xslots * sizeof(double));
    if (rd == NULL) {
        fprintf(stderr, "Error: memory allocation failed\n");
        exit(1);
    }
    eval_c(p, concentration,rd) ; 
   
    // Rewrite concentration with the concentration of the new time step 
    for(i=1;i<n_xslots-1;i++){
        // Advance in space except for the boundaries 
        // calculate the function using the previous step
        concentration[i] = concentration[i] + delta_t * rd[i];    
        }
        // Neumann Boundary Conditions set to 0 
        concentration[0]=(4*concentration[1]-concentration[2])/3.0;
        concentration[n_xslots-1]=(4*concentration[n_xslots-2]-concentration[n_xslots-3])/3.0;
        // Save data of each time iteration

    // Update number of cells for the next iteration
   
    free(rd); 


      }

void PDE_Handler( struct all_parameters *params,
                double *concentration){
                   
    struct all_parameters *p;
    p= (struct all_parameters *) params;
    double delta_t; 
    delta_t = (p->delta_t); 
    euler(p,concentration);
}