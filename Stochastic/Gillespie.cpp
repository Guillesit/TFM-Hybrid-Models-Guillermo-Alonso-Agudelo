

#include "Header.h"





//SUBROUTINES

void Initialize_propensities(struct all_parameters *params,
                            struct age_structure **StochasticAges,
                            double *oxygen_concentration
                            ){

	long i, ref_index;
  double *division_threshold;
  struct all_parameters *p;

  //double death_rate;
  //double tau_p;
  //double p6_over_p3;
  //double aplus;
  double diff_coef_pop;
  double r;
  double K; 
  //double survival_rate;
  long n_xslots;

  p= (struct all_parameters *) params;

  //death_rate=(p->death_rate_hat);

  //death_rate=0.1;

  //tau_p=(p->tau_p);
  //aplus=(p->aplus);
  //p6_over_p3=(p->p6_over_p3);
  diff_coef_pop=(p->D)/((p->delta_x)*(p->delta_x));
  r=(p->r);
  K=(p->K);
  //survival_rate=(p->survival_rate);
  n_xslots=(p->n_xslots);

  if((division_threshold= (double *) malloc(sizeof(double)*n_xslots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
  };

       //loop through the structure array
	for (i = 0; i < n_xslots; i++){


        //Set all propensities
        //0=birth, 1=death, 2=left diffusion, 3=right diffusion
        (StochasticAges[i]->propensity_vector)[BIRTH]=(StochasticAges[i]->ncells)*r;//tau_p;
        (StochasticAges[i]->propensity_vector)[DEATH]=(StochasticAges[i]->ncells)*((StochasticAges[i]->ncells)-1)*r/K;
        if(i==0){ //cannot diffuse to the left at the left end
            (StochasticAges[i]->propensity_vector)[LDIFF]=0.0;
        }
        else{
            (StochasticAges[i]->propensity_vector)[LDIFF]=(StochasticAges[i]->ncells)*diff_coef_pop;
        };
        if(i==n_xslots-1){ //cannot diffuse to the right at the right end
            (StochasticAges[i]->propensity_vector)[RDIFF]=0.0;
        }
        else{
            (StochasticAges[i]->propensity_vector)[RDIFF]=(StochasticAges[i]->ncells)*diff_coef_pop;
        };


	};//end spatial for
  free(division_threshold);

};//end Initialize_propensities



double Compute_total_propensity(struct age_structure **StochasticAges,
                                long n_xslots){

    long i=0;
    double total_propensity=0.0;

    for (i = 0; i < n_xslots; i++){
        total_propensity=total_propensity
                          +(StochasticAges[i]->propensity_vector)[BIRTH] //birth
                          +(StochasticAges[i]->propensity_vector)[DEATH] //death
                          +(StochasticAges[i]->propensity_vector)[LDIFF] //left_diffusion
                          +(StochasticAges[i]->propensity_vector)[RDIFF]; //right_diffusion
    };

    if(total_propensity<=0.0){
      fprintf(stderr,"Error: nonpositive propensity computed at time t=...\n");
      fprintf(stderr,"Aborting execution\n");
      exit(1);
    };

    return(total_propensity);

};// end Compute_total_propensity



void Determine_event(struct all_parameters *params,
                      long *pfiring_slot,
                      int *pfiring_reaction,
                      long *pfiring_cell,
                      struct age_structure ** StochasticAges,
                      double fraction_of_total_propensity
                      ){

    long i=0; //spatial location
    int j=0; //reaction tag, 0=birth, 1=death, 2=left diffusion, 3=right diffusion
    long k=0; //index for firing cell
    double sum=0.0;


    struct all_parameters *p;

    //double death_rate;
    //double tau_p;
    double diff_coef_pop;
    double r;
    double K; 
    //long n_xslots;

    p= (struct all_parameters *) params;
    //death_rate=(p->death_rate_hat);
    //tau_p=(p->tau_p);
    diff_coef_pop=(p->D)/((p->delta_x)*(p->delta_x));
    r=(p->r);
    K=(p->K);
    //n_xslots=(p->n_xslots);


    //First pass: determine firing_slot and firing_reaction
     
    sum=(StochasticAges[0]->propensity_vector)[0];
    
    while(sum<fraction_of_total_propensity){
        j++;
        if(j==NUMBER_OF_REACTIONS){
            i=i+1;
            j=0;
        };
        sum=sum+(StochasticAges[i]->propensity_vector)[j];
        
    };//end while
    
    *pfiring_slot=i; // compartment where reaction takes place
    *pfiring_reaction=j; //reaction tag


    //////////////
    //Second pass: determine firing_cell

    sum=sum-(StochasticAges[i]->propensity_vector)[j];
    
    switch(j){ //doing things like case '0': causes it to fail
            //removing the literals, no itoa or the like are needed

      case BIRTH: //birth 0
        k=0;
        sum=sum+r;// /tau_p;
        while(sum<fraction_of_total_propensity){
          k++; 
          sum=sum+r;// /tau_p;                                          //REVISAR
        };
        if(k>=(StochasticAges[i]->ncells)){ //hopefully not needed
          fprintf(stderr,"Something went wrong when determining firing_cell (birth)\n");
          exit(1);
        };
        *pfiring_cell=k;
        break;


      case DEATH: //death 1
        k=0;
        
        sum=sum+((StochasticAges[i]->ncells)-1)*r/K;
        while(sum<fraction_of_total_propensity){
          k++; 
          sum=sum+((StochasticAges[i]->ncells)-1)*r/K;
        };

        if(k>=(StochasticAges[i]->ncells)){ //hopefully not needed
          fprintf(stderr,"Something went wrong when determining firing_cell (death)\n");
          exit(1);
        };
        *pfiring_cell=k;
        break;
        

      case LDIFF: //left_diffusion 2
        k=0;
        sum=sum+diff_coef_pop;
        while(sum<fraction_of_total_propensity){
          k++; 
          sum=sum+diff_coef_pop;
        };

        if(k>=(StochasticAges[i]->ncells)){ //hopefully not needed
          fprintf(stderr,"Something went wrong when determining firing_cell (left_diff)\n");
          exit(1);
        };
        *pfiring_cell=k;
        break;
        

      case RDIFF: //right_diffusion 3
        k=0;
        sum=sum+diff_coef_pop;
        while(sum<fraction_of_total_propensity){
          k++; 
          sum=sum+diff_coef_pop;
        };

        if(k>=(StochasticAges[i]->ncells)){ //hopefully not needed
          fprintf(stderr,"Something went wrong when determining firing_cell (right_diff)\n");
          exit(1);
        };
        *pfiring_cell=k;
        break;
        

      default: 
        printf("We should not have entered here when determining firing_cell\n");
        exit(1);     

    }; //end switch j

}; //end Determine_event




void Advance_ages(struct age_structure **StochasticAges,
                    long n_xslots,
                    double tau){

      long i, j;

      for (i = 0; i < n_xslots;i++){
          for (j = 0; j < (StochasticAges[i]->ncells);j++)
          {
              (StochasticAges[i]->age_distribution)[j]+=tau;
          };
      };

};//end Advance_ages



void Remove_cell(double *vector,
                long location,
                long n_elements){

  long i;

  for (i = location; i < n_elements-1; i++){
      vector[i]=vector[i+1];
  };
  vector[n_elements-1]=-1.0;//hope this failsafe works for the age threshold recomputation issues
  //Anyhow, this should be controlled somewhere else (uncontrolled access to garbage)

  //NOTE: we do not update n_elements here. This should be done within the calling routine

};// end Remove_cell




void Insert_cell(double *vector,
                double age,
                long n_elements){

  long i;

  i = n_elements;
  while((age>vector[i-1])&&(i>=1)){	
    //printf("oki ");
      vector[i]=vector[i-1];
      //printf("doki\n");
      i--;
  };
  vector[i]=age;

  //NOTE: we do not update n_elements here. This should be done within the calling routine

};// end Insert_cell




//Birth propensities and thresholds are not updated here, but in the specific subroutine
void Handle_event(struct all_parameters *params,
                  struct age_structure **StochasticAges,
                  long firing_slot,    
                  int firing_reaction,
                  long firing_cell,
                  double *ncells_per_slot,
                  double r3){

    double age_aux;

    struct all_parameters *p;

    //double death_rate;
    double diff_coef_pop;
    double r;
    double K; 
    //double survival_rate;
    long n_xslots;

    p= (struct all_parameters *) params;
    //death_rate=(p->death_rate_hat);
  
    diff_coef_pop=(p->D)/((p->delta_x)*(p->delta_x));
    r=(p->r);
    K=(p->K);
    //survival_rate=(p->survival_rate);
    n_xslots=(p->n_xslots);

    switch(firing_reaction){ 

      case BIRTH: //birth 0

      /*
      if(r3>survival_rate){ //dividing cell is killed by therapy

        //Implementación alternativa: olvidar esto
        //y multiplicar la ratio de nacimiento por F_S
          
          //remove dying cell 
          Remove_cell((StochasticAges[firing_slot]->age_distribution),
                        firing_cell,
                        (StochasticAges[firing_slot]->ncells));

          (StochasticAges[firing_slot]->ncells)--;
          (ncells_per_slot[firing_slot])--;

          //update 3 propensities
          (StochasticAges[firing_slot]->propensity_vector)[DEATH]-=death_rate;
          if(firing_slot>0){(StochasticAges[firing_slot]->propensity_vector)[LDIFF]-=diff_coef_pop;};
          if(firing_slot<n_xslots-1){(StochasticAges[firing_slot]->propensity_vector)[RDIFF]-=diff_coef_pop;};
      
      }//end if (case of death by therapy)
      else{ //dividing cell proceeds as usual
      
          //remove dividing cell (no cambiar ncells aqui)
          Remove_cell((StochasticAges[firing_slot]->age_distribution),
                    firing_cell,
                    (StochasticAges[firing_slot]->ncells));

          //insert two new cells with zero age
          //AQUI NOS DA UNA VIOLACION DE SEGMENTO
          //LA AMPLIACION DE MEMORIA NO DEBE ESTAR FUNCIONANDO BIEN
          //Pero si el buffer es muy grande no da problema (no hay necesidad de ampliar memoria)
          //ALTERNATIVA: HACER LAS COMPROBACIONES DEL REALLOC AQUI EN VEZ DE EN LA RUTINA PRINCIPAL
          //Y TAMBIÉN EN LAS PARTES SIMILARES CORRESPONDIENTES A DIFUSION

          (StochasticAges[firing_slot]->age_distribution)[(StochasticAges[firing_slot]->ncells)-1]=0.0;
          (StochasticAges[firing_slot]->age_distribution)[(StochasticAges[firing_slot]->ncells)]=0.0;
          */
          //Es supuestamente equivalente a...
          //Insert_cell((StochasticAges[firing_slot]->age_distribution),
            //        0.0,
              //      (StochasticAges[firing_slot]->ncells-1));
          //Insert_cell((StochasticAges[firing_slot]->age_distribution),
            //        0.0,
              //      (StochasticAges[firing_slot]->ncells));

          (StochasticAges[firing_slot]->ncells)++;
          (ncells_per_slot[firing_slot])++;

          //update 3(4) propensities
          (StochasticAges[firing_slot]->propensity_vector)[BIRTH]+=r;
          (StochasticAges[firing_slot]->propensity_vector)[DEATH]=(StochasticAges[firing_slot]->ncells)*((StochasticAges[firing_slot]->ncells)-1)*r/K;//Recompute not as easy because of the n-1 inside the prop formula
          if(firing_slot>0){(StochasticAges[firing_slot]->propensity_vector)[LDIFF]+=diff_coef_pop;};
          if(firing_slot<n_xslots-1){(StochasticAges[firing_slot]->propensity_vector)[RDIFF]+=diff_coef_pop;};

      //};//end else (usual birth)
      
      break;


      //////
      ///////////////
      case DEATH: //death 1

      //remove dying cell
      /*
      Remove_cell((StochasticAges[firing_slot]->age_distribution),
                    firing_cell,
                    (StochasticAges[firing_slot]->ncells));
      */

      (StochasticAges[firing_slot]->ncells)--;
      (ncells_per_slot[firing_slot])--;

      //update 3 propensities
      (StochasticAges[firing_slot]->propensity_vector)[BIRTH]-=r;
      (StochasticAges[firing_slot]->propensity_vector)[DEATH]=(StochasticAges[firing_slot]->ncells)*((StochasticAges[firing_slot]->ncells)-1)*r/K;
      if(firing_slot>0){(StochasticAges[firing_slot]->propensity_vector)[LDIFF]-=diff_coef_pop;};
      if(firing_slot<n_xslots-1){(StochasticAges[firing_slot]->propensity_vector)[RDIFF]-=diff_coef_pop;};

      break;




      //////
      ////////////////////////
      case LDIFF: //left_diffusion 2
      
      //remove moving cell 
      age_aux=(StochasticAges[firing_slot]->age_distribution)[firing_cell];
      /*
      Remove_cell((StochasticAges[firing_slot]->age_distribution),
                    firing_cell,
                    (StochasticAges[firing_slot]->ncells));
      */
      (StochasticAges[firing_slot]->ncells)--;
      (ncells_per_slot[firing_slot])--;

      //insert at new destination 
      /*
      Insert_cell((StochasticAges[firing_slot-1]->age_distribution),
                    age_aux,
                    (StochasticAges[firing_slot-1]->ncells));
      */
      (StochasticAges[firing_slot-1]->ncells)++;
      (ncells_per_slot[firing_slot-1])++;

      //update 3+3 propensities (careful with boundaries)
      (StochasticAges[firing_slot]->propensity_vector)[BIRTH]-=r;
      (StochasticAges[firing_slot]->propensity_vector)[DEATH]=(StochasticAges[firing_slot]->ncells)*((StochasticAges[firing_slot]->ncells)-1)*r/K;
      (StochasticAges[firing_slot]->propensity_vector)[LDIFF]-=diff_coef_pop;
      if(firing_slot<n_xslots-1){(StochasticAges[firing_slot]->propensity_vector)[RDIFF]-=diff_coef_pop;};//TO CORRECT (why??)

      (StochasticAges[firing_slot-1]->propensity_vector)[BIRTH]+=r;
      (StochasticAges[firing_slot-1]->propensity_vector)[DEATH]=(StochasticAges[firing_slot-1]->ncells)*((StochasticAges[firing_slot-1]->ncells)-1)*r/K;
      if(firing_slot>1){(StochasticAges[firing_slot-1]->propensity_vector)[LDIFF]+=diff_coef_pop;};
      (StochasticAges[firing_slot-1]->propensity_vector)[RDIFF]+=diff_coef_pop;
      
      break;



      //////
      /////////////////////////
      case RDIFF: //right_diffusion 3
      
      //remove moving cell 
      age_aux=(StochasticAges[firing_slot]->age_distribution)[firing_cell];
      /*
      Remove_cell((StochasticAges[firing_slot]->age_distribution),
                    firing_cell,
                    (StochasticAges[firing_slot]->ncells));
      */
      (StochasticAges[firing_slot]->ncells)--;
      (ncells_per_slot[firing_slot])--;

      //insert at new destination 
      Insert_cell((StochasticAges[firing_slot+1]->age_distribution),
                    age_aux,
                    (StochasticAges[firing_slot+1]->ncells));

      (StochasticAges[firing_slot+1]->ncells)++;
      (ncells_per_slot[firing_slot+1])++;


      //update 3+3 propensities (careful with boundaries)
      (StochasticAges[firing_slot]->propensity_vector)[BIRTH]-=r;
      (StochasticAges[firing_slot]->propensity_vector)[DEATH]=(StochasticAges[firing_slot]->ncells)*((StochasticAges[firing_slot]->ncells)-1)*r/K;
      if(firing_slot>0){(StochasticAges[firing_slot]->propensity_vector)[LDIFF]-=diff_coef_pop;};
      (StochasticAges[firing_slot]->propensity_vector)[RDIFF]-=diff_coef_pop;

      (StochasticAges[firing_slot+1]->propensity_vector)[BIRTH]+=r;
      (StochasticAges[firing_slot+1]->propensity_vector)[DEATH]=(StochasticAges[firing_slot+1]->ncells)*((StochasticAges[firing_slot+1]->ncells)-1)*r/K;
      (StochasticAges[firing_slot+1]->propensity_vector)[LDIFF]+=diff_coef_pop;
      if(firing_slot<n_xslots-2){(StochasticAges[firing_slot+1]->propensity_vector)[RDIFF]+=diff_coef_pop;};
      
      break;

      ////////
      default:
      printf("We should not have entered here when handling firing event\n");
      exit(1); 

    };//end switch  

};//end Handle_event



void Update_birth_issues(struct all_parameters *params,
                        struct age_structure **StochasticAges,
                        double *oxygen_concentration
                        ){

  long i, ref_index;
  double *division_threshold;

  struct all_parameters *p;

  //double tau_p;
  //double p6_over_p3;
  //double aplus;
  long n_xslots;

  p= (struct all_parameters *) params;
  //tau_p=(p->tau_p);
  //aplus=(p->aplus);
  //p6_over_p3=(p->p6_over_p3);
  n_xslots=(p->n_xslots);



  if((division_threshold= (double *) malloc(sizeof(double)*n_xslots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
  };

  /* 
  //Recompute division threshold  
  Compute_division_threshold(p,
                              division_threshold,
                                oxygen_concentration,
                                n_xslots-1);

  */ 

  //loop through the structure array
  for (i = 0; i < n_xslots; i++){

  //printf("loc=%ld, ag1s_dyn=%lf\n",i,division_threshold[i]);

    //determining division index
    ref_index=0;//means no cell committed to divide

    while(((StochasticAges[i]->age_distribution)[ref_index]>=division_threshold[i])
        &&(ref_index<(StochasticAges[i]->ncells))){ //to avoid slipping out of range into the garbage
      ref_index++;
    };

    (StochasticAges[i]->ncells_ready_to_divide)=ref_index;
              //substract 1 to get the vector index of the youngest cell ready to divide
              //that would be -1 for no cells commited to divide

    //Reset birth propensities
    (StochasticAges[i]->propensity_vector)[BIRTH]=(StochasticAges[i]->ncells_ready_to_divide);// /tau_p;

  };//end spatial for

  free(division_threshold);

};//end Update_birth_issues

