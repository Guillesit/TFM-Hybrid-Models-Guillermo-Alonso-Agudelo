#include "Header.h"						
// #include <gsl/gsl_odeiv.h>
// #include <gsl/gsl_errno.h>




double Reactions_invader( //sum of propensities for the population at the left
						int n_spatial_slots,
                        double **Propensities_invader,
                        long *J_inv,
                        long **StochasticNcells_inv,
                        double **b_inv,
                        long interface_invader,
                        double death_rate_invader,
                        double diffusion_invader,
                        double survival_rate  //associated with a therapy, affects both populations
                        ){
		//0=death, 1=birth, 2=right diffusion, 3=left diffusion
		
		int i,j;
		double suma=0.0;
		//The host is at the right, the invader is at the left.
		//Recall that, at the interface, the diffusion rate towards the deterministic domain is capped


		for(j=0;j<J_inv[interface_invader];j++){
				Propensities_invader[interface_invader][j]=
							StochasticNcells_inv[interface_invader][j]
                                            *(death_rate_invader+(1.0-survival_rate)*b_inv[interface_invader][j]);
				Propensities_invader[interface_invader][j+J_inv[interface_invader]]
								=StochasticNcells_inv[interface_invader][j]*b_inv[interface_invader][j]*survival_rate;
				Propensities_invader[interface_invader][j+2*J_inv[interface_invader]]
								=StochasticNcells_inv[interface_invader][j]*diffusion_invader;
				Propensities_invader[interface_invader][j+3*J_inv[interface_invader]]=0.0;
				suma=suma+StochasticNcells_inv[interface_invader][j]*(death_rate_invader+b_inv[interface_invader][j]+diffusion_invader);
		};


		//CAjas interiores
		for(i=1+interface_invader;i<(n_spatial_slots-1);i++)
		{

			// Save invader propensity and sum
			for(j=0;j<J_inv[i];j++){
				Propensities_invader[i][j]=StochasticNcells_inv[i][j]
                                *(death_rate_invader+(1.0-survival_rate)*b_inv[i][j]);
				Propensities_invader[i][j+J_inv[i]]=StochasticNcells_inv[i][j]*b_inv[i][j]*survival_rate;
				Propensities_invader[i][j+2*J_inv[i]]=StochasticNcells_inv[i][j]*diffusion_invader;
				Propensities_invader[i][j+3*J_inv[i]]=StochasticNcells_inv[i][j]*diffusion_invader;
				suma=suma+StochasticNcells_inv[i][j]*(death_rate_invader+b_inv[i][j]+2*diffusion_invader);
			}
			
		};

		//actually i=n_spatial_slots-1 when we reach this part
		for(j=0;j<J_inv[n_spatial_slots-1];j++){
				Propensities_invader[n_spatial_slots-1][j]=
                    StochasticNcells_inv[n_spatial_slots-1][j]
                        *(death_rate_invader+(1.0-survival_rate)*b_inv[n_spatial_slots-1][j]);
				Propensities_invader[n_spatial_slots-1][j+J_inv[n_spatial_slots-1]]=
							StochasticNcells_inv[n_spatial_slots-1][j]*b_inv[n_spatial_slots-1][j]*survival_rate;
				Propensities_invader[n_spatial_slots-1][j+2*J_inv[n_spatial_slots-1]]=0.0;
				Propensities_invader[n_spatial_slots-1][j+3*J_inv[n_spatial_slots-1]]=
							StochasticNcells_inv[n_spatial_slots-1][j]*diffusion_invader;
				suma=suma+StochasticNcells_inv[n_spatial_slots-1][j]*(death_rate_invader+b_inv[n_spatial_slots-1][j]+diffusion_invader);
		};

	return suma;

}
//Notes: (i) Includes therapy effects!!!!
//(ii) En verdad hacer una actualizacion tan exhaustiva solo parece merecer la pena la primera vez que se entra en el bucle
//Para las siguientes deberia bastar actualizar solo las casillas afectadas por el canal que disparo



/*******/
double Reactions_host( //sum of propensities for the population at the right
						int n_spatial_slots,
                        double **Propensities_host,
                        long *J_host,
                        long **StochasticNcells_host,
                        double **b_host,
                        long interface_host,
                        double death_rate_host,
                        double diffusion_host,
                        double survival_rate //associated with a therapy, affects both populations
                        ){

		//0=death, 1=birth, 2=right diffusion, 3=left diffusion
		
		int i,j;
		double suma=0.0;
		//The host is at the right, the invader is at the left.
		//Recall that, at the interface, the diffusion rate towards the deterministic domain is capped

		// Leftmost box diffuse to the right
		// Save host propensity and sum
		for(j=0;j<J_host[0];j++){
				Propensities_host[0][j]=
						StochasticNcells_host[0][j]
                                *(death_rate_host+(1-survival_rate)*b_host[0][j]);
				Propensities_host[0][j+J_host[0]]=
						StochasticNcells_host[0][j]*b_host[0][j]*survival_rate;
				Propensities_host[0][j+2*J_host[0]]=
						StochasticNcells_host[0][j]*diffusion_host;
				Propensities_host[0][j+3*J_host[0]]=0.0;
				suma=suma+StochasticNcells_host[0][j]*(death_rate_host+b_host[0][j]+diffusion_host);
		}; //end for leftmost location

		for(i=1;i<interface_host;i++)
		{
			// Save host propensity and sum
			for(j=0;j<J_host[i];j++){
				Propensities_host[i][j]=StochasticNcells_host[i][j]
                                *(death_rate_host+(1-survival_rate)*b_host[i][j]);
				Propensities_host[i][j+J_host[i]]=StochasticNcells_host[i][j]*b_host[i][j]*survival_rate;
				Propensities_host[i][j+2*J_host[i]]=StochasticNcells_host[i][j]*diffusion_host;
				Propensities_host[i][j+3*J_host[i]]=StochasticNcells_host[i][j]*diffusion_host;
				suma=suma+StochasticNcells_host[i][j]*(death_rate_host+b_host[i][j]+diffusion_host+diffusion_host);
			};
		}; // end for inner locations	

		//Rightmost box only can diffuse to the left
		// Save host propensity and sum		
		
		for(j=0;j<J_host[interface_host];j++){
				Propensities_host[interface_host][j]=StochasticNcells_host[interface_host][j]
                    *(death_rate_host+(1-survival_rate)*b_host[interface_host][j]);
				Propensities_host[interface_host][j+J_host[interface_host]]=StochasticNcells_host[interface_host][j]
						*b_host[interface_host][j]*survival_rate;
				Propensities_host[interface_host][j+2*J_host[interface_host]]=0.0;
				Propensities_host[interface_host][j+3*J_host[interface_host]]=StochasticNcells_host[interface_host][j]*diffusion_host;
				suma=suma+StochasticNcells_host[interface_host][j]*(death_rate_host+b_host[interface_host][j]+diffusion_host);
		};

	return suma;
}
//Notes: (i) Includes therapy effects!!!!
//(ii) En verdad hacer una actualizacion tan exhaustiva solo parece merecer la pena la primera vez que se entra en el bucle
//Para las siguientes deberia bastar actualizar solo las casillas afectadas por el canal que disparo

/***********/
double Reactions_two_pop(int n_spatial_slots,
				double **gilAux,  
				long *J_host,
				long **StochasticNcells_host,
				double **b_host,
				long *J_inv,
				long **StochasticNcells_inv,
				double **b_inv, 
				long interface_host,
				long interface_invader,
				double death_rate_host,
                double death_rate_invader,
				double diffusion_host, 
				double diffusion_invader,
                double survival_rate //associated with a therapy, affects both populations
				) {
	//0=death, 1=birth, 2=right diffusion, 3=left diffusion
		
		int i,j;
		double suma=0.0;
		//The host is at the right, the invader is at the left.
		//Recall that, at the interface, the diffusion rate towards the deterministic domain is capped


		//MAnaging the host (DONE)

 		// Leftmost box diffuse to the right
		// Save host propensity and sum
		for(j=0;j<J_host[0];j++){
				gilAux[0][j]=
						StochasticNcells_host[0][j]
                                *(death_rate_host+(1-survival_rate)*b_host[0][j]);
				gilAux[0][j+J_host[0]]=
						StochasticNcells_host[0][j]*b_host[0][j]*survival_rate;
				gilAux[0][j+2*J_host[0]]=
						StochasticNcells_host[0][j]*diffusion_host;
				gilAux[0][j+3*J_host[0]]=0.0;
				suma=suma+StochasticNcells_host[0][j]*(death_rate_host+b_host[0][j]+diffusion_host);
		}; //end for leftmost location

		for(i=1;i<interface_host;i++)
		{
			// Save host propensity and sum
			for(j=0;j<J_host[i];j++){
				gilAux[i][j]=StochasticNcells_host[i][j]
                                *(death_rate_host+(1-survival_rate)*b_host[i][j]);
				gilAux[i][j+J_host[i]]=StochasticNcells_host[i][j]*b_host[i][j]*survival_rate;
				gilAux[i][j+2*J_host[i]]=StochasticNcells_host[i][j]*diffusion_host;
				gilAux[i][j+3*J_host[i]]=StochasticNcells_host[i][j]*diffusion_host;
				suma=suma+StochasticNcells_host[i][j]*(death_rate_host+b_host[i][j]+diffusion_host+diffusion_host);
			};
		}; // end for inner locations	

		//Rightmost box only can diffuse to the left
		// Save host propensity and sum		
		
		for(j=0;j<J_host[interface_host];j++){
				gilAux[interface_host][j]=StochasticNcells_host[interface_host][j]
                    *(death_rate_host+(1-survival_rate)*b_host[interface_host][j]);
				gilAux[interface_host][j+J_host[interface_host]]=StochasticNcells_host[interface_host][j]
						*b_host[interface_host][j]*survival_rate;
				gilAux[interface_host][j+2*J_host[interface_host]]=0.0;
				gilAux[interface_host][j+3*J_host[interface_host]]=StochasticNcells_host[interface_host][j]*diffusion_host;
				suma=suma+StochasticNcells_host[interface_host][j]*(death_rate_host+b_host[interface_host][j]+diffusion_host);
		};
			
		//printf("sumahost=%lf\n",suma);


		/********/
		//Manage invader (DONE)
		
		for(j=0;j<J_inv[interface_invader];j++){
				gilAux[interface_invader][j+4*J_host[interface_invader]]=StochasticNcells_inv[interface_invader][j]
                                            *(death_rate_invader+(1-survival_rate)*b_inv[interface_invader][j]);
				gilAux[interface_invader][j+J_inv[interface_invader]+4*J_host[interface_invader]]
								=StochasticNcells_inv[interface_invader][j]*b_inv[interface_invader][j]*survival_rate;
				gilAux[interface_invader][j+2*J_inv[0]+4*J_host[interface_invader]]
								=StochasticNcells_inv[interface_invader][j]*diffusion_invader;
				gilAux[interface_invader][j+3*J_inv[0]+4*J_host[interface_invader]]=0.0;
				suma=suma+StochasticNcells_inv[interface_invader][j]*(death_rate_invader+b_inv[interface_invader][j]+diffusion_invader);
		};
		
    
    //Tambien lo puedo hacer asi:
    /*
    for(j=4*J_host[0];j<J_inv[0]+4*J_host[0];j++){ //Pues si lo hago asi funciona (?)
        gilAux[0][j]=StochasticNcells_inv[0][j-4*J_host[0]]*death_rate;
        gilAux[0][j+J_inv[0]]=StochasticNcells_inv[0][j-4*J_host[0]]*b_inv[0][j-4*J_host[0]];
        gilAux[0][j+2*J_inv[0]]=StochasticNcells_inv[0][j-4*J_host[0]]*diffusion_R;
        gilAux[0][j+3*J_inv[0]]=0.0;
        suma=suma+StochasticNcells_inv[0][j-4*J_host[0]]*(death_rate+b_inv[0][j-4*J_host[0]]+diffusion_R);
    }
     */
    //La clave es que el indice del gilaux no va correlativo al de la segunda poblacion
    //(que empieza de cero)
						
    
    
     //CAjas interiores
		for(i=1+interface_invader;i<(n_spatial_slots-1);i++)
		{

			// Save invader propensity and sum
			for(j=0;j<J_inv[i];j++){
				gilAux[i][j+4*J_host[i]]=StochasticNcells_inv[i][j]
                                *(death_rate_invader+(1-survival_rate)*b_inv[i][j]);
				gilAux[i][j+J_inv[i]+4*J_host[i]]=StochasticNcells_inv[i][j]*b_inv[i][j]*survival_rate;
				gilAux[i][j+2*J_inv[i]+4*J_host[i]]=StochasticNcells_inv[i][j]*diffusion_invader;
				gilAux[i][j+3*J_inv[i]+4*J_host[i]]=StochasticNcells_inv[i][j]*diffusion_invader;
				suma=suma+StochasticNcells_inv[i][j]*(death_rate_invader+b_inv[i][j]+diffusion_invader+diffusion_invader);
			}
			
		};


		// Save invader propensity and sum
        //actually i=n_spatial_slots-1 when we reach this part
		for(j=0;j<J_inv[n_spatial_slots-1];j++){
				gilAux[n_spatial_slots-1][j+4*J_host[i]]=
                    StochasticNcells_inv[n_spatial_slots-1][j]
                        *(death_rate_invader+(1-survival_rate)*b_inv[n_spatial_slots-1][j]);
				gilAux[n_spatial_slots-1][j+J_inv[i]+4*J_host[i]]=StochasticNcells_inv[n_spatial_slots-1][j]*b_inv[n_spatial_slots-1][j]*survival_rate;
				gilAux[n_spatial_slots-1][j+2*J_inv[i]+4*J_host[i]]=0.0;
				gilAux[n_spatial_slots-1][j+3*J_inv[i]+4*J_host[i]]=StochasticNcells_inv[n_spatial_slots-1][j]*diffusion_invader;
				suma=suma+StochasticNcells_inv[i][j]*(death_rate_invader+b_inv[i][j]+diffusion_invader);
		};

		//printf("sumatot=%lf\n",suma);

	return suma;
}


/***********/


void Gillespie(int i,
			int j, 
			long **StochasticNcells,
			long *J,
			double **StochasticAge, 
			double **b
			){
	int k;
	
	if(j<J[i]) //death case
		{
				//printf("death\n");
			StochasticNcells[i][j]=StochasticNcells[i][j]-1;
			if(StochasticNcells[i][j]==0) //remove that slot-age pair (case that no cells with the given age remain)
				{	
					for(k=j;k<J[i];k++)
					{
						StochasticNcells[i][k]=StochasticNcells[i][k+1];
						StochasticAge[i][k]=StochasticAge[i][k+1];
						b[i][k]=b[i][k+1];
					}
					J[i]=J[i]-1;
				}
		}
		else if(j<2*J[i]) //birth case
		{
			//printf("birth\n");
			StochasticNcells[i][j-J[i]]=StochasticNcells[i][j-J[i]]-1;
			if(StochasticNcells[i][j-J[i]]==0)  //remove that slot-age pair (case that no cells with the given age remain)
				{	
					for(k=(j-J[i]);k<J[i];k++)
					{
						StochasticNcells[i][k]=StochasticNcells[i][k+1];
						StochasticAge[i][k]=StochasticAge[i][k+1];
						b[i][k]=b[i][k+1];
					}
					J[i]=J[i]-1;
			}
			J[i]=J[i]+1; //enlarge the structure to accomodate newborns
			for(k=J[i]-1;k>0;k=k-1)
			{
				StochasticNcells[i][k]=StochasticNcells[i][k-1];
				StochasticAge[i][k]=StochasticAge[i][k-1];
				b[i][k]=b[i][k-1];
			}
			StochasticAge[i][0]=0.0;
			StochasticNcells[i][0]=2.0;  //OPTIMIZE (insert new cells at the end of the present vector)
			b[i][0]=0.0;	
		}
		else if(j<3*J[i]) //diffusion to the right
		{
			//printf("rdiff\n");
			StochasticNcells[i][j-2*J[i]]=StochasticNcells[i][j-2*J[i]]-1;
			J[i+1]=J[i+1]+1;  //enlarge the structure to accomodate the coming particle
			for(k=J[i+1]-1;k>0;k=k-1)
			{
				StochasticNcells[i+1][k]=StochasticNcells[i+1][k-1];
				StochasticAge[i+1][k]=StochasticAge[i+1][k-1];
				b[i+1][k]=b[i+1][k-1];
			}
			StochasticAge[i+1][0]=StochasticAge[i][j-2*J[i]];  //OPTIMIZE (insert new cells at the end of the present vector)
			StochasticNcells[i+1][0]=1.0;
			b[i+1][0]=b[i][j-2*J[i]];
			if(StochasticNcells[i][j-2*J[i]]==0) //remove that slot-age pair (case that no cells with the given age remain)
				{	
					for(k=(j-2*J[i]);k<J[i];k++)
					{
						StochasticNcells[i][k]=StochasticNcells[i][k+1];
						StochasticAge[i][k]=StochasticAge[i][k+1];
						b[i][k]=b[i][k+1];
					}
					J[i]=J[i]-1;
			}
		}
		else //diffusion to the left
		{
			//printf("ldiff\n");
			StochasticNcells[i][j-3*J[i]]=StochasticNcells[i][j-3*J[i]]-1;
			J[i-1]=J[i-1]+1;
			for(k=J[i-1]-1;k>0;k=k-1)
			{
				StochasticNcells[i-1][k]=StochasticNcells[i-1][k-1];
				StochasticAge[i-1][k]=StochasticAge[i-1][k-1];
				b[i-1][k]=b[i-1][k-1];
			}
			StochasticAge[i-1][0]=StochasticAge[i][j-3*J[i]];  //OPTIMIZE (insert new cells at the end of the present vector)
			StochasticNcells[i-1][0]=1.0;
			b[i-1][0]=b[i][j-3*J[i]];
			if(StochasticNcells[i][j-3*J[i]]==0)  //remove that slot-age pair (case that no cells with the given age remain)
				{	
					for(k=(j-3*J[i]);k<J[i];k++)
					{
						StochasticNcells[i][k]=StochasticNcells[i][k+1];
						StochasticAge[i][k]=StochasticAge[i][k+1];
						b[i][k]=b[i][k+1];
					}
					J[i]=J[i]-1;
			}
		}
}






//This one is superseded by the stuff in "Utilities"
double birth(double ox){//para (p6/p3)=1.0
	
	double a_=8250;
	//double Ccr=0.02266;
	double Ccr;
	Ccr=1-0.4*log(1/0.017*(2.02-1/(1-0.45185)));
	if(ox>Ccr)
	{
		return a_*pow(ox/Ccr-1,-0.2);
	}
	else
	{
		return 1000000000;//9
	}
}






