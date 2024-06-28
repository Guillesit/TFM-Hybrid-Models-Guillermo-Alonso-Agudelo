#include "Header.h"						
// #include <gsl/gsl_odeiv.h>
// #include <gsl/gsl_errno.h>

double Reactions(int n_spatial_slots,
				double **gilAux,  
				long *J,
				long **StochasticNcells,
				double **b, 
				double death_rate,
				double diffusion_R, 
				double diffusion_L
				) {
		//supongo que first es el primer índice (sería 0 si consideramos todo el vector)
		int i,j;
		double suma=0.0;
 
		for(j=0;j<J[0];j++){
				gilAux[0][j]=StochasticNcells[0][j]*death_rate;
				gilAux[0][j+J[0]]=StochasticNcells[0][j]*b[0][j];
				gilAux[0][j+2*J[0]]=StochasticNcells[0][j]*diffusion_R;
				gilAux[0][j+3*J[0]]=0.0;
				suma=suma+StochasticNcells[0][j]*(death_rate+b[0][j]+diffusion_R);
		}
						
		for(i=1;i<(n_spatial_slots-1);i++)
		{
			for(j=0;j<J[i];j++){
				gilAux[i][j]=StochasticNcells[i][j]*death_rate;
				gilAux[i][j+J[i]]=StochasticNcells[i][j]*b[i][j];
				gilAux[i][j+2*J[i]]=StochasticNcells[i][j]*diffusion_R;
				gilAux[i][j+3*J[i]]=StochasticNcells[i][j]*diffusion_L;
				suma=suma+StochasticNcells[i][j]*(death_rate+b[i][j]+diffusion_R+diffusion_L);
			}
		}
		
		for(j=0;j<J[n_spatial_slots-1];j++){
				gilAux[n_spatial_slots-1][j]=StochasticNcells[n_spatial_slots-1][j]*death_rate;
				gilAux[n_spatial_slots-1][j+J[i]]=StochasticNcells[n_spatial_slots-1][j]*b[n_spatial_slots-1][j];
				gilAux[n_spatial_slots-1][j+2*J[i]]=0.0;
				gilAux[n_spatial_slots-1][j+3*J[i]]=StochasticNcells[n_spatial_slots-1][j]*diffusion_L;
				suma=suma+StochasticNcells[i][j]*(death_rate+b[i][j]+diffusion_L);
		}

	return suma;
}


//Includes therapy effects!!!!
double Reactions_two_pop(int n_spatial_slots,
				double **gilAux,  
				long *J_host,
				long **StochasticNcells_host,
				double **b_host,
				long *J_inv,
				long **StochasticNcells_inv,
				double **b_inv, 
				double death_rate_host,
                double death_rate_invader,
				double diffusion_R, 
				double diffusion_L,
                double survival_rate //associated with a therapy, affects both populations
				) {
		//supongo que first es el primer índice (sería 0 si consideramos todo el vector)
		int i,j;
		double suma=0.0;
 // First box only can diffuse to the right
		// Save host propensity and sum
		for(j=0;j<J_host[0];j++){
				gilAux[0][j]=StochasticNcells_host[0][j]
                                *(death_rate_host+(1-survival_rate)*b_host[0][j]);
				gilAux[0][j+J_host[0]]=StochasticNcells_host[0][j]*b_host[0][j]*survival_rate;
				gilAux[0][j+2*J_host[0]]=StochasticNcells_host[0][j]*diffusion_R;
				gilAux[0][j+3*J_host[0]]=0.0;
				suma=suma+StochasticNcells_host[0][j]*(death_rate_host+b_host[0][j]+diffusion_R);
		}
    
		// Save invader propensity and sum
		for(j=0;j<J_inv[0];j++){
				gilAux[0][j+4*J_host[0]]=StochasticNcells_inv[0][j]
                                            *(death_rate_invader+(1-survival_rate)*b_inv[0][j]);
				gilAux[0][j+J_inv[0]+4*J_host[0]]=StochasticNcells_inv[0][j]*b_inv[0][j]*survival_rate;
				gilAux[0][j+2*J_inv[0]+4*J_host[0]]=StochasticNcells_inv[0][j]*diffusion_R;
				gilAux[0][j+3*J_inv[0]+4*J_host[0]]=0.0;
				suma=suma+StochasticNcells_inv[0][j]*(death_rate_invader+b_inv[0][j]+diffusion_R);
		}
    
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
		for(i=1;i<(n_spatial_slots-1);i++)
		{
			// Save host propensity and sum
			for(j=0;j<J_host[i];j++){
				gilAux[i][j]=StochasticNcells_host[i][j]
                                *(death_rate_host+(1-survival_rate)*b_host[i][j]);
				gilAux[i][j+J_host[i]]=StochasticNcells_host[i][j]*b_host[i][j]*survival_rate;
				gilAux[i][j+2*J_host[i]]=StochasticNcells_host[i][j]*diffusion_R;
				gilAux[i][j+3*J_host[i]]=StochasticNcells_host[i][j]*diffusion_L;
				suma=suma+StochasticNcells_host[i][j]*(death_rate_host+b_host[i][j]+diffusion_R+diffusion_L);
			}

			// Save invader propensity and sum
			for(j=0;j<J_inv[i];j++){
				gilAux[i][j+4*J_host[i]]=StochasticNcells_inv[i][j]
                                *(death_rate_invader+(1-survival_rate)*b_inv[i][j]);
				gilAux[i][j+J_inv[i]+4*J_host[i]]=StochasticNcells_inv[i][j]*b_inv[i][j]*survival_rate;
				gilAux[i][j+2*J_inv[i]+4*J_host[i]]=StochasticNcells_inv[i][j]*diffusion_R;
				gilAux[i][j+3*J_inv[i]+4*J_host[i]]=StochasticNcells_inv[i][j]*diffusion_L;
				suma=suma+StochasticNcells_inv[i][j]*(death_rate_invader+b_inv[i][j]+diffusion_R+diffusion_L);
			}
			
		}
// Last box only can diffuse to the Left
		// Save host propensity and sum		
		for(j=0;j<J_host[n_spatial_slots-1];j++){
				gilAux[n_spatial_slots-1][j]=StochasticNcells_host[n_spatial_slots-1][j]
                    *(death_rate_host+(1-survival_rate)*b_host[n_spatial_slots-1][j]);
				gilAux[n_spatial_slots-1][j+J_host[i]]=StochasticNcells_host[n_spatial_slots-1][j]*b_host[n_spatial_slots-1][j]*survival_rate;
				gilAux[n_spatial_slots-1][j+2*J_host[i]]=0.0;
				gilAux[n_spatial_slots-1][j+3*J_host[i]]=StochasticNcells_host[n_spatial_slots-1][j]*diffusion_L;
				suma=suma+StochasticNcells_host[i][j]*(death_rate_host+b_host[i][j]+diffusion_L);
		}

		// Save invader propensity and sum
        //actually i=n_spatial_slots-1 when we reach this part
		for(j=0;j<J_inv[n_spatial_slots-1];j++){
				gilAux[n_spatial_slots-1][j+4*J_host[i]]=
                    StochasticNcells_inv[n_spatial_slots-1][j]
                        *(death_rate_invader+(1-survival_rate)*b_inv[n_spatial_slots-1][j]);
				gilAux[n_spatial_slots-1][j+J_inv[i]+4*J_host[i]]=StochasticNcells_inv[n_spatial_slots-1][j]*b_inv[n_spatial_slots-1][j]*survival_rate;
				gilAux[n_spatial_slots-1][j+2*J_inv[i]+4*J_host[i]]=0.0;
				gilAux[n_spatial_slots-1][j+3*J_inv[i]+4*J_host[i]]=StochasticNcells_inv[n_spatial_slots-1][j]*diffusion_L;
				suma=suma+StochasticNcells_inv[i][j]*(death_rate_invader+b_inv[i][j]+diffusion_L);
		}

	return suma;
}


void Gillespie(int i,
			int j, 
			long **StochasticNcells,
			long *J,
			double **StochasticAge, 
			double **b
			){
	int k;
	
	if(j<J[i])
		{
			StochasticNcells[i][j]=StochasticNcells[i][j]-1;
			if(StochasticNcells[i][j]==0)
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
		else if(j<2*J[i])
		{
			StochasticNcells[i][j-J[i]]=StochasticNcells[i][j-J[i]]-1;
			if(StochasticNcells[i][j-J[i]]==0)
				{	
					for(k=(j-J[i]);k<J[i];k++)
					{
						StochasticNcells[i][k]=StochasticNcells[i][k+1];
						StochasticAge[i][k]=StochasticAge[i][k+1];
						b[i][k]=b[i][k+1];
					}
					J[i]=J[i]-1;
			}
			J[i]=J[i]+1;
			for(k=J[i]-1;k>0;k=k-1)
			{
				StochasticNcells[i][k]=StochasticNcells[i][k-1];
				StochasticAge[i][k]=StochasticAge[i][k-1];
				b[i][k]=b[i][k-1];
			}
			StochasticAge[i][0]=0.0;
			StochasticNcells[i][0]=2.0;
			b[i][0]=0.0;	
		}
		else if(j<3*J[i]) //R
		{
			StochasticNcells[i][j-2*J[i]]=StochasticNcells[i][j-2*J[i]]-1;
			J[i+1]=J[i+1]+1;
			for(k=J[i+1]-1;k>0;k=k-1)
			{
				StochasticNcells[i+1][k]=StochasticNcells[i+1][k-1];
				StochasticAge[i+1][k]=StochasticAge[i+1][k-1];
				b[i+1][k]=b[i+1][k-1];
			}
			StochasticAge[i+1][0]=StochasticAge[i][j-2*J[i]];
			StochasticNcells[i+1][0]=1.0;
			b[i+1][0]=b[i][j-2*J[i]];
			if(StochasticNcells[i][j-2*J[i]]==0)
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
		else //L
		{
			StochasticNcells[i][j-3*J[i]]=StochasticNcells[i][j-3*J[i]]-1;
			J[i-1]=J[i-1]+1;
			for(k=J[i-1]-1;k>0;k=k-1)
			{
				StochasticNcells[i-1][k]=StochasticNcells[i-1][k-1];
				StochasticAge[i-1][k]=StochasticAge[i-1][k-1];
				b[i-1][k]=b[i-1][k-1];
			}
			StochasticAge[i-1][0]=StochasticAge[i][j-3*J[i]];
			StochasticNcells[i-1][0]=1.0;
			b[i-1][0]=b[i][j-3*J[i]];
			if(StochasticNcells[i][j-3*J[i]]==0)
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






