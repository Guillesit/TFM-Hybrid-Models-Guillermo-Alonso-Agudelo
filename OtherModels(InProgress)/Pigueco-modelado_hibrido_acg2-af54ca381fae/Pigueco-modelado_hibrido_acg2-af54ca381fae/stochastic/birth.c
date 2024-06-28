#include <math.h>							
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

double reactions(int end, double **gilAux,  int *J, double **age, double **b, double nu, double R, double L) {
		//supongo que first es el primer índice (sería 0 si consideramos todo el vector)
		int i,j;
		double suma=0.0;
 
		for(j=0;j<J[0];j++){
				gilAux[0][j]=age[0][j]*nu;
				gilAux[0][j+J[0]]=age[0][j]*b[0][j];
				gilAux[0][j+2*J[0]]=age[0][j]*R;
				gilAux[0][j+3*J[0]]=0.0;
				suma=suma+age[0][j]*(nu+b[0][j]+R);
		}
						
		for(i=1;i<(end-1);i++)
		{
			for(j=0;j<J[i];j++){
				gilAux[i][j]=age[i][j]*nu;
				gilAux[i][j+J[i]]=age[i][j]*b[i][j];
				gilAux[i][j+2*J[i]]=age[i][j]*R;
				gilAux[i][j+3*J[i]]=age[i][j]*L;
				suma=suma+age[i][j]*(nu+b[i][j]+R+L);
			}
		}
		
		for(j=0;j<J[end-1];j++){
				gilAux[end-1][j]=age[end-1][j]*nu;
				gilAux[end-1][j+J[i]]=age[end-1][j]*b[end-1][j];
				gilAux[end-1][j+2*J[i]]=0.0;
				gilAux[end-1][j+3*J[i]]=age[end-1][j]*L;
				suma=suma+age[i][j]*(nu+b[i][j]+L);
		}

	return suma;
}


void gillespie(int i, int j, double **age, int *J, double **ag, double **b){
	int k;
	
	if(j<J[i])
		{
			age[i][j]=age[i][j]-1;
			if(age[i][j]==0)
				{	
					for(k=j;k<J[i];k++)
					{
						age[i][k]=age[i][k+1];
						ag[i][k]=ag[i][k+1];
						b[i][k]=b[i][k+1];
					}
					J[i]=J[i]-1;
				}
		}
		else if(j<2*J[i])
		{
			age[i][j-J[i]]=age[i][j-J[i]]-1;
			if(age[i][j-J[i]]==0)
				{	
					for(k=(j-J[i]);k<J[i];k++)
					{
						age[i][k]=age[i][k+1];
						ag[i][k]=ag[i][k+1];
						b[i][k]=b[i][k+1];
					}
					J[i]=J[i]-1;
			}
			J[i]=J[i]+1;
			for(k=J[i]-1;k>0;k=k-1)
			{
				age[i][k]=age[i][k-1];
				ag[i][k]=ag[i][k-1];
				b[i][k]=b[i][k-1];
			}
			ag[i][0]=0.0;
			age[i][0]=2.0;
			b[i][0]=0.0;	
		}
		else if(j<3*J[i]) //R
		{
			age[i][j-2*J[i]]=age[i][j-2*J[i]]-1;
			J[i+1]=J[i+1]+1;
			for(k=J[i+1]-1;k>0;k=k-1)
			{
				age[i+1][k]=age[i+1][k-1];
				ag[i+1][k]=ag[i+1][k-1];
				b[i+1][k]=b[i+1][k-1];
			}
			ag[i+1][0]=ag[i][j-2*J[i]];
			age[i+1][0]=1.0;
			b[i+1][0]=b[i][j-2*J[i]];
			if(age[i][j-2*J[i]]==0)
				{	
					for(k=(j-2*J[i]);k<J[i];k++)
					{
						age[i][k]=age[i][k+1];
						ag[i][k]=ag[i][k+1];
						b[i][k]=b[i][k+1];
					}
					J[i]=J[i]-1;
			}
		}
		else //L
		{
			age[i][j-3*J[i]]=age[i][j-3*J[i]]-1;
			J[i-1]=J[i-1]+1;
			for(k=J[i-1]-1;k>0;k=k-1)
			{
				age[i-1][k]=age[i-1][k-1];
				ag[i-1][k]=ag[i-1][k-1];
				b[i-1][k]=b[i-1][k-1];
			}
			ag[i-1][0]=ag[i][j-3*J[i]];
			age[i-1][0]=1.0;
			b[i-1][0]=b[i][j-3*J[i]];
			if(age[i][j-3*J[i]]==0)
				{	
					for(k=(j-3*J[i]);k<J[i];k++)
					{
						age[i][k]=age[i][k+1];
						ag[i][k]=ag[i][k+1];
						b[i][k]=b[i][k+1];
					}
					J[i]=J[i]-1;
			}
		}
}


int jac(){
	return 0;
	}



/*int odeOx (double t,const double y[], double f[], void *params){
	double  k=1.57*1e-4;
	double *data=(double *)params;
	double Omega, densidad;
	Omega=data[0];
	densidad=data[1];
	f[0]=k*Omega-k*densidad*y[0];
	return GSL_SUCCESS;
	}*/

/*int Oxygen(double tini, double tfin,double y[],double datos[] ){
	const gsl_odeiv_step_type * T1 = gsl_odeiv_step_rk2;  //
	gsl_odeiv_step * s1 = gsl_odeiv_step_alloc (T1, 1); //(T y n número de edo-s)
	gsl_odeiv_control * c = gsl_odeiv_control_y_new (1e-6, 0.0);	//error control
	gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (1);
	gsl_odeiv_system sys = {odeOx, jac, 1, &*datos};	//add equation system
	  		
	double h = (tfin-tini)/2.0;			//step

	
	while (tini < tfin){
		int status = gsl_odeiv_evolve_apply (e, c, s1, &sys, &tini, tfin, &h, y);
		if (status != GSL_SUCCESS)
		break;
		}
	//fprintf(stderr,"Despues %f %f \n",tini, y[0]);
	//free memory
	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s1);

	
	return 0;
	}*/

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
/*double birth2(double ox){//para (p6/p3)=0.989
		
		double a_=8250;
		double Ccr=0.0996697;
		if(ox>Ccr)
		{
			return a_*pow(ox/Ccr-1,-0.2);
		}
		else
		{	
			return 1000000000;//9
		}
	}*/



