#include "Header.h"


double reactions(long first, long end, double **gilAux,  long *J, double **age, double **b, double nu, double delta_x) {
		//supongo que first es el primer índice (sería 0 si consideramos todo el vector)

		long i,j;
		double suma=0.0;
        double R,L;
    
    R=DIFFUSION_COEF_CELLS/(delta_x*delta_x);
    L=R;
		
		for(j=0;j<J[first];j++){
				gilAux[first][j]=age[first][j]*nu;
				gilAux[first][j+J[first]]=age[first][j]*b[first][j];
				gilAux[first][j+2*J[first]]=age[first][j]*R;
				gilAux[first][j+3*J[first]]=0.0;
				suma=suma+age[first][j]*(nu+b[first][j]+R);
		}
    /*for(j=0;j<J[first+1];j++){
        gilAux[first+1][j]=age[first+1][j]*nu;
        gilAux[first+1][j+J[first+1]]=age[first+1][j]*b[first+1][j];
        gilAux[first+1][j+2*J[first+1]]=age[first+1][j]*R;
        gilAux[first+1][j+3*J[first+1]]=0.0;
        suma=suma+age[first+1][j]*(nu+b[first+1][j]+R);
    }*/
						
		for(i=(first+1);i<(end-1);i++)
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
//----------------------------------------------------------------------------------

//double  //Now we compute tau in the main body, then we print and only then we update
void gillespie(double tau, 
                long first,    
                long end, 
                double **gilAux, 
                double ** ageAux,
                double **age, 
                long *J, 
                double **b, 
                double **ag, 
                double suma, 
                long *pseed, 
                double *number_of_cells, 
                double *division_threshold, 
                double *oxygen_level, 
                double p6_over_p3, 
                double tau_p, 
                long *ptam
                ){
    
	long i,j,k;
	double r2, rtot;
	double auxsuma=0.0;
	
	//r1=ran2(pseed);//drand48();
	r2=ran2(pseed);//drand48();
	//tau=log(1/r1)/suma;
	rtot=suma*r2;
	
	i=first-1;
	while(auxsuma<=rtot)
	{
		i=i+1;
		for(j=0;j<4*J[i];j++)
        {
            auxsuma=auxsuma+gilAux[i][j];
            if(auxsuma>rtot)
            {break;}
        }
	}
	if(j<J[i]) //Muerte
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
        number_of_cells[i]=number_of_cells[i]-1;
	}
	else if(j<2*J[i]) //Nacimiento
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
        number_of_cells[i]=number_of_cells[i]+1;
        
        if(J[i]>=*ptam){
            *ptam=*ptam+TAM;
            
            if((b[i]= (double *) realloc(b[i],sizeof(double)*(*ptam)))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            
            if((age[i]= (double *) realloc(age[i],sizeof(double)*(*ptam)))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            
            if((ag[i]= (double *) realloc(ag[i],sizeof(double)*(*ptam)))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            
            if((gilAux[i]= (double *) realloc(gilAux[i],sizeof(double)*(*ptam)*4))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            
            if((ageAux[i]= (double *) realloc(ageAux[i],sizeof(double)*(*ptam)))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            
        };//end if
        
        
	}
	else if(j<3*J[i]) //Difusion derecha
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
        number_of_cells[i+1]=number_of_cells[i+1]+1;
        number_of_cells[i]=number_of_cells[i]-1;
        
        if(J[i+1]>=*ptam){
            *ptam=*ptam+TAM;
            
            if((b[i+1]= (double *) realloc(b[i+1],sizeof(double)*(*ptam)))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            
            if((age[i+1]= (double *) realloc(age[i+1],sizeof(double)*(*ptam)))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            
            if((ag[i+1]= (double *) realloc(ag[i+1],sizeof(double)*(*ptam)))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            
            if((gilAux[i+1]= (double *) realloc(gilAux[i+1],sizeof(double)*(*ptam)*4))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            
            if((ageAux[i+1]= (double *) realloc(ageAux[i+1],sizeof(double)*(*ptam)))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            
        };//end if
	}
	else //difusion izquierda
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
        number_of_cells[i-1]=number_of_cells[i-1]+1;
        number_of_cells[i]=number_of_cells[i]-1;
        
        if(J[i-1]>=*ptam){
            *ptam=*ptam+TAM;
            
            if((b[i-1]= (double *) realloc(b[i-1],sizeof(double)*(*ptam)))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            
            if((age[i-1]= (double *) realloc(age[i-1],sizeof(double)*(*ptam)))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            
            if((ag[i-1]= (double *) realloc(ag[i-1],sizeof(double)*(*ptam)))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            
            if((gilAux[i-1]= (double *) realloc(gilAux[i-1],sizeof(double)*(*ptam)*4))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            
            if((ageAux[i-1]= (double *) realloc(ageAux[i-1],sizeof(double)*(*ptam)))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            
        };//end if
        
	}
    
    for(i=first;i<end;i++)
	{
		for(j=0;j<J[i];j++)
		{
			if(oxygen_level[i]>c_cr(p6_over_p3))// este es el valor límite para la relación p6/p3=1.0
			{
				ag[i][j]=ag[i][j]+tau;
			}
			if(ag[i][j]>=division_threshold[i])
			{
				b[i][j]=1/tau_p;
			}
		}
	}
    
	//return tau;

}
