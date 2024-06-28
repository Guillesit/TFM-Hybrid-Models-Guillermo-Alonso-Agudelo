

//We don't include here a damping term at oxigen's equation so far

//Template for hybrid code, a single cellular line
//One spatial dimension

//As far as the code is developed, initial conditions are supposed to flow from left to right across the spatial domain

//The interface is now composed of two slots. Interface_location traces the leftmost

#include "Header.h"




int main(int argc, char **argv)
{
    FILE *INITIAL_DATUM, *POPULATION, *PARAMETERS, *OUTPUT_DATA, *MESSAGES, *ROBERTO_FILE;
    
    long i_iteration,i,j,k,x_slot, n_xslots=0;  //counters and tracers
    
    double t=0.0, tau, tstop; //current time, Gillespie's time step (measured in funny adimensional time units)
    double delta_x=1.0; //user's value for this quantity, in mm.
    long seed; //To handle ran2()
    double r1; //first random number, associated with tau
    
    //QUANTITIES TO MONITOR AND SYSTEM PARAMETERS
    double death_rate_inv; //death rates^-1
    double tau_p; //parameter related to stepwise birth rates
    double *division_threshold; //threshold parameter at birth rates
    double source_oxygen, k_oxygen; //parameters in the equation for the oxygen
    
    double *oxygen_level; //depending on time and space
    double p6_over_p3; //ratio of p's
    //double aplus; //switch at zero oxigen, if applicable
    double *number_of_cells,*density_of_cells; //number of cells at each spatial slot, regardless of their age (resp. density of cells at each spatial slot)
    //The mean field equation for the population uses densities, the global equation for the oxygen uses number_of_cells
    
    double total_number_of_cells; //overall number of cells in the system
    long Interface_location, aux_interface; //To handle interface location and displacements
        //Interface_location traces the leftmost slot in the interface compartment. Note that the rightmost slot may well be below threshold
    
	//Stuff to handle output files
	long n_files=1; //desired number of output files (typically 100)
    //long sampling_period=1; //sampling period
    //long number_of_iterations; //100*sampling_period as designed
    char output_path1[200]="OutputValuesPopulation"; //NUEVO (17/09/2016)
    char output_path2[200]="OutputValuesOxygen"; //NUEVO (17/09/2016)
    char usertag[100]=""; //NUEVO (16/09/2016)
    char label1[200]="OutputValuesPopulation"; //NUEVO (17/09/2016)
    char label2[200]="OutputValuesOxygen"; //NUEVO (17/09/2016)
    char label3[200]="messages"; //NUEVO (17/09/2016)
    char label4[200]="total_population_vs_time"; //NUEVO (17/09/2016)
    int counter=1; //MAY cause trouble if we want a large number of output files
    long max_n_iter=0; //max number of events allowed
    double sampling_time_window=0.0;
    char create_folder[300]="mkdir "; //NUEVO (17/09/2016)
    
    
    double Diff_Coef_Pop,Diff_Coef_O2;
    
    
    //Estructuras de datos de la parte de Roberto:
    double **age,**b,**ag;
    long *J;
    long tam=TAM; //500000; //arreglar este manejo en cuanto el codigo funcione
    //Tener mucho cuidado porque aqui no hay garbage collector!!!!!
    double **ageAux;
    double **gilAux;
    double suma;
    
    //Variables para procesar la distribucion de edad en equilibrio que nos pasa Roberto
    int MAGIC_NUMBER=10; //To deal with the specific boundary condition required for comparison with Roberto
    int c; //to sweep and count lines in Roberto file
    double basurilla; //to store Roberto's auxiliary info
    double RobNcells, RobAge;
    double *RobertoNcells, *RobertoAges;
    long sparse_slots,Roberto_lines=0;
    long age_slots=2; //Number of nonempty age slots at equilibrium...
    
    
    double newAux; //, newAux2; //variables necesarias para corregir la renormalizacion. Corresponde al número 
    // de partículas del estócastico. Hay que pasarlas al Renormalize_ porque las necesita (ver codigo renormalize)
    /************************/
    
    /*****************/
    //OPENING FILES
    
    if(argc!=4){
		printf("Error, incorrect argument number\n");
		printf("Use as: \n exec_name  \n initial_data_file \n parameters_file \n output_files_tag \n");
		return(1);
	};
    
    if((INITIAL_DATUM=fopen(argv[1],"rt"))==NULL){
		printf("Error: initial_data_file could not be opened \n");
		return(1) ;
	};
    
    if((PARAMETERS=fopen(argv[2],"rt"))==NULL){
		printf("Error: parameters_file could not be opened \n");
		return(1) ;
	};
    
    strcpy(usertag,argv[3]);
    
    //create "total_population_vs_time" (plus user tag)
    strcat(label4,usertag);
    if((POPULATION=fopen(label4,"w"))==NULL){
		printf("Error: output_population_file could not be opened \n");
		return(1) ;
	};

   
    
    /*********************/
    //GATHERING ADDITIONAL INITIAL INFORMATION AND INITIALIZING
    
    n_xslots=Read_Init_Space(INITIAL_DATUM,PARAMETERS,&death_rate_inv,&tau_p,&source_oxygen,&k_oxygen,&oxygen_level,&p6_over_p3,&number_of_cells);
    //We assume that the numbers in the initial condition file represent number of cells per slot, not densities
    
    //ask for sampling period
    printf("Enter stopping time (adimensional units...)\n");
    scanf( "%lf",&tstop);
    printf("Enter the desired number of output files\n");
    scanf("%ld",&n_files);
    sampling_time_window=tstop/n_files;   
    //printf("%lf\n",sampling_time_window);
    
    //ask for slot size
    //Once this is finely checked we should provide this in the parameters file
    //printf("Enter slot size (in mm)\n");
    //scanf("%lf",&delta_x);
    delta_x=1;
    
    Diff_Coef_Pop=DIFFUSION_COEF_CELLS/(delta_x*delta_x);
    Diff_Coef_O2=DIFFUSION_COEF_O2/(delta_x*delta_x);
    
    
    //COMPUTE SYSTEM's initial INTERFACE
    Interface_location=0;
    while((number_of_cells[Interface_location]>=POPULATION_THRESHOLD)&&(Interface_location<n_xslots-2)){
        Interface_location++;
    }; //Vector index, nor the natural one
    
    Interface_location=Interface_location-1; //Sets the leftmost interface slot at the last slot above threshold
    //book density_of_cells
    if((density_of_cells= (double *) malloc(sizeof(double)*n_xslots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    
    //Initializing random number generator
    //NOTE THAT ROBERTO WILL BE USING A DIFFERENT ONE
    srand(time(NULL));
    seed= -rand();
    
    
    /********************************/
    //Handling Roberto's equilibrium age structure file (part I...)
    
    //Open Roberto file and compute number of lines
    if((ROBERTO_FILE=fopen("edadesJuan.dat","rt"))==NULL){
		printf("Error: Roberto_file could not be opened \n");
		return(1) ;
	};
    while((c=getc(ROBERTO_FILE))!=EOF){if(c=='\n'){Roberto_lines++;};};
    //Roberto_lines++; //We do not jump from last line (according to Roberto)
    //Esto ha demostrado ser falso
    rewind(ROBERTO_FILE);
    
    //Copy Roberto's data (we assume that there is no redundant data)
    if((RobertoNcells= (double *) malloc(sizeof(double)*Roberto_lines))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    if((RobertoAges= (double *) malloc(sizeof(double)*Roberto_lines))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    for(i=0;i<Roberto_lines;i++){
        fscanf(ROBERTO_FILE,"%lf %lf",RobertoNcells+i, RobertoAges+i);//do not trust this so much, do some trials
    };
    
    fclose(ROBERTO_FILE);
    //Apparently after that we have removed redundancies (...)
    //ASK TO HAVE A FILE WITHOUT REDUNDANCIES
    //printf("%lf %lf \n", RobertoNcells[Roberto_lines-1],RobertoAges[Roberto_lines-1]);
   
   
    
    
    
    /*******************************/
    //INITIALIZING DIVISION RATES
    
    //We need to compute the division rate at zero oxygen for the sake of formulas
    //Just call the ODE solver once
    //As it stands, it works only for a single cellular line
    
    //book division_threshold
    if((division_threshold= (double *) malloc(sizeof(double)*n_xslots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    //Proceed to initialize the stuff
    Compute_division_threshold(division_threshold,oxygen_level,p6_over_p3,n_xslots-1);
    

    
    /**********************/
    //INITIALIZING ROBERTO'S STUFF...
    //n_xslots plays the role of box
    
    //Handling Roberto's equilibrium age structure file (part II) while we are at it
    //Copy Roberto's age structure into Roberto's funny arrays
    //Just for slots in the interface, given the very special initial condition
    // we are dealing with. Well, we changed our minds, see below
    ag=(double **)malloc(n_xslots*sizeof(double *));
	for(i=0;i<n_xslots;i++){ag[i]= (double*)calloc(TAM,sizeof(double));};
    
    J=(long *)calloc(n_xslots,sizeof(long));
	for (i=0;i<n_xslots;i++){
        if(number_of_cells[i]>0){J[i]=Roberto_lines;/*1;*/}
        //pienso ahora que vamos a meter la estructura computacional para tener edades en equilibrio en toda parte poblada, 
        //no solo en la interfaz, total poco nos cuesta y puede prevenir (o crear?) side effects
        else{J[i]=0;};
    };
    
	age=(double **)malloc(n_xslots*sizeof(double *));
	
	for(i=0;i<n_xslots;i++){age[i]= (double*)calloc(TAM,sizeof(double));};
    
	b=(double **)malloc(n_xslots*sizeof(double *)); 	//birth rate
	for(i=0;i<n_xslots;i++){b[i]= (double*)calloc(TAM,sizeof(double));};
    
    for (i=0;i<n_xslots;i++){
		for(j=0;j<J[i];j++){
			age[i][j]=RobertoNcells[j]; //I do not trust the +1, we slip out of range
            //(recall that Roberto_lines-1 is the last meaningful index)
            ag[i][j]=RobertoAges[j];
            if(ag[i][j]>=division_threshold[i]){b[i][j]=1/tau_p;};
		};
	};
    
    ageAux=(double **)malloc(n_xslots*sizeof(double *));
	for(i=0;i<n_xslots;i++){ageAux[i]= (double*)calloc(TAM,sizeof(double));};
    
    gilAux=(double **)malloc(n_xslots*sizeof(double *));
	for(i=0;i<n_xslots;i++){gilAux[i]= (double*)calloc(4*TAM,sizeof(double));};
   //Se inicializa en el bucle del Gillespie
   
    
    
    
    /*******************/
    //PRINTING INITIAL CONDITION
    
    //computing and printing total number of cells
    //initializing cell density while we are at it
    total_number_of_cells=0.0;
    for(x_slot=0;x_slot<n_xslots;x_slot++){
        density_of_cells[x_slot]=number_of_cells[x_slot]/delta_x;
        total_number_of_cells=total_number_of_cells+number_of_cells[x_slot];
    };
    fprintf(POPULATION,"%lf %lf \n",t,total_number_of_cells);
    
    //create "OutputValuesPopulation" and "OutputValuesOxygen" (plus user tag)
    strcat(label1,usertag);
    strcat(label2,usertag);
    
    strcat(create_folder,label1);
    system(create_folder);
    strcpy(create_folder,"mkdir ");
    strcat(create_folder,label2);
    system(create_folder);
    
    
    //Printing spatial DENSITY of cells regardless of their age
    strcat(label1,"/Out"); //Final label for population folder
    strcpy(output_path1,label1);
    //strcpy(output_path1,"OutputValuesPopulation/Out");
    Print_Vector(OUTPUT_DATA,output_path1,0,n_xslots,density_of_cells);

    //Printing oxygen spatial distribution
    strcat(label2,"/Out"); //Final label for oxygen folder
    strcpy(output_path2,label2);
    //strcpy(output_path,"OutputValuesOxygen/Out");
    Print_Vector(OUTPUT_DATA,output_path2,0,n_xslots,oxygen_level);
    
    
    
    
    /*******************************/
    //PRINTING SIMULATION INFOS
    
    
    //create "messages" (plus user tag)
    strcat(label3,usertag);
    if((MESSAGES=fopen(label3,"w"))==NULL){ //Master output file
		printf("Error: output messages file could not be created \n");
		return(1) ;
	};
    
    fprintf(MESSAGES,"#Debug info: \n \n");
    fprintf(MESSAGES,"#Data generated with the following .exe file:%s\n",argv[0]);
    fprintf(MESSAGES,"#Initial condition file:%s\n",argv[1]);
    fprintf(MESSAGES,"#Parameters file:%s\n",argv[2]);
    fprintf(MESSAGES,"#Stopping time: %lf \n \n",tstop); 
    fprintf(MESSAGES,"#Specific infos: no oxygen damping, eigenvalue shortcut\n");
    fprintf(MESSAGES,"#Specific infos: hybrid method, a single cellular line\n");
    
    fprintf(MESSAGES,"#Death rate: %lf \n",1.0/death_rate_inv);
    fprintf(MESSAGES,"#tau_p: %lf\n",tau_p);
    fprintf(MESSAGES,"#source_oxygen: %lf\n",source_oxygen);
    fprintf(MESSAGES,"#k_oxygen: %lf\n",k_oxygen);
    fprintf(MESSAGES,"#p6/p3: %lf\n",p6_over_p3);
    fprintf(MESSAGES,"#slot size: %lf mm \n",delta_x);
    fprintf(MESSAGES,"#Diff coef cells: %.15lf (user units) \n",Diff_Coef_Pop);
    fprintf(MESSAGES,"#Diff coef oxygen: %.15lf (user units)\n \n \n",Diff_Coef_O2);

    fprintf(MESSAGES,"#threshold value: %lf\n",(double) POPULATION_THRESHOLD);
    fprintf(MESSAGES,"#initial memory size: %ld\n",(long) TAM);
     
    fprintf(MESSAGES,"#List of issues:\n");
    //So far "issues" are printed in stderr, to avoid passing *MESSAGES to subroutines
    //Think about this...
    
    
     
    
    /******************************/
    
    
    /******************************/
    //MAIN LOOP
    
    while(t<tstop){
        
        
        //DEAL WITH THE STOCHASTIC PART
        suma=reactions(Interface_location,n_xslots,gilAux,J,age,b,death_rate_inv,delta_x);

        r1=ran2(&seed);//drand48();
		tau=log(1/r1)/suma;
        //printf("tau=%lf\n",tau);
        
        

       	  /******************/
       //PRINTING DATA ON SELECTED ITERATIONS
        //again the issue numbers vs densities

       	if(t+tau>=counter*sampling_time_window){ //print the last state prior surpassing the time mark

       		//printf("Interface set at %ld\n",Interface_location);
            
            //printing total number of cells
            //We take the opportunity to compute densities even in the stochastic part,
            //to use this for the output
            total_number_of_cells=0.0;
            for(x_slot=0;x_slot<n_xslots;x_slot++){
                density_of_cells[x_slot]=number_of_cells[x_slot]/delta_x;
                total_number_of_cells=total_number_of_cells+number_of_cells[x_slot];
            };
            fprintf(POPULATION,"%lf %lf \n",t,total_number_of_cells);
            
            //Printing spatial density of cells
            strcpy(output_path1,label1);
            //strcat(output_path,"/Out");
            //strcpy(output_path,"OutputValuesPopulation/Out");
            Print_Vector(OUTPUT_DATA,output_path1,counter,n_xslots,density_of_cells);
            
            //Printing oxigen spatial distribution
            strcpy(output_path2,label2);
            //strcat(output_path,"/Out");
            //strcpy(output_path,"OutputValuesOxygen/Out");
            Print_Vector(OUTPUT_DATA,output_path2,counter,n_xslots,oxygen_level);

            printf("t=%lf\n",counter*sampling_time_window);
        	
        	counter++;
        };
        //print stuff, update right after

        
        //NOW advance intracellular part from t to t+tau
        gillespie(tau,Interface_location,n_xslots,gilAux,ageAux,age,J,b,ag,suma,&seed,number_of_cells,division_threshold,oxygen_level,p6_over_p3,tau_p,&tam);
       	// number_of_cells is updated inside this module


        
        //COMPUTE "FLUX" FROM STOCHASTIC TO DETERMINISTIC REGION
        density_of_cells[Interface_location]=number_of_cells[Interface_location]/delta_x;
        //density_of_cells[Interface_location+1]=number_of_cells[Interface_location+1]/delta_x; Fuera porque solo queremos un hueco de la interfaz
        
        
        //Add 14/11/2016
        newAux=density_of_cells[Interface_location];
        //newAux2=density_of_cells[Interface_location+1];
        //fprintf(stderr,"t-%e \n",number_of_cells[Interface_location]);
        /*******************/
        //EVOLVE MEAN FIELD REGION
        
        //Determine suitable refinement of the time step tau before we loop through population-oxygen as many times as needed (splitting 1/2 1/2 once the time step is short enough)
        
        Global_Mean_Field_Handler(tau,number_of_cells,density_of_cells,division_threshold,oxygen_level,Diff_Coef_Pop,Diff_Coef_O2,k_oxygen,source_oxygen, delta_x,tau_p,death_rate_inv,p6_over_p3,Interface_location,n_xslots);
        //WE SHOULD USE THE DAMPED VERSIN HERE (11/11/16)
       
        
        /************/
        
        //fprintf(stderr,"u-%e \n",number_of_cells[Interface_location]);
		//1/0;
        //COMPUTE FLUX FROM  DETERMINISTIC To STOCHASTIC REGION
        
        //Queremos poner solo un espacio en la interfase, basta con quitar sólo esta instruccion?: (y todo lo demas en interface_location+1)
		//Renormalize_center(number_of_cells,Interface_location+1,&seed,division_threshold[Interface_location+1], death_rate_inv, tau_p,age,J,b,ag,gilAux,ageAux,&tam,newAux2);
        
        Renormalize_center(number_of_cells,Interface_location,&seed,division_threshold[Interface_location], death_rate_inv, tau_p,age,J,b,ag,gilAux,ageAux,&tam,newAux);
        
        //Now we recompute again cell densities (as masses may have been renormalized)
        //and when we sweep we DO include the interface (two compartments now)
        //for(x_slot=0;x_slot<=Interface_location+1;x_slot++){
		for(x_slot=0;x_slot<=Interface_location;x_slot++){	
            density_of_cells[x_slot]=number_of_cells[x_slot]/delta_x;
        };
        
        
        /********************/
        //RECOMPUTE SETTINGS FOR NEXT ITERATION
    
        t=t+tau;
        
       //RECOMPUTE SYSTEM INTERFACE
        //Determine where should it be
        aux_interface=Realocate_interface(number_of_cells,n_xslots);
        
        //Renormalization procedure if the interface moves left
        if(aux_interface<Interface_location){
            Renormalize_left(Interface_location,aux_interface,number_of_cells,&seed,division_threshold,death_rate_inv,tau_p,age,J,b,ag,gilAux,ageAux,&tam,newAux);
            printf("We shifted left at time t=%lf\n",t);
            fprintf(MESSAGES,"We shifted left at time t=%lf\n",t);
        };
        
        //Now we recompute again cell densities (as masses may have been renormalized)
        for(x_slot=0;x_slot<=aux_interface;x_slot++){
            density_of_cells[x_slot]=number_of_cells[x_slot]/delta_x;
        };
        
        
        //Renormalization procedure if the interface moves right
        if(aux_interface>Interface_location){
           // for(x_slot=Interface_location+1;x_slot<=aux_interface+1;x_slot++){
			  for(x_slot=Interface_location+1;x_slot<=aux_interface;x_slot++){
                density_of_cells[x_slot]=number_of_cells[x_slot]/delta_x;
            };
            printf("We shifted right at time t=%lf\n",t);
            fprintf(MESSAGES,"We shifted right at time t=%lf\n",t);
        };
        
        
        //Actualize interface tag for next iteration
        Interface_location=aux_interface;
        
        //recompute division_threshold
        Compute_division_threshold(division_threshold,oxygen_level,p6_over_p3,n_xslots-1);
      
        
    }; //end main loop

    
    /***************************************/
     
    /***************************************/
    
    
    //CLOSING FILES AND RETRIEVING MEMORY
    
    free(number_of_cells);
    free(density_of_cells);
    free(oxygen_level);
    free(division_threshold);
    
    fclose(POPULATION);
    fclose(MESSAGES);
    
    //Roberto's stuff:
    
    free(J);
	//free(ag1s);
	
	for(k=0;k<n_xslots;k++){
		free(age[k]);
		free(ageAux[k]);
		free(b[k]);
		free(ag[k]);
		free(gilAux[k]);
	}
	free(b);
	free(age);
	free(ageAux);
    free(gilAux);
    free(ag);
    
    free(RobertoNcells);
    free(RobertoAges);

    /***************************************/
    //Probably this
    //won't be used anymore (17/09/2016)
    //if(user_flag1){
    //command="mv OutputValuesPopulation OutputValuesPopulation";
    //strcat(command,usertag);
    //system(command);
    //system(mkdir OutputValuesPopulation);
    //command="mv OutputValuesOxygen OutputValuesOxygen";
    //strcat(command,usertag);
    //system(command);
    //system(mkdir OutputValuesOxygen);
    //command="mv messages messages";
    //strcat(command,usertag);
    //system(command);
    //command="mv messages messages";
    //strcat(command,usertag);
    //system(command);
    //command="mv total_population_vs_time total_population_vs_time";
    //strcat(command,usertag);
    //system(command);
	//}
    
    
	return(0);  
	
}      //end main




