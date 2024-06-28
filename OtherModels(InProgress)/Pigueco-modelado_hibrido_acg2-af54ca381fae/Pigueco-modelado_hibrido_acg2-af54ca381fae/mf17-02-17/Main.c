
//We changed the IO slightly to account for initial conditions with well mixed age structure
//This additional info is provided by a "Roberto file"

//One spatial dimension plus age structure



#include "Header.h"




int main(int argc, char **argv)
{
    FILE *INITIAL_DATUM, *POPULATION, *PARAMETERS, *OUTPUT_DATA, *MESSAGES, *ROBERTO_FILE;
    
    long age_slots=2; //Number of age slots (we will handle this with malloc-realloc)
    long i_iteration, i_slot,i,j,x_slot,a_slot, n_xslots=0;  //counters and tracers
    
    double t=0.0;
    double delta_t=1.0;
    double delta_x=1.0; //user's value for this quantity, in mm.
    
    //QUANTITIES TO MONITOR AND SYSTEM PARAMETERS
    
    double * *Actual_population; //population vector(s) as a function of age (array of dynamic vectors)
    //This tracks spatial dependence (outer index) and age dependence (inner index)
    //SO FAR the age grid will be taylored to account for an age resolution of 0.5 time units
    //But be very careful about this, if we refine the numerical solver then we must change the overall structure
    
    double death_rate_inv; //death rates^-1
    double tau_p; //parameter related to stepwise birth rates
    double *division_threshold; //threshold parameter at birth rates
    double source_oxygen, k_oxygen; //parameters in the equation for the oxygen
    
    
    double *oxygen_level; //depending on time and space
    double p6_over_p3; //ratio of p's
    double aplus; //switch at zero oxygen, if applicable
    double *number_of_cells; //number of cells at each spatial slot, regardless of their age
    double total_number_of_cells; //overall number of cells in the system
    double peak_number_of_cells;
    
    
	//Stuff to handle output files
    long sampling_period=1; //sampling period
    long number_of_iterations;
    char output_path[70]="OutputValuesPopulation/Out";
    int output_iter=0;

    int MAGIC_NUMBER=10; //To deal with the specific initial condition required for comparison with Roberto
    
    int c; //to sweep and count lines in Roberto file
    double basurilla; //to store Roberto's auxiliary info
    double RobNcells, RobAge;
    double *RobertoNcells, *RobertoAges;
    long sparse_slots,Roberto_lines=0;
   
   
    /************************/
    
    /*****************/
    //OPENING FILES
    
    if(argc!=3){
		printf("Error, incorrect argument number\n");
		printf("Use as: \n exec_name  \n initial_data_file \n parameters_file \n");
		return(1);
	};
    
    if((INITIAL_DATUM=fopen(argv[1],"rt"))==NULL){		
		fprintf(stderr,"Error: initial_data_file could not be opened \n");
		return(1) ;
	};
    
    if((PARAMETERS=fopen(argv[2],"rt"))==NULL){
		fprintf(stderr,"Error: parameters_file could not be opened \n");
		return(1) ;
	};
    
    if((POPULATION=fopen("total_population_vs_time","w"))==NULL){
		fprintf(stderr,"Error: output_population_file could not be opened \n");
		return(1) ;
	};
	

   
    
    /*********************/
    //GATHERING ADDITIONAL INITIAL INFORMATION AND INITIALIZING (I)
    
    n_xslots=Read_Init_Space(INITIAL_DATUM,PARAMETERS,&death_rate_inv,&tau_p,&source_oxygen,&k_oxygen,&oxygen_level,&p6_over_p3,&number_of_cells);
    
    printf("Enter the sampling period\n");
    scanf("%ld",&sampling_period);
    //sampling_period=5000;
    number_of_iterations=N_OUTPUTS*sampling_period; //tune this at the header file
    
    delta_x=1;
    
    
    /**********************/
    //Processing Roberto's stuff...
    
    //Open Roberto file and compute number of lines
    if((ROBERTO_FILE=fopen("edadesJuan.dat","rt"))==NULL){
		printf("Error: Roberto_file could not be opened \n");
		return(1) ;
	};
    
    while((c=getc(ROBERTO_FILE))!=EOF){if(c=='\n'){Roberto_lines++;};};
    Roberto_lines++; //We do not jump from last line (according to Roberto)
    rewind(ROBERTO_FILE);
    
    //Book memory to store Roberto's data
    if((RobertoNcells= (double *) malloc(sizeof(double)))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    if((RobertoAges= (double *) malloc(sizeof(double)))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    //Copy and distill Roberto's data
    fscanf(ROBERTO_FILE,"%lf %lf",&RobertoNcells[0], &RobertoAges[0]);
    RobertoAges[0]=round(2*RobertoAges[0]); //Here and afterwards this is to take into account the halved age grid induced by the splitting method
    j=0;
    
    for(i=1;i<Roberto_lines-1;i++){ //The "-1" should not be correct, but it is...
        fscanf(ROBERTO_FILE,"%lf %lf",&RobNcells, &RobAge);
        RobAge=round(2*RobAge);
        if(fabs(RobAge-RobertoAges[j])<0.5){ //in fact it would be the same if we test equality
            RobertoNcells[j]=RobertoNcells[j]+RobNcells;
        }
        else{
            j++;
            if((RobertoNcells= (double *) realloc(RobertoNcells,sizeof(double)*(j+1)))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            if((RobertoAges= (double *) realloc(RobertoAges,sizeof(double)*(j+1)))==NULL){
                fprintf(stderr,"Error, memory could not be assigned \n");
                exit(1);
            };
            RobertoNcells[j]=RobNcells;
            RobertoAges[j]=RobAge;
        };
    };
    
    fclose(ROBERTO_FILE);
    
    age_slots=RobAge+1; //last age that was read, converted to the halved scale (+1 to account for zero age)
    //j stands now for the non-empty number of age bins (minus one)
    sparse_slots=j;
    //for(i=0;i<=j;i++){printf("%lf %lf \n",RobertoNcells[i],RobertoAges[i]);};
    
    
    
    /**********************/
   //Initialization (array of dynamic vectors) for an age distribution given by
    //the data file provided by Roberto
    //WE WILL NOT EVEN BOTHER TO MOVE THIS TO A SUBROUTINE FOR NOW
    
    //Book the outermost structure
    if((Actual_population= (double **) malloc(sizeof(double *)*n_xslots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
   
    //Copy Roberto's age structure
    if((Actual_population[0]= (double *) malloc(sizeof(double)*age_slots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    j=0;
    for(i=0;i<age_slots;i++){
    
        if(i<RobertoAges[j]){
            Actual_population[0][i]=0.0;
        }
        else{
            Actual_population[0][i]=RobertoNcells[j];
            j++;
        };
    };
    
    //Copy this again to the other spatial slots
    for(x_slot=1;x_slot<MAGIC_NUMBER;x_slot++){
        if((Actual_population[x_slot]= (double *) malloc(sizeof(double)*age_slots))==NULL){
            fprintf(stderr,"Error, memory could not be assigned \n");
            exit(1);
        };
        //Copy from template
        for(i=0;i<age_slots;i++){
            Actual_population[x_slot][i]=Actual_population[0][i];
        };
    };
    
    //Unpopulated slots
    for(x_slot=MAGIC_NUMBER;x_slot<n_xslots;x_slot++){
        if((Actual_population[x_slot]= (double *) malloc(sizeof(double)*age_slots))==NULL){
            fprintf(stderr,"Error, memory could not be assigned \n");
            exit(1);
        };
        Actual_population[x_slot][0]=number_of_cells[x_slot];
        Actual_population[x_slot][1]=0.0;
    };
    
    
    
    /*
     //OLDER CODE, YOU NEVER KNOW...
     for(j=0;j<=sparse_slots;j++){};
     while(i<RobAge[j]){ //Esto funcionara bien si no hay redundancias en el fichero de Roberto (aunque visto lo visto es problematico)
     Actual_population[0][i]=0.0;
     i++;
     };
     Actual_population[0][i]= RobNcells;
     
     
     j=0;
     for(i=0;i<Roberto_lines;i++){
     fscanf(ROBERTO_FILE,"%lf %lf",&RobNcells, &RobAge);
     //printf("%lf %lf\n",RobNcells, RobAge);
     RobAge=round(2*RobAge); //compute closest age in the "0.5-grid" and rank the position in that grid (e.g. an age of 0.8 should go to the third position -labeled with 2, as usual with arrays in C-, after 0,0.5,1,...)
     
     while(j<RobAge){ //Esto funcionara bien si no hay redundancias en el fichero de Roberto (aunque visto lo visto es problematico)
     Actual_population[0][j]=0.0;
     j++;
     };
     Actual_population[0][j]= RobNcells;
     j++;
     };
     
     //Aqui hay un asunto gordo, NO SOLO IMPORTAN EL NUMERO DE CASILLAS SINO TAMBIEN LAS EDADES?
     //Sí, va a haber que hacer un if para decidir que slots se llenan y cuales se dejan vacios; age_slots sera la ultima lectura, muchos slots van a quedar vacios, hay que redondear las edades a enteros -si usamos una malla de enteros para las edades. De hecho no a enteros sino a MULTIPLOS DE 0.5 (provided that delta_t=1 and then we halve the time step when doing a splitting to olve the age transport). If we were to refine delta_t then everything would be much more involved.
     //By the way, this requires an EXTREMELY large amount of memory, roughly 100*2*18800 entries (the structure is quite sparse; for now we will do it by brute force, but this is something to think of in the future).
     */
    
    
  
    /*******************/
    //PRINTING INITIAL CONDITION
    
    //Printing spatial distribution of the number of cells regardless their age
    strcpy(output_path,"OutputValuesPopulation/Out");
    Print_Vector(OUTPUT_DATA,output_path,output_iter,n_xslots,number_of_cells);

    //computing and printing total number of cells
    total_number_of_cells=0.0;
    for(x_slot=0;x_slot<n_xslots;x_slot++){
        total_number_of_cells=total_number_of_cells+number_of_cells[x_slot];
    };
    fprintf(POPULATION,"%lf %lf \n",t,total_number_of_cells);

    //Printing oxygen spatial distribution
    strcpy(output_path,"OutputValuesOxygen/Out");
    Print_Vector(OUTPUT_DATA,output_path,output_iter,n_xslots,oxygen_level);
 
    
    /*******************************/
    //PRINTING SIMULATION INFOS
    
    if((MESSAGES=fopen("messages","w"))==NULL){ //Master output file
		printf("Error: output messages file could not be created \n");
		return(1) ;
	};
    
    fprintf(MESSAGES,"#Debug info: \n \n");
    fprintf(MESSAGES,"#Data generated with the following .exe file:%s\n",argv[0]);
    fprintf(MESSAGES,"#Initial condition file:%s\n",argv[1]);
    fprintf(MESSAGES,"#Parameters file:%s\n",argv[2]);
    fprintf(MESSAGES,"#Event number: %ld \n \n",number_of_iterations);
    fprintf(MESSAGES,"#Specific infos: oxygen damping, finite differences, a single cellular line\n");
    
    fprintf(MESSAGES,"#Death rate: %lf \n",1.0/death_rate_inv);
    fprintf(MESSAGES,"#tau_p: %lf\n",tau_p);
    fprintf(MESSAGES,"#source_oxygen: %lf\n",source_oxygen);
    fprintf(MESSAGES,"#k_oxygen: %lf\n",k_oxygen);
    fprintf(MESSAGES,"#p6/p3: %lf\n",p6_over_p3);
    fprintf(MESSAGES,"#slot size: %lf mm\n \n \n",delta_x);
     
    fprintf(MESSAGES,"#List of issues:\n");
    //So far "issues" are printed in stderr, to avoid passing *MESSAGES to subroutines
   
    
    
    /*******************************/
    //INITIALIZING DIVISION RATES
   
    //We need to compute the division at zero oxygen for the sake of formulas
    //book division_threshold
    if((division_threshold= (double *) malloc(sizeof(double)*n_xslots))==NULL){
        fprintf(stderr,"Error, memory could not be assigned \n");
        exit(1);
    };
    
    Compute_division_threshold(division_threshold,oxygen_level,p6_over_p3,n_xslots-1);
    
    /******************************/
    
    
    /******************************/
    //MAIN LOOP
    
     delta_t=1;
    //The OVERALL age structure is adapted to the fact that age transport solvers will advance delta_t/2 at a time (this is our age quanta), as the code is now designed.
    //Be EXTREMELY CAREFUL with this
    
    for(i_iteration=1;i_iteration<=number_of_iterations;i_iteration++){
        
        //Book additional memory for the population vector if needed
        //UNTIL we find a more clever structure, if we are to enlarge one row then we must enlarge all (to be able to accomodate spatial diffusion of old cells into a spatial slot populated by younger ones)
        Handle_memory_extension(Actual_population,&age_slots,n_xslots);
        
        
        //ADVANCE THE POPULATION
        //Calling the solver (boundary conditions are determined in this module)
        //A particular SPLITTING PROPOSAL: transport first, diffusion next
        
        for(x_slot=0;x_slot<n_xslots;x_slot++){ //loop in space
            FakeTransportSolver(delta_t/2,Actual_population[x_slot],age_slots,death_rate_inv,division_threshold[x_slot],tau_p);
        };
     //PERO ESTO NO ES UN BUG, LO DEL 1/2, QUE DeSCUADRA LA RIGIDEZ DEL TRANSPORT SOLVER???
     //Al solver no le afecta, pero la interpretacion del mallado en edad cambia -y por tanto los division threshold se han de recalcular de forma acorde, CAGADA? No, porque la forma de manejar los division thresholds en el modulo de transporte tiene en cuenta que el delta_t no necesariamente ha de ser uno. PERO, el mallado en edad ha de estar ajustado globalmente a delta_t/2 cuando leemos de Roberto (que es lo que se pasa al transport_solver) o la cagamos.
        
        for(a_slot=0;a_slot<age_slots;a_slot++){ //loop in age
            DiffusionSolver(delta_t/2,Actual_population,a_slot,n_xslots,delta_x);
        };
       
        //repeat for a total (virtual) advancement of 2*delta_t
        //to be honest it seems fair to recompute division_threshold in between
        //but we will disregard it due to time computing issues...
        for(x_slot=0;x_slot<n_xslots;x_slot++){ //loop in space
            FakeTransportSolver(delta_t/2,Actual_population[x_slot],age_slots,death_rate_inv,division_threshold[x_slot],tau_p);
        };
        for(a_slot=0;a_slot<age_slots;a_slot++){ //loop in age
            DiffusionSolver(delta_t/2,Actual_population,a_slot,n_xslots,delta_x);
        };
        //Afterwards we take care of this and regard this whole round as an advancement over just delta_t.
        /************/
        
        
        //ADVANCE THE OXYGEN
        
        //compute number_of_cells at each spatial slot
    
        peak_number_of_cells=Cummulative_Cells(number_of_cells, age_slots, Actual_population,n_xslots);
        
        OxygenSolverDamp(oxygen_level,n_xslots,number_of_cells,delta_t,k_oxygen,source_oxygen,peak_number_of_cells,delta_x);
       
        
        
        /********************/
        //RECOMPUTE SETTINGS FOR NEXT ITERATION
    
        t=t+delta_t;
        
        //recompute division_threshold
        Compute_division_threshold(division_threshold,oxygen_level,p6_over_p3,n_xslots-1);
               
        /******************/
       //PRINTING DATA ON SELECTED ITERATIONS
        
        if((i_iteration%sampling_period)==0){
            
            output_iter=i_iteration/sampling_period;
            
            //printing total number of cells
            total_number_of_cells=0.0;
            for(x_slot=0;x_slot<n_xslots;x_slot++){
                total_number_of_cells=total_number_of_cells+number_of_cells[x_slot];
            };
            fprintf(POPULATION,"%lf %lf \n",t,total_number_of_cells);
            
            //Printing spatial distribution of cells
            strcpy(output_path,"OutputValuesPopulation/Out");
            Print_Vector(OUTPUT_DATA,output_path,output_iter,n_xslots,number_of_cells);
            
            //Printing oxigen spatial distribution
            strcpy(output_path,"OutputValuesOxygen/Out");
            Print_Vector(OUTPUT_DATA,output_path,output_iter,n_xslots,oxygen_level);
            
        
        }; //end if printing outputs
        
    }; //end main loop

    
    /***************************************/
    

     
    /***************************************/
    
    
    //CLOSING FILES AND RETRIEVING MEMORY
    
    for(x_slot=0;x_slot<n_xslots;x_slot++){
        free(Actual_population[x_slot]);
    };
    free(Actual_population);
    free(number_of_cells);
    free(oxygen_level);
    free(division_threshold);
    
    fclose(POPULATION);
    fclose(MESSAGES);
    
	return(0);  
	
}      //end main




