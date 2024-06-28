#########################################################################PARALELL################################################
#use multiprocessing to do this computation in parallel in multiple processe
# Additionally the script is able to create a pool of processes so different calls
# for different input files can be run in parallel, one in each processor
from multiprocessing import Process,Pool,Lock  #parallel processing
import multiprocessing as mp
#########################################################################PARALELL################

import os

try:
    file = open('tiempos', 'r')
except IOError:
    file = open('tiempos', 'w')

# tag for the differents realisations
Number_simulations=range(1,50)
def hybrid_ejecutable(i):
    # outputdir="OutputValuesHostPopulation_2pop_p3p61d001_therapy0d6_100_%s"%i
    # if not os.path.exists(outputdir):
    print i
	#executable and files used for it
    os.system("./EXE host100derechaEdad invader100izdaEdad Init_oxy100 Param_host Param_host Parsim _2pop_p3p6equal1_%s >> tiempos &"%i)
    #os.system("./EXE host100derechaEdad invader100izdaEdad Init_oxy100 Param_host Param_invader Parsim _2pop_p3p61d001_%s >> tiempos"%i)
    # os.system("./EXE host100derechaEdad invader100izdaEdad Init_oxy100 Param_host Param_invader Parsim_0d6 _2pop_p3p61d001_therapy0d6_100_%s >> tiempos &"%i)
	

# # # ########  RUN IN PARALLEL #############################################
if __name__ == '__main__' or True: # This is necessary... (don't ask me why). Note: Actually it wasn't as the True option states
    cpunum=mp.cpu_count()
    pool = Pool(processes=cpunum) # creating a pool with processors equal to the number of processors
    print("Number of processors "+str(cpunum))
    # print 'Start at %s' %datetime.datetime.now()
    pool.map(hybrid_ejecutable,Number_simulations)  
 	# Obs: change the parameter pair choosen inside the function 
    pool.close()
    pool.join()
    # print 'All subprocesses done at %s' %datetime.datetime.now()
#####  RUN IN PARALLEL #############################################