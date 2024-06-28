
# coding: utf-8


# Setting the terminal postscript with the options
import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import math as m

t_end = 2000011
n_time_save = 100
Number_simulation = 100
threshold = 2000.0

###### take files data name, every file is a snashop in time
foldername='DatosSimulacion/1/OutputValuesPopulation/'
name=[]
for i in range(0,10):
    for file in os.listdir(foldername):
        if file[:5]=="Out%s."%i:
#             print file
            name.append(file)
for i in range(10,n_time_save):
    for file in os.listdir(foldername):
        if file[:6]=="Out%s."%i:
#             print file
            name.append(file)
for file in os.listdir(foldername):
        if file[:7]=="Out%s."%n_time_save:
#             print file
            name.append(file)
 



###### Stochastich files information
# column[0] = space position (box)
# column[1] = birth events in this box
# column[2] = population (number of cells) in box


Pop_stoc={}
for file in name:
    j=1
    foldername='DatosSimulacion/%s/OutputValuesPopulation/'%j
    aux=[]
    crs = open(foldername+file, "r")
    for columns in ( raw.strip().split() for raw in crs ):
        aux.append(columns[2])
    aux=np.array(map(float,aux))
    Pop_stoc[file]=aux
    
    for j in range(1,(Number_simulation+1)):
        foldername='DatosSimulacion/%s/OutputValuesPopulation/'%j
        aux=[]
        
        crs = open(foldername+file, "r")
        for columns in ( raw.strip().split() for raw in crs ):
            aux.append(columns[2])
        aux=np.array(map(float,aux))
        Pop_stoc[file]=aux+Pop_stoc[file]
    Pop_stoc[file]=Pop_stoc[file]/Number_simulation # divide by number of simulations


# # Empirical birth rate

# Values used in the simulations
# tp= 480 
# nu_inv = 0.000041667

"""Save create the folder to save the data"""
outputdir="images"
if not os.path.exists(outputdir): # if the folder doesn't exist create it
    os.makedirs(outputdir)


plt.figure()
Pop_staux_mean=[]
for l in name:
    Pop_staux=np.array(Pop_stoc[l])
    Pop_staux_mean.append(np.mean(np.where((Pop_staux<(Pop_staux[0]-200)) & (Pop_staux>0.01))[0]))
plt.plot(Pop_staux_mean,'-',linewidth=3,label="$n$")
plt.xticks(range(0,100,100/5),range(0,2000000/n_time_save,2000000/(5*n_time_save)), fontsize = 15) # work on current fig
plt.legend()
plt.yticks( fontsize = 15) # work on current fig
plt.xticks( fontsize = 15)
plt.xlabel('time', fontsize=15)
plt.ylabel(r'Mean X axis TW pop', fontsize=15)
plt.savefig(outputdir+'/mean_position_TW_h_cg2.eps', dpi=100, bbox_inches='tight')
    

plt.figure()
velocity_st=[]
t_step = t_end/n_time_save
aux=1
for i in range(0,len(name)-aux):
    l=name[i]
    m=name[i+aux] 
    
    Pop_staux_1=np.array(Pop_stoc[l])
    Pop_staux_2=np.array(Pop_stoc[m])

    flag1=np.min([np.max(Pop_staux_2), np.max(Pop_staux_1)])
    flag2=0
    flag=2000#(flag1-flag2)/2
    y1=np.max(Pop_staux_1[(Pop_staux_1==flag2)])  
    y2=np.min(Pop_staux_1[(Pop_staux_1>=flag1)]) 
    y1_=np.max(Pop_staux_2[(Pop_staux_2==flag2)])  
    y2_=np.min(Pop_staux_2[(Pop_staux_2>=flag1)]) 
    x1 = np.min(np.where((Pop_staux_1==flag2)))
    x2 = np.max(np.where((Pop_staux_1>=flag1)))
    x1_= np.min(np.where((Pop_staux_2==flag2)))  
    x2_= np.max(np.where((Pop_staux_2>=flag1))) 
    ct= (y2*x1-y1*x2)/(x1-x2)
    mt=(y1-y2)/(x1-x2)
    ct1= (y2_*x1_-y1_*x2_)/(x1_-x2_)
    mt1=(y1_-y2_)/(x1_-x2_)
    xt= (flag- ct)/mt 
    xt1= (flag- ct1)/mt1
    
    vel_st= abs(xt1-xt)/(t_step*aux)
    
    
    velocity_st.append(vel_st)
    

pop_sto,=plt.plot(velocity_st,'o-',linewidth=3)

plt.legend([pop_sto],['n'])

# plt.xticks(range(0,(100-aux)+1,(100-aux)/5),range(0,500000+1,500000/5), fontsize = 15) # work on current fig
plt.yticks( fontsize = 15) # work on current fig
plt.xticks( fontsize = 15)
plt.xlabel('t', fontsize=15)
plt.ylabel(r'TW Velocity', fontsize=15)
plt.tight_layout()
plt.savefig(outputdir+'/Velocity_TW_ng_st.eps')





##### Save velocity for all the different realizations to do the histogram
Pop_stoc_realization = {}
for file in name:
    for j in range(1,(Number_simulation+1)):
        foldername='DatosSimulacion/%s/OutputValuesPopulation/'%j
        aux=[]
        crs = open(foldername+file, "r")
        for columns in ( raw.strip().split() for raw in crs ):
            aux.append(columns[2])
        aux=np.array(map(float,aux))
        save_name = file +'%s'%(j-1)
        Pop_stoc_realization[save_name]=aux

Velocity_stoc_hist=[]
aux=1
for j in range(Number_simulation):
    for i in range(90,len(name)-aux):
        l=name[i]+'%s'%(j)
        m=name[i+aux] +'%s'%(j)
        Pop_staux_1=np.array(Pop_stoc_realization[l])
        Pop_staux_2=np.array(Pop_stoc_realization[m])

        flag1=np.min([np.max(Pop_staux_2), np.max(Pop_staux_1)])
        flag2=0
        flag=2000#(flag1-flag2)/2
        y1=np.max(Pop_staux_1[(Pop_staux_1==flag2)])  
        y2=np.min(Pop_staux_1[(Pop_staux_1>=flag1)]) 
        y1_=np.max(Pop_staux_2[(Pop_staux_2==flag2)])  
        y2_=np.min(Pop_staux_2[(Pop_staux_2>=flag1)]) 
        x1 = np.min(np.where((Pop_staux_1==flag2)))
        x2 = np.max(np.where((Pop_staux_1>=flag1)))
        x1_= np.min(np.where((Pop_staux_2==flag2)))  
        x2_= np.max(np.where((Pop_staux_2>=flag1))) 
        ct= (y2*x1-y1*x2)/(x1-x2)
        mt=(y1-y2)/(x1-x2)
        ct1= (y2_*x1_-y1_*x2_)/(x1_-x2_)
        mt1=(y1_-y2_)/(x1_-x2_)
        xt= (flag- ct)/mt 
        xt1= (flag- ct1)/mt1
        
        vel_st= abs(xt1-xt)/(t_step*aux)
        Velocity_stoc_hist.append(vel_st)
plt.figure()
nbin=100
plt.hist(Velocity_stoc_hist,bins=nbin,color='g',alpha=0.25,histtype='stepfilled',label='Velocity TW')
plt.hist(Velocity_stoc_hist, bins=nbin,color='g',histtype='step',lw=2)
plt.legend()

# plt.xticks(range(0,(100-aux)+1,(100-aux)/5),range(0,500000+1,500000/5), fontsize = 15) # work on current fig
plt.yticks( fontsize = 15) # work on current fig
plt.xticks( fontsize = 15)
plt.xlabel('Velocity', fontsize=15)
plt.tight_layout()
plt.savefig(outputdir+'/Velocity_TW_st_hist.eps')









