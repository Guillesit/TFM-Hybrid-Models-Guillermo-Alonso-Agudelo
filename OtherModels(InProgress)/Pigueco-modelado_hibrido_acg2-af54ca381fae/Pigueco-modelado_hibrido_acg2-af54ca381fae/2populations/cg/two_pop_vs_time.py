#####################################################################
# Scrip to plot two population vs time  (cg with i.c. in equilibrio)#
#####################################################################


# Setting the terminal postscript with the seaborn options
import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import math as m
from glob import glob
import seaborn as sns

###COLORBAR 
# Using contourf to provide my colorbar info, then clearing the figure
plot_step = 5
num_plots=(100/plot_step+1)
colormap =sns.cubehelix_palette( start=1,rot=1, as_cmap=True)#plt.cm.plasma#gist_ncar
sc=plt.gca().set_color_cycle([colormap(j) for j in np.linspace(0, .9, num_plots)])
min, max = (0, 500000)
step = 500000/(num_plots)
Z = [[0,100],[0,4000]]
levels = range(min,max+step,step)
CS3 = plt.contourf(Z, levels, cmap=colormap)
plt.clf()
plt.close()
# cbar = plt.colorbar(CS3, ticks=range(0,550000,500000/5))
# cbar.set_label('time', rotation=90)

#read files in order of name vector
name=[]
for i in range(0,101):
	file="Out%s"%i+".txt"
	name.append(file)
	# print file
output_files_tag = '7'
foldername_host = 'OutputValuesPopulationHost'+output_files_tag+'/'
foldername_invader = 'OutputValuesPopulationInvader'+output_files_tag+'/'
foldername_oxygen = 'OutputValuesOxygen'+output_files_tag+'/'

#population information
Pop_Host={}
Pop_Invader={}
#Oxygen data
Oxygen={}

for file in os.listdir(foldername_host):
	aux=[]
	if file[:3]=="Out":
		#print file
		crs = open(foldername_host+file, "r")
		for columns in ( raw.strip().split() for raw in crs ):
		    aux.append(columns[1])
		aux=np.array(map(float,aux))
		Pop_Host[file]=aux


for file in os.listdir(foldername_invader):
	aux=[]
	if file[:3]=="Out":
		#print file
		crs = open(foldername_invader+file, "r")
		for columns in ( raw.strip().split() for raw in crs ):
		    aux.append(columns[1])
		aux=np.array(map(float,aux))
		Pop_Invader[file]=aux

for file in os.listdir(foldername_oxygen):
	aux=[]
	if file[:3]=="Out":
		#print file
		crs = open(foldername_oxygen+file, "r")
		for columns in ( raw.strip().split() for raw in crs ):
		    aux.append(columns[1])
		aux=np.array(map(float,aux))
		Oxygen[file]=aux


## Pictures of $P_{host}$ and $P_{invader}$ at same time (Population and Oxigen)

#Same time snapshot of Population both populations
plt.figure()
sns.set_style("whitegrid")
colormap =sns.cubehelix_palette( start=1,rot=1, as_cmap=True)#plt.cm.plasma#gist_ncar
plt.gca().set_color_cycle([colormap(j) for j in np.linspace(0, .9, 2*num_plots)])
for l in name[0::plot_step]:
    Pop_aux=np.array(Pop_Host[l])
    Pop_maux=np.array(Pop_Invader[l])
    pop_host,=plt.plot(Pop_aux,'--',linewidth=3)
    pop_inv,=plt.plot(Pop_maux,'-',linewidth=3)
plt.title('Host and Invader')
plt.legend([pop_host, pop_inv], ['$n_{host}$', '$n_{invader}$'],frameon=True, fancybox=True, shadow=True,fontsize=15)
# plt.xticks(range(0,110,100/5),range(0,550000,500000/5), fontsize = 15) # work on current fig
plt.yticks( fontsize = 15) # work on current fig
plt.xticks( fontsize = 15)
#plt.xlim(0,100)
# plt.ylim(0,4500)
plt.xlabel('space (X)', fontsize=15)
plt.ylabel('Population', fontsize=15)
# plt.ylim(0,100000)
cbar = plt.colorbar(CS3, ticks=range(0,550000,500000/5))
cbar.set_label('time', rotation=90, fontsize=15)
plt.savefig('two_pop'+output_files_tag+'.eps', dpi=100, bbox_inches='tight')
plt.show()

#Host population
plt.figure()
sns.set_style("whitegrid")
colormap =sns.cubehelix_palette( start=1,rot=1, as_cmap=True)#plt.cm.plasma#gist_ncar
plt.gca().set_color_cycle([colormap(j) for j in np.linspace(0, .9, num_plots)])
for l in name[0::plot_step]:
	Pop_aux=np.array(Pop_Host[l])
	pop_host,=plt.plot(Pop_aux,'--',linewidth=3)
plt.title('Host')
plt.legend([pop_host], ['$n_{host}$'],frameon=True, fancybox=True, shadow=True,fontsize=15)
# plt.xticks(range(0,110,100/5),range(0,550000,500000/5), fontsize = 15) # work on current fig
plt.yticks( fontsize = 15) # work on current fig
plt.xticks( fontsize = 15)
#plt.xlim(0,100)
# plt.ylim(0,4500)
plt.xlabel('space (X)', fontsize=15)
plt.ylabel('Population', fontsize=15)
cbar = plt.colorbar(CS3, ticks=range(0,550000,500000/5))
cbar.set_label('time', rotation=90, fontsize=15)
plt.savefig('host_pop'+output_files_tag+'.eps', dpi=100, bbox_inches='tight')
plt.show()

#Invader population
plt.figure()
sns.set_style("whitegrid")
colormap =sns.cubehelix_palette( start=1,rot=1, as_cmap=True)#plt.cm.plasma#gist_ncar
plt.gca().set_color_cycle([colormap(j) for j in np.linspace(0, .9, num_plots)])
for l in name[0::plot_step]:
	Pop_maux=np.array(Pop_Invader[l])
	pop_inv,=plt.plot(Pop_maux,'-',linewidth=3, label='%s'%l)
	plt.legend()
plt.title('Invader population')
plt.legend([pop_inv], ['$n_{invader}$'],frameon=True, fancybox=True, shadow=True,fontsize=15)
# plt.xticks(range(0,110,100/5),range(0,550000,500000/num_plots), fontsize = 15) # work on current fig
plt.yticks( fontsize = 15) # work on current fig
plt.xticks( fontsize = 15)
# #plt.xlim(0,5)
# plt.ylim(0,4500)
# plt.ylim(50000,120000000) #para hacer zoom
# plt.ylim(4000,700000)
plt.xlabel('space (X)', fontsize=15)
plt.ylabel('Population', fontsize=15)
cbar = plt.colorbar(CS3, ticks=range(0,550000,500000/5))
cbar.set_label('time', rotation=90, fontsize=15)
plt.savefig('invader_pop'+output_files_tag+'.eps', dpi=100, bbox_inches='tight')
plt.show()

#Oxigen
plt.figure()
sns.set_style("whitegrid")
colormap =sns.cubehelix_palette( start=1,rot=1, as_cmap=True)#plt.cm.plasma#gist_ncar
plt.gca().set_color_cycle([colormap(j) for j in np.linspace(0, .9, num_plots)])
for l in name[0::plot_step]:
    Oxy_aux=np.array(Oxygen[l])
    oxy,=plt.plot(Oxy_aux,'-',linewidth=3)
plt.title('Oxygen')
# plt.legend([oxy], ['Oxygen'],frameon=True, fancybox=True, shadow=True,fontsize=15)
# plt.xticks(range(0,110,100/5),range(0,550000,500000/5), fontsize = 15) # work on current fig
plt.yticks( fontsize = 15) # work on current fig
plt.xticks( fontsize = 15)
#plt.xlim(0,100)
plt.xlabel('space (X)', fontsize=15)
plt.ylabel('Oxygen', fontsize=15)
cbar = plt.colorbar(CS3, ticks=range(0,550000,500000/5))
cbar.set_label('time', rotation=90, fontsize=15)
plt.savefig('oxygen_'+output_files_tag+'.eps', dpi=100, bbox_inches='tight')
plt.show()
