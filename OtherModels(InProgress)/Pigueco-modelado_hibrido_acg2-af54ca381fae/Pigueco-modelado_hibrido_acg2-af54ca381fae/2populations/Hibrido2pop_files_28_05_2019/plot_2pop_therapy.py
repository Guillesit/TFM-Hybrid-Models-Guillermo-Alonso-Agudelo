
# coding: utf-8


# Setting the terminal postscript with the options
import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import math as m
import seaborn as sns

extension = '_00'#'_2pop_p3p61d001_100_1'
#case = 'PruebaA2'
fromDropbox=''
n_space = 100
###### take files data name
foldername_Host = fromDropbox +'OutputValuesHostPopulation'+extension+'/'
foldername_Invader = fromDropbox +'OutputValuesInvaderPopulation'+extension+'/'
foldername_Oxygen = fromDropbox +'OutputValuesOxygen'+extension+'/'
name=[]
for i in range(0,10):
    for file in os.listdir(foldername_Oxygen):
        if file[:5]=="Out%s."%i:
#             print file
            name.append(file)
# for i in range(10,100):
#     for file in os.listdir(foldername_Oxygen):
#         if file[:6]=="Out%s."%i:
# #             print file
#             name.append(file)
# for file in os.listdir(foldername_Oxygen):
#         if file[:7]=="Out%s."%100:
# #             print file
#             name.append(file)
 



###### Stochastich files information
# column[0] = population in this box

#Population data
Pop_Host={}
Pop_Invader={}


#Oxygen data
Pop_ox={}

for file in os.listdir(foldername_Host):
    aux=[]
    if file[:3]=="Out":
        crs = open(foldername_Host+file, "r")
        for columns in ( raw.strip().split() for raw in crs ):
            aux.append(columns[0])
        aux=np.array(map(float,aux))
        Pop_Host[file]=aux

for file in os.listdir(foldername_Invader):
    aux=[]
    if file[:3]=="Out":
        crs = open(foldername_Invader+file, "r")
        for columns in ( raw.strip().split() for raw in crs ):
            aux.append(columns[0])
        aux=np.array(map(float,aux))
        Pop_Invader[file]=aux

for file in os.listdir(foldername_Oxygen):
    aux=[]
    if file[:3]=="Out":
        crs = open(foldername_Oxygen+file, "r")
        for columns in ( raw.strip().split() for raw in crs ):
            aux.append(columns[0])
        aux=np.array(map(float,aux))
        Pop_ox[file]=aux



i=0
frames=[]
"""Save images in a folder before video"""
outputdir="images"
if not os.path.exists(outputdir): # if the folder doesn't exist create it
    os.makedirs(outputdir)
for m in name:
#     print m
    fig=plt.figure()
    #sns.set_style("whitegrid")
    Pop_Host_time=np.array(Pop_Host[m])
    Pop_Invader_time=np.array(Pop_Invader[m])
    Pop_ox_time = np.array(Pop_ox[m])
    #x=np.array(range(0,len(Pop_time)))
    plt.plot(Pop_Host_time,'r-',label='Host',linewidth=3)
    plt.plot(Pop_Invader_time,'b--',label='Invader',linewidth=3)
    # plt.plot(Pop_ox_time*4000,'g-',label='Oxygen',linewidth=3)
    # plt.title(case)
    plt.legend()#'Time %s'%i)
    plt.xlabel('space (X)', fontsize=15)
    plt.ylabel(r'Population', fontsize=15)
    plt.xlim(0,n_space)
    plt.ylim(0,6000)
    #plt.show()
    i=i+1
    frame="images/image%03i.png" % i
    fig.savefig(frame,dpi=100)
    frames.append(frame)  

#os.system("mencoder 'mf://images/image"+case+"*.png' -mf type=png:fps=5 -ovc lavc -lavcopts vcodec=wmv2 -oac copy  -o " + "prueba"+case+".mpg")        
#for frame in frames: os.remove(frame)


# make a video in linux
# os.system("mencoder 'mf://images/image_right*.png' -mf type=png:fps=5 -ovc lavc -lavcopts vcodec=wmv2 -oac copy  -o " + "left.mpg")
# make a video in Mac
os.system("ffmpeg -framerate 5/1 -i images/image%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p pop_"+extension+".mp4") #for Mac computer

