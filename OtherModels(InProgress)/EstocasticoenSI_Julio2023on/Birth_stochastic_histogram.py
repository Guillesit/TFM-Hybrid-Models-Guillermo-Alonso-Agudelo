
# coding: utf-8


# Setting the terminal postscript with the options
import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import math as m


###### take files data name
foldername='DatosSimulacion/1/OutputValuesPopulation/'
name=[]
for i in range(0,10):
    for file in os.listdir(foldername):
        if file[:5]=="Out%s."%i:
#             print file
            name.append(file)
for i in range(10,100):
    for file in os.listdir(foldername):
        if file[:6]=="Out%s."%i:
#             print file
            name.append(file)
for file in os.listdir(foldername):
        if file[:7]=="Out%s."%100:
#             print file
            name.append(file)
 



###### Stochastich files information
# column[0] = space position (box)
# column[1] = birth events in this box
# column[2] = population in this box

threshold = 2000.0
Pop_stoc={}
Birth_stoc = {}
Birth_left = {}
Birth_right = {}
for file in name:
    j=1
    foldername='DatosSimulacion/%s/OutputValuesPopulation/'%j
    aux=[]
    aux_birth_left = []
    aux_birth_right = []
    birth = []
    crs = open(foldername+file, "r")
    for columns in ( raw.strip().split() for raw in crs ):
        aux.append(columns[2])
        if (float(columns[2])<threshold):
            if (float(columns[2])>0.0):
	            	# ratio of birth events between total number of cells, tp*nu*B
                    aux_birth_left.append(float(columns[1])/float(columns[2]))
                    birth.append(columns[1])
            else:
                birth.append(0.0)
        if (float(columns[2])>=threshold):
            aux_birth_right.append(float(columns[1])/float(columns[2]))
            birth.append(0.0)
    aux=np.array(map(float,aux))
    birth= np.array(map(float,birth))
    Pop_stoc[file]=aux
    Birth_stoc[file]=birth
    # Birth_left[file] = aux_birth_left
    # Birth_right[file] = aux_birth_right
    # aux_birth_left = []
    # aux_birth_right = []
    for j in range(1,101):
        foldername='DatosSimulacion/%s/OutputValuesPopulation/'%j
        aux=[]
        birth=[]
        crs = open(foldername+file, "r")
        for columns in ( raw.strip().split() for raw in crs ):
            aux.append(columns[2])
            if (float(columns[2])<threshold): 
            	if (float(columns[2])>0.0):
	            	# ratio of birth events between total number of cells, tp*nu*B
                    aux_birth_left.append(float(columns[1])/float(columns[2]))
                    birth.append(columns[1])
                else:
                    birth.append(0.0)
            if (float(columns[2])>=threshold):
                aux_birth_right.append(float(columns[1])/float(columns[2]))
                birth.append(0.0)
        aux=np.array(map(float,aux))
        birth= np.array(map(float,birth))
        Pop_stoc[file]=aux+Pop_stoc[file]
        Birth_stoc[file]=birth+Birth_stoc[file]
        Birth_left[file] = aux_birth_left
        Birth_right[file] = aux_birth_right
    Pop_stoc[file]=Pop_stoc[file]/100.0 # divide by number of simulations
    Birth_stoc[file]=Birth_stoc[file]/100.0


# # Empirical birth rate

# Values used in the simulations
# tp= 480 
# nu_inv = 0.000041667

"""Save create the folder to save the data"""
outputdir="images"
if not os.path.exists(outputdir): # if the folder doesn't exist create it
    os.makedirs(outputdir)

# i=0
# frames = []

# for m in name:
#     fig=plt.figure()
#     B_=Birth_right[m] #right
#     x_max = 0.05
#     binwidth = x_max/20.0
#     nbin=np.linspace(0,x_max,21)
#     plt.hist(B_, bins=nbin,color='g',alpha=0.25,histtype='stepfilled',label='Empirical birth rate')
#     plt.hist(B_, bins=nbin,color='g',histtype='step',lw=2)
#     plt.xlim(0, x_max)
#     plt.ylim(0,1400)
#     plt.xlabel(r'$\tau_p \nu B$')
#     lgd=plt.legend()
#     i = i+1
#     frame = "images/images_right%03d.png"%i
#     fig.savefig(frame, dpi=100)
#     frames.append(frame)
# plt.close('all')
# # make a video in linux
# # os.system("mencoder 'mf://images/image_left*.png' -mf type=png:fps=5 -ovc lavc -lavcopts vcodec=wmv2 -oac copy  -o " + "left.mpg")
# # make a video in Mac
# os.system("ffmpeg -framerate 5/1 -i images/images_right%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p birth_left.mp4") #for Mac computer

# # remove files in images
# for frame in frames: os.remove(frame)



i=0
frames = []

for m in name:
    fig=plt.figure()
    B_=Birth_left[m] #left
    x_max = 1.0
    binwidth = x_max/20.0
    nbin = np.linspace(0,x_max,210)
    plt.hist(B_, bins=nbin,color='g',alpha=0.25,histtype='stepfilled',label='Empirical birth rate')
    plt.hist(B_, bins=nbin,color='g',histtype='step',lw=2)
    plt.xlim(0, x_max)
    plt.ylim(0,70)
    plt.xlabel(r'$\tau_p \nu B$')
    lgd=plt.legend()
    i = i+1
    frame = "images/images_left%03d.png"%i
    fig.savefig(frame, dpi=100)
    frames.append(frame)
plt.close('all')
# make a video in linux
# os.system("mencoder 'mf://images/image_left*.png' -mf type=png:fps=5 -ovc lavc -lavcopts vcodec=wmv2 -oac copy  -o " + "left.mpg")
# make a video in Mac
os.system("ffmpeg -framerate 5/1 -i images/images_left%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p birth_right.mp4") #for Mac computer

# remove files in images
for frame in frames: os.remove(frame)


i=0
frames = []

for m in name:
    fig=plt.figure()
    B_=Birth_stoc[m] #left
    plt.plot(B_, 'g-',label='Number of births')
    plt.xlim(0, 100)
    # plt.ylim(0,70)
    # plt.xlabel(r'$\tau_p \nu B$')
    lgd=plt.legend()
    i = i+1
    frame = "images/images_births%03d.png"%i
    fig.savefig(frame, dpi=100)
    frames.append(frame)
plt.close('all')
# make a video in linux
# os.system("mencoder 'mf://images/image_left*.png' -mf type=png:fps=5 -ovc lavc -lavcopts vcodec=wmv2 -oac copy  -o " + "left.mpg")
# make a video in Mac
os.system("ffmpeg -framerate 5/1 -i images/images_births%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p births_stoc.mp4") #for Mac computer

# remove files in images
for frame in frames: os.remove(frame)

i=0
frames = []

for m in name:
    fig=plt.figure()
    B_=Pop_stoc[m] #left
    plt.plot(B_, 'g-',label='Stoc')
    plt.xlim(0, 100)
    # plt.ylim(0,70)
    # plt.xlabel(r'$\tau_p \nu B$')
    lgd=plt.legend()
    i = i+1
    frame = "images/images_pop%03d.png"%i
    fig.savefig(frame, dpi=100)
    frames.append(frame)
plt.close('all')
# make a video in linux
# os.system("mencoder 'mf://images/image_left*.png' -mf type=png:fps=5 -ovc lavc -lavcopts vcodec=wmv2 -oac copy  -o " + "left.mpg")
# make a video in Mac
os.system("ffmpeg -framerate 5/1 -i images/images_pop%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p pop_stoc.mp4") #for Mac computer

# remove files in images
for frame in frames: os.remove(frame)


