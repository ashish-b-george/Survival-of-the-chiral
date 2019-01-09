#!/usr/bin/python
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import glob
import os
import math
import numpy as np
import sys
sys.path.append('/usr3/graduate/ashishge/off_lattice/')
from off_lattice_population_class import population_Cpp

#lvl1foldertype="boxL"
lvl1foldertype="Ms"

path = '/projectnb/qedksk/ashishge/data/off_lattice/chirality_test/try1/'
destpath = '/projectnb/qedksk/ashishge/data/off_lattice/chirality_test/try1_avg/'
file_ctr_max=4
def plot_points(xval,yval,xlabel,ylabel,file_name):
    fig = plt.figure(figsize=(6,6))  
    ax = fig.add_subplot(111)
    plt.plot(xval,yval,'o',markeredgecolor='none') 
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)   
    #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)) 
    #plt.tight_layout()
    plt.savefig(file_name) 
    
def plot_points_and_theory(xval,yval,ytheory,xlabel,ylabel,file_name):
    fig = plt.figure(figsize=(6,6))  
    ax = fig.add_subplot(111)
    plt.plot(xval,yval,'bo',markeredgecolor='none') 
    plt.plot(xval,ytheory,'r-',markeredgecolor='none') 
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)   
    #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0)) 
    #plt.tight_layout()
    plt.savefig(file_name)     
    

velocity=[]
velocity_theoretical=[]
box_length=[]
for lvl1fold in glob.glob(path+'*/'):
    foldercounter=0
    print lvl1fold
    for fold in glob.glob(lvl1fold+'*/'):
        filectr=0
        for params_file in glob.glob(fold+"params*.txt"):
            filectr+=1
        if file_ctr_max!=0:
            filectr=file_ctr_max ## manual override if the runs quit from max population size
        
        if foldercounter==0:  
            meany=[]
            simulation_time=[]
            pop1=population_Cpp()
            for i in range (filectr):
                filename="pos_"+str(i)+".txt"
                pop1.read_data(fold,filename)
                meany.append(np.mean(pop1.ypos))
                simulation_time.append(pop1.simulation_time)
        else:
            for i in range (filectr):
                filename="pos_"+str(i)+".txt"
                pop1.read_data(fold,filename)
                meany[i]+=np.mean(pop1.ypos)
                simulation_time[i]+=pop1.simulation_time
        
        destfold=fold.replace(path, destpath)
        if not os.path.exists(destfold): os.makedirs(destfold)
        pop1.plot_histograms(destfold,"histogram"+str(filectr-1),5) ##plots for the last value of i only.
        pop1.plot_scatter(destfold,"scatter_labelled")
        pop1.plot_image(destfold,"concentration")
        foldercounter+=1 
        
    destfold=lvl1fold.replace(path, destpath)
    if not os.path.exists(destfold): os.makedirs(destfold)
    meany=np.array(meany)/foldercounter
    simulation_time=np.array(simulation_time)/foldercounter
    np.savetxt(destfold+"meany.txt",meany)
    np.savetxt(destfold+"simulation_time.txt",simulation_time)
    plot_points(simulation_time,meany,"simulation time","meany",destfold+"meany_t.png")
    Diffusion= (pop1.params.loc["migration_step_size"].values[0])**2. /4.  ### diffusion = (step^2/4)
    growth_rate=pop1.params.loc["growth_rate"].values[0]
    velocity_theoretical.append( 2.* math.sqrt(Diffusion * growth_rate) )    
    plt.close("all")
    velocity.append(2.* (meany[-1]-meany[len(meany)/2])/(simulation_time[-1]-simulation_time[len(meany)/2]) )
    bl=lvl1fold.replace(path,"")
    bl=bl.replace(lvl1foldertype,"")
    bl=bl.replace("/","")
    box_length.append(float(bl))
    
destfold=destpath
if not os.path.exists(destfold): os.makedirs(destfold)

plot_points(box_length,velocity,str(lvl1foldertype),"velocity",destfold+"vel_"+str(lvl1foldertype)+".png")

plot_points_and_theory(box_length,velocity,velocity_theoretical,str(lvl1foldertype),"velocity",destfold+"velPredict_"+str(lvl1foldertype)+".png")

plt.close("all")