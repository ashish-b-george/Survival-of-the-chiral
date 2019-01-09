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
import argparse

path = '/projectnb/qedksk/ashishge/data/off_lattice/oppchiral_destination/try1/'
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="base_path")
args = parser.parse_args()
if args.i:
    print "input path given"
    path =args.i

destpath =path [:-1]+'_plots/'






## change number of */ in first for loop to number of loop levels as required.
file_ctr_max=0
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
    

for lvl1fold in glob.glob(path+'*/*/'):
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
    #Diffusion= (pop1.params.loc["migration_step_size"].values[0])**2. /4.  ### diffusion = (step^2/4)
    #growth_rate=pop1.params.loc["growth_rate"].values[0]   
    plt.close("all")

plt.close("all")
