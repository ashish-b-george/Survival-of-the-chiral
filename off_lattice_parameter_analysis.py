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
from copy import deepcopy
sys.path.append('/usr3/graduate/ashishge/off_lattice/')
#sys.path.append('/Users/ashish/Dropbox/research/off lattice simulations/python codes/')

from off_lattice_population_class import population_Cpp
import argparse
import shutil

path = '/projectnb/qedksk/ashishge/data/off_lattice/systematic_analysis/oppch_inv/2try1smallK/'
parameter_type="ch"
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="base_path")
parser.add_argument("-p", help="parameter_type")
args = parser.parse_args()
if args.i:
    print "input path given"
    path =args.i
    print path
if args.p:
    print "parameter type given"
    parameter_type =args.p
    print parameter_type

destpath = path[:-1]+'_analysis/'
if not os.path.exists(destpath): os.makedirs(destpath)





color_cycler=['b','r','g','y','black','c']
def plot_points(xarr,yarr,xlabel,ylabel,file_name,label_arr=None):
    fig = plt.figure(figsize=(6,6))  
    ax = fig.add_subplot(111)
    if len(np.shape(xarr))==1:
        plt.plot(xarr,yarr,'o',markeredgecolor='none') 
    elif np.shape(xarr)[0]==1: ## if there is an extra dummy dimension!
        plt.plot(xarr[0],yarr[0],'o',markeredgecolor='none') 
    else:        
        for i in range(len(xarr)):
            print i
            if label_arr!=None:
                plt.plot(xarr[i],yarr[i],'o',markeredgecolor='none',markerfacecolor=color_cycler[i],label=label_arr[i%len(color_cycler)])  
            else:   
                plt.plot(xarr[i],yarr[i],'o',markeredgecolor='none',markerfacecolor=color_cycler[i%len(color_cycler)])  
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)  
    if label_arr!=None:
        legend = plt.legend(loc="best",fontsize=10)  
    plt.savefig(file_name+".png") 


#/projectnb/qedksk/ashishge/data/off_lattice/systematic_analysis/ch/try1/Mfunc2/ms0.5/ch0.0/1



for lvl1_fold in glob.glob(path+'*/*/'):
    
    lvl1_suffix=lvl1_fold.replace(path,"_")
    lvl1_suffix=lvl1_suffix.replace("/","_")
    
    vel_p=[]
    globalhet_p=[]
    param_list=[]
    roughness_p=[]
    MacroscopicChirality_p=[]
    for param_fold in glob.glob(lvl1_fold+'*/'): ## folder that is different for each parameter.
        print param_fold
        
        param_val=param_fold.replace(lvl1_fold,"")
        param_val=param_val.replace("/","")
        param_val=param_val.replace("_","")
        param_val=param_val.replace(parameter_type,"")
        param_list.append(float(param_val))
        
        foldercounter=0
        meany=[]
        simulation_time=[]
        fraction_a=[]
        global_het=[]
        roughness=[]
        boundary1_list=[]
        boundary1flag_list=[]
        boundary2_list=[]
        boundary2flag_list=[]
        avg_height=np.zeros(1)
        for fold in glob.glob(param_fold+'*/'):
            filectr=0
            file_suffix_list=[]
            meany.append([])
            simulation_time.append([])
            fraction_a.append([])
            global_het.append([])
            roughness.append([])
            boundary1_list.append([])
            boundary1flag_list.append([])
            boundary2_list.append([])
            boundary2flag_list.append([])
            
            destfold=fold.replace(path, destpath)
            if not os.path.exists(destfold): os.makedirs(destfold)
            for params_file in glob.glob(fold+"params_*.txt"):
                file_suffix=params_file.replace(fold+"params_","")
                file_suffix=file_suffix.replace(".txt","")
                file_suffix_list.append(file_suffix)
                filectr+=1
            for i in range (filectr):
                pop1=population_Cpp()
                filename="pos_"+file_suffix_list[i]+".txt"
                pop1.read_data(fold,filename)
                meany[foldercounter].append(np.mean(pop1.ypos))
                simulation_time[foldercounter].append(pop1.simulation_time)
                
                
                #fraction_a[foldercounter].append( np.mean(pop1.species_label) )
                #global_het[foldercounter].append( np.mean(pop1.species_label)*np.mean(1-pop1.species_label) )
                fatemp,hettemp=pop1.calculate_Fraction_Globalhet_at_front()
                fraction_a[foldercounter].append(fatemp) ## fraction of species in the movibgbox
                global_het[foldercounter].append( hettemp ) ## global het at front              
                pop1.calculate_height_and_roughness() ## to calculate roughness
                roughness[foldercounter].append( pop1.roughness )
                
                
                
                pop1.plot_image(destfold,"concentration_"+file_suffix_list[i])
                pop1.plot_scatter(destfold,"scatter_labelled_"+file_suffix_list[i])
                pop1.plot_height(destfold,"height_"+file_suffix_list[i])
                
                plt.close("all")
            
            if "final" not in file_suffix_list:
                max_suffix=str(0)
                for i in file_suffix_list :
                    if float(max_suffix)<float(i): 
                        max_suffix=i 
                filename="pos_"+ max_suffix+".txt"
                pop1.read_data(fold,filename)
                pop1.plot_sensedC(fold,"sensedC_"+max_suffix,destfold)
                if len(avg_height)==1:
                    avg_height=deepcopy(np.ravel(pop1.Height_x) )
                else:
                    avg_height+=np.ravel(pop1.Height_x)
                
            else:
                
                filename="pos_final.txt"
                pop1.read_data(fold,filename) 
                pop1.plot_sensedC(fold,"sensedC_"+"final",destfold) 
                
                
            pop1.plot_histograms(destfold,"histogram"+str(filectr-1),5) ##plots for the last value of i only.
            pop1.plot_scatter(destfold,"scatter_labelled")
            pop1.plot_image(destfold,"concentration")
            pop1.plot_height(destfold,"height")
            pop1.plot_domain_boundaries(destfold,"boundary")
            

            
            
            
            boundary1,boundary2,boundary_flag1,boundary_flag2 =pop1.calculate_boundaries()
            
            boundary1_list[foldercounter].append(boundary1)
            boundary1flag_list[foldercounter].append(boundary_flag1)
            boundary2_list[foldercounter].append(boundary2)
            boundary2flag_list[foldercounter].append(boundary_flag2)
            
            
            
            
            

            plt.close("all")
            
            foldercounter+=1 
        
            
        destfold=param_fold.replace(path, destpath)
        if not os.path.exists(destfold): os.makedirs(destfold)
        print "shape of simulation time",np.shape(simulation_time)
        plot_points(simulation_time,meany,"time","mean Y position ",destfold+"mean_y")
        plot_points(simulation_time,fraction_a,"time","fa ",destfold+"fraction_front")       
        plot_points(simulation_time[0],np.mean(fraction_a,axis=0),"time","mean_fa ",destfold+"mean_fraction_front")        
        plot_points(simulation_time,global_het,"time","global het ",destfold+"global_het_at_front")
        plot_points(simulation_time,roughness,"time","roughness ",destfold+"roughness")      
        
        np.savetxt(destfold+"simulation_time.txt",simulation_time)
        np.savetxt(destfold+"mean_ypos.txt",meany)
        np.savetxt(destfold+"fraction_front.txt",fraction_a)
        np.savetxt(destfold+"mean_fa.txt",np.mean(fraction_a,axis=0) )
        np.savetxt(destfold+"global_het.txt",global_het)
                
        if len(avg_height)!=1:
           
            plot_points(np.arange(len(avg_height)),avg_height,"x","height_height ",destfold+"height_sum")     
        
        
        if foldercounter>1:
            v=0
            gh=0
            rough=0
            macro_ch=0
            for i in range(foldercounter):
                v+= 2.*np.polyfit(simulation_time[i], meany[i], 1)[0]  
                gh+=  global_het[i][np.argmax(simulation_time[i])]   
                rough+= roughness[i][np.argmax(simulation_time[i])] 
                
                
                xval=np.ravel(np.where(boundary1flag_list[i][0]==1))
                yval=np.ravel(boundary1_list[i])[np.where(boundary1flag_list[i][0]==1)]
                
               # print "len(xval)",len(xval), boundary1flag_list[i]
                if len(xval) >2: 
                    macro_ch+=np.polyfit(xval, yval, 1)[0]
                else:
                    print " not enough yvals to fit!!!1", len(yval),fold
                xval=np.ravel(np.where(boundary2flag_list[i][0]==1) )
                yval=np.ravel(boundary2_list[i])[np.where(boundary2flag_list[i][0]==1)]
                if len(xval) >2:
                    macro_ch+=np.polyfit(xval, yval, 1)[0]
                else:
                    print " not enough yvals to fit!!!1", len(yval),fold
                
                
                
                
                
            vel_p.append(v/foldercounter)
            globalhet_p.append(gh/foldercounter)
            roughness_p.append(rough/foldercounter)
            MacroscopicChirality_p.append(macro_ch/(2.*foldercounter))
        else:
            macro_ch=0
            vel_p.append(2.*np.polyfit(simulation_time[0], meany[0], 1)[0]  )
            globalhet_p.append( global_het[0][np.argmax(simulation_time[0])] ) 
            roughness_p.append( roughness[0][np.argmax(simulation_time[0])] )
            xval=np.ravel(np.where(boundary1flag_list[0][0]==1))
            yval=np.ravel(boundary1_list)[np.where(boundary1flag_list[0][0]==1)]
            if len(xval) >2:
                macro_ch+=np.polyfit(xval, yval, 1)[0]
            else:
                print " not enough yvals to fit!!!1", len(yval),fold
            xval=np.ravel(np.where(boundary2flag_list[0][0]==1))
            yval=np.ravel(boundary2_list)[np.where(boundary2flag_list[0][0]==1)]
            if len(xval) >2:
                macro_ch+=np.polyfit(xval, yval, 1)[0]
            else:
                print " not enough yvals to fit!!!1", len(yval),fold
            
            MacroscopicChirality_p.append( macro_ch/2.)
            
           
        #destfold=destpath
        #
        #
        #
        #
        #if not os.path.exists(destfold): os.makedirs(destfold)
        
        
        plot_points(param_list,vel_p,parameter_type,"velocity ",destfold+"vel"+lvl1_suffix)
        plot_points(param_list,globalhet_p,parameter_type,"global heterozygosity f(1-f) ",destfold+"globalHet_front"+lvl1_suffix)
        plot_points(param_list,roughness_p,parameter_type,"roughness ",destfold+"roughness"+lvl1_suffix)        
        plot_points(param_list,MacroscopicChirality_p,parameter_type,"macroscopic chirality ",destfold+"macroscopic_chirality"+lvl1_suffix)
        np.savetxt(destfold+"param_list.txt",param_list)
        np.savetxt(destfold+"MacroscopicChirality.txt",MacroscopicChirality_p)
        np.savetxt(destfold+"velocity.txt",vel_p)
        np.savetxt(destfold+"globalhet.txt",globalhet_p)
  
plt.close("all")


for submission_script in glob.glob(path+'*.sh'):
    print "copying submission script ", submission_script
    shutil.copy(submission_script, submission_script.replace(path, destpath))