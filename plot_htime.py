#!/usr/bin/env python

import numpy as np
import glob
import os
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

chirality=[]
slopes=[]


path='/projectnb/qedksk/ashishge/data/new_sim/effect_of_N/samech_gdnoise/try6/'
destpath=path[:-1]+"_Hplots/"
if not os.path.exists(destpath): os.makedirs(destpath)


plot_from_average_folder=0 ### this can plot from average folder but then can't estimate error bars.

lvl1_prefix='ch'
lvl2_prefix='N'
r_IC=20
vel_list=[]
vel_stderr_list=[]
het_mean=[]
het_mean2=[]
het_stderr=[]
het_stderr2=[]
chlist=[]
Nlist=[]
foldercounter=0


def plot_Ht(h,t,destfile_name):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.set_title('$h(x) ')
    plt.plot(t,h,'b-')
    plt.savefig(destfile_name)
if plot_from_average_folder==1:
    for lvl1 in glob.glob(path+'*/'):  ## ch value
        print lvl1
        #chval=lvl1.replace(path,'')
        #chval=chval.replace(lvl1_prefix,'')
        #chval=float(chval.replace('/',''))
        #chlist.append(chval)
        run_counter=0
        het_mean.append([])       
        for fold in glob.glob(lvl1+'*/'):  
            print fold
            times=np.loadtxt(fold+"time.txt")
            het=np.loadtxt(fold+"HET_time.txt")
    
            tval=np.linspace(0,np.min(times),len(het))
            
            destfold=fold.replace(path,destpath)
            if not os.path.exists(destfold): os.makedirs(destfold)
            plot_Ht(het,tval,destfold+"Het_t.png")
            
            print "50% and 90%",np.mean(het[:len(het)/2]),np.mean(het[:9*len(het)/10])

        foldercounter+=1


elif plot_from_average_folder==0:             ### plot_from_average_folder=0!
    lvl1foldercounter=0
    for lvl1 in glob.glob(path+'*/'):  ## ch value
        #print lvl1
        chval=lvl1.replace(path,'')
        chval=chval.replace(lvl1_prefix,'')
        chval=float(chval.replace('/',''))        
        
        v_templist=[]
        het_mean.append([])
        het_mean2.append([])
        het_stderr.append([])
        het_stderr2.append([])
        chlist.append([])
        Nlist.append([])
        for lvl2 in glob.glob(lvl1+'*/'):  ## N value
            Nval=lvl2.replace(lvl1,'')
            Nval=Nval.replace(lvl2_prefix,'')
            Nval=float(Nval.replace('/',''))
            
            
            run_counter=0
            mean1=mean2=var1=var2=0
            for fold in glob.glob(lvl2+'*/'):  
                #print fold
                times=np.loadtxt(fold+"time.txt")
                het=np.loadtxt(fold+"heterozygosity_time.txt")
        
                tval=np.linspace(0,np.min(times),len(het))
                
                destfold=fold.replace(path,destpath)
                if not os.path.exists(destfold): os.makedirs(destfold)
                plot_Ht(het,tval,destfold+"Het_t.png")
                
                print "50% and 90%",np.mean(het[:len(het)/2]),np.mean(het[:9*len(het)/10])
                mean1+=np.mean(het[:len(het)/2])
                mean2+=np.mean(het[:9*len(het)/10])
                var1+=np.mean(het[:len(het)/2])**2
                var2+=np.mean(het[:9*len(het)/10])**2

                run_counter+=1
                
                
            mean1=mean1/ run_counter 
            mean2=mean2/ run_counter 
            var1=var1/run_counter -mean1**2
            var2=var2/run_counter -mean2**2
            
            het_mean[ lvl1foldercounter].append(mean1)
            het_mean2[ lvl1foldercounter].append(mean2)
            het_stderr[ lvl1foldercounter].append(np.sqrt(var1/run_counter))
            het_stderr2[ lvl1foldercounter].append(np.sqrt(var2/run_counter))
            
            Nlist[ lvl1foldercounter].append(Nval)
            chlist[ lvl1foldercounter].append(chval)
            
            
            print "ch=",chval,Nval,   mean1,"+-", np.sqrt(var1/run_counter),  mean2,"+-",  np.sqrt(var2/run_counter)       
        lvl1foldercounter+=1
        

        
    chlist=np.array(  chlist)  
    Nlist=np.array(  Nlist)  
    het_mean=np.array(het_mean)
    het_mean2=np.array(het_mean2)
    het_stderr=np.array(het_stderr)
    het_stderr2=np.array(het_stderr2)
    
    
    np.savetxt(destpath+"het_mean1.txt",het_mean)
    np.savetxt(destpath+"het_mean2.txt",het_mean2)
    np.savetxt(destpath+"het_stderr1.txt",het_stderr)
    np.savetxt(destpath+"het_stderr2.txt",het_stderr2)
    np.savetxt(destpath+"chlist.txt",chlist)
    np.savetxt(destpath+"Nlist.txt",Nlist)
    
    
    print chlist[0], Nlist[0], het_mean[0],het_stderr[0]
    print np.shape(Nlist)
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111) 
    for i in range(len(Nlist[0])):
        plt.errorbar(chlist[:,i],het_mean[:,i], yerr=  het_stderr[:,i], fmt='o')
    destfile_name=destpath+"width1_"+lvl1_prefix+"_"+lvl2_prefix+str(Nlist[0][i])
    plt.savefig(destfile_name+".png")     
    
    for i in range(len(Nlist[0])):
        plt.errorbar(chlist[:,i],het_mean2[:,i], yerr=  het_stderr2[:,i], fmt='o')
    destfile_name=destpath+"width2_"+lvl1_prefix+"_"+lvl2_prefix+str(Nlist[0][i])
    plt.savefig(destfile_name+".png")     
    
    

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111) 
    for i in range(len(chlist.T[0])):
        plt.errorbar(Nlist[i,:],het_mean[i,:], yerr=  het_stderr[i,:], fmt='o')
    destfile_name=destpath+"width1_"+lvl2_prefix+"_"+lvl1_prefix+str(chlist.T[0][i])
    plt.savefig(destfile_name+".png")      
            
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111) 
    for i in range(len(chlist.T[0])):
        plt.errorbar(Nlist[i,:],het_mean2[i,:], yerr=  het_stderr2[i,:], fmt='o')
    destfile_name=destpath+"width2_"+lvl2_prefix+"_"+lvl1_prefix+str(chlist.T[0][i])
    plt.savefig(destfile_name+".png")  