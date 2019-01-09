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


path='/projectnb/qedksk/ashishge/data/new_sim/velocity_with_oppch/gdnoise/try1/'
destpath=path[:-1]+"_vplots/"
if not os.path.exists(destpath): os.makedirs(destpath)

lvl1_prefix='ch'
lvl2_foldername='h0.0'
r_IC=20
vel_list=[]
vel_stderr_list=[]
xb_time_list=[]
timeArr_list=[]
chlist=[]

foldercounter=0
for lvl1 in glob.glob(path+'*/'):  ## ch value

    print lvl1
    chval=lvl1.replace(path,'')
    chval=chval.replace(lvl1_prefix,'')
    chval=float(chval.replace('/',''))
    chlist.append(chval)
    run_counter=0
    v_templist=[]
    xb_time_list.append([])
    timeArr_list.append([])
    for fold in glob.glob(lvl1+lvl2_foldername+'/*/'):  
        #with open(fold+'ca.txt') as file:
        #    ca = [[float(digit) for digit in line.split()] for line in file]

        with open(fold+'c.txt') as file:       
            c = [[float(digit) for digit in line.split()] for line in file] 
            
        with open(fold+'time.txt') as file:
            time_temp = [[float(digit) for digit in line.split()] for line in file]
        time=time_temp[0][0]
        c=np.array(c)              
        l= len(c)
        b=len(c.T)
        
        with open(fold+'c2.txt') as file:       
            c2 = [[float(digit) for digit in line.split()] for line in file]             
        with open(fold+'time2.txt') as file:
            time_temp = [[float(digit) for digit in line.split()] for line in file]
        time2=time_temp[0][0]        
        velocity=(np.sum(c)/l-np.sum(c2)/l)  /  (time-time2)

        #### if c2 and time 2 dont exist we can use the following method to compute velocity
        #velocity=(np.sum(c)/l-r_IC)  /time
        
        
        xbtime=np.loadtxt(fold+"sumc_array.txt")/l
        time_arr=np.loadtxt(fold+"time_array.txt")
        
        xb_time_list[foldercounter].append(xbtime)
        timeArr_list[foldercounter].append(time_arr)
        v_templist.append(velocity)
        run_counter+=1
    #print time-time2, time , time2
        
    vel_list.append(np.mean(v_templist)) 
    stderr=np.sqrt(np.var(v_templist))/np.sqrt(run_counter)
    vel_stderr_list.append( stderr)
    foldercounter+=1
 
    
chirality_of_strain1 =2*(np.array(chlist)-0.05)        
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111) 
ax.plot(chirality_of_strain1,vel_list,'bo',markersize=5,markeredgecolor=None)           
destfile_name=destpath+"velocity_"+lvl1_prefix+"_"+str(lvl2_foldername)
plt.savefig(destfile_name+".png")
    
    
print vel_stderr_list     



fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111) 
ax.plot(chirality_of_strain1,vel_list,'bo',markersize=5,markeredgecolor=None) 
ax.errorbar(chirality_of_strain1,vel_list,yerr=  vel_stderr_list , fmt='o',markersize=5)      
destfile_name=destpath+"velocity_"+lvl1_prefix+"_"+str(lvl2_foldername)
plt.savefig(destfile_name+"_err.png")

color_cycler=['b','y','g','r','cyan','brown','orange','black']
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111) 
for i in range(len(chirality_of_strain1)):
    for j in range(len(xb_time_list[i])):
        if j==0: ## with label
            ax.plot(timeArr_list[i][j],xb_time_list[i][j],'o',color=color_cycler[i%len(color_cycler)],label="ch="+str(chirality_of_strain1[i]),markeredgecolor=None,markevery=100  ) 
        else: ## no label
            ax.plot(timeArr_list[i][j],xb_time_list[i][j],'o',color=color_cycler[i%len(color_cycler)],markeredgecolor=None,markevery=100) 
plt.legend(loc="best"  ) 
destfile_name=destpath+"boundary_in_time_"+lvl1_prefix+"_"+str(lvl2_foldername)
plt.savefig(destfile_name+".png")



np.savetxt(destpath+"chlist"+lvl2_foldername+".txt",chlist)
np.savetxt(destpath+"chirality_of_strain1"+lvl2_foldername+".txt",chirality_of_strain1)
np.savetxt(destpath+"velocity_"+lvl1_prefix+"_"+str(lvl2_foldername)+".txt",vel_list)
np.savetxt(destpath+"velocity_stderr_"+lvl1_prefix+"_"+str(lvl2_foldername)+".txt",vel_stderr_list)

np.savetxt(destpath+"runlist"+lvl2_foldername+".txt",np.arange( len(xb_time_list[0]) ) )

for i in range(len(chirality_of_strain1)):
    for j in range(len(xb_time_list[i])):
        np.savetxt(destpath+"xbtime"+lvl2_foldername+"_"+str(i)+"_"+str(j)+".txt",xb_time_list[i][j])
        np.savetxt(destpath+"timeArr"+lvl2_foldername+"_"+str(i)+"_"+str(j)+".txt",timeArr_list[i][j])

            