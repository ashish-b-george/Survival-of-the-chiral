#!/usr/bin/env python
import numpy as np
import glob
import os
import math
from sys import exit
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

chirality=[]
slopes=[]


path='/projectnb/qedksk/ashishge/data/new_sim/bulge_slope_from_track_height/gdnoise/try3/'
destpath=path[:-1]+"_wallExtract/"
if not os.path.exists(destpath): os.makedirs(destpath)

lvl1_prefix='ch'

r_IC=100
beta_list=[]
vel_list=[]
chlist=[]


def plot_figure(ca,cb,destfile_name,color_scheme,tick_flag=None):
    l= len(ca)
    b=len(ca.T)
    
    fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(b/100.0,l/100.0))   
    img = np.zeros((l,b,3))

    if color_scheme==0: #red green
        img[:,:,0]=ca[:,:]
        img[:,:,1]=cb[:,:]
    img3 = ax.imshow(img,aspect="equal")  
    if tick_flag==None or tick_flag==0:   
        ax.set_xticks([])
        ax.set_yticks([])
    
    plt.savefig(destfile_name) 



def find_wall(li,lf,topb,botb):
    wall=np.zeros(lf)
    wall_flag=0
    for i in range(li, lf):
        dcslice=dc[:,i]
        errormin=float(4.0)
        ls=botb           
        location =-1000
        break_condition=1
        while (ls>topb) :
#                    if ls>0:
            address=np.arange(ls,ls+21)
            #if(ls+21>l): # if the wall hasnt reached the edges of the domain
            #    if (ls<l):
            #         address=np.arange(ls,ls+21)                        
            #error= np.sum(dcslice[address])
            error=np.sum(dcslice.take(address, mode='wrap'))

            if (abs(error)<abs(errormin)):
                location =ls
                errormin=error

            ls=ls-1
            
            if (ls==topb+1 and location==-1000): # a solution to minimise was not found!!
                #print "using this condition"
                if break_condition==1: # we have not tried seraching in an expanded domain!
                    topb=int(np.mean(wall[li:i]))-int(3*gap)
                    botb=int(np.mean(wall[li:i]))+int(3*gap)
                    
                    ls=botb
                    #print ls
                    #print topb
                    break_condition=0
                else:
                    print "still not found.. wtf,"
                    lf=i-10
                    wall_flag=1
                    break
                    #location =-1000             
        if wall_flag==1:
            print " breaking, wall not found"
            break
        topb=location-6 # we just search 20 points on either side of the previous wall point.
        botb=location+6     
        wall[i]=(location+10)
    return wall[li:lf],li,lf,wall_flag



for lvl1 in glob.glob(path+'*/*/'):  ## ch value

    run_counter=0
    v_templist=[]
    beta_avg=0
    for fold in glob.glob(lvl1+'/*/'):
        print fold  
        with open(fold+'ca.txt') as file:
            ca = [[float(digit) for digit in line.split()] for line in file]

        with open(fold+'c.txt') as file:       
            c = [[float(digit) for digit in line.split()] for line in file] 
            
        with open(fold+'time.txt') as file:
            time_temp = [[float(digit) for digit in line.split()] for line in file]
        time=time_temp[0][0]
        
        
        c=np.array(c) 
        ca=np.array(ca)              
        l= len(c)
        b=len(c.T)
        ## doing a roll to extract both walls that move so that zero crossings sees two walls 
        ca=np.roll(ca,30,axis=0)
        print np.mean(c[0])
        c=np.roll(c,30,axis=0) 
        print np.mean(c[0])
        cb=c-ca
        
        
        
        destfold=fold.replace(path,destpath)   
        if not os.path.exists(destfold): os.makedirs(destfold)     
        image_file_name=destfold+"ca.png"
        
        plot_figure(ca,cb,image_file_name,0)

        
        lexp=np.where (c[l/2,:]<0.8)[0][0]-10  
        print lexp
        
        dc=2*ca-c
        zero_crossings = np.where(np.diff(np.sign(dc[:,0])))[0] # all the zero crossing gives you number of initial walls 
        print  zero_crossings
        
        lf =lexp# final length for domain wall after this value, the domain walls coalesce
        li =r_IC
        
        ftime=np.sum(cb[:,:lexp],axis=0)/np.sum(c[:,:lexp],axis=0)
        np.savetxt(destfold+"/ftime.txt", ftime)
        
        for k in range(len(zero_crossings)): # this is the number of walls we want to track??
        
            wall_flag=0
            
            if k==0:
                topb=(zero_crossings[k])/2
                ## if there is only one boundary detected then k+1  is invalid:/
                ## so in those cases issues can arise
                botb=(zero_crossings[k]+zero_crossings[k+1])/2  
            elif k==len(zero_crossings)-1:
                topb=(zero_crossings[k]+zero_crossings[k-1])/2
                botb=(zero_crossings[k]+l)/2
            else:
                topb=(zero_crossings[k]+zero_crossings[k-1])/2
                botb=(zero_crossings[k]+zero_crossings[k+1])/2 
            #topb=(zero_crossings[k]+min(zero_crossings[k-1],0))/2
            #botb=(zero_crossings[k]+max(zero_crossings.take(k+1, mode='wrap'),l))/2  #this is to get the really chiral boundarries that move past 3/4 and so on
            gap=(botb-topb)/10
            print gap
            
            wall,li,lf,wall_flag=find_wall(li,lf,topb,botb)

            
            
            destfold=fold.replace(path,destpath)   
            if not os.path.exists(destfold): os.makedirs(destfold)          
            
            np.savetxt(destfold+"/wall"+str(k)+".txt", wall)
            np.savetxt(destfold+"/wallFlag"+str(k)+".txt", [wall_flag,])
            lilf=np.array((li,lf,lexp))
            np.savetxt(destfold+"/lilf"+str(k)+".txt", lilf)
    
            boundary_arg=wall[li:lf]        
            xval=np.arange(len(boundary_arg))
            slope,intercept=np.polyfit(xval, boundary_arg,1)
            fig = plt.figure(figsize=(6,6))
            ax = fig.add_subplot(111) 
            ax.plot(xval,boundary_arg,'bo',markersize=5,markeredgecolor=None)      
            ax.plot(xval,slope*xval+intercept,'r-',)           
            destfile_name=destfold+"wall"+str(k)+".png"
            plt.savefig(destfile_name)

        
        #beta_avg+=slope#*xval[-1]/time
        #print slope, beta_avg
        #
        #velocity=(np.sum(c)/l-r_IC)  /time
        
        
        
        run_counter+=1
        


