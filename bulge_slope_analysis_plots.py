#!/usr/bin/env python

import numpy as np
import glob
import os
import math
from sys import exit
#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm


path='/Users/ashish/Downloads/bulge_slope_from_track_height/try3_slope_analysis/'
destpath='/Users/ashish/Downloads/bulge_slope_from_track_height/try3_slope_analysis_plots/'


  
if not os.path.exists(destpath): os.makedirs(destpath) 


lvl1_prefix="chsum"
lvl2_prefix='h'

chlist=[]
hlist=[]
lvl1_foldercounter=0  

vParallel_list=[]
alphaD_list=[]
vParallel_stderr_list=[]
alphaD_stderr_list=[]         
                           
for lvl1 in glob.glob(path+'*/'):  ## chsum value
    chsumval=lvl1.replace(path,'')            
    chsumval=chsumval.replace(lvl1_prefix,'')
    chsumval=float(chsumval.replace('/',''))
    print "chsumval is",chsumval        
    
        
            
    slope1_list=[] 
    slope2_list=[] 
    std_err1_list=[] 
    std_err2_list=[] 
    r_value1_list=[] 
    r_value2_list=[] 
    tlist_list=[]
    
    chlist.append([] )
    hlist.append([] )              
    lvl2_foldercounter=0                                             
    for lvl2 in glob.glob(lvl1+'/*/'):
        
        hval=lvl2.replace(lvl1,'')            
        hval=hval.replace(lvl2_prefix,'')
        hval=float(hval.replace('/',''))
        print "hval is",hval        
        
        
        
        slope1_list.append([] )
        slope2_list.append([] ) 
        std_err1_list.append([] ) 
        std_err2_list.append([] ) 
        r_value1_list.append([] ) 
        r_value2_list.append([] ) 
        tlist_list.append([] ) 
        
        
        run_counter=0
        for fold in glob.glob(lvl2+'/*/'):
            
            
            
            slope1=np.loadtxt(fold+"slope1.txt")
            slope2=np.loadtxt(fold+"slope2.txt")
            std_err1=np.loadtxt(fold+"std_err1.txt")
            std_err2=np.loadtxt(fold+"std_err2.txt")           
            r_value1=np.loadtxt(fold+"R1.txt")
            r_value2=np.loadtxt(fold+"R2.txt")
            tlist=np.loadtxt(fold+"tlist.txt")
            
            
            
            slope1_list[lvl2_foldercounter].append(slope1 )
            slope2_list[lvl2_foldercounter].append(slope2 ) 
            std_err1_list[lvl2_foldercounter].append(std_err1 ) 
            std_err2_list[lvl2_foldercounter].append(std_err2 ) 
            r_value1_list[lvl2_foldercounter].append(r_value1 ) 
            r_value2_list[lvl2_foldercounter].append(r_value2 ) 
            tlist_list[lvl2_foldercounter].append(tlist)
            run_counter+=1
            
        chlist[lvl1_foldercounter].append(chsumval)
        hlist[lvl1_foldercounter].append(hval)
        lvl2_foldercounter+=1
        
    normalized_arr=(np.array(hlist)-np.min(hlist) )/ (np.max(hlist)-np.min(hlist))
    line_colors = cm.rainbow(normalized_arr)   
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)     
    for i in range(lvl2_foldercounter)   :
        for j in range(run_counter):
            if j==0: ## label is required
                plt.errorbar( tlist_list[i][j],slope1_list[i][j],label=lvl2_prefix+"="+str(hlist[lvl1_foldercounter][i]),yerr=std_err1_list[i][j],color=line_colors[lvl1_foldercounter][i] )
            else: ## label is not required
                plt.errorbar( tlist_list[i][j],slope1_list[i][j],yerr=std_err1_list[i][j],color=line_colors[lvl1_foldercounter][i] )
    plt.legend(loc="best")
    destfile_name=destpath+"slope1_"+lvl1_prefix+str(chlist[lvl1_foldercounter][i])+".png"
    plt.savefig(destfile_name)    
    
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)     
    for i in range(lvl2_foldercounter)   :
        for j in range(run_counter):
            if j==0: ## label is required
                plt.errorbar( tlist_list[i][j],slope2_list[i][j],label=lvl2_prefix+"="+str(hlist[lvl1_foldercounter][i]),yerr=std_err2_list[i][j],color=line_colors[lvl1_foldercounter][i] )
            else: ## label is not required
                plt.errorbar( tlist_list[i][j],slope2_list[i][j],yerr=std_err2_list[i][j],color=line_colors[lvl1_foldercounter][i] )
    plt.legend(loc="best")
    destfile_name=destpath+"slope2_"+lvl1_prefix+str(chlist[lvl1_foldercounter][i])+".png"
    plt.savefig(destfile_name)
    

    slopes1=[]
    slopes2=[]
    for i in range(lvl2_foldercounter)   :
        for j in range(run_counter):
          slopes1.append(  slope1_list[i][j][5:])
          slopes2.append(  slope2_list[i][j][5:])
            
    slopes1=np.ravel(slopes1) ## we only use the fits towards the end!
    slopes2=np.ravel(slopes2)
    
    
    
    ### (slope1+slope2 )/2 = vparallel/v0
    ### (slope1-slope2 ) = alpha/D
    ### beta=vparallel/(0.5-f*)
    
    slope1_mean=np.mean(slopes1)
    slope2_mean=np.mean(slopes2)
    slope1_var=np.var(slopes1)
    slope2_var=np.var(slopes2)
    
    alphaD=slope1_mean-slope2_mean
    alphaD_stderr=np.sqrt(  slope1_var/run_counter + slope2_var/run_counter) ## we assume that pnly the different replicas are taken into account to divide by,
                                                                            ## the fact that we used multiple values from each replica is ignored!
    
    
    vParallel=(slope1_mean+slope2_mean )/2
    vParallel_stderr=np.sqrt(  0.5*slope1_var/run_counter + 0.5*slope2_var/run_counter) ## 0.5 cos vparallel formula has a factor of 1/2
    
    
    
    vParallel_list.append(vParallel)
    vParallel_stderr_list.append(vParallel_stderr)
    alphaD_list.append(alphaD)
    alphaD_stderr_list.append(alphaD_stderr)
    
    
    print np.shape(slope1_list),np.shape(slope1_list[0]),np.shape(slope1_list[0][0])   
    lvl1_foldercounter+=1
    
    
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)            
plt.errorbar(np.array(chlist)[:,0], vParallel_list,yerr=vParallel_stderr_list,fmt='o' )
destfile_name=destpath+"vParallel.png"
plt.savefig(destfile_name)  


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)            
plt.errorbar(np.array(chlist)[:,0], alphaD_list,yerr=alphaD_stderr_list ,fmt='o')
destfile_name=destpath+"alphaD.png"
plt.savefig(destfile_name) 


charr=np.array(chlist)[:,0]

fstar=charr[np.where(charr!=1.)]/2  
beta= np.array(vParallel_list)[np.where(charr!=1.)] /(0.5-fstar)


xaxis=charr[np.where(charr!=1.)]
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)            
plt.errorbar(xaxis,beta,yerr=np.array(2*vParallel_stderr_list)[np.where(charr!=1.)] ,fmt='o')
destfile_name=destpath+"betaV.png"
plt.savefig(destfile_name) 


fbar=0.5+ (  beta/ np.array(alphaD_list)[np.where(charr!=1.) ]  ) *  (fstar-0.5)
  
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)            
plt.errorbar(xaxis,fbar,yerr= np.array(alphaD_stderr_list +vParallel_stderr_list*2)[np.where(charr!=1.)] ,fmt='o')
destfile_name=destpath+"Fbar.png"
plt.savefig(destfile_name)   
  
  
np.savetxt(destpath+"vParallel.txt", vParallel_list)
np.savetxt(destpath+"vParallel_stderr.txt", vParallel_stderr_list)
np.savetxt(destpath+"alphaD.txt", alphaD_list)
np.savetxt(destpath+"alphaD_stderr.txt", alphaD_stderr_list)
np.savetxt(destpath+"chlist.txt", chlist)
