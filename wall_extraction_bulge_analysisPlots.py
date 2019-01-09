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
from scipy import stats


path='/Users/ashish/Downloads/bulge_slope_from_track_height/try3_wallExtract/'
destpath='/Users/ashish/Downloads/bulge_slope_from_track_height/try3_wallExtract_plots/'
destpath2='/Users/ashish/Downloads/bulge_slope_from_track_height/try3_wallExtract_ftimePlots/'


  
if not os.path.exists(destpath): os.makedirs(destpath) 
if not os.path.exists(destpath2): os.makedirs(destpath2) 


def plot_wall_and_fit(xval,yval,fit,filename):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111) 
    ax.plot(xval,yval,'bo',markersize=5,markeredgecolor=None)      
    ax.plot(xval,fit,'r-',)           
    plt.savefig(filename)



lvl1_prefix="chsum"
lvl2_prefix='h'

chlist=[]
hlist=[]
fbar_list=[]
fbar_stderr_list=[]

vParallel_list=[]
vParallel_stderr=[]

lvl1_foldercounter=0           
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
    wall_list=[]
    ftime_list=[]
    
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
        wall_list.append([] ) 
        ftime_list.append([] ) 
        
        
        run_counter=0
        for fold in glob.glob(lvl2+'/*/'):
            
            
            ftime=np.loadtxt(fold+"ftime.txt")
            
            wall1=np.loadtxt(fold+"wall1.txt")
            [li,lf,lexp]=np.loadtxt(fold+"lilf1.txt")
            destfold=fold.replace(path,destpath)   
            if not os.path.exists(destfold): os.makedirs(destfold)       
            
            xval=np.arange(li,lf)
            slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress( xval,wall1  )          
            fit= slope1*xval+intercept1            
            fit_file_name=destfold+"wall1.png"
            plot_wall_and_fit(xval,wall1,fit,fit_file_name)
            
            slope1_list[lvl2_foldercounter].append(slope1 )
            std_err1_list[lvl2_foldercounter].append(std_err1 ) 
            r_value1_list[lvl2_foldercounter].append(r_value1 ) 
            wall_list[lvl2_foldercounter].append(wall1 )
            ftime_list[lvl2_foldercounter].append(ftime )
            
            
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
                plt.plot(np.arange(li,li+len(wall_list[i][j]) ), wall_list[i][j]-wall_list[i][j][0],label=lvl2_prefix+"="+str(hlist[lvl1_foldercounter][i]),color=line_colors[lvl1_foldercounter][i] )
            else: ## label is not required
                plt.plot(np.arange(li,li+len(wall_list[i][j]) ), wall_list[i][j]-wall_list[i][j][0],color=line_colors[lvl1_foldercounter][i] )
    plt.legend(loc="best")
    destfile_name=destpath+"wall1_"+lvl1_prefix+str(chlist[lvl1_foldercounter][i])+".png"
    plt.savefig(destfile_name)    
    
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)     
    for i in range(lvl2_foldercounter)   :        
        plt.errorbar(np.arange(run_counter), slope1_list[i],label=lvl2_prefix+"="+str(hlist[lvl1_foldercounter][i]),yerr=std_err1_list[i],color=line_colors[lvl1_foldercounter][i] )
    plt.legend(loc="best")
    destfile_name=destpath+"slope1_"+lvl1_prefix+str(chlist[lvl1_foldercounter][i])+".png"
    plt.savefig(destfile_name)  
    
    
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)     
    for i in range(lvl2_foldercounter)   :
        for j in range(run_counter):
            if j==0: ## label is required
                plt.plot(np.arange(len(ftime_list[i][j]) ), ftime_list[i][j],label=lvl2_prefix+"="+str(hlist[lvl1_foldercounter][i]),color=line_colors[lvl1_foldercounter][i] )
            else: ## label is not required
                plt.plot(np.arange(len(ftime_list[i][j]) ), ftime_list[i][j],color=line_colors[lvl1_foldercounter][i] )
    plt.legend(loc="best")
    destfile_name=destpath2+"ftime_"+lvl1_prefix+str(chlist[lvl1_foldercounter][i])+".png"
    plt.savefig(destfile_name)    
    
    
    slopes=np.ravel(  slope1_list)  
    vParallel_list.append(np.mean(slopes))
    print len(slopes)
    stderror_slopes=np.sqrt( np.var(slopes)/len(slopes) )
    vParallel_stderr.append( stderror_slopes)
    
    
    ftimes=[]
    for i in range(lvl2_foldercounter):
        for j in range(run_counter):
          ftimes.append(  np.mean(ftime_list[i][j][  int(len(ftime_list[i][j])*0.85) :] )  )
         
    
    fbar_list.append(np.mean(ftimes))
    fbar_stderr_list.append(  np.sqrt( np.var(ftimes)/len(ftimes) )  )
    
    lvl1_foldercounter+=1
    
    
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)            
plt.errorbar(np.array(chlist)[:,0], vParallel_list,yerr=vParallel_stderr,fmt='o' )
destfile_name=destpath+"vParallel.png"
plt.savefig(destfile_name)  
np.savetxt(destpath+"vParallel.txt", vParallel_list)
np.savetxt(destpath+"vParallel_stderr.txt", vParallel_stderr)
np.savetxt(destpath+"chlist.txt", chlist)
  
    
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)            
plt.errorbar(np.array(chlist)[:,0], fbar_list,yerr=fbar_stderr_list,fmt='o' )
destfile_name=destpath2+"Fbar.png"
plt.savefig(destfile_name)  


np.savetxt(destpath+"fbar_stderr.txt", fbar_stderr_list)
np.savetxt(destpath+"fbar.txt", fbar_list)
    
    
    
plt.close("all")