# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 23:23:19 2015

@author: ashish
"""

# Reads from ca to plot
import numpy as np
#import matplotlib.pyplot as plt
import glob
import os
import math
import matplotlib.pyplot as plt


from scipy.interpolate import griddata

#from scipy.optimize import curve_fit
#doesnt average over runs!! Just the run1 value is used!!!
ms=8
fs=16
lw=2
chirality=[]
slopes=[]


alpha_beta_diff=1 ## alpha-beta which is fixed needs to be given!
alpha_beta_flag=1


#####this script can be used for fstar runs where chsum is used to vary chirality 
path = '/Volumes/Ashish_backup/research data/new_figures/fstar/gdnoise_many runs useful/gdnoise9/'
path = '/Users/ashish/Downloads/fstar_avg/gdnoise9/'
Nvalue='/N400/' #
lvl1foldertype='h'
foldertype='chsum'



#path = '/Volumes/Ashish_backup/research data/new_figures/neutralHdecay/N/gdnoise/ftime_4_avg/'
#Nvalue='/ch0.075/' #
#lvl1foldertype='h'
#foldertype='N'
#

#path = '/projectnb/qedksk/ashishge/data/new_sim/fstar_avg/gdnoise10Steep/'
#Nvalue='/N200/' #
#lvl1foldertype='h'
#foldertype='chsum'



#####,for fstarb runs runs use dfferent script fig_fstarb.py


### for plotting coexistence average data also,( notfstar or fstarb runs)####
#path ='/Volumes/Ashish_backup/research data/to be backed up/new_figures/PAPER FIGURES/oppch_coexistence_fig/moran_avg/try1_shuffle/alpha/'
#Nvalue='/N400'
#lvl1foldertype='h'
#foldertype='chsum'

fractionlist=[]

cdratiolist=[]
fss_conditional=[]
chirality=[]
timelist=[]
#lenpercentlist=[]
ftimelist=[]
frac_steady_state=[]



foldercount=0

for lvl1 in glob.glob(path+'*/'):
    #print lvl1
    frac_steady_state.append([])
    fractionlist.append([])
    cdratiolist.append([])
    fss_conditional.append([])
    #lenpercentlist.append([])
    timelist.append([])
    ftimelist.append([])
    #slopes.append([]) #for each value of lvl1 folders,a new slope array is made with the slopes in it.
    chirality.append([])
    
    for ch in glob.glob(lvl1+'*'):    # for finding chirality value
        
        if os.path.exists(ch+Nvalue):
           
            #FOR FINDING A COMMON LEXP FOR ALL THE RUNS!!
            fold=ch+Nvalue+'/'
            
            
            with open(fold+'fraction_time.txt') as file:
                    ftime = [[float(digit) for digit in line.split()] for line in file]
            ftime=np.array(ftime)
            #ftime is already the mean over the runs unlike cbpops and ca pops but is an array in time

            
            ftimelist[foldercount].append(ftime)#ftime is already the mean over the runs unlike cbpops and ca pops but is an array in time
            #lenpercentlist[foldercount].append(np.mean(lenpercent))
            with open(fold+'fss_runs.txt') as file:
                    fss_runs = [[float(digit) for digit in line.split()] for line in file]
            fss_runs=np.array(fss_runs)
            
            frac_steady_state[foldercount].append(np.mean(ftime[int(len(ftime)*.90):int(len(ftime)*.99)])) # assume at .7*len, steady state is reached
            
            ratiovalue=lvl1.replace(path+lvl1foldertype,'')
            ratiovalue=ratiovalue.replace('/','')
            
            cdratiolist[foldercount].append(float(ratiovalue))
            

            chiralvalue=ch.replace(lvl1+foldertype,'')
            print chiralvalue
            chirality[foldercount].append(float(chiralvalue))
            
            
            ensembles=len(fss_runs)
            if np.sum(fss_runs<0.001) >0 and np.sum(fss_runs<0.001) <ensembles/6: # less than half the runs have fixed accidentally
                fss_conditional[foldercount].append(np.mean(fss_runs[fss_runs>0.001]))
                print "bad fixations at 0"
                print np.mean(fss_runs[fss_runs>0.001])
                print np.mean(ftime[int(len(ftime)*.90):int(len(ftime)*.99)])
            elif np.sum(fss_runs>1-0.001) >0 and np.sum(fss_runs>1-0.001) <ensembles/6:
                fss_conditional[foldercount].append(np.mean(fss_runs[fss_runs<1-0.001])) 
                print "bad fixations at 1"
                print np.mean(fss_runs[fss_runs<1-0.001])
                print np.mean(ftime[int(len(ftime)*.90):int(len(ftime)*.99)])
                
            else:
                fss_conditional[foldercount].append(np.mean(fss_runs)) 
                #print "no bad fixations i think"
                #print np.mean(fss_runs)
                #print np.mean(ftime[int(len(ftime)*.90):int(len(ftime)*.99)])

    
    if foldertype =='ch':
        chirality[foldercount]=(chirality[foldercount]-np.min(np.asarray(chirality[foldercount])))*2 

    #np.savetxt(lvl1+"/frac_steady_state.txt", frac_steady_state[foldercount])
    #
    #np.savetxt(lvl1+"/eff_fitness.txt", np.abs(np.log(1+np.asarray(frac_steady_state[foldercount])-0.5)/lexp)*100)
    #np.savetxt(lvl1+"/chirality.txt", chirality[foldercount])
    

    
    foldercount=foldercount+1

X=np.ravel(np.asarray(cdratiolist))
Y=np.ravel(np.asarray(chirality))

Z=np.ravel(np.asarray(frac_steady_state))

time_unit=36.0
cutoff=1.0 ## we plot ftime till point total time/cutoff


for chidx,chval in enumerate(np.asarray(chirality[0])):
    fig = plt.figure(figsize=(6,6))  
    ax = fig.add_subplot(111)
    limit_line=0.0
    minlength=len(np.asarray(ftimelist)[0][chidx])
    for idx,val in enumerate(np.asarray(cdratiolist).T[chidx]):
        if minlength>len(np.asarray(ftimelist)[idx][chidx]):
            minlength=len(np.asarray(ftimelist)[idx][chidx])-1
            
    yval=[]
    xval=np.arange(minlength)
    for idx,val in enumerate(np.asarray(cdratiolist).T[chidx]):
        xaxis=np.arange(len(np.asarray(ftimelist)[idx][chidx]))*time_unit/len(np.asarray(ftimelist)[idx][chidx])
        ax.plot(xaxis[:int( len(xaxis)/cutoff)],np.asarray(ftimelist)[idx][chidx][:int( len(xaxis)/cutoff)],linestyle='none',marker='o', markeredgecolor='blue',markerfacecolor='none',label="$F_c$="+str(val),markevery=10,markersize=ms)  
        limit_line= limit_line+np.mean(np.asarray(ftimelist)[idx][chidx][int( len(xaxis)/cutoff) - 20:int( len(xaxis)/cutoff)])
        
        #xval.append(np.arange(minlength))## not converted to arbitrary time units, call runs chopped off to same length
        yval.append(np.ravel(np.asarray(ftimelist)[idx][chidx])[:minlength])
        
        #print len(np.asarray(ftimelist)[idx][chidx])*1.0
    #plt.legend( loc=1, borderaxespad=0.,prop={'size':fs,'weight':'bold'})
    #ax.plot([0, time_unit/cutoff], [0.5, 0.50],'black', linestyle='--',linewidth =lw)
    plt.xticks(np.arange(0, time_unit/cutoff +1, 3.0))
    #ax.plot([0, time_unit/cutoff], [limit_line/(1+idx), limit_line/(1+idx)] ,'black', linestyle='--',linewidth =lw)
    ax.plot([0, time_unit/cutoff], [1.0, 1.0] ,'black', linestyle='--',linewidth =lw)
    ax.set_ylim(0.0,1.1)
    plt.xlabel('time (a.u.)',fontsize = fs,fontweight='bold')
    plt.ylabel('fraction of left-handed strain',fontsize = fs, fontweight='bold')
    plt.xticks([0,2,4,6])
    plt.tick_params(axis='both', which='major', labelsize=fs)   
    plt.tight_layout() 
    plt.savefig(path+'ftime_vs_IC_ch'+str(chval)+'.png')
    
    
    np.savetxt(path+'time_ch'+str(chval)+'.txt',np.array(xval))
    np.savetxt(path+'fraction_ch'+str(chval)+'.txt',np.asarray(yval))
for chidx,chval in enumerate(np.asarray(chirality[0])):
    fig = plt.figure(figsize=(6,6))  
    ax = fig.add_subplot(111)
    for idx,val in enumerate(np.asarray(cdratiolist).T[chidx]):
        ax.semilogx(np.asarray(ftimelist)[idx][chidx],linestyle='none',marker='o',markeredgecolor='blue',markerfacecolor='none',label="$F_c$="+str(val),markevery=5)  
    #plt.legend( loc=1, borderaxespad=0.,prop={'size':fs,'weight':'bold'})
    ax.set_ylim(0.0,1.02)
    #ax.semilogx([0, len(np.asarray(ftimelist)[idx][chidx])], [0.5, 0.50],'black', linestyle='--',linewidth =lw)
    plt.xlabel('time (a.u.)',fontsize = fs,fontweight='bold')
    plt.ylabel('mutant frequency',fontsize = fs, fontweight='bold')
    plt.tick_params(axis='both', which='major', labelsize=fs)   
    plt.tight_layout() 
    plt.savefig(path+'flogtime_vs_IC_ch'+str(chval)+'.png')

#### DONT REMEMBER WHAT THIS STUFF IS ABOUT, SO COMMENTED
#if alpha_beta_flag==1:
#    xaxis=np.array(chirality[0])-1
#    xaxis=(alpha_beta_diff+xaxis)*0.5/alpha_beta_diff
#else:
#    xaxis=np.array(np.asarray(chirality[0])-(chirality[0][0]+chirality[0][len(chirality[0])-1])/2)
if foldertype =='chsum':
        xaxis=1.0 - np.asarray(chirality[0])/2



frac_steady_state=np.array(frac_steady_state)
fig = plt.figure(figsize=(6,6))  
ax = fig.add_subplot(111)
plt.plot(xaxis,np.mean(frac_steady_state, axis=0),'o')
ax.axvspan(-0.02, 0.5, alpha=0.25, color='red')
ax.axvspan(0.5, 1.02, alpha=0.25, color='green')
ax.set_ylim(-0.02,1.02)
ax.set_xlim(-0.02,1.02)
ax.plot([-1, 1], [0.5, 0.50],'black', linestyle='--',linewidth =lw)
#ax.plot([-1,1], [1.0, 1.0],'black', linestyle='-',linewidth =lw)
#ax.plot([-1,1], [0.0, 0.0],'black', linestyle='-',linewidth =lw)
#plt.xlabel('chirality difference, $\\frac{|v_a|-|v_b|}{|v_a-v_b|}$',fontsize = fs,fontweight='bold')
plt.xlabel('relative chirality $\\mathbf{f^* }$',fontsize = fs, fontweight="bold")
plt.ylabel('steady state fraction, $ \\mathbf{\\bar{f}_{\\mathrm{\\mathbf{eq}}}}}$',fontsize = fs, fontweight='bold') 
plt.tick_params(axis='both', which='major', labelsize=fs)
plt.tight_layout()
#plt.savefig(path+'poster_fss_vs_chirality.png')
plt.savefig(path+'fbar_fstar.pdf')




fig = plt.figure(figsize=(9,9))  
ax = fig.add_subplot(111)
plt.plot(xaxis,np.mean(fss_conditional, axis=0),'o')
#ax.fill_between(0, 0.0, 1.0, facecolor='yellow', alpha=0.5,label='1 sigma range')
ax.axvspan(np.min(xaxis)*1.1, 0.5, alpha=0.25, color='red')
ax.axvspan(0.5, np.max(xaxis)*1.1, alpha=0.25, color='green')
ax.set_ylim(-0.01,1.01)
ax.set_xlim(np.min(xaxis)*1.1,np.max(xaxis)*1.1)
ax.plot([-1, 1], [0.5, 0.50],'black', linestyle='--',linewidth =lw)
#ax.plot([-1, 1], [1.0, 1.0],'black', linestyle='-',linewidth =lw)
#ax.plot([-1,1], [0.0, 0.0],'black', linestyle='-',linewidth =lw)
#plt.xlabel('chirality difference, $\\frac{\\chi_a}{\\chi_a-\\chi_b}$',fontsize = fs,fontweight='bold')
plt.xlabel('$\\mathbf{f^*}$, $\\mathbf{ \\frac{A^{(1)} }{A^{(1)}-A^{(2)}}  }$',fontsize = fs)
plt.ylabel('steady state fraction, $\\bar{f}$',fontsize = fs, fontweight='bold') 
plt.tick_params(axis='both', which='major', labelsize=fs) 
plt.tight_layout()
plt.savefig(path+'fbar_fstar_conditional_vs_chirality.pdf')

#np.savetxt(path+"/scaled_ch.txt", xaxis)
#np.savetxt(path+"/fstar.txt", np.mean(frac_steady_state, axis=0))
print path

np.savetxt(path+'fbar_fstar.txt',np.mean(frac_steady_state, axis=0))
np.savetxt(path+'fbar_fstar_conditional.txt',np.mean(fss_conditional, axis=0))

np.savetxt(path+'fbar_fstar_xaxis.txt',xaxis)


plt.close("all")

