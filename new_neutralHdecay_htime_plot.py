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
import matplotlib.ticker as mtick

from scipy.interpolate import griddata

#from scipy.optimize import curve_fit
#doesnt average over runs!! Just the run1 value is used!!!

fs=20

lw=3
chirality=[]

slopes=[]
slopesloglog=[]
path = '/Volumes/Ashish_backup/research data/new_figures/neutralHdecay/N/gdnoise/try2_rescale_avg/'

Nvalue='ch0.75/'
lvl1foldertype='a'
foldertype='N'


het_spatial_flag=1



fs=20
sfs=12
lw=2
ms=10



import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = lw
mpl.rcParams['font.size'] = sfs
mpl.rcParams['axes.labelweight'] = "bold"
mpl.rcParams['axes.labelsize'] = fs



smlist=[]
#velalist=[]
#velblist=[]
calist=[]
cblist=[]
chirality=[]
timelist=[]

het_all_list=[] #het time of all folders together
hetlist=[]
het_all_notlog_list=[]
het_ss_list=[]
het_spatial_ss_list=[]
lenpercentlist=[]

foldercount=0
for lvl1 in glob.glob(path+'*/'):
    print lvl1
    
    #slope1list.append([])
    #slope2list.append([])
    #velalist.append([])
    #velblist.append([])
    smlist.append([])
    hetlist.append([])
    calist.append([])
    cblist.append([])
    lenpercentlist.append([])
    timelist.append([])
    slopes.append([]) #for each value of lvl1 folders,a new slope array is made with the slopes in it.
    slopesloglog.append([])
    chirality.append([])
    het_ss_list.append([])
    het_spatial_ss_list.append([])
    het_all_list.append([])
    het_all_notlog_list.append([])
    #velocities.append([])
    for ch in glob.glob(lvl1+'*/'):    # for finding chirality value
        if os.path.exists(ch+Nvalue):
            print ch
            #FOR FINDING A COMMON LEXP FOR ALL THE RUNS!!
            fold=ch+Nvalue+''
            with open(fold+'capops.txt') as file:
                ca = [[float(digit) for digit in line.split()] for line in file]
            
            ca=np.array(ca)
            #print ca
            with open(fold+'cbpops.txt') as file:       
                cb = [[float(digit) for digit in line.split()] for line in file]
            with open(fold+'time.txt') as file:
                    time = [[float(digit) for digit in line.split()] for line in file]
            time=np.array(time)
            #with open(fold+'lenpercent.txt') as file:
            #        lenpercent = [[float(digit) for digit in line.split()] for line in file]
            #lenpercent=np.array(lenpercent)
            #with open(fold+'heterozygosity.txt') as file:
            #    het = [[float(digit) for digit in line.split()] for line in file]         
            #het=np.array(het)
            #with open(fold+'het_exp_avg.txt') as file:
            #    het_exp = [[float(digit) for digit in line.split()] for line in file]
            
            if het_spatial_flag==1:
                with open(fold+'HET_spatial_time.txt') as file:
                   het_spatial_time= [[float(digit) for digit in line.split()] for line in file]
                het_spatial_time=np.ravel(het_spatial_time)
                
            with open(fold+'HET_time.txt') as file:
                HET_time = [[float(digit) for digit in line.split()] for line in file]
            HET_time=np.ravel(HET_time)
            if os.path.isfile(fold+'sumclist.txt'):
                 with open(fold+'sumclist.txt') as file:
                    sumc = [[float(digit) for digit in line.split()] for line in file]
                 sumc=np.ravel(np.asarray(sumc))
                 fig = plt.figure(figsize=(6,6))
                 ax = fig.add_subplot(111)
                 ax.set_title('sumc')
                 plt.plot(sumc,'o')
                 plt.savefig(fold+'sumc.png')
            
            c=cb+ca
            c=np.array(c)  
                        
            #fig = plt.figure(figsize=(6,6))  
            #ax = fig.add_subplot(111)
            #plt.plot(het_exp,'o') 
            #plt.savefig(fold+'het_exp_avg.png')
            #
            #fig = plt.figure(figsize=(6,6))  
            #ax = fig.add_subplot(111)
            #plt.semilogy(het_exp,'o') 
            #plt.savefig(fold+'HET_exp_semilog.png')
            #
            #fig = plt.figure(figsize=(6,6))  
            #ax = fig.add_subplot(111)
            #plt.loglog(het_exp,'o') 
            #plt.savefig(fold+'HET_exp_loglog.png')

            
            fig = plt.figure(figsize=(6,6))  
            ax = fig.add_subplot(111)
            plt.plot(HET_time,'o') 
            plt.savefig(fold+'HET_time.png')
            
            
            #lexp=len(HET_time)
            zero_element=np.where(np.asarray(HET_time[30:])==0) # starts from 30 in case 1 boundary, where het at beginning is 0
            if len(zero_element[0]) :
                lexp= zero_element[0][0]
            else :
                lexp=len(HET_time)
            
            het_all_list[foldercount].append(np.log10(np.ravel(np.asarray(HET_time[:lexp]))))
            het_all_notlog_list[foldercount].append(np.ravel(np.asarray(HET_time[:lexp])))
            SS=2
            lin=400
            #lin=lexp/SS+1
            #if lin<500:
            #    lin=500
            lst=lexp-lin
            xdata=np.arange(lin,lst+lin)
            ydata=np.asarray(HET_time[lin:lexp])
            ydata = ydata.reshape((lst))
            
            
            logx=np.log(xdata)
            logy=np.log(ydata)


            fig = plt.figure(figsize=(6,6))  
            ax = fig.add_subplot(111)
            plt.semilogy(HET_time,'o') 
            plt.savefig(fold+'HET_time_semilog.png')
            
            
            if het_spatial_flag==1:
                if len(np.where(het_spatial_time>0)[0])>10:
                    fig = plt.figure(figsize=(6,6))  
                    ax = fig.add_subplot(111)
                    plt.semilogy(het_spatial_time,'o') 
                    plt.savefig(fold+'HET__spatial_semilog.png')
            
            z=np.polyfit(xdata, logy, 1)
            p=np.poly1d(z)
            slopes[foldercount].append(z[0])
                            
            k=np.zeros(lst)
        
            for i in range(0, lst):
                k[i]=p(xdata[i]) 
                
            fig = plt.figure(figsize=(6,6))  
            ax = fig.add_subplot(111)
            #plt.plot(xdata,logy,'o',xdata,k)
            plt.semilogy(xdata,ydata,'o',xdata,np.exp(k), linewidth=lw) 
            #plt.semilogy(xdata,k, linewidth=lw) 
            plt.xlabel("time", fontsize=fs, fontweight ='bold')
            plt.ylabel("heterozygosity",fontsize=fs, fontweight ='bold')
            plt.savefig(fold+'HET_time_semilog_fit.png')
            
            
              
                                          
            fig = plt.figure(figsize=(6,6))  
            ax = fig.add_subplot(111)
            plt.loglog(HET_time,'o') 
            plt.savefig(fold+'HET_time_loglog.png')
            
            
            z=np.polyfit(logx, logy, 1)
            p=np.poly1d(z)
            slopesloglog[foldercount].append(z[0])        
            k=np.zeros(lst)
        
            for i in range(0, lst):
                k[i]=p(logx[i]) 
            
            fig = plt.figure(figsize=(6,6))  
            ax = fig.add_subplot(111)
            #plt.plot(xdata,logy,'o',xdata,k)
            plt.loglog(xdata,ydata,'o',xdata,np.exp(k), linewidth=lw) 
            #plt.semilogy(xdata,k, linewidth=lw) 
            plt.xlabel("time", fontsize=fs, fontweight ='bold')
            plt.ylabel("heterozygosity",fontsize=fs, fontweight ='bold')
            plt.savefig(fold+'HET_time_loglog_fit.png')
            
            
            #fig = plt.figure(figsize=(6,6))  
            #ax = fig.add_subplot(111)
            #plt.semilogy(xdata,HET_time,'o') 
            #plt.savefig(fold+'HET_time_loglog.png')
            
            fig = plt.figure(figsize=(6,6))  
            ax = fig.add_subplot(111)
            plt.plot(cb/c,'*') 
            plt.savefig(fold+'cb_fraction.png')
            
            #lenpercentlist[foldercount].append(np.mean(lenpercent))            
            # the velocities are found by meanfield simulations so no np.mean reqd
            #velalist[foldercount].append(vela)
            #velblist[foldercount].append(velb)            
            #not sure whether to average slopes or tanthetas so no averaging done. only meanfield really makes sense
            #slope1list[foldercount].append(slope1)
            #slope2list[foldercount].append(slope2)
            
            #averaging over identical runs if noise exists
            
            het_ss_list[foldercount].append(np.mean(HET_time[len(HET_time)-4000:len(HET_time)-1000]))
            if het_spatial_flag==1:
                het_spatial_ss_list[foldercount].append(np.mean(het_spatial_time[len(het_spatial_time)-400:len(het_spatial_time)-100]))
            
            #hetlist[foldercount].append(np.mean(het))
            cblist[foldercount].append(np.mean(cb))
            calist[foldercount].append(np.mean(ca))
            timelist[foldercount].append(np.mean(time))
            chiralvalue=ch.replace(lvl1+foldertype,'')
            chiralvalue=chiralvalue.replace('/','')
            chirality[foldercount].append(float(chiralvalue))
            sm=lvl1.replace(path+lvl1foldertype,'')
            sm=sm.replace('/','')
            smlist[foldercount].append(float(sm))
    
    fraction=np.asarray(cblist[foldercount])/(np.asarray(calist[foldercount])+np.asarray(cblist[foldercount]))
    arealvelocity=(np.asarray(calist[foldercount])+np.asarray(cblist[foldercount]))/np.asarray(timelist[foldercount])
    arealvelocityb=(np.asarray(cblist[foldercount]))/np.asarray(timelist[foldercount])
    
    if foldertype == 'ch':
        chirality[foldercount]=(chirality[foldercount]-np.min(np.asarray(chirality[foldercount])))*2
     
     
     
     # to get actual al -ar
    
    #fig = plt.figure(figsize=(6,6))  
    #ax = fig.add_subplot(111)
    #plt.plot(chirality[foldercount],np.abs(slope1list[foldercount])+np.abs(slope2list[foldercount]),'o') 
    #plt.savefig(lvl1+'tangents vs chirality.png')
    #
    #fig = plt.figure(figsize=(6,6))  
    #ax = fig.add_subplot(111)
    #plt.plot(chirality[foldercount],np.arctan(np.abs(slope1list[foldercount]))+np.arctan(np.abs(slope1list[foldercount])),'o') 
    #plt.savefig(lvl1+'thetas vs chirality.png')
    

    
    fig = plt.figure(figsize=(6,6))  
    ax = fig.add_subplot(111)
    plt.plot(chirality[foldercount],np.abs(slopes[foldercount]),'o') 
    plt.xlabel("microscopic chirality, $a_l-a_r$", fontsize=fs, fontweight='bold')
    plt.ylabel("decay timescale$\\tau$", fontsize=fs, fontweight='bold')
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    plt.tight_layout()
    plt.savefig(lvl1+'Hetdecay_vs_chirality.png')
    
    
    
    fig = plt.figure(figsize=(6,6))  
    ax = fig.add_subplot(111)
    plt.plot(chirality[foldercount],np.abs(slopesloglog[foldercount]),'o') 
    plt.xlabel("microscopic chirality, $a_l-a_r$", fontsize=fs, fontweight='bold')
    plt.ylabel("decay timescale$\\tau$", fontsize=fs, fontweight='bold')
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    plt.tight_layout()
    plt.savefig(lvl1+'Hetdecay_loglog_vs_chirality.png')
    
    
    fig = plt.figure(figsize=(6,6))  
    ax = fig.add_subplot(111)
    plt.plot(chirality[foldercount],het_ss_list[foldercount],'o') 
    #plt.xlabel("microscopic chirality, $a_l-a_r$", fontsize=fs, fontweight='bold')
    plt.xlabel("constant diffusion, $a_0$", fontsize=fs, fontweight='bold')
    plt.ylabel("final heterozygosity", fontsize=fs, fontweight='bold')
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    plt.tight_layout()
    plt.savefig(lvl1+'Het_final_vs_chirality.png')
    
    np.savetxt(lvl1+"het_ss.txt",het_ss_list[foldercount])
    np.savetxt(lvl1+"chirality.txt",chirality[foldercount])
    
    fig = plt.figure(figsize=(6,6))  
    ax = fig.add_subplot(111)
    plt.plot(chirality[foldercount],het_ss_list[foldercount],'o')
    plt.xlabel('microscopic chirality $a_l-a_r$', fontsize=fs, fontweight='bold')
    plt.ylabel('productivity, mixing',fontsize=fs, fontweight='bold') 
    plt.tight_layout()
    plt.savefig(lvl1+'het_vs_chirality.png')
    
    if het_spatial_flag==1:
        fig = plt.figure(figsize=(6,6))  
        ax = fig.add_subplot(111)
        plt.plot(chirality[foldercount],het_spatial_ss_list[foldercount],'o')
        plt.xlabel('microscopic chirality $a_l-a_r$', fontsize=fs, fontweight='bold')
        plt.ylabel('productivity, mixing',fontsize=fs, fontweight='bold') 
        plt.tight_layout()
        plt.savefig(lvl1+'hetspatial_vs_chirality.png')
    
    #fig = plt.figure(figsize=(6,6))  
    #ax = fig.add_subplot(111)
    #plt.plot(chirality[foldercount],fraction,'o') 
    #plt.savefig(lvl1+'fraction_vs_chirality.png')
    #
    #fig = plt.figure(figsize=(6,6))  
    #ax = fig.add_subplot(111)
    #plt.plot(chirality[foldercount],hetlist[foldercount],'o')
    #plt.xlabel('microscopic chirality $a_l-a_r$', fontsize=fs, fontweight='bold')
    #plt.ylabel('productivity, mixing',fontsize=fs, fontweight='bold') 
    #plt.tight_layout()
    #plt.savefig(lvl1+'het_vs_chirality.png')
    #
    #fig = plt.figure(figsize=(6,6))  
    #ax = fig.add_subplot(111)
    #plt.plot(chirality[foldercount],cblist[foldercount],'o') 
    #plt.savefig(lvl1+'cb population.png')
    
    #fig = plt.figure(figsize=(6,6))  
    #ax = fig.add_subplot(111)
    #plt.plot(chirality[foldercount],velblist[foldercount],'o') 
    #plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    #plt.savefig(lvl1+'vel_vs_chirality.png')
    #
    #
    #fig = plt.figure(figsize=(6,6))  
    #ax = fig.add_subplot(111)
    #plt.plot(velblist[foldercount],cblist[foldercount],'o') 
    #plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    #plt.savefig(lvl1+'cb vs velb.png')
    #
    #fig = plt.figure(figsize=(6,6))  
    #ax = fig.add_subplot(111)
    #plt.plot(chirality[foldercount],velblist[foldercount]-np.mean(velalist[foldercount]),'o') 
    #plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    #plt.savefig(lvl1+'veldiff_vs_chirality.png')
    #
    #fig = plt.figure(figsize=(6,6))  
    #ax = fig.add_subplot(111)
    #plt.plot(chirality[foldercount],arealvelocity,'o') 
    #plt.savefig(lvl1+'areal_velocity_vs_chirality.png')
    #
    #fig = plt.figure(figsize=(6,6))  
    #ax = fig.add_subplot(111)
    #plt.plot(chirality[foldercount],arealvelocityb,'o') 
    #plt.savefig(lvl1+'areal_velocity_vs_chiralityb.png')
    #
    #fig = plt.figure(figsize=(6,6))  
    #ax = fig.add_subplot(111)
    #plt.plot(chirality[foldercount],lenpercentlist[foldercount],'o') 
    #plt.savefig(lvl1+'lenpercent_vs_chirality.png')
    
    foldercount=foldercount+1

time_unit=36.0


fig = plt.figure(figsize=(8,8))  
ax = fig.add_subplot(111)

for idx,val in enumerate(chirality[0]):
    xaxis=np.arange(len(het_all_list[0][idx]))*time_unit/len(het_all_list[0][idx])
    
    if foldertype =='ch':
        ax.plot(xaxis,het_all_list[0][idx],'o',label="$a_l-a_r$ ="+str(val), markevery=200)
    else: 
        ax.plot(xaxis,het_all_list[0][idx],'o',label="$a_0$ ="+str(val), markevery=200)
plt.legend(bbox_to_anchor=(0.0, 0.0), loc=3, borderaxespad=0.1)
#plt.xlabel('Diffusion constant, $a_0$',fontsize = fs,fontweight='bold')


plt.xlabel('time (a.u.)',fontsize = fs*8/6,fontweight='bold')

plt.ylabel('log Heterozygosity, $\\log H $',fontsize = fs*8/6,fontweight='bold')
plt.tick_params(axis='both', which='major', labelsize=fs*8/6, pad =0.04) 
plt.tight_layout()
plt.savefig(path+'logH_all.png')


fig = plt.figure(figsize=(8,8))  
ax = fig.add_subplot(111)

for idx,val in enumerate(chirality[0]):
    xaxis=np.arange(len(het_all_list[0][idx]))*time_unit/len(het_all_list[0][idx])
    
    if foldertype =='ch':
        ax.plot(xaxis,np.power(10,het_all_list[0][idx]),'o',label="$a_l-a_r$ ="+str(val), markevery=200)
    else: 
        ax.plot(xaxis,np.power(10,het_all_list[0][idx]),'o',label="$a_0$ ="+str(val), markevery=200)
plt.legend(bbox_to_anchor=(1.0, 1.0), loc=1, borderaxespad=0.1)
#plt.xlabel('Diffusion constant, $a_0$',fontsize = fs,fontweight='bold')

plt.ylim(0.0,0.25)
plt.xlabel('time (a.u.)',fontsize = fs*8/6,fontweight='bold')

plt.ylabel('Heterozygosity, $H$',fontsize = fs*8/6,fontweight='bold')
plt.tick_params(axis='both', which='major', labelsize=fs*8/6, pad =0.04) 
plt.tight_layout()
plt.savefig(path+'H_all.png')



#fig = plt.figure(figsize=(8,8))  
#ax = fig.add_subplot(111)
#
#for idx,val in enumerate(chirality[0]):
#    xaxis=np.arange(len(het_all_notlog_list[0][idx]))*time_unit/len(het_all_notlog_list[0][idx])
#    
#    if foldertype =='ch':
#        ax.plot(xaxis,np.exp(het_all_notlog_list[0][idx]),'o',label="$a_l-a_r$ ="+str(val), markevery=200)
#    else: 
#        ax.plot(xaxis,het_all_notlog_list[0][idx],'o',label="$a_0$ ="+str(val), markevery=200)
#plt.legend(bbox_to_anchor=(0.0, 0.0), loc=3, borderaxespad=0.1)
##plt.xlabel('Diffusion constant, $a_0$',fontsize = fs,fontweight='bold')
#
#
#plt.xlabel('time (a.u.)',fontsize = fs*8/6,fontweight='bold')
#
#plt.ylabel('Heterozygosity, $H$',fontsize = fs*8/6,fontweight='bold')
#plt.tick_params(axis='both', which='major', labelsize=fs*8/6, pad =0.04) 
#plt.tight_layout()
#plt.savefig(path+'H_notlog_all.png')

















fig = plt.figure(figsize=(8,8))  
ax = fig.add_subplot(111)

for idx,val in enumerate(chirality[0]):
    
    
    #ax.plot(het_all_list[0][idx][:len(het_all_list[0][idx])/5],'o',label= foldertype +"="+str(val), markevery=50)
     ax.plot(het_all_list[0][idx][20:10000],'o',label= foldertype +"="+str(val), markevery=50)  
 
plt.legend(bbox_to_anchor=(0.0, 0.0), loc=3, borderaxespad=0.1)
#plt.xlabel('Diffusion constant, $a_0$',fontsize = fs,fontweight='bold')

if foldertype =='ch':
    plt.xlabel('chirality, $a_l-a_r$',fontsize = fs*8/6,fontweight='bold')
elif foldertype=='D':
    plt.xlabel('Diffusion constant, $a_0$',fontsize =fs*8/6,fontweight='bold')
plt.ylabel('log Heterozygosity, $\\log H $',fontsize = fs*8/6,fontweight='bold')
plt.tick_params(axis='both', which='major', labelsize=fs*8/6, pad =0.04) 
plt.tight_layout()
plt.savefig(path+'logH_all_initialstage.png')

#
#fig = plt.figure(figsize=(10,10))  
#ax = fig.add_subplot(111)
#
#for idx,val in enumerate(chirality[0]):
#    
#    
#    #ax.plot(het_all_list[0][idx][:len(het_all_list[0][idx])/5],'o',label= foldertype +"="+str(val), markevery=50)
#     ax.plot(het_all_list[0][idx][:200],'o',label= foldertype +"="+str(val), markevery=50)  
# 
#plt.legend(bbox_to_anchor=(0.0, 0.0), loc=3, borderaxespad=0.1)
##plt.xlabel('Diffusion constant, $a_0$',fontsize = fs,fontweight='bold')
#
#if foldertype =='ch':
#    plt.xlabel('chirality, $a_l-a_r$',fontsize = fs,fontweight='bold')
#elif foldertype=='D':
#    plt.xlabel('Diffusion constant, $a_0$',fontsize = fs,fontweight='bold')
#plt.ylabel('log Heterozygosity, $\\log H $',fontsize = fs,fontweight='bold')
#plt.tick_params(axis='both', which='major', labelsize=fs, pad =0.04) 
#plt.tight_layout()
#plt.savefig(path+'logH_all_EXTREMELY_initialstage.png')


#
plt.close("all")