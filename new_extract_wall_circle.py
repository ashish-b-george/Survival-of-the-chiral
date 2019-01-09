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
#from scipy.optimize import curve_fit
#doesnt average over runs!! Just the run1 value is used!!!

chirality=[]
slopes=[]
path = '/Users/ashish/Downloads/meanfield_fit/'
#path = '/Users/ashish/Downloads/figures/avgvel/chn'
Nvalue=''
destpath= '/Users/ashish/Downloads/logspiral/meanfield_fit/'
#dest= '/usr3/graduate/ashishge/saved data'
for ch in glob.glob(path+'*/'):    # for finding chirality value
    if os.path.exists(ch+Nvalue):
        #FOR FINDING A COMMON LEXP FOR ALL THE RUNS!!
        fold1=ch+Nvalue+'1/'
        with open(fold1+'ca.txt') as file:
            ca = [[float(digit) for digit in line.split()] for line in file]
        
        ca=np.array(ca)
        print ca
        with open(fold1+'c.txt') as file:       
            c = [[float(digit) for digit in line.split()] for line in file]
        
        
        c=np.array(c)  
        cb=c-ca              
        print  'gkfk'
        l=len(ca)
        b=len(ca.T)
        for i in range(0,b/2):
            if c[l/2,b/2+i]<.85:   #condition for the colony to have expanded
                lexp=i-10
                break
        
        print lexp #length of expansion
        # LEXP IS i-`10!! just to be safe
       
        cumwall=np.zeros(lexp)
        wall=np.zeros(lexp)
        cumca=np.zeros((l,b))
        cumc=np.zeros((l,b))
        
        runcounter=0
        for fold in glob.glob(ch+ Nvalue+'*/'):
            runcounter=runcounter+1
            with open(fold+'ca.txt') as file:
                ca = [[float(digit) for digit in line.split()] for line in file]
            
            ca=np.array(ca)
            print ca
            with open(fold+'c.txt') as file:       
                c = [[float(digit) for digit in line.split()] for line in file]
            
            
            c=np.array(c)  
            cb=c-ca              
            l= len(ca)
            b=len(ca.T)
            cumca=cumca+ca
            cumc=cumc+c
            
#            fig = plt.figure(figsize=(6,6))
#            ax = fig.add_subplot(111)
#            ax.set_title('species a diffusion')
#            sca=plt.imshow(ca)
#            plt.xlabel('space')
#            ax.set_aspect('equal')
#            cbar = fig.colorbar(sca);
#            plt.savefig(fold+'ca.png')
#            
#            dc=2*ca-c
#            fig = plt.figure(figsize=(6,6))
#            ax = fig.add_subplot(111)
#            ax.set_title('ca-cb')
#            sca=plt.imshow(dc)
#            plt.xlabel('space')
#            ax.set_aspect('equal')
#            cbar = fig.colorbar(sca);
#            plt.savefig(fold+'dc.png')



            
            SS=6 #fraction of lexp to exclude to reach steady state
            
            dc=2*ca-c
            for r in range(lexp/SS, lexp):
                mf=6.0
                # 6r points because 6r>4r number of points on perimeter of square
                thetaray=np.zeros(mf*r)
                for theta in range(0,int(mf*r)):
                    x= r* math.cos(2.0* math.pi * float(theta)/(mf*r))
                    y= r* math.sin(2.0* math.pi * float(theta)/(mf*r))
                    x1=math.floor(x)
                    y1=math.floor(y)
                    x2=math.ceil(x)
                    y2=math.ceil(y)
#                    x2-x1 and y2-y1 are both POSITIVE ONE so first fraction is not needed
                    if x2==x1:
                        x2=x1+1
                    if y2==y1:
                        y2=y1+1
                    thetaray[theta]= (x2-x)*(y2-y)*dc[l/2+x1,b/2+y1]+(x-x1)*(y2-y)*dc[l/2+x2,b/2+y1]+(x2-x)*(y-y1)*dc[l/2+x1,b/2+y2]+(x-x1)*(y-y1)*dc[l/2+x2,b/2+y2]
              
                errormin=float(8.0)
                ls=int(0)
                location =10000.0
                while (ls>-mf*r/2) :
#                    if ls>0:
                    address=np.arange(ls,ls+21)
                    error= np.sum(thetaray[address])
#                    if ls<0:
#                        error=np.sum(thetaray[-20+ls:21])
                    if (abs(error)<abs(errormin)):
                        location =ls
                        errormin=error
#                        print "errormin"
#                        print errormin
#                        print error
                    ls=ls-1
#                print errormin
#                print location
#                print thetaray[location:location+21]

                
                wall[r]=(location+10)*2*math.pi/(mf*r)
 
#SOMETHING IS GOIGN WRONG IN THE CUMWALLSTEPS;/               
            cumwall=cumwall+wall
        cumwall=cumwall/runcounter
        cumca=cumca/runcounter
        cumc=cumc/runcounter
        fold=ch+ Nvalue        
#        fig = plt.figure(figsize=(6,10))  
#        ax = fig.add_subplot(111)
#        ax.set_title('total wall')
#        plt.plot(cumwall)
#
#            #linis the length to which we assume the expansion has not reached steady state
#        lin=lexp/SS+1
#        lst=lexp-lin
#        xdata=np.arange(lin,lst+lin)
#        ydata=cumwall[lin:]
        
        
        
#        CURVEFIT ISNT WORKING AT ALL FOR SOME REASON.
#        def func(z, p1,p2):
#            return p1*np.log(z)+p2
#        popt, pcov = curve_fit(func, xdata, ydata)
#        
#        k=np.zeros(lst)
#        for i in range(0,lst):
#            k[i]=func(i,popt[0],popt[1])
#
#        fig = plt.figure(figsize=(6,10))
#        plt.plot(xdata,ydata,xdata,k)
#        plt.savefig(fold+'wall.png') 
#        print popt
        
        
#        SO FITTING POLY FIT WITH LOG FOR LOG FIT
#        logx=np.log(xdata)
#        
#        z=np.polyfit(logx, ydata, 1)
#        p=np.poly1d(z)
#
#        k=np.zeros(lst)
#            
#        for i in range(0, lst):
#            k[i]=p(logx[i])
#            
#        print z
#        fig = plt.figure(figsize=(6,10))
#        plt.plot(logx,ydata,logx,k)
#        plt.title("log plot for logfit")
#        plt.savefig(fold+'logwall.png') 
## FOR 1/R FIT I AM FITTING IN A LOG LOG HOPEFULLY GETTIN ALMOST ONE FOR EXPONENT
#        logy=np.log(np.fabs(ydata))
#        z2=np.polyfit(logx, logy, 1)
#        p2=np.poly1d(z2)
#
#        k=np.zeros(lst)
#            
#        for i in range(0, lst):
#            k[i]=p2(logx[i])
#            
#        print z2
#        fig = plt.figure(figsize=(6,10))
#        plt.plot(logx,logy,logx,k)
#        plt.title("log log plot for 1/r fit")
#        plt.savefig(fold+'loglogwall.png')
        chiralvalue=ch.replace(path+'ch','')
        print chiralvalue
        chiralvalue2=ch.replace(path,'')
        newpath = destpath+Nvalue+chiralvalue2
        if not os.path.exists(newpath): os.makedirs(newpath)
#        np.savetxt(newpath+"/fit.txt",z)
        np.savetxt(newpath+"/wall.txt", cumwall)
        np.savetxt(newpath+"/ca.txt", ca)
        np.savetxt(newpath+"/c.txt", c)
        
        
#   
#        chirality.append(float(chiralvalue))
#        slopes.append(z[[0]])
#    plt.close("all")
#newpath2= destpath+Nvalue
#if not os.path.exists(newpath2): os.makedirs(newpath2)
#np.savetxt(newpath2+'chirality.txt',chirality)
#np.savetxt(newpath2+'slopes.txt',slopes)

#fig = plt.figure(figsize=(6,10))  
#ax = fig.add_subplot(111)
#ax.set_title('chirality variation')
#plt.plot(chirality,slopes,'*') 
#plt.savefig(newpath2+'ch_vs_slope.png')
