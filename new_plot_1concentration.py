# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:32:44 2015

@author: ashish
"""

import numpy as np
from matplotlib.colors import colorConverter
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import os
import math
import pandas as pd
from copy import deepcopy
SS=6
r=40
#path = '/Volumes/Ashish_backup/research data/new_figures/PAPER FIGURES/oppch_coexistence_fig/gdnoise/try1_paper/N200/ch0.09/h0.05/1/'
suffix=""
path = '/Users/ashish/Downloads/N100/ch0.75/6/'
destpath =path

gz_flag=0

#Figure parameters
fs=24 #fontzize
sfs=14#small fontsize
lw=2#linewidth
ms=12#markersize


fold=path






if gz_flag==0:
    if not os.path.isfile(fold+'ca'+suffix+'.txt'):
        print "file not found"

    with open(fold+'ca'+suffix+'.txt') as file:
        ca = [[float(digit) for digit in line.split()] for line in file]

    with open(fold+'c'+suffix+'.txt') as file:       
        c = [[float(digit) for digit in line.split()] for line in file]
else:
    if not os.path.isfile(fold+'ca'+suffix+'.gz'):
        print " file not found"

    ca= pd.read_csv(fold+'ca'+suffix+'.gz',sep=' ').values
    c= pd.read_csv(fold+'c'+suffix+'.gz',sep=' ').values
    #ca=np.loadtxt(fold+'ca.gz')
    #c=np.loadtxt(fold+'c.gz')

ca=np.array(ca)
c=np.array(c)  
cb=c-ca 
        
l= len(ca)
b=len(ca.T)


#ca=np.roll(ca,l/6,axis=0)
#c=np.roll(c, l/6,axis=0)






#ca=ca[l/2-400:l/2+400,:400]
#c=c[l/2-400:l/2+400,:400]
##ca=np.roll(ca,l/3,axis=0)
##c=np.roll(c, l/3,axis=0)
#l= len(ca)
#b=len(ca.T)

cb=c-ca


#ctemp=deepcopy(ca)
#ca=cb
#cb=ctemp
#wallx=np.zeros(lexp)
#wally=np.zeros(lexp)        
#for idx, val in enumerate(wall):
#     wallx[idx]= idx
#     wally[idx]= val
#

fold = destpath   
if not os.path.exists(fold): os.makedirs(fold)

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
ax.set_title('species a diffusion')
sca=plt.imshow(ca)
plt.xlabel('space')
ax.set_aspect('equal')
cbar = fig.colorbar(sca);
#plt.plot(wallx[lexp/SS:],wally[lexp/SS:],'black',linestyle='--')
plt.savefig(fold+'ca'+suffix+'.png')






#fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(b/100,l/100))

fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))

img = np.zeros((l,b,3))
#img2 = np.ones((l,b,3))
#img[:,:,2]=cb[:,:]
#img[:,:,0]=ca[:,:]
#img[:,:,1]=ca[:,:]

img[:,:,0]=ca[:,:]+cb[:,:]
img[:,:,1]=ca[:,:]

            
ax.set_xticks([])
ax.set_yticks([])
img3 = ax.imshow(img,aspect='auto')

plt.savefig(fold+'fig3'+suffix+'.png') 




#fig4, ax = plt.subplots(nrows=1, ncols=1, figsize=(b/100,l/100))
fig4, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
img = np.ones((l,b,3))
#img[:,:,0]=1-cb[:,:]
#img[:,:,1]=1-ca[:,:]
#img[:,:,2]=1-cb[:,:]-ca[:,:]

img[:,:,0]=ca[:,:]
img[:,:,1]=cb[:,:]
img[:,:,2]=0

#
#circle=plt.Circle((l/2,b/2),r,color='black', fill=False)
#fig4.gca().add_artist(circle)

#ax.annotate('Initial droplet', xy=(l/2, b/2-r), xytext=(l/2, b/2-3*r),
#            arrowprops=dict(facecolor='black', shrink=0.05,width =1, headwidth=4),
#            )
#
ax.set_xticks([])
ax.set_yticks([])
img3 = ax.imshow(img)
#ax.plot(wallx[lexp/SS:],wally[lexp/SS:],'black',linestyle='--')
#ax.annotate('domain boundary', xy=(wallx[lexp/2], wally[lexp/2]), xytext=(wallx[lexp/2], wally[lexp/2]+2*r),ha ='center',
#            arrowprops=dict(facecolor='black', shrink=0.05,width =1, headwidth=5),
#            )     
plt.tight_layout()
plt.savefig(fold+'fig4'+suffix+'.png')  


sumc=np.sum(c,axis=0)

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
ax.set_title('sumc')
plt.plot(sumc,'o')

plt.savefig(fold+'sumc'+suffix+'.png')

#
#
# 
#beg=np.where(sumc<l-1)[0][0]
#end=np.where(sumc>0)[0][-1]+2
#
#
#if beg<0:     #just conditions for when we have diagonal expansion!
#    beg=np.where(sumc<10)[0][0]-1
#    
#    
#if (beg<end)and beg>0:
#    fig = plt.figure(figsize=(6,6))
#    ax = fig.add_subplot(111)
#    ax.set_title('sumc')
#    plt.plot(np.arange(beg,end),sumc[beg:end],'o')
#    
#    plt.savefig(fold+'sumc_zoom.png')
#    fig = plt.figure(figsize=(6,6))
#    ax = fig.add_subplot(111)
#    ax.set_title('sumc')
#    plt.semilogy(np.arange(beg,end),sumc[beg:end],'o')
#    
#    plt.savefig(fold+'sumc_semilog_zoom.png')    
#    print end-beg
#sumca=np.sum(ca,axis=0)
#fig = plt.figure(figsize=(6,6))
#ax = fig.add_subplot(111)
#ax.set_title('sumca')
#plt.plot(sumca,'o')
#
#plt.savefig(fold+'sumca.png') 





plt.show()


plt.close("all")
