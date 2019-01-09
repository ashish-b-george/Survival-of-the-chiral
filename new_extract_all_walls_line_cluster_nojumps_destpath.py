import numpy as np
import glob
import os
import math
import matplotlib
from copy import deepcopy
matplotlib.use('GTK')
import matplotlib.pyplot as plt

#SS=10
r=20
path = '/projectnb/qedksk/ashishge/data/new_sim/boundary_motion/moran/try2/'
destpath='/projectnb/qedksk/ashishge/data/new_sim/boundary_motion/moran_walls_only/try2_test/'
#path ='/projectnb/qedksk/ashishge/data/new_sim/alpha_beta/moran/try2_N200_invasion/alpha/'

chirality=[]
slopes=[]


for ch in glob.glob(path+'*/*/*/1/'):	# for finding chirality value
#for ch in glob.glob(path+''):
    print ch
    fold=ch
    with open(fold+'ca.txt') as file:
        ca = [[float(digit) for digit in line.split()] for line in file]
    
    ca=np.array(ca)
    #print ca
    with open(fold+'c.txt') as file:       
        c = [[float(digit) for digit in line.split()] for line in file]
    
    
    c=np.array(c) 
              
    print  'gkfk'
    l=len(ca)
    b=len(ca.T)
    ## doing a roll to extract both walls that move (they wrap around otherwise)
    ca=np.roll(ca,l/3,axis=0)
    c=np.roll(c,l/3,axis=0) 
    cb=c-ca
    
    ## doing a roll to keep walls from wrapping around
    #c= np.roll(c,l/2,axis=0)
    #ca= np.roll(ca,l/2,axis=0)
    
    for i in range(0,b):
        if np.sum(c[:,i])<l*0.8:   #condition for the colony to have expanded
            lexp=i-10
            break
    
    print lexp
    
    wall=np.zeros(lexp)
    

    dc=2*ca-c
    
    zero_crossings = np.where(np.diff(np.sign(dc[:,0])))[0] # all the zero crossing gives you number of initial walls 
    print  zero_crossings
    # that you want to track and  their starting positions! 
    wall_list=[]
    zero_crossings_mid = np.where(np.diff(np.sign(dc[:,lexp/2])))[0]
    print "mid", zero_crossings_mid
    
    zero_crossings_end = np.where(np.diff(np.sign(dc[:,lexp-20])))[0]
    print "end", zero_crossings_end
    
    zero_crossings_34 = np.where(np.diff(np.sign(dc[:,round(lexp*0.75)])))[0]
    # you have to find walls in two ways, where there are equal species on both sides and 
    #where heterozygosity is equal on both sides
    
    for k in range(len(zero_crossings)): # this is the number of walls we want to track??

        lf =lexp -20 # final length for domain wall after this value, the domain walls coalesce
        li =r
        wall_flag=0
        
        if k==0:
             topb=(zero_crossings[k])/2
             ## if there is only one boundary detected then k+1  is invalid:/
             ## so in those cases issues can arise
             botb=(zero_crossings[k]+zero_crossings[k+1])/2  #this is to get the really chiral boundarries that move past 3/4 and so on
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
        
        for i in range(li, lf):
            dcslice=dc[:,i]
            #caslice=ca[:,i]
                            
            errormin=float(4.0)
            ls=botb
            
            location =-1000
            break_condition=1
            while (ls>topb) :
    #                    if ls>0:
                address=np.arange(ls,ls+20)
                #if(ls+21>l): # if the wall hasnt reached the edges of the domain
                #    if (ls<l):
                #         address=np.arange(ls,ls+21)
                        
                #error= np.sum(dcslice[address])
                error=np.sum(dcslice.take(address, mode='wrap'))

                if (abs(error)<abs(errormin)):
                    location =ls
                    errormin=error
                    #if (abs(error)<abs(errorminlim)):
                    #    #to find the first wall
                    #    break
    #                        print "errormin"
    #                        print errormin
    #                        print error
                ls=ls-1
                
                if (ls==topb+1 and location==-1000): # a solution to minimise was not found!!
                    #print "using this condition"
                    if break_condition==1: # we have not tried seraching in an expanded domain!
                        topb=int(np.mean(wall[li:i]))-int(2*gap)
                        botb=int(np.mean(wall[li:i]))+int(2*gap)
                        
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
                print k
                #print fold
                break
            topb=location-6 # we just search 20 points on either side of the previous wall point.
            botb=location+6
            

            
            wall[i]=(location+10)
            
            #topb=int(wall[i]-gap)
            #botb=int(wall[i]+gap) #this is the new area the search will run over
        
        destfold=fold.replace(path,destpath)    
        
        if not os.path.exists(destfold): os.makedirs(destfold)
        np.savetxt(destfold+"/wall"+str(k)+".txt", wall[li:lf])
        temp_wall=deepcopy(wall[li:lf])
        wall_list.append(temp_wall)
        

        #prints 1 if wall extraction failed at some value
        with open(destfold+"flag_wall"+str(k)+".txt", "w") as text_file:
            text_file.write(str(wall_flag))
        
        #prints li and lf of the wall extraction
        lilf=np.array((li,lf,lexp))
        np.savetxt(destfold+"/lilf"+str(k)+".txt", lilf)
        #fig = plt.figure(figsize=(6,6))
        #ax1 = fig.add_subplot(1, 1, 1)
        #plt.plot(wall)
        #ax1.set_xlabel('time')
        ##plt.tight_layout()
        #plt.savefig(destfold+'wall'+str(k)+'.png') 
    for i in range(len(zero_crossings)):
        print wall_list[i]
        print i,wall_list[i][0]
        
        
    wallx_coordinates=[]
    wally_coordinates=[]
    for i in range(len(zero_crossings)):
        wallx=np.zeros(lexp)
        wally=np.zeros(lexp)       
        for idx, val in enumerate(wall_list[i]):
            
            wallx[idx]= li+ idx
            wally[idx]= val
        wallx_coordinates.append(np.array(wallx))
        wally_coordinates.append(np.array(wally))
    print np.shape(  wallx_coordinates)              
    fig3 = plt.figure(figsize=(6,6))  
    ax = fig3.add_subplot(1, 1, 1)
    #fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
    img = np.zeros((l,b,3))
    img2 = np.ones((l,b,3))
    img[:,:,0]=ca[:,:]
    img[:,:,1]=cb[:,:]
    ax.set_xticks([])
    ax.set_yticks([])
    img3 = ax.imshow(img)
    for i in range(len(zero_crossings)):
        print i
        ax.plot(wallx_coordinates[i][li:],wally_coordinates[i][li:],'black',linestyle='--')
    plt.savefig(destfold+'fig3.png') 
    fig4 = plt.figure(figsize=(b/200,l/200))  
    ax = fig4.add_subplot(1, 1, 1)
    #fig4, ax = plt.subplots(nrows=1, ncols=1, figsize=(b/100,l/100))
    img = np.ones((l,b,3))
    img[:,:,0]=1-cb[:,:]
    img[:,:,1]=1-ca[:,:]
    img[:,:,2]=1-cb[:,:]-ca[:,:]
    ax.set_xticks([])
    ax.set_yticks([])
    img3 = ax.imshow(img)
    for i in range(len(zero_crossings)):  
        ax.plot(wallx_coordinates[i][li:],wally_coordinates[i][li:],'black',linestyle='--')
    plt.savefig(destfold+'fig4.png') 
        
        
        
        
    newdomaincounter=0
    if len(zero_crossings)!=len(zero_crossings_end):
        newdomaincounter=newdomaincounter+1
    
    if len(zero_crossings)!=len(zero_crossings_mid):
        newdomaincounter=newdomaincounter+1
    if len(zero_crossings)!=len(zero_crossings_34):
        newdomaincounter=newdomaincounter+1
        ## 0 if no new domains, 1 , 2 etc if either at end or middle or 3/4ths there were new domains, 
        ##3 if there were new domains throughout that lived from the beginning to the end
    with open(destfold+"domain_birth_flag.txt", "w") as text_file:
        text_file.write(str(newdomaincounter))
    
    print "new domain counter", newdomaincounter 
    

plt.close("all")                                     
    #het=ca*cb                  
    
                             
                                                      
                                                                   
                                                                                                        
                                                                                                                                 
                                                                                                                                                          
                                                                                                                                                                                                            
