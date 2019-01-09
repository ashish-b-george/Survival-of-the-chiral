import numpy as np
import glob
import os
import math
#SS=10
r=80
#path = '/projectnb/qedksk/ashishge/data/new_sim/boundary_motion/gdnoise/ch/try2_rescale_epsilon/'

path='/projectnb/qedksk/ashishge/data/new_sim/bulge_slope_from_track_height/gdnoise/try2/'
destpath=path[:-1]+"_wallExtract/"
chirality=[]
slopes=[]

## activate np.roll to extract walls that wrap  around

for ch in glob.glob(path+'*/*/*/'):	# for finding chirality value
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
    
    ## doing a roll to extract both walls that move so that zero crossings sees two walls 
    ca=np.roll(ca,10,axis=0)
    c=np.roll(c,10,axis=0) 
    cb=c-ca
    
    ## doing a roll to keep walls from wrapping around
    #c= np.roll(c,l/2,axis=0)
    #ca= np.roll(ca,l/2,axis=0)
    
    for i in range(0,b):
        if np.sum(c[:,i])<l*0.8:   #condition for the colony to have expanded
            lexp=int(i-10)
            break
    
    print lexp
    
    wall=np.zeros(lexp)
    

    dc=2*ca-c
    
    zero_crossings = np.where(np.diff(np.sign(dc[:,0])))[0] # all the zero crossing gives you number of initial walls 
    print  zero_crossings
    # that you want to track and  their starting positions! 
    
    zero_crossings_mid = np.where(np.diff(np.sign(dc[:,lexp/2])))[0]
    print "mid", zero_crossings_mid
    
    zero_crossings_end = np.where(np.diff(np.sign(dc[:,lexp-20])))[0]
    print "end", zero_crossings_end
    
    zero_crossings_34 = np.where(np.diff(np.sign(dc[:,int(round(lexp*0.75))])))[0]
    # you have to find walls in two ways, where there are equal species on both sides and 
    #where heterozygosity is equal on both sides
    
    for k in range(len(zero_crossings)): # this is the number of walls we want to track??
        
        lf =lexp# final length for domain wall after this value, the domain walls coalesce
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
                print k
                #print fold
                break
            topb=location-6 # we just search 20 points on either side of the previous wall point.
            botb=location+6
            

            
            wall[i]=(location+10)
            
            #topb=int(wall[i]-gap)
            #botb=int(wall[i]+gap) #this is the new area the search will run over
        destfold=fold.replace(path,destpath)
        np.savetxt(destfold+"/wall"+str(k)+".txt", wall[li:lf])
        
        #prints 1 if wall extraction failed at some value
        with open(destfold+"flag_wall"+str(k)+".txt", "w") as text_file:
            text_file.write(str(wall_flag))
        
        #prints li and lf of the wall extraction
        lilf=np.array((li,lf,lexp))
        np.savetxt(destfold+"/lilf"+str(k)+".txt", lilf)
        
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
                                     
    #het=ca*cb                  
    
                             
                                                      
                                                                   
                                                                                                        
                                                                                                                                 
                                                                                                                                                          
                                                                                                                                                                                                            
    #for k in range(len(zero_crossings)): # this is the number of walls we want to track??
    #    
    #    if k==0:
    #         topb=(zero_crossings[k])/2
    #         botb=(zero_crossings[k]+zero_crossings[k+1])/2  #this is to get the really chiral boundarries that move past 3/4 and so on
    #    elif k==len(zero_crossings)-1:
    #        topb=(zero_crossings[k]+zero_crossings[k-1])/2
    #        botb=(zero_crossings[k]+l)/2
    #    else:
    #        topb=(zero_crossings[k]+zero_crossings[k-1])/2
    #        botb=(zero_crossings[k]+zero_crossings[k+1])/2 
    #    #topb=(zero_crossings[k]+min(zero_crossings[k-1],0))/2
    #    #botb=(zero_crossings[k]+max(zero_crossings.take(k+1, mode='wrap'),l))/2  #this is to get the really chiral boundarries that move past 3/4 and so on
    #    print "botb"
    #    print "topb"
    #    print botb
    #    print topb
    #    for i in range(li, lexp):
    #        hetslice=het[:,i]
    #        #icaslice=ca[:,i]
    #                        
    #        errormin=float(2.0)
    #        ls=botb
    #        
    #        location =10000.0
    #        while (ls>topb) :
    ##                    if ls>0:
    #            address1=np.arange(ls,ls+10)
    #            address2=np.arange(ls+10,ls+20)
    #            error= np.sum(hetslice[address1])-np.sum(hetslice[address2])
    #            #print ls
    #            if (abs(error)<abs(errormin)):
    #                location =ls
    #                errormin=error
    #                print location
    #                #if (abs(error)<abs(errorminlim)):
    #                #    #to find the first wall
    #                #    break
    #                #print "errormin"
    #                #print errormin
    #                #print error
    #            ls=ls-1
    #            
    #        #print errormin
    #        #print location
    #        wall[i]=(location+10)
    #        
    #        topb=int(wall[i]-10)
    #        botb=int(wall[i]+10) #this is the new area the search will run over
    #
    #    np.savetxt(fold+"/wall_het"+str(k)+".txt", wall[li:lf])       