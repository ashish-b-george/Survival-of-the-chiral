#!/usr/bin/python
import numpy as np
#import matplotlib.pyplot as plt
import glob
import os
import math
#import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
#doesnt average over runs!! Just the run1 value is used!!!

chirality=[]
slopes=[]
#path = '/projectnb/qedksk/ashishge/data/new_sim/neutral_Hdecay/N/gdnoise/try2_rescale/'
#Nvalue='/ch0.75/'
#destpath= '/projectnb/qedksk/ashishge/data/new_sim/neutral_Hdecay/N/gdnoise/try2_rescale_avg/'




path = '/projectnb/qedksk/ashishge/data/new_sim/effect_of_N/samech_gdnoise/try5/'
destpath= '/projectnb/qedksk/ashishge/data/new_sim/effect_of_N/samech_gdnoise/try5_avg/'

single_run_no_avg_flag=0

print "runnning"
lvl2foldertype='ch'
SS=6
#*/ can be removed in lvl1 if rquried!


for lvl1 in glob.glob(path+''):  
    print lvl1	 
    for ch in glob.glob(lvl1+'*/'):    # for finding chirality value
        for Nfold in glob.glob(ch +'*/'):
            #if os.path.exists(ch+Nvalue):
            #fold1=ch+Nvalue+'1/'
            fold1=Nfold+'1/'
            #print fold1
            with open(fold1+'ca.txt') as file:
                ca = [[float(digit) for digit in line.split()] for line in file]
            
            ca=np.array(ca)
            #print ca
            with open(fold1+'c.txt') as file:       
                c = [[float(digit) for digit in line.split()] for line in file]
            with open(fold1+'time.txt') as file:
                    time = [[float(digit) for digit in line.split()] for line in file]
            print "time is ", time
            time=int(time[0][0])
            c=np.array(c)  
            cb=c-ca              

            l=len(ca)
            b=len(ca.T)

            cumca=np.zeros((l,b))
            cumc=np.zeros((l,b))
            cumcb=np.zeros((l,b))
            oldca=np.zeros((l,b))
            sumc=np.zeros(b)
            cumHET=np.zeros(int(time))            
            mintime=int(time)-1
            min_lexp=b-1
            cumHET_spatial=np.zeros(b)            
            cblist=[]
            calist=[]            
            cumtimelist=[]
            lenpercent=[]
            
            hetlist=[]
            het_timelist=[]

            
            runcounter=0
            
            
            for fold in glob.glob(Nfold+'*/'):
                if single_run_no_avg_flag!=0:
                    fold=Nfold+str(single_run_no_avg_flag)+'/'
                
                runcounter=runcounter+1
                with open(fold+'ca.txt') as file:
                    ca = [[float(digit) for digit in line.split()] for line in file]
                with open(fold+'time.txt') as file:
                    time = [[float(digit) for digit in line.split()] for line in file]
                time=int(time[0][0])
                ca=np.array(ca)
                #print ca
                with open(fold+'c.txt') as file:       
                    c = [[float(digit) for digit in line.split()] for line in file]
                with open(fold+'heterozygosity_time.txt') as file:       
                    HET = [[float(digit) for digit in line.split()] for line in file]
                HET=np.array(HET)
                HET=np.ravel(HET)
                print fold
                c=np.array(c)  
                cb=c-ca              
                l= len(ca)
                b=len(ca.T)
                              
                #for i in range(0,b):
                #    if c[l/2,i]<.8:   #condition for the colony to have expanded
                #        lexp=i-10
                #        break
            
                c_mean=np.mean(c,axis=0)
                lexp=np.where(c_mean<0.9)[0][0]-20
                
                
 
                cumtimelist.append(time)
                cblist.append(np.sum(cb))
                calist.append(np.sum(ca))

                if mintime>time-1:
                    mintime=time-1
                    cumHET[:mintime]=cumHET[:mintime]+HET[:mintime]
                else:
                    cumHET[:mintime]=cumHET[:mintime]+HET[:mintime]
                
                #het_spatial=np.sum(2.*ca[:,lexp-20:lexp]*cb[:,lexp-20:lexp],axis=0)/(20.0*l)
                
                
                het_spatial=np.sum(2.* ( ca[:,:lexp]/c[:,:lexp] )  * (cb[:,:lexp]/c[:,:lexp]),axis=0)/(l*1.0)
                ## we can calculate the heterozygosity in space as well
                if min_lexp>lexp:
                    min_lexp=lexp
                    cumHET_spatial[:min_lexp]=cumHET_spatial[:min_lexp]+het_spatial[:min_lexp]
                else:
                    cumHET_spatial[:min_lexp]=cumHET_spatial[:min_lexp]+het_spatial[:min_lexp]
                


                if single_run_no_avg_flag!=0:
                    print "just single  run with fold= single_run_no_avg_flag chosen, no averaging!"
                    break


                #het_list.append(np.sum(2*ca[:,lexp-20:lexp]*cb[:,lexp-20:lexp])/(20*l))
                #sumc=sumc+np.sum(c,axis=0) # sum of c in slices along direction of expansion
                #het_exp=np.sum(2*ca*cb/l,axis=0)
                #sumhet_exp=sumhet_exp+het_exp
       
            #print hetlist
            #hetlist=np.array(hetlist)
            #slope1list=np.array(slope1list)
            #slope2list=np.array(slope2list)
            #velalist=np.array(velalist)
            #velblist=np.array(velblist)
            #sumc=np.array(sumc)/runcounter
            #sumhet_exp=sumhet_exp/runcounter
            calist=np.array(calist)
            cblist=np.array(cblist)
            
            #lenpercent=np.array(lenpercent)
            #chiralvalue=ch.replace(lvl1+lvl2foldertype,'')
            #print chiralvalue
            #chiralvalue2=ch.replace(path,'')
            

            cumHET=cumHET/runcounter
            cumHET_spatial=cumHET_spatial/runcounter
            newpath = Nfold.replace(path,destpath)
            if not os.path.exists(newpath): os.makedirs(newpath)

            print newpath, runcounter

            np.savetxt(newpath+"/capops.txt", calist)
            np.savetxt(newpath+"/cbpops.txt", cblist)
            #np.savetxt(newpath+"/sumc.txt", sumc)
            np.savetxt(newpath+"/time.txt", cumtimelist)
            np.savetxt(newpath+"/HET_time.txt", cumHET[:mintime])
            np.savetxt(newpath+"/HET_spatial_time.txt", cumHET_spatial[:min_lexp])
            
            #np.savetxt(newpath+"/slope1.txt", slope1list)
            #np.savetxt(newpath+"/slope2.txt", slope2list)
            #np.savetxt(newpath+"/vela.txt", velalist)
            #np.savetxt(newpath+"/velb.txt", velblist) 
            #np.savetxt(newpath+"/heterozygosity.txt", hetlist)
            #np.savetxt(newpath+"/lenpercent.txt", lenpercent)
            #np.savetxt(newpath+"/het_exp_avg.txt", sumhet_exp)
