#!/usr/bin/python
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from copy import deepcopy
'''
This is a code to implement a prelimiary off lattice chiral simulation 
'''
class population_Cpp:
    def __init__(self):
        #self.xpos=np.zeros(max_population_size)
        #self.ypos=np.zeros(max_population_size) 
        self.max_population_size=0
        self.coarseningFactor=0
        #self.xpos=np.zeros(10)
        #self.pos =np.zeros((max_population_size,2)) ## [0] is x and [1] is y
    def read_data(self,path,pos_filename):
        positions=np.loadtxt(path+pos_filename)
        self.xpos=np.ravel(positions[0])
        self.ypos=np.ravel(positions[1])
        label_filename=pos_filename.replace("pos","label")
        self.species_label=np.loadtxt(path+label_filename)     
        time_filename=pos_filename.replace("pos","time")
        self.simulation_time=np.loadtxt(path+time_filename)   
        params_filename=pos_filename.replace("pos","params")
        self.params=pd.read_table(path+params_filename,header=None, sep=' ',lineterminator="\n") ## header =none means 0 and 1 are headers used
                        ## otherwise the first parameter would become the header
        self.params.set_index(0, inplace=True)  ### now we can address each row by  kk.loc["f_IC"] or kk.loc["f_IC"].values for the value

    def plot_histograms(self,path,file_name,coarsening_factor):
        
        bin_size=self.params.loc["spread_in_concentration_sensing"].values[0]
        binValues=np.arange(0,round(np.max(self.ypos)+10*bin_size),bin_size)
        #w_x_K=0.5*self.params.loc["per_deme_carrying_capacity"].values[0]*self.params.loc["domain_width"].values[0] ## width times K to normalise to1

        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        #plt.hist(self.ypos,bins=binValues ,facecolor="blue",edgecolor="none",weights=1./w_x_K*np.ones(len(self.ypos)) )
        plt.hist(self.ypos,bins=binValues ,facecolor="blue",edgecolor="none" )
        plt.hist(self.ypos[self.species_label==1],bins=binValues,facecolor="red",edgecolor="none")
        plt.savefig(path+file_name+".png")
 
        binValues=np.arange(0,round(np.max(self.ypos)+10*coarsening_factor*bin_size),coarsening_factor*bin_size)
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        plt.hist(self.ypos,bins=binValues ,facecolor="blue",edgecolor="none" )
        plt.hist(self.ypos[self.species_label==1],bins=binValues,facecolor="red",edgecolor="none")
        plt.savefig(path+file_name+"_"+str(coarsening_factor)+"x_coarser.png")
                
    def plot_scatter(self,path,filename):
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        plt.scatter(self.xpos, self.ypos, c=self.species_label,s=2, cmap='Accent',edgecolor="none")
        plt.savefig(path+filename+".png")



    def plot_image(self,path,filename,coarsening_factor=None): 
        
        #if coarsening_factor==None:
        #    coarsening_factor=1.
        #dx= self.params.loc["spread_in_concentration_sensing"].values[0]
        #dx=dx*coarsening_factor
        #xgrid=np.arange(0,round(np.max(self.xpos)),dx)
        #ygrid=np.arange(0,round(np.max(self.ypos)+10*dx),dx) 
        #
        #concentration=np.zeros( (len(xgrid)+1,len(ygrid)+1) )
        #concentration_a=np.zeros( (len(xgrid)+1,len(ygrid)+1) )
        ###normalise so that just rounding gets you the closest lattice point
        #xpos_normed=self.xpos/dx
        #ypos_normed=self.ypos/dx
        #for i in range(len(self.xpos)):
        #    xval=int(round(xpos_normed[i]))
        #    yval=int(round(ypos_normed[i]))
        #    concentration[xval,yval]+=1
        #    if self.species_label[i]==1:
        #        concentration_a[xval,yval]+=1
        #concentration_b= concentration- concentration_a   
        
        self.calculate_concentrations(coarsening_factor) 
        img = np.zeros( (len(self.xgrid),len(self.ygrid),3) )
        img[:,:,0]=self.concentration_b[:,:]
        img[:,:,1]= self.concentration_a[:,:]
        
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))  
        ax.imshow(img,aspect="equal") 
        ax.set_xticks([])
        ax.set_yticks([])
        plt.savefig(path+filename+".png")
        
        
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        plt.imshow(self.concentration, cmap='plasma' ,aspect="equal",interpolation='nearest')
        plt.colorbar()
        plt.savefig(path+filename+"_totalc.png")
        
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        plt.imshow(self.local_heterozygosity, cmap='plasma' ,aspect="equal",interpolation='nearest')
        plt.colorbar()
        plt.savefig(path+filename+"_het.png")
        
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        plt.imshow(self.concentration_a, cmap='plasma' ,aspect="equal",interpolation='nearest')
        plt.colorbar()
        plt.savefig(path+filename+"_a.png")
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        plt.imshow(self.concentration_b, cmap='plasma' ,aspect="equal",interpolation='nearest')
        plt.colorbar()
        plt.savefig(path+filename+"_b.png")
        
        
        
        
    def calculate_concentrations(self,coarsening_factor=None):          
        if coarsening_factor==None:
            self.coarseningFactor=1.
        else :
            self.coarseningFactor=coarsening_factor
        dx= self.params.loc["spread_in_concentration_sensing"].values[0]
        dx=dx*self.coarseningFactor
        K=self.params.loc["per_deme_carrying_capacity"].values[0]
        self.xgrid=np.arange(0,round(np.max(self.xpos)),dx)
        self.ygrid=np.arange(0,round(np.max(self.ypos)+10*dx),dx)       
        self.concentration=np.zeros( (len(self.xgrid),len(self.ygrid)) )
        self.concentration_a=np.zeros( (len(self.xgrid),len(self.ygrid)) )
        ##normalise so that just rounding gets you the closest lattice point
        xpos_normed=self.xpos/dx
        ypos_normed=self.ypos/dx
        xpos_normed_max=int(round (round(np.max(self.xpos) )/dx))
        for i in range(len(self.xpos)):
            xval=int(   int(round(xpos_normed[i]))%xpos_normed_max  )### this implements periodic boundary conditions!
            yval=int(round(ypos_normed[i]))
            self.concentration[xval,yval]+=1
            if self.species_label[i]==1:
                self.concentration_a[xval,yval]+=1
                
        self.concentration=self.concentration/K
        self.concentration_a=self.concentration_a/K
        self.concentration_b= self.concentration- self.concentration_a 
        self.local_heterozygosity=2*self.concentration_b*self.concentration_a
        
        
    def calculate_fraction_y(self,fraction_coarsening_factor=None):    
        if fraction_coarsening_factor==None:
            fraction_coarsening_factor=1.
        self.calculate_concentrations(fraction_coarsening_factor)        
        lexp=np.where( np.sum(self.concentration,axis=0) <0.1*len(self.xgrid)  )[0][0] -2
        self.fraction_y=  np.sum(self.concentration_a,axis=0)[:lexp]/np.sum(self.concentration,axis=0)[:lexp]
        return self.fraction_y, self.ygrid[:lexp]
        
        
    def calculate_height_and_roughness(self):
        if self.coarseningFactor!=1:    ## will be 1 if previously calcualted with coarsenign factor of 1, in whcih case no need to recalculate
            self.calculate_concentrations()
        self.Height_x=np.sum(self.concentration,axis=1) 
        self.roughness=np.var(self.Height_x)
        #self.Roughness=np.var(self.concentration,axis=1)
        
    def plot_height(self,path,height_filename):
        self.calculate_height_and_roughness()
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        plt.plot(self.xgrid, self.Height_x )
        plt.savefig(path+height_filename+".png")
  
    ###### fraction of species in the moving box where box length is the box length used in the simulations.. #####
    def calculate_Fraction_Globalhet_at_front(self):
        y_max=np.max(self.ypos)
        box_length=self.params.loc["box_length"].values[0]       
        idx_cells_at_front=np.where(self.ypos>y_max-box_length)        
        fraction_at_front=np.mean(self.species_label[idx_cells_at_front])
        globalhet_at_front=np.mean( np.ravel(self.species_label[idx_cells_at_front])*np.ravel(self.species_label[idx_cells_at_front]) )
        return fraction_at_front,globalhet_at_front
        

        
    def calculate_boundaries2(self):      
        
        domain_a=np.zeros( len(self.ygrid) ) ## locations of max ca and cb
        domain_b=np.zeros( len(self.ygrid) )
        
        boundary1=np.zeros( len(self.ygrid) )
        boundary2=np.zeros( len(self.ygrid) )
        
        boundary_flag=np.zeros(len(self.ygrid) ).astype(int)
        
        coarsening_factor=20 ## the size of the window we average over 
        for y_idx,y in enumerate(self.ygrid):
             if np.mean(self.concentration[:,y_idx]>0.0):
                wrapped_c=np.pad(self.concentration[:,y_idx],(0,coarsening_factor-1),mode="wrap") ## wrapped by adding at the end only!!!
                wrapped_ca=np.pad(self.concentration_a[:,y_idx],(0,coarsening_factor-1),mode="wrap")
                wrapped_cb=np.pad(self.concentration_b[:,y_idx],(0,coarsening_factor-1),mode="wrap")
                
                coarsened_c= np.convolve(wrapped_c, np.ones(coarsening_factor), mode="valid") ## length =self.conc and <wrapped_c
                coarsened_ca= np.convolve(wrapped_ca, np.ones(coarsening_factor), mode="valid")
                coarsened_cb= np.convolve(wrapped_cb, np.ones(coarsening_factor), mode="valid")
                domain_a[y_idx]=np.argmax(coarsened_ca)
                domain_b[y_idx]=np.argmax(coarsened_cb)
                
                het=self.local_heterozygosity[:,y_idx]
#                if domain_a[y_idx] <domain_b[y_idx]:
#                    
#                    boundary1[y_idx] =( np.argmax( np.take( het,np.arange(domain_a[y_idx], domain_b[y_idx]).astype(int),mode="wrap" ) ) +domain_a[y_idx] )%len(self.xgrid)
#                    boundary2[y_idx] =( np.argmax( np.take( het,np.arange(domain_b[y_idx] ,domain_a[y_idx] +len(self.xgrid) ).astype(int),mode="wrap" )  )  +domain_b[y_idx] )%len(self.xgrid) 
#
#                    boundary_flag[y_idx]=1
#                    
#                    #print  np.arange(domain_a[y_idx], domain_b[y_idx]).astype(int),het[np.arange(domain_a[y_idx], domain_b[y_idx]).astype(int)]
#                    #print "max",np.max(het), np.argmax(het), boundary1[y_idx] ,boundary2[y_idx] 
#                    #print het[boundary1[y_idx]],het[boundary2[y_idx]]
#                elif domain_a[y_idx] >domain_b[y_idx]:
#                    
#                    boundary1[y_idx] =(np.argmax( np.take( het,np.arange(  domain_a[y_idx], domain_b[y_idx]+len(self.xgrid) ).astype(int),mode="wrap" ) ) +domain_a[y_idx] )%len(self.xgrid)
#                    boundary2[y_idx] =(np.argmax( np.take( het,np.arange( domain_b[y_idx] ,domain_a[y_idx]  ).astype(int),mode="wrap" ) ) +domain_b[y_idx] )%len(self.xgrid)   
#                    
#                    boundary_flag[y_idx]=1

                    
                    
                
                boundary1[y_idx]=np.argmin( np.roll(coarsened_ca,10)-np.roll(coarsened_cb,-10) )
                boundary2[y_idx]=np.argmin( np.roll(coarsened_cb,10)-np.roll(coarsened_ca,-10) )
                
                boundary_flag[y_idx]=1
                
        return (boundary1,boundary2, boundary_flag,domain_a, domain_b)
    
        
            
    def calculate_boundaries(self):   ####### adapted from new_extract_all_walls_line_cluster_nojumps_destpath.py#####
        
        dc=self.concentration_a-self.concentration_b
        r=5
        coarsening_factor=10
        lexp=0
        wrapped_c=np.pad(self.concentration[:,r],(0,coarsening_factor-1),mode="wrap") ## wrapped by adding at the end only!!!
        wrapped_ca=np.pad(self.concentration_a[:,r],(0,coarsening_factor-1),mode="wrap")
        wrapped_cb=np.pad(self.concentration_b[:,r],(0,coarsening_factor-1),mode="wrap")
        coarsened_c= np.convolve(wrapped_c, np.ones(coarsening_factor), mode="valid") ## length =self.conc and <wrapped_c
        coarsened_ca= np.convolve(wrapped_ca, np.ones(coarsening_factor), mode="valid")
        coarsened_cb= np.convolve(wrapped_cb, np.ones(coarsening_factor), mode="valid")
        domain_a=np.argmax(coarsened_ca)
        domain_b=np.argmax(coarsened_cb)
        
        for y_idx in range( len(self.ygrid) ):
             if np.mean(self.concentration[:,y_idx]>0.0):
                 lexp=y_idx

        lx=len(self.xgrid)
        ly=len(self.ygrid)
        
        lf =lexp -5 # final length for domain wall after this value, the domain walls coalesce
        li =r
        
        
        #  find walls as where there are equal species on both sides and 
        
        
        for k in range(2): # this is the number of walls we want to track
            wall=np.zeros(ly)
            boundary_flag=np.zeros(ly)
            wall_flag=0  
            if k==0:
                topb=min(domain_a,domain_b)
                botb=max(domain_a,domain_b)
            elif k==1:
                topb=max(domain_a,domain_b)
                botb=min(domain_a,domain_b)+lx
            #else:
            #    topb=(zero_crossings[k]+zero_crossings[k-1])/2
            #    botb=(zero_crossings[k]+zero_crossings[k+1])/2 
            gap=(botb-topb)/2
            #print "gap",gap
            #print li,lf
            for i in range(li, lf):
                dcslice=dc[:,i]
                #caslice=ca[:,i]
                                
                errormin=float(4.0)
                ls=botb
                
                location =-1000
                break_condition=1
                while (ls>topb) :
                    address=np.arange(ls,ls+20)
                    error=np.sum(dcslice.take(address, mode='wrap'))
    
                    if (abs(error)<abs(errormin)):
                        location =ls
                        errormin=error
                        boundary_flag[i]=1
                        #if (abs(error)<abs(errorminlim)):
                        #    #to find the first wall
                        #    break
                        
        #                        print errormin
        #                        print error
                    ls=ls-1
                    
                    if (ls==topb+1 and location==-1000): # a solution to minimise was not found!!
                        print "using this condition"
                        if break_condition==1: # we have not tried seraching in an expanded domain!
                            if i>li:
                                topb=int(np.mean(wall[li:i]))-int(0.6*gap)
                                botb=int(np.mean(wall[li:i]))+int(0.6*gap)
                            else :
                                botb=botb+10
                                topb=topb-10
                            
                            ls=botb
                            #print ls
                            #print topb
                            break_condition=0
                        else:
                            print "still not found.. wtf,"
                            #lf=i-10
                            wall_flag+=1
                            break
                            #location =-1000
                
                
                if wall_flag>10:
                    print " breaking, wall not found atleast for 5 values of ls"
                    print k
                    break
                if boundary_flag[i]==1:    
                    topb=location-10 # we just search 20 points on either side of the previous wall point.
                    botb=location+10
                
                wall[i]=(location+10)
                
                #topb=int(wall[i]-gap)
                #botb=int(wall[i]+gap) #this is the new area the search will run over
            

            if k==0:
                boundary1=deepcopy(wall)
                boundary_flag1=deepcopy(boundary_flag)
            elif k==1:
                boundary2=deepcopy(wall)
                boundary_flag2=deepcopy(boundary_flag)

        print len(boundary_flag1)  ,np.sum(boundary_flag1)
        return boundary1,boundary2, boundary_flag1 ,boundary_flag2               
                                                  
                        
    def plot_domain_boundaries(self,path,filename): 
        boundary1,boundary2,boundary_flag, domain_a, domain_b =self.calculate_boundaries2()
        #print boundary1[boundary1!=boundary2]
        #print boundary1[boundary1==boundary2]
        #print boundary1!=boundary2
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        plt.plot(self.ygrid[boundary_flag==1], boundary1[boundary_flag==1] ,'o')
        plt.savefig(path+filename+"_1alt.png")
        
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        plt.plot(self.ygrid[boundary_flag==1], boundary2[boundary_flag==1] ,'o')
        plt.savefig(path+filename+"_2alt.png")
        
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        plt.plot(self.ygrid, domain_a,'o')
        plt.savefig(path+"domain_a.png")
        
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        plt.plot(self.ygrid, domain_b,'o')
        plt.savefig(path+"domain_b.png")        
        boundary1,boundary2,boundary_flag1,boundary_flag2 =self.calculate_boundaries()
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        plt.plot(self.ygrid[boundary_flag1==1], boundary1[boundary_flag1==1] ,'o')
        plt.savefig(path+filename+"_1.png")
        
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        plt.plot(self.ygrid[boundary_flag2==1], boundary2[boundary_flag2==1] ,'o')
        plt.savefig(path+filename+"_2.png")
    
    


    def plot_sensedC(self,path,sensedC_filename,destfold=None):
        sensed_C_all=np.loadtxt(path+sensedC_filename+".txt")
        if destfold==None:
            destfold=path
        self.sensed_C=np.ravel(sensed_C_all[0])
        self.sensed_Cx=np.ravel(sensed_C_all[1])
        self.sensed_Cy=np.ravel(sensed_C_all[2])
        #print "mean max min ",np.mean(self.sensed_C),np.max(self.sensed_C),np.min(self.sensed_C)
        #print "X mean max min ",np.mean(self.sensed_Cx),np.max(self.sensed_Cx),np.min(self.sensed_Cx)
        #print "Y mean max min ",np.mean(self.sensed_Cy),np.max(self.sensed_Cy),np.min(self.sensed_Cy)
        #print "address of min" ,np.argmin(self.sensed_C), len(self.sensed_C),self.sensed_C[len(self.sensed_C)-10:]
        #print "shape(self.sensed_C),", np.shape(self.sensed_C)
        #print "shape(self.xpos),", np.shape(self.xpos)
        #print "popsize", self.params.loc["current_population_size"].values[0]
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        plt.scatter(self.xpos, self.ypos, c=self.sensed_C,s=8, cmap='Accent',edgecolor="none")       
        #plt.text(0.0,0.95, "mean max min "+str(np.mean(self.sensed_C) )+", "+ str(np.max(self.sensed_C))+", "+ str(np.min(self.sensed_C)) ,color="r")
        plt.colorbar()
        plt.savefig(destfold+"sensed_C.png")
        
        
        #fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        #plt.scatter(self.ygrid, np.mean(self.sensed_C,axis=0),'o',markeredgecolor="none")
        #plt.ylabel("mean concentration")
        #plt.xlabel("yposition")
        #plt.savefig(destfold+"sensed_meanC_along_y.png")
        
        
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        plt.scatter(self.xpos, self.ypos, c=self.sensed_Cx,s=8, cmap='Accent',edgecolor="none")
        plt.colorbar()
        plt.savefig(destfold+"sensed_Cx.png")
        
        fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))   
        plt.scatter(self.xpos, self.ypos, c=self.sensed_Cy,s=8, cmap='Accent',edgecolor="none")
        plt.colorbar()
        plt.savefig(destfold+"sensed_Cy.png")
        
        
        #sensed_C=np.ravel(sensed_C_all[0])