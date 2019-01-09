#!/usr/bin/env python
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
#from matplotlib.colors import colorConverter
import glob
import os
import time
import pandas as pd
from scipy import stats
option2=0
SS=6
r=24



path = '/projectnb/qedksk/ashishge/data/new_sim/bulge_slope_from_track_height/gdnoise/try3/'
specific_path="*/*/*/"
destpath=path[:-1]+"_slope_analysis/"

#path = '/Users/ashish/Downloads/bulge_trackheight_test/chsum1 h0.5/'
#specific_path=""
#destpath='/Users/ashish/Downloads/bulge_trackheight_test_slope_analysis/chsum1 h0.5/'
if not os.path.exists(destpath): os.makedirs(destpath) 

gz_flag=0
img_shift_flag=0
print_txt_flag=1


color_scheme=0


def plot_figure(ca,cb,destfile_name,color_scheme,tick_flag=None):
    l= len(ca)
    b=len(ca.T)
    
    fig3, ax = plt.subplots(nrows=1, ncols=1, figsize=(b/100.0,l/100.0))   
    img = np.zeros((l,b,3))

    if color_scheme==3: #red_blue or blue red is next option
                            ## since red chiral and blue  nonchiral, youm might need to use other color scheme
        img[:,:,0]=cb[:,:]
        img[:,:,2]=ca[:,:]     
            
            
    if color_scheme==2: #blue_red
        img[:,:,0]=ca[:,:]
        img[:,:,2]=cb[:,:]     
        
    elif color_scheme==1: #blue yellow
        img[:,:,1]=0.9*cb[:,:]
        img[:,:,0]=0.9*cb[:,:]
        img[:,:,2]=ca[:,:] 
        
    elif color_scheme==0: #red green
        img[:,:,0]=ca[:,:]
        img[:,:,1]=cb[:,:]
    img3 = ax.imshow(img,aspect="equal")  
    if tick_flag==None or tick_flag==0:   
        ax.set_xticks([])
        ax.set_yticks([])
    
    plt.savefig(destfile_name) 


def plot_fractions(ca,c,destfile_name):
    
    mean_c=np.mean(c,axis=0)
    mean_ca=np.mean(ca,axis=0)
    fa=mean_ca[np.where(mean_c>0.1)]/mean_c[np.where(mean_c>0.1)]
    
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.set_title('fraction of a ')
    plt.plot(fa,'o')
    plt.savefig(destfile_name)
    
    
def plot_height(h,destfile_name):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.set_title('h(x) ')
    plt.plot(h,'b-')
    plt.savefig(destfile_name)

def plot_data_errorbar(xval,yval,xlabel,ylabel,destfile_name,yerror=0):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    if yerror==0:
        plt.plot(xval,yval,'o')
    else:
        plt.errorbar(xval,yval,yerr=yerror)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(xlabel)
    plt.savefig(destfile_name)

def plot_height_and_fit(h,fit1,x1,fit2,x2,destfile_name):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.set_title('h(x) ')
    plt.plot(h,'b-')
    plt.plot(x1,fit1,'r-')
    plt.plot(x2,fit2,'r-')
    plt.savefig(destfile_name)
    
def fit_bulge_slopes(harr,hfit_file_name):
        ########### this method finds slopes by shifting hfront pos such that xmin is at 0  ################
        xmin=np.argmin(harr)       
        harr=np.roll(harr,-xmin)       
        xmax=np.argmax(harr)
        

        xval1=np.arange(10,xmax)
        xval2=np.arange(xmax,l-10)
        
        if min(len(xval1),len(xval2))>10:
        
            slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress( xval1,harr[xval1]  )
            slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress( xval2,harr[xval2]   )
            print slope1,slope2
            fit1=slope1*xval1+intercept1
            fit2=slope2*xval2+intercept2
            
            if hfit_file_name!='n':
                plot_height_and_fit(harr,fit1,xval1%l,fit2,xval2%l,hfit_file_name)
            return (slope1, slope2, r_value1,r_value2,std_err1,std_err2)
        else :
            print "too early", xmax,l-xmax
        
        ########### this method finds slopes wihtout shifting hfront pos and wrapping every time, not as easy to read ################                    
        #xmin=np.argmin(hfrontpos[lexp])
        #xmax=np.argmax(hfrontpos[lexp])
        #if xmin<xmax:
        #    xval1=np.arange(xmin,xmax)## positive slop
        #    xval2=np.arange(xmax,xmin+l) ## -ve slope
        #    
        #else:
        #    xval1=np.arange(xmin,xmax+l) ##
        #    xval2=np.arange(xmax,xmin) ##
        #
        #slope1, intercept1, r_value1, p_value, std_err = stats.linregress( xval1,np.take(hfrontpos[lexp],xval1,mode='wrap')   )
        #slope2, intercept2, r_value1, p_value, std_err = stats.linregress( xval2,np.take(hfrontpos[lexp],xval2,mode='wrap')   )
        #
        #fit1=slope1*xval1+intercept1
        #fit2=slope2*xval2+intercept2
        #hfit_file_name=image_file_name[:-7]+image_file_name_subset.replace("ca","hfit")
        #plot_height_and_fit(hfrontpos[lexp],fit1,xval1%l,fit2,xval2%l,hfit_file_name)



for ch in glob.glob(path+specific_path):	
    print ch
    fold=ch
    
    
    #if gz_flag==0:      
    for ca_file  in glob.glob(ch+'ca.txt'):   
        
        ca = np.array(pd.read_csv(ca_file,sep=' ').values   ) 
        
        ca_file_name=ca_file.replace(ch,"")        
        c_file=ch+ca_file_name.replace("a","")
        c = np.array(pd.read_csv(c_file,sep=' ').values   ) 
    #else:
    #    for ca_file  in glob.glob(ch+'ca*.gz'):
    #   
    #        ca= pd.read_csv(ca_file,sep=' ').values
    #        ca_file_name=ca_file_name.replace(ch ,"")
    #        c_file=ch+ca_file_name.replace("a","")
    #        c= pd.read_csv(c_file,sep=' ').values
    
        cb=c-ca             
        l= len(ca)
        b=len(ca.T)
        frontpos_file=ch+ca_file_name.replace("ca","frontpos")
        hfrontpos=np.array(pd.read_csv(frontpos_file,sep=' ').values   ) 
        

        destfold=fold.replace(path,destpath)   
        if not os.path.exists(destfold): os.makedirs(destfold)     
        image_file_name=ca_file.replace(path,destpath)
        #if gz_flag==0:  
        image_file_name=image_file_name.replace(".txt",".png")
        #else:
        #    image_file_name=image_file_name.replace(".gz",".png")
     
        if img_shift_flag==1:
            print "shifting image!!!"
            img_shift=l/2
            ca=np.roll(ca,img_shift,axis=0)
            c=np.roll(c,img_shift,axis=0)
            cb=np.roll(cb,img_shift,axis=0)
        plot_figure(ca,cb,image_file_name,color_scheme)
        
        image_file_name_subset=image_file_name[-7:]
        hx_file_name=image_file_name[:-7]+image_file_name_subset.replace("ca","height_x")
        h_x=np.sum(c,axis=1)
        plot_height(h_x,hx_file_name)
        lbeg=np.where(hfrontpos[:,0]>0)[0][0]+1
        print lbeg
        lexp=np.where(hfrontpos[:,0]>0)[0][-1]-1
        print lexp
        hfrontpos_file_name=image_file_name[:-7]+image_file_name_subset.replace("ca","height_frontpos")
        plot_height(hfrontpos[lexp],hfrontpos_file_name)
                
        start_time = time.time()        
        n_fits=10
        
        slope1_list=[]
        slope2_list=[]
        r_value1_list=[]
        r_value2_list=[]
        std_err1_list=[]
        std_err2_list=[]
                
        tlist=np.linspace(lbeg+lexp/2,lexp,n_fits,dtype=int)
        for idx,t in enumerate(tlist):
            hfit_file_name=image_file_name[:-7]+image_file_name_subset.replace("ca","hfit"+str(idx))
            slope1, slope2, r_value1,r_value2,std_err1,std_err2=fit_bulge_slopes(hfrontpos[t],hfit_file_name)
            
            slope1_list.append(slope1)
            slope2_list.append(slope2)
            r_value1_list.append(r_value1)
            r_value2_list.append(r_value2)
            std_err1_list.append(std_err1)
            std_err2_list.append(std_err2)

        end_time = time.time()
        print "time taken for all the loops in secs", (end_time-start_time)   
                
        slope1_file_name=image_file_name[:-7]+image_file_name_subset.replace("ca","slope1")      
        plot_data_errorbar(tlist,slope1_list,"time","slope1",slope1_file_name,yerror=std_err1_list)
        slope2_file_name=image_file_name[:-7]+image_file_name_subset.replace("ca","slope2")
        plot_data_errorbar(tlist,slope2_list,"time","slope2",slope2_file_name,yerror=std_err2_list)
        
        r1_file_name=image_file_name[:-7]+image_file_name_subset.replace("ca","R1")
        plot_data_errorbar(tlist,np.square(r_value1_list),"time","$R^2$",r1_file_name)
        r2_file_name=image_file_name[:-7]+image_file_name_subset.replace("ca","R2")
        plot_data_errorbar(tlist,np.square(r_value2_list),"time","$R^2$",r2_file_name)
         
                                          
                                                                                                                                      
        if print_txt_flag==1:
            
            np.savetxt(slope1_file_name.replace(".png",".txt"),slope1_list) 
            np.savetxt(slope2_file_name.replace(".png",".txt"),slope2_list) 
            np.savetxt(slope1_file_name.replace("slope1.png","std_err1.txt"),std_err1_list) 
            np.savetxt(slope2_file_name.replace("slope2.png","std_err2.txt"),std_err2_list)        
            np.savetxt(r1_file_name.replace(".png",".txt"),r_value1_list) 
            np.savetxt(r2_file_name.replace(".png",".txt"),r_value2_list) 
            
            np.savetxt(slope1_file_name.replace("slope1.png","tlist.txt"),tlist) 
            np.savetxt(hx_file_name.replace(".png",".txt"),h_x)
    
    
    plt.close("all")
    
    


        
        




plt.close("all")
