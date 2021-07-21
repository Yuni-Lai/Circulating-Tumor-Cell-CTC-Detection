#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 20:32:22 2020

@author: Yuni
"""

import sys
import os
import pickle
import matplotlib.pyplot as plt
from skimage import io
import fnmatch
import pickle
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.stats import *
from scipy.stats import gaussian_kde
import seaborn as sns
import itertools
import scipy
import math
from skimage import measure
import cv2
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
Blue_dia_all=[]
Cell_dia_all=[]
Yel_int_all=[]
Blue_int_all=[]

result_Path="/dataserver145/image/Yuni/CTC/01.1_patients_samples/01_batch/02_segmented/"
postfix=result_Path.split('/dataserver145/image/Yuni/CTC/')[-1]
postfix=postfix.split('02_segmented/')[0]
files=os.listdir(result_Path)

for ff in files:
    if fnmatch.fnmatch(ff, '*.pkl'):
        f=open(f'{result_Path}/{ff}','rb')
        Blue_dia,Cell_dia,Yel_int,Blue_int,seg_all1,All_num,CTC_num,seg_CTC,nmask1,cmask1,img_c,Added_image=pickle.load(f)
        f.close()
        Blue_dia_all.append(Blue_dia)
        Cell_dia_all.append(Cell_dia)
        Yel_int_all.append(Yel_int)
        Blue_int_all.append(Blue_int)

PATH='/dataserver145/image/Yuni/CTC/02_dsitribution_files/'+postfix
if os.path.exists(PATH):
    print('file exit: ',PATH)
if not os.path.exists(PATH):
        os.makedirs(PATH)
os.chdir(PATH)
f = open('Blue_dia_all.pkl', 'wb')
pickle.dump(Blue_dia_all, f)
f.close()
f = open('Cell_dia_all.pkl', 'wb')
pickle.dump(Cell_dia_all, f)
f.close()
f = open('Yel_int_all.pkl', 'wb')
pickle.dump(Yel_int_all, f)
f.close()
f = open('Blue_int_all.pkl', 'wb')
pickle.dump(Blue_int_all, f)
f.close()


# os.chdir('/home/Yuni/CTC/variables/healthy_01/')
# f = open('Blue_dia_all.pckl', 'rb')
# Blue_dia_all = pickle.load(f)
# f.close()
# f = open('Cell_dia_all.pckl', 'rb')
# Cell_dia_all = pickle.load(f)
# f.close()
# f = open('Blue_int_all.pckl', 'rb')
# Blue_int_all = pickle.load(f)
# f.close()
# #os.chdir('/home/Yuni/CTC/variables/34')
# f = open('Yel_int_all.pckl', 'rb')
# Yel_int_all = pickle.load(f)
# f.close()
matplotlib.use('Qt5Agg',force=True)
plt.figure()
plt.title("Blue diameter")
sns.set_style('darkgrid')
Blue_dia_all1= list(itertools.chain.from_iterable(Blue_dia_all))
Blue_dia_all1=[x for x in Blue_dia_all1 if x<20]
sns.distplot(Blue_dia_all1)#,hist_kws={'color':'green'}
plt.figure()
plt.title("Blue intensity")
sns.set_style('darkgrid')
Blue_int_all1= list(itertools.chain.from_iterable(Blue_int_all))
sns.distplot(Blue_int_all1,hist_kws={'color':'green'})#
plt.figure()
plt.title("Yellow intensity")
sns.set_style('darkgrid')
Yel_int_all1= list(itertools.chain.from_iterable(Yel_int_all))
sns.distplot(Yel_int_all1,hist_kws={'color':'yellow'})#
plt.figure()
plt.title("Cell diameter")
sns.set_style('darkgrid')
Cell_dia_all1= list(itertools.chain.from_iterable(Cell_dia_all))
Cell_dia_all1=[x for x in Cell_dia_all1 if x<25]
sns.distplot(Cell_dia_all1,hist_kws={'color':'yellow'})#


#=======cutoff=============================================
def drawcutoff(Blue_dia_all1,direction):#direction={0,1},0=left,1=right
    Blue_dia_array=np.array(Blue_dia_all1)
    p1, p99 = np.percentile(Blue_dia_array, (1, 99))
    Blue_dia_all_part=Blue_dia_array[(p1 < Blue_dia_array) & (Blue_dia_array < p99)]
    ax=sns.distplot(Blue_dia_all_part)#,fit=norm
    
    #Get the data from the KDE line
    xdata, ydata = ax.get_lines()[0].get_data()
    
    cdf = scipy.integrate.cumtrapz(ydata, xdata, dx=1, initial=0)
    if direction==1:
        index=(np.where(cdf>=0.95))[0][0]
        ydata[index-1]
        cut=xdata[index-1]
    if direction==0:
        index=(np.where(cdf<=0.05))[0][-1]
        ydata[index-1]
        cut=xdata[index-1]
    value=cut
    #plt.plot(xdata,cdf)
    #Find the closest point on the curve
    idx = (np.abs(xdata-value)).argmin()
    #Interpolate to get a better estimate
    p = np.interp(value,xdata[idx:idx+2],ydata[idx:idx+2])
    print('Point on PDF for X = {} is: {}'.format(value,p))
    #Plot the line
    ax.vlines(value, 0, p ,colors='r')
    xy=(value,p)
    ax.annotate("P<0.05: (%s)" % round(xy[0], 2), xy=xy, xytext=(-20, 10), textcoords='offset points')
    #Find the probability for an interval of one (e.g. between 20 and 100)
    if direction==1:
        ecart = 100
        idx = (np.abs(xdata-value)).argmin()
        idx2 = (np.abs(xdata-(value+ecart))).argmin()
        pr = cdf[idx2] - cdf[idx] # Error here see old code, need to define idx_ 
        print('Probability of X <{},{}> is: {}'.format(value,value+ecart,pr))
        
        # Fill the area 
        plt.fill_between(xdata,ydata, where = (xdata>=value) & (xdata<=value+ecart), color='g')
        #plt.pause(60)
        #plt.show(block=True)
        #plt.close()
    if direction==0:
        ecart = -100
        idx = (np.abs(xdata-value)).argmin()
        idx2 = (np.abs(xdata-(value+ecart))).argmin()
        pr = cdf[idx2] - cdf[idx] # Error here see old code, need to define idx_ 
        print('Probability of X <{},{}> is: {}'.format(value,value+ecart,pr))
        
        # Fill the area 
        plt.fill_between(xdata,ydata, where = (xdata<=value) & (xdata>=value+ecart), color='g')
        #plt.pause(60)
        
        #plt.close()
    return(cut,xdata, ydata,cdf)

plt.figure()
plt.title("Blue diameter")
Blue_dia_cut,Blue_dia_xdata, Blue_dia_ydata,Blue_dia_cdf=drawcutoff(Blue_dia_all1,1)

plt.figure()
plt.title("Blue intensity")
Blue_int_cut,Blue_int_xdata, Blue_int_ydata,Blue_int_cdf=drawcutoff(Blue_int_all1,0)

plt.figure()
plt.title("Yellow intensity")
Yel_int_cut,Yel_int_xdata, Yel_int_ydata,Yel_int_cdf=drawcutoff(Yel_int_all1,0)

plt.figure()
plt.title("Cell diameter")
Cell_dia_cut,Cell_dia_xdata, Cell_dia_ydata,Cell_dia_cdf=drawcutoff(Cell_dia_all1,1)

#plt.show(block=True)

#--------????------------------------
os.chdir('/home/Yuni/CTC/02_dsitribution_files/Tumor02/')
f = open('Blue_dia_all.pckl', 'rb')
Blue_dia_all = pickle.load(f)
f.close()
f = open('Cell_dia_all.pckl', 'rb')
Cell_dia_all = pickle.load(f)
f.close()
f = open('Yel_int_all.pckl', 'rb')
Yel_int_all = pickle.load(f)
f.close()
f = open('Blue_int_all.pckl', 'rb')
Blue_int_all = pickle.load(f)
f.close()


#-----------------??????-------------------------
def Normal(x,mu,sigma):#????????????
    
    return np.exp(-(x-mu)**2/(2*sigma**2))/(np.sqrt(2*np.pi)*sigma)

def GMM2(data,n,Mu,SigmaSquare,Alpha):
    N=data.shape[0]#?????
    #n:????   
    np.random.seed(1)
    
    i=0#????
    
    while(i<=n):#?EM?????????
        i+=1
        #Expectation--------
        if SigmaSquare[0][0]<=0.001:
            SigmaSquare[0][0]=SigmaSquare[0][0]+0.1
        if SigmaSquare[0][1]<=0.001:
            SigmaSquare[0][1]=SigmaSquare[0][1]+0.1
        gauss1=Normal(data,Mu[0][0],np.sqrt(SigmaSquare[0][0]))#?????
        gauss2=Normal(data,Mu[0][1],np.sqrt(SigmaSquare[0][1]))#?????
        
        Gamma1=Alpha[0][0]*gauss1
        Gamma2=Alpha[0][1]*gauss2
        M=Gamma1+Gamma2+0.00001

    
        #Gamma=np.concatenate((Gamma1/m,Gamma2/m),axis=1) ??(j,k)??j??????k??????,???????????
    
        #Maximization--------
        #??SigmaSquare
        SigmaSquare[0][0]=np.dot((Gamma1/M).T,(data-Mu[0][0])**2)/np.sum(Gamma1/M)
        SigmaSquare[0][1]=np.dot((Gamma2/M).T,(data-Mu[0][1])**2)/np.sum(Gamma2/M)
    
        #??Mu  
        Mu[0][0]=np.dot((Gamma1/M).T,data)/np.sum(Gamma1/M)
        Mu[0][1]=np.dot((Gamma2/M).T,data)/np.sum(Gamma2/M)
    
        #??Alpha
        Alpha[0][0]=np.sum(Gamma1/M)/N
        Alpha[0][1]=np.sum(Gamma2/M)/N
        
        if(i%10==0):
            print ("The ",i," iteration:")
            print ("Mu:",Mu)
            print ("Sigma:",np.sqrt(SigmaSquare))
            print ("Alpha",Alpha)
    return(Mu,SigmaSquare,Alpha)


def GMM3(data,n,Mu,SigmaSquare,Alpha):

    N=data.shape[0]#?????
    #n:????
    np.random.seed(1)

    i=0#????
    
    while(i<=n):#?EM?????????
        
        i+=1
        
        #Expectation
        if SigmaSquare[0][0]<=0.001:
            SigmaSquare[0][0]=SigmaSquare[0][0]+0.03
        if SigmaSquare[0][1]<=0.001:
            SigmaSquare[0][1]=SigmaSquare[0][1]+0.03
        if SigmaSquare[0][2]<=0.001:
            SigmaSquare[0][2]=SigmaSquare[0][2]+0.03
        gauss1=Normal(data,Mu[0][0],np.sqrt(SigmaSquare[0][0]))#?????
        gauss2=Normal(data,Mu[0][1],np.sqrt(SigmaSquare[0][1]))#?????
        gauss3=Normal(data,Mu[0][2],np.sqrt(SigmaSquare[0][2]))#?????
        
        Gamma1=Alpha[0][0]*gauss1
        Gamma2=Alpha[0][1]*gauss2
        Gamma3=Alpha[0][2]*gauss3
    
        M=Gamma1+Gamma2+Gamma3+0.00001
    
        #Gamma=np.concatenate((Gamma1/m,Gamma2/m),axis=1) ??(j,k)??j??????k??????,???????????
    
        #Maximization
        
        #??SigmaSquare
        
        SigmaSquare[0][0]=np.dot((Gamma1/M).T,(data-Mu[0][0])**2)/np.sum(Gamma1/M)
        
        SigmaSquare[0][1]=np.dot((Gamma2/M).T,(data-Mu[0][1])**2)/np.sum(Gamma2/M)
        SigmaSquare[0][2]=np.dot((Gamma3/M).T,(data-Mu[0][2])**2)/np.sum(Gamma3/M)
    
        #??Mu       
    
        Mu[0][0]=np.dot((Gamma1/M).T,data)/np.sum(Gamma1/M)
        Mu[0][1]=np.dot((Gamma2/M).T,data)/np.sum(Gamma2/M)
        Mu[0][2]=np.dot((Gamma3/M).T,data)/np.sum(Gamma3/M)
    
        #??Alpha
    
        Alpha[0][0]=np.sum(Gamma1/M)/N  
        Alpha[0][1]=np.sum(Gamma2/M)/N
        Alpha[0][2]=np.sum(Gamma3/M)/N
        
        if(i%10==0):
            print ("The ",i," iteration:")
            print ("Mu:",Mu)
            print ("Sigma:",np.sqrt(SigmaSquare))
            print ("Alpha",Alpha)
    return(Mu,SigmaSquare,Alpha)


def gaussian(sigma, x, u):
	y = np.exp(-(x - u) ** 2 / (2 * sigma ** 2)) / (sigma * math.sqrt(2 * math.pi))
	return y

#???????-----------
sns.set_style('darkgrid')
Blue_dia_all1= list(itertools.chain.from_iterable(Blue_dia_all))
Blue_dia_all_array=np.array(Blue_dia_all1)
Blue_dia_all_array.shape=Blue_dia_all_array.shape[0],1
data=Blue_dia_all_array#??????,N?1?
data=data[np.where(data < 25)]
data=data[np.where(data > 5)]
Mu=np.array([[5.5,12.5]])#?????
SigmaSquare=np.array([[2,8]]) #?????Sigma??
Alpha=np.array([[0.5,0.5]])  
n=300
Mu,SigmaSquare,Alpha=GMM2(data,n,Mu,SigmaSquare,Alpha)
x = np.linspace(data.min(), data.max(), 1000)
plt.figure()
plt.title("Blue diameter_GMM2")
sns.set_style('darkgrid')
sns.distplot(list(data))#,hist_kws={'color':'green'}
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][0]), x, Mu[0][0])*Alpha[0][0], "g-", linewidth=1)
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][1]), x, Mu[0][1])*Alpha[0][1], "r-", linewidth=1)
a=Mu[0][1]
b=gaussian(np.sqrt(SigmaSquare[0][1]), Mu[0][1], Mu[0][1])*Alpha[0][1]
plt.text(a,b , (a,b),ha='center', va='bottom', fontsize=10)
plt.plot([a, a,], [0, b,], 'k--', linewidth=2.5)
plt.show()
#plt.show(block=True)
Blue_dia_Mu,Blue_dia_SigmaSquare,Blue_dia_Alpha=Mu,SigmaSquare,Alpha

##????-----------
#plt.figure()
#plt.title("Blue intensity")
#sns.set_style('darkgrid')
Blue_int_all1= list(itertools.chain.from_iterable(Blue_int_all))
#sns.distplot(Blue_int_all1,hist_kws={'color':'green'})#
all_array=np.array(Blue_int_all1)
all_array.shape=all_array.shape[0],1
data=all_array#??????,N?1?
#data=data[np.where(data < 25)]
#data=data[np.where(data > 5)]
Mu=np.array([[10,230]])#?????
SigmaSquare=np.array([[2,8]]) #?????Sigma??
Alpha=np.array([[0.5,0.5]])  
n=300
Mu,SigmaSquare,Alpha=GMM2(data,n,Mu,SigmaSquare,Alpha)
x = np.linspace(data.min()-50, data.max()+50, 1000)
plt.figure()
plt.title("Blue intensity_GMM2")
sns.set_style('darkgrid')
sns.distplot(list(data))#,hist_kws={'color':'green'}
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][0]), x, Mu[0][0])*Alpha[0][0], "g-", linewidth=1,label='G1')
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][1]), x, Mu[0][1])*Alpha[0][1], "r-", linewidth=1,label='G2')
a=Mu[0][1]
b=gaussian(np.sqrt(SigmaSquare[0][1]), Mu[0][1], Mu[0][1])*Alpha[0][1]
plt.text(a,b , (a,b),ha='center', va='bottom', fontsize=10)
plt.plot([a, a,], [0, b,], 'k--', linewidth=2.5)
plt.legend(loc='upper right')
plt.show()
#plt.show(block=True)
Blue_int_Mu,Blue_int_SigmaSquare,Blue_int_Alpha=Mu,SigmaSquare,Alpha

#????-----------
#plt.figure()
#plt.title("Yellow intensity")
#sns.set_style('darkgrid')
Yel_int_all1= list(itertools.chain.from_iterable(Yel_int_all))
#sns.distplot(Yel_int_all1,hist_kws={'color':'yellow'})#
all_array=np.array(Yel_int_all1)
all_array.shape=all_array.shape[0],1
data=all_array#??????,N?1?
#data=data[np.where(data < 25)]
#data=data[np.where(data > 5)]
Mu=np.array([[0.1,0.5]])#?????
SigmaSquare=np.array([[10.5,5.5]]) #?????Sigma??
Alpha=np.array([[0.6,0.4]])  
n=100
Mu,SigmaSquare,Alpha=GMM2(data,n,Mu,SigmaSquare,Alpha)
x = np.linspace(data.min(), data.max(), 1000)
plt.figure()
plt.title("Yellow intensity_GMM2")
sns.set_style('darkgrid')
sns.distplot(list(data))#,hist_kws={'color':'green'}
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][0]), x, Mu[0][0])*Alpha[0][0], "g-", linewidth=1,label='G1')
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][1]), x, Mu[0][1])*Alpha[0][1], "r-", linewidth=1,label='G2')
a=Mu[0][0]
b=gaussian(np.sqrt(SigmaSquare[0][0]), Mu[0][0], Mu[0][0])*Alpha[0][0]
plt.text(a,b , (a,b),ha='center', va='bottom', fontsize=10)
plt.plot([a, a,], [0, b,], 'k--', linewidth=2.5)
plt.legend(loc='upper right')
plt.show()
#plt.show(block=True)
Yel_int_Mu,Yel_int_SigmaSquare,Yel_int_Alpha=Mu,SigmaSquare,Alpha


#??????-----------
#plt.figure()
#plt.title("Cell diameter")
#sns.set_style('darkgrid')
Cell_dia_all1= list(itertools.chain.from_iterable(Cell_dia_all))
#sns.distplot(Cell_dia_all1,hist_kws={'color':'yellow'})#
all_array=np.array(Cell_dia_all1)
all_array.shape=all_array.shape[0],1
data=all_array#??????,N?1?
data=data[np.where(data < 25)]
data=data[np.where(data > 3)]

Mu=np.array([[5.5,10.5,12.5]])#?????
SigmaSquare=np.array([[2.1,6.5,8.1]]) #?????Sigma??
Alpha=np.array([[0.33,0.33,0.34]])#????????????(????0,???1)
n=500
Mu,SigmaSquare,Alpha=GMM3(data,n,Mu,SigmaSquare,Alpha)#data,iterate number
x = np.linspace(data.min()-10, data.max()+10, 1000)
plt.figure()
plt.title("Cell diameter_GMM3")
sns.set_style('darkgrid')
sns.distplot(list(data),kde=True)#,hist_kws={'color':'green'}
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][0]), x, Mu[0][0])*Alpha[0][0], "g-", linewidth=1,label='G1')
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][1]), x, Mu[0][1])*Alpha[0][1], "r-", linewidth=1,label='G2')
plt.plot(x, gaussian(np.sqrt(SigmaSquare[0][2]), x, Mu[0][2])*Alpha[0][2], color='gold',linestyle='-', linewidth=1,label='G3')
a=Mu[0][2]
b=gaussian(np.sqrt(SigmaSquare[0][2]), Mu[0][2], Mu[0][2])*Alpha[0][2]
plt.text(a,b , (a,b),ha='center', va='bottom', fontsize=10)
plt.plot([a, a,], [0, b,], 'k--', linewidth=2.5)
plt.legend(loc='upper right')
plt.show()
plt.show(block=True)
Cell_dia_Mu,Cell_dia_SigmaSquare,Cell_dia_Alpha=Mu,SigmaSquare,Alpha



os.makedirs(f'/home/Yuni/CTC/variables/GMMparameters')
os.chdir(f'/home/Yuni/CTC/variables/GMMparameters')
f = open('GMMparameters.pckl', 'wb')
pickle.dump([Cell_dia_Mu,Cell_dia_SigmaSquare,Cell_dia_Alpha,Yel_int_Mu,Yel_int_SigmaSquare,Yel_int_Alpha,Blue_int_Mu,Blue_int_SigmaSquare,Blue_int_Alpha,Blue_dia_Mu,Blue_dia_SigmaSquare,Blue_dia_Alpha], f)
f.close()

os.makedirs(f'/home/Yuni/CTC/variables/P_parameters')
os.chdir(f'/home/Yuni/CTC/variables/P_parameters')
f = open('P_parameters.pckl', 'wb')
pickle.dump([Blue_dia_cut,Blue_dia_xdata, Blue_dia_ydata,Blue_dia_cdf,Blue_int_cut,Blue_int_xdata, Blue_int_ydata,Blue_int_cdf,Yel_int_cut,Yel_int_xdata, Yel_int_ydata,Yel_int_cdf,Cell_dia_cut,Cell_dia_xdata, Cell_dia_ydata,Cell_dia_cdf], f)
f.close()

#--------write excel--------------------------------
import xlrd
import xlwt
from xlutils.copy import copy
 
def write_excel_xls(path, sheet_name, value):
    index = len(value)  # ???????????
    workbook = xlwt.Workbook()  # ???????
    sheet = workbook.add_sheet(sheet_name)  # ???????????
    for i in range(0, index):
        for j in range(0, len(value[i])):
            sheet.write(i, j, value[i][j])  # ????????(??????)
    workbook.save(path)  # ?????
    print("write in xls success!")
 
 
def write_excel_xls_append(path, value):
    index = len(value)  # ???????????
    workbook = xlrd.open_workbook(path)  # ?????
    sheets = workbook.sheet_names()  # ???????????
    worksheet = workbook.sheet_by_name(sheets[0])  # ??????????????????
    rows_old = worksheet.nrows  # ??????????????
    new_workbook = copy(workbook)  # ?xlrd???????xlwt??
    new_worksheet = new_workbook.get_sheet(0)  # ???????????????
    for i in range(0, index):
        for j in range(0, len(value[i])):
            new_worksheet.write(i+rows_old, j, value[i][j])  # ??????,????i+rows_old?????
    new_workbook.save(path)  # ?????
    print("append xls success!")
 
 
def read_excel_xls(path):
    workbook = xlrd.open_workbook(path)  # ?????
    sheets = workbook.sheet_names()  # ???????????
    worksheet = workbook.sheet_by_name(sheets[0])  # ??????????????????
    for i in range(0, worksheet.nrows):
        for j in range(0, worksheet.ncols):
            print(worksheet.cell_value(i, j), "\t", end="")  # ????????
        print()




#--CTC_estimate-------------
def CTC_estimate(cmask2,nmask2,info_yellow,info_yellow_backg,seg_all2,img_nucle_rescale,Added_image,cell_minus,n,num,xls_path):
    info_cell=measure.regionprops(cmask2,intensity_image=cell_minus)
    info_nuccc=measure.regionprops(nmask2,intensity_image=img_nucle_rescale)
    Blue_int=np.zeros(len(info_nuccc))
    Cell_int=np.zeros(len(info_nuccc))
    Blue_dia=np.zeros(len(info_nuccc))
    Cell_dia=np.zeros(len(info_cell))
    Yel_int=np.zeros(len(info_yellow))
    #intensity_Y=np.zeros(len(info_yellow))
    for i in range(0,len(info_nuccc)-1):
        Blue_int[i]=info_nuccc[i]['mean_intensity']
        Yel_int[i]=info_yellow[i]['mean_intensity']
        Cell_int[i]=info_cell[i]['mean_intensity']
        Blue_dia[i]=info_nuccc[i]['equivalent_diameter']*0.172
        Cell_dia[i]=info_cell[i]['equivalent_diameter']*0.172
        #intensity_Y[i]= info_yellow[i]['mean_intensity']
    #p70_B = np.percentile(Blue_int, 70)
    #p40_Y,p50_Y = np.percentile(Yel_int, (40,50))
    #p10_C=np.percentile(Cell_int, 30)
   # p10_C
#    plt.figure()
#    plt.title("Blue intensity")
#    sns.set_style('darkgrid')
#    sns.distplot(Blue_int)#,hist_kws={'color':'green'}
#    plt.figure()
#    plt.title("Blue diameter")
#    sns.set_style('darkgrid')
#    sns.distplot(Blue_dia,hist_kws={'color':'green'})#
#    plt.figure()
#    plt.title("Yellow intensity")
#    sns.set_style('darkgrid')
#    sns.distplot(Yel_int,hist_kws={'color':'yellow'})#
#    plt.figure()
#    plt.title("Yellow relevent intensity")
#    sns.set_style('darkgrid')
#    sns.distplot(intensity_Y,hist_kws={'color':'yellow'})#
    CTC=[]
    #plt.figure()
    #plt.imshow(seg_all1)
    CTC_mask_strong=np.zeros(nmask2.shape, dtype=np.uint8) 
    CTC_mask_very_strong=np.zeros(nmask2.shape, dtype=np.uint8) 
    CTC_mask_extreme=np.zeros(nmask2.shape, dtype=np.uint8) 
    for i in range(0,len(info_cell)): 
        #i=152 
        result_a_cell=[]
        
        info_a_cell=[]
        info_a_cell.append(i)
        info_a_cell.append(Cell_dia[i])
        info_a_cell.append(Blue_dia[i])
        info_a_cell.append(Yel_int[i])
        
        #Blue_dia-----------------
        M2_Blue_dia=gaussian(np.sqrt(Blue_dia_SigmaSquare[0][1]), Blue_dia[i], Blue_dia_Mu[0][1])*Blue_dia_Alpha[0][1]
        M1_Blue_dia=gaussian(np.sqrt(Blue_dia_SigmaSquare[0][0]), Blue_dia[i], Blue_dia_Mu[0][0])*Blue_dia_Alpha[0][0]
        
        if M1_Blue_dia==0:
            M1_Blue_dia=M1_Blue_dia+0.00001
        BayF_Blue_dia=math.exp(math.log(M2_Blue_dia)-math.log(M1_Blue_dia))
        #Cell_dia-----------
        M2_Cell_dia=gaussian(np.sqrt(Cell_dia_SigmaSquare[0][2]), Cell_dia[i], Cell_dia_Mu[0][2])*Cell_dia_Alpha[0][2]
        M1_Cell_dia=gaussian(np.sqrt(Cell_dia_SigmaSquare[0][0]), Cell_dia[i], Cell_dia_Mu[0][0])*Cell_dia_Alpha[0][0]+gaussian(np.sqrt(Cell_dia_SigmaSquare[0][1]), Cell_dia[i], Cell_dia_Mu[0][1])*Cell_dia_Alpha[0][1]
        if M1_Cell_dia==0:
            M1_Cell_dia=M1_Cell_dia+0.00001
        BayF_Cell_dia=math.exp(math.log(M2_Cell_dia)-math.log(M1_Cell_dia))
        #Yel_int-------------
        M2_Yel_int=gaussian(np.sqrt(Yel_int_SigmaSquare[0][0]), Yel_int[i], Yel_int_Mu[0][0])*Yel_int_Alpha[0][0]
        M1_Yel_int=gaussian(np.sqrt(Yel_int_SigmaSquare[0][1]), Yel_int[i], Yel_int_Mu[0][1])*Yel_int_Alpha[0][1]
        
        if M1_Yel_int==0:
            M1_Yel_int=M1_Yel_int+0.00001
        BayF_Yel_int=math.exp(math.log(M2_Yel_int)-math.log(M1_Yel_int))
        #ALl-------------------
        BayF_ALL=BayF_Blue_dia*BayF_Cell_dia*BayF_Yel_int#(BayF_Blue_dia+BayF_Cell_dia+BayF_Yel_int)/3
        #-------------
        BayF_a_cell=[]
        BayF_a_cell.append("Bayes_Factor")
        BayF_a_cell.append(BayF_Blue_dia)
        BayF_a_cell.append(BayF_Cell_dia)
        BayF_a_cell.append(BayF_Yel_int)
        BayF_a_cell.append(BayF_ALL)
        
        
        Blue_dia_p=1-Blue_dia_cdf[(np.abs(Blue_dia_xdata-Blue_dia[i])).argmin()]
        Cell_dia_p=1-Cell_dia_cdf[(np.abs(Cell_dia_xdata-Cell_dia[i])).argmin()]
        Yel_int_p=Yel_int_cdf[(np.abs(Yel_int_xdata-Yel_int[i])).argmin()]#
        #print("Yel_int_cdf",Yel_int_cdf)
        #print("Yel_int_xdata",Yel_int_xdata)
        All_p=Blue_dia_p*Cell_dia_p*Yel_int_p#Blue_dia_p*Cell_dia_p*Yel_int_p   (Blue_dia_p+Cell_dia_p+Yel_int_p)/3
        P_value_a_cell=[]
        P_value_a_cell.append("P_value")
        P_value_a_cell.append(Blue_dia_p)
        P_value_a_cell.append(Cell_dia_p)
        P_value_a_cell.append(Yel_int_p)
        P_value_a_cell.append(All_p)
        
        
        result_a_cell.append(info_a_cell)
        result_a_cell.append(BayF_a_cell)
        result_a_cell.append(P_value_a_cell)
        write_excel_xls_append(xls_path, result_a_cell)
    
        
#        if BayF_ALL>=10 and BayF_ALL<30:
#            label=info_nuccc[i]['label']
#            CTC.append(i)
#            temp=(nmask2==label)
#            CTC_mask_strong=CTC_mask_strong+temp    
#            
#        if BayF_ALL>=30 and BayF_ALL<100:
#            label=info_nuccc[i]['label']
#            CTC.append(i)
#            temp=(nmask2==label)
#            CTC_mask_very_strong=CTC_mask_very_strong+temp 
#        if BayF_ALL>=100:
#            label=info_nuccc[i]['label']
#            CTC.append(i)
#            temp=(nmask2==label)
#            CTC_mask_extreme=CTC_mask_extreme+temp   
        print("All_p:",All_p)
        if All_p>=0.05 and All_p<0.3:
            label=info_nuccc[i]['label']
            CTC.append(i)
            temp=(nmask2==label)
            CTC_mask_strong=CTC_mask_strong+temp    
            
        if All_p>=0.01 and All_p<0.05:
            label=info_nuccc[i]['label']
            CTC.append(i)
            temp=(nmask2==label)
            CTC_mask_very_strong=CTC_mask_very_strong+temp 
        if All_p<=0.01 :
            label=info_nuccc[i]['label']
            CTC.append(i)
            temp=(nmask2==label)
            CTC_mask_extreme=CTC_mask_extreme+temp                    
    #plt.imshow(CTC_mask_strong)
    contours1 , _ = cv2.findContours (CTC_mask_strong , cv2.RETR_EXTERNAL , cv2.CHAIN_APPROX_SIMPLE )
    contours2 , _ = cv2.findContours (CTC_mask_very_strong , cv2.RETR_EXTERNAL , cv2.CHAIN_APPROX_SIMPLE )
    contours3 , _ = cv2.findContours (CTC_mask_extreme , cv2.RETR_EXTERNAL , cv2.CHAIN_APPROX_SIMPLE )
    import copy
    seg_CTC=copy.deepcopy(Added_image)
    #seg_CTC=cv2.cvtColor(seg_CTC, cv2.COLOR_GRAY2BGR)
    cv2.drawContours(seg_CTC,contours1,-1,(255,174,201),8) 
    cv2.drawContours(seg_CTC,contours2,-1,(255,234,0),8)
    cv2.drawContours(seg_CTC,contours3,-1,(255,0,0),8)  
#    plt.figure()
#    plt.imshow(cv2.cvtColor(seg_CTC, cv2.COLOR_BGR2RGB))
    #plt.imshow(cmask2)
    CTC_num=len(CTC)  
    All_num=len(info_cell)
    
    plt.figure()
    plt.style.use('classic')
    plt.imshow(seg_CTC)
    #savefile='/home/Yuni/CTC/03_CTC_detect/02_CTC_visualize_images/'+n+'_'+num+'_result.png'
    #plt.savefig(savefile,dpi=600)
    plt.show(block=True)
    return(CTC_num,All_num)

# import pathlib

# xls_path="/home/Yuni/CTC/03_CTC_detect/01_cell_info_estimate/estimate_result.xls"
# file = pathlib.Path("xls_path")
# if file.exists ():
#     print ("File exist")
# else:
#     print ("Create ",xls_path)
#     sheet_name_xls="cell_all_info"
#     value = [["Cell", "Cell_diameter", "Nucleus_diameter", "Yellow_intensity", "All"],]
#     write_excel_xls(xls_path, sheet_name_xls, value)

 
# #processing-------------------       
# #os.makedirs(f'/home/Yuni/CTC/result/')
# input_path="/home/Yuni/CTC/01.2_50tumor_samples/02_segmentation_result/"
# os.chdir(input_path)
# resultlist = os.listdir(input_path)
# print("list:",resultlist)
# #leng=0
# for files in resultlist:
#     if fnmatch.fnmatch(files, '02*'):
#         print("processing",files)
#         #leng=leng+1
#         name_num=files.split(".")[1]
#         n='02'
#         #print(name_num,"of ",n)
#         f=open(f'{input_path}{n}.{name_num}.result.pkl','rb')
#         seg_all2,nmask2,cmask2,img_c2,Blue_dia,Cell_dia,Yel_int,Blue_int,All_num,info_yellow,info_yellow_backg,img_nucle_rescale,Added_image,cell_minus=pickle.load(f)
#         f.close()
#         #plt.imshow(seg_all2)
#         CTC_estimate(cmask2,nmask2,info_yellow,info_yellow_backg,seg_all2,img_nucle_rescale,Added_image,cell_minus,n,name_num,xls_path)
#         print(f"finished {name_num} image of {n}")          


                   