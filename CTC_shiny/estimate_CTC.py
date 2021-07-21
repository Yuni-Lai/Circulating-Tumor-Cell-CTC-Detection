# -*- coding: utf-8 -*-
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
import pathlib

os.chdir('C:/Users/10158/OneDrive - City University of Hong Kong/ProgramCode/CTC_Image/02_shiny/CTC/CTC_shiny/patientes_78/')
f = open('Blue_dia_all.pckl', 'rb')
Blue_dia_all = pickle.load(f)
f.close()
f = open('Cell_dia_all.pckl', 'rb')
Cell_dia_all = pickle.load(f)
f.close()
f = open('Blue_int_all.pckl', 'rb')
Blue_int_all = pickle.load(f)
f.close()

os.chdir('C:/Users/10158/OneDrive - City University of Hong Kong/ProgramCode/CTC_Image/02_shiny/CTC/CTC_shiny/34')
f = open('Yel_int_all.pckl', 'rb')
Yel_int_all = pickle.load(f)
f.close()

os.chdir('C:/Users/10158/OneDrive - City University of Hong Kong/ProgramCode/CTC_Image/02_shiny/CTC/CTC_shiny/')
f = open('GMMparameters.pckl', 'rb')
Cell_dia_Mu,Cell_dia_SigmaSquare,Cell_dia_Alpha,Yel_int_Mu,Yel_int_SigmaSquare,Yel_int_Alpha,Blue_int_Mu,Blue_int_SigmaSquare,Blue_int_Alpha,Blue_dia_Mu,Blue_dia_SigmaSquare,Blue_dia_Alpha = pickle.load(f)
f.close()

f = open('P_parameters.pckl', 'rb')
Blue_dia_cut,Blue_dia_xdata, Blue_dia_ydata,Blue_dia_cdf,Blue_int_cut,Blue_int_xdata, Blue_int_ydata,Blue_int_cdf,Yel_int_cut,Yel_int_xdata, Yel_int_ydata,Yel_int_cdf,Cell_dia_cut,Cell_dia_xdata, Cell_dia_ydata,Cell_dia_cdf = pickle.load(f)
f.close()


Blue_dia_all1= list(itertools.chain.from_iterable(Blue_dia_all))

Blue_int_all1= list(itertools.chain.from_iterable(Blue_int_all))

Yel_int_all1= list(itertools.chain.from_iterable(Yel_int_all))

Cell_dia_all1= list(itertools.chain.from_iterable(Cell_dia_all))


def gaussian(sigma, x, u):
	y = np.exp(-(x - u) ** 2 / (2 * sigma ** 2)) / (sigma * math.sqrt(2 * math.pi))
	return y


#--------write excel--------------------------------
import xlrd
import xlwt
#from xlutils.copy import copy
 
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
    from xlutils.copy import copy
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
def CTC_estimate(info_cell,info_nuccc,cmask1,nmask1,info_yellow,img_yellow_backg,img_c,Added_image):
    xls_path="C:/Users/10158/OneDrive - City University of Hong Kong/ProgramCode/CTC_Image/02_shiny/CTC/CTC_shiny/estimate_result.xls"
    file = pathlib.Path("xls_path")
    if file.exists ():
        print ("File exist")
        os.remove(xls_path)

    print ("Create ",xls_path)
    sheet_name_xls="cell_all_info"
    value = [["Cell", "Cell_diameter", "Nucleus_diameter", "Yellow_intensity", "All"],]
    write_excel_xls(xls_path, sheet_name_xls, value)
    #aaa=measure.label(nmask1)
    #info_nuccc=measure.regionprops(aaa,intensity_image=img_nucle_rescale)
    Blue_int=np.zeros(len(info_nuccc))
    #Cell_int=np.zeros(len(info_nuccc))
    Blue_dia=np.zeros(len(info_nuccc))
    Cell_dia=np.zeros(len(info_cell))
    Yel_int=np.zeros(len(info_yellow))
    #intensity_Y=np.zeros(len(info_yellow))
    for i in range(0,len(info_nuccc)-1):
        Blue_int[i]=info_nuccc[i]['mean_intensity']
        Yel_int[i]=info_yellow[i]['mean_intensity']
        #Cell_int[i]=info_cell[i]['mean_intensity']
        Blue_dia[i]=info_nuccc[i]['equivalent_diameter']*0.172
        Cell_dia[i]=info_cell[i]['equivalent_diameter']*0.172
    #plt.imshow()
    CTC=[]
    #plt.figure()
    #plt.imshow(seg_all1)
    CTC_mask_strong=np.zeros(nmask1.shape, dtype=np.uint8) 
    CTC_mask_very_strong=np.zeros(nmask1.shape, dtype=np.uint8) 
    CTC_mask_extreme=np.zeros(nmask1.shape, dtype=np.uint8) 
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
        if M2_Blue_dia==0:
            M2_Blue_dia=M2_Blue_dia+0.00001
        BayF_Blue_dia=math.exp(math.log(M2_Blue_dia)-math.log(M1_Blue_dia))
        #Cell_dia-----------
        M2_Cell_dia=gaussian(np.sqrt(Cell_dia_SigmaSquare[0][2]), Cell_dia[i], Cell_dia_Mu[0][2])*Cell_dia_Alpha[0][2]
        M1_Cell_dia=gaussian(np.sqrt(Cell_dia_SigmaSquare[0][0]), Cell_dia[i], Cell_dia_Mu[0][0])*Cell_dia_Alpha[0][0]+gaussian(np.sqrt(Cell_dia_SigmaSquare[0][1]), Cell_dia[i], Cell_dia_Mu[0][1])*Cell_dia_Alpha[0][1]
        if M1_Cell_dia==0:
            M1_Cell_dia=M1_Cell_dia+0.00001
        if M2_Cell_dia==0:
            M2_Cell_dia=M2_Cell_dia+0.00001
        BayF_Cell_dia=math.exp(math.log(M2_Cell_dia)-math.log(M1_Cell_dia))
        #Yel_int-------------
        M2_Yel_int=gaussian(np.sqrt(Yel_int_SigmaSquare[0][0]), Yel_int[i], Yel_int_Mu[0][0])*Yel_int_Alpha[0][0]
        M1_Yel_int=gaussian(np.sqrt(Yel_int_SigmaSquare[0][1]), Yel_int[i], Yel_int_Mu[0][1])*Yel_int_Alpha[0][1]
        print(i)
        if M1_Yel_int==0:
            M1_Yel_int=M1_Yel_int+0.00001
        if M2_Yel_int==0:
            M2_Yel_int=M2_Yel_int+0.00001
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
        All_p=max(Blue_dia_p,Cell_dia_p,Yel_int_p)#Blue_dia_p*Cell_dia_p*Yel_int_p   (Blue_dia_p+Cell_dia_p+Yel_int_p)/3
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
#            temp=(nmask1==label)
#            CTC_mask_strong=CTC_mask_strong+temp    
#            
#        if BayF_ALL>=30 and BayF_ALL<100:
#            label=info_nuccc[i]['label']
#            CTC.append(i)
#            temp=(nmask1==label)
#            CTC_mask_very_strong=CTC_mask_very_strong+temp 
#        if BayF_ALL>=100:
#            label=info_nuccc[i]['label']
#            CTC.append(i)
#            temp=(nmask1==label)
#            CTC_mask_extreme=CTC_mask_extreme+temp   
        print("All_p:",All_p)
        if All_p<=0.05 and All_p<0.1:
            label=info_nuccc[i]['label']
            CTC.append(i)
            temp=(nmask1==label)
            CTC_mask_strong=CTC_mask_strong+temp    
            
        if All_p>=0.01 and All_p<0.05:
            label=info_nuccc[i]['label']
            CTC.append(i)
            temp=(nmask1==label)
            CTC_mask_very_strong=CTC_mask_very_strong+temp 
        if All_p<=0.01 :
            label=info_nuccc[i]['label']
            CTC.append(i)
            temp=(nmask1==label)
            CTC_mask_extreme=CTC_mask_extreme+temp                    
    #    matplotlib.use('Qt5Agg',force=True)
    #plt.imshow(nmask1)
    #plt.imshow(CTC_mask_extreme)
    #plt.imshow(CTC_mask_strong)
    #plt.imshow(CTC_mask_very_strong)
    contours1 , _ = cv2.findContours (CTC_mask_strong , cv2.RETR_EXTERNAL , cv2.CHAIN_APPROX_SIMPLE )
    contours2 , _ = cv2.findContours (CTC_mask_very_strong , cv2.RETR_EXTERNAL , cv2.CHAIN_APPROX_SIMPLE )
    contours3 , _ = cv2.findContours (CTC_mask_extreme , cv2.RETR_EXTERNAL , cv2.CHAIN_APPROX_SIMPLE )
    import copy
    seg_CTC=copy.deepcopy(Added_image)
    #seg_CTC=cv2.cvtColor(seg_CTC, cv2.COLOR_GRAY2BGR)
    cv2.drawContours(seg_CTC,contours1,-1,(255,174,201),8) 
    cv2.drawContours(seg_CTC,contours2,-1,(255,234,0),8)
    cv2.drawContours(seg_CTC,contours3,-1,(255,0,0),8)
    #plt.imshow(seg_CTC)
#    plt.figure()
#    plt.imshow(cv2.cvtColor(seg_CTC, cv2.COLOR_BGR2RGB))
    #fig.savefig('test.jpg', format='jpg', transparent=True, dpi=1000, pad_inches = 0)
    #plt.imshow(cmask2)
    CTC_num=len(CTC)
    All_num=len(info_cell)
    return(CTC_num,All_num,seg_CTC)

