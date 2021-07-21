# -*- coding: utf-8 -*-
import matplotlib
#matplotlib.use('TkAgg')
import cv2
import numpy as np
import os
import pickle
from skimage import data,color
import matplotlib.pyplot as plt
from skimage.morphology import disk
from skimage import measure
import skimage.filters.rank as sfr
from PIL import ImageFilter
from PIL import Image
from PIL import ImageOps
from numpy import asarray
from skimage.filters.rank import median
from skimage.filters.rank import maximum, minimum, gradient
from skimage.filters.rank import entropy
from skimage.filters.rank import autolevel_percentile
from skimage.filters.rank import enhance_contrast
from skimage import morphology
#from skimage.measure import label
import copy
import random
from skimage import transform as tf
from skimage import io
from skimage import img_as_ubyte
from skimage import exposure
from scipy import ndimage as ndi
from skimage import morphology,feature
import pylab
import mahotas as mh
from skimage.segmentation import slic
from skimage.segmentation import mark_boundaries
import fnmatch
from PIL import Image
import matplotlib
from scipy.spatial import Voronoi, voronoi_plot_2d
from skimage import filters

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


rootPath='/dataserver145/image/Yuni/CTC/01.3_healthy_samples/'


os.chdir('/home/Yuni/CTC_shiny/patientes_78/')
f = open('Blue_dia_all.pckl', 'rb')
Blue_dia_all = pickle.load(f)
f.close()
f = open('Cell_dia_all.pckl', 'rb')
Cell_dia_all = pickle.load(f)
f.close()
f = open('Blue_int_all.pckl', 'rb')
Blue_int_all = pickle.load(f)
f.close()

os.chdir('/home/Yuni/CTC_shiny/34/')
f = open('Yel_int_all.pckl', 'rb')
Yel_int_all = pickle.load(f)
f.close()

os.chdir('/home/Yuni/CTC_shiny/')
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
def CTC_estimate(info_cell,info_nuccc,cmask1,nmask1,info_yellow,img_yellow_backg,img_c,Added_image,label_cell):
    xls_path=rootPath+"estimate_result.xls"
    file = pathlib.Path("xls_path")
    if file.exists():
        #print ("File exist")
        pass
        #os.remove(xls_path)
    if file.exists()==False:   
        print ("Create ",xls_path)
        sheet_name_xls="cell_all_info"
        value = [["Cell", "Cell_diameter", "Nucleus_diameter", "Yellow_intensity", "All"],]
        write_excel_xls(xls_path, sheet_name_xls, value)
    Blue_int=np.zeros(len(info_nuccc))
    #Cell_int=np.zeros(len(info_nuccc))
    Blue_dia=np.zeros(len(info_nuccc))
    Cell_dia=np.zeros(len(info_cell))
    Yel_int=np.zeros(len(info_yellow))
    #intensity_Y=np.zeros(len(info_yellow))
    for i in range(0,len(info_cell)):
        #i=10
        label=label_cell[i]
        index=label-1
        Blue_int[index]=info_nuccc[index]['mean_intensity']
        Yel_int[index]=info_yellow[index]['mean_intensity']
        #Cell_int[i]=info_cell[i]['mean_intensity']
        Blue_dia[index]=info_nuccc[index]['equivalent_diameter']*0.172
        Cell_dia[index]=info_cell[i]['equivalent_diameter']*0.172
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
        #print(i)
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
        #write_excel_xls_append(xls_path, result_a_cell)
    
        
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
        #print("All_p:",All_p)
        if All_p<=0.05 and All_p<0.1:
            label=info_nuccc[i]['label']
            #CTC.append(i)
            temp=(nmask1==label)
            CTC_mask_strong=CTC_mask_strong+temp    
            
        if All_p>=0.001 and All_p<0.05:
            label=info_nuccc[i]['label']
            #CTC.append(i)
            temp=(nmask1==label)
            CTC_mask_very_strong=CTC_mask_very_strong+temp 
        if All_p<=0.001 :
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
    print('----------1--------------------------------')
    #seg_CTC=cv2.cvtColor(seg_CTC, cv2.COLOR_GRAY2BGR)
    cv2.drawContours(seg_CTC,contours1,-1,(255,174,201),8) 
    cv2.drawContours(seg_CTC,contours2,-1,(255,234,0),8)
    cv2.drawContours(seg_CTC,contours3,-1,(255,0,0),8)
    print('------------2-----------------------')
    #plt.imshow(seg_CTC)
    plt.figure()
    plt.style.use('classic')
    plt.imshow(seg_CTC)
    savefile=rootPath+'03_1_CTC_result_images/'+name+"_seg_CTC.jpg"
    plt.savefig(savefile,dpi=600)
    #cv2.imwrite('/home/Yuni/CTC/01.3_healthy_patients/03_1_CTC_result_images/'+name+"_seg_CTC.jpg",seg_CTC)
    CTC_num=len(CTC)
    All_num=len(info_cell)
    #plt.close()
    return(Blue_dia,Cell_dia,Yel_int,Blue_int,CTC_num,All_num,seg_CTC)


matplotlib.use('Agg')
def preprocessing(num,name):  
    #num='1'
    #name="1_BF.tif"
    #======nucleus================
    nucle_n=num+"_DAPI.tif"
    img_nucle=cv2.imread(nucle_n)
    img_nucle_g = cv2.cvtColor(img_nucle, cv2.COLOR_BGR2GRAY)
    img_nucle_g =img_as_ubyte(img_nucle_g)
    p2, p98 = np.percentile(img_nucle_g, (1, 99))
    img_nucle_rescale = exposure.rescale_intensity(img_nucle_g, in_range=(p2, p98))
    del p2, p98
    edges = cv2.Canny(img_nucle_rescale,50,255)
    #plt.imshow(edges)
    #threshold------
    th2 = cv2.adaptiveThreshold(img_nucle_rescale,200,cv2.ADAPTIVE_THRESH_MEAN_C,cv2.THRESH_BINARY,79,-20)#79
    th2=th2-edges
    #plt.imshow(th2)
    th2=(th2>0)
    th2 = morphology.remove_small_holes(th2, 2500)
    th2=morphology.remove_small_objects(th2, min_size=300,connectivity=2)
    th2=morphology.opening(th2, disk(5))
    nuc_adp_th=np.uint8(th2)*255
        #nuc_adp_th=maximum(minimum(th2, disk(5)), disk(5))
    #plt.imshow(nuc_adp_th)
    del th2
    #watershade--------
    #nucle_median_g = cv2.medianBlur(img_nucle_rescale, 21)
    img_nucle_Gauss = cv2.GaussianBlur(img_nucle_rescale, ( 75,75), 0)
    rmax = mh.regmax(img_nucle_Gauss)
    #pylab.imshow(mh.overlay(img_nucle_g, rmax))
    seeds,nr_nuclei = mh.label(rmax)
    dist = mh.distance(nuc_adp_th)
    dist = 255 - mh.stretch(dist)  
    #pylab.imshow(dist)
    nuclei = mh.cwatershed(dist, seeds)
    nuclei *= nuc_adp_th
    del dist
    #plt.imshow(nuclei)
    info_nuccc=measure.regionprops(nuclei,intensity_image=img_nucle_g)

    #info_nuccc_list=list(range(0,len(info_nuccc)))
    nuc=np.zeros(img_nucle_g.shape, dtype=np.uint8) 
    global area
    print(len(info_nuccc))
    area=np.zeros(len(info_nuccc))
    for i in range(0,len(info_nuccc)): 
        #label=info_nuccc[i]['label']
        area[i]=info_nuccc[i]['area']
        min_row, min_col, max_row, max_col=info_nuccc[i]['bbox']
        temp=info_nuccc[i]['image']        
        #temp=morphology.closing(temp, disk(5))
        #temp=morphology.closing(temp, disk(5))
        temp=morphology.opening(temp, disk(10))
        temp=morphology.opening(temp, disk(10))
        temp=morphology.opening(temp, disk(10))
        temp=morphology.opening(temp, disk(10))
        temp=morphology.opening(temp, disk(10))
        #temp=morphology.closing(temp, disk(3))
        temp=morphology.binary_erosion(temp)
        nuc[min_row:max_row,min_col:max_col]=nuc[min_row:max_row,min_col:max_col]+temp
    area=np.median(area)

    nuc=nuc>0
    nuc_wat_th = morphology.remove_small_holes(nuc, area)
    nuc_wat_th=morphology.remove_small_objects(nuc_wat_th, min_size=area*0.3,connectivity=2)
    nuc_wat_th=np.uint8(nuc_wat_th)*255

    cv2.imwrite(name+"_nucle_mask.jpg",nuc_wat_th)
    
    #==========DIC==img====================================================== 
    name=num+'_BF.tif'  
    img = cv2.imread(name)
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    # Adaptive Equalization
    gray_adapteq = exposure.equalize_adapthist(gray,kernel_size=(70,70), clip_limit=0.05)  # Ã¨â¡ÂªÃ©â¬âÃ¥ÂºâÃ¥Ââ¡Ã¨Â¡Â?#    plt.imshow(gray_adapteq)
    gray_adapteq=img_as_ubyte(gray_adapteq)
    median_g = cv2.medianBlur(gray_adapteq, 21)
    enh = enhance_contrast(median_g, disk(10))
    backg =sfr.mean(median_g, disk(80)) 
    return(img,enh,backg,median_g,img_nucle_g,gray_adapteq,nuc_wat_th,img_nucle_rescale,img_nucle)


def segmentation_gradient(num,img,enh,img_nucle_rescale):
    #======nuc======
    grad_nuc = gradient(img_nucle_rescale, disk(8))
    grad_nuc = cv2.GaussianBlur(grad_nuc, ( 95,95), 0)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(5,5)) #MORPH_CROSS#MORPH_ELLIPSE
    #plt.imshow(grad_nuc)
    edges = cv2.Canny(img_nucle_rescale,10,200)
    edges = cv2.morphologyEx(edges, cv2.MORPH_CLOSE, kernel,iterations=5)
    #plt.imshow(edges)
   # median(img_nucle_rescale, morphology.disk(5))
    grad_nuc=grad_nuc>35
    summ=(grad_nuc==True).sum()
    if summ>img.shape[0]*img.shape[1]*0.7:#70%
        grad_nuc=grad_nuc>45
    if summ>img.shape[0]*img.shape[1]*0.85:#85%
        grad_nuc=grad_nuc>55
    grad_nuc = np.uint8(grad_nuc)*255
    #plt.imshow(grad_nuc)
    
    grad_nuc = cv2.morphologyEx(grad_nuc, cv2.MORPH_CLOSE, kernel,iterations=3)
    nuc_th_c = grad_nuc > 0
    nuc_th_c = morphology.remove_small_holes(nuc_th_c, 3500)
    nuc_th_c = morphology.remove_small_objects(nuc_th_c, min_size=area*0.4,connectivity=2)
    #plt.imshow(nuc_th_c)
    nuc_th_c=np.uint8(nuc_th_c)*255
    
    #====cell image===========
    grad = gradient(enh, disk(8))

    grad_Gauss = cv2.GaussianBlur(grad, ( 75,75), 0)

    grad_mask=grad_Gauss>25
    grad_mask = np.uint8(grad_mask)*255
    #plt.imshow(grad_mask)
    closed_g = cv2.morphologyEx(grad_mask, cv2.MORPH_CLOSE, kernel,iterations=5)

    arr = closed_g > 0
    arr1 = morphology.remove_small_holes(arr, 2000)
    arr2 = morphology.remove_small_objects(arr1, min_size=area*0.4,connectivity=2)
    #plt.imshow(arr2)
    img_th=np.uint8(arr2)*255
    return(img_th,nuc_th_c)
    


def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)

  
def re_segmentation(img_c,cellmask,nuc_wat_th,img_nucle_rescale):   
    nmask=(nuc_wat_th>0)
    nmask=morphology.remove_small_objects(nmask, min_size=500,connectivity=2)
    #plt.imshow(nmask)
    #plt.savefig("/home/Yuni/CTC/01.3_healthy_patients/test.png")
    nmask=measure.label(nmask)
    info_nuccc=measure.regionprops(nmask,intensity_image=img_nucle_rescale)
    center=[]
    for i in range(0,len(info_nuccc)): 
        center.append(info_nuccc[i]['centroid'])
    points = np.array(center)
    vor = Voronoi(points)
#    regions, vertices = voronoi_finite_polygons_2d(vor)
#    plt.figure()
#    cellmask222=cv2.flip(np.rot90(cellmask),0)
#    plt.imshow(cv2.cvtColor(cellmask222, cv2.COLOR_BGR2RGB))
#    #plt.plot(points[:, 0], points[:, 1], 'o')
#    for region in regions:
#        polygon = vertices[region]
#        plt.fill(*zip(*polygon), alpha=0.4 ,color = "k",edgecolor="k",linewidth=3)
    #plt.plot(points[:,0], points[:,1], 'ko')
#    plt.xlim(vor.min_bound[0], vor.max_bound[0])
#    plt.ylim(vor.min_bound[1], vor.max_bound[1])
    #plt.show()
#    img_1 = np.zeros(cellmask.shape, dtype=np.uint8) 
    #cv2.flip(image, 0)
    #-------------------
    #matplotlib.use('Agg',force=True)

    plt.figure()
    plt.imshow(cv2.cvtColor(cellmask, cv2.COLOR_BGR2RGB))
    #plt.plot(points[:, 1], points[:, 0], 'o')
    for simplex in vor.ridge_vertices:
         simplex = np.asarray(simplex)
         if np.all(simplex >= 0):
             plt.plot(vor.vertices[simplex, 1], vor.vertices[simplex, 0], 'k-',linewidth=0.5)
    
    #-------------------------
    plt.xlim(0,cellmask.shape[1]); plt.ylim(cellmask.shape[0],0)
    plt.axis('off')
    #cv2.flip(np.rot90(np.rot90(np.rot90(seg_all))),1)  
    fig = plt.gcf()
    fig.set_size_inches(cellmask.shape[1]/1000,cellmask.shape[0]/1000)
    
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
    plt.margins(0,0)
    fig.savefig('cmask.jpg', format='jpg', transparent=True, dpi=1000, pad_inches = 0)
    #plt.close()
    cmask=cv2.imread('cmask.jpg')
    #cmask.shape
    cmask=cmask[:,:,0]
    cmask[cmask<=220]=0
    cmask[cmask>220]=1
    #nmask[nmask>0]=1
    #cmask=cmask+nmask
    #combine=np.multiply(cmask,nmask)
    cmask=(cmask>=1)
    cmask=morphology.remove_small_objects(cmask, min_size=500,connectivity=2)
    cmask=morphology.remove_small_holes(cmask, 500)
    #cmask[cmask>0]=255
    #matplotlib.use('Qt5Agg',force=True)
    cmask=measure.label(cmask)
    
  #np.uint8(arr2)*255
#    matplotlib.use('Qt5Agg',force=True)
#    plt.imshow(cmask)
    #plt.imshow(nmask)
    from scipy import stats
    import copy
    cmask_c = copy.deepcopy(cmask)
    info_cell=measure.regionprops(cmask,intensity_image=nmask)

    for i in range(0,len(info_cell)):
        #i=len(info_cell)-1
        label=info_cell[i]['label']
        if info_cell[i]['mean_intensity']==0:
            #print(info_cell[i]['mean_intensity'])
            cmask[cmask_c == label] = 0
        else:
            #print(info_cell[i]['max_intensity'])
            temp=info_cell[i]['intensity_image']
            temp=temp[temp>0]
            temp=stats.mode(temp)[0][0]
            cmask[cmask_c == label]=temp
    cmask[(cmask==0) & (nmask>0)]=cmask[(cmask==0) & (nmask>0)]+nmask[(cmask==0) & (nmask>0)]
    info_cell=measure.regionprops(cmask,intensity_image=nmask)
    
    # plt.figure()
    # Temp_c=cmask[min_row-100:max_row+100,min_col-100:max_col+100]
    # plt.imshow(Temp_c)
    # Temp_n=nmask[min_row-100:max_row+100,min_col-100:max_col+100]
    # #1850:1980,2100:2250
    # plt.figure()
    # plt.imshow(Temp_n)
    # Temp=np.zeros(cmask.shape, dtype=np.uint8)
    # Temp[(cmask==0)&(nmask>0)]=1
    # plt.imshow(Temp)
    label_nuc=[]
    for i in range(0,len(info_nuccc)): 
        label_nuc.append(info_nuccc[i]['label'])
    label_cell=[]
    for i in range(0,len(info_cell)): 
        label_cell.append(info_cell[i]['label'])
    for number in label_nuc:
        if number not in label_cell:
            print(number)#1657
            #number=1657
            info_nuc_temp=measure.regionprops(nmask,intensity_image=cmask)
            if info_nuc_temp[number-1]['label']==number:
                temp=info_nuc_temp[number-1]['intensity_image']
                temp=temp[temp>0]
                temp=stats.mode(temp)[0][0]
                cmask[cmask == temp]=0
    cmask[(cmask==0) & (nmask>0)]=cmask[(cmask==0) & (nmask>0)]+nmask[(cmask==0) & (nmask>0)]
    info_cell=measure.regionprops(cmask,intensity_image=nmask)
    label_cell=[]
    for i in range(0,len(info_cell)): 
        label_cell.append(info_cell[i]['label'])
    #info_nuccc[1656]['label']
    #info_nuccc[1656]['bbox']
    #min_row, min_col, max_row, max_col=info_nuccc[1656]['bbox']
    #temp=nmask[min_row-50:max_row+50,min_col-50:max_col+50]
  # plt.imshow(temp)
    # plt.imshow(color.label2rgb(cmask))
    # plt.imshow(nmask)
    
    # nmask_temp=(nmask==1)
    # cmask_temp=(cmask==1)
    # plt.figure()
    # plt.imshow(cmask_temp)
    # plt.figure()
    # plt.imshow(nmask_temp)
 
    ccc=color.label2rgb(cmask,img_c,bg_label=0,colors=[(0,50,0)],alpha=0.008)
#    plt.figure()
#    plt.imshow(ccc)
    cccc=color.label2rgb(nuc_wat_th,bg_label=0,colors=[(0,0,10)])
    #plt.imshow(cccc)
    seg_all=cv2.addWeighted(ccc, 0.7, cccc, 0.3, 0)
    #plt.imshow(seg_all)
    return (cmask,nmask,seg_all,info_cell,info_nuccc,label_cell)

'''
    info_nuc_cell=measure.regionprops(nmask,intensity_image=cmask)
    for i in range(0,len(info_nuc_cell)):
       if info_nuc_cell[i]['mean_intensity']==0:
           label=info_nuc_cell[i]['label']
           nmask[nmask == label] = 0
    for i in range(0,len(info_nuccc)): 
        #label=info_nuccc[i]['label']
        area[i]=info_nuccc[i]['area']
        min_row, min_col, max_row, max_col=info_nuccc[i]['bbox']
        temp=info_nuccc[i]['image']        
        nuc[min_row:max_row,min_col:max_col]=nuc[min_row:max_row,min_col:max_col]+temp
  
'''
def segmentation_diff(num,img_c,enh_c,backg_c,nuc_wat_th,name):               
    cell_minus=cv2.absdiff(enh_c,backg_c) 
    _,minus_th = cv2.threshold(cell_minus,45,255,cv2.THRESH_BINARY)
    overlap1=(minus_th*nuc_wat_th>0)
    overlap_ratio=((overlap1==True).sum())/(((nuc_wat_th>0)==True).sum())
#    plt.imshow(nuc_wat_th)
#    plt.imshow(minus_th)
    if overlap_ratio<0.25:
        _,minus_th = cv2.threshold(cell_minus,20,255,cv2.THRESH_BINARY)
    if overlap_ratio>0.6:
        _,minus_th = cv2.threshold(cell_minus,65,255,cv2.THRESH_BINARY)
    minus_th=(minus_th>0)
    minus_th = morphology.remove_small_holes(minus_th, area)
    cell=np.uint8(minus_th)*255
    #minus_th=maximum(minimum(minus_th, disk(5)), disk(5))

#    plt.imshow(cell)
    comb=cell+nuc_wat_th
#    plt.imshow(comb)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(5,5)) #MORPH_CROSS#MORPH_ELLIPSE
    closed = cv2.morphologyEx(comb, cv2.MORPH_CLOSE, kernel,iterations=3)
#    plt.imshow(closed) 
    opened1 = cv2.morphologyEx(closed, cv2.MORPH_OPEN, kernel,iterations=3) 
#    plt.imshow(opened1)  
    arr = opened1 > 0
    arr1 = morphology.remove_small_holes(arr, area)
    arr2 = morphology.remove_small_objects(arr1, min_size=area*0.5,connectivity=2)
    #plt.imshow(arr2)
    cellmask=np.uint8(arr2)*255
    cv2.imwrite(name+"_cell_mask_c2.jpg",cellmask)
    return (cellmask,cell_minus)


def draw_approx_hull_polygon(img, cnts):
    img1 = np.zeros(img.shape, dtype=np.uint8) 
    contour_smooth=[]
    for cnt in cnts: 
        cnt=np.squeeze(cnt, axis=None)
        row_array=range(0,cnt.shape[0],30)
        cnt = cnt[row_array]
        contour_smooth.append(np.expand_dims(np.array(measure.subdivide_polygon(cnt , degree=2)),1))
    contour_smooth=[np.around(c) for c in contour_smooth]
    contour_smooth=[np.int32(c) for c in contour_smooth]
    cv2.drawContours(img1,contour_smooth,-1, (255, 255, 255), thickness=-1)
    img1=cv2.cvtColor(img1, cv2.COLOR_BGR2GRAY)
    return img1

#def compute_centroids(X):
#    n,_ = X.shape
#    centroids=[ (sum(X[:,0]) / n),(sum(X[:,1]) / n)]
#    return(centroids)
def check( cmask, nmask, rowMove,colMove):
    tform = tf.SimilarityTransform(translation=(-(rowMove),(colMove)))
    import copy
    cmask_c = copy.deepcopy(cmask)
    #cmask_c.transform=np.array([[1,0,rowMove],[0,1,colMove],[0,0,1]])
    cmask_c = tf.warp(cmask,tform)
    overlap_t=(cmask_c*nmask>0)
    overlap_count=(overlap_t==True).sum()
    return (overlap_count)
def adjust_img_ovelap(img,enh,img_th,nuc_th_c,name,img_nucle_g,backg,nuc_wat_th):
    overlap=(img_th*nuc_th_c>0)
    overlap_count=(overlap==True).sum()
    maxRow = 0
    maxCol=0  
    for rowMove in range(2,25):
        overlap_count_t=check(img_th, nuc_th_c, rowMove, maxCol)
        if overlap_count_t>overlap_count:
            maxRow=rowMove
            overlap_count=overlap_count_t
    tform = tf.SimilarityTransform(translation=(-(maxRow+3),(maxCol)))
    for colMove in range(1,25):
        overlap_count_t=check(img_th, nuc_th_c, maxRow, colMove)
        if overlap_count_t>overlap_count:
            maxCol=colMove
            overlap_count=overlap_count_t
        #maxRow = max(overlap_count, overlap_count_t)#????????????
    tform = tf.SimilarityTransform(translation=(-(maxRow+3),(maxCol+3)))
    img_th_c = copy.deepcopy(img_th)
    #cmask_c.transform=np.array([[1,0,rowMove],[0,1,colMove],[0,0,1]])
    img_th_c = tf.warp(img_th,tform)
    img_th_c=img_th_c*255
    img_c= img_as_ubyte(tf.warp(img,tform))
#    plt.figure()
#    plt.imshow(img_c)
    io.imsave(name+'_c.jpg', img_c)
    contours4 , _ = cv2.findContours (nuc_wat_th , cv2.RETR_EXTERNAL , cv2.CHAIN_APPROX_SIMPLE )
    img_cont3=copy.deepcopy(img_c)
    cv2.drawContours(img_cont3,contours4,-1,(0,0,255),3) 
#    plt.imshow(img_cont3)   
    enh_c = img_as_ubyte(tf.warp(enh,tform))
    backg_c=img_as_ubyte(tf.warp(backg,tform))
    return(img_c,img_th_c,enh_c,backg_c,tform)


def compute_centroids(X):
    n,_ = X.shape
    centroids=[ (sum(X[:,0]) / n),(sum(X[:,1]) / n)]
    return(centroids)  


def adjust_img_center(img,enh,cmask,info_nuc,info_cell,name,backg):
    #info_cell=measure.regionprops(cmask)
    #info_nuc=measure.regionprops(nmask)
    list_nuc=list(range(0,len(info_nuc)))
    random.seed(10)
    
    sample_N=len(info_nuc)
    random_num = random.sample(list_nuc, sample_N) 
    
    nuc_X=np.zeros((sample_N,2))
    cell_X=np.zeros((sample_N,2))
    for i,label in enumerate(random_num):
        nuc_X[i,]=info_nuc[label]["centroid"]
        cell_X[i,]=info_cell[label]["centroid"]
    nuc_X[1:30,]
    cell_X[1:30,]
    cent_nuc=compute_centroids(nuc_X)
    cent_cell=compute_centroids(cell_X)
    
    move_x=cent_cell[0]-cent_nuc[0]
    move_x
    move_y=cent_cell[1]-cent_nuc[1]
    move_y
    tform = tf.SimilarityTransform(translation=(-(move_x),-(move_y)))
    img_c= img_as_ubyte(tf.warp(img,tform))
    io.imsave(name+'_c2.jpg', img_c)
    cmask_c = img_as_ubyte(tf.warp(cmask,tform))
    enh_c = img_as_ubyte(tf.warp(enh,tform))
    backg_c=img_as_ubyte(tf.warp(backg,tform))
    return(img_c,cmask_c,enh_c,backg_c,tform)

def yellow(num,cmask1,nmask1,img_c,img_nucle_rescale):
    yel_n=num+"_TRITC.tif"
    img_yellow=cv2.imread(yel_n)
    #img_yellow_g = cv2.cvtColor(img_yellow, cv2.COLOR_BGR2GRAY)
    img_yellow_e = exposure.equalize_adapthist(img_yellow,kernel_size=(10,10), clip_limit=0.01)
    img_yellow_e=img_as_ubyte(img_yellow_e)
    img_yellow_g=cv2.cvtColor(img_yellow_e, cv2.COLOR_BGR2GRAY)
    info_yellow=measure.regionprops(nmask1,intensity_image=img_yellow_g/255)
    cmask22=morphology.dilation(cmask1,morphology.square(5)) 
    cmask22=morphology.dilation(cmask22,morphology.square(5)) 
    neg_mask=(cmask22==0)
    neg_mask=np.uint8(neg_mask)*1
    #plt.imshow(neg_mask)
    img_yellow_backg =(neg_mask*img_yellow_g).sum()/(neg_mask==1).sum()
    #gray_adapteq_C=cv2.cvtColor(img_c2, cv2.COLOR_GRAY2BGR)
    #gray_adapteq_C = exposure.equalize_adapthist(img_c,kernel_size=(5,5), clip_limit=0.05)  # Ã¨â¡ÂªÃ©â¬âÃ¥ÂºâÃ¥Ââ¡Ã¨Â¡Â?#    plt.imshow(gray_adapteq)
    gray_adapteq_C=img_c
    a,b=img_nucle_rescale.shape
    img_nucle_B = np.zeros((a, b,3)).astype("uint8")
    #img_nucle_rescale=cv2.cvtColor(img_nucle_rescale, cv2.COLOR_BGR2RGB)
    img_nucle_B[:,:,2] = img_nucle_rescale
    img_yellow_e=cv2.cvtColor(img_yellow_e, cv2.COLOR_BGR2RGB)
    #gray_adapteq_C=img_as_ubyte(tf.warp(img_as_ubyte(tf.warp(gray_adapteq_C,tform)),tform2))
    #plt.imshow(img_yellow_e)
    #plt.imshow(img_yellow_new)
    #plt.imshow(img_nucle_B)
    #plt.imshow(Added_image)
    A_image=cv2.addWeighted(gray_adapteq_C, 0.8, img_nucle_B, 0.2, 0)
    Added_image=cv2.addWeighted(A_image, 0.7, img_yellow_e, 0.3, 0)
    return(info_yellow,Added_image,img_yellow_backg)


def cell_detection(inputt):
    global name
    
    num,imageP=inputt[0],inputt[1]
    #num,imageP='25','/home/Yuni/CTC/01.3_healthy_patients/01_1_data_healthy/01_0902_healthy_data/474748_1'
    n=imageP.split('/')[-1]
    os.chdir(imageP)
    name=num+'_BF.tif'
    img,enh,backg,median_g,img_nucle_g,gray_adapteq,nuc_wat_th,img_nucle_rescale,img_nucle=preprocessing(num,name)
    img_th,nuc_th_c=segmentation_gradient(num,img,enh,img_nucle_rescale)
    img_c,img_th_c,enh_c,backg_c,tform=adjust_img_ovelap(img,enh,img_th,nuc_th_c,name,img_nucle_g,backg,nuc_wat_th)
    #plt.imshow(img)
    cellmask,cell_minus=segmentation_diff(num,img_c,enh_c,backg_c,nuc_wat_th,name)
    cmask1,nmask1,seg_all1,info_cell,info_nuccc,label_cell=re_segmentation(img_c,cellmask,nuc_wat_th,img_nucle_rescale)
    #img_c2,cmask_c2,enh_c2,backg_c2,tform2=adjust_img_center(img_c,enh_c,cmask1,info_nuccc,info_cell,name,backg_c)
    #cellmask2,cell_minus=segmentation_diff(num,img_c2,enh_c2,backg_c2,nuc_wat_th)
    #cmask2,nmask2,seg_all2,info_cell,info_nuccc=re_segmentation(img_c,cellmask2,nuc_wat_th,img_nucle_rescale)
    #plt.imshow(seg_all1)
    info_yellow,Added_image,img_yellow_backg=yellow(num,cmask1,nmask1,img_c,img_nucle_rescale)
    #seg_CTC,CTC_mask,Blue_dia,Cell_dia,Yel_int,Blue_int,CTC_num,All_num=CTC(cmask2,nuc_wat_th,info_yellow,info_yellow_backg,seg_all2,img_nucle_rescale,Added_image,cell_minus)
    Blue_dia,Cell_dia,Yel_int,Blue_int,CTC_num,All_num,seg_CTC = CTC_estimate(info_cell,info_nuccc,cmask1,nmask1,info_yellow,img_yellow_backg,img_c,Added_image,label_cell)
    f = open(f"{rootPath}02_1_segmented/{n}.{num}.result.pkl", 'wb')
    pickle.dump([Blue_dia,Cell_dia,Yel_int,Blue_int,seg_all1,All_num,CTC_num,seg_CTC,nmask1,cmask1,img_c,Added_image],f)
    f.close()
    plt.figure()
    plt.style.use('classic')
    plt.imshow(seg_all1)
    savefile=rootPath+'02_1_segmented/'+name+"_seg_image.jpg"
    plt.savefig(savefile,dpi=600)
    #plt.close()
    return(Blue_dia,Cell_dia,Yel_int,Blue_int)


#=====================================
def equalize_preprocessing(num,name):  
    nucle_n=num+"_DAPI.tif"
    img_nucle=cv2.imread(nucle_n)
    imgnucle_adapteq = exposure.equalize_adapthist(img_nucle, clip_limit=0.03) 
    name=num+'_BF.tif'  
    img = cv2.imread(name)
    img_adapteq = exposure.equalize_adapthist(img, clip_limit=0.03)  # èªéåºåè¡¡
    yel_n=num+"_TRITC.tif"
    img_yellow=cv2.imread(yel_n)
    yellow_adapteq = exposure.equalize_adapthist(img_yellow, clip_limit=0.03) 
    img_adapteq<-img_as_ubyte(img_adapteq)
    imgnucle_adapteq<-img_as_ubyte(imgnucle_adapteq)
    yellow_adapteq<-img_as_ubyte(yellow_adapteq)
    # cv2.imwrite(name+"_img_adapteq.jpg",img_as_ubyte(img_adapteq))
    # cv2.imwrite(name+"_nucle_adapteq.jpg",imgnucle_adapteq)
    # cv2.imwrite(name+"_yellow_adapteq.jpg",yellow_adapteq)
    return(img_adapteq,imgnucle_adapteq,yellow_adapteq)


def original(num,name):  
    nucle_n=num+"_DAPI.tif"
    img_nucle=cv2.cvtColor(cv2.imread(nucle_n),cv2.COLOR_BGR2RGB)
    name=num+'_BF.tif'  
    img = cv2.imread(name)
    yel_n=num+"_TRITC.tif"
    img_yellow=cv2.cvtColor(cv2.imread(yel_n),cv2.COLOR_BGR2RGB)
    return(img,img_nucle,img_yellow)
#os.chdir('C:/Users/10158/OneDrive - City University of Hong Kong/ProgramCode/CTC_Image/02_shiny/CTC/CTC_shiny')
#seg_all2,All_num,CTC_num,seg_CTC,Added_image=cell_detection('1','C:/Users/10158/OneDrive - City University of Hong Kong/ProgramCode/CTC_Image/02_shiny/CTC/CTC_shiny')
#matplotlib.use('Qt5Agg',force=True)
#plt.imshow(seg_all2)
#plt.imshow(seg_CTC)
#plt.imshow(Added_image)

Path="/dataserver145/image/Yuni/CTC/01.3_healthy_patients/01_1_data_healthy/01_0902_healthy_data/"
subfolder=os.listdir(Path)
imagePAll=[]
for sub in subfolder:
    imagePAll.append(Path+sub)
    
Blue_dia_all=[]
Cell_dia_all=[]
Yel_int_all=[]
Blue_int_all=[]

list_process=[]
for imageP in imagePAll:
    #imageP=imagePAll[0]
    os.chdir(imageP)
    n=imageP.split('/')[-1]
    imagelist = os.listdir(imageP)
    #leng=0
    for files in imagelist:
        if fnmatch.fnmatch(files, '*_BF.tif'):
            #leng=leng+1
            name_num=files.split("_")[0]
            list_process.append([name_num,imageP])
            print(name_num,"of ",imageP)


#cell_detection(list_process[0])

from tqdm import tqdm
import concurrent.futures
import time

def run(function, inputs):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(function, inputs), total=len(inputs)))
    return results


'''
start = time.perf_counter()
results = run(cell_detection, list_process)
for result in results:
    print(result)
    #Blue_dia,Cell_dia,Yel_int,Blue_int
    Blue_dia_all.append(result[0])
    Cell_dia_all.append(result[1])
    Yel_int_all.append(result[2])
    Blue_int_all.append(result[3])
    finish = time.perf_counter()
    print(f'Finished in {round(finish-start, 2)} second(s)')  
'''


def main():
    start = time.perf_counter()
    p = concurrent.futures.ProcessPoolExecutor(max_workers=20)
    results = p.map(cell_detection, list_process)
    p.shutdown(wait=True)
    for result in results:
        print(result[0])
        # #Blue_dia,Cell_dia,Yel_int,Blue_int
        # Blue_dia_all.append(result[0])
        # Cell_dia_all.append(result[1])
        # Yel_int_all.append(result[2])
        # Blue_int_all.append(result[3])
    finish = time.perf_counter()
    print(f'Finished in {round(finish-start, 2)} second(s)')
    # PATH='/home/Yuni/CTC/variables/healthy_01/'
    # if os.path.exists(PATH):
    #    print('file exit: ',PATH)
    # if not os.path.exists(PATH):
    #        os.makedirs(PATH)
    # os.chdir(PATH)
    # f = open('Blue_dia_all.pckl', 'wb')
    # pickle.dump(Blue_dia_all, f)
    # f.close()
    # f = open('Cell_dia_all.pckl', 'wb')
    # pickle.dump(Cell_dia_all, f)
    # f.close()
    # f = open('Yel_int_all.pckl', 'wb')
    # pickle.dump(Yel_int_all, f)
    # f.close()
    # f = open('Blue_int_all.pckl', 'wb')
    # pickle.dump(Blue_int_all, f)
    # f.close()
    
#args = list_process
if __name__ == "__main__":
    main()
    #results = executor.map(do_something, secs)


  
'''    
    # start = time.perf_counter()
    # p = ProcessPoolExecutor(max_workers=20)
    # results = p.map(cell_detection, list_process[21:40])
    # p.shutdown(wait=True)
    # for result in results:
    #     print(result)
    #     Blue_dia_all.append(result[0])
    #     Cell_dia_all.append(result[1])
    #     Yel_int_all.append(result[2])
    #     Blue_int_all.append(result[3])
    # finish = time.perf_counter()
    # print(f'Finished in {round(finish-start, 2)} second(s)')
    
    # start = time.perf_counter()
    # p = ProcessPoolExecutor(max_workers=20)
    # results = p.map(cell_detection, list_process[41:60])
    # p.shutdown(wait=True)
    # for result in results:
    #     print(result)
    #     Blue_dia_all.append(result[0])
    #     Cell_dia_all.append(result[1])
    #     Yel_int_all.append(result[2])
    #     Blue_int_all.append(result[3])
    # finish = time.perf_counter()
    # print(f'Finished in {round(finish-start, 2)} second(s)')
    
    
    # start = time.perf_counter()
    # p = ProcessPoolExecutor(max_workers=20)
    # results = p.map(cell_detection, list_process[61:])
    # p.shutdown(wait=True)
    # for result in results:
    #     print(result)
    #     Blue_dia_all.append(result[0])
    #     Cell_dia_all.append(result[1])
    #     Yel_int_all.append(result[2])
    #     Blue_int_all.append(result[3])
    # finish = time.perf_counter()
    # print(f'Finished in {round(finish-start, 2)} second(s)')
    
    os.makedirs(f'/home/Yuni/CTC/variables/{n}')
    os.chdir(f'/home/Yuni/CTC/variables/{n}')
    f = open('Blue_dia_all.pckl', 'wb')
    pickle.dump(Blue_dia_all, f)
    f.close()
    f = open('Cell_dia_all.pckl', 'wb')
    pickle.dump(Cell_dia_all, f)
    f.close()
    f = open('Yel_int_all.pckl', 'wb')
    pickle.dump(Yel_int_all, f)
    f.close()
    f = open('Blue_int_all.pckl', 'wb')
    pickle.dump(Blue_int_all, f)
    f.close()
    '''
