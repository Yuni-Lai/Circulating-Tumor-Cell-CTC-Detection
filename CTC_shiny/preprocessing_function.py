# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 12:06:38 2020
@author: Yuni
"""
import cv2
import numpy as np
import os
import PIL.Image as Image
from skimage import color,data
import matplotlib.pyplot as plt
from skimage.morphology import disk
from skimage import measure
import skimage.filters.rank as sfr
# from PIL import ImageFilter
# from PIL import Image
# from PIL import ImageOps
from numpy import asarray
from skimage.filters.rank import median
from skimage.filters.rank import maximum, minimum, gradient
from skimage.filters.rank import entropy
from skimage.filters.rank import autolevel_percentile
from skimage.filters.rank import enhance_contrast
from skimage import morphology
import copy
import random
#from rpy2.robjects import pandas2ri
#import rpy2.robjects as robjects
#from rpy2.robjects.packages import importr
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
#from PIL import Image
#import seaborn as sns
#for num in (1:4):
#def whole_hist(image):
#	'''
#	绘制整幅图像的直方图
#	'''
#	plt.hist(image.ravel(), 256, [0, 256]) #numpy的ravel函数功能是将多维数组降为一维数组
#	plt.show()
#
#def channel_hist(image):
#	'''
#	画三通道图像的直方图
#	'''
#	color = ('b', 'g', 'r')   #这里画笔颜色的值可以为大写或小写或只写首字母或大小写混合
#	for i , color in enumerate(color):
#		hist = cv2.calcHist([image], [i], None, [256], [0, 256])  #计算直方图
#		plt.plot(hist, color)
#		plt.xlim([0, 256])
#	plt.show()
#
#image = cv2.imread('/home/liubo/tmp/sample/00001.jpg')
#cv2.imshow('image', image)
#cv2.waitKey(0)


def equalize_preprocessing(num,name):  
    nucle_n=num+"_DAPI.tif"
    img_nucle=cv2.imread(nucle_n)
    imgnucle_adapteq = exposure.equalize_adapthist(img_nucle, clip_limit=0.03) 
    name=num+'_BF.tif'  
    img = cv2.imread(name)
    img_adapteq = exposure.equalize_adapthist(img, clip_limit=0.03)  # 自适应均衡
    yel_n=num+"_TRITC.tif"
    img_yellow=cv2.imread(yel_n)
    yellow_adapteq = exposure.equalize_adapthist(img_yellow, clip_limit=0.03) 
    img_adapteq<-img_as_ubyte(img_adapteq)
    imgnucle_adapteq<-img_as_ubyte(imgnucle_adapteq)
    yellow_adapteq<-img_as_ubyte(yellow_adapteq)
    cv2.imwrite(name+"_img_adapteq.jpg",img_as_ubyte(img_adapteq))
    cv2.imwrite(name+"_nucle_adapteq.jpg",imgnucle_adapteq)
    cv2.imwrite(name+"_yellow_adapteq.jpg",yellow_adapteq)
    return(img_adapteq,imgnucle_adapteq,yellow_adapteq)
    
def preprocessing(num,name):  
    #======nucleus================
    nucle_n=num+"_DAPI.tif"
    img_nucle=cv2.imread(nucle_n)
    #whole_hist(img_nucle)
    #channel_hist(img_nucle)

    #plt.figure()
    #plt.title("original DAPI image")
    #plt.imshow(img_nucle)
    #plt.imshow(cv2.cvtColor(img_nucle, cv2.COLOR_BGR2RGB))
    img_nucle_g = cv2.cvtColor(img_nucle, cv2.COLOR_BGR2GRAY)
    #img_nucle_g = exposure.equalize_adapthist(img_nucle_g, clip_limit=0.03)
    img_nucle_g =img_as_ubyte(img_nucle_g)
    #plt.imshow(cv2.cvtColor(img_nucle_g, cv2.COLOR_BGR2RGB))
    # Contrast stretching    #对比度拉伸-----
    p2, p98 = np.percentile(img_nucle_g, (1, 99))
    img_nucle_rescale = exposure.rescale_intensity(img_nucle_g, in_range=(p2, p98))
    #plt.imshow(cv2.cvtColor(img_nucle_rescale, cv2.COLOR_BGR2RGB))
    del p2, p98
    #edges = cv2.Canny(img_nucle_rescale,50,255)
    #plt.imshow(edges)
    #threshold------
    th2 = cv2.adaptiveThreshold(img_nucle_rescale,255,cv2.ADAPTIVE_THRESH_MEAN_C,cv2.THRESH_BINARY,79,0)
    #plt.imshow(th2)
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
    #plt.figure()
    #plt.imshow(nuc)
    nuc=nuc>0
    nuc_wat_th = morphology.remove_small_holes(nuc, area)
    nuc_wat_th=morphology.remove_small_objects(nuc_wat_th, min_size=area*0.3,connectivity=2)
    nuc_wat_th=np.uint8(nuc_wat_th)*255
#    plt.figure()
#    plt.title("mask of DAPI image")
    #plt.imshow(nuc_wat_th) 
#    contours_tt , _ = cv2.findContours(nuc_wat_th, 2, 1)
#    nuc2=copy.deepcopy(img_nucle_rescale)
#    cv2.drawContours(nuc2,contours_tt,-1,(0,0,255),3) 
#    plt.imshow(nuc2)
    cv2.imwrite(name+"_nucle_mask.jpg",nuc_wat_th)
    
#    img_nucle_g_m = median(img_nucle_rescale, morphology.disk(5))
#    _,nuc_th_c = cv2.threshold(img_nucle_g_m,5,255,cv2.THRESH_BINARY)
#    nuc_th_c=maximum(minimum(nuc_th_c, disk(5)), disk(5))
#    nuc_th_c = morphology.remove_small_holes(nuc_th_c, area*0.7)
#    nuc_th_c=np.uint8(nuc_th_c)*255
    #plt.imshow(nuc_th_c)
#    del img_nucle_g_m
    #==========DIC==img====================================================== 
    name=num+'_BF.tif'  
    img = cv2.imread(name)
#    plt.figure()
#    plt.title("original DIC image")
#    plt.imshow(img) 
    #print(img.shape)
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    # Adaptive Equalization
    gray_adapteq = exposure.equalize_adapthist(gray, clip_limit=0.03)  # 自适应均衡
#    plt.imshow(gray_adapteq)
    gray_adapteq=img_as_ubyte(gray_adapteq)
    median_g = cv2.medianBlur(gray_adapteq, 21)
    #plt.imshow(median_g)

    #median_g=cv2.cvtColor(median_img, cv2.COLOR_BGR2GRAY)
#    gauss = cv2.GaussianBlur(median_g, ( 51,51), 0)
    #plt.imshow(gauss)
    enh = enhance_contrast(median_g, disk(10))
    #plt.imshow(enh)
    #plt.imshow(cv2.cvtColor(enh, cv2.COLOR_BGR2RGB))
    backg =sfr.mean(median_g, disk(80)) 
    #plt.imshow(backg)
    #plt.imshow(cv2.cvtColor(backg, cv2.COLOR_BGR2RGB))
    
    return(img,enh,backg,median_g,img_nucle_g,gray_adapteq,nuc_wat_th,img_nucle_rescale,img_nucle)

def segmentation_gradient(num,img,enh,img_nucle_rescale):
    #======nuc======
    grad_nuc = gradient(img_nucle_rescale, disk(8))
    grad_nuc = cv2.GaussianBlur(grad_nuc, ( 95,95), 0)
   # median(img_nucle_rescale, morphology.disk(5))
    grad_nuc=grad_nuc>35
    summ=(grad_nuc==True).sum()
    if summ>3509452:#70%
        grad_nuc=grad_nuc>45
    if summ>4261478:#85%
        grad_nuc=grad_nuc>55
    grad_nuc = np.uint8(grad_nuc)*255
    #plt.imshow(grad_nuc)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(5,5)) #MORPH_CROSS#MORPH_ELLIPSE
    grad_nuc = cv2.morphologyEx(grad_nuc, cv2.MORPH_CLOSE, kernel,iterations=3)
    nuc_th_c = grad_nuc > 0
    nuc_th_c = morphology.remove_small_holes(nuc_th_c, 3500)
    nuc_th_c = morphology.remove_small_objects(nuc_th_c, min_size=area*0.4,connectivity=2)
    #plt.imshow(nuc_th_c)
    nuc_th_c=np.uint8(nuc_th_c)*255
    
    #====cell image===========
    grad = gradient(enh, disk(8))
#    grad2 = gradient(gauss, disk(8))
#    Added=cv2.addWeighted(gauss, 0.7, grad2, 0.3, 0)
    grad_Gauss = cv2.GaussianBlur(grad, ( 75,75), 0)
    #plt.imshow(Added)
#    plt.imshow(grad2)
   # median(img_nucle_rescale, morphology.disk(5))
    grad_mask=grad_Gauss>25
    grad_mask = np.uint8(grad_mask)*255
    #plt.imshow(grad_mask)
    closed_g = cv2.morphologyEx(grad_mask, cv2.MORPH_CLOSE, kernel,iterations=5)
    #plt.imshow(closed_g) 
    #opened1_g = cv2.morphologyEx(closed_g, cv2.MORPH_OPEN, kernel,iterations=3) 
#    plt.imshow(opened1_g) 
    arr = closed_g > 0
    arr1 = morphology.remove_small_holes(arr, 2000)
    arr2 = morphology.remove_small_objects(arr1, min_size=area*0.4,connectivity=2)
    #plt.imshow(arr2)
    img_th=np.uint8(arr2)*255
    #contours2 , _ = cv2.findContours(img_th, 2, 1)
    #img_cont2=copy.deepcopy(img)
    #cv2.drawContours(img_cont2,contours2,-1,(0,0,255),3) 
    #plt.imshow(img_cont2) 
    #cv2.imwrite(name+"_cell_mask.jpg",img_th)
    return(img_th,nuc_th_c)
# def re_segmentation(imgN,cmaskN,nmaskN,imageP):   
#     r=robjects.r
#     EBImage= importr("EBImage", lib_loc = "C:/Program Files/R/R-3.6.2/library")
#     rscript ='''
#     function(imgN,cmaskN,nmaskN,imageP){
#       setwd(imageP)
#       img = readImage(imgN)
#       r = channel(img,"r")
#       g = channel(img,"g")
#       b = channel(img,"b")
#       gray = 0.21*r+0.71*g+0.07*b
#       #display(gray) 
#       cell_m = readImage(cmaskN)
#       #display(cell_m)
#       nucle_m = readImage(nmaskN)
#       nucle_m=opening(nucle_m, makeBrush(5, shape='disc'))
#       nucle_m=opening(nucle_m, makeBrush(5, shape='disc'))
#       nucle_m = bwlabel(nucle_m)
#       #display(nucle_m,all=TRUE)
#       cells = rgbImage(green=gray+0.3*cell_m, blue=gray+5*nucle_m,red=gray)
#       # cells=img
#       #display(cells, all = TRUE)
#       cell_m = opening(cell_m>0.1, makeBrush(5, shape='disc'))
#       cmask = propagate(gray, seeds=nucle_m, mask=cell_m)
#       #display(cmask, all = TRUE)
#       segmented1 = paintObjects(cmask, cells, col='#ff00ff')
#       segmented = paintObjects(nucle_m, segmented1, col='#ffff00')
#       #display(segmented, all=TRUE)
#       res<-list(cmask=cmask,nmask=nucle_m, seg_cell=segmented1,seg_all=segmented )
#       return(res)
#     }
#     '''
#     Voronoi_segmentation=r(rscript)
#     res=Voronoi_segmentation(imgN,cmaskN,nmaskN,imageP)
#     #print(res[0])
#     res=list(res)
#     pandas2ri.activate()
#     cmask=pandas2ri.ri2py(res[0])
#     nmask=pandas2ri.ri2py(res[1])
#     seg_cell=pandas2ri.ri2py(res[2])
#     seg_all=pandas2ri.ri2py(res[3])
#     
#     cmask = cv2.flip(np.rot90(np.rot90(np.rot90(cmask))),1)
#     nmask = cv2.flip(np.rot90(np.rot90(np.rot90(nmask))),1)
#     seg_cell = cv2.flip(np.rot90(np.rot90(np.rot90(seg_cell))),1)
#     seg_all =cv2.flip(np.rot90(np.rot90(np.rot90(seg_all))),1)
#     #print(seg_all.max())
# #    plt.figure()
# #    plt.imshow(cmask) 
# #    plt.figure()
# #    plt.imshow(nmask) 
# #    plt.figure()
# #    plt.imshow(seg_cell) 
# #    plt.figure()
# #    plt.imshow(seg_all) 
#     return (cmask,nmask,seg_all)

def segmentation_diff(num,img_c2,enh_c2,backg_c2,nuc_wat_th):               
    cell_minus=cv2.absdiff(enh_c2,backg_c2) 
    #hist = cv2.calcHist([cell_minus], [0], None, [256], [0, 256])  #计算直方图
#		plt.plot(hist)
    #plt.imshow(cell_minus)
    _,minus_th = cv2.threshold(cell_minus,45,255,cv2.THRESH_BINARY)
    #minus_th = cv2.adaptiveThreshold(cell_minus,255,cv2.ADAPTIVE_THRESH_MEAN_C,cv2.THRESH_BINARY,79,2)
    overlap1=(minus_th*nuc_wat_th>0)
    overlap_ratio=((overlap1==True).sum())/(((nuc_wat_th>0)==True).sum())
#    plt.imshow(nuc_wat_th)
#    plt.imshow(minus_th)
    if overlap_ratio<0.25:
        _,minus_th = cv2.threshold(cell_minus,20,255,cv2.THRESH_BINARY)
    if overlap_ratio>0.6:
        _,minus_th = cv2.threshold(cell_minus,65,255,cv2.THRESH_BINARY)
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
    cellmask2=np.uint8(arr2)*255
#    plt.imshow(opened2) 
    #contours , _ = cv2.findContours(opened2, 2, 1)
    # B样条   
    
    #cellmask2 = draw_approx_hull_polygon(img_c2, contours)#----------------------------------
    #cellmask2=cellmask2+nuc_wat_th
    #plt.imshow(cellmask2) 
#    fig, axes = plt.subplots(1,2,figsize=(8,8))
#    ax0, ax1= axes.ravel()
#    ax0.imshow(opened2,plt.cm.gray)
#    ax0.set_title('last image')
#    
#    ax1.imshow(cellmask2,plt.cm.gray)
#    ax1.set_title('B-spline image')
    cv2.imwrite(name+"_cell_mask_c2.jpg",cellmask2)
#    plt.imshow(cellmask2)
    #contours2 , _ = cv2.findContours(cellmask2, 2, 1)
    #img_cont1=copy.deepcopy(img_c2)
    #cv2.drawContours(img_cont1,contours2,-1,(0,0,255),3) 
    #plt.imshow(img_cont1) 
    #cv2.imwrite(name+"_cell_mask_c2222.jpg",img_cont1)
    return (cellmask2,cell_minus)

def draw_approx_hull_polygon(img, cnts):
    # img = np.copy(img)
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
    #plt.figure()
    #plt.title('subdivide_polygon',fontsize='large',fontweight='bold') 
   # plt.imshow(img1) 
    return img1

#def compute_centroids(X):
#    n,_ = X.shape
#    centroids=[ (sum(X[:,0]) / n),(sum(X[:,1]) / n)]
#    return(centroids)
def check( cmask, nmask, rowMove,colMove):
    tform = tf.SimilarityTransform(translation=(-(rowMove),(colMove)))
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
    #img_th_c = copy.deepcopy(img_th)
    #img_th_c = tf.warp(img_th,tform)
    #plt.imshow(img_th_c)
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
    
    #plt.imshow(img_th_c)
    img_c= img_as_ubyte(tf.warp(img,tform))
#    plt.figure()
#    plt.imshow(img_c)
    io.imsave(name+'_c.jpg', img_c)
    # contours4 , _ = cv2.findContours (nuc_wat_th , cv2.RETR_EXTERNAL , cv2.CHAIN_APPROX_SIMPLE )
    # img_cont3=copy.deepcopy(img_c)
    # cv2.drawContours(img_cont3,contours4,-1,(0,0,255),3) 
#    plt.imshow(img_cont3)   
    enh_c = img_as_ubyte(tf.warp(enh,tform))
    backg_c=img_as_ubyte(tf.warp(backg,tform))
    return(img_c,img_th_c,enh_c,backg_c,tform)
def compute_centroids(X):
    n,_ = X.shape
    centroids=[ (sum(X[:,0]) / n),(sum(X[:,1]) / n)]
    return(centroids)   
def adjust_img_center(img,enh,cmask,nmask,name,backg):
    #label_cell_rgb = color.label2rgb(cmask)
    #plt.figure()
    #plt.imshow(label_cell_rgb)
    info_cell=measure.regionprops(cmask)
    
    #label_nuc_rgb = color.label2rgb(nmask)
    #plt.figure()
    #plt.imshow(label_nuc_rgb)
    info_nuc=measure.regionprops(nmask)
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
    tform = tf.SimilarityTransform(translation=(-(move_x+5),-(move_y-2)))
    img_c= img_as_ubyte(tf.warp(img,tform))
#    plt.figure()
    #plt.imshow(img_c)
    io.imsave(name+'_c2.jpg', img_c)
#    contours4 , _ = cv2.findContours (nuc_wat_th , cv2.RETR_EXTERNAL , cv2.CHAIN_APPROX_SIMPLE )
#    img_cont3=copy.deepcopy(img_c)
#    cv2.drawContours(img_cont3,contours4,-1,(0,0,255),3) 
#    plt.imshow(img_cont3)   
    cmask_c = img_as_ubyte(tf.warp(cmask,tform))
    enh_c = img_as_ubyte(tf.warp(enh,tform))
    backg_c=img_as_ubyte(tf.warp(backg,tform))
    return(img_c,cmask_c,enh_c,backg_c,tform)

def yellow(num,nmask2,gray_adapteq,img_nucle_rescale,tform):
    yel_n=num+"_TRITC.tif"
    img_yellow=cv2.imread(yel_n)
#    plt.figure()
#    plt.title("original TRITC image")
    #plt.imshow(img_yellow)
#    plt.imshow(cv2.cvtColor(img_yellow, cv2.COLOR_BGR2RGB))
    img_yellow_g = cv2.cvtColor(img_yellow, cv2.COLOR_BGR2GRAY)
    info_yellow=measure.regionprops(nmask2,intensity_image=img_yellow_g)
    img_yellow_backg =sfr.mean(img_yellow_g, disk(300))
#    plt.imshow(img_yellow_backg)
    info_yellow_backg=measure.regionprops(nmask2,intensity_image=img_yellow_backg)

    #gray_adapteq_C=gray_adapteq#cv2.cvtColor(gray_adapteq, cv2.COLOR_GRAY2BGR)
    #img_nucle_B = np.zeros((2048, 2448,3)).astype("uint8")
    #img_nucle_B[:,:,0] = img_nucle_rescale
#plt.imshow(cv2.cvtColor(img_nucle_B, cv2.COLOR_BGR2RGB))
    Added_image=gray_adapteq
    #plt.figure()
    #plt.imshow(cv2.cvtColor(Added_image, cv2.COLOR_BGR2RGB))
    return(info_yellow,info_yellow_backg,Added_image)
    
    
def gaussian(sigma, x, u):
	y = np.exp(-(x - u) ** 2 / (2 * sigma ** 2)) / (sigma * math.sqrt(2 * math.pi))
	return y
	
	
#--CTC_estimate-------------
def CTC_estimate(cmask2,nmask2,info_yellow,info_yellow_backg,seg_all2,img_nucle_rescale,Added_image,num,xls_path):
    info_cell=measure.regionprops(cmask2)
    info_nuccc=measure.regionprops(nmask2,intensity_image=img_nucle_rescale)
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
    
    #savefile=num+'_result.png'
    #plt.savefig(savefile,dpi=600)
    #plt.show(block=True)
    return(seg_CTC,CTC_num,All_num)



def cell_detection(num,imageP):
    global name
    num='1'
    name=num+'_BF.tif'
    img,enh,backg,median_g,img_nucle_g,gray_adapteq,nuc_wat_th,img_nucle_rescale,img_nucle=preprocessing(num,name)
    img_th,nuc_th_c=segmentation_gradient(num,img,enh,img_nucle_rescale)
    img_c,img_th_c,enh_c,backg_c,tform=adjust_img_ovelap(img,enh,img_th,nuc_th_c,name,img_nucle_g,backg,nuc_wat_th)
    segmentation_diff(num,img_c,enh_c,backg_c,nuc_wat_th)
    return(img_c,enh_c,backg_c,nuc_wat_th,tform,gray_adapteq,img_nucle_rescale)

def CTC_e(seg_all,num,cmask,nmask,gray_adapteq,img_nucle_rescale,tform,xls_path):
    cmask = cv2.flip(np.rot90(np.rot90(np.rot90(cmask))),1)
    nmask = cv2.flip(np.rot90(np.rot90(np.rot90(nmask))),1)
    #seg_cell = cv2.flip(np.rot90(np.rot90(np.rot90(seg_cell))),1)
    seg_all =cv2.flip(np.rot90(np.rot90(np.rot90(seg_all))),1)
    info_yellow,info_yellow_backg,Added_image=yellow(num,nmask,gray_adapteq,img_nucle_rescale,tform)
    seg_CTC,CTC_num,All_num=CTC_estimate(cmask,nmask,info_yellow,info_yellow_backg,seg_all,img_nucle_rescale,Added_image,num,xls_path)
    return(seg_CTC,CTC_num,All_num)
    

#os.chdir('C:/Users/10158/OneDrive - City University of Hong Kong/ProgramCode/CTC_Image/02_shiny/CTC/CTC_shiny')
#cell_detection('1','C:/Users/10158/OneDrive - City University of Hong Kong/ProgramCode/CTC_Image/02_shiny/CTC/CTC_shiny')




#     imgN1=name+'_c.jpg'
#     cmaskN1=name+"_cell_mask_c2.jpg"
#     nmaskN=name+"_nucle_mask.jpg"
#     cmask1,nmask1,seg_all1=re_segmentation(imgN1,cmaskN1,nmaskN,imageP)
#     #plt.imshow(seg_all1)


# def cell_detection2(img_c,enh_c,cmask1,nmask1,name,backg_c):
#     img_c2,cmask_c2,enh_c2,backg_c2,tform2=adjust_img_center(img_c,enh_c,cmask1,nmask1,name,backg_c)
#     _,cell_minus=segmentation_diff(num,img_c2,enh_c2,backg_c2,nuc_wat_th)
#     return(img_c2,cmask_c2,enh_c2,backg_c2,tform2,cell_minus)
#     
    
#     imgN2=name+'_c2.jpg'
#     cmaskN2=name+"_cell_mask_c2.jpg"
#     cmask2,nmask2,seg_all2=re_segmentation(imgN2,cmaskN2,nmaskN,imageP)
#     #plt.imshow(seg_all2)

#     info_yellow,info_yellow_backg,Added_image=yellow(num,nmask2,gray_adapteq,img_nucle_rescale,tform,tform2)
#     seg_CTC,CTC_mask,Blue_dia,Cell_dia,Yel_int,Blue_int=CTC(cmask2,nmask2,info_yellow,info_yellow_backg,seg_all2,img_nucle_rescale,Added_image,cell_minus)
#     return(seg_all1,seg_all2,nmask2,cmask2,img_c2,seg_CTC,CTC_mask,Blue_dia,Cell_dia,Yel_int,Blue_int)


# def CTC_detection():
#     global name
#     num='1'
#     name=num+'_BF.tif'
#     print("processing")
#     _,seg_all2,_,_,_,seg_CTC,_,Blue_dia,Cell_dia,Yel_int,Blue_int=cell_detection(num,imageP)
#     print("finish CTC detection")
#     return(seg_all2,Blue_dia,Cell_dia,Yel_int,Blue_int)


