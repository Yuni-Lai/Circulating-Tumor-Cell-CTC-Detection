# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 12:06:38 2020

@author: Yuni
"""

import cv2
from skimage import exposure

def equalize_preprocessing(num,name):  
    #======nucleus================
    nucle_n=num+"_DAPI.tif"
    img_nucle=cv2.imread(nucle_n)
    imgnucle_adapteq = exposure.equalize_adapthist(img_nucle, clip_limit=0.03) 
    #==========DIC==img====================================================== 
    name=num+'_BF.tif'  
    img = cv2.imread(name)
    img_adapteq = exposure.equalize_adapthist(img, clip_limit=0.03)  # 自适应均衡
    yel_n=num+"_TRITC.tif"
    img_yellow=cv2.imread(yel_n)
    yellow_adapteq = exposure.equalize_adapthist(img_yellow, clip_limit=0.03)  
    cv2.imwrite("C:/Users/10158/OneDrive - City University of Hong Kong/ProgramCode/CTC_Image/02_shiny/uplaod/"+name+"_img_adapteq.jpg",img_adapteq)
    cv2.imwrite("C:/Users/10158/OneDrive - City University of Hong Kong/ProgramCode/CTC_Image/02_shiny/uplaod/"+name+"_nucle_adapteq.jpg",imgnucle_adapteq)
    cv2.imwrite("C:/Users/10158/OneDrive - City University of Hong Kong/ProgramCode/CTC_Image/02_shiny/uplaod/"+name+"_yellow_adapteq.jpg",yellow_adapteq)
    return(img_adapteq,imgnucle_adapteq,yellow_adapteq)


