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
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
Blue_dia_all=[]
Cell_dia_all=[]
Yel_int_all=[]
Blue_int_all=[]
# white_blood
#result_Path="/dataserver145/image/Yuni/CTC/01.2_cellline_samples/01_50tumor_50whiteBlood/02.1_whiteBlood_samples/02_segmented/"
#50WBC+50CMF
result_Path="/dataserver145/image/Yuni/CTC/01.2_cellline_samples/01_50tumor_50whiteBlood/02.2_50CMF+50whiteBlood_samples/02_segmented/"

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
