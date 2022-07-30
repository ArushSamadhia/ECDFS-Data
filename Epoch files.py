# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 16:29:11 2022

@author: samad
"""

from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
from astropy import units as u
import pandas as pd
import numpy as np
from astropy.table import join, Table, hstack

def mapping(a):
    print('Epoch'+str(a))
    print()
    if a==1:
        d = 2
        e = 1
    elif a == 9:
        d = 2
        e = 3
    else:
        d = 1
        e = 1
    length = 200
    for i in range(d,7,e):
        path = "C:/Users/samad/OneDrive/Desktop/KMooley Internship/Data/PNT"+str(i)+"_"+str(a)+".ecsv"
    
        ref_catalog = ascii.read(path,delimiter=' ')
        if len(ref_catalog)<length: 
            ref = ref_catalog
            ref_path = path
            length = len(ref_catalog)
        
    #print(ref)
    ref_cat = ref['col4','col5','col10']
    ref_cat_err = ref['col4','col5','col11']
    ref_cat.rename_column('col10', 'peak1')
    ref_cat_err.rename_column('col11', 'peak1_err')
    
    ra1 = ref_cat['col4']
    dec1 = ref_cat['col5']
    c1 = SkyCoord(ra1, dec1, frame='icrs', unit='deg')

    #ref_cat['peak2'] = None
    #print(ref_cat)
    radius = 1/3600
    n=2
    for i in range(d,7,e):
        
        path1 = "C:/Users/samad/OneDrive/Desktop/KMooley Internship/Data/PNT"+str(i)+"_"+str(a)+".ecsv"
        cat1 = ascii.read(path1, delimiter=' ')
        cat_1 = cat1['col4','col5','col10']
        cat_1_err = cat1['col4','col5','col11']
        
        ra2 = cat1['col4']
        dec2 = cat1['col5']
        c2 = SkyCoord(ra2,dec2, frame='icrs',unit= 'deg')
        idx,d2d,d3d = c1.match_to_catalog_sky(c2)
        selection = (d2d > radius*u.deg)
        match_index = idx
            
        match_index[selection] = -1
            #print(match_index)
        s2 = (match_index >= 0)
            
            
        j = 'col10'
        ind = (s2)
            #print(ind)
            
        ind2 = ((cat_1[j][match_index][ind]>0) & (ref_cat['peak1'][ind]>0))
        ind2_err = ((cat_1_err['col11'][match_index][ind]>0) & (ref_cat_err['peak1_err'][ind]>0))
        b1 = cat_1[match_index][ind][ind2]
        b1_err = cat_1_err[match_index][ind][ind2_err]
        b1.rename_column('col10', 'peak'+str(n))
        b1_err.rename_column('col11', 'peak'+str(n))
            #print(b1)
        b11 = b1['peak'+str(n)]
        b11_err = b1_err['peak'+str(n)]
        b2 = ref_cat[ind][ind2]
        b2_err = ref_cat_err[ind][ind2_err]
        if i == 1 or (i==2 and (a==1 or a==9)) :
                if path1 == ref_path:
                    final= Table()
                    final_err = Table()
                    n-=1
                else:
                    final = hstack([b2,b11])
                    final_err = hstack([b2_err, b11_err])
                    
        else:
                if len(final) == 0:
                    final = hstack([b2,b11])
                    final = hstack([b2_err, b11_err])
                else:
                    if path1 != ref_path:
                        final = hstack([final,b11])
                        final_err = hstack([final_err, b11_err])
        n+=1
            #print(final)
    final.rename_column('col4', 'RA')
    final.rename_column('col5', 'DEC')
    print(final_err)
    
    #final.write("C:/Users/Arush Samadhia/Desktop/ECDFS Data/Final_Output_"+str(a)+".ecsv",overwrite= True)
        
for m in range(1,10):        
   mapping(m)