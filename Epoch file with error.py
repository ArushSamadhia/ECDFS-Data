# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 15:06:09 2022

@author: samad
"""

from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
from astropy import units as u
import pandas as pd
import numpy as np
from astropy.table import join, QTable, hstack, vstack

path = "C:/Users/samad/OneDrive/Desktop/KMooley Internship/Dates.txt"
date = ascii.read(path, delimiter = ' ')

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
            ref_i = i
            ref_a = a
            length = len(ref_catalog)
    for j in range(len(date)):
        if date['col2'][j] == str(ref_i)+'_'+str(ref_a):
            ref_date = date['col1'][j]
    print(str(ref_a)+'_'+str(ref_i))      
    #print(ref_path)    
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
        print(i)
        
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
        for k in range(len(date)):
            if date['col2'][k] == str(i)+'_'+str(a):
                d_o = date['col1'][k]
        
        ind2 = ((cat_1[j][match_index][ind]>0) & (ref_cat['peak1'][ind]>0))
        ind2_err = ((cat_1_err['col11'][match_index][ind]>0) & (ref_cat_err['peak1_err'][ind]>0))
        b1 = cat_1[match_index][ind][ind2]
        b1_err = cat_1_err[match_index][ind][ind2_err]
        b1.rename_column('col10', d_o)
        b1_err.rename_column('col11', d_o)
            #print(b1)
        b11 = b1[d_o]
        b11_err = b1_err[d_o]
        b2 = ref_cat[ind][ind2]
        b2_err = ref_cat_err[ind][ind2_err]
        b2.rename_column('peak1', ref_date)
        b2_err.rename_column('peak1_err', ref_date)
        if (i ==2 or i==3 or i==4) and (a == 6 or a == 5) :
            
            print(b1)
            print(b2)
        b1c = SkyCoord(b1['col4'],b1['col5'], frame='icrs',unit= 'deg')
        if i == d:
            final = QTable()
            col = []
            for t in b2:
               # print(t)
                b2c = SkyCoord(t['col4'],t['col5'], frame='icrs',unit= 'deg')
                idx_2,d2d_2,d3d_2 = b2c.match_to_catalog_sky(b1c)
                #Id = idx_2[(d2d_2 < (radius/30)*u.deg)]
                #print(idx_2)
                #print(Id)
                #print(len(b1))
                #print(len(b2))
                col.append(b1[d_o][idx_2])
            final = b2
            final[d_o] = col
        else:
            col = []
            for t in final:
                final_c = SkyCoord(t['col4'],t['col5'], frame='icrs',unit= 'deg')
                idx_2,d2d_2,d3d_2 = final_c.match_to_catalog_sky(b1c)
                col.append(b1[d_o][idx_2])
            final[d_o] = col
    print(final)
    final.rename_column('col4', 'RA')
    final.rename_column('col5', 'DEC')
    final.write("C:/Users/samad/OneDrive/Desktop/KMooley Internship/Updated Epochs/Epoch"+str(a)+".ecsv",overwrite= True)
        
'''
            if i == 1 or (i==2 and (a==1 or a==9)) :
                    if path1 == ref_path:
                        final= Table()
                    else:
                        final = hstack([b2[t], b1['col10'][idx_2]])
            '''
            
'''    
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
                    final_err = hstack([b2_err, b11_err])
                else:
                    if path1 != ref_path:
                        final = hstack([final,b11])
                        final_err = hstack([final_err, b11_err])
                n+=1
          
        if a == 4:
            #print(b2_err)
            #print(final_err)
            
            #print(final)
    
    final.rename_column('col4', 'RA')
    final.rename_column('col5', 'DEC')
    final_err.rename_column('col4', 'RA')
    final_err.rename_column('col5', 'DEC')
  
    final_err.write("C:/Users/samad/OneDrive/Desktop/KMooley Internship/Epochs_error/Epoch_err"+str(a)+".ecsv",overwrite= True)
    '''   
for m in range(1,10):        
   mapping(m)