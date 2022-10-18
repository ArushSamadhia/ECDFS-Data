# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 15:43:14 2022

@author: samad
"""

from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
from astropy import units as u
import pandas as pd
import numpy as np
from astropy.table import join, Table, hstack, vstack

p1 = "C:/Users/samad/OneDrive/Desktop/KMooley Internship/Updated Epochs/Epoch1.ecsv"
e1 = ascii.read(p1, delimiter = ' ')
ra1 = e1['RA']
dec1 = e1['DEC']
c1 = SkyCoord(ra1, dec1, frame='icrs',unit= 'deg')
RA = []
DEC = []
co = []
CO = []
def unique(x):
    p2 = "C:/Users/samad/OneDrive/Desktop/KMooley Internship/Updated Epochs/Epoch"+str(x)+".ecsv"
    e2 = ascii.read(p2, delimiter = ' ')
    
    
    ra2 = e2['RA']
    dec2 = e2['DEC']
    
    c2 = SkyCoord(ra2, dec2, frame='icrs',unit= 'deg')
    
    idx,d2d,d3d = c1.match_to_catalog_sky(c2)
    radius = (1/3600)*u.deg
    sep = d2d < radius
    sep1 = d2d >radius
    c1_match = c1[sep]
    c2_match = c2[idx[sep]]
    c2_match_1 = c2[idx[sep1]]
    
    index1 = []
    for i in idx:
        if i not in index1:
            index1.append(int(i))
    index1.sort()
    #print(index1)
    
    '''
    index2 = []
    for j in range(len(c2)):
        if j not in index1:
            index2.append(j)
    print(index2)
    '''
    idx1,d2d1,d3d1 = c2.match_to_catalog_sky(c1)
    sep2 = d2d1 > radius
    c1_match_1 = c1[idx1[sep2]]
    
    for l in c2_match_1:
        co.append(l)
    
    for k in c2_match:
        co.append(k)
        
    if x == 2:
        for m in c1_match_1:
            co.append(m)
            
    '''
    for m in index2:
        RA.append(c2[m].ra.degree)
        DEC.append(c2[m].dec.degree)
    '''
for q in range(2, 10):
    unique(q)
for z in co:
    if z not in CO:
        CO.append(z)
for y in CO:
    RA.append(y.ra.degree)
    DEC.append(y.dec.degree)

#print(RA)
#print(DEC)
coordinates = Table()
coordinates = hstack([coordinates, RA])
coordinates = hstack([coordinates, DEC])
coordinates.rename_column('col0_1', 'RA')
coordinates.rename_column('col0_2', 'DEC')
coordinates.write("C:/Users/samad/OneDrive/Desktop/KMooley Internship/Unique_Sources.ecsv", overwrite = True)