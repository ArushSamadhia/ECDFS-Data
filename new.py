# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 18:50:17 2022

@author: samad
"""
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from astropy.table import Table, hstack, QTable
import matplotlib.pyplot as plt

def unique():
    path = "C:/Users/samad/OneDrive/Desktop/KMooley Internship/Epochs_Dates/Epoch1.ecsv"
    path_err = "C:/Users/samad/OneDrive/Desktop/KMooley Internship/Epochs_error/Epoch_err1.ecsv"
    d_path = "C:/Users/samad/OneDrive/Desktop/KMooley Internship/Dates.txt"
    date = ascii.read(d_path, delimiter = ' ')
    file = ascii.read(path, delimiter = ' ')
    file_err = ascii.read(path_err, delimiter = ' ')
    sample = file[0]
    #print(sample)
    sample_err = file_err[0]
    ra1 = sample['RA']
    dec1 = sample['DEC']
    c1 = SkyCoord(ra1, dec1, frame='icrs', unit='deg')
    radius= 1/3600
    RA1 = []
    DEC1= []
    for r in range(len(sample)-2):
        RA1.append(sample['RA'])
        DEC1.append(sample['DEC'])
    final = Table()
    for i in range(1,10):
        p = "C:/Users/samad/OneDrive/Desktop/KMooley Internship/Epochs_Dates/Epoch"+str(i)+".ecsv"
        p_err = "C:/Users/samad/OneDrive/Desktop/KMooley Internship/Epochs_error/Epoch_err"+str(i)+".ecsv"
        if p == path:
            continue
        ep= ascii.read(p, delimiter = ' ')
        ep_err = ascii.read(p_err, delimiter = ' ')
        ra2 = ep['RA']
        dec2 = ep['DEC']
        c2 = SkyCoord(ra2,dec2, frame='icrs',unit= 'deg')
        idx,d2d,d3d = c1.match_to_catalog_sky(c2)
        #print(d2d<radius*u.deg)
        match_index = idx
        #print(match_index)
        RA = ep['RA']
        DEC = ep['DEC']
        
        cat = ep[match_index]
        cat_err = ep_err[match_index]
        #print(i)
        #print(cat)
        #print(len(cat))
        for q in (d2d>radius*u.deg):
            if q == True:
                for r in range(len(cat)-2):
                    RA1.append(cat['RA'])
                    DEC1.append(cat['DEC'])
        del ep['RA']
        del ep['DEC']
        del ep_err['RA']
        del ep_err['DEC']
        for q in (d2d>radius*u.deg):
            if q == True:
                sample = hstack([sample, cat])
                sample_err = hstack([sample_err, cat_err])
    print(sample)
    del sample['RA']
    del sample['DEC']
    del sample_err['RA']
    del sample_err['DEC']
    #print(RA1)
    #print(DEC1)
    date_j = []
    peaks = []
    peaks_err = []
    pnt = []
    epoch = []
    for j in sample.colnames:
        #print(j)
        date_j.append(int(Time(j).mjd))
        peaks.append(float(sample[j][0]))
        peaks_err.append(float(sample_err[j][0]))
        for k in range(len(date)):
            if date['col1'][k] == j:
                pnt.append(int(date['col2'][k][0]))
                epoch.append(int(date['col2'][k][2]))

    #print(pnt)
    #print(epoch)
    #print(date_j)
    #print(peaks_err)
    curve = QTable()
    curve = hstack([curve, RA1])
    curve = hstack([curve, DEC1])
