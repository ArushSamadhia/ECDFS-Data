# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 12:50:00 2022

@author: samad
"""

from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from astropy.table import Table, hstack, QTable
import matplotlib.pyplot as plt

d_path = "C:/Users/samad/OneDrive/Desktop/KMooley Internship/Dates.txt"
date = ascii.read(d_path, delimiter = ' ')
coord = ascii.read("C:/Users/samad/OneDrive/Desktop/KMooley Internship/PNT_Coord.txt", delimiter= ' ')

path = "C:/Users/samad/OneDrive/Desktop/KMooley Internship/Epochs_Dates/Epoch1.ecsv"
path_err = "C:/Users/samad/OneDrive/Desktop/KMooley Internship/Epochs_error/Epoch_err1.ecsv"
file = ascii.read(path, delimiter = ' ')
file_err = ascii.read(path_err, delimiter = ' ')
sample = file[0]
sample_err = file_err[0]
ra1 = sample['RA']
dec1 = sample['DEC']
c1 = SkyCoord(ra1, dec1, frame='icrs', unit='deg')
radius= 1/3600

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
    del ep['RA']
    del ep['DEC']
    del ep_err['RA']
    del ep_err['DEC']
    cat = ep[match_index]
    cat_err = ep_err[match_index]
    #print(cat)
    
    for i in (d2d<radius*u.deg):
        if i == True:
            sample = hstack([sample, cat])
            sample_err = hstack([sample_err, cat_err])
#print(sample_err)
del sample['RA']
del sample['DEC']
del sample_err['RA']
del sample_err['DEC']
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
curve['PNT'] = pnt
curve['Epoch'] = epoch
curve['MJD'] = date_j

curve = hstack([curve, RA])
curve = hstack([curve, DEC])
CF = []
for o in range(len(curve)):
    if curve['PNT'][o] == 1:
        rac = coord['RA'][0]
        decc = coord['DEC'][0]
    elif curve['PNT'][o] == 2:
        rac = coord['RA'][1]
        decc = coord['DEC'][1]
    elif curve['PNT'][o] == 3:
        rac = coord['RA'][2]
        decc = coord['DEC'][2]
    elif curve['PNT'][o] == 4:
        rac = coord['RA'][3]
        decc = coord['DEC'][3]
    elif curve['PNT'][o] == 5:
        rac = coord['RA'][4]
        decc = coord['DEC'][4]
    elif curve['PNT'][o] == 6:
        rac = coord['RA'][5]
        decc = coord['DEC'][5]
    C1 = SkyCoord(curve['RA'][o], curve['DEC'][o], frame = 'icrs', unit = 'deg')
    C2 = SkyCoord(rac, decc, frame = 'icrs', unit = 'deg')
    off = C2.position_angle(C1).to(u.deg) 
    offset = off.arcminute
    F_GHz=1.4
    G1,G2,G3 = -1.343e-3, 6.579e-7, -1.186e-10 # VLA at 1.465 GHz
    X = offset * F_GHz
    beam = 1 + G1*X**2 + G2*X**4 + G3*X**6
    CF.append(beam)
for p in range(len(curve)):
    peaks[p]/=CF[p]
    peaks_err[p]/=CF[p]
curve['Flux D'] = peaks
curve['Error'] = peaks_err
#print(CF)
print(curve)
curve.write("C:/Users/samad/OneDrive/Desktop/KMooley Internship/Light_Curve.ecsv", overwrite = True)
plt.scatter(date_j, peaks)
plt.errorbar(date_j, peaks, yerr=peaks_err, fmt="o")
    



    