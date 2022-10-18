# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 18:15:12 2022

@author: samad
"""
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from astropy.table import Table, hstack, QTable

RA = []
DEC = []

def sep(x,y):
    if len(RA) == 0:
        RA.append(x)
        DEC.append(y)
    else:
        C1 = SkyCoord(x, y, unit=(u.hourangle, u.deg))
        for i in range(len(RA)):
            C2 = SkyCoord(RA[i], DEC[i], unit=(u.hourangle, u.deg))
            if C1.separation(C2).arcsecond > 1:
                RA.append(x)
                DEC.append(y)
                

p = "C:/Users/samad/OneDrive/Desktop/KMooley Internship/Epochs_Dates/Epoch1.ecsv"
file = ascii.read(p, delimiter = ' ')
print(len(file))
for k in range(len(file)):
        print(k)
        print(file['RA'][k])
        sep(file['RA'][k], file['DEC'][k])
    #print(RA)
    #print(DEC)
        
print(RA)
print(DEC)
unique = QTable()
unique = hstack([unique, RA])
unique = hstack([unique, DEC])

unique.write("C:/Users/samad/OneDrive/Desktop/KMooley Internship/Unique_Sources_EP1.ecsv", overwrite = True)