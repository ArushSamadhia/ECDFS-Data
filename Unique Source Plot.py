# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 23:32:44 2022

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
import matplotlib.pyplot as plt

path = "C:/Users/samad/OneDrive/Desktop/KMooley Internship/Unique_Sources.ecsv"

file = ascii.read(path)

plt.scatter(file['RA'], file['DEC'])
