# ECDFS-Data
## Introduction
Our universe is filled with a large number of strange objects one of them being transients. Transients refer to the sources whose intensity or the flux emission varies as a function of time. A few of the transients include supernovae, Active Galactic Nuclei (AGN), gamma-ray bursts, etc. 
The project deals with the identification of supernovae by plotting a light curve, i.e., a plot between the flux/intensity and time. If we observe a variable intensity of the source as a function of time, we can say that we have found a transient. But we will need further verifications to conclude whether the observed source is a supernova.
## Dataset
The data on which I am working is the Extended Chandra Deep Field South (ECDFS) data obtained from the Chandra X-Ray observatory. The plotting of the light curve begins by identifying the same sources from different catalogs. For this, I ran a script based on Astropy for catalog cross-matching. Then the flux from the same sources was arranged in a table as a function of time.  
## Script
The script involves first locating the closest source to a particular source in some other catalog. Then a threshold radius in units of arcseconds is set depending on the distance of the source from Earth. If the detected closest source lies within the threshold radius, it is the same source and we add its flux to the table, else we check for the other source matches. This is required as the coordinates the telescope returns for a particular source might vary by a few arcseconds due to the radius of the observed source.
