# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 16:31:04 2022

@author: samad
"""

import requests

index1 = [2,3,4,5,6,7,8]
index2 = [1,2,3,4,5,6,7,8]
index3 = [1,2,3,4,5,6]
index4 = [1,2,3,4,5,6,7,8,9]

def download(a,b):
        for i in b:
            url = 'http://www.tauceti.caltech.edu/kunal/uploads/ecdfs/PNT'+str(a)+'_'+str(i)+'.aegean'
            
            response = requests.get(url)
            
            open("C:/Users/samad/OneDrive/Desktop/KMooley Internship/PNT"+str(a)+"_"+str(i)+".ecsv","w+")
            open("C:/Users/samad/OneDrive/Desktop/KMooley Internship/PNT"+str(a)+"_"+str(i)+".ecsv", "wb").write(response.content)

download(1,index1)
download(2,index4)
download(3,index2)
download(4,index2)
download(5,index4)
download(6,index2)
