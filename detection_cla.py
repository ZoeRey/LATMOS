#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 15:38:42 2019

@author: rey
"""
import numpy as np
import os.path
import matplotlib.pyplot as plt

#%%Lecture des fichiers de la journée

annee=["2019"]
dates=["20190329"]

if annee[0]=="2018" :
    trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/ascii/2018/'+dates[0]+'/mli2MOY_'+dates[0]+'.000000.JUS'
else :
    trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/ascii/2019/'+dates[0]+'/mli2MOY_'+dates[0]+'.000000.JUS'

file=np.loadtxt(trajet,skiprows=36)
altitude=file[:,0]
altitude=altitude[::-1]  
#zone aveugelement jusqu'à 75m
#altitude=altitude[6:]

num1=["0","1","2"]
num2=["0","1","2","3","4","5"]
num3=["0"]
num4=["0","1","2"]
num=["0","1","2","3","4","5","6","7","8","9"]
hours=[]
minutes=[]
secondes=[]

var808paral=[[0]*len(altitude) for k in range(0,143)]
var532paral=[[0]*len(altitude) for k in range(0,143)]
i=0
for a in num1:
    for b in num:
        for c in num2:
            for d in num3:
                for e in num4:
                    for f in num:
                        trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/ascii/'+annee[0]+'/'+dates[0]+'/'+'mli2MOY_'+dates[0]+'.'+a+b+c+d+e+f+'.JUS'
                        if os.path.exists(trajet):
                            #print(i)
                            hours.append(a+b)
                            minutes.append(c+d)
                            secondes.append(e+f)
                            file=np.loadtxt(trajet,skiprows=36)
                            col7=file[:,7]
                            col5=file[:,5]
                            col7reverse=col7[::-1]
                            col5reverse=col5[::-1]
                            var808paral[i]=col7reverse[:]
                            var532paral[i]=col5reverse[:]
                            i=i+1
  
#%%
                      
zmaxvar808paral=[] 
zmaxvar532paral=[]            
for i in range(len(var808paral)):
    var=list(var808paral[i][30:])
    var2=list(var532paral[i][30:])
    zmaxvar808paral.append(altitude[var.index(max(var))])
    zmaxvar532paral.append(altitude[var2.index(max(var2))])

time=[]    
for i in range(len(var808paral)):
    time.append(float(hours[i])+((100./60.)*float(minutes[i])+float(secondes[i]))/100.)

#%%

plt.figure()
plt.title('var532// date')
plt.plot(var532paral[60],altitude,label='11h')
plt.plot(var532paral[61],altitude,label='11h')
plt.plot(var532paral[62],altitude,label='11h')
plt.plot(var532paral[63],altitude,label='11h')
"""
plt.plot(var532paral[78],altitude,label='13h')
plt.plot(var532paral[84],altitude,label='14h')
plt.plot(var532paral[90],altitude,label='15h')
plt.plot(var532paral[96],altitude,label='16h')
"""
plt.ylabel("altitude (km)")
plt.xlabel("var532//")
plt.ylim(0,2)
plt.legend()
plt.show()

 
plt.figure()
plt.title('zmax var532// date')
plt.plot(time,zmaxvar532paral)
plt.xlabel("time")
plt.ylabel("zmax var532//")
#plt.xlim(11,17)
plt.legend()
plt.show()
