#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 15:38:42 2019

@author: rey
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.signal

#%%Lecture des fichiers de la journée

annee=["2019"]
dates=["20190429"]

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


data532paral=[[0]*len(altitude) for k in range(0,143)]
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
                            col1=file[:,1]
                            col1reverse=col1[::-1]                    
                            data532paral[i]=col1reverse[:]
                            i=i+1

#%%

time=[]    
for i in range(len(data532paral)):
    time.append(float(hours[i])+(float(minutes[i])*60+float(secondes[i]))/3600)
#%%read fonction recouvrement 532 + data corrigé

trajet_geom='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/fgeom/milan_20190121_fgeom_v2.txt'
file=np.loadtxt(trajet_geom,skiprows=13)
fgeom532=file[:,1]
fgeom808=file[:,3]

geom532_fit=list(fgeom532[::-1])
geom808_fit=list(fgeom808[::-1])

for i in range(len(altitude)-len(geom532_fit)):
    geom532_fit.append(1)
    geom808_fit.append(1)

    
#%%Moyenne sur 30min
"""
time_moy=[0]
for i in range(1,48):
    time_moy.append(time_moy[i-1]+0.5)

data532moy=[[0]*len(altitude) for i in range(0,48)]

for i in range(0,40):
    for k in range(0,3):
        for j in range(len(altitude)):
            data532moy[i][j]=data532moy[i][j]+data532paral[k+3*i][j]
                
for i in range(len(data532moy)):
    for j in range(len(altitude)):
        data532moy[i][j]=data532moy[i][j]/3

plt.figure()
plt.title("Verification moyenne")
plt.plot(data532paral[8*6+2][:300],altitude[:300],label="normal")   
plt.plot(data532moy[16][:300],altitude[:300],label="moyenne")
plt.grid()
plt.legend()
plt.show()

"""        

#%% 
"""
for i in range(len(time_moy)):
    for j in range(len(altitude)):
        data532moy[i][j]=data532moy[i][j]/geom532[j]
"""        
        
plt.figure()
plt.plot(data532paral[102][:150],altitude[:150])
plt.grid()
plt.show()

#%%recherche du max

hcla=[]
iymax=[]
time_hcla=time[42:114]

for i in range(42,114): #de 6h à 17h
    y=list(data532paral[i][:150])
    iymax.append(y.index(max(y)))
    hcla.append(altitude[iymax[i-42]])


#%% plot hcla
    
plt.figure()
plt.plot(time_hcla[:],hcla)
plt.title(dates[0])
plt.xlabel("time (UTC)")
plt.ylabel("hcla (km)")
plt.grid()
plt.show()

#%%Save data
hcla_string=[]
time_hcla_string=[]

fichier = open("/net/nfs/home/rey/Documents/code_Zoe/detection_CLA/stat/dev_cla_"+dates[0]+".txt", "w")
fichier.write("time(UTC)|hcla(km)\n")
fichier.write("---------------\n")

for i in range(len(time_hcla)):
    hcla_string.append("%.4f" %hcla[i])
    time_hcla_string.append("%.4f" %time_hcla[i])
    fichier.write(time_hcla_string[i]+"\t"+hcla_string[i]+"\n")
fichier.close()

"""
file=np.loadtxt("/net/nfs/home/rey/Documents/code_Zoe/detection_CLA/stat/dev_cla_"+dates[0]+".txt",skiprows=2)
time2=file[:,0]
hcla2=file[:,1]
"""

#%% dérivée
deriv_hcla=[]

for k in range(len(time_hcla)-1):
        deriv_hcla.append((hcla[k+1]-hcla[k])/time[1])
        
plt.figure()
plt.plot(time_hcla[10:len(time_hcla)-1],deriv_hcla[10:])
plt.title(dates[0])
plt.xlabel("time (UTC)")
plt.ylabel("dérivée hcla (km/UTC)")
plt.grid()
plt.show()

#%%Approximation hcla par un polynôme (désolé François, encore un !)


x=time_hcla[:30]
y=hcla[:30]
coeff_poly_hcla=np.polyfit(x,y,1)
hcla_poly=[]
hcla_poly_deriv=[]

for i in range(len(x)):
    hcla_poly.append(coeff_poly_hcla[0]*x[i]+coeff_poly_hcla[1])
    hcla_poly_deriv.append(coeff_poly_hcla[0])

plt.figure()
plt.plot(time_hcla,hcla)
plt.plot(x,hcla_poly,label="Approximation polynomiale")
plt.xlabel("time (h)")
plt.ylabel("hcla (km)")
plt.legend()
plt.grid()
plt.show()

"""
plt.figure()
plt.plot(x,hcla_poly_deriv)
plt.xlabel("time (h)")
plt.ylabel("Vitesse développement cla (km/h)")
plt.legend()
plt.grid()
plt.show()
"""
print("Vitesse moyenne de developpement de la couche limite (m/s)=",np.mean(hcla_poly_deriv)*1000/3600)