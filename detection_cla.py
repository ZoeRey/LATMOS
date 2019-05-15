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
dates=["20190328"]

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
                            col7=file[:,7]
                            col5=file[:,5]
                            col1=file[:,1]
                            col1reverse=col1[::-1]                    
                            col7reverse=col7[::-1]
                            col5reverse=col5[::-1]
                            var808paral[i]=col7reverse[:]
                            var532paral[i]=col5reverse[:]
                            data532paral[i]=col1reverse[:]
                            i=i+1

#%%
time=[]    
for i in range(len(var808paral)):
    time.append(float(hours[i])+((100./60.)*float(minutes[i])+float(secondes[i]))/100.)

#%%
    
#plt.plot(data532paral[60][:200],altitude[:200])
#plt.plot(data532paral[72][:200],altitude[:200])
#plt.plot(data532paral[84][:200],altitude[:200])
#plt.plot(data532paral[96][:200],altitude[:200])
plt.plot(data532paral[90][:200],altitude[:200])
plt.grid()
plt.show()
    
#%%Détection couche limite avec la recherche du point d'inflexion
iymax=[]
iymax4=[]
hcla=[]
#data532paralfit=[[0]*60 for k in range(len(time))]
#data532paralfitderiv2=[[0]*60 for k in range(len(time))]

for i in range(len(time)):
    y=list(data532paral[i][:200])
    iymax.append(y.index(max(y)))
    
    diff1=100.
    for l in range(iymax[i],200):
        diff=abs(data532paral[i][l]-data532paral[i][iymax[i]]/4)
        if(diff<diff1):
            diff1=diff
            i4=l
    iymax4.append(i4)
    y=data532paral[i][iymax[i]-4:iymax4[i]+10]
    x=altitude[iymax[i]-4:iymax4[i]+10]
    coeff_poly=np.polyfit(x,y,5)
    datafit=[]

    for j in range(len(x)):
        datafit.append(coeff_poly[0]*x[j]*x[j]*x[j]*x[j]*x[j]+coeff_poly[1]*x[j]*x[j]*x[j]*x[j]+coeff_poly[2]*x[j]*x[j]*x[j]+coeff_poly[3]*x[j]*x[j]+coeff_poly[4]*x[j]+coeff_poly[5])
    """plt.figure()
    plt.plot(datafit,x)
    plt.plot(y,x)
    plt.grid()
    plt.show()
    datafitderiv2=[]"""
    for j in range(len(x)):
        datafitderiv2.append(20*coeff_poly[0]*x[j]*x[j]*x[j]+12*coeff_poly[1]*x[j]*x[j]+6*coeff_poly[2]*x[j]+2*coeff_poly[3])


    for k in range(len(x)):
        if (datafitderiv2[k]>=0):
            izero=k
            break

    hcla.append(x[izero])

#plt.plot(time,hcla)
#plt.plot(datafitderiv2[90],altitude[iymax[90]-10:iymax[90]+50])
#plt.grid()    
#plt.show()
    
plt.plot(data532paral[90][iymax[90]-4:iymax4[90]+10],altitude[iymax[90]-4:iymax4[90]+10],label="data532paral")
#plt.plot(data532paralfit[90],altitude[iymax[90]-10:iymax4[90]+10],label="fit")
plt.legend()
plt.grid()
plt.show()     


#plt.plot(data532paral[90][:200],altitude[:200])
#plt.grid()
#plt.show()


plt.plot(time,hcla)
plt.xlabel("time")
plt.ylabel("hcla")
plt.grid()
plt.show()

#%% Marche pas bien, une polynôme de degré 3 est pas assez pentu
hcla_fit=[]
time_hcla=time[42:84]
coeff_poly_hcla=np.polyfit(time_hcla,hcla[42:84],3)
for i in range(len(time_hcla)):
    hcla_fit.append(coeff_poly_hcla[0]*time_hcla[i]*time_hcla[i]*time_hcla[i]+coeff_poly_hcla[1]*time_hcla[i]*time_hcla[i]+coeff_poly_hcla[2]*time_hcla[i]+coeff_poly_hcla[3])

plt.plot(time_hcla,hcla_fit)
plt.plot(time[42:84],hcla[42:84])
plt.show()




#%% Détection couche limite avec la recherche du max de variance
                      
zmaxvar808paral=[] 
zmaxvar532paral=[]            
for i in range(len(var808paral)):
    var=list(var808paral[i][30:])
    var2=list(var532paral[i][30:])
    zmaxvar808paral.append(altitude[var.index(max(var))])
    zmaxvar532paral.append(altitude[var2.index(max(var2))])



#%%

plt.figure()
plt.title('var532// date')
plt.plot(var532paral[60],altitude,label='10h')
plt.plot(var532paral[66],altitude,label='11h')
plt.plot(var532paral[69],altitude,label='11h30')
plt.plot(var532paral[72],altitude,label='12h')
"""
plt.plot(var532paral[78],altitude,label='13h')
plt.plot(var532paral[84],altitude,label='14h')
plt.plot(var532paral[90],altitude,label='15h')
plt.plot(var532paral[96],altitude,label='16h')
"""
plt.ylabel("altitude (km)")
plt.xlabel("var532//")
plt.ylim(0,2)
#plt.xlim(0,10000)
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
