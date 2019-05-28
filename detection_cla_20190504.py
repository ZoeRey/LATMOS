#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 15:38:42 2019

@author: rey
"""
import numpy as np
import os.path
import matplotlib.pyplot as plt
import scipy.signal

#%%Lecture des fichiers de la journée

annee=["2019"]
dates=["20190514"]

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
    time.append(float(hours[i])+((100./60.)*float(minutes[i])+float(secondes[i]))/100.)


#%%read fonction recouvrement 532 + data corrigé

trajet_geom='/net/nfs/home/rey/Documents/code_Zoe/MILAN/fct_geom532_JUS.txt'
geom532=list(np.loadtxt(trajet_geom))
geom532[0]=1.

if annee[0]=="2018" :
    for i in range(len(altitude)-len(geom532)):
        geom532.append(1.)
#%% recherche du max et du max/4

plt.plot(data532paral[49][:200],altitude[:200])
plt.grid()
plt.show()

iymax=[]

for i in range(len(time)):
    for j in range(len(altitude)):
        data532paral[i][j]=data532paral[i][j]/geom532[j]

for i in range(49,90):
    diff1=100.
    if i<=60:
        y=list(data532paral[i][:60])
        iymax.append(y.index(max(y)))  
    elif i>60 and i <78:
        y=list(data532paral[i][10:100])
        iymax.append(y.index(max(y))+10)
    elif i>=78:
        y=list(data532paral[i][65:110])
        iymax.append(y.index(max(y))+65)

    
#%%filtre médian
    
time_hcla=time[49:90]
data532_medianfilter=[[] for i in range(len(time_hcla))]

for i in range(len(iymax)):
    if i<=22:
        data532_medianfilter[i]=data532_medianfilter[i]+list(scipy.signal.medfilt(data532paral[i+49][iymax[i]-10:iymax[i]+100]))
    else:
        data532_medianfilter[i]=data532_medianfilter[i]+list(scipy.signal.medfilt(data532paral[i+49][iymax[i]-10:iymax[i]+70]))



"""
plt.plot(data532_medianfilter[80-49],altitude[iymax[80-49]-4:iymax4[80-49]+12])
#plt.plot(data532paral[80][iymax[80-49]-4:iymax4[80-49]+12],altitude[iymax[80-49]-4:iymax4[80-49]+12])
plt.grid()
plt.show()
"""
#%% Approximation polynôme de degré 5 et calcul dérivé 2nd de ce polynome

datafit=[[] for i in range(len(time_hcla))]
datafitderiv2=[[] for i in range(len(time_hcla))]

for i in range(len(time_hcla)):
    if i<=22:
        x=altitude[iymax[i]-10:iymax[i]+100]
    else:
        x=altitude[iymax[i]-10:iymax[i]+70]
    y=data532_medianfilter[i]
    coeff_poly=np.polyfit(x,y,5)

    for j in range(len(x)):
        datafit[i].append(coeff_poly[0]*x[j]*x[j]*x[j]*x[j]*x[j]+coeff_poly[1]*x[j]*x[j]*x[j]*x[j]+coeff_poly[2]*x[j]*x[j]*x[j]+coeff_poly[3]*x[j]*x[j]+coeff_poly[4]*x[j]+coeff_poly[5])
        datafitderiv2[i].append(20*coeff_poly[0]*x[j]*x[j]*x[j]+12*coeff_poly[1]*x[j]*x[j]+6*coeff_poly[2]*x[j]+2*coeff_poly[3])
        
    plt.figure()
    plt.plot(datafit[i],x,color='black',linewidth=2.5,linestyle='--',label="approximation")
    #plt.plot(data532paral[i+49][iymax[i]-4:iymax4[i]+12],altitude[iymax[i]-4:iymax4[i]+12],label="S 532//")
    plt.plot(y,x)
    plt.xlabel("signal rétrodiffusé, 532// nm")
    plt.ylabel("altitude (km)")
    plt.legend()
    #plt.grid()
    plt.show()
    
    
#%% Recherche point d'inflexion
#il peut y avoir plusieurs points d'inflexion

pt_inflexion=[[] for i in range(len(time_hcla))]

for i in range(len(time_hcla)):
    if i<=22:
        x=altitude[iymax[i]-10:iymax[i]+100]
    else:
        x=altitude[iymax[i]-10:iymax[i]+70]    
    for k in range(len(x)-1):
        if(datafitderiv2[i][k+1]>=0 and datafitderiv2[i][k]<=0) or (datafitderiv2[i][k+1]<=0 and datafitderiv2[i][k]>=0):
            pt_inflexion[i].append(x[k])
            
            
#%%Vérif s'il y a des endroits sans points d'inflexion
                
for i in range(len(pt_inflexion)):
    if len(pt_inflexion[i])==0:
        pt_inflexion[i]=pt_inflexion[i-1]       

           
#%%
hcla=[]
#poser la question !!!!
print("A 15h, le sommet de la couche limite était plutôt à :")
print(pt_inflexion[len(pt_inflexion)-1])
hcla_15h=input("Choisissez: ")
hcla_15h=int(hcla_15h)
hcla.append(pt_inflexion[len(pt_inflexion)-1][hcla_15h])
for i in range(len(pt_inflexion)-1,0,-1):
    diff1=100.
    for j in range(len(pt_inflexion[i])):
        diff=abs(pt_inflexion[i][j]-hcla[abs(i-len(pt_inflexion)+1)])
        if diff<diff1:
            diff1=diff
            jgood=j
    hcla.append(pt_inflexion[i][jgood])  
hcla=hcla[::-1]
"""   
plt.figure()
plt.plot(time_hcla,hcla)
plt.xlabel("time")
plt.ylabel("hcla")
plt.grid()
plt.show()
"""
#%%filtre médian sur hcla
    
hcla_medianfilter=list(scipy.signal.medfilt(hcla))

plt.figure()
plt.plot(time_hcla,hcla_medianfilter)
plt.xlabel("time (h)")
plt.ylabel("hcla (km)")
plt.title(dates[0])
plt.legend()
plt.grid()
plt.show()

#%%Approximation hcla par un polynôme (désolé François, encore un !)

x=time_hcla[16:]
y=hcla_medianfilter[16:]
coeff_poly_hcla=np.polyfit(x,y,1)
hcla_poly=[]
hcla_poly_deriv=[]

for i in range(len(x)):
    hcla_poly.append(coeff_poly_hcla[0]*x[i]+coeff_poly_hcla[1])
    hcla_poly_deriv.append(coeff_poly_hcla[0])

plt.figure()
plt.plot(time_hcla,hcla_medianfilter)
plt.plot(x,hcla_poly,label="Approximation polynomiale")
plt.xlabel("time (h)")
plt.ylabel("hcla (km)")
plt.title(dates[0])
plt.legend()
plt.grid()
plt.show()



print("Vitesse moyenne de developpement de la couche limite (m/s)=",np.mean(hcla_poly_deriv)*1000/3600)


