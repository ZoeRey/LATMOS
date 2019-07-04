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

annee=["2018"]
dates=["20180929"]

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


data532paral=[[0]*len(altitude) for k in range(0,120)]
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
 

#%%
plt.figure()
plt.plot(data532paral[75][:300],altitude[:300])
plt.show()


#%% recherche du max
iymax=[]
   


for i in range(40,75): #de 8h à 16h
    if i<=49: #de 8h à 10h
        y=list(data532paral[i][:70])
        iymax.append(y.index(max(y)))
    else :
        y=list(data532paral[i][:110])
        iymax.append(y.index(max(y)))

#%%filtre médian
    
time_hcla=time[40:75]
data532_medianfilter=[[] for i in range(len(time_hcla))]
xbis=[[] for i in range(len(time_hcla))]

for i in range(len(iymax)):
    if i<=2:
        data532_medianfilter[i]=data532_medianfilter[i]+list(data532paral[i+40][:30])
        xbis[i]=xbis[i]+list(altitude[:30])      
    elif i>2 and i<=9:
        data532_medianfilter[i]=data532_medianfilter[i]+list(data532paral[i+40][:60])
        xbis[i]=xbis[i]+list(altitude[:60])       
    else:
        data532_medianfilter[i]=data532_medianfilter[i]+list(scipy.signal.medfilt(data532paral[i+40][iymax[i]-15:iymax[i]+50]))
        xbis[i]=xbis[i]+list(altitude[iymax[i]-15:iymax[i]+50])
        

#%% Approximation polynôme de degré 5 et calcul dérivé 2nd de ce polynome

datafit=[[] for i in range(len(time_hcla))]
datafitderiv2=[[] for i in range(len(time_hcla))]
datafitderiv1=[[] for i in range(len(time_hcla))]

for i in range(len(time_hcla)):
    print(i)
    x=xbis[i]
    y=data532_medianfilter[i]
    coeff_poly=np.polyfit(x,y,5)

    for j in range(len(x)):
        datafit[i].append(coeff_poly[0]*x[j]*x[j]*x[j]*x[j]*x[j]+coeff_poly[1]*x[j]*x[j]*x[j]*x[j]+coeff_poly[2]*x[j]*x[j]*x[j]+coeff_poly[3]*x[j]*x[j]+coeff_poly[4]*x[j]+coeff_poly[5])
        datafitderiv1[i].append(5*coeff_poly[0]*x[j]*x[j]*x[j]*x[j]+4*coeff_poly[1]*x[j]*x[j]*x[j]+3*coeff_poly[2]*x[j]*x[j]+2*coeff_poly[3]*x[j]+coeff_poly[4])
        datafitderiv2[i].append(20*coeff_poly[0]*x[j]*x[j]*x[j]+12*coeff_poly[1]*x[j]*x[j]+6*coeff_poly[2]*x[j]+2*coeff_poly[3])
        
    plt.figure()
    plt.plot(datafit[i],x,color='black',linewidth=2.5,linestyle='--',label="approximation")
    #plt.plot(data532paral[i+49][iymax[i]-4:iymax4[i]+12],altitude[iymax[i]-4:iymax4[i]+12],label="S 532//")
    plt.plot(y,x)
    plt.xlabel("signal rétrodiffusé, 532// nm")
    plt.ylabel("altitude (km)")
    plt.legend()
    plt.grid()
    plt.show()
    

#%% Recherche point d'inflexion
#il peut y avoir plusieurs points d'inflexion

pt_inflexion=[[] for i in range(len(time_hcla))]

for i in range(len(time_hcla)):
    x=xbis[i]
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
print("A "+str(time[75])+" le sommet de la couche limite était plutôt à :")
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
 
plt.figure()
plt.plot(time_hcla,hcla)
plt.xlabel("time (UTC)")
plt.ylabel("hcla (km)")
plt.title(dates[0])
plt.grid()
plt.show()

#%%filtre médian sur hcla
    
hcla_medianfilter=list(scipy.signal.medfilt(hcla))

plt.figure()
plt.plot(time_hcla,hcla_medianfilter)
plt.xlabel("time (h)")
plt.ylabel("hcla (km)")
plt.grid()
plt.show()

#%% test en partant du début

hcla2=[]
#poser la question !!!!
print("A "+str(time[40])+", le sommet de la couche limite était plutôt à :")
print(pt_inflexion[0])
hcla_8h=input("Choisissez: ")
hcla_8h=int(hcla_8h)
hcla2.append(pt_inflexion[0][hcla_8h])
for i in range(1,len(pt_inflexion)):
    diff1=100.
    for j in range(len(pt_inflexion[i])):
        diff=abs(pt_inflexion[i][j]-hcla2[i-1])
        if diff<diff1:
            diff1=diff
            jgood=j
    hcla2.append(pt_inflexion[i][jgood])  

    
hcla2_medianfilter=list(scipy.signal.medfilt(hcla2))

plt.figure()
plt.plot(time_hcla,hcla2)
plt.xlabel("time (UTC)")
plt.ylabel("hcla (km)")
plt.title(dates[0])
plt.legend()
plt.grid()
plt.show()


#%%Save data
hcla_string=[]
time_hcla_string=[]

fichier = open("/net/nfs/home/rey/Documents/code_Zoe/detection_CLA/stat/dev_cla_"+dates[0]+".txt", "w")
fichier.write("time(UTC)|hcla(km)\n")
fichier.write("---------------\n")

for i in range(len(time_hcla)):
    hcla_string.append("%.4f" %hcla2[i])
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
plt.plot(time_hcla[:len(time_hcla)-1],deriv_hcla[:])
plt.title(dates[0])
plt.xlabel("time (UTC)")
plt.ylabel("dérivée hcla (km/UTC)")
plt.grid()
plt.show()

#%%Approximation hcla par un polynôme (désolé François, encore un !)


x=time_hcla[:20]
y=hcla[:20]
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

