#Version Zoe
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 16:03:36 2019

Calcul du parametre geometrique et de la constante systeme pour le lidar Milan.

"""
import statistics
import jutilities as juti
import matplotlib.pyplot as plt
import numpy as np
import outils_OHP
from scipy.optimize import curve_fit


#%% Read file

annee=["2018","2018","2018","2018","2018"]
dates=["20180902","20180902","20180902","20180902","20180902"]
#hour=["050000","050100","050201","050301","050402","050502"]
hour=["010000","011006","012010","013014","014018"]
moytmps=["2MOY_","2_","MOY"]



trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/ICOS_OHP/ascii/'+annee[0]+'/'+'CE376_'+dates[0]+'/mli'+moytmps[0]+dates[0]+'.'+hour[0]+'.OHP'
if annee[0]=="2018" :
    skip=36
if annee[0]=="2017":
    skip=33
    
file=np.loadtxt(trajet,skiprows=skip)
altitude=file[:,0]
altitude=altitude[::-1]   

data532=[[0]*len(altitude) for i in range(len(dates))]
data532_dp=[[0]*len(altitude) for i in range(len(dates))]

for i in range(len(dates)):
    trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/ICOS_OHP/ascii/'+annee[i]+'/'+'CE376_'+dates[i]+'/mli'+moytmps[0]+dates[i]+'.'+hour[i]+'.OHP'
    file=np.loadtxt(trajet,skiprows=skip) 
    voie532=file[:,1]
    voie532_dp=file[:,2]
    data532[i]=voie532[::-1]
    data532_dp[i]=voie532_dp[::-1]


data532moy=[0]*len(altitude)
data532dpmoy=[0]*len(altitude)


for i in range(len(altitude)):
    for n in range(len(dates)):
        data532moy[i]=data532moy[i]+data532[n][i]
        data532dpmoy[i]=data532dpmoy[i]+data532_dp[n][i]
    data532moy[i]=data532moy[i]/(len(dates))
    data532dpmoy[i]=data532dpmoy[i]/(len(dates))


plt.figure()
plt.title("Vérification moyenne, 532//")
plt.plot(altitude,data532[1])
plt.plot(altitude,data532moy)
plt.show()

plt.figure()
plt.title("Vérification moyenne, 532ppd")
plt.plot(altitude,data532_dp[1])
plt.plot(altitude,data532dpmoy)
plt.show()

#%%read fonction geom
ll=["","2","3"]
geom_som=[0]*len(altitude)
geom_moy=[]
for i in range(len(ll)):
    trajet_geom='/net/nfs/home/rey/Documents/code_Zoe/OHP/fct_geom'+ll[i]+'_OHP.txt'
    geom=np.loadtxt(trajet_geom)
    for i in range(len(altitude)):
        geom_som[i]=geom_som[i]+geom[i]
for i in range(len(altitude)):
    geom_moy.append(geom_som[i]/len(ll))

plt.figure()
plt.plot(altitude[:200],geom_moy[:200])
plt.xlabel("altitude (km)")
plt.ylabel("fonction recouvrement 532//")
plt.legend()
plt.grid()
plt.show()

geom532_fit_file=geom_moy


#%% calcul profil moleculaire et rayleigh
source_mol='ISA'
rayleigh532=juti.calc_molecular_profile(source_mol,altitude,532)

#%% calcul constante
idx1=np.where(altitude==6.)
id1=idx1[0][0]
idx2=np.where(altitude==8.01)
id2=idx2[0][0]

idx3=np.where(altitude==2.01)
id3=idx3[0][0]
idx4=np.where(altitude==3)
id4=idx4[0][0]


cste532=outils_OHP.calc_cste(data532moy,rayleigh532,id1,id2)
cste532_dp=outils_OHP.calc_cste(data532dpmoy,rayleigh532,id3,id4)

plt.figure()
plt.title('Vérification superpostition courbes, date+heure, 532//nm')
plt.plot(data532moy[:700],altitude[:700],linewidth=1.,label="PR2")
plt.plot(cste532*rayleigh532['beta_m'][:700]*rayleigh532['tm'][:700],altitude[:700],linestyle='--',color='k',label="Rayleigh*constante")
plt.ylabel("altitude (km)")
plt.grid()
plt.xlabel("voie1")
#plt.ylim(0,10)
#plt.xlim(-2.5,2.5)
plt.legend()
plt.show()

plt.figure()
plt.title('Vérification superpostition courbes, date+heure, 532//nm')
plt.plot(data532moy[:],altitude[:],linewidth=1.,label="PR2")
plt.plot(cste532*rayleigh532['beta_m'][:]*rayleigh532['tm'][:],altitude[:],linestyle='--',color='k',label="Rayleigh*constante")
plt.ylabel("altitude (km)")
plt.grid()
plt.xlabel("voie1")
#plt.ylim(0,10)
#plt.xlim(-2.5,2.5)
plt.legend()
plt.show()

"""
plt.figure()
plt.title('Vérification superpostition courbes, date+heure, 532ppd nm')
plt.plot(data532dpmoy[:300],altitude[:300],linewidth=1.,label="PR2")
plt.plot(cste532_dp*rayleigh532['beta_m'][:300]*rayleigh532['tm'][:300],altitude[:300],linestyle='--',color='k',label="Rayleigh*constante")
plt.ylabel("altitude (km)")
plt.xlabel("voie1")
#plt.ylim(0,10)
#plt.xlim(-2.5,2.5)
plt.legend()
plt.show()
"""

#%%Calcul fonction de recouvrement calculée avec la voie 532dp pour le 20170811 19h

"""
geom532=outils_OHP.calc_geom(data532dpmoy,cste532_dp,rayleigh532,altitude)
geom532_fit=outils_OHP.calc_geom_fit_moy(geom532,altitude)

plt.figure()
plt.title("Fonction de recouvrement, 532dpnm")
plt.xlabel("altitude (km)")
plt.ylabel("fonction de recouvrement")
plt.plot(altitude[:200],geom532[:200],label="532dp")
plt.plot(altitude[:200],geom532_fit_file[:200],label="532dp fit_file")
plt.plot(altitude[:200],geom532_fit[:200],label="532dp fit")
plt.grid()
plt.legend()
plt.show()


plt.figure()
plt.title("Fonction de recouvrement, 532dpnm")
plt.xlabel("altitude (km)")
plt.ylabel("fonction de recouvrement")
plt.plot(altitude[:600],geom532_fit[:600],label="532dp fit")
plt.plot(altitude[:600],geom532_fit_file[:600],label="532dp fit_file")
plt.grid()
plt.legend()
plt.show()


fichier = open("fct_geomtest_OHP.txt", "a")
for i in range(len(altitude)):
    fichier.write(str(geom532_fit[i])+"\n")
fichier.close()
"""
geom532_fit=geom532_fit_file


#%%Calcul data corrigé de constante

data532cor=[]
data532dpcor=[]
for i in range(len(altitude)):
    data532cor.append(data532moy[i]/cste532/geom532_fit[i])
    data532dpcor.append(data532dpmoy[i]/cste532_dp/geom532_fit[i])

plt.figure()
plt.title('Vérification data corrigé, date+heure, 532//nm')
plt.plot(data532cor[20:600],altitude[20:600],linewidth=1.,label="PR2")
plt.plot(rayleigh532['beta_m'][20:600]*rayleigh532['tm'][20:600],altitude[20:600],linestyle='--',color='k',label="Rayleigh*constante")
plt.ylabel("altitude (km)")
plt.grid()
plt.xlabel("voie1")
#plt.ylim(0,10)
#plt.xlim(-2.5,2.5)
plt.legend()
plt.show()
 
#%%Estimation transmission aérosols à zmax

LRp=50.
LRm=(8*np.pi)/3.
T_p532=outils_OHP.calc_transmission_aer(LRp,altitude,data532cor,rayleigh532)

#T_p532=outils_OHP.calc_transmission_aer2(LRp,altitude,data532cor,rayleigh532)

AOT=0.06
T_p532=np.exp(-1*AOT)

#%%Nouvelle estimation de la constante instrumentale

cste532=outils_OHP.calc_cste2(data532moy,rayleigh532,id1,id2,T_p532)
#cste532_dp=outils_OHP.calc_cste2(data532dpmoy,rayleigh532,id3,id4,1)

#%%Nouveau Calcul data corrigé de constante

data532cor=[]
data532dpcor=[]
for i in range(len(altitude)):
    data532cor.append(data532moy[i]/cste532/geom532_fit[i])
    data532dpcor.append(data532dpmoy[i]/cste532_dp/geom532_fit[i])

plt.figure()
plt.title('Vérification superpostition courbes, date+heure, 532//nm')
plt.plot(data532cor[10:600],altitude[10:600],linewidth=1.,label="PR2")
plt.plot(rayleigh532['beta_m'][10:600]*rayleigh532['tm'][10:600],altitude[10:600],linestyle='--',color='k',label="Rayleigh*constante")
plt.ylabel("altitude (km)")
plt.grid()
plt.xlabel("voie1")
plt.legend()
plt.show()




plt.figure()
plt.title('Vérification superpostition courbes, date+heure, 532dp nm')
plt.plot(data532dpcor[:600],altitude[:600],linewidth=1.,label="PR2")
plt.plot(rayleigh532['beta_m'][:600]*rayleigh532['tm'][:600],altitude[:600],linestyle='--',color='k',label="Rayleigh*constante")
plt.ylabel("altitude (km)")
plt.grid()
plt.xlabel("voie1")
#plt.ylim(0,10)
#plt.xlim(-2.5,2.5)
plt.legend()
plt.show()
#%% plot


plt.figure()
plt.title('Vérification superpostition courbes, date+heure, 532nm')
plt.plot(data532moy,altitude,linewidth=1.,label="PR2")
plt.plot(cste532[2]*rayleigh532['beta_m']*rayleigh532['tm']*T_p532*T_p532,altitude,linestyle='--',color='k',label="Rayleigh*constante")
#plt.yscale('Log')
plt.ylim(0,30)
#plt.ylim(-100,100)
plt.ylabel("altitude (km)")
plt.xlabel("log(voie1)")
plt.legend()
plt.show()


plt.figure()
plt.title("Fonction de recouvrement, 532nm")
plt.xlabel("altitude (km)")
plt.ylabel("fonction de recouvrement")
plt.plot(altitude,geom532_fit)
plt.xlim(0,10)
plt.legend()
plt.show()

plt.figure()
plt.title('Correction avec K, 2019 03 08 8h 808nm')
plt.plot(data808cor[id0:id9],altitude[id0:id9],label="PR2 corrigé de 0(z) et K")
plt.plot(rayleigh808['beta_m'][id0:id9]*rayleigh808['tm'][id0:id9]*T_p808*T_p808,altitude[id0:id9],label="Rayleigh")
plt.xlabel("PR2")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()

plt.figure()
plt.title('Correction avec K, 2019 03 08 8h 532nm')
plt.plot(data532cor[id0:id9],altitude[id0:id9],label="PR2 corrigé de 0(z) et K")
plt.plot(rayleigh532['beta_m'][id0:id9]*rayleigh532['tm'][id0:id9]*T_p532*T_p532,altitude[id0:id9],label="Rayleigh")
plt.xlabel("PR2")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()

print("Constante 808nm égale à", statistics.mean(cste808_liste),"à plus ou moins ",np.std(cste808_liste))
print("Constante 532nm égale à", statistics.mean(cste532_liste),"à plus ou moins ",np.std(cste532_liste))
