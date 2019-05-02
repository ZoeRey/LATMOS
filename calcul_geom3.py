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
import outils
from scipy.optimize import curve_fit


#%% Read file

#les deux premiers fichiers sont utilisés pour calculer une moyenne de la fonction de recouvremenent (couche limite très peu développée)
#les autres sont utilisés pour le calcul de la constante instrumentale
#annee=["2019","2018","2018","2018","2018","2018","2018"]
annee=["2019","2018","2019","2019","2019","2019","2019"]
dates=["20190215","20181226","20190121","20190121","20190121","20190121"]
hour=["000000","020000","010000","011002","012004","013007","014009"]


trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/ascii/'+annee[0]+'/'+dates[0]+'/mli2MOY_'+dates[0]+'.'+hour[0]+'.JUS'
file=np.loadtxt(trajet,skiprows=36)
altitude=file[:,0]
altitude=altitude[::-1]   

data808=[[0]*len(altitude) for i in range(len(dates))]
data532=[[0]*len(altitude) for i in range(len(dates))]
data808_dp=[[0]*len(altitude) for i in range(len(dates))]
data532_dp=[[0]*len(altitude) for i in range(len(dates))]

for i in range(len(dates)):
    trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/ascii/'+annee[i]+'/'+dates[i]+'/mli2MOY_'+dates[i]+'.'+hour[i]+'.JUS'
    file=np.loadtxt(trajet,skiprows=36) 
    voie808=file[:,3]
    voie532=file[:,1]
    voie808_dp=file[:,4]
    voie532_dp=file[:,2]
    data808[i]=voie808[::-1]
    data532[i]=voie532[::-1]
    data808_dp[i]=voie808_dp[::-1]
    data532_dp[i]=voie532_dp[::-1]

data532moy=[0]*len(altitude)
data808moy=[0]*len(altitude)
data532moy_dp=[0]*len(altitude)
data808moy_dp=[0]*len(altitude)

for i in range(len(altitude)):
    for n in range(2,len(dates)):
        data808moy[i]=data808moy[i]+data808[n][i]
        data532moy[i]=data532moy[i]+data532[n][i]
        data808moy_dp[i]=data808moy_dp[i]+data808_dp[n][i]
        data532moy_dp[i]=data532moy_dp[i]+data532_dp[n][i]
    data808moy[i]=data808moy[i]/(len(dates)-2)
    data532moy[i]=data532moy[i]/(len(dates)-2)
    data808moy_dp[i]=data808moy_dp[i]/(len(dates)-2)
    data532moy_dp[i]=data532moy_dp[i]/(len(dates)-2)


"""
plt.figure()
plt.title('Verif data + data_dp')
plt.plot(altitude,data808tot)
plt.plot(altitude,data808moy_dp)
plt.show()

plt.plot(altitude,data808moy)
plt.plot(altitude,data808moy_dp)

plt.figure()
plt.title("Vérification moyenne")
plt.plot(altitude,data808_dp[3])
plt.plot(altitude,data808moy_dp)
plt.xlim(0,14)
plt.ylim(-100,150)
plt.show()
"""

#%% calcul profil moleculaire et rayleigh
source_mol='ISA'
rayleigh808=juti.calc_molecular_profile(source_mol,altitude,808)
rayleigh532=juti.calc_molecular_profile(source_mol,altitude,532)

"""
plt.plot(altitude,rayleigh532['beta_m']*rayleigh532['tm'], label='Rayleigh 532')
plt.plot(altitude,rayleigh808['beta_m']*rayleigh808['tm'],label='Rayleigh 808')
"""
#%% calcul constante
idx_3000m=np.where(altitude==0.6)
id3=idx_3000m[0][0]
idx_4000m=np.where(altitude==1.005)
id4=idx_4000m[0][0]
idx_8000m=np.where(altitude==6)
id8=idx_8000m[0][0]
idx_9000m=np.where(altitude==8.01)
id9=idx_9000m[0][0]
idx_0m=np.where(altitude==0.)
id0=idx_0m[0][0]
idx_5000m=np.where(altitude==5.01)
id5=idx_5000m[0][0]


cste808=outils.calc_cste(data808,data808moy,rayleigh808,id3,id4,id8,id9)
cste532=outils.calc_cste(data532,data532moy,rayleigh532,id3,id4,id8,id9)


#%%Calcul fonction de recouvrement

geom808=outils.calc_geom(data808,cste808,rayleigh808,altitude)
geom532=outils.calc_geom(data532,cste532,rayleigh532,altitude)

geom808_fit=outils.calc_geom_fit_moy(geom808,altitude,id0,id4)
geom532_fit=outils.calc_geom_fit_moy(geom808,altitude,id0,id4)


#%%Calcul data corrigé de géométrie et constante

data808cor=[]
data532cor=[]
for i in range(len(altitude)):
    data808cor.append(data808moy[i]/cste808[2]/geom808_fit[i])
    data532cor.append(data532moy[i]/cste532[2]/geom532_fit[i])
 
#%%Estimation transmission aérosols à zmax

LRp=50.
LRm=(8*np.pi)/3.

T_p808=outils.calc_transmission_aer(LRp,LRm,altitude,data808moy,data808cor,rayleigh808)
T_p532=outils.calc_transmission_aer(LRp,LRm,altitude,data532moy,data532cor,rayleigh532)


#%%Nouvelle estimation de la constante instrumentale

cste808=outils.calc_cste2(data808,data808moy,rayleigh808,id3,id4,id8,id9,T_p808)
cste532=outils.calc_cste2(data532,data532moy,rayleigh532,id3,id4,id8,id9,T_p532)

#%%Nouveau calcul de data corrigé
data808cor=[]
data532cor=[]
for i in range(len(altitude)):
    data808cor.append(data808moy[i]/cste808[2]/geom808_fit[i])
    data532cor.append(data532moy[i]/cste532[2]/geom532_fit[i])

#%%Calcul coefficient rétrodiffusion particule tot, méthode de Fernald, avec des intégrales de 0 à z, trop facile pour être juste

beta808_aer=outils.calc_beta_aer(rayleigh808,altitude,data808moy,LRp,LRm,cste808[2])
beta532_aer=outils.calc_beta_aer(rayleigh532,altitude,data532moy,LRp,LRm,cste532[2])

"""
plt.figure()
plt.title("coefficient retrodiffusion, methode Fernald 808nm")
plt.plot(altitude,beta808_aer,label="beta aerosol")
plt.plot(altitude,rayleigh808['beta_m'],label="beta moléculaire")
plt.xlabel("altitude (km)")
plt.xlim(0,30)
plt.legend()
plt.show()
"""

#%%Calcul coefficient rétrodiffusion particule, Fernald,plus stable
"""
izrkm=np.where(altitude==8.01)
izr=izrkm[0][0]
data808moy_zr=statistics.mean(data808moy[izr-10:izr+10])


beta808_aer2=outils.calc_beta_aer_stable(izr,data808moy_zr,altitude,rayleigh808,data808moy,LRp,LRm)


plt.figure()
plt.title("coefficient retrodiffusion aérosols, 808nm")
plt.plot(altitude,beta808_aer2/cste808[2],label="beta aerosol")
plt.plot(altitude,rayleigh808['beta_m']/cste808[2],label="beta moléculaire")
plt.legend()
plt.show()    

beta808_p_fit=beta808_p[:id8]
coeff_beta808_p_fit=np.polyfit(altitude[id8:],beta808_p[id8:],2)
for i in range(len(altitude[id8:])):
    beta808_p_fit.append(coeff_beta808_p_fit[0]*altitude[i]*altitude[i]+coeff_beta808_p_fit[1]*altitude[i]+coeff_beta808_p_fit[2])

plt.figure()
plt.title("coefficient retrodiffusion aérosols, 808nm")
plt.plot(altitude,beta808_p,label="beta aerosol")
plt.plot(altitude,beta808_p_fit,label="beta aérosol fit")
plt.legend()
plt.show()    
"""

 #%%Moyenne cste
cste808_liste=[36689419.37807519,33214615.779274035,33606423.776109919,31568627.592203494,33948537.277566925,41044935.117830724,33606423.776109919,34789707.198917709]
cste532_liste=[13716012.962917635,14174571.356842197,14607243.276845891,15351104.91029889,12607659.628663693,12345142.824067345,14607243.276845891,12332518.529443894]


#%% plot

plt.figure()
plt.title('Vérification superpostition courbes, date+heure, 808nm')
plt.plot(altitude,data808moy,linewidth=1.,label="PR2")
plt.plot(altitude,cste808[2]*rayleigh808['beta_m']*rayleigh808['tm']*T_p808*T_p808,linestyle='--',color='k',label="Rayleigh*constante")
#plt.yscale('Log')
plt.xlim(0,30)
#plt.ylim(-100,150)
plt.xlabel("altitude (km)")
plt.ylabel("log(voie1)")
plt.legend()
plt.show()

plt.figure()
plt.title('Vérification superpostition courbes, date+heure, 532nm')
plt.plot(altitude,data532moy,linewidth=1.,label="PR2")
plt.plot(altitude,cste532[2]*rayleigh532['beta_m']*rayleigh532['tm']*T_p532*T_p532,linestyle='--',color='k',label="Rayleigh*constante")
#plt.yscale('Log')
plt.xlim(0,30)
#plt.ylim(-100,100)
plt.xlabel("altitude (km)")
plt.ylabel("log(voie1)")
plt.legend()
plt.show()


plt.figure()
plt.title("Fonction de recouvrement, 808nm")
plt.xlabel("altitude (km)")
plt.ylabel("fonction de recouvrement")
plt.plot(altitude,geom808_fit)
plt.xlim(0,2)
plt.legend()
plt.show()

plt.figure()
plt.title("Fonction de recouvrement, 532nm")
plt.xlabel("altitude (km)")
plt.ylabel("fonction de recouvrement")
plt.plot(altitude,geom532_fit)
plt.xlim(0,2)
plt.legend()
plt.show()

plt.figure()
plt.title('Correction avec K, 2019 03 28 17h 808nm')
plt.plot(altitude[id0:id9],data808cor[id0:id9],label="PR2 corrigé de 0(z) et K")
plt.plot(altitude[id0:id9],rayleigh808['beta_m'][id0:id9]*rayleigh808['tm'][id0:id9]*T_p808*T_p808,label="Rayleigh")
plt.ylabel("PR2")
plt.xlabel("altitude (km)")
plt.legend()
plt.show()

plt.figure()
plt.title('Correction avec K, 2019 03 28 17h 532nm')
plt.plot(altitude[id0:id9],data532cor[id0:id9],label="PR2 corrigé de 0(z) et K")
plt.plot(altitude[id0:id9],rayleigh532['beta_m'][id0:id9]*rayleigh532['tm'][id0:id9]*T_p532*T_p532,label="Rayleigh")
plt.ylabel("PR2")
plt.xlabel("altitude (km)")
plt.legend()
plt.show()

print("Constante 808nm égale à", statistics.mean(cste808_liste),"à plus ou moins ",np.std(cste808_liste))
print("Constante 532nm égale à", statistics.mean(cste532_liste),"à plus ou moins ",np.std(cste532_liste))
