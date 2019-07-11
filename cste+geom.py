#Version Zoe
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 16:03:36 2019

Calcul du parametre geometrique et de la constante systeme pour le lidar Milan.

"""
#import statistics
import jutilities as juti
import matplotlib.pyplot as plt
import numpy as np
import outils
#from scipy.optimize import curve_fit


#%% Read file

#annee=["2018","2018","2018","2018","2018"]
annee=["2019","2019","2019","2019","2019"]
dates=["20190216","20190216","20190216","20190216","20190216"]
hour=["010000","011002","012004","013005","014007"]


trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/ascii/'+annee[2]+'/'+dates[2]+'/mli2MOY_'+dates[2]+'.'+hour[2]+'.JUS'
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
    for n in range(len(dates)):
        data808moy[i]=data808moy[i]+data808[n][i]
        data532moy[i]=data532moy[i]+data532[n][i]
        data808moy_dp[i]=data808moy_dp[i]+data808_dp[n][i]
        data532moy_dp[i]=data532moy_dp[i]+data532_dp[n][i]
    data808moy[i]=data808moy[i]/(len(dates))
    data532moy[i]=data532moy[i]/(len(dates))
    data808moy_dp[i]=data808moy_dp[i]/(len(dates))
    data532moy_dp[i]=data532moy_dp[i]/(len(dates))



plt.figure()
plt.title('Verif data moy')
plt.plot(data532_dp[3],altitude)
plt.plot(data532moy_dp,altitude)
plt.show()


plt.figure()
plt.title(dates[2]+', moyenne de 8h à 8h40, 532nm')
plt.plot(data532moy[:500],altitude[:500])
#plt.xscale('log')
plt.ylabel("altitude (km)")
plt.xlabel("intensité du signal rétrodiffusé (UA) ")
plt.legend()
plt.show()


#%% calcul profil moleculaire et rayleigh
source_mol='ISA'
rayleigh808=juti.calc_molecular_profile(source_mol,altitude,808)
rayleigh532=juti.calc_molecular_profile(source_mol,altitude,532)




plt.figure()
plt.title(dates[2]+', moyenne de 1h à 1h40, 808nm')
plt.plot(data808moy[:300],altitude[:300],label="PR2")
plt.plot(rayleigh808['beta_m'][:300]*rayleigh808['tm'][:300],altitude[:300],label="Profil Rayleigh")
plt.ylabel("altitude (km)")
plt.xscale('log')
plt.xlabel("intensité du signal rétrodiffusé (UA) ")
plt.legend()
plt.show()




#%% limite de bruit étendue
z0=0.08

data1=[]
data2=[]
data3=[]
data4=[]
for i in range(len(altitude)):
    data1.append(data532moy[i]/(altitude[i]-z0)**2)
    data2.append(data532moy_dp[i]/(altitude[i]-z0)**2)
    data3.append(data808moy[i]/(altitude[i]-z0)**2)
    data4.append(data808moy_dp[i]/(altitude[i]-z0)**2)

bruit532=[]
bruit532_dp=[]
bruit808=[]
bruit808_dp=[]

ectype1=np.std(data1[1300:])
ectype2=np.std(data2[1300:])
ectype3=np.std(data3[1300:])
ectype4=np.std(data4[1300:])


for i in range(len(altitude)):
    bruit532.append(2.*ectype1*(altitude[i]-z0)*(altitude[i]-z0))
    bruit532_dp.append(2.*ectype2*(altitude[i]-z0)*(altitude[i]-z0))
    bruit808.append(2.*ectype3*(altitude[i]-z0)*(altitude[i]-z0))
    bruit808_dp.append(2.*ectype4*(altitude[i]-z0)*(altitude[i]-z0))

#%% calcul constante

idx_8000m=np.where(altitude==6)
id8=idx_8000m[0][0]
idx_9000m=np.where(altitude==7.005)
id9=idx_9000m[0][0]

idx_3000m=np.where(altitude==2.01)
id3=idx_3000m[0][0]
idx_4000m=np.where(altitude==3)
id4=idx_4000m[0][0]



print("808")
cste808=outils.calc_cste(data808,data808moy,rayleigh808,id8,id9)
print("532")
cste532=outils.calc_cste(data532,data532moy,rayleigh532,id8,id9)
print("808_dp")
cste808_dp=outils.calc_cste(data808,data808moy_dp,rayleigh808,id3,id4)
print("532_dp")
cste532_dp=outils.calc_cste(data532,data532moy_dp,rayleigh532,id3,id4)

plt.figure()
plt.title('Vérification superpostition courbes, '+dates[2]+' '+hour[2]+', 532 dp')
plt.plot(data532moy_dp,altitude,linewidth=1.,label="PR2")
plt.plot(cste532_dp*rayleigh532['beta_m']*rayleigh532['tm'],altitude,linestyle='--',color='k',label="Rayleigh*constante")
plt.plot(bruit532_dp,altitude,label="limite de bruit")
#plt.xscale('Log')
plt.ylim(0,10)
plt.xlim(-5,5)
plt.ylabel("altitude (km)")
plt.xlabel("voie 2")
plt.grid()
plt.legend()
plt.show()

plt.figure()
plt.title('Vérification superpostition courbes, '+dates[2]+' '+hour[2]+', 532')
plt.plot(data532moy,altitude,linewidth=1.,label="PR2")
plt.plot(cste532*rayleigh532['beta_m']*rayleigh532['tm'],altitude,linestyle='--',color='k',label="Rayleigh*constante")
#plt.plot(bruit,altitude)
#plt.yscale('Log')
#plt.plot(bruit532_dp,altitude,label="limite de bruit")
plt.ylim(0,10)
plt.xlim(-5,75)
plt.ylabel("altitude (km)")
plt.xlabel("voie 1")
plt.grid()
plt.legend()
plt.show()

plt.figure()
plt.title('Vérification superpostition courbes, '+dates[2]+' '+hour[2]+', 808 dp')
plt.plot(data808moy_dp,altitude,linewidth=1.,label="PR2")
plt.plot(cste808_dp*rayleigh808['beta_m']*rayleigh808['tm'],altitude,linestyle='--',color='k',label="Rayleigh*constante")
#plt.yscale('Log')
plt.plot(bruit532_dp,altitude,label="limite de bruit")
plt.ylim(0,10)
plt.xlim(-5,5)
plt.ylabel("altitude (km)")
plt.xlabel("voie 4")
plt.grid()
plt.legend()
plt.show()

plt.figure()
plt.title(dates[2]+' '+hour[2]+', 808')
plt.plot(data808moy,altitude,linewidth=1.,label="PR2")
plt.plot(cste808*rayleigh808['beta_m']*rayleigh808['tm'],altitude,linestyle='--',color='k',label="Rayleigh*constante")
plt.plot(rayleigh808['beta_m']*rayleigh808['tm'],altitude,linestyle='--',color='g',label="Rayleigh")
#plt.yscale('Log')
#plt.plot(bruit532_dp,altitude,label="limite de bruit")
plt.ylim(0,6)
plt.xscale('log')
plt.ylabel("altitude (km)")
plt.xlabel("voie 3")
#plt.grid()
plt.legend()
plt.show()

#%%Calcul fonction de recouvrement

trajet_geom='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/fgeom/milan_20190121_fgeom_v6.txt'
file=np.loadtxt(trajet_geom,skiprows=13)
fgeom532=file[:,1]
fgeom808=file[:,3]

geom532_fit=list(fgeom532[::-1])
geom808_fit=list(fgeom808[::-1])

for i in range(len(altitude)-len(geom532_fit)):
    geom532_fit.append(1)
    geom808_fit.append(1)

"""
geom808=outils.calc_geom(data808,cste808,rayleigh808,altitude)
geom532=outils.calc_geom(data532,cste532,rayleigh532,altitude)

geom808_fit=outils.calc_geom_fit_moy(geom808,altitude,id4)
geom532_fit=outils.calc_geom_fit_moy(geom532,altitude,id4)
"""
"""
fichier = open("fct_geom532_JUS.txt", "a")
for i in range(len(altitude)):
    fichier.write(str(geom532_fit[i])+"\n")
fichier.close()

fichier = open("fct_geom808_JUS.txt", "a")
for i in range(len(altitude)):
    fichier.write(str(geom808_fit[i])+"\n")
fichier.close()
"""

plt.figure()
plt.ylabel("altitude (km)")
plt.xlabel("fonction de géométrie, 808//nm")
plt.plot(geom808_fit,altitude)
plt.plot(geom532_fit,altitude)
plt.ylim(0,2)
plt.legend()
plt.show()


                        
#%%Calcul data corrigé de constante

data808cor=[]
data532cor=[]
data532cor_dp=[]
data808cor_dp=[]
for i in range(len(altitude)):
    data808cor.append(data808moy[i]/cste808)
    data808cor_dp.append(data808moy_dp[i]/cste808_dp)
    data532cor.append(data532moy[i]/cste532)
    data532cor_dp.append(data532moy_dp[i]/cste532_dp)
    


plt.figure()
plt.title('Vérification superpostition courbes, date+heure, 532')
plt.plot(data532cor[:300],altitude[:300],linewidth=1.,label="PR2")
plt.plot(rayleigh532['beta_m'][:300]*rayleigh532['tm'][:300],altitude[:300],linestyle='--',color='k',label="Rayleigh*constante")
plt.ylabel("altitude (km)")
plt.xlabel("voie1")
plt.grid()
plt.legend()
plt.show()

plt.figure()
plt.title('Vérification superpostition courbes, date+heure, 808')
plt.plot(data808cor[:300],altitude[:300],linewidth=1.,label="PR2")
plt.plot(rayleigh808['beta_m'][:300]*rayleigh808['tm'][:300],altitude[:300],linestyle='--',color='k',label="Rayleigh*constante")
plt.ylabel("altitude (km)")
plt.xlabel("vioe3")
plt.grid()
plt.legend()
plt.show()



 
#%%Estimation transmission aérosols à zmax

LRp=50.
LRm=(8*np.pi)/3.

print("808")
T_p808,imax808=outils.calc_transmission_aer(LRp,LRm,altitude,data808cor,rayleigh808)
print("532")
T_p532,imax532=outils.calc_transmission_aer(LRp,LRm,altitude,data532cor,rayleigh532)
print("808_dp")
T_p808_dp,imax808=outils.calc_transmission_aer(LRp,LRm,altitude,data808cor_dp,rayleigh808)
print("532_dp")
T_p532_dp,imax532=outils.calc_transmission_aer(LRp,LRm,altitude,data532cor_dp,rayleigh532)

#%%Calcul couleur ratio à zmax

betap808=(data808cor[imax808]/(rayleigh808['tm'][imax808]*T_p808)-rayleigh808['beta_m'][imax808])
betap532=data532cor[imax532]/(rayleigh532['tm'][imax532]*T_p532)-rayleigh532['beta_m'][imax532]

print("color ratio=",betap808/betap532)

#%%Nouvelle estimation de la constante instrumentale

print("808")
cste808=outils.calc_cste2(data808,data808moy,rayleigh808,id8,id9,T_p808)
print("532")
cste532=outils.calc_cste2(data532,data532moy,rayleigh532,id8,id9,T_p532)
print("808_dp")
cste808_dp=outils.calc_cste2(data808,data808moy_dp,rayleigh808,id3,id4,T_p808_dp)
print("532_dp")
cste532_dp=outils.calc_cste2(data532,data532moy_dp,rayleigh532,id3,id4,T_p532_dp)

#%%Zone aveuglement

i_avg=5

altitude1=altitude
altitude=altitude[i_avg:]
data808moy1=data808moy
data532moy=data532moy[i_avg:]
data808moy=data808moy[i_avg:]
data532moy_dp=data532moy_dp[i_avg:]
data808moy_dp=data808moy_dp[i_avg:]

rayleigh532['beta_m']=rayleigh532['beta_m'][i_avg:]
rayleigh532['tm']=rayleigh532['tm'][i_avg:]
rayleigh808['beta_m']=rayleigh808['beta_m'][i_avg:]
rayleigh808['tm']=rayleigh808['tm'][i_avg:]
geom532_fit=geom532_fit[i_avg:]
geom808_fit=geom808_fit[i_avg:]
                        

#%%Nouveau calcul de data corrigé

data808cor=[]
data532cor=[]
data532cor_dp=[]
data808cor_dp=[]
for i in range(len(altitude)):
    data808cor.append(data808moy[i]/cste808/geom808_fit[i])
    data808cor_dp.append(data808moy_dp[i]/cste808_dp/geom808_fit[i])
    data532cor.append(data532moy[i]/cste532/geom532_fit[i])
    data532cor_dp.append(data532moy_dp[i]/cste532_dp/geom532_fit[i])
   
plt.figure()
plt.title('Correction du signal, '+dates[2]+" "+hour[2]+", 532nm")
plt.plot(data532cor[:300],altitude[:300],linewidth=1.,label="PR2 corrigé")
plt.plot(rayleigh532['beta_m'][:300]*rayleigh532['tm'][:300],altitude[:300],linestyle='--',color='k',label="Rayleigh*constante")
plt.ylabel("altitude (km)")
plt.xlabel("voie1")
plt.grid()
plt.legend()
plt.show()
    
    
plt.figure()
plt.title('Vérification superpostition courbes, date+heure, 532dp')
plt.plot(data532cor_dp[:300],altitude[:300],linewidth=1.,label="PR2")
plt.plot(rayleigh532['beta_m'][:300]*rayleigh532['tm'][:300],altitude[:300],linestyle='--',color='k',label="Rayleigh*constante")
plt.ylabel("altitude (km)")
plt.xlabel("voie2")
plt.grid()
plt.legend()
plt.show()

plt.figure()
plt.title('Vérification superpostition courbes, date+heure, 808dp')
plt.plot(data808cor_dp[:300],altitude[:300],linewidth=1.,label="PR2")
plt.plot(rayleigh808['beta_m'][:300]*rayleigh808['tm'][:300],altitude[:300],linestyle='--',color='k',label="Rayleigh*constante")
plt.ylabel("altitude (km)")
plt.xlabel("voie4")
plt.grid()
plt.legend()
plt.show()

plt.figure()
plt.title('Vérification superpostition courbes, date+heure, 808')
plt.plot(data808cor[:300],altitude[:300],linewidth=1.,label="PR2")
plt.plot(rayleigh808['beta_m'][:300]*rayleigh808['tm'][:300],altitude[:300],linestyle='--',color='k',label="Rayleigh*constante")
plt.ylabel("altitude (km)")
plt.xlabel("voie4")
plt.grid()
plt.legend()
plt.show()


#%%plot avant après

plt.figure()
plt.title(dates[2]+" "+hour[2]+", 808nm")
plt.subplot(1,2,1)

plt.subplot(121)
plt.title('PR2')
plt.plot(data808moy1[:68],altitude1[:68],linewidth=1.,label="PR2")
plt.ylim(0,1)
plt.xscale('log')
plt.ylabel("altitude (km)")
plt.xlabel("voie 808// (UA)")


plt.subplot(122)
plt.title('PR2 corrigé')
plt.plot(data808cor[:68],altitude[:68],linewidth=1.)
plt.xlabel("voie 808// (km⁻¹ sr⁻¹)")
#plt.grid()
plt.ylim(0,1)
plt.xscale('log')
plt.legend()
plt.show()

