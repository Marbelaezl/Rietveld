# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 21:24:02 2022

@author: migue
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import Wyckoff
import tkinter as tk
from tkinter import filedialog
import tablas
import time
from numpy.linalg import norm, solve
import matplotlib.animation as animation
#mpl.use('TkAgg')
#IMPORTANTE PARA EL EJECUTABLE AAAAAAAAAAAAAAA :)


minang=10
maxang=90

np.set_printoptions(suppress=True)
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 100
#Condiciones iniciales
#Longitudes en A, ángulos en radianes
#Estructura actual: Ti hexagonal
k=1
a=2.9500
b=2.9500
c=4.6800
al=90*np.pi/180
be=90*np.pi/180
ga=120*np.pi/180
la=1.54060
la2=1.544398
delta=0.00001
cel1=np.array([
    #[22,0,0,0.25,0,1],
    [22,1/3,2/3,0.25,0,1]
    ])
ks1=[k,a,b,c,al,be,ga]
#TiO2
ks2=[100,4.593,4.593,2.9590,90.*np.pi/180,90*np.pi/180,90*np.pi/180]
#ks2=[100,100,50,40,90*np.pi/180,90*np.pi/180,90*np.pi/180]



n=1
cel=np.array([
      [22,0,0,0,0,1],
      [8,0.306,0.306,0,0,n],
      ])
cel2=Wyckoff.Wyckoff(136, cel)
ks3=[100,3.7850 ,3.7850 ,9.520,90*np.pi/180,90*np.pi/180,90*np.pi/180]
cel3=Wyckoff.Wyckoff(141, np.array([[22,0,0-1/4,0-3/8,0,1],[8,0,0-1/4,0.20806-3/8,0,1]]))
testvecs=np.array([
    [0.33898305,0,0],
    [0,0.33898305,0],
    [0,0,0.2177226214]
    ])

filename= "C:\\Users\migue\Downloads\HF 2.csv"
arr=np.loadtxt(filename,skiprows=27,delimiter=",")

data=np.genfromtxt("HF2.csv",delimiter=",")
ba2sm=np.genfromtxt("Ba2SbSmO6.ASC",delimiter=" ")
#------------------------------------------------------------------------------------------------------------
#Regresión polinómica para el fondo.
class datos():
    
    def PowBackground(self):
        "Se espera una entrada de 2 x n; en el que la primera columna son ángulos y la segunda son intensidades. Determina el fondo"
        #Primero, se escribe la matriz de datos
        matriz=np.vstack((
            self.angs**2,
            self.angs,
            np.ones_like(self.angs),
            1/(self.angs),
            1/(self.angs**2),
            1/(self.angs**3),
            1/(self.angs**4),
            1/(self.angs**5),
            1/(self.angs**6))).T
        #Vector de intensidades observadas
        y=np.array([self.ints]).T
        #Fórmula de regresión polinómica: (X^T X)^-1 X^T
        MatFin= np.linalg.inv((matriz.T @ matriz)) @ matriz.T
        #La estimación de los coeficientes es MatFin @ Intensidades.
        #La interpretación de los coeficientes es tal que el fondo es fondo[0] + fondo[1]/theta + fondo[2]/theta^2 + ...
        return MatFin @ y
    def __init__(self,entrada):
        self.angs=entrada[:,0]
        self.ints=entrada[:,1]
        self.bgcoeff= self.PowBackground()
        self.bgints = np.zeros_like(self.angs)
        for i in range(0,len(self.bgcoeff)):
            self.bgints = self.bgints + (self.bgcoeff[i]/(self.angs**(i-2)))
        sd= np.std(self.ints-self.bgints)
        
        #Se guardan los puntos poco intensos en datos_ruido y los picos en intensos. 
        #Al usar np.where, se dan como listas de indices, no como los puntos en sí
        datos_ruido=np.where(self.ints-self.bgints <= sd)
        intensos=np.where(np.abs(self.ints-self.bgints > sd))
        
        #Se vuelve a calcular el fondo después de eliminar los datos intensos (picos)
        self.angs=entrada[datos_ruido,0].flatten()
        self.ints=entrada[datos_ruido,1].flatten()
        self.bgcoeff= self.PowBackground()
        #Habiendo hallado los coeficientes, se restituye la información comlpeta
        self.angs=entrada[:,0]
        self.ints=entrada[:,1]
        self.bgints = np.zeros_like(self.angs)
        for i in range(0,len(self.bgcoeff)):
            self.bgints = self.bgints + (self.bgcoeff[i]/(self.angs**(i-2)))

ASCBa =datos(ba2sm)
HF2 = datos(data)

fig,ax=plt.subplots()

plt.xlabel(r'$2\theta($°)')
plt.ylabel("Intensidad (counts)")
fig,ax=plt.subplots()
ax.plot(ASCBa.angs,ASCBa.ints-ASCBa.bgints)


# Funciones requeridas para el patrón:
    
def disf(vecs,hkl):
    "Retorna la distancia interplanar del plano con índices h,k,l. vecs es la matriz de vectores EN ESPACIO RECÍPROCO"
    h=hkl[0]
    k=hkl[1]
    l=hkl[2]
    ac=vecs[0]
    bc=vecs[1]
    cc=vecs[2]
    try:
        return 1/(sum((h*ac+k*bc+l*cc)**2)**(1/2))
    except:
        return[-1]

def angf(vecs,hkl,la=la):
    h=hkl[0]
    k=hkl[1]
    l=hkl[2]
    if (la/(2*disf(vecs,[h,k,l])))<1:
        theta=np.arcsin(la/(2*disf(vecs,[h,k,l])))
        return ([360*theta/np.pi,h,k,l,1])
#Si el ángulo no existe (la codnición de Bragg da sin(theta)>1, la función devuelve [-1], negativo para que nunca sea igual a minang)
    return [-1]
    
def vecdisf(vecs,hkl):
    "disf vectorizada para tomar un array hkl más grande (bidimensional)"
    return (1/np.linalg.norm(hkl@vecs,axis=1))

def vecangf(vecs,hkl,la=la):

    prov=la/(2*vecdisf(vecs,hkl))
    indexes=np.where(prov<1)
    angs=(360*np.arcsin(prov[indexes]))/np.pi
    prov=np.hstack([np.array([angs]).T,hkl[indexes],np.array([np.ones_like(angs)]).T])
    indexes=np.where((prov[:,0]<maxang))
    return prov[indexes]
    
def fastpicoshelp(ind,flag=True):
    #Función ayudante para el indexado de picos. Calcula todas las combinaciones de enteros entre
    #0 e ind si flag=True, o entre -ind e ind si flag=False
    if flag:
        prov=np.arange(0,(ind[0]+1)*(ind[1]+1)*(ind[2]+1),1,dtype='i4')
        res=np.fix(np.vstack(np.array([prov/((ind[1]+1)*(ind[2]+1)),(prov/(ind[2]+1))%(ind[1]+1),prov%(ind[2]+1)]).T))
        
    else:
        prov=np.arange(0,(2*ind[0]+1)*(2*ind[1]+1)*(2*ind[2]+1),1,dtype='i4')
        res=np.floor(np.vstack(np.array([prov/((2*ind[1]+1)*(2*ind[2]+1)),(prov/(2*ind[2]+1))%(2*ind[1]+1),prov%(2*ind[2]+1)]).T))
        res[:]=res[:]-ind
    res=np.delete(res, (np.where(np.all(res==0,axis=-1))),axis=0)
    return res


def vecpicosf(vecs,indmax,la=la):
    logd=6 #logaritmo de delta
    delta=10**-logd
    prov=vecangf(vecs,fastpicoshelp(indmax),la)
    prov=prov[prov[:, 0].argsort()]
    prov = np.delete(prov, np.argwhere(np.ediff1d(prov[:,0].round(logd)) < delta/10) + 1,axis=0)
    prov2=vecangf(vecs,fastpicoshelp(indmax,flag=False),la=la)
    minii=np.searchsorted(prov[:,0],minang,side='left')
    maxii=np.searchsorted(prov[:,0],maxang,side='right')
    prov=prov[minii:maxii]
    prov2=prov2[prov2[:, 0].argsort()]
    minii=np.searchsorted(prov2[:,0],minang,side='left')
    maxii=np.searchsorted(prov2[:,0],maxang,side='right')
    prov2=prov2[minii:maxii]
    unique, counts = np.unique(np.round(prov2[:,0],logd), return_counts=True)
    prov[:,-1]+=counts-1
    return prov

Cro_Mann=tablas.Cro_Mann

def atom_strucf(n,par,ion=0):
    "Devuelve el factor de estructura del átomo con número n según el parámetro par = sin(theta)/lambda"
    "Solo considera la parte real; se omiten efectos de dispersión anómala"
    return (Cro_Mann[ion,n,0]*np.e**(-Cro_Mann[ion,n,4]*(par**2))+
          Cro_Mann[ion,n,1]*np.e**(-Cro_Mann[ion,n,5]*(par**2))+
          Cro_Mann[ion,n,2]*np.e**(-Cro_Mann[ion,n,6]*(par**2))+
          Cro_Mann[ion,n,3]*np.e**(-Cro_Mann[ion,n,7]*(par**2))+
          Cro_Mann[ion,n,8])

def matrizstrucf(multip,celda, la=la):

    m1=multip[:,1:4]
    m2=celda[:,1:4].T
    fases= m1@m2
    factores=np.ones((3,len(multip),len(celda)))
    factores[0]= (factores[0] * [celda[:,0]])
    factores[1]= (factores[1] * [celda[:,4]])
    factores[2]= (factores[2] * [celda[:,5]])
    par=np.ones((len(celda),len(multip)))
    par=par * np.sin(multip[:,0]*np.pi/360)/la
    struct=(atom_strucf(factores[0].astype(int), par.T,factores[1].astype(int)))
    res= np.exp((0+2j)*np.pi*fases)*struct
    res= res * factores[2]
    return(np.sum(res,axis=1))
    

def fastintensidadesf(multip,struct):
    #como la entrada a[i,0] es 2theta, se multiplica por pi/360 en lugar de pi/180
    ang=multip[:,0]*np.pi/360
    #Factor de Lorentz - Geometría Bragg-Brentano 
    lo= (1+ (np.cos(2*ang)**2))/(np.sin(2*ang)*np.sin(ang))
    multip[:,4]=np.abs(struct)**2 *lo*multip[:,4]
    multip[:,4]=multip[:,4]/max(multip[:,4])
    return multip
    
def ord2f(intens):
    "Ordena los picos según el ángulo 2theta"
    intens=np.delete(intens,np.where(intens[:,-1]<0.001),axis=0)
    intens=intens[np.argsort(intens[:,0])]
    return intens


def MinInfoPat(ctes=ks1,celda=cel1,ForceIndMax=9999, ReturnMultip=False, la=la):
    start=time.time()
    "Regresa la lista ordenada de picos e intensidades relativas, dadas las condiciones"
    #Primero, se determinan los vectores unitarios en el sistema de ejes perpendiculares con 
    #x en dirección a, y sobre el plano definido por a y b, y el componente de c en z positivo.
    #Esto se hace porque la distancia interplanar es 1/|ha*+kb*+lc*|, con a*,b* y c* vectores recíprocos
    avec=np.array([ctes[1],0,0])
    #La fórmula de b se puede hallar usando a . b = abcos(gamma)
    bvec=np.array([ctes[2]*np.cos(ctes[6]),ctes[2]*np.sin(ctes[6]),0])
    #Para hallar el componente en x de c (cprov[0]), se usa a . c = accos(beta), y lo mismo para el componente y
    cprov=np.array([ctes[3]*np.cos(be),ctes[3]*((np.cos(ctes[4])-((np.cos(ctes[6]))*np.cos(ctes[5])))/np.sin(ctes[6]))])
    #El componente z se halla usando cx y cy, y usando que la magnitud de cvec es c
    cvec=np.array([cprov[0],cprov[1],-(ctes[3]**2-cprov[0]**2-cprov[1]**2)**(1/2)])*-1
    #el orden en la celda unitaria es átomo, ion, x, y, z
    #V es el volumen de la celda unidad
    V=np.cross(bvec,cvec)@avec
    # ac es el complemento de a, bc el de b y cc el de c
    ac=np.cross(bvec,cvec)/V
    bc=np.cross(cvec,avec)/V
    cc=np.cross(avec,bvec)/V
    Vecs=[ac,bc,cc]
    #Matriz de vectores recíprocos
    #El índice máximo se obtiene con la máxima reflexión del tipo x 0 0. Como Lsin/2d = ha* + kb* + lc*,
    #entonces indmax = int(max(Lsin/2a*, Lsin/2b*, Lsin/2c*))
    Lsin= 2*np.sin(maxang/2)/(la)
    indMax= np.ceil([
        min([Lsin/(sum(ac**2)**(1/2))+1,ForceIndMax]),
        min([Lsin/(sum(bc**2)**(1/2))+1,ForceIndMax]),
        min([Lsin/(sum(cc**2)**(1/2))+1,ForceIndMax])
        ])
    print("indmax =" + str(indMax))

    multip=vecpicosf(Vecs, indMax,la)
    prov=np.copy(multip)
    end1=time.time()
    struct=matrizstrucf(multip, celda,la)
    end2=time.time()
    intensidades=fastintensidadesf(multip,struct)
    end3=time.time()
    patron=ord2f(intensidades)
    end=time.time()
    print("Tiempo de cálculo total: " + str(end-start))
    if ReturnMultip:
        return patron, prov
    else:
        return patron


def MultipPat(ctes,celda,multip,CheckAngs=False, la=la):
    start=time.time()
    "Regresa la lista ordenada de picos e intensidades relativas, dadas las condiciones"
    #Primero, se determinan los vectores unitarios en el sistema de ejes perpendiculares con 
    #x en dirección a, y sobre el plano definido por a y b, y el componente de c en z positivo.
    #Esto se hace porque la distancia interplanar es 1/|ha*+kb*+lc*|, con a*,b* y c* vectores recíprocos
    avec=np.array([ctes[1],0,0])
    #La fórmula de b se puede hallar usando a . b = abcos(gamma)
    bvec=np.array([ctes[2]*np.cos(ctes[6]),ctes[2]*np.sin(ctes[6]),0])
    #Para hallar el componente en x de c (cprov[0]), se usa a . c = accos(beta), y lo mismo para el componente y
    cprov=np.array([ctes[3]*np.cos(be),ctes[3]*((np.cos(ctes[4])-((np.cos(ctes[6]))*np.cos(ctes[5])))/np.sin(ctes[6]))])
    #El componente z se halla usando cx y cy, y usando que la magnitud de cvec es c
    cvec=np.array([cprov[0],cprov[1],-(ctes[3]**2-cprov[0]**2-cprov[1]**2)**(1/2)])*-1
    #el orden en la celda unitaria es átomo, ion, x, y, z
    #V es el volumen de la celda unidad
    V=np.cross(bvec,cvec)@avec
    # ac es el complemento de a, bc el de b y cc el de c
    ac=np.cross(bvec,cvec)/V
    bc=np.cross(cvec,avec)/V
    cc=np.cross(avec,bvec)/V
    Vecs=[ac,bc,cc]
    #Matriz de vectores recíprocos
    #El índice máximo se obtiene con la máxima reflexión del tipo x 0 0. Como Lsin/2d = ha* + kb* + lc*,
    #entonces indmax = int(max(Lsin/2a*, Lsin/2b*, Lsin/2c*))
    Lsin= 2*np.sin(maxang/2)/(la)
    if CheckAngs:
        multip[:,0]=vecangf(Vecs,multip[:,1:4])
    prov=np.copy(multip)
    end1=time.time()
    struct=matrizstrucf(prov, celda,la)
    intensidades=fastintensidadesf(prov,struct)
    patron=ord2f(intensidades)
    end=time.time()
    return patron

 
U=0.00
V=0.00
W=0.055

def Caglioti(ang,UVWg):
    "Fórmula de Caglioti. El ángulo de la entrada es en radianes"
    a=np.tan(2*np.pi * ang/360)
    return np.sqrt((UVWg[0]*(a**2))+(UVWg[1]*a)+UVWg[2])

#IMPORTANTE: no confundir gamma (fracción de caracter lorentziano) 
#            con ga (ángulo que forman a y b()

def PseVoi(pico,UVWg,ang):
    "Halla el valor de la distribución Pseudo-Voigt centrada en pico en el valor ang"
    #Se usa ang como 2theta
    ang2= np.tile(ang,(len(pico),1))
    dif=(ang2.T-pico.T).T/Caglioti(ang/2,UVWg)
    lor= (2/(np.pi*Caglioti(ang/2,UVWg)))/(1+((2*dif)**2))
    gau= (np.sqrt(4*np.log(2)/np.pi)/Caglioti(ang/2,UVWg)) * np.exp(-(dif**2)*4*np.log(2))
    return (UVWg[3]*lor) + ((1-UVWg[3])*gau)

def Intensidad (structs,ang,CheckAngs=False):
    "Función que devuelve la intensidad predicha en ang=2theta"
    "fracKa es el coeficiente que multiplica a Ka, no Ka/kb. Usar fracKa = (Ka/Kb)(1+ Ka/Kb)"
    partial = np.zeros_like(ang)
    for i in structs:
        patron =i.patron()
        partial = partial + (np.sum ((PseVoi(patron[:,0], i.UVWg, ang) * patron[:,4,np.newaxis]), axis=0 )/i.ks[0])
    return partial

group_cutoffs= np.array([2.1,15.1,74.1,142.1,167.1,194.1,230.1])

class estructura():
    def check(self, action):
        assert(type(self.gnum) is int)
        if action not in ["warn","max","min","first"]:
            raise ValueError("Acción inválida en check. opciones: warn, max, min, first")
        group_cutoffs= np.array([2.1,15.1,74.1,142.1,167.1,194.1,230.1])
        caso = np.searchsorted(group_cutoffs,self.gnum,side='left')

        if caso==0:
            self.refinableks=np.array([True,True,True,True,True,True,True])
        elif caso==1:
            if action=="warn":
                if np.array([np.abs(self.ks[4] - 90 *np.pi/180) < delta, 
                             np.abs(self.ks[6] - 90 *np.pi/180) < delta]).any():
                    raise AssertionError("La estructura no corresponde con un sistema monoclínico. Recuerde que el ángulo != 90° debe estar en la posición beta")
            self.ks[4]=90*np.pi/180
            self.ks[6]=90*np.pi/180
            self.refinableks=np.array([True,True,True,True,False,True,False])
        elif caso==2:
            if action=="warn":
                if np.array([np.abs(self.ks[4] - 90 *np.pi/180) > delta, 
                             np.abs(self.ks[5] - 90 *np.pi/180) > delta, 
                             np.abs(self.ks[6] - 90 *np.pi/180) > delta]).any():
                    raise AssertionError("La estructura no corresponde con un sistema ortorrómbico")
            self.ks[4]=90*np.pi/180
            self.ks[5]=90*np.pi/180
            self.ks[6]=90*np.pi/180
            self.refinableks=np.array([True,True,True,True,False,False,False])
        elif caso==3:
            if action=="warn":
                if np.array([np.abs(self.ks[4] - 90 *np.pi/180) > delta, 
                             np.abs(self.ks[5] - 90 *np.pi/180) > delta, 
                             np.abs(self.ks[6] - 90 *np.pi/180) > delta,
                             self.ks[1]!=self.ks[2]]).any():
                    raise ValueError("La estructura no corresponde con un sistema tetragonal. Recuerde que las longitudes iguales son a y b")
            elif action=="max":
                self.ks[1]=max(self.ks[1],self.ks[2])
                self.ks[2]=max(self.ks[1],self.ks[2])
            elif action=="min":
                self.ks[1]=min(self.ks[1],self.ks[2])
                self.ks[2]=min(self.ks[1],self.ks[2])
            elif action=="first":
                self.ks[2]=self.ks[1]
                
            self.ks[4]=90*np.pi/180
            self.ks[5]=90*np.pi/180
            self.ks[6]=90*np.pi/180
            self.refinableks=np.array([True,True,False,True,False,False,False])
            
        elif caso in (4,5): #Las limitaciones para los sistemas trigonales y hexagonales son los mismos
            if action=="warn":
                if np.array([np.abs(self.ks[4] - 90 *np.pi/180) > delta, 
                             np.abs(self.ks[5] - 90 *np.pi/180) > delta, 
                             np.abs(self.ks[6] - 120 *np.pi/180) > delta,
                             self.ks[1]!=self.ks[2]]).any():
                    raise AssertionError("La estructura no corresponde con un sistema hexagonal o trigonal")
            elif action=="max":
                self.ks[1]=max(self.ks[1],self.ks[2])
                self.ks[2]=max(self.ks[1],self.ks[2])
            elif action=="min":
                self.ks[1]=min(self.ks[1],self.ks[2])
                self.ks[2]=min(self.ks[1],self.ks[2])
            elif action=="first":
                self.ks[2]=self.ks[1]
                
            self.ks[4]=90*np.pi/180
            self.ks[5]=90*np.pi/180
            self.ks[6]=120*np.pi/180
            self.refinableks=np.array([True,True,False,True,False,False,False])
        elif caso ==6:
            if action=="warn":
                if np.array([np.abs(self.ks[4] - 90 *np.pi/180) > delta, 
                             np.abs(self.ks[5] - 90 *np.pi/180) > delta, 
                             np.abs(self.ks[6] - 90 *np.pi/180) > delta,
                             self.ks[1]!=self.ks[2],self.ks[2]!=self.ks[3]]).any():
                    raise AssertionError("La estructura no corresponde con un sistema cúbico")
            elif action=="max":
                self.ks[1]=max([self.ks[1],self.ks[2],self.ks[3]])
                self.ks[2]=max([self.ks[1],self.ks[2],self.ks[3]])
                self.ks[3]=max([self.ks[1],self.ks[2],self.ks[3]])
            elif action=="min":
                self.ks[1]=min([self.ks[1],self.ks[2],self.ks[3]])
                self.ks[2]=min([self.ks[1],self.ks[2],self.ks[3]])
                self.ks[3]=min([self.ks[1],self.ks[2],self.ks[3]])
            elif action=="first":
                self.ks[2]=self.ks[1]
                self.ks[3]=self.ks[1]
            self.ks[4]=90*np.pi/180
            self.ks[5]=90*np.pi/180
            self.ks[6]=90*np.pi/180
            self.refinableks=np.array([True,True,False,False,False,False,False])
        else:
            raise ValueError("El número de grupo es inválido")
            
    def __init__(self,ks,celda,gnum,grupo=None,UVWg=[0,0,0.05,0.7]):
        self.ks=ks
        self.celda=celda
        self.gnum=gnum
        self.grupo=grupo
        if grupo==None:
            self.grupo=tablas.Stdgroups[gnum]
        else:
            self.grupo=grupo
        self.FullCell= Wyckoff.custom_Wyckoff(self.grupo, celda)
        self.multip=None
        self.check("warn")
        self.refinableCel=None
        self.UVWg=UVWg
        self.refinableUVWg=np.array([0,0,1,0],dtype=bool)
        self.la =la
    def Setla (self, x):
        self.la=x
    def moveatom(self,varnum,cambio):
        "Mueve el átomo la cantidad correspondiente y retorna el movimiento como array de 6 componentes"
        a=np.zeros(6)
        if varnum%6 in (1,2,3):
            deltas=cambio*(tablas.DeltaCombinations[np.where(tablas.DeltaCombinations[:,varnum%6]!=0)])
            flag=True
            base=Wyckoff.custom_group(self.grupo, self.celda[int(varnum/6)])
            for i in deltas:
              if len(base) == len(Wyckoff.custom_group(self.grupo, self.celda[int(varnum/6)]+i)):
                  a=i
                  self.celda[int(varnum/6)] += a
                  self.FullCell = Wyckoff.custom_Wyckoff(self.grupo, self.celda)
                  flag=False
                  break
        elif varnum %6==5:
            a=np.array([0,0,0,0,0,cambio])
            self.celda[int(varnum/6),5]+=cambio;
        return a
         
    def patron(self,ForceIndMax=9999,CheckAngs=False):
        if self.multip is None:
            out, self.multip = MinInfoPat(self.ks,self.FullCell,ForceIndMax=9999,ReturnMultip=True, la = self.la)
            return out
        else:
            return MultipPat(self.ks, self.FullCell, self.multip, CheckAngs, la = self.la)
        
    def deriv(self, varnum, CellOrAtom, angs, constrained=True):
        "Encuentra la derivada parcial del patrón de difracción sobre los ángulos dados respecto a cierta variable"
        varnum= int(varnum)
        
        if CellOrAtom not in ["A","C","U"]:
            raise ValueError("Opción inválida. ingrese C para parámetro de celda o A para posición atómica")
        ProvStruct = estructura(np.copy(self.ks),np.copy(self.celda),self.gnum,self.grupo,np.copy(self.UVWg))
        #print(ProvStruct.ks)
        #print(ProvStruct.celda)
        if CellOrAtom=="C":
            locdelta = 1e-10
            ProvStruct.ks[varnum] = ProvStruct.ks[varnum] + locdelta
            if constrained:
                ProvStruct.check("max")
            res=(Intensidad([ProvStruct], angs,CheckAngs=True)-Intensidad([self], angs))/locdelta
            ProvStruct.ks=np.copy(self.ks)
            return res
                
        elif CellOrAtom == "A":
            locdelta = delta
            if (varnum % 6 == 0) or (varnum % 6 == 4):
                raise ValueError ("Error. No se puede derivar respecto al tipo de átomo o el ion")
            a=ProvStruct.moveatom(varnum,locdelta)
            if constrained and np.all(a=0):
                print("Error: No se encontró ninguna posición cercana con la misma multiplicidad. Intente nuevamente con la opción constrained=False")
                return np.zeros_like(angs)
            ProvStruct.celda[int(varnum/6)] -=a
            #Por alguna razón esta línea es necesaria. Parece que python optimiza al no borrar la variable local ProvStruc, así que hay que restar a celda.
            #Como FullCell no se actualiza, la derivada funciona como se espera
        elif CellOrAtom == "U":
            locdelta = 1e-10
            ProvStruct.UVWg[varnum] += locdelta
        return (Intensidad([ProvStruct],angs)-Intensidad([self],angs))/locdelta
    def idref(self):
        "Identifica qué parámetros de posición atómica son refinables. No está vectorizada completamente"
        if self.refinableCel is None:
            prov = np.zeros_like(self.celda,dtype=bool)
            prov[:,5] = (self.celda[:,5]!=1)
            for i in range(0,len(self.celda[:,0])):
                movimientos=np.array([[0,0,0,0,0,0.]])
                for j in range(1,4):
                    a=self.moveatom(6*i+j,delta)
                    self.moveatom(6*i+j,-delta)
                    #movimientos=np.vstack([movimientos,np.array([a])])
                    if np.all((np.any(a != movimientos, axis=1 ))):
                        movimientos=np.vstack([movimientos,np.array([a])])
                        prov[i,j]=True
            self.refinableCel=prov
            self.refinableUVWg = np.where(self.UVWg!=0)
    def iteration(self, datasource,flag=False):
        data_ints = datasource.ints-datasource.bgints
        weights= 1/datasource.ints
        self.idref()
        derivadasC=[]
        derivadasU=[]
        derivadasA=[]
        for i in np.where(self.refinableks)[0]:
            derivadasC.append(self.deriv(i,"C",datasource.angs)*weights)
        for i in np.where(self.refinableUVWg)[0]:
            derivadasU.append(self.deriv(i,"U",datasource.angs)*weights)
        
        for i in np.where(self.refinableCel.flatten())[0]:
            derivadasA.append(self.deriv(i,"A",datasource.angs)*weights)
        derivadas = derivadasC+derivadasU+derivadasA
        
        #Se dividen todas las derivadas entre la intensidad. Para restaurar que el peso es 1/I en lugar de
        #1/I^2, se multiplica por I al armar la matriz
        matriz = np.ones((len(derivadas),len(derivadas)))
        
        for i in range (0,len(derivadas)):
            for j in range(0,len(derivadas)):
                matriz[i,j]= np.sum(derivadas[i]*derivadas[j]/weights)
        delta = (Intensidad([self],datasource.angs) - data_ints)
        deltavec=np.ones((len(derivadas),1))
        #print(matriz)
        for i in range(0,len(derivadas)):
            deltavec[i,0]= np.sum(delta * derivadas[i])
        #print(deltavec)
        print(matriz)
        if flag:
            fig,ax=plt.subplots()
            ax.plot(datasource.angs,delta/datasource.ints)
            plt.show()
        res =  solve(matriz,deltavec)
        print(res)
        return res
Ti=estructura(ks1, cel1, 194)

prov=Wyckoff.semi_readCIF("Ba2SbSmO6.cif")

print (prov[1])
Ba2SbSmO6= estructura(prov[1],prov[2],prov[0],prov[3])
#Ba2SbSmO6.refinableks[0]=False
Ba2SbSmO6.ks[0]=1.02e-4

r=[]
As=[]
Xs=[]
k=[]
u=[]
Ba2SbSmO6.refinableCel=np.array([[0,0,0,0,0,0],
                           [0,0,0,0,0,0],
                           [0,0,0,0,0,0],
                           [0,1,0,0,0,0]],dtype=bool)
Ba2SbSmO6.refinableks = np.array([1,1,0,0,0,0,0])
suma = np.sum(ASCBa.ints)
Ba2SbSmO6.ks[[1,2,3]]=8.493
#Ba2SbSmO6.ks[[1,2,3]]=8.503
Ba2SbSmO6.celda[3,1]=0.211

Ba2SbSmO6.UVWg[2]=0.01
#Ba2SbSmO6.celda[3,1]=0.235

#Algo extraño está pasando cuando f1 no es 0 o 1. Revisar en dónde aparece el error
f1=1
ax.scatter(ASCBa.angs, Intensidad([Ba2SbSmO6],ASCBa.angs),label="Inicial",s=4,color="black")
r.append(np.sum(((np.abs(ASCBa.ints - ASCBa.bgints - Intensidad([Ba2SbSmO6],ASCBa.angs)))/suma)))
As.append(Ba2SbSmO6.ks[1])
Xs.append(Ba2SbSmO6.celda[3,1])
k.append(Ba2SbSmO6.ks[0])
u.append(Ba2SbSmO6.UVWg[2])

Ba2SbSmO6.FullCell = Wyckoff.custom_Wyckoff(Ba2SbSmO6.grupo, Ba2SbSmO6.celda)

NumIters=30
Ba2SbSmO6.iteration(ASCBa,True)

Ba2SbSmO6.refinableCel=np.array([[0,0,0,0,0,0],
                           [0,0,0,0,0,0],
                           [0,0,0,0,0,0],
                           [0,0,0,0,0,0]],dtype=bool)
Ba2SbSmO6.refinableks = np.array([1,1,0,0,0,0,0])
Ba2SbSmO6.refinableUVWg=np.array([0,0,0,0])
for i in range(0,30):
      print(i)
      vector = Ba2SbSmO6.iteration(ASCBa)
      Ba2SbSmO6.ks[0] -= vector[0,0]
      Ba2SbSmO6.ks[[1,2,3]] -= vector[1,0]
#      Ba2SbSmO6.UVWg[2] = np.abs(Ba2SbSmO6.UVWg[2] - vector[1,0])
#      Ba2SbSmO6.moveatom(19,-vector[1,0])
      Ba2SbSmO6.FullCell = Wyckoff.custom_Wyckoff(Ba2SbSmO6.grupo, Ba2SbSmO6.celda)
      r.append(np.sum(((np.abs(ASCBa.ints - ASCBa.bgints - Intensidad([Ba2SbSmO6],ASCBa.angs)))/suma)))
      As.append(Ba2SbSmO6.ks[1])
      Xs.append(Ba2SbSmO6.celda[3,1])
      k.append(Ba2SbSmO6.ks[0])
      u.append(Ba2SbSmO6.UVWg[2])

Ba2SbSmO6.refinableCel=np.array([[0,0,0,0,0,0],
                           [0,0,0,0,0,0],
                           [0,0,0,0,0,0],
                           [0,1,0,0,0,0]],dtype=bool)
Ba2SbSmO6.refinableks = np.array([1,1,0,0,0,0,0])
Ba2SbSmO6.refinableUVWg=np.array([0,0,1,0])
ax.scatter(ASCBa.angs, Intensidad([Ba2SbSmO6],ASCBa.angs),label="Intermedio",s=4,color="green")
for i in range(0,NumIters):
      print(i)
      vector = Ba2SbSmO6.iteration(ASCBa)
      Ba2SbSmO6.ks[0] -= vector[0,0]
      Ba2SbSmO6.ks[[1,2,3]] -= vector[1,0]
      Ba2SbSmO6.UVWg[2] = np.abs(Ba2SbSmO6.UVWg[2] - vector[2,0])
      Ba2SbSmO6.moveatom(19,-vector[3,0])
      Ba2SbSmO6.FullCell = Wyckoff.custom_Wyckoff(Ba2SbSmO6.grupo, Ba2SbSmO6.celda)
      r.append(np.sum(((np.abs(ASCBa.ints - ASCBa.bgints - Intensidad([Ba2SbSmO6],ASCBa.angs)))/suma)))
      As.append(Ba2SbSmO6.ks[1])
      Xs.append(Ba2SbSmO6.celda[3,1])
      k.append(Ba2SbSmO6.ks[0])
      u.append(Ba2SbSmO6.UVWg[2])

Ba2SbSmO6.refinableCel=np.array([[0,0,0,0,0,0],
                           [0,0,0,0,0,0],
                           [0,0,0,0,0,0],
                           [0,1,0,0,0,0]],dtype=bool)
Ba2SbSmO6.refinableks = np.array([1,1,0,0,0,0,0])
Ba2SbSmO6.refinableUVWg=np.array([0,0,1,0])
ax.scatter(ASCBa.angs, Intensidad([Ba2SbSmO6],ASCBa.angs),label="Intermedio2",s=4,color="orange")
for i in range(0,NumIters):
      print(i)
      vector = Ba2SbSmO6.iteration(ASCBa)
      Ba2SbSmO6.ks[0] -= vector[0,0]
      Ba2SbSmO6.ks[[1,2,3]] -= vector[1,0]
      Ba2SbSmO6.UVWg[2] = np.abs(Ba2SbSmO6.UVWg[2] - vector[1,0])
      Ba2SbSmO6.moveatom(19,-vector[2,0])
      Ba2SbSmO6.FullCell = Wyckoff.custom_Wyckoff(Ba2SbSmO6.grupo, Ba2SbSmO6.celda)
      r.append(np.sum(((np.abs(ASCBa.ints - ASCBa.bgints -Intensidad([Ba2SbSmO6],ASCBa.angs)))/suma)))
      As.append(Ba2SbSmO6.ks[1])
      Xs.append(Ba2SbSmO6.celda[3,1])
      k.append(Ba2SbSmO6.ks[0])
      u.append(Ba2SbSmO6.UVWg[2])



Ba2SbSmO6.iteration(ASCBa,True)   
  
ax.scatter(ASCBa.angs, Intensidad([Ba2SbSmO6],ASCBa.angs),label="Final",s=4,color="red")
ax.legend()
print(Ba2SbSmO6.ks)
print(Ba2SbSmO6.celda)
plt.legend()
# plt.title(r'$w_i = I^{-1}$')
fig,ax=plt.subplots()
# plt.title(r'$w_i = I^{-1}$')
ax.plot(np.arange(len(r)),np.array(r))
plt.xlabel("Iteraciones")
plt.ylabel("R")

fig,ax=plt.subplots()
ax.plot(np.arange(len(As)),np.array(As))
#plt.title(r'$w_i = I^{-1}$')
plt.xlabel("Iteraciones")
plt.ylabel("a (A)")
ax.plot(np.array([0,NumIters]),np.array([8.503,8.503]))

fig,ax=plt.subplots()
ax.plot(np.arange(len(Xs)),np.array(Xs))
#plt.title(r'$w_i = I^{-1}$')
plt.xlabel("Iteraciones")
plt.ylabel("x (rel. units)")
ax.plot(np.array([0,NumIters]),np.array([0.235,0.235]))

fig,ax=plt.subplots()
ax.plot(np.arange(len(k)),np.array(k))
#plt.title(r'$w_i = I^{-1}$')
plt.xlabel("Iteraciones")
plt.ylabel("k")

fig,ax=plt.subplots()
ax.plot(np.arange(len(u)),np.array(u))
#plt.title(r'$w_i = I^{-1}$')
plt.xlabel("Iteraciones")
plt.ylabel("u")



        



   