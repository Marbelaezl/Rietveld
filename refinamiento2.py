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
#mpl.use('TkAgg')
#IMPORTANTE PARA EL EJECUTABLE AAAAAAAAAAAAAAA :)

indmax=8
minang=10
maxang=80

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
delta=0.0001
cel1=np.array([
    #[22,0,0,0.25,0,1],
    [22,1/3,2/3,0.25,0,1]
    ])
ks1=[k,a,b,c,al,be,ga]
# Ninguna de las estructuras del TiO2 estaba dando los resultados correctos
# debido a un error en el cálculo de cvec. Buscar una nueva estructura e intentar de nuevo
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
#print(cel3)
#Predicciones correctas: Ti(alpha), diamante, NaCl, TiO2
#Resultados = open("Resultados.dat",mode="w")
#Datos_ruido=open("Datos_ruido.dat",mode="w")
# fig, ax=plt.subplots()
# ax.set_xlabel(r'ángulo $(2 \theta)$ ')
# ax.set_ylabel("Intensidad relativa")
filename= "C:\\Users\migue\Downloads\HF 2.csv"
arr=np.loadtxt(filename,skiprows=27,delimiter=",")
#data=np.array([arr[:,0]+(arr[:,1]/(10**9)),arr[:,2]]).T
data=np.genfromtxt("HF2.csv",delimiter=",")
ba2sm=np.genfromtxt("Ba2SbSmO6.ASC",delimiter=" ")
#------------------------------------------------------------------------------------------------------------
#Regresión polinómica para el fondo. Esta sección  está comentada para la versión con interfaz gráfica
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
ax.plot(HF2.angs,HF2.ints,linewidth=0.5)
plt.xlabel(r'$2\theta($°)')
plt.ylabel("Intensidad (counts)")
fig,ax=plt.subplots()
ax.plot(ASCBa.angs,ASCBa.ints-ASCBa.bgints)
#ax.plot(ASCBa.angs,ASCBa.bgints)





#resfondo= background(datos)
#fondo=np.zeros_like(datos[:,0])
#La interpretación de los coeficientes es tal que el fondo es fondo[0] + fondo[1]/theta + fondo[2]/theta^2 + ...

# for i in range(0,len(resfondo)):
#     fondo=fondo + (resfondo[i]/(datos[:,0]**(i-2)))
# print(np.std(datos_ruido))
# sd=np.std(datos[:,1]-fondo)
# iterador=0
# intensos=[]
# while iterador<len(datos_ruido):
#     if (datos_ruido[:,1]-fondo)[iterador] > sd:
#         intensos.append(datos_ruido[iterador])
#         datos_ruido=np.delete(datos_ruido,iterador,axis=0)
#         fondo=np.delete(fondo,iterador,axis=0)
#     else:
#         iterador=iterador+1
# #with np.printoptions(threshold=np.inf):
# #    print(datos_ruido,file=Datos_ruido)

# print (np.array(intensos))

# resfondo=background(datos_ruido.reshape((-1,2)))
# fondo=np.zeros_like(datos[:,0])
# fondo1=np.zeros_like(datos_ruido[:,0])
# for i in range(0,len(resfondo)):
#     fondo=fondo + (resfondo[i]/(datos[:,0]**(i-2)))
#     fondo1=fondo1+(resfondo[i]/(datos_ruido[:,0]**(i-2)))

# SDRuido=np.std(datos_ruido[:,1]-fondo1)
# print("Desviación estándar sin los picos: " + str(SDRuido))
# #with np.printoptions(threshold=np.inf):
# #    print(np.vstack([datos[:,0],datos[:,1]-fondo]).T,file=Resultados)
# #print (np.std(datos[:,1]-fondo))
# # ax.plot(datos[:,0],datos[:,1])
# # ax.plot(datos[:,0],fondo)
# print(datos[:,1]-fondo)
#ax.plot(datos[:,0],datos[:,1]-fondo,linewidth=0.3)
#ax.plot(datos[:,0],datos[:,1],linewidth=0.3)
# #Parece que ya funciona, pero no lo he implementado completamente porque no tengo los datos completos
#------------------------------------------------------------------------------------------------------------------------------------------------------------
#Función de peso, basada en la desviación estándar
#Implementación actual: sin pesos

def peso(datos):
    # prov=np.array([datos[:,0],datos[:,1]-fondo]).T
    # for i in range(0,len(datos)):
    #     prov[i,1]= max(1,(abs(prov[i,1])-(3*SDRuido)))
    prov=np.array([datos[:,0],1/datos[:,1]])
    prov=np.array([datos[:,0],np.ones_like(datos[:,1])])
    return prov.T
#prov=peso(datos)
#ax.plot(prov[:,0],prov[:,1],color="black")


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


def angf(vecs,hkl):
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

def vecangf(vecs,hkl):
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

def fastpicosf(vecs,indmax):
    #Implementación relativamente rápida, pero no completamente vectorizada. vecpicosf debe hacer lo mismo que
    #esta función, no al contrario.
    delta=1e-5
    res=np.array([[0.0,0.0,0.0,0.0,1]])
    for i in fastpicoshelp(indmax):
        prov=angf(vecs,i)
        ii=np.searchsorted(res[:,0],prov[0],side='right')
        if prov[0]!=-1:
            if np.abs(prov[0]-res[ii-1,0])>delta:
                res=np.insert(res,ii,prov,axis=0)
            else:
                res[ii-1,-1]+=1
    minii=np.searchsorted(res[:,0],minang,side='left')
    maxii=np.searchsorted(res[:,0],maxang,side='right')
    res=res[minii:maxii]
    res[:,-1]=0
    for i in fastpicoshelp(indmax,flag=False):
        prov=angf(vecs,i)
        prov[0]+=(delta/100)
        ii=np.searchsorted(res[:,0],prov[0],side='right')
        if np.abs(prov[0]-res[ii-1,0])<delta:
            res[ii-1,-1]+=1
    return (res)

def vecpicosf(vecs,indmax):
    logd=6 #logaritmo de delta
    delta=10**-logd
    prov=vecangf(vecs,fastpicoshelp(indmax))
    prov=prov[prov[:, 0].argsort()]
    prov = np.delete(prov, np.argwhere(np.ediff1d(prov[:,0].round(logd)) < delta/10) + 1,axis=0)
    prov2=vecangf(vecs,fastpicoshelp(indmax,flag=False))
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
def structf(multip,celda):
    "Halla el factor F_hkl de los planos"
    b=[]
    res=0
    #Esto usa un for,  probablemente ineficiente para celdas grandes. Estoy casi seguro de que el problema
    #puede reducirse a una multiplicación de matrices
    for i in multip:
        for j in celda:
            # Se halla hx + ky + lz
            fase= (i[1]*j[1])+(i[2]*j[2])+(i[3]*j[3])
            #Factor de estructura = suma a lo largo de j: f_j * (e^2(pi)j*fase)
            prov = (atom_strucf(int(j[0]),(np.sin(i[0]*np.pi/360)/la),int(j[4]))*(np.exp((0+2j)*np.pi*fase)))
            prov=prov*j[5]
            res=res+ prov
        b.append(res)

        res=0
    return np.array(b)

def matrizstrucf(multip,celda):
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
    
def intensidadesf(multip,struct):
    "Convierte amplitudes complejas a intensidades: I proporcional a F_hkl^2 *" 
    "factor de Lorentz * factor de multiplicidad"
    #Las predicciones de intensidad a partir de la amplitud compleja son calculadas
    #de la manera esperada. El factor de multiplicidad se añade en el siguiente paso,
    #cuando se suman las intensidades de los picos con el mismo ángulo
    for i in range(0,len(multip)):
        try:
        #como la entrada a[i,0] es 2theta, se multiplica por pi/360 en lugar de pi/180
            ang=multip[i,0]*np.pi/360
            #lo: factor de lorentz y polarización. Se deja 1 antes de asignar el valor 
            # correcto para facilitar la corrección del programa
            lo=1
            #Geometría Bragg-Brentano
            lo= (1+ (np.cos(2*ang)**2))/(np.sin(2*ang)*np.sin(ang))
            multip[i,4]=np.abs(struct[i])**2 *lo*multip[i,4]
        except:
            pass
        if max(multip[:,4]>1000):
            multip[:,4]=multip[:,4]/max(multip[:,4])
    multip[:,4]=multip[:,4]/max(multip[:,4])
    return multip

def fastintensidadesf(multip,struct):
    #como la entrada a[i,0] es 2theta, se multiplica por pi/360 en lugar de pi/180
    ang=multip[:,0]*np.pi/360
    #Factor de Lorentz - Geometría Bragg-Brentano 
    lo= (1+ (np.cos(2*ang)**2))/(np.sin(2*ang)*np.sin(ang))
    multip[:,4]=np.abs(struct)**2 *lo*multip[:,4]
    multip[:,4]=multip[:,4]/max(multip[:,4])
    return multip
    
    

def ord1f(intens):
    "Ordena los picos según el ángulo 2theta. Implementación recursiva y no muy buena"
    cambios=0
    for i in range(1,len(intens)):
        try:
            if intens[i,4]<0.01:
                intens=np.delete(intens,i,axis=0)
                cambios=cambios+1
            if intens[i,0]< intens[i-1,0]:
                #a, a+b
                intens[i-1]=intens[i-1]+intens[i]
                #b, a+b
                intens[i]=intens[i-1]-intens[i]
                intens[i-1]=intens[i-1]-intens[i]
                cambios=cambios+1
        except:
            pass
    if cambios==0:
        return intens
    return ord1f(intens)
    
def ord2f(intens):
    "Ordena los picos según el ángulo 2theta"
    intens=np.delete(intens,np.where(intens[:,-1]<0.001),axis=0)
    intens=intens[np.argsort(intens[:,0])]
    return intens


def MinInfoPat(ctes=ks1,celda=cel1,ForceIndMax=9999, ReturnMultip=False):
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

    multip=vecpicosf(Vecs, indMax)
    #print("Factor de multiplicidad encontrado. Picos obtenidos : "+str(len(multip)))
    prov=np.copy(multip)
    end1=time.time()
    #print("Tiempo de cálculo picos y multiplicidad: " + str(end1-start))
    struct=matrizstrucf(multip, celda)
    #print(struct)
    #print(matrizstrucf(prov, celda))
    end2=time.time()
    #print("Tiempo de cálculo structf: " + str(end2-end1))
    intensidades=fastintensidadesf(multip,struct)
    end3=time.time()
    #print("Tiempo de cálculo intensidadesf: " + str(end3-end2))
    
    #print("Intensidades sin ordenar:")

    #end4=time.time()
    #print("Tiempo de cálculo: " + str(end4-end3))
    
    patron=ord2f(intensidades)
    end=time.time()
    print("Tiempo de cálculo total: " + str(end-start))
    if ReturnMultip:
        return patron, prov
    else:
        return patron


def MultipPat(ctes,celda,multip,CheckAngs=False):
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
        multip[:,0]=vecdisf(Vecs,multip[:,1:4])
    prov=np.copy(multip)
    end1=time.time()
    struct=matrizstrucf(prov, celda)
    end2=time.time()
    #print("Tiempo de cálculo structf: " + str(end2-end1))
    intensidades=fastintensidadesf(prov,struct)
    end3=time.time()
    #print("Tiempo de cálculo intensidadesf: " + str(end3-end2))
    patron=ord2f(intensidades)
    end=time.time()
    #print("Tiempo de cálculo total: " + str(end-start))
    return patron

#pat1, multip1 = MinInfoPat(ks2,cel2,ReturnMultip=True)
#pat2 = MultipPat(ks2, cel2, multip1)

    
U=0.00
V=0.00
W=0.01


def Caglioti(ang):
    "Fórmula de Caglioti. El ángulo de la entrada es en radianes"
    a=np.tan(ang)
    return np.sqrt((U*(a**2))+(V*a)+W)

gamma=0.5
#IMPORTANTE: no confundir gamma (fracción de caracter lorentziano) 
#            con ga (ángulo que forman a y b()

def PseVoi(pico,ang):
    "Halla el valor de la distribución Pseudo-Voigt centrada en pico en el valor ang"
    #Se usa ang como 2theta
    dif=(ang-pico)/Caglioti(ang/2)
    lor= (2/(np.pi*Caglioti(ang/2)))/(1+((2*dif)**2))
    gau= (np.sqrt(4*np.log(2)/np.pi)/Caglioti(ang/2)) * np.exp(-(dif**2)*4*np.log(2))
    return (gamma*lor) + ((1-gamma)*gau)

def Intensidad (patron,k,ang):
    "Función que devuelve la intensidad predicha en ang=2theta"
    suma=0
    for i in patron:
        suma=suma+(i[4]*PseVoi(i[0],ang))
    return k*suma
#IMPORTANTE: FALTA AÑADIR LA INTENSIDAD DE FONDO AL CÁLCULO
#IMPORTANTE: EL PROGRAMA NO CONSIDERA K_ALPHA Y K_BETA
#print (MinInfoPat(ks2,cel2))
entrada=np.array([
    [27.44,1,1,0,1],
    [36.079,1,0,1,.471],
    [39.197,2,0,0,.0472],
    [41.240,1,1,1,.183],
    [44.050,2,1,0,0.065],
    [54.325,2,1,1,0.583],
    [56.635,2,2,0,0.163],
    [62.752,0,0,2,0.083],
    [64.058,3,1,0,0.080],
#    [65.517,2,2,1,0.006],
    [69.012,3,0,1,0.195],
    [69.797,1,1,2,0.104],
    [72.427,3,1,1,0.01],
#    [74.414,3,2,0,0.002],
    [76.536,2,0,2,0.024],
    [79.830,2,1,2,0.012]
    
    ])

group_cutoffs= np.array([2.1,15.1,74.1,142.1,167.1,194.1,230.1])
print(np.searchsorted(group_cutoffs, 210))


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
            
    def __init__(self,ks,celda,gnum,grupo=None):
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
            out, self.multip = MinInfoPat(self.ks,self.FullCell,ForceIndMax=9999,ReturnMultip=True)
            return out
        else:
            return MultipPat(self.ks, self.FullCell, self.multip, CheckAngs)
        
    def deriv(self, varnum, CellOrAtom, angs, constrained=True):
        "Encuentra la derivada parcial del patrón de difracción sobre los ángulos dados respecto a cierta variable"
        varnum= int(varnum)
        if CellOrAtom not in ["A","C"]:
            raise ValueError("Opción inválida. ingrese C para parámetro de celda o A para posición atómica")
        ProvStruct = estructura(np.copy(self.ks),np.copy(self.celda),self.gnum,self.grupo)
        #print(ProvStruct.ks)
        #print(ProvStruct.celda)
        if CellOrAtom=="C":
            ProvStruct.ks[varnum] = ProvStruct.ks[varnum] + delta
            if varnum==0:
                return Intensidad(self.patron(), 1, angs)
            elif constrained:
                ProvStruct.check("max")
            res=(Intensidad(ProvStruct.patron(CheckAngs=True), ProvStruct.ks[0], angs)-Intensidad(self.patron(), self.ks[0], angs))/delta
            ProvStruct.ks=np.copy(self.ks)
            return res
                
        elif CellOrAtom == "A":
            if (varnum % 6 == 0) or (varnum % 6 == 4):
                raise ValueError ("Error. No se puede derivar respecto al tipo de átomo o el ion")
            a=ProvStruct.moveatom(varnum,delta)
            if constrained and np.all(a=0):
                print("Error: No se encontró ninguna posición cercana con la misma multiplicidad. Intente nuevamente con la opción constrained=False")
                return np.zeros_like(angs)
            ProvStruct.celda[int(varnum/6)] -=a
            #Por alguna razón esta línea es necesaria. Parece que python optimiza al no borrar la variable local ProvStruc, así que hay que restar a celda.
            #Como FullCell no se actualiza, la derivada funciona como se espera
            
            return (Intensidad(ProvStruct.patron(), ProvStruct.ks[0], angs)-Intensidad(self.patron(), self.ks[0], angs))/delta
        else:
            return (Intensidad(ProvStruct.patron(), ProvStruct.ks[0], angs)-Intensidad(self.patron(), self.ks[0], angs))/delta
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
    def iteration(self, datasource):
        data_ints = datasource.ints-datasource.bgints
        weights= 1/data_ints
        self.idref()
        derivadasC=[]
        derivadasA=[]
        for i in np.where(self.refinableks)[0]:
            derivadasC.append(self.deriv(i,"C",datasource.angs)*weights)
        for i in np.where(self.refinableCel.flatten())[0]:
            derivadasA.append(self.deriv(i,"A",datasource.angs)*weights)
        derivadas = derivadasC+derivadasA
        
        #Se dividen todas las derivadas entre la intensidad. Para restaurar que el peso es 1/I en lugar de
        #1/I^2, se multiplica por I al armar la matriz
        matriz = np.ones((len(derivadas),len(derivadas)))
        
        for i in range (0,len(derivadas)):
            for j in range(0,len(derivadas)):
                matriz[i,j]= np.sum(derivadas[i]*derivadas[j]/weights)
        delta = (Intensidad(self.patron(), self.ks[0], datasource.angs) - data_ints)
        deltavec=np.ones((len(derivadas),1))
        #print(matriz)
        for i in range(0,len(derivadas)):
            deltavec[i,0]= np.sum(delta * derivadas[i])
        #print(deltavec)
        res =  solve(matriz,deltavec)
        print(res)
        return res
Ti=estructura(ks1, cel1, 194)
#print(Ti.FullCell)
#print(Ti.patron())

prov=Wyckoff.semi_readCIF("Ba2SbSmO6.cif")

print (prov[1])
Ba2SbSmO6= estructura(prov[1],prov[2],prov[0],prov[3])
#Ba2SbSmO6.refinableks[0]=False
Ba2SbSmO6.ks[0]=4800

r=[]
As=[]
Xs=[]
suma = np.sum(ASCBa.ints)
Ba2SbSmO6.ks[[1,2,3]]=8.463
#Ba2SbSmO6.ks[[1,2,3]]=8.503
Ba2SbSmO6.celda[3,1]=0.211
#Ba2SbSmO6.celda[3,1]=0.235
ax.scatter(ASCBa.angs, Intensidad(Ba2SbSmO6.patron(), Ba2SbSmO6.ks[0], ASCBa.angs),label="Inicial",s=4,color="black")
r.append(np.sum(((np.abs(ASCBa.ints - ASCBa.bgints - Intensidad(Ba2SbSmO6.patron(), Ba2SbSmO6.ks[0], ASCBa.angs)))/suma)))
As.append(Ba2SbSmO6.ks[1])
Xs.append(Ba2SbSmO6.celda[3,1])

Ba2SbSmO6.FullCell = Wyckoff.custom_Wyckoff(Ba2SbSmO6.grupo, Ba2SbSmO6.celda)

NumIters=300

for i in range(0,NumIters):
      print(i)
      vector = Ba2SbSmO6.iteration(ASCBa)
      Ba2SbSmO6.ks[0] -= vector[0,0] #REVISAR EL SIGNO DE ESTA LÍNEA AAAAAAAAAAAAAAAAA!!!!!!!!!!!!!!!!
      Ba2SbSmO6.ks[[1,2,3]] -= vector[1,0]
      Ba2SbSmO6.moveatom(19,-vector[2,0])
      Ba2SbSmO6.FullCell = Wyckoff.custom_Wyckoff(Ba2SbSmO6.grupo, Ba2SbSmO6.celda)
      r.append(np.sum(((np.abs(ASCBa.ints - ASCBa.bgints - Intensidad(Ba2SbSmO6.patron(), Ba2SbSmO6.ks[0], ASCBa.angs)))/suma)))
      As.append(Ba2SbSmO6.ks[1])
      Xs.append(Ba2SbSmO6.celda[3,1])
      if i==500:
          constant500=[Ba2SbSmO6.ks[0],Ba2SbSmO6.ks[1],Ba2SbSmO6.celda[3,1]]
      
ax.scatter(ASCBa.angs, Intensidad(Ba2SbSmO6.patron(), Ba2SbSmO6.ks[0], ASCBa.angs),label="Final",s=4,color="red")
    
    
ax.scatter(ASCBa.angs, Intensidad(Ba2SbSmO6.patron(), Ba2SbSmO6.ks[0], ASCBa.angs),label="Final",s=4,marker="s")
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

#print(constant500)
# TiO2= estructura(ks2,cel,136)
# TiO2.idref()
# print(TiO2.refinableCel)
# print(TiO2.patron())

# print (TiO2.iteration(HF2))
    
# ax.legend()
# fig,ax=plt.subplots()
# #ax.plot(np.arange(0,90,0.1),Intensidad(TiO2.patron(), 1, np.arange(0,90,0.1)))
# print(TiO2.moveatom(8, delta))
# ax.plot(np.arange(25,100,0.2),TiO2.deriv(8,"A", np.arange(25,100,0.2)))
# # plt.show()
# print (TiO2.FullCell)
        


        



   