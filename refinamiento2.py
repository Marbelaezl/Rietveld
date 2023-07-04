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
from numpy.linalg import norm
mpl.use('TkAgg')

indmax=8
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
delta=0.0001
cel1=np.array([
    [22,0,0,0,0,1],
    [22,1/3,2/3,0.50,0,1]
    ])
ks1=[k,a,b,c,al,be,ga]
# Ninguna de las estructuras del TiO2 estaba dando los resultados correctos
# debido a un error en el cálculo de cvec. Buscar una nueva estructura e intentar de nuevo
ks2=[100,4.593,4.593,2.9590,90*np.pi/180,90*np.pi/180,90*np.pi/180]
ks2=[100,100,300,40,90*np.pi/180,90*np.pi/180,90*np.pi/180]
#La estructura del rutilo depende de un parámetro interno, que era el que aparecía
#como "x" en las posiciones de Wyckoff.
#El rutilo (ahora simulado correctamente) no parece poder explicar el pico súbito
#cerca a 70

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
# filename= "C:\\Users\migue\Downloads\HF 2.csv"
# arr=np.loadtxt(filename,skiprows=27,delimiter=",")
# datos=np.array([arr[:,0]+(arr[:,1]/(10**9)),arr[:,2]]).T
# datos=np.genfromtxt("HF2.csv",delimiter=",")
#------------------------------------------------------------------------------------------------------------
#Regresión polinómica para el fondo. Esta sección  está comentada para la versión con interfaz gráfica
# def background(data):
#     "Se espera una entrada de 2 x n; en el que la primera columna son ángulos y la segunda son intensidades. Determina el fondo"
#     #Primero, se escribe la matriz de datos
#     matriz=np.vstack((
#         data[:,0]**2,
#         data[:,0],
#         np.ones_like(data[:,0]),
#         1/(data[:,0]),
#         1/(data[:,0]**2),
#         1/(data[:,0]**3),
#         1/(data[:,0]**4),
#         1/(data[:,0]**5),
#         1/(data[:,0]**6))).T
#     #Vector de intensidades observadas
#     y=np.array([data[:,1]]).T
#     #Fórmula de regresión polinómica: (X^T X)^-1 X^T
#     MatFin= np.linalg.inv((matriz.T @ matriz)) @ matriz.T
#     #La estimación de los coeficientes es MatFin @ Intensidades.
#     return MatFin @ y
# datos_ruido=datos
# resfondo= background(datos)
# fondo=np.zeros_like(datos[:,0])
# #La interpretación de los coeficientes es tal que el fondo es fondo[0] + fondo[1]/theta + fondo[2]/theta^2 + ...

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
# #ax.plot(datos[:,0],datos[:,1]-fondo,linewidth=0.3)
# #ax.plot(datos[:,0],datos[:,1],linewidth=0.3)
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
    angs=(360*np.arcsin(la/(2*vecdisf(vecs,hkl))))/np.pi
    prov=np.hstack([np.array([angs]).T,hkl,np.array([np.ones_like(angs)]).T])
    indexes=np.where((prov[:,0]<maxang))
    return prov[indexes]
    



def fastpicoshelp(ind,flag=True):
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

    print(prov)
    
    unique, counts = np.unique(np.round(prov2[:,0],logd), return_counts=True)
    print(unique)
    prov[:,-1]+=counts-1
    return prov

def longpicosf(vecs,indmax):
    delta=1e-5
    res1=[]
    flag=True
    for l in np.arange(0,indmax[2]+1):
        if l%20==1:
            print ("índice actual:" +str(l-1))
        for k in np.arange(0,indmax[1]+1):
            for h in np.arange(0,indmax[0]+1):
                if [h,k,l]!=[0,0,0]:
#Esta condición controla los ángulos máximos y mínimos del patrón predicho
                    prov=angf(vecs,[h,k,l]) #Se guarda para no tener que calcularlo 3 veces
                    if prov[0]<=maxang and prov[0]>=minang:
                        for i in res1:
                            if np.abs(i[0]-prov[0])<delta:
                                flag=False
                        if flag:
                            res1.append(prov)
                        flag=True    
    #Factor de multiplicidad, no está funcionando. La implementación comentada funciona pero es lenta                                                                       
    cont=0
    res1=np.array(res1)
    res1[:,-1]=0
    res1=res1[res1[:, 0].argsort()]
    for l in np.arange(-indmax[2],indmax[2]+1):
        for k in np.arange(-indmax[1],indmax[1]+1):
            for h in np.arange(-indmax[0],indmax[0]+1):
                cont+=1
                prov=angf(vecs,[h,k,l])
                # if prov[0]<=maxang and prov[0]>=minang:
                #     for i in res1:
                #         if np.abs(i[0]-prov[0])<0.001:
                #             i[-1]+=1
                #             break
                
                prov[0]+=(delta/100)
                if prov[0]<=maxang and prov[0]>=minang:
                    index=np.searchsorted(res1[:,0], prov[0],side='right')
                    res1[index-1,-1]+=1
    print("Picos considerados por longpicosf:")
    print(cont)
    return res1


Cro_Mann=tablas.Cro_Mann
def structf(multip,celda):
    "Halla el factor F_hkl de los planos"
    def atom_strucf(n,par,ion=0):
        "Devuelve el factor de estructura del átomo con número n según el parámetro par = sin(theta)/lambda"
        "Solo considera la parte real; se omiten efectos de dispersión anómala"
        return (Cro_Mann[ion,n,0]*np.e**(-Cro_Mann[ion,n,4]*(par**2))+
              Cro_Mann[ion,n,1]*np.e**(-Cro_Mann[ion,n,5]*(par**2))+
              Cro_Mann[ion,n,2]*np.e**(-Cro_Mann[ion,n,6]*(par**2))+
              Cro_Mann[ion,n,3]*np.e**(-Cro_Mann[ion,n,7]*(par**2))+
              Cro_Mann[ion,n,8]
        )
#Las predicciones de factores de estructura son correctas y respetan las reglas de
#selección de picos
    b=[]
    res=0
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
    "Ordena los picos según el ángulo 2theta"
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
    
def MinInfoPat(ctes=ks1,celda=cel1,ForceIndMax=9999):
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
    Lsin= 2*np.sin(maxang)/(la)
    #print(Lsin)
    #print(sum(ac**2)**(1/2))
    indMax= np.ceil([
        min([Lsin/(sum(ac**2)**(1/2))+1,ForceIndMax]),
        min([Lsin/(sum(bc**2)**(1/2))+1,ForceIndMax]),
        min([Lsin/(sum(cc**2)**(1/2))+1,ForceIndMax])
        ])
    print("indmax =" + str(indMax))
    #multip=fastpicosf(Vecs,indMax)
    multip=vecpicosf(Vecs, indMax)
     #   print("Factor de multiplicidad encontrado. Picos obtenidos ("+str(len(multip))+"):")
     #   print("picos calculados correctamente. Picos encontrados:" + str(len(picos)))
    print("Factor de multiplicidad encontrado. Picos obtenidos ("+str(len(multip))+"):")
    print(multip)
    end1=time.time()
    print("Tiempo de cálculo: " + str(end1-start))
    struct=structf(multip,celda)
    end2=time.time()
    print("Tiempo de cálculo structf: " + str(end2-end1))
    
    #print("Factores de estructura encontrados:")
    #end3=time.time()
    #print("Tiempo de cálculo: " + str(end3-end2))
    
    intensidades=fastintensidadesf(multip,struct)
    end3=time.time()
    print("Tiempo de cálculo intensidadesf: " + str(end3-end2))
    
    #print("Intensidades sin ordenar:")

    #end4=time.time()
    #print("Tiempo de cálculo: " + str(end4-end3))
    
    patron=ord1f(intensidades)
    
    end=time.time()
    print("Tiempo de cálculo: " + str(end-start))
    return patron
#patron=MinInfoPat()
MinInfoPat(ks2)
    
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
#Entrada 9004141, con a y c diferentes e incompleta

print(" ")
#print (MinInfoPat(ks3,cel3))
#funcion=Intensidad(MinInfoPat(ks1,cel1),ks1[0],datos[:,0])
#funcion2=Intensidad(MinInfoPat(ks2,cel2),8000/50,datos[:,0])
#funcion3=Intensidad(entrada,65/500,datos[:,0])
#funcion=funcion+ fondo
#funcion4=Intensidad(MinInfoPat(ks3,cel3),ks3[0],datos[:,0])
#datos[:,1]=datos[:,1]-fondo
#diferencia= datos[:,1]-funcion
#print (diferencia)
#ax.plot(prov[:,0],diferencia* prov[:,1],linewidth=0.5)

#R=np.sum(np.abs(diferencia))/np.sum(np.abs(datos[:,1]))
#print ("R = " +str(R))
# prov=0
# prov2=0
#a=MinInfoPat(ks2,cel2)
# for i in range(0,min(len(a)-1,14)):
#     prov= prov+ np.abs(a[i+1,4]-entrada[i,4])
#     prov2=prov2 + entrada[i,4]
# print ("R_hkl = " + str(prov/prov2))

        
#ax.plot(datos[:,0],datos[:,1],label="datos",linewidth=0.2)
#ax.plot(datos[:,0],funcion,label="predicción (Ti)",linewidth=0.5)
#ax.plot(datos[:,0],funcion2,label="predicción (TiO2)",linewidth=1)
#ax.plot(datos[:,0],funcion4,label="predicción (Anatasa)",linewidth=0.5)
#ax.scatter(datos[:,0],funcion3,label="entrada 96-900-4142 COD",s=5,color="red")
#ax.legend()
#pyplot.annotate("R = " + str(R),xy=(10,1300))
#plt.savefig("filename.png")

#print(Intensidad(ks1,cel1,datos[:,0]))
#plt.show()

##Interfaz gráfica
#--------------------------------------------------------------------------------

#a=input()

   