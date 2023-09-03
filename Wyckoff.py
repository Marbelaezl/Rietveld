# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 18:34:05 2023

@author: migue
"""
import numpy as np
import tablas

a=223
celda= np.array([
    [22,0.1,0.2,0.3,0],
    [8,0.306,0.306,0,0]
    ])

#IMPORTANTE: LA SALIDA ES DE 3 VALORES, NO DE 6. 


n=4
def standard(arr):
    for i in [0,1,2]:
        while np.abs(arr[i])>=1 or arr[i]<0:
            if arr[i]>=1:
                arr[i]=round(arr[i]-1,n)
            if arr[i]<0:
                arr[i]=round(arr[i]+1,n)
        arr[i]=round(arr[i],n)
    return arr



    
def stol(algo,delim=","):
    #string to list
    res=[]
    prov=""
    for i in algo:
        if i==delim:
            res.append(prov)
            prov=""
        else:

            prov=prov+i
    res.append(prov)
    return res

Stdgroups=tablas.Stdgroups
    
    
def custom_group(grupo,entrada):
    x=str(entrada[1])
    y=str(entrada[2])
    z=str(entrada[3])
    res=[]
    prov=[]
    g2=[]
    for i in grupo:
        for j in stol(i,","):
            prov.append((float(eval(j.replace("x",x).replace("y",y).replace("z",z)))))
        g2.append(prov)
        prov=[]
    flag=False
    delta=0.0001
    for i in g2:
        for j in res:
            if np.linalg.norm(np.array(standard(i))-np.array(j))<delta:
                flag=True
        if flag==False:
            #print(str(i)+ "->" +str(standard(i)))
            res.append(standard(i))
        flag=False
    prov=[]
    if len(entrada)==6:
        for i in res:
            prov.append([entrada[0],i[0],i[1],i[2],entrada[4],entrada[5]])
        return prov
    elif len(entrada)==5:
        for i in res:
            prov.append([entrada[0],i[0],i[1],i[2],entrada[4]])
        return prov


def atom(num, entrada):
    return custom_group(Stdgroups[num], entrada)


def Wyckoff(n, celda):
    "Construye la celda unitaria a partir de la entrada equivalente y el grupo espacial"
    c=[]
    for j in celda:
        c=c+ atom(n,j)
    return np.array(c)
#print (Wyckoff(136,celda))

def custom_Wyckoff(grupo,celda):
    "Equivalente a Wyckoff para grupos personalizados"
    c=[]
    for i in range(0,len(celda)):
        c=c+ custom_group(grupo, celda[i])
    return np.array(c)


elementos=tablas.elementos

#print(standard([0.5,0.3,0.7]))
#print(Wyckoff(2,[[53,0.1,0.2,0.3,0,1]]))

def Standard_ion(non_standard,mute=False):
    "Convierte un ion a forma estándar. La entrada son strings!!!"
    assert(type(non_standard)==str)
    #Por si acaso, porque el usuario me cae mal >:C
    elemento=""
    ion=""
    for i in non_standard:
        if i==" ":
            pass
        elif i in "1234567890":
            ion=ion+i
        elif i in "+-":
            ion= i + ion
        else:
            elemento=elemento+i
    try:
        elementos[elemento+ion]
    except:
        try:
            elementos[elemento]
            if mute==False:
                print((elemento+ion)+ " no es un ion reconocido. Iones disponibles para el " + elemento + ":")
                for i in range(0,5):
                    try:
                        if mute==False:
                            print(list(elementos.keys())[list(elementos.values()).index([elementos[elemento][0],i])])
                    except:
                        return None
        except:
            if mute==False:
                print (elemento + " no es un elemento reconocido.")
            return None
    return elemento+ion



def semi_readCIF(nombre):
    try:
        prov1=[]
        prov2=[]
        grupo=[]
        flag=0
        constantes=[1,0,0,0,0,0,0]
        celda=[]
        g=1
        file1=open(nombre,"r")
        Lines = file1.readlines()
        for i in Lines:
            
            if "_symmetry_Int_Tables_number" in i:
                g=int(i.replace("_symmetry_Int_Tables_number"," ").strip())
                print("El archivo fue leído, pero corresponde a un estándar cif antiguo o tiene opciones de retrocompatibilidad. Esto puede afectar ligeramente el funcionamiento")
    
                
            if "_space_group_IT_number" in i:
                g=int(i.replace("_space_group_IT_number"," ").strip())
    
                
            if "_cell_angle_alpha" in i:
                constantes[4]=float(i.replace("_cell_angle_alpha "," "))*np.pi/180
    
            
            if "_cell_angle_beta" in i:
                constantes[5]=float(i.replace("_cell_angle_beta "," "))*np.pi/180
    
                
            if "_cell_angle_gamma" in i:
                constantes[6]=float(i.replace("_cell_angle_gamma "," "))*np.pi/180
    
              
            if "_cell_length_a" in i:
                constantes[1]=float(i.replace("_cell_length_a"," "))
                
            if "_cell_length_b" in i:
                constantes[2]=float(i.replace("_cell_length_b"," "))
                
            if "_cell_length_c" in i:
                constantes[3]=float(i.replace("_cell_length_c"," "))
                
            if i.strip()=="loop_":
    
                if prov1!=[]:
                    prov2.append(prov1)
                    prov1=[]
                flag=1
                
            if i.strip()==";":
                if prov1!=[]:
                    prov2.append(prov1)
                    prov1=[]
                flag=0
                
            if flag:
                prov1.append(i.strip())
        prov2.append(prov1)
        prov1=[0,0,0,0,0]   
        for j in prov2:
            if j[1]=="_symmetry_equiv_pos_as_xyz" or j[1]=="_space_group_symop_operation_xyz":
                grupo=j[2:]
                #print(grupo)
            if "_atom_site" in j[1]:
                for i in range(0,len(j)):
                    if j[i]=="_atom_site_type_symbol":
                        prov1[0]=i-1
                    if j[i]=="_atom_site_fract_x":
                        prov1[1]=i-1
                    if j[i]=="_atom_site_fract_y":
                        prov1[2]=i-1
                    if j[i]=="_atom_site_fract_z":
                        prov1[3]=i-1
                    if j[i]=="_atom_site_occupancy":
                        prov1[4]=i-1
                    if ("_atom_site" not in j[i]) and j[i]!="loop_":
                        lista=stol(j[i]," ")
                        celda.append([
                            elementos[Standard_ion(lista[prov1[0]])][0],
                            float(lista[prov1[1]]),
                            float(lista[prov1[2]]),
                            float(lista[prov1[3]]),
                            elementos[Standard_ion(lista[prov1[0]])][1],
                            float(lista[prov1[4]])
                                     ])
        #print(custom_group(grupo, celda[0]))
        return([g,np.array(constantes),np.array(celda),grupo])
    except:
        raise ValueError ("Error al leer el archivo")
