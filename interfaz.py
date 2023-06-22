# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 08:51:27 2023

@author: migue
"""
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import Wyckoff
import tkinter as tk
from tkinter import filedialog
import tablas
import time
import refinamiento2 as ref

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
    [22,0,0,0,0],
    [22,1/3,2/3,0.50,0]
    ])
ks1=[k,a,b,c,al,be,ga]


ks2=[100,4.593,4.593,2.9590,90*np.pi/180,90*np.pi/180,90*np.pi/180]

n=1
cel=np.array([
      [22,0,0,0,0,1],
      [8,0.306,0.306,0,0,n],
      ])
cel2=Wyckoff.Wyckoff(136, cel)
ks3=[100,3.7850 ,3.7850 ,9.520,90*np.pi/180,90*np.pi/180,90*np.pi/180]
cel3=Wyckoff.Wyckoff(141, np.array([[22,0,0-1/4,0-3/8,0,1],[8,0,0-1/4,0.20806-3/8,0,1]]))


root = tk.Tk()
# Labels
title = tk.Label(root,text="Programa en desarrollo")
grouplabel=tk.Label(root,text="Grupo espacial")

alabel=tk.Label(root,text="a(Å)")
blabel=tk.Label(root,text="b(Å)")
clabel=tk.Label(root,text="c(Å)")

alphalabel=tk.Label(root,text="α(°)")
betalabel=tk.Label(root,text="β(°)")
gammalabel=tk.Label(root,text="γ(°)")

cellabel=tk.Label(root,text="Celda unitaria")
Zlabel=tk.Label(root,text="Símbolo")
xlabel=tk.Label(root,text="x")
ylabel=tk.Label(root,text="y")
zlabel=tk.Label(root,text="z")
ulabel=tk.Label(root,text="U")

cellabel2=tk.Label(root,text="")

title.grid(row=0,column=3)

grouplabel.grid(row=1,column=0)

alabel.grid(row=3,column=0)
blabel.grid(row=3,column=1)
clabel.grid(row=3,column=2)

alphalabel.grid(row=3,column=3)
betalabel.grid(row=3,column=4)
gammalabel.grid(row=3,column=5)

cellabel.grid(row=5,column=0)
Zlabel.grid(row=6,column=0)
xlabel.grid(row=6,column=1)
ylabel.grid(row=6,column=2)
zlabel.grid(row=6,column=3)
ulabel.grid(row=6,column=4)
cellabel2.grid(row=5,column=5)

#Textboxes
groupentry=tk.Entry(root)
aentry=tk.Entry(root)
bentry=tk.Entry(root)
centry=tk.Entry(root)
alphaentry=tk.Entry(root)
betaentry=tk.Entry(root)
gammaentry=tk.Entry(root)
Params=[aentry,bentry,centry,alphaentry,betaentry,gammaentry]
def switch_params(state):
    for i in Params:
        if state==0:
          i.config(state="disabled")
        elif state==1:
            i.config(state="normal")
switch_params(0)
Zentry=tk.Entry(root)
xentry=tk.Entry(root)
yentry=tk.Entry(root)
zentry=tk.Entry(root)
uentry=tk.Entry(root)
CellParams=[Zentry,xentry,yentry,zentry,uentry]
def switch_cparams(state):
    for i in CellParams:
        if state==0:
          i.config(state="disabled")
        elif state==1:
            i.config(state="normal")
switch_cparams(0)
groupentry.grid(row=2,column=0)
aentry.grid(row=4,column=0)
bentry.grid(row=4,column=1)
centry.grid(row=4,column=2)
alphaentry.grid(row=4,column=3)
betaentry.grid(row=4,column=4)
gammaentry.grid(row=4,column=5)
Zentry.grid(row=7,column=0)
xentry.grid(row=7,column=1)
yentry.grid(row=7,column=2)
zentry.grid(row=7,column=3)
uentry.grid(row=7,column=4)
#Buttons
g=1
def groupbut():
    global g
    global ks1
    g=groupentry.get()
    try:
        g=int(g)
        gl2= tk.Label(root, text="Grupo seleccionado: "+str(g)).grid(row=2,column=1)
    except:
        print("Error. Por favor ingrese el número del grupo espacial")
        switch_params(0)
    if g in range(1,3):
        #Monoclínicos
        switch_params(1)
        switch_cparams(1)
    elif g in range(3,16):
        #Triclínicos
        ks1[4]=90*np.pi/180
        ks1[6]=90*np.pi/180
        switch_params(1)
        switch_cparams(1)
        alphaentry.config(state="disabled")
        gammaentry.config(state="disabled")
    elif g in range(16,75):
        #Ortorrómbicos
        ks1[4]=90*np.pi/180
        ks1[5]=90*np.pi/180
        ks1[6]=90*np.pi/180
        switch_params(1)
        switch_cparams(1)
        alphaentry.config(state="disabled")
        betaentry.config(state="disabled")
        gammaentry.config(state="disabled")
    elif g in range(75,143):
        #Tetragonales
        ks1[4]=90*np.pi/180
        ks1[5]=90*np.pi/180
        ks1[6]=90*np.pi/180
        switch_params(0)
        switch_cparams(1)
        aentry.config(state="normal")
        centry.config(state="normal")
    elif g in range(143,195):
        #Hexagonales y trigonales
        ks1[4]=90*np.pi/180
        ks1[5]=90*np.pi/180
        ks1[6]=120*np.pi/180
        switch_params(0)
        switch_cparams(1)
        aentry.config(state="normal")
        centry.config(state="normal")
    elif g in range(195,231):
        #Cúbicos
        ks1[4]=90*np.pi/180
        ks1[5]=90*np.pi/180
        ks1[6]=90*np.pi/180
        switch_params(0)
        switch_cparams(1)
        aentry.config(state="normal")
    else:
        print("Grupo no reconocido. Por favor ingrese un número entero entre 1 y 230")
    rmr()
        
groupbutton=tk.Button(root,text="Confirmar grupo espacial",command=groupbut)


cel=[]
cel1=[]
def celadd():
    global cel
    global cel1
    try:
        a=[Wyckoff.elementos[Wyckoff.Standard_ion(Zentry.get())][0],float(eval(xentry.get())),float(eval(yentry.get())),float(eval(zentry.get())),Wyckoff.elementos[Wyckoff.Standard_ion(Zentry.get())][1],float(eval(uentry.get()))]
        cel.append(a)
    except:
        try:
            a=[Wyckoff.elementos[Wyckoff.Standard_ion(Zentry.get(),mute=True)][0],float(eval(xentry.get())),float(eval(yentry.get())),float(eval(zentry.get())),Wyckoff.elementos[Wyckoff.Standard_ion(Zentry.get(),mute=True)][1],1]
            cel.append(a)
        except:
            print("Entrada inválida")
    cel1=Wyckoff.Wyckoff(g,cel)
    cellabel2.config(text=cel1)
addbutton=tk.Button(root,text="Añadir",command=celadd)
def celrm():
    global cel
    global cel1
    try:
        cel[-1]
        cel.pop(-1)
    except:
        print("Nada que remover")
    cel1=Wyckoff.Wyckoff(g,cel)
    cellabel2.config(text=cel1)
rmbutton=tk.Button(root,text="Quitar último",command=celrm)
def rmr():
    global cel
    global cel1
    cel=[]
    cel1=[]
    cellabel2.config(text=cel1)
rmrbutton=tk.Button(root,text="Limpiar celda",command=rmr)

def calc():
    ks=ks1
    cel2=cel1
    if g in range(1,3):
        #Monoclínicos
        ks[1]=eval(aentry.get())
        ks[2]=eval(bentry.get())
        ks[3]=eval(centry.get())
        ks[4]=eval(alphaentry.get())*np.pi/180
        ks[5]=eval(betaentry.get())*np.pi/180
        ks[6]=eval(gammaentry.get())*np.pi/180
    elif g in range(3,16):
        #Triclínicos
        ks[1]=eval(aentry.get())
        ks[2]=eval(bentry.get())
        ks[3]=eval(centry.get())
        ks[5]=eval(betaentry.get())*np.pi/180
    elif g in range(16,75):
        #Ortorrómbicos
        ks[1]=eval(aentry.get())
        ks[2]=eval(bentry.get())
        ks[3]=eval(centry.get())
    elif g in range(75,143):
        #Tetragonales
        ks[1]=eval(aentry.get())
        ks[2]=eval(aentry.get())
        ks[3]=eval(centry.get())
    elif g in range(143,195):
        #Hexagonales y trigonales
        ks[1]=eval(aentry.get())
        ks[2]=eval(aentry.get())
        ks[3]=eval(centry.get())
    elif g in range(195,231):
        #Cúbicos
        ks[1]=eval(aentry.get())
        ks[2]=eval(aentry.get())
        ks[3]=eval(aentry.get())
    for i in range(0,len(ks)):
        ks[i]=float(ks[i])
    prov=ref.MinInfoPat(ks,cel2)
    print(prov)
    patron=ref.Intensidad(prov, 1, np.linspace(ref.minang-1,ref.maxang+1,1000))
    patron=patron*100/np.max(patron)
    fig,ax=plt.subplots()
    ax.plot(np.linspace(ref.minang-1,ref.maxang+1,1000),patron)
    plt.show()
calcbutton=tk.Button(root,text="Calcular",command=calc)

def Opencif():
    global cel
    global cel1
    global g
    filename = filedialog.askopenfilename()
    print('Archivo cargado', filename)
    provisional=Wyckoff.semi_readCIF(filename)
    #Todas estas líneas tienen sentido al ver la salida de Wyckoff.semi_readCIF.
    g=provisional[0]
    ks=provisional[1]
    cel=provisional[2]
    grupo=provisional[3]
    cel1=Wyckoff.custom_Wyckoff(grupo, cel)
    cellabel2.config(text=cel1)
    prov=ref.MinInfoPat(ks,cel1)
    print(prov)
    patron=ref.Intensidad(prov, 1, np.linspace(ref.minang-1,ref.maxang+1,1000))
    patron=patron*100/np.max(patron)
    fig,ax=plt.subplots()
    ax.plot(np.linspace(ref.minang-1,ref.maxang+1,1000),patron)
    plt.show()
    
cifbutton = tk.Button(root, text="Leer desde archivo cif", command=Opencif)


cifbutton.grid(row=2,column=4)
groupbutton.grid(row=2,column=3)
addbutton.grid(row=8,column=0)
rmbutton.grid(row=8, column=2)
rmrbutton.grid(row=8,column=3)
calcbutton.grid(row=8,column=5)






root.mainloop()