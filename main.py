from testdata import testdatarr
import testdata as t
import numpy as np
import matplotlib as plt
from refinamiento2 import *

lam=1.54060
la2=1.5443
kb=[la2,0.5]



r=[]
As=[]
Xs=[]
k=[]
k2=[]
u=[]

prov=Wyckoff.semi_readCIF("Ba2SbSmO6.cif")

Ba2SbSmO6= estructura(prov[1],prov[2],prov[0],prov[3])
#Ba2SbSmO6.celda[3,1] = 0.225
Ba2SbSmO6.FullCell=Wyckoff.custom_Wyckoff(Ba2SbSmO6.grupo, Ba2SbSmO6.celda)
Ba2SbSmO6.ks[0]=2.05e-4
Ba2SbSmO6.ks[[1,2,3]]=8.48
#Ba2SbSmO6.moveatom(19, +0.1)
Ba2SbSmO6.refinableCel=np.array([[0,0,0,0,0,0],
                            [0,0,0,0,0,0],
                            [0,0,0,0,0,0],
                            [0,1,0,0,0,0]],dtype=bool)
Ba2SbSmO6.refinableks = np.array([1,1,0,0,0,0,0])
Ba2SbSmO6.UVWg[0]= 0.025
Ba2SbSmO6.UVWg[1]= 0
Ba2SbSmO6.UVWg[2]= 0.05
Ba2SbSmO6.refinableUVWg=np.array([0,0,1,0])
ASCBa=np.genfromtxt("Ba2SbSmO6.ASC",delimiter=" ")
res,ax2=plt.subplots()

#ASCBa=ASCBa[np.where(ASCBa[:,0] < 84)]

ASCBa = datos(ASCBa)
ASCBa.angs -= 0.03
ax2.plot(ASCBa.angs,1000 + ASCBa.ints,color="black",label="Experimental")
ax2.scatter(ASCBa.angs, 1000 + ASCBa.bgints + intkb(Ba2SbSmO6,ASCBa.angs,kb),color="blue",label="inicial",s=4)



suma = np.sum(ASCBa.ints)

NumIters=20

# #IMPORTANTE: ACELERAR CUSTOM_WYCKOFF Y MOVEATOM PARA TENER UN IDREF FUNCIONAL
# t1=time.time()
# for i in range(0,100):
#     prov = Wyckoff.custom_Wyckoff(Ba2SbSmO6.grupo, Ba2SbSmO6.celda)
# t2=time.time()
# print(t2-t1)


# Ba2SbSmO6.refinableCel=None
# t1= time.time()
# Ba2SbSmO6.idref()
# t2=time.time()
# print(t2-t1)
# print(Ba2SbSmO6.refinableCel)



r.append(np.sum(np.abs(ASCBa.ints - ASCBa.bgints - intkb(Ba2SbSmO6,ASCBa.angs,kb,CheckAngs=True)))/suma)
As.append(Ba2SbSmO6.ks[1])
Xs.append(Ba2SbSmO6.celda[3,1])
k.append(Ba2SbSmO6.ks[0])
u.append(Ba2SbSmO6.UVWg[2])

Ba2SbSmO6.refinableCel=np.array([[0,0,0,0,0,0],
                            [0,0,0,0,0,0],
                            [0,0,0,0,0,0],
                            [0,0,0,0,0,0]],dtype=bool)
Ba2SbSmO6.refinableks = np.array([1,1,0,0,0,0,0])
Ba2SbSmO6.refinableUVWg=np.array([0,0,1,0])

for i in range(0,NumIters):
      print(i)
      vector = Ba2SbSmO6.iteration(ASCBa,kb=kb)
      Ba2SbSmO6.ks[0] -= vector[0,0]
      Ba2SbSmO6.ks[[1,2,3]] += vector[1,0]
      Ba2SbSmO6.UVWg[2] = np.abs(Ba2SbSmO6.UVWg[2] + vector[2,0])
#      Ba2SbSmO6.moveatom(19,vector[3,0])
      Ba2SbSmO6.FullCell= Wyckoff.custom_Wyckoff(Ba2SbSmO6.grupo, Ba2SbSmO6.celda)
#      Ba2SbSmO6.FullCell = Wyckoff.custom_Wyckoff(Ba2SbSmO6.grupo, Ba2SbSmO6.celda)
      r.append(np.sum(np.abs(ASCBa.ints - ASCBa.bgints - intkb(Ba2SbSmO6,ASCBa.angs,kb,CheckAngs=True)))/suma)
      As.append(Ba2SbSmO6.ks[1])
      Xs.append(Ba2SbSmO6.celda[3,1])
      k.append(Ba2SbSmO6.ks[0])
      u.append(Ba2SbSmO6.UVWg[2])


#vector = Ba2SbSmO6.iteration(ASCBa,kb=kb,flag=True)





if NumIters !=0:
    print(Ba2SbSmO6.ks)
    print(Ba2SbSmO6.celda)
    # # plt.legend()
    # # # plt.title(r'$w_i = I^{-1}$')
    fig,ax=plt.subplots()
    # plt.title(r'$w_i = I^{-1}$')
    ax.plot(np.arange(len(r)),np.array(r),color="blue")
    plt.xlabel("Iteraciones")
    plt.ylabel(r'$R_p$')
    ax.set_xticks(np.arange(0,NumIters+1,2))
    
    fig,ax=plt.subplots()
    ax.plot(np.arange(len(As)),np.array(As),color="black")
    #plt.title(r'$w_i = I^{-1}$')
    plt.xlabel("Iteraciones")
    plt.ylabel(r'$a (\AA)$')
    ax.plot(np.array([0,NumIters+1]),np.array([8.503,8.503]),linestyle="--",color="red")
    ax.set_xticks(np.arange(0,NumIters,2))

    
    fig,ax=plt.subplots()
    ax.plot(np.arange(len(k)),np.array(k))
    #plt.title(r'$w_i = I^{-1}$')
    plt.xlabel("Iteraciones")
    plt.ylabel("k")
    
    fig,ax=plt.subplots()
    ax.plot(np.arange(len(u)),np.array(u))
    #plt.title(r'$w_i = I^{-1}$')
    plt.xlabel("Iteraciones")
    plt.ylabel("W")

ax2.scatter(ASCBa.angs, 1000+ ASCBa.bgints + intkb(Ba2SbSmO6,ASCBa.angs,kb),color="red",label="final",s=4)

ax2.scatter(ASCBa.angs, ASCBa.ints - ASCBa.bgints - intkb(Ba2SbSmO6,ASCBa.angs,kb),color="green",s=4,label="diferencia")
# # #fig, ax=plt.subplots()
# # # ax.plot(ASCBa.angs, ASCBa.ints - ASCBa.bgints - Intensidad([Ba2SbSmO6,Copia], ASCBa.angs))   
ax2.legend()

