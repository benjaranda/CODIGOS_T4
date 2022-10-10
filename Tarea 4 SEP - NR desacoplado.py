import numpy as np
import sympy as sp
import sympy
import math
import cmath
import time

t0 = time.time()

#DATOS PREGUNTA
y_12 = 1/(0.12+0.2j)
y_14 = 1/(0.05+0.2j)
y_15 = 1/(0.08+0.3j)
y_23 = 1/(0.06+0.25j)
y_24 = 1/(0.05+0.1j)
y_25 = 1/(0.1+0.3j)
y_26 = 1/(0.09+0.2j)
y_35 = 1/(0.12+0.24j)
y_36 = 1/(0.03+0.1j)
y_45 = 1/(0.2 + 0.4j)
y_56 = 1/(0.14 + 0.3j)

Y = np.array([[y_12+y_14+y_15,-y_12,0,-y_14,-y_15,0],
     [-y_12,y_12+y_24+y_25+y_26+y_23,-y_23,-y_24,-y_25,-y_26],
     [0,-y_23,y_36+y_35+y_23,0,-y_35,-y_36],
     [-y_14,-y_24,0,y_14+y_24+y_45,-y_45,0],
     [-y_15,-y_25,-y_35,-y_45,y_15+y_25+y_35+y_45+y_56,-y_56],
     [0,-y_26,-y_36,0,-y_56,y_26+y_36+y_56]])

P = np.array([0.5, 0.6, -0.7, -0.7, -0.7]) #2,3,4,5,6
Q = np.array([-0.7, -0.7, -0.7]) #4,5,6

#Estado inicial
Inc0 = [0, 0, 0, 0, 0, 1, 1, 1] #Incógnitas a2, a3, a4, a5, a6, v4, v5, v6
M = np.array([1.05, 1.05, 1.07, 1, 1, 1]) #magnitud voltaje inicial
A = np.array([0, 0, 0, 0, 0, 0]) #fase voltaje inicial

def VaS(Y,Mi,Ai): #entra M y A y sale S para todos los nodos
    V = [] #voltaje 6 nodos
    for i in range (6):
        V.append(cmath.rect(Mi[i],Ai[i])) 
    S = []
    for i in range (6):
        Ys = []
        for j in range(6):
            Ys.append(Y[i][j]*V[j])
            Yss = sum(Ys)
        S.append(V[i] * Yss.conjugate())
    S = np.array(S)
    return S #S para todo nodo

def PQ(s,i,o,k): #s= P o Q, i= nodo i+1, o= m o a, k nodo k+1
    Pa = 0
    Qa = 0
    for j in range(6): #j va de 0 a 5, es decir, nodo 1 a 6
        if o == "a": #deriva respecto ángulo k
            v = A[k] 
            if j != k:
                Pa += M[i]*M[j]*abs(Y[i][j])*sympy.cos(cmath.phase(Y[i][j])-A[i]+A[j])
                Qa -= M[i]*M[j]*abs(Y[i][j])*sympy.sin(cmath.phase(Y[i][j])-A[i]+A[j])
            else:
                x = sp.Symbol('x')
                Pa += M[i]*M[j]*abs(Y[i][j])*sympy.cos(cmath.phase(Y[i][j]) -A[i] + x)
                Qa -= M[i]*M[j]*abs(Y[i][j])*sympy.sin(cmath.phase(Y[i][j]) -A[i] + x)
        if o == "v": #deriva respecto magnitud k
            v = M[k]
            if j != k:
                Pa += M[i]*M[j]*abs(Y[i][j])*sympy.cos(cmath.phase(Y[i][j])-A[i]+A[j])
                Qa -= M[i]*M[j]*abs(Y[i][j])*sympy.sin(cmath.phase(Y[i][j])-A[i]+A[j])
            else:
                x = sp.Symbol('x')
                Pa += M[i]*x*abs(Y[i][j])*sympy.cos(cmath.phase(Y[i][j])-A[i]+ A[j])
                Qa -= M[i]*x*abs(Y[i][j])*sympy.sin(cmath.phase(Y[i][j])-A[i]+ A[j])
    if s == "P":
        dP = sp.diff(Pa, x).subs(x, v) #derivada de P respecto a x evaluada en v
        ret = dP
    else:
        dQ = sp.diff(Qa, x).subs(x, v) #derivada de Q respecto a x evaluada en v
        ret = dQ
    return ret

def R(S): #vector residuos (variación P y variación Q)
    VP = P - np.array([S.real[1], S.real[2], S.real[3], S.real[4], S.real[5]]) #2,3,4,5,6
    VQ = Q - np.array([S.imag[3], S.imag[4], S.imag[5]]) #4,5,6
    R = np.concatenate((VP, VQ))
    return R

def J(): #matriz jacobiana
    J = np.array([[PQ("P",1,"a",1), PQ("P",1,"a",2), PQ("P",1,"a",3), PQ("P",1,"a",4), PQ("P",1,"a",5), 0, 0, 0],
        [PQ("P",2,"a",1), PQ("P",2,"a",2), PQ("P",2,"a",3), PQ("P",2,"a",4), PQ("P",2,"a",5), 0, 0, 0],
        [PQ("P",3,"a",1), PQ("P",3,"a",2), PQ("P",3,"a",3), PQ("P",3,"a",4), PQ("P",3,"a",5), 0, 0, 0],
        [PQ("P",4,"a",1), PQ("P",4,"a",2), PQ("P",4,"a",3), PQ("P",4,"a",4), PQ("P",4,"a",5), 0, 0, 0],
        [PQ("P",5,"a",1), PQ("P",5,"a",2), PQ("P",5,"a",3), PQ("P",5,"a",4), PQ("P",5,"a",5), 0, 0, 0],
        [0, 0, 0, 0, 0, PQ("Q",3,"v",3), PQ("Q",3,"v",4), PQ("Q",3,"v",5)],
        [0, 0, 0, 0, 0, PQ("Q",4,"v",3), PQ("Q",4,"v",4), PQ("Q",4,"v",5)],
        [0, 0, 0, 0, 0, PQ("Q",5,"v",3), PQ("Q",5,"v",4), PQ("Q",5,"v",5)]
        ], dtype = float)
    return J

#Se realiza la primera iteración, donde se calcula la matriz jacobiana en 0 para usarla en las siguientes
S0 = VaS(Y, M, A)
J0 = J()
Inc = np.dot(np.linalg.inv(J0), R(S0)) + Inc0 #a2, a3, a4, a5, a6, v4, v5, v6 incógnitas
M =  np.array([1.05, 1.05, 1.07, Inc[5], Inc[6], Inc[7]])
A =  np.array([0, Inc[0], Inc[1], Inc[2], Inc[3], Inc[4]])
S = VaS(Y, M, A)

i=1
while i >= 0: #OTRAS ITERACIONES CON J0 FIJO
    i+=1
    Inc = np.dot(np.linalg.inv(J()), R(S)) + Inc #a2, a3, a4, a5, a6, v4, v5, v6 incógnitas #PARA DESHONESTO IRÍA J0 Y NO J()
    M =  np.array([1.05, 1.05, 1.07, Inc[5], Inc[6], Inc[7]])
    A =  np.array([0, Inc[0], Inc[1], Inc[2], Inc[3], Inc[4]])
    S = VaS(Y, M, A)

    t = 0
    for k in range(5): #condición de salida
        if R(S)[k] >= 0.00001:
            t = 1
    if t != 1:
        ni = i
        i = -1  

tf = time.time()

print("DESACOPLADO")

for i in range(6):
    print("Voltaje nodo", i+1, ":", np.round(M[i],4), "<", np.round(A[i]*180/np.pi,4))

print("Calculado en",np.round(tf - t0,4),"segundos")
print("Total de iteraciones:", ni)
print(np.round((tf - t0)/ni,4), "segundos por iteración")
