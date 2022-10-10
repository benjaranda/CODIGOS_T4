import math
import cmath
import numpy as np
import time

t0 = time.time()

########## Definición tipo de barra ##########
########## k Type Pgk  Qgk Pdk   Qdk   V0#####
Barras = [[1, 0 ,  0,   0,   0,    0,   1.05],
          [2, 1 , 0.5,  0,   0,    0,   1.05],
          [3, 1,  0.5,  0,   0,    0,   1.07],
          [4, 2,  0,    0,  0.7,  0.7,    1],
          [5, 2,  0,    0,  0.7,  0.7,    1],
          [6, 2,  0,    0,  0.7,  0.7,    1]]

##################################
#Generación Matriz de admitancias#
##################################
y_12 = 1/(0.12+0.2j); y_14 = 1/(0.05+0.2j); y_15 = 1/(0.08+0.3j)
y_23 = 1/(0.06+0.25j); y_24 = 1/(0.05+0.1j); y_25 = 1/(0.1+0.3j); y_26 = 1/(0.09+0.2j)
y_35 = 1/(0.12+0.24j); y_36 = 1/(0.03+0.1j)
y_45 = 1/(0.2 + 0.4j)
y_56 = 1/(0.14 + 0.3j)

Y = [[y_12+y_14+y_15,-y_12,0,-y_14,-y_15,0],
     [-y_12,y_12+y_24+y_25+y_26+y_23,-y_23,-y_24,-y_25,-y_26],
     [0,-y_23,y_36+y_35+y_23,0,-y_35,-y_36],
     [-y_14,-y_24,0,y_14+y_24+y_45,-y_45,0],
     [-y_15,-y_25,-y_35,-y_45,y_15+y_25+y_35+y_45+y_56,-y_56],
     [0,-y_26,-y_36,0,-y_56,y_26+y_36+y_56]]

NB = 6 #Numero de barras

##################################
######Método Gauss-Seidel#########
##################################

error = (0.0001)
Iteraciones = 300000
x = 6
V_0 = [fila[x] for fila in Barras]; Vt = V_0
V_k = [1.05, 1.05, 1.07,0,0,0];
for t in range(Iteraciones):
    for k in range(NB):
        if Barras[k][1] == 0:
           Vt[k] = Barras[k][6]
           V_0[k] = Vt[k]

        elif Barras[k][1] == 1:
            Sdk = Barras[k][4]+Barras[k][5]*1j
            Qgk = -1*(Sdk.conjugate()+sum(np.multiply((np.multiply(Vt[k].conjugate(),Y[k])),Vt))).imag
            Sgk = Barras[k][2]+1j*Qgk
            Sk=Sgk-Sdk
            V_pv = (1/(Y[k][k]))*((Sk.conjugate())/Vt[k].conjugate()-sum((np.multiply(Y[k],Vt)))+Y[k][k]*Vt[k])
            V_angle = math.atan(V_pv.imag/V_pv.real)
            Vt[k] = abs(Vt[k])*(math.cos(V_angle)+1j*math.sin(V_angle))
            V_0[k] = Vt[k]

        else:
            Sdk = Barras[k][4] + Barras[k][5] * 1j
            Sgk = Barras[k][2] + 1j * Barras[k][3]
            Sk = Sgk - Sdk
            Vt[k] = (1 / (Y[k][k])) * ((Sk.conjugate()) / Vt[k].conjugate() - sum((np.multiply(Y[k], Vt))) + Y[k][k] * Vt[k])
            V_0[k] = Vt[k]

    if max(abs(abs(np.array(Vt))-abs(np.array(V_k)))) < error:
        break
    if max(abs(abs(np.array(Vt)) - abs(np.array(V_k)))) > error:
        V_k = [Vt[0],Vt[1],Vt[2],Vt[3],Vt[4],Vt[5]]
        V_0 = Vt

tf = time.time()
print("GAUSS - SEIDEL")

for i in range(6):
    print("Voltaje nodo", i+1, ":", np.round(abs(V_0[i]), 4), "<", np.round((180*math.atan(V_0[i].imag/V_0[i].real)/(math.pi)),4))

print("Calculado en",np.round(tf - t0,4),"segundos")
print("Iteraciones: " + str(t+1))
