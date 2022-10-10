import numpy as np
import sympy as sp
import sympy
import math
import time

t0 = time.time()

#DATOS PREGUNTA

b_12 = 1/(0.2)
b_14 = 1/(0.2)
b_15 = 1/(0.3)
b_23 = 1/(0.25)
b_24 = 1/(0.1)
b_25 = 1/(0.3)
b_26 = 1/(0.2)
b_35 = 1/(0.24)
b_36 = 1/(0.1)
b_45 = 1/(0.4)
b_56 = 1/(0.3)

B = np.array([[b_12+b_24+b_25+b_26+b_23, -b_23, -b_24, -b_25, -b_26],
    [-b_23, b_23+b_35+b_36, 0, -b_35, -b_36],
    [-b_24, 0, b_14+b_24+b_45, -b_45, 0],
    [-b_25, -b_35, -b_45, b_25+b_35+b_45+b_56, -b_56],
    [-b_26, -b_36, 0, -b_56, b_26+b_36+b_56]])


P = np.array([0.5, 0.6, -0.7, -0.7, -0.7])

Binv = np.linalg.inv(B)

A = np.dot(Binv,P)

tf = time.time()

print("APROXIMACIÓN DC")
for i in range (5):
    if i==0:
        print("Voltaje nodo", 1, ": 1 < 0")
    print("Voltaje nodo", i+2, ": 1 <", np.round(A[i]*180/np.pi,4))
print("Calculado en",np.round(tf - t0,8),"segundos")
