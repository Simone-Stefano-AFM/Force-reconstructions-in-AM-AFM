# -*- coding: utf-8 -*-
"""
@author: Simone Benaglia
@author: Stefano Chiodini

If you use the following code please cite: "Quantification of solvation forces with amplitude modulation AFM"
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
from scipy.signal import savgol_filter
# import scipy.signal as signal



# data from "Holscher_FRM.py file"
k_Holscher = np.loadtxt("k_Holscher.txt")[1, :]
Dmin = np.loadtxt("k_Holscher.txt")[0, :]
amplitude = np.loadtxt("Amp.txt")[1, :]
dforce = np.loadtxt("force.txt")[1, :]
D1 = np.loadtxt("force.txt")[0, :]


# interpolation of data
n = 7000 # at least 1000. For high values better 64 bit than 32 bit versions of Python

Dmin_interp = np.linspace(np.min(Dmin), np.max(Dmin), n)

interp_k = si.interp1d(Dmin, k_Holscher, kind = "linear")
k_interp = interp_k(Dmin_interp)

interp_A = si.interp1d(Dmin, amplitude, kind = "linear")
A_interp = interp_A(Dmin_interp) 

a1 = Dmin_interp
b1 = a1 + 2.*A_interp

# initialization of matrix elements to zero
mat = np.zeros([n, n])

# bin of Dmin to compute matrix
step = (np.max(a1) - np.min(a1))/(n-1)


# choice of the kernel: the real kernel (line 51) or the Holscher approximated one (line 52) 
# see paper Holscher, APL 86, 123109, 2006

def ker(x, t, A):
    
    return ((t-x-A)/(A**2-(t-x-A)**2)**0.5)
#   return 1.0/(t - x)**0.5
    

#Construction of the matrix coefficients (correct)

for i in range (0, n):
    
    print i
    
    if b1[i] < np.max(a1):

        jmax = np.max(np.where(a1 < b1[i])[0]) 
        
        if i == 0:
            
            mat[0][0] = 0.5*step*ker(a1[0], a1[0] + step/2, A_interp[0]) 
            
            mat[0][jmax] = 0.5*step*ker(a1[0], a1[jmax] + step/2, A_interp[0])
            
            for j in range (1, jmax):
                
                mat[0][j] = step*ker(a1[0], a1[j] + step/2, A_interp[0])
        
        else: 
            
            for j in range (0, i): 
                
                mat[i][j] = 0.
    
            for j in range (i, jmax):
                
                if j == i:
                    
                   mat[i][j] = 0.5*step*ker(a1[i], a1[j] + step/2, A_interp[i]) 
                   
                elif j > i:
    
                    mat[i][j] = step*ker(a1[i], a1[j] + step/2, A_interp[i])
                
            
            for j in range (jmax, n):
                
                if j == jmax:
                    
                    mat[i][j] = 0.5*step*ker(a1[i], a1[j] + step/2, A_interp[i])
                    
                else:    
            
                    mat[i][j] = 0.
        

    elif b1[i] > np.max(a1):

        for j in range (0, i):

            if i == 0:
                
                mat[i][j] = 0.5*step*ker(a1[0], a1[0] + step/2, A_interp[0]) 
        
            else: 
                
                mat[i][j] = 0.       
        
        
        for j in range(i, n):
            
            if j == n-1:
                
                mat[i][j] = 0.5*step*ker(a1[i], a1[j] + step/2, A_interp[i])

            elif j == i and i != 0 : 
                
                mat[i][j] = 0.5*step*ker(a1[i], a1[j] + step/2, A_interp[i])
                
            elif j > i and j != n-1:
                
                mat[i][j] = step*ker(a1[i], a1[j] + step/2, A_interp[i])

          
mat[np.isnan(mat)] = 0.

# invert matrix
inv_mat = np.linalg.inv(mat)

# find the solution (the force)
force_matrix = inv_mat.dot(-np.pi*np.sqrt(A_interp/2)*k_interp) 

smooth_force_matrix = savgol_filter(force_matrix, 251, 3) # to be checked every time

# final plots
plt.plot(a1, smooth_force_matrix, label = 'Matrix', color = 'blue', linewidth = 3)
plt.plot(D1, dforce, label = 'Force', color = 'purple', linewidth = 3)
plt.legend()



# -------------------------------------------------------------
# same inversion for the Holscher approximated kernel (line 52).
# if used comment line 51 and from 132 to 145
# -------------------------------------------------------------

#force_approx = inv_mat.dot(np.pi*I_interp) 
#smooth_force_approx = savgol_filter(force_approx, 51, 3)
#D = np.linspace (np.min(a1), np.max(a1), n)
#plt.plot(D, smooth_force_approx, label = 'Approx. Holscher')
 
