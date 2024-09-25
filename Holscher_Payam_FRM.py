# -*- coding: utf-8 -*-
"""
@author: Simone Benaglia
@author: Stefano Chiodini

If you use the following code please cite: "Quantification of solvation forces with AFM"
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
import scipy.signal as signal
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline
from numpy import inf
from scipy.signal import savgol_filter


#-----------------------------------------------------------------------------
# Code required to obtain Holscher and Payam force
# reconstruction from amplitude-, phase- and deflection-distance experimental 
# AFM data
#------------------------------------------------------------------------------

# AFM parameters
k = 3.0 # spring constant
Q = 6.0 # Q-factor
f = 800000. # resonance frequency
w = 2.0*np.pi*f # angular frequency

# experimental data 
# i=0 close, i=N far
amplitude = np.loadtxt("Amplitude.txt")
phase_grad = np.loadtxt("Phase.txt")
deflection = np.loadtxt("deflection.txt") # if deflection is required use this. If not comment and modify line 85

# reference force. If not needed, comment
D1 = np.loadtxt("F-d.txt")[:,0]*1.0e-9
f_d = np.loadtxt("F-d.txt")[:,1]*1.0e-12

# interpolation of reference force 
D = np.linspace(np.min(D1), np.max(D1), 1000)
interp_f_d = si.interp1d(D1, f_d, kind = "slinear")
f_d = interp_f_d(D)


# data
zc1 = np.flipud(amplitude[:, 0]*1.0e-9)

Amp = np.flipud(amplitude[:, 1]*1.0e-9)
A0 = Amp[-1]

# for phase data we follow the Asylum Research convention where free phase = 90° and
# the attractive (repulsive) regime corresponds to phase > 90° (<90°)
phi_grad = np.flipud(phase_grad[:, 1]) 
phi_grad = phi_grad + 90.0 - phi_grad[-1]
phi_rad1 = phi_grad*np.pi/180.0

defl1 = np.flipud(deflection[:, 1]*1.0e-9)


# interpolation of data
N = 1000
zc2 = np.linspace (np.min(zc1), np.max(zc1), N)

interpolation1 = si.interp1d(zc1, Amp, kind = "slinear")
A = interpolation1(zc2)

interpolation2 = si.interp1d(zc1, phi_rad1, kind = "slinear")
phi_rad = interpolation2(zc2)

interpolation4 = si.interp1d(zc1, defl1, kind = "slinear")
defl = interpolation4(zc2)

# resetting of name zc
zc = zc2


# definition of variables
k_Holscher = np.zeros(N)
integrand = np.zeros(N)
Dmin = np.zeros(N)
Dmax = np.zeros(N)

Dmin = zc - A + defl
Dmax = Dmin + 2.0*A

# setting amplitude and phase as function of Dmin
interpolation_A = si.interp1d(Dmin, A, kind = "slinear")
interpolation_phi = si.interp1d(Dmin, phi_rad, kind = "slinear")

# definition of the k function of Holscher, see paper Holscher, APL 86, 123109, 2006
k_Holscher = (k*A0/Q)*np.sqrt(interpolation_A(Dmin)/2.0)*np.cos(interpolation_phi(Dmin))

# definition of variables
integral = np.zeros(N) 
error = np.zeros(N)
integral_payam1 = np.zeros(N) 
error_payam1 = np.zeros(N)
integral_payam2 = np.zeros(N) 
error_payam2 = np.zeros(N)

index1 = 0
index2 = N - 2


for i in range (index1, index2): 
    
    print i
    
    
#   ******HOLSCHER****** see paper Holscher, APL 86, 123109, 2006
    integrand1 = (k_Holscher[(i+1):N]/(np.sqrt(Dmin[i+1:N] - Dmin[i]))) 
    integrand1[integrand1 == -inf] = 0
    interpolation_H = si.interp1d(Dmin[i+1:N], integrand1, kind = "slinear")
    
    
#   ******PAYAM******** see paper Payam, Nanotechnology 26, 185706, 2015
    X = ((A0)/(2.0*Q*interpolation_A(Dmin[i+1:N])))*np.cos(interpolation_phi(Dmin[i+1:N]))  
    Y = ((A0)/(2.0*Q*interpolation_A(Dmin[i+1:N])*w))*np.sin(interpolation_phi(Dmin[i+1:N])) - ((1.0)/(2.0*Q*w))  
    payam1 = 2.0*k*X 
    interp_payam1 = si.interp1d(Dmin[i+1:N], payam1, kind = "slinear") 
    payam2 = -2.0*k*((X*(interpolation_A(Dmin[i+1:N]))**1.5)/(np.sqrt(2.0*(Dmin[i+1:N] - Dmin[i]))))
    interp_payam2 = si.interp1d(Dmin[i+1:N], payam2, kind = "slinear")        
    omega = np.sqrt(1.0 + ((A0/Q)/(interpolation_A(Dmin[i+1:N])))*np.cos(interpolation_phi(Dmin[i+1:N]))) - 1.0
    diff_omega = np.diff(omega)/np.diff(Dmin[i+1:N]) 
    diff_omega = np.append(diff_omega, diff_omega[-1])

# ******COMMON TO BOTH****** 
# consider that Holscher integrates between Dmin and Dmax, Payam between Dmin and infinity
    if (Dmax[i+1] <= Dmin[-1]):
        integral[i], error[i] = integrate.quad(interpolation_H, Dmin[i+1], Dmax[i+1])
        integral_payam1[i], error_payam1[i] = integrate.quad(interp_payam1, Dmin[i+1], Dmin[-1])
#        integral_payam1[i], error_payam1[i] = integrate.quad(interp_payam1, Dmin[i+1], Dmax[i+1])
        integral_payam2[i], error_payam2[i] = integrate.quad(interp_payam2, Dmin[i+1], Dmin[-1])
#        integral_payam2[i], error_payam2[i] = integrate.quad(interp_payam2, Dmin[i+1], Dmax[i+1])
    
        
    else:
        integral[i], error[i] = integrate.quad(interpolation_H, Dmin[i+1], Dmin[-1])
        integral_payam1[i], error_payam1[i] = integrate.quad(interp_payam1, Dmin[i+1], Dmin[-1])
        integral_payam2[i], error_payam2[i] = integrate.quad(interp_payam2, Dmin[i+1], Dmin[-1]) 

                                     
integral[np.isnan(integral)] = 0.
smooth_integral = savgol_filter(integral, 7, 1) # to be checked everytime

payam2 = np.diff(integral_payam2)/np.diff(Dmin)

N1  = 2    # Filter order
Wn = 0.09 # Cutoff frequency (to be checked everytime)
b, a = signal.butter(N1, Wn, output='ba')
smooth_payam2 = signal.filtfilt(b, a, payam2)
smooth_payam2 = np.asarray(smooth_payam2)

smooth_payam2 = np.append(smooth_payam2, smooth_payam2[-1])

payam = integral_payam1 + smooth_payam2


# "InterpolationUnivariateSpline" requires a monotonically increasing x-function
# for this reason, sometimes, better to smooth Dmin, in order to  remove the non-increasing parts
smooth_Dmin = savgol_filter(Dmin, 127, 1)

interpolazione = InterpolatedUnivariateSpline(smooth_Dmin, smooth_integral, k = 1)
force_smooth = interpolazione.derivative()


# final plotting
plt.plot(smooth_Dmin, -force_smooth(smooth_Dmin), color = 'red', label = "Holscher", linewidth = 3)
plt.plot(D, f_d, label = "Force", color = 'purple', linewidth = 3)
plt.plot(Dmin, payam, color = 'blue', label = "Payam", linewidth = 3)

   