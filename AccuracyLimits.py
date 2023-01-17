################################################################################################################
# Name: AccuracyLimits.py
#
# Purpose1: Plot Relative error of \gamma vs. \delta\Theta
# Purpose2: Animate Diffraction against aperture position
#			
# Things to fix: Comments
#
#################################################################################################################

import matplotlib.pyplot as plt

import UndulatorDiffractionOneDimensional_more_U1U2 as UnDi
import InterferenceFunction_more as IFunc
import AmplitudeFormula as AUndu

from numpy import transpose as nptp
import numpy as np

import cv2
import os , sys
import re

import csv

from natsort import natsorted
from scipy.optimize import curve_fit
from lmfit import minimize, Parameters, report_fit, Minimizer
from matplotlib.animation import FuncAnimation, PillowWriter 
import random

import imageio

ratio = 1000. / 606#606. #how many pixel are aperture
UnDi.apy = .004


zetaprime2 = np.linspace(-ratio*UnDi.apy,ratio*UnDi.apy,1000)
fixedsize = .245


NumOfStripes = 1
NumOfSeries = 100#825   
LengthOfInterval = 25#825

dispersion = 0.004808759252942803e-9 #Former 0.0048085643228767475e-9
pixelof404nm = 1315.3024745554976 #Former: 1243.8551272671261
CalibWl = 404.6565e-9

#(407.7837-404.6565)

chosenpx = 1300

global listOfboots
listOfboots = []

Yscreenfit = .3e-4#-.5e-4
LambdaRfit = 400e-9

def diffraction_dataset(params,i,yy):
    Amplfit = params['Amplfit_%i' % (i)].value
    Ytran   = params['Ytran_%i' % (i)].value
    gammafit  = params['gamma_%i' % (i)].value
    delta0fit  = params['delta0_%i' % (i)].value
    bfit  = params['bfit_%i' % (i)].value
    zpfit  = params['zpfit_%i' % (i)].value

    #Yscreenfit  = params['Yscreenfit_%i' % (i+1)].value
    #LambdaRfit  = params['LambdaRfit_%i' % (i+1)].value
    #print(delta0fit+i*.001)

    Diff_intensity = UnDi.diffractionintensity(yy,Ytran,i,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+i*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit)


    return Diff_intensity

def objective_diffraction(params, yy, data,weights):
    ndata,nx = data.shape
    
    resid = 0.0*data[:]
    
    for Strpidx in range(ndata):
        if Strpidx > 0:
            UnDi.SingleAmplitudeDone = 1
        resid[Strpidx] = (diffraction_dataset(params,listOfboots[Strpidx],yy) - data[Strpidx])*weights[Strpidx]#*(1.+50000.*yy**2.)#**2.
        
    return np.array(resid.flatten())

maybeBadGamma = 0

badk = []

DoNotSave = 0

Amplfit = 0.95769674
Ytran = 1.5618e-04
gammafit = 352.408492
delta0fit = 0.20647915
bfit = 0.58404239

#Amplfit_709:  0.95769674 +/- 0.00151847 (0.16%) == 'Amplfit_596'
#    Ytran_709:    1.5618e-04 +/- 7.2134e-07 (0.46%) == 'Ytran_596'
#    gamma_709:    352.408492 +/- 0.04331662 (0.01%) == 'gamma_596'
#    delta0_709:   0.20647915 +/- 2.0594e-04 (0.10%) == 'delta0_596'
#    bfit_709:     0.58404239 +/- 0.01084113 (1.86%) == 'bfit_596'

kstep = 1

WLofPixel = (chosenpx - pixelof404nm)*dispersion + CalibWl

k0 = 000
k = k0

checkstart = 1

i = 100

Posi = 16#156
allNpics = 5
mm=.001

def dgovergnolambda(gamma,dTheta):
    return gamma**2.*dTheta**2./(1+gamma**2.*dTheta**2.)

print(dgovergnolambda(352.4,7.3*10**-7))

dThetaLst=np.linspace(0,.16,100)

dgammaLst = dgovergnolambda(352.4,dThetaLst/1000.)

fntsze = 17

fig, ax = plt.subplots(figsize=(8,6))
ax.set(title=r'Relative error of $\gamma$ vs. $\delta\Theta$')
ax.plot(dThetaLst, 1000.*dgammaLst)
ax.grid()
ax.set_xlabel(r'$\delta\Theta$ [mrad]',fontsize=fntsze)
ax.set_ylabel(r'$\delta\gamma/\gamma $ e-3',fontsize=fntsze)
ax.tick_params(axis='x', labelsize=fntsze)
ax.tick_params(axis='y', labelsize=fntsze)
plt.show()

#"""
def plot_for_offset(power):
    # Data for plotting
    t = np.arange(0.0, 100, 1)

    s = 1.*UnDi.diffractionintensitymodY(zetaprime2,Ytran,1.*(-.1*mm+power*mm/100.),allNpics*Posi,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+allNpics*Posi*.001,gammafit,Amplfit,bfit,5.3,Yscreenfit)

    fig, ax = plt.subplots(figsize=(10,5))
    ax.set(title='Simulated Diffraction')
    ax.plot(zetaprime2*1000., s)
    ax.grid()
    ax.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
    ax.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
    #axA.tick_params(axis='both', which='major', labelsize=fntsze)
    #ax.set(xlabel='X', ylabel='x^{}'.format(power),title='Powers of x')

    # IMPORTANT ANIMATION CODE HERE
    # Used to keep the limits constant
    #ax.set_ylim(0, 80)

    # Used to return the plot as an image rray
    fig.canvas.draw()       # draw the canvas, cache the renderer
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return image

kwargs_write = {'fps':1.0, 'quantizer':'nq'}
imageio.mimsave('./powers.gif', [plot_for_offset(i) for i in range(40)], fps=8)
#"""
"""
power = 42
s=UnDi.diffractionintensitymodY(zetaprime2,Ytran,1.*(-4*mm+power*mm/10.),allNpics*Posi,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+allNpics*Posi*.001,gammafit,Amplfit,bfit,5.3,Yscreenfit)

fig, ax = plt.subplots(figsize=(10,5))
ax.plot(zetaprime2, s)
ax.grid()
plt.show()


#axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
#axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
#axA.tick_params(axis='both', which='major', labelsize=fntsze)
#ax.set(xlabel='X', ylabel='x^{}'.format(power),title='Powers of x')


# IMPORTANT ANIMATION CODE HERE
# Used to keep the limits constant
ax.set_ylim(0, 80)
#"""