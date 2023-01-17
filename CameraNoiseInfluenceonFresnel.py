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

from scipy.optimize import curve_fit

from matplotlib.animation import FuncAnimation, PillowWriter 
import random

ratio = 1000. / 606#606. #how many pixel are aperture
UnDi.apy = .004

zetaprime2 = np.linspace(-ratio*UnDi.apy,ratio*UnDi.apy,1000)

dispersion = 0.004808759252942803e-9 #Former 0.0048085643228767475e-9
pixelof404nm = 1315.3024745554976 #Former: 1243.8551272671261
CalibWl = 404.6565e-9

#(407.7837-404.6565)

def IF(r,d):
    return (1/r+1/(r-d))**2.

dLst = np.linspace(.52,.52+1.,1000)

fntsze = 17

"""
fig, ax = plt.subplots(figsize=(8,6))
ax.set(title=r'Relative Intensity vs. undulator position')
ax.plot(dLst, IF(AUndu.z0,dLst)/IF(AUndu.z0,dLst[0]),"g")
ax.grid()
ax.set_xlabel(r'Undulator position d [m]',fontsize=fntsze)
ax.set_ylabel(r'Relative Intensity',fontsize=fntsze)
ax.tick_params(axis='x', labelsize=fntsze)
ax.tick_params(axis='y', labelsize=fntsze)
plt.show()
"""

Blst = [
129.859374,
129.932969,
129.734838,
129.942522,
129.764601,
129.810968,
129.999594,
129.884920,
129.706112,
129.781073,
129.805769
]



stdB = np.std(np.array(Blst))
print(stdB)


#lRlist = np.array([400.,404.,408.,412.,416.])*10**-9
#gammalist = np.array([352.404,352.397,352.392,352.387,352.383])

lRlist = np.array([404.,408.,412.,416.])#*10**-9
gammalist = 1000.*(np.array([352.397,352.392,352.387,352.383])*.511-180.)
chisqr = np.array([296546.,294782.,294511.,295534.])/294500

fig, ax1 = plt.subplots(figsize=(8,6))
ax1.set(title=r'Relative chisqr and energy vs. $\lambda_R$')
ax1.plot(lRlist, gammalist,"g.")
ax1.grid()

ax1.set_xlabel(r'Resonance wavelength $\lambda_R$ [nm]',fontsize=fntsze)
ax1.set_ylabel(r'energy [keV] + 180 MeV', color = 'green',fontsize=fntsze)
ax1.tick_params(axis='x', labelsize=fntsze)
ax1.tick_params(axis='y', labelcolor = 'green', labelsize=fntsze)
#ax.set_ylim(352.38,352.41)
ax2 = ax1.twinx() 
ax2.set_ylabel('Relative chisqr', color = 'blue', fontsize=fntsze)
ax2.plot(lRlist, chisqr, color = 'blue')
#ax2.tick_params(axis='y', color = 'blue', labelsize=fntsze)
ax2.tick_params(axis ='y', labelcolor = 'blue', labelsize=fntsze)
ax2.set_ylim(.999,1.01)

plt.show()

chosenpx = 1300

Yscreenfit = .3e-4#-.5e-4
LambdaRfit = 400e-9
Amplfit = 0.95769674
Ytran = 1.5618e-04
gammafit = 352.408492
delta0fit = 0.20647915
bfit = 1.5#0.58404239

WLofPixel = (chosenpx - pixelof404nm)*dispersion + CalibWl

Posi = 14.5#156
allNpics = 5
mm=.001

def dgovergnolambda(gamma,dTheta):
    return gamma**2.*dTheta**2./(1+gamma**2.*dTheta**2.)

print(dgovergnolambda(352.4,7.3*10**-7))

dThetaLst=np.linspace(0,.16,100)

dgammaLst = dgovergnolambda(352.4,dThetaLst/1000.)

Powaa = -0
"""
bfit = 0
sb0 = UnDi.diffractionintensity(zetaprime2,Ytran+1.*(-.1*mm+Powaa*mm/100.),allNpics*Posi,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+allNpics*Posi*.001,gammafit,Amplfit,bfit,5.2,Yscreenfit)
bfit = .8
sb2 = UnDi.diffractionintensity(zetaprime2,Ytran+1.*(-.1*mm+Powaa*mm/100.),allNpics*Posi,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+allNpics*Posi*.001,gammafit,Amplfit,bfit,5.2,Yscreenfit)

fig, ax = plt.subplots(figsize=(8,6))
ax.set(title=r'Intensity vs. angle')
ax.plot(range(len(zetaprime2)), sb0,"r")
#ax.plot(zetaprime2, sb1,"g")
ax.plot(range(len(zetaprime2)), sb2,"b")
ax.grid()
ax.set_xlabel(r'Position on screen [pixel]',fontsize=fntsze)
ax.set_ylabel(r'Intensity [a.u.]',fontsize=fntsze)
ax.tick_params(axis='x', labelsize=fntsze)
ax.tick_params(axis='y', labelsize=fntsze)
ax.legend([r'P$_b=0$',r'P$_b=0.8$'])
plt.show()
#"""
"""
#AUndu.z0 =100000.
print(AUndu.z0)
Posi = 22
zpfit = 5.23
sb1 = UnDi.diffractionintensity(zetaprime2,Ytran+1.*(-.1*mm+Powaa*mm/100.),allNpics*Posi,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+allNpics*Posi*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit)
zpfit = 5.45
AUndu.z0 =10.8
sb2 = UnDi.diffractionintensity(zetaprime2,Ytran+1.*(-.1*mm+Powaa*mm/100.),allNpics*Posi,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+allNpics*Posi*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit)

fig, ax = plt.subplots(figsize=(8,6))
ax.set(title=r'Intensity vs. angle')
#ax.plot(zetaprime2, sb0,"r")
ax.plot(range(len(zetaprime2)), sb1,"g")
ax.plot(range(len(zetaprime2)), sb2,"b")
ax.grid()
ax.set_xlabel(r'Position on screen [pixel]',fontsize=fntsze)
ax.set_ylabel(r'Intensity [a.u.]',fontsize=fntsze)
ax.tick_params(axis='x', labelsize=fntsze)
ax.tick_params(axis='y', labelsize=fntsze)
ax.legend([r'P$_z=5.23$m, with $z_0=11.45$m',r'P$_z=5.45$m, with $z_0=10.8$m'])
plt.show()
#"""

Posi = 42
zpfit = 5.23
LambdaRfit = 400*10**-9
sb1 = UnDi.diffractionintensity(zetaprime2,Ytran+1.*(-.1*mm+Powaa*mm/100.),allNpics*Posi,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+allNpics*Posi*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit)
LambdaRfit = 404*10**-9

sb2 = .7*UnDi.diffractionintensity(zetaprime2,Ytran+1.*(-.1*mm+Powaa*mm/100.),allNpics*Posi,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+allNpics*Posi*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit)

fig, ax = plt.subplots(figsize=(8,6))
ax.set(title=r'Intensity vs. angle')
#ax.plot(zetaprime2, sb0,"r")
ax.plot(range(len(zetaprime2)), sb1,"g")
ax.plot(range(len(zetaprime2)), sb2,"b")
ax.grid()
ax.set_xlabel(r'Position on screen [pixel]',fontsize=fntsze)
ax.set_ylabel(r'Intensity [a.u.]',fontsize=fntsze)
ax.tick_params(axis='x', labelsize=fntsze)
ax.tick_params(axis='y', labelsize=fntsze)
ax.legend([r'P$_z=5.23$m, with $z_0=11.45$m',r'P$_z=5.45$m, with $z_0=10.8$m'])
plt.show()

exit()

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