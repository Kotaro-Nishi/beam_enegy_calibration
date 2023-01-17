################################################################################################################
# Name: AmplitudeFormula.py
# 
# Purpose: return x and y Amplitude of one undulator
#
# Things to fix:...
#
################################################################################################################

import scipy.special as special
import numpy as np

global LambU
global n_U
global Lu

LambU,n_U,Lu = .08,6,.52

global z0
global z1
global z2
global zc
global f
global zp

z0,z1,z2,zc,f,zp = 11.45,3.03,.32,.68,1.,3.#11.45,3.03,.32,.68,1.,3.


def LambdaR(LambU,Ku,gamma,Theta):
    return LambU/(2.*gamma**2)*(1+Ku**2./2.+Theta**2*gamma**2/2.)

def AA(Ku,gamma,Theta):
    return 1.+Ku**2/2.+gamma**2*Theta**2

def XX(Ku,gamma,AA,Theta,phi):
    return (2*gamma*Theta*Ku*np.cos(phi))/AA

def YY(Ku,gamma,AA,Theta):
    return Ku**2/(4*AA)

def Sq(q,xx,yy):
    SqList = []
    for pp in range(-5,5):
        SqList.append(special.jv(pp,yy)*special.jv(1+2*pp+q,xx))
    return sum(SqList)

def Ax(xx,yy,Ku,gamma,Theta,phi):
    return 2.*gamma*Theta*np.cos(phi)*Sq(0,xx,yy)-Ku*(Sq(1,xx,yy)+Sq(-1,xx,yy))

def Ay(xx,yy,Ku,gamma,Theta,phi):
    return 2.*gamma*Theta*np.sin(phi)*Sq(0,xx,yy)

def Amplitude(gamma,WLofPixel,LambdaRfit,Ku,Theta,phi):

                                                         # LambdaRfit is guessed Resonant Wavelength

    #if LambdaRfit == WLofPixel:
    #    LambdaRfit = LambdaRfit + 1.e-11

    #LR = LambdaRfit#LambdaR(LambU,Ku,gamma,Theta)                   # Resonant wavelength at angle Theta
    LR = LambdaR(LambU,Ku,gamma,Theta)                   # Resonant wavelength at angle Theta
    LP = WLofPixel                                       # Wavelength of a pixel

    aa = AA(Ku,gamma,Theta)                              # Walker called this A

    xx = XX(Ku,gamma,aa,Theta,phi)                       # Walker called this X
    yy = YY(Ku,gamma,aa,Theta)                           # Walker called this Y

    ax = Ax(xx,yy,Ku,gamma,Theta,phi)                    # x-component of Vector potential
    ay = Ay(xx,yy,Ku,gamma,Theta,phi)                    # y-component of Vector potential
    
    return [ax/aa * np.sin(n_U*np.pi*(LR-LP)/LR)/np.sin(np.pi*(LR-LP)/LR), ay/aa * np.sin(n_U*np.pi*(LR-LP)/LR)/np.sin(np.pi*(LR-LP)/LR)]