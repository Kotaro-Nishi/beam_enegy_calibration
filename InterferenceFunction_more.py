import matplotlib.pyplot as plt
import AmplitudeFormula as AUndu
from scipy import ndimage

from numpy import transpose as nptp
import numpy as np

def convolveplain(Sourcelst1D1000,ElectronSource):

    Sourcelst1D1000real = np.array(Sourcelst1D1000).real
    Sourcelst1D1000imag = np.array(Sourcelst1D1000).imag

    Sourcelst1D1000real = ndimage.convolve(Sourcelst1D1000real,ElectronSource)
    Sourcelst1D1000imag = ndimage.convolve(Sourcelst1D1000imag,ElectronSource)

    Sourcelst1D1000 = Sourcelst1D1000real + 1.j*Sourcelst1D1000imag

    return Sourcelst1D1000

def convolvefresnel(Sourcelst1D1000,fresnelforward,ElectronSource):

    Sourcelst1D1000real = np.array(Sourcelst1D1000).real
    Sourcelst1D1000imag = np.array(Sourcelst1D1000).imag

    Sourcelst1D1000real = ndimage.convolve(Sourcelst1D1000real,ElectronSource)
    Sourcelst1D1000imag = ndimage.convolve(Sourcelst1D1000imag,ElectronSource)

    Sourcelst1D1000 = Sourcelst1D1000real + 1.j*Sourcelst1D1000imag

    fresnelforward_f = np.fft.fft(fresnelforward)

    orig_f = np.fft.fft(np.array(Sourcelst1D1000))
    origprime_f = nptp(orig_f * nptp(fresnelforward_f))
    origprime = np.fft.ifft(origprime_f)

    return np.fft.fftshift(origprime)
    
def ElectronSourcefunction(yy):
    Esize = .00009#.00009
    sourcesamples = 51
    sourcesampleshalf = int((sourcesamples-1)/2)

    ElectronSource = []
    for i in range(sourcesamples):
        ElectronSource.append(np.exp(-(((yy[i-int((len(yy))/2)-sourcesampleshalf])**2)/Esize**2.)/2.)/2000.*Esize) #denominator is not necessary for making unit
    return ElectronSource

def ElectronSourcefunctionSet(yy,Esize):
    sourcesamples = 51
    sourcesampleshalf = int((sourcesamples-1)/2)

    ElectronSource = []
    for i in range(sourcesamples):
        ElectronSource.append(np.exp(-(((yy[i-int((len(yy))/2)-sourcesampleshalf])**2)/Esize**2.)/2.)/2000.*Esize) #denominator is not necessary for making unit
    return ElectronSource

def squaredphase(yy,Ytran,deltad,zp,WLofPixel):
    fresnelforward = []
    for i in range(len(yy)):
        thetaradius = abs(yy[i]-Ytran) #add pythagoras if 2D
        fresnelforward.append(np.exp(1.j * 2.*np.pi/WLofPixel * (((yy[i]-Ytran)**2)/(2.*zp)))*np.exp(1.j * 2.*np.pi/WLofPixel * (((yy[i]-Ytran)**2)/(2.*(AUndu.z0-deltad/1.)))))#change 20210512 changed deltad/2. to deltad/1.
        #fresnelforward.append(np.exp(1.j * 2.*np.pi/WLofPixel * (((yy[i]-Ytran)**2)/(2.*zp))))#change 20210512 changed deltad/2. to deltad/1.
        
    return fresnelforward