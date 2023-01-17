################################################################################################################
# Name: AvarageWavelengthsToMakeFresnelVis.py
#
# Purpose: Avaragingscript which makes the fresnel diffraction more visible
#			
# Things to fix: *.*
#
#################################################################################################################

import matplotlib.pyplot as plt

import TransverseUndulator_function_Walker as Undu

from numpy import transpose as nptp
import cv2
import os, sys
import csv
import re
import numpy as np
from natsort import natsorted
from scipy.optimize import curve_fit
#from lmfit import minimize, Parameters, report_fit

#from lmfit import Minimizer

lambdaU=.08
rat4 = 1.65
a2 = .004#.0040
failed = False

import random

def secmax(lst):
    maxlst = max(lst)
    return [c for c in lst if c != maxlst]

def secmin(lst):
    minlst = min(lst)
    return [c for c in lst if c != minlst]

def linefunc(xxx,a1,b1):
    return b1*xxx+a1

def per_iteration(fit_params, iter, resid, *args, **kws):

    global deter
    print("ampl:      " + str(fit_params['ampl_1'].value))   
    if 'asymmetry_1' not in deter:
        print("asymmetry: " + str(fit_params['asymmetry_1'].value))
    print('offset: ' + str(fit_params['offset_1'].value))
    print("sFactor: " + str(fit_params['sFactor_1'].value))
    print("deltad: " + str(fit_params['deltad_1'].value))
    #print("dispersion: " + str(fit_params['dispersion_1'].value))
    print("center0: " + str(fit_params['center0_1'].value))

    global asymmetryfirst

    lineresiduals = 1.

    NumOfstripe = int(len(resid)/1000)

    residN = resid.reshape(NumOfstripe, 1000)
    residNLst = []
    for res in range(len(residN)):
        print(sum(residN[res]**2.))
        residNLst.append(sum(residN[res]**2.))

    print(iter)
    print(sum(resid**2.))

    if iter >90:
        print("Please quit")
        return True
    #"""
    if (max(residNLst)<.4): #160
        return True
        
    if iter > 150 and max(secmax(residNLst))<40.:
        print("max after fit: " + str(max(residNLst)))
        print("all resids: " + str(residNLst))
        return True
    
    if iter > 190 and min(secmin(secmin(residNLst)))<6.:
        print("min after fit")
        print("all resids: " + str(residNLst))
        return True

    if iter > 210 and max(residNLst)>1000.:
        return True

    if iter > 230:
        return True
    #"""
    else:
        return False

def per_iterationfinal(fit_params, iter, resid, *args, **kws):
    #print("asymmetry: " + str(fit_params['asymmetry_1'].value))
    #print('offset: ' + str(fit_params['offset_1'].value))
    #print("sFactor: " + str(fit_params['sFactor_1'].value))
    #print("deltad: " + str(fit_params['deltad_1'].value))
    #print("dispersion: " + str(fit_params['dispersion_1'].value))    

    print("per_iterationfinal")

    lineresiduals = 1.

    NumOfstripe = int(len(resid)/1000)

    residN = resid.reshape(NumOfstripe, 1000)
    residNLst = []
    for res in range(len(residN)):
        print(sum(residN[res]**2.))
        residNLst.append(sum(residN[res]**2.))
        
    print(iter)
    print(sum(resid**2.))

    if iter >444:
        return True

    if iter > 100:
        if (min(residNLst)<.9): #160
            return True
        
    if iter > 150 and max(secmax(residNLst))<40.:
        print("max after fit: " + str(max(residNLst)))
        print("all resids: " + str(residNLst))
        return True
    
    if iter > 190 and min(secmin(secmin(residNLst)))<6.:
        print("min after fit")
        print("all resids: " + str(residNLst))
        return True

    if iter > 210 and max(residNLst)>1000.:
        return True

    if iter > 230:
        return True
        
    else:
        return False

deter = 0

asymmetrydeter = 0

def diffraction_dataset(params, i, x,y):
    global deter
    global asymmetrydeter
    ampl      = params['ampl_%i' % (i+1)].value
    if 'asymmetry_1' not in deter:
        asymmetry = params['asymmetry_%i' % (i+1)].value
    else:
        asymmetry = asymmetrydeter
    sFactord   = params['sFactor_%i' % (i+1)].value
    center0     = params['center0_%i' % (i+1)].value
    offset    = params['offset_%i' % (i+1)].value
    #dispersion    = params['dispersion_%i' % (i+1)].value
    deltad    = params['deltad_%i' % (i+1)].value
    
    global amplfixed
    global asymmetryfixed
    global sFactordfixed
    global center0fixed
    global offsetfixed
    #global dispersionfixed
    global delta0fixed

    global corr

    amplfixed = ampl
    asymmetryfixed = asymmetry
    sFactordfixed = sFactord
    center0fixed = center0
    offsetfixed = offset
    #dispersionfixed = dispersion
    delta0fixed = deltad
    
    Strpidx = i
    
    return diffraction(x,y,Strpidx,ampl,asymmetry,sFactord,center0,offset,deltad)

def diffraction_datasetfinal(params, i, x,y):
    if 'ampl_1' not in corr:
        ampl = params['ampl_%i' % (i+1)].value
    else:
        ampl = amplfixed

    if 'asymmetry_1' not in corr:
        asymmetry = params['asymmetry_%i' % (i+1)].value
    else:
        asymmetry = asymmetryfixed

    if 'sFactor_1' not in corr:
        sFactord = params['sFactor_%i' % (i+1)].value
    else:
        sFactord = sFactordfixed

    #if 'dispersion_1' not in corr:
    #    dispersion = params['dispersion_%i' % (i+1)].value
    #else:
    #    dispersion = dispersionfixed

    if 'center0_1' not in corr:
        center0 = params['center0_%i' % (i+1)].value
    else:
        center0 = center0fixed
        
    if 'offset_1' not in corr:
        offset = params['offset_%i' % (i+1)].value
    else:
        offset = offsetfixed

    if 'deltad_1' not in corr:
        deltad = params['deltad_%i' % (i+1)].value
    else:
        deltad = delta0fixed
    
    Strpidx = i
    return diffraction(x,y,Strpidx,ampl,asymmetry,sFactord,center0,offset,deltad)

def objective_diffraction(params, x, data,weights):
    """ calculate total residual for fits to several data sets held
    in a 2-D array, and modeled by Gaussian functions"""
    ndata, nx = data.shape #spli tupel

    resid = 0.0*data[:]

    # make residual per data set
    for i in range(ndata):
        resid[i] = (data[i] - diffraction_dataset(params, i, x))*weights[i]**2.
    # now flatten this to a 1D array, as minimize() needs
    return np.array(resid.flatten())
    
def objective_diffractionfinal(params, x,y, data,weights):
    """ calculate total residual for fits to several data sets held
    in a 2-D array, and modeled by Gaussian functions"""
    ndata, nx = data.shape #spli tupel

    resid = 0.0*data[:]

    # make residual per data set
    for i in range(ndata):
        resid[i] = (data[i] - diffraction_datasetfinal(params, i, x))*weights[i]**2.
    # now flatten this to a 1D array, as minimize() needs
    return np.array(resid.flatten())

zetaprime1=np.linspace( -(rat4)*a2,(rat4)*a2,100)
zetaprime2=np.linspace( -(rat4)*a2,(rat4)*a2,100)

####################

def diffraction(xx,yy,Strpidx,ampl,asymmetry,sFactord,center0,offset,deltad):
    distance=3.444

    WLofPixel = 409.140e-9 +Strpidx*-.478e-9*1.

    k = 2.*np.pi/WLofPixel 

    ElectronSource = []
    Esize = .00007###
    sourcesamples = 51
    sourcesampleshalf = int((sourcesamples-1)/2)

    for i in range(sourcesamples):
        ElectronSource.append(np.exp(-(((xx[i-int((len(xx))/2)-sourcesampleshalf])**2)/Esize**2.)/2.)/(2000.*Esize))

    #fig1 = plt.figure()
    #ax2 = fig1.add_subplot(111)
    #ax2.plot(range(len(ElectronSource)),ElectronSource,'r')
    #plt.show()

    ElectronSource = np.array(ElectronSource)

    fresnelforward2D = []
    for i in range(len(xx)):
        fresnelforward = []
        for ii in range(len(yy)):
            fresnelforward.append(np.exp(1.j * k*(xx[i]**2 + yy[ii]**2)/(2.*distance)))
        fresnelforward2D.append(fresnelforward)

    #fresnelforward = np.convolve(fresnelforward, ElectronSource)[sourcesampleshalf:1000+sourcesampleshalf]/36.

    fresnelforward = np.transpose(fresnelforward2D)

    fresnelforward_f = np.fft.fft2(np.array(fresnelforward))#[0:100]

    #fresnelforward_f = np.transpose(fresnelforward_f[0:100])

    #fresnelforward_f = np.transpose(fresnelforward_f)

    """
    fig2 = plt.figure()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Angle [mm]')
    #plt.plot(1000.*zetaprime2[0:100], fresnelforward2D[47], '-')
    #plt.imshow((np.array(fresnelforward).real))
    plt.imshow((np.array(fresnelforward_f).imag))
    plt.show()
    #"""


    Sourcelst2D=[]

    for ii in range(len(yy)):
        Sourcelst = []
        for i in range(len(xx)):
            if (xx[i]-asymmetry > -a2*1. and xx[i]-asymmetry < a2*1.) and (yy[ii]-asymmetry > -a2*.5 and yy[ii]-asymmetry < a2*.5):
                lambdaNull = Undu.lambda0fit(center0,(((xx[i])**2 + yy[ii]**2)**.5)/sFactord,lambdaU)
                deltaL = WLofPixel-lambdaNull
                #Sourcelst.append(Undu.analyticalspecPhase((xx[i]-asymmetry*10**-5-offset)/sFactord,deltaL,center0,WLofPixel,deltad))
                Sourcelst.append(Undu.analyticalspecPhase((((xx[i])**2 + yy[ii]**2)**.5)/sFactord,deltaL,center0,WLofPixel,deltad))
            else:
                Sourcelst.append(0.)
        Sourcelst2D.append(Sourcelst)

    Sourcelst = np.transpose(Sourcelst2D)

    orig_f = np.fft.fft2(np.array(Sourcelst))
    #print("after Source FFT")

    origprime_f = orig_f * fresnelforward_f

    origprime = ampl*np.fft.ifft2(origprime_f)

    result = np.conjugate(origprime)*origprime
    #print(result[51].real)
    #result = np.conjugate(Sourcelst)*Sourcelst

    result = np.fft.fftshift(result)
    gfdsgfds = sum(np.transpose(result)).real
    gfdsgfds2 = sum(np.transpose(np.transpose(np.conjugate(Sourcelst2D)*Sourcelst2D))).real

    fig2 = plt.figure()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Angle [mm]')
    plt.plot(range(len(gfdsgfds2)),gfdsgfds2)

    fig2 = plt.figure()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Angle [mm]')
    plt.plot(range(len(gfdsgfds)),gfdsgfds)#,aspect=20,extent=(wlLstplot[0],wlLstplot[-1]+1,Thetaplot[0],Thetaplot[-1]))
    #plt.show()

    fig2 = plt.figure()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Angle [mm]')
    plt.imshow(result.real)#,aspect=20,extent=(wlLstplot[0],wlLstplot[-1]+1,Thetaplot[0],Thetaplot[-1]))
    #plt.show()

    fig2 = plt.figure()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Angle [mm]')
    plt.imshow(np.transpose(np.conjugate(Sourcelst2D)*Sourcelst2D).real)#,aspect=20,extent=(wlLstplot[0],wlLstplot[-1]+1,Thetaplot[0],Thetaplot[-1]))
    plt.show()



    return(abs(np.array(result)))

####################

"""
arg = sys.argv

tiflist = []

filelist = os.listdir(arg[1])
for i in range(len(filelist)):
	if filelist[i][-4:] == "tiff" or filelist[i][-3:] == "tif":
		tiflist.append(filelist[i])
filelist=natsorted(tiflist)

Commonpath = str(os.path.dirname(os.path.abspath(__file__)))

print(Commonpath)

DataPath = arg[1]

print(DataPath)

slashes = [i for i, ltr in enumerate(DataPath) if ltr == "\\"]
DataPathBase = DataPath[0:slashes[-2]]

print(DataPathBase)
"""
datarangelow=1804

#Guessvalues=[.122,.0004,9.5,408.3e-9,0.0002,delta0+k*.001]

#print diffraction(zetaprime2,zetaprime2,2,ampl,asymmetry,sFactord,center0,offset,deltad)

#print(diffraction(zetaprime1,zetaprime2,2,.122,.0004,9.5,408.3e-9,.0002,.2))
for i in range(200):
    print(.09+.1*i/100.)
    UnduDiffractionImage = diffraction(zetaprime1,zetaprime2,2,.122,.000,9.5,408.3e-9,.000,.09+.1*i/100.)

fig2 = plt.figure()
plt.xlabel('Wavelength [nm]')
plt.ylabel('Angle [mm]')
plt.imshow(UnduDiffractionImage,aspect=1)#,aspect=20,extent=(wlLstplot[0],wlLstplot[-1]+1,Thetaplot[0],Thetaplot[-1]))
plt.show()
exit()

Backgroundimage = cv2.imread(DataPathBase + "/BG/raw_shortExp/File1.tif",cv2.IMREAD_UNCHANGED)
ListOfAllMatrices=[]
ListOfAllSigmaMatrices=[]
NumOfstripes=5 #because data looks flat

WidthOfStripe = int((2304-datarangelow)/NumOfstripes)
print(WidthOfStripe)

MatrixOfSlit = []

# create 5 sets of parameters, one per data set

deltaON=True

delta0 = .207
delta0min = delta0-.008
delta0max = delta0+.008

fixedsize = .245#.075
k0=0#28#restart from k0=30 to find optimum
k=k0#+40#+120

maybeBadGamma = 0

targetfiles = ["Undupos.csv",
               "Amplitudes.csv",
               "center0s.csv",
               "sFactors.csv",
               "dispersions.csv",
               "Bsymmetries.csv",
               "Offsets.csv",
               "deltads.csv",
               "AmploffsetCorrel.csv",
               "AmpldeltaCorrel.csv",
               "sFactordeltaCorrel.csv"]

badk = []

determined = ["asymmetry_1"]



while k < int(len(filelist)/4):
    Guessvalues=[.122,.0004,9.5,408.3e-9,0.0002,delta0+k*.001]
    asymmetrydeter = .00036
    correlpair = []
    NoneCases = 0
    while True:
        print(badk)
        if len(badk) > 0:
            if badk[0] < k - 100:
                del badk[0]
                maybeBadGamma -= 1

        if k==k0:
            Undu.gamma = Undu.gamma-.05
            checklist = os.listdir(DataPath + "\SharedP")
            maybeBadGamma = 0
            for ckid in range(len(targetfiles)):
                if targetfiles[ckid] in checklist:
                    os.remove(DataPath + "\SharedP" + "//" + targetfiles[ckid])

        GuessWLth = Guessvalues[3]

        OneLineEachFromTenPicturesAtOnePosition=[]
        OneLineEachFromTenPicturesAtOnePositionForSigma=[]
        
        SigmaOneLineEachFromTenPicturesAtOnePosition=[]
        for i in range(4):
            noti = i
            PictureMinusBG = cv2.imread(DataPath + "/" + filelist[noti+k*4],cv2.IMREAD_UNCHANGED).astype(int)-Backgroundimage*.99
            PictureMinusBGForSigma = cv2.imread(DataPath + "/" + filelist[i+k*4],cv2.IMREAD_UNCHANGED).astype(int)-Backgroundimage*.99
            
            PictureMinusBG = np.float32(PictureMinusBG)
            PictureMinusBGForSigma = np.float32(PictureMinusBGForSigma)
            
            PictureMinusBG=np.transpose(np.transpose(PictureMinusBG)[(datarangelow-100):])
            PictureMinusBGForSigma=np.transpose(np.transpose(PictureMinusBGForSigma)[(datarangelow-100):])
            
            OneLineEachFromTenPicturesAtOnePosition.append(PictureMinusBG)#.tolist())
            OneLineEachFromTenPicturesAtOnePositionForSigma.append(PictureMinusBGForSigma)#.tolist())

        TenPixelFromTenPicturesTimesLineAtOnePosition = sum(OneLineEachFromTenPicturesAtOnePosition)/4.
        TenPixelFromTenPicturesTimesLineAtOnePositionForSigma = sum(OneLineEachFromTenPicturesAtOnePositionForSigma)/4.
        
        print(filelist[i+k*4])

        #SigmaOneLineEachFromTenPicturesAtOnePosition = np.std(OneLineEachFromTenPicturesAtOnePosition,0,ddof=1).tolist()#ddof=1 because we divide by (N-1)

        #cv2.imshow("imahe",TenPixelFromTenPicturesTimesLineAtOnePosition)
        #cv2.waitKey(0)
        interpolatedpicture=[]
        interpolatedsigmapicture=[]
        for ii in range(NumOfstripes):
            stripe = np.transpose(np.transpose(TenPixelFromTenPicturesTimesLineAtOnePosition)[WidthOfStripe*ii:WidthOfStripe+WidthOfStripe*ii])
            stripesigma = np.transpose(np.transpose(TenPixelFromTenPicturesTimesLineAtOnePositionForSigma)[WidthOfStripe*ii:WidthOfStripe+WidthOfStripe*ii])#/.5

            stripe = cv2.resize(stripe,(1,len(TenPixelFromTenPicturesTimesLineAtOnePosition)),interpolation=cv2.INTER_AREA) #1 zu 1000
            #stripesigma = cv2.resize(stripesigma,(1,len(TenPixelFromTenPicturesTimesLineAtOnePosition)),interpolation=cv2.INTER_AREA) #1 zu 1000

            print(len(stripe))
            print(len(stripe[0]))

            stripesigma = np.std(np.transpose(stripesigma),0,ddof=1).tolist()

            interpolatedpicture.append(np.transpose(stripe)[0])
            interpolatedsigmapicture.append(np.transpose(stripesigma))
            
            data2D = np.transpose(stripe)[0]/200.
            datasigma2D = np.transpose(stripesigma)[0]/200.
            
        data = np.array(interpolatedpicture)
        datasigma = np.array(interpolatedsigmapicture)


        #assert(data.shape) == (5, 151)
        print(k)
        print((NumOfstripes/5.))

        deter = determined

        fit_params = Parameters()####[1.,.00001,.13,-4.5*10**-9.,0.0001,.00009]
        for iy, y in enumerate(data):#,7.14716241e-03, 5.55463305e-04, 9.97300981e-02, 7.52119302e-09, 1.18235062e-03)
            fit_params.add( 'ampl_%i' % (iy+1), value=Guessvalues[0], min=.02, max=.9)
            if 'asymmetry_1' not in determined:
                fit_params.add( 'asymmetry_%i' % (iy+1), value=Guessvalues[1], min=0.0003, max=0.0007)#1.55463305e-05)
                
            fit_params.add( 'sFactor_%i' % (iy+1), value=Guessvalues[2], min=8.5,max=113.)
            #fit_params.add( 'dispersion_%i' % (iy+1), value=.478e-9/(NumOfstripes/5.), min=.477e-9/(NumOfstripes/5.),max=.479e-9/(NumOfstripes/5.))
            fit_params.add( 'center0_%i' % (iy+1), value=GuessWLth, min=GuessWLth-4.5*10**-9, max=GuessWLth+3.5*10**-9)
            fit_params.add( 'offset_%i' % (iy+1), value=Guessvalues[4], min=.0000, max=.0008)#1.18235062e-03)
            #fit_params.add( 'offset_%i' % (iy+1), value=0., min=-.00001, max=.00001)#1.18235062e-03)
            fit_params.add( 'deltad_%i' % (iy+1), value=Guessvalues[5], min=delta0min+k*.001, max=delta0max+k*.001)

        if 'asymmetry_1' not in determined:
            for iy in tuple(range(NumOfstripes+1)[2:]):
                fit_params['asymmetry_%i' % iy].expr='asymmetry_1'

        for iy in tuple(range(NumOfstripes+1)[2:]):
            fit_params['ampl_%i' % iy].expr='ampl_1'

        for iy in tuple(range(NumOfstripes+1)[2:]):
            fit_params['center0_%i' % iy].expr='center0_1'

        #for iy in tuple(range(NumOfstripes+1)[2:]):
        #    fit_params['dispersion_%i' % iy].expr='dispersion_1'

        for iy in tuple(range(NumOfstripes+1)[2:]):
            fit_params['sFactor_%i' % iy].expr='sFactor_1'

        for iy in tuple(range(NumOfstripes+1)[2:]):
            fit_params['offset_%i' % iy].expr='offset_1'

        for iy in tuple(range(NumOfstripes+1)[2:]):
            fit_params['deltad_%i' % iy].expr='deltad_1'

        #while True:

        #minimize(objective_diffraction, fit_params,iter_cb=per_iteration, args=(zetaprime2[0:1000], data, 1./datasigma),)
        #print(datasigma[0])
        for ghg in range(len(datasigma)):
            print(min(datasigma[ghg]))
            print(max(datasigma[ghg]))

        mini = Minimizer(objective_diffraction, fit_params,iter_cb=per_iteration, nan_policy = "raise",scale_covar=False , fcn_args=(zetaprime2[0:1000], data, 1./datasigma),)
        
        result = mini.leastsq(ftol=1.e-3)

        residNresult = (result.residual).reshape(NumOfstripes, 1000)
        
        residNLstresult = []
        for res in range(len(residNresult)):
            residNLstresult.append(sum(residNresult[res]**2.))

        report_fit(result.params)
        
        for key, par in result.params.items():
            print('{} {:.3f} +/- {}'.format(key, par.value, par.stderr))

        if result.params["ampl_1"].stderr!=None and min(residNLstresult)<1.5:
            DoNotSave = 0
            amplitudefit = result.params['ampl_1'].value
            
            if 'asymmetry_1' not in deter:
                asymmetryfit = result.params['asymmetry_1'].value
            else:
                asymmetryfit = asymmetrydeter
            
            sizefactor = result.params['sFactor_1'].value
            #dispersionfit = result.params['dispersion_1'].value
            center0fit = result.params['center0_1'].value
            offsetfit = result.params['offset_1'].value
            deltadfit = result.params['deltad_1'].value

            plainparameters = ["ampl_1","asymmetry_1","sFactor_1","center0_1","offset_1","deltad_1"]
            corrcoeff = 0.
            canceledlist = [1,1,1,1,1,1,1]

            """
            for cancelindex1 in range(6):
                for cancelindex2 in range(1+cancelindex1,6):
                    #print(plainparameters[cancelindex1],plainparameters[cancelindex2])
                    #print(result.params[plainparameters[cancelindex1]].correl[plainparameters[cancelindex2]])
                    if corrcoeff < abs(result.params[plainparameters[cancelindex1]].correl[plainparameters[cancelindex2]]):
                        corrcoeff = abs(result.params[plainparameters[cancelindex1]].correl[plainparameters[cancelindex2]])
                        #print(corrcoeff)
                        #print(plainparameters[cancelindex1],plainparameters[cancelindex2])
            if corrcoeff >.4:
            """

            correlpair.append(plainparameters[3])#,plainparameters[cancelindex2]]
            #correlpair.append(plainparameters[4])
            correlpair.append(plainparameters[1])
            break
            
        else:
            GuessvaluesAll=[[.095,.00035,9.0,408.3e-9,0.0003,delta0+k*.001-.0005],
             [.10,.00037,9.2,408.3e-9,0.00028,delta0+k*.001-.0003],
             [.11,.00039,9.4,408.3e-9,0.00026,delta0+k*.001-.0001],
             [.12,.00041,9.6,408.3e-9,0.00024,delta0+k*.001+.0001],
             [.13,.00043,9.8,408.3e-9,0.00022,delta0+k*.001+.0003],
             [.14,.00045,10.0,408.3e-9,0.0002,delta0+k*.001+.0005]]
            Guessvalues = GuessvaluesAll[NoneCases]
            NoneCases += 1
            if NoneCases > 5:
                DoNotSave = 1
                break
            else:
                DoNotSave = 0

        if DoNotSave > 0:
            break

    corr = correlpair
    #"""
    fit_params = Parameters()####[1.,.00001,.13,-4.5*10**-9.,0.0001,.00009]
    for iy, y in enumerate(data):#,7.14716241e-03, 5.55463305e-04, 9.97300981e-02, 7.52119302e-09, 1.18235062e-03)

        if 'ampl_1' not in correlpair:
            fit_params.add( 'ampl_%i' % (iy+1), value=.12, min=.02, max=.9)

        if 'asymmetry_1' not in correlpair:
            fit_params.add( 'asymmetry_%i' % (iy+1), value=.0004, min=0.0003, max=0.0007)#1.55463305e-05)

        if 'sFactor_1' not in correlpair:
            fit_params.add( 'sFactor_%i' % (iy+1), value=sizefactor, min=8.5,max=15.)

        #if 'dispersion_1' not in correlpair:
        #    fit_params.add( 'dispersion_%i' % (iy+1), value=dispersionfit, min=.477e-9/(NumOfstripes/5.),max=.479e-9/(NumOfstripes/5.))

        if 'center0_1' not in correlpair:
            fit_params.add( 'center0_%i' % (iy+1), value=center0fit, min=GuessWLth-4.5*10**-9, max=GuessWLth+3.5*10**-9)

        if 'offset_1' not in correlpair:
            fit_params.add( 'offset_%i' % (iy+1), value=offsetfit, min=.0000, max=.0006)#1.18235062e-03)

        if 'deltad_1' not in correlpair:
            fit_params.add( 'deltad_%i' % (iy+1), value=delta0+k*.001, min=delta0min+k*.001, max=delta0max+k*.001)

    #"""

    print(fit_params)

    if 'ampl_1' not in correlpair:
        for iy in tuple(range(NumOfstripes+1)[2:]):
            fit_params['ampl_%i' % iy].expr='ampl_1'

    if 'asymmetry_1' not in correlpair:
        for iy in tuple(range(NumOfstripes+1)[2:]):        
            fit_params['asymmetry_%i' % iy].expr='asymmetry_1'

    if 'sFactor_1' not in correlpair:
        for iy in tuple(range(NumOfstripes+1)[2:]):
            fit_params['sFactor_%i' % iy].expr='sFactor_1'

    #if 'dispersion_1' not in correlpair:
    #    for iy in tuple(range(NumOfstripes+1)[2:]):
    #        fit_params['dispersion_%i' % iy].expr='dispersion_1'

    if 'center0_1' not in correlpair:
        for iy in tuple(range(NumOfstripes+1)[2:]):
            fit_params['center0_%i' % iy].expr='center0_1'

    if 'offset_1' not in correlpair:
        for iy in tuple(range(NumOfstripes+1)[2:]):
            fit_params['offset_%i' % iy].expr='offset_1'

    if 'deltad_1' not in correlpair:
        for iy in tuple(range(NumOfstripes+1)[2:]):
            fit_params['deltad_%i' % iy].expr='deltad_1'



    print(fit_params)

    #for ghg in range(len(datasigma)):
    #    print(min(datasigma[ghg]))
    #    print(max(datasigma[ghg]))

    mini = Minimizer(objective_diffractionfinal, fit_params,iter_cb=per_iterationfinal, nan_policy = "raise",scale_covar=False , fcn_args=(zetaprime2[0:1000], data, 1./datasigma),)

    result = mini.leastsq(ftol=1.e-3)

    residNresult = (result.residual).reshape(NumOfstripes, 1000)
    
    if result.params["ampl_1"].stderr!=None:
        
        print(result.covar)
        
        #print((result.covar[0][0])**.5)
        #print((result.covar[1][1])**.5)
        #print((result.covar[2][2])**.5)
        #print((result.covar[3][3])**.5)
        
        ampldeltaCorrel = result.covar[0][3] / (result.params["ampl_1"].stderr*result.params["deltad_1"].stderr)
        sFactordeltaCorrel = result.covar[1][3] / (result.params["sFactor_1"].stderr*result.params["deltad_1"].stderr)
        amploffsetCorrel = result.covar[0][2] / (result.params["ampl_1"].stderr*result.params["offset_1"].stderr)
    
    else:
        ampldeltaCorrel = 0.
        sFactordeltaCorrel = 0.
        amploffsetCorrel = 0.
    
    #exit()
    
    residNLstresult = []
    for res in range(len(residNresult)):
        residNLstresult.append(sum(residNresult[res]**2.))

    report_fit(result.params)
    for key, par in result.params.items():
        print('{} {:.3f} +/- {}'.format(key, par.value, par.stderr))
    
    failed = False

    if 'ampl_1' not in correlpair:
        amplitudefit = result.params['ampl_1'].value
    else:
        amplitudefit = amplfixed


    if 'asymmetry_1' not in correlpair:
        asymmetryfit = result.params['asymmetry_1'].value
    else:
        asymmetryfit = asymmetryfixed


    if 'sFactor_1' not in correlpair:
        sizefactor = result.params['sFactor_1'].value
    else:
        sizefactor = amplfixed


    if 'center0_1' not in correlpair:
        center0fit = result.params['center0_1'].value
    else:
        center0fit = center0fixed


    if 'offset_1' not in correlpair:
        offsetfit = result.params['offset_1'].value
    else:
        offsetfit = offsetfixed


    if 'deltad_1' not in correlpair:
        deltadfit = result.params['deltad_1'].value
    else:
        deltadfit = delta0fixed
        

    UndWlLst = []
    Undposlst = []
    
    for i in range(NumOfstripes):
        UndWlLst.append(408.184e-9-.478e-9*1.*(i-(int(NumOfstripes/2))))
        Undposlst.append(k)


    if max(residNLstresult) > 12.:
        print("Gamma was maybe bad")
        maybeBadGamma += 1
        badk.append(k)
    
    # plot the data sets and fits
    """

    ThetaLst = np.linspace(-2000.*10**-6,2000.*10**-6,101)
    wlLst = np.linspace((400.-41.5)*10**-9,(400.+38.5)*10**-9,100)
    #wlLst = np.linspace((400.-81.5)*10**-9,(400.+78.5)*10**-9,100)
    Undu2Dspectrum = []
    p1=[]
    for j in range(len(ThetaLst)):
        Unduspectrum = []

        for i in range(len(wlLst)):
            lambdaNull = Undu.lambda0fit(center0fit,(ThetaLst[j]*sizefactor-offsetfit)/sizefactor,lambdaU)
            deltaL = wlLst[i]-lambdaNull
            WLofPixel = wlLst[i]

            Unduspectrum.append(amplitudefit*Undu.analyticalspecPhase((ThetaLst[j]*sizefactor-offsetfit)/sizefactor,deltaL,center0fit,WLofPixel,deltadfit))

            if (wlLst[i])>min(UndWlLst) and (wlLst[i])<max(UndWlLst):
                if ThetaLst[j]>(-.004+asymmetryfit-offsetfit)/sizefactor and ThetaLst[j]<(+.004+asymmetryfit-offsetfit)/sizefactor:
                    #p1.append([i,j])
                    p1.append([wlLst[i]*10**9,ThetaLst[j]*10**3])

        Unduspectrum = np.transpose(Unduspectrum)

        SourcelstIntensity = abs(np.array(Unduspectrum)**2.)

        Undu2Dspectrum.append(SourcelstIntensity)

    wlLstplot = []
    for i in range(6):
        wlLstplot.append(int(10**9*(i*(wlLst[-1] - wlLst[0])/5+wlLst[0])))
        
    Thetaplot = []
    for i in range(6):
        Thetaplot.append(10**3*(i*(ThetaLst[-1] - ThetaLst[0])/5+ThetaLst[0]))

    fig2 = plt.figure()
    plt.plot([p1[0][0],p1[-1][0],p1[-1][0],p1[0][0],p1[0][0]],[p1[0][1],p1[0][1],p1[-1][1],p1[-1][1],p1[0][1]],'r')
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Angle [mm]')
    plt.imshow(Undu2Dspectrum,aspect=20,extent=(wlLstplot[0],wlLstplot[-1]+1,Thetaplot[0],Thetaplot[-1]))

    #plt.colorbar()

    figA = plt.figure("Diffraction")
    axA = figA.add_subplot(111)
    figB = plt.figure("Sources")
    axB = figB.add_subplot(111)
    figC = plt.figure("FFT")
    axC = figC.add_subplot(111)
    
    for Strpidx in range(NumOfstripes):
        #Strpidx = 2
        

        y_fit = diffraction_datasetfinal(result.params, Strpidx, zetaprime2[0:1000])

       
        WLofPixel = 408.184e-9 -.478e-9*1.*(Strpidx-(int(NumOfstripes/2)))
        Sourcelst=[]
        for ii in range(len(zetaprime2)):

            lambdaNull = Undu.lambda0fit(center0fit,(zetaprime2[ii]-offsetfit)/sizefactor,lambdaU)
            deltaL = WLofPixel-lambdaNull
            #print("deltaL Stripe: ",deltaL)
            if (zetaprime2[ii]-asymmetryfit > -a2 and zetaprime2[ii]-asymmetryfit < a2):
                if deltaON:
                    Sourcelst.append(amplitudefit*Undu.analyticalspecPhase((zetaprime2[ii]-offsetfit)/sizefactor,deltaL,center0fit,WLofPixel,deltadfit))
                else:
                    Sourcelst.append(amplitudefit*Undu.analyticalspec((zetaprime2[ii]-offsetfit)/sizefactor,UndWlLst[Strpidx]))
            else:
                Sourcelst.append(0.)

        Sourcelst = np.transpose(Sourcelst)

        SourcelstIntensity1D = 55000.*abs(np.array(Sourcelst)**2.)

        residualbeforeplot = sum(((np.array(y_fit-data[Strpidx, :])[0:1000])**2.)/(1000.*datasigma[Strpidx][0:1000]**2.))
        print(residualbeforeplot)

        innerresidual = sum(((np.array(y_fit-data[Strpidx, :])[400:600])**2.)/(200.*datasigma[Strpidx][400:600]**2.))
        print(innerresidual)

        maxofdata = np.around(max(data[Strpidx, :]), decimals=-1)
        maxofSourcelstIntensity1D = np.around(max(SourcelstIntensity1D), decimals=-1)
        maxofplt = max([maxofdata,maxofSourcelstIntensity1D])

        fntsze = 13


        axA.annotate(r'$\Delta d + d = %d $' % k, xy=(0, maxofplt), fontsize=fntsze)
        axA.plot(1000.*zetaprime2[0:1000], data[Strpidx, :], '-', 1000.*zetaprime2[0:1000], y_fit, '-')
        axA.set_xlim(-7, 7)
        axA.set_ylim(-10, maxofplt+10)
        axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
        axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
        axA.tick_params(axis='both', which='major', labelsize=fntsze)
        
        axB.annotate(r'$\Delta d + d = %d $' % k, xy=(0, maxofplt), fontsize=fntsze)
        axB.plot(1000.*zetaprime2[0:1000], SourcelstIntensity1D, '-')
        axB.set_xlim(-7, 7)
        axB.set_ylim(-10, maxofplt+10)
        axB.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
        axB.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
        axB.tick_params(axis='both', which='major', labelsize=fntsze)

    plt.show()
    #exit()
    #"""

    k = k + 1# + 78

    interpolatedpicture = np.array(np.transpose(interpolatedpicture))
    interpolatedsigmapicture = np.array(np.transpose(interpolatedsigmapicture))

    #"""
    if DoNotSave == 0:
        with open(DataPath + "\SharedP\\Undupos.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow(Undposlst)

        with open(DataPath + "\SharedP\AmpldeltaCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([ampldeltaCorrel])

        with open(DataPath + "\SharedP\sFactordeltaCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([sFactordeltaCorrel])

        with open(DataPath + "\SharedP\AmploffsetCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([amploffsetCorrel])

        with open(DataPath + "\SharedP\Amplitudes.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([amplitudefit])

        with open(DataPath + "\SharedP\center0s.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([center0fit])

        with open(DataPath + "\SharedP\sFactors.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([sizefactor])

        with open(DataPath + "\SharedP\Bsymmetries.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([asymmetryfit])

        with open(DataPath + "\SharedP\Offsets.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([offsetfit])

        with open(DataPath + "\SharedP\deltads.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([deltadfit])
    #"""