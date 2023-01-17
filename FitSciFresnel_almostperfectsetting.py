import matplotlib.pyplot as plt

import UndulatorDiffractionOneDimensional as UnDi

from numpy import transpose as nptp
import numpy as np

import cv2
import os , sys
import re

import csv

from natsort import natsorted
from scipy.optimize import curve_fit
from lmfit import minimize, Parameters, report_fit, Minimizer

Lu =.52
LambU=.08

z0,z1,z2,zc,f = 10.,3.03,.32,.68,1.

zp = 3.0/1. # should be z1+z2
zpp = zc*f/(f-zc) + zp
zppp = z0 + zpp

ratio = 1000. / 600.#606. #how many pixel are aperture
apy = .004

zetaprime2 = np.linspace(-ratio*apy,ratio*apy,1000)
fixedsize = .245
k0 = 40
k = k0


def diffraction_dataset(params,i,yy):
    Amplfit = params['Amplfit_%i' % (i+1)].value
    Ytran   = params['Ytran_%i' % (i+1)].value
    gammafit  = params['gamma_%i' % (i+1)].value
    delta0fit  = params['delta0_%i' % (i+1)].value
    print(delta0fit+i*.001)

    return UnDi.diffractionintensity(yy,Ytran,i,Amplfit,z0,408.3e-9,delta0fit+i*.001,gammafit)

def objective_diffraction(params, yy, data,weights):
    ndata,nx = data.shape
    
    resid = 0.0*data[:]
    
    for Strpidx in range(ndata):
        resid[Strpidx] = (data[Strpidx] - diffraction_dataset(params,Strpidx,yy))*weights[Strpidx]**2.
        
    return np.array(resid.flatten())

arg = sys.argv

tiflist = []

filelist = os.listdir(arg[1])
for i in range(len(filelist)):
    if filelist[i][-4:] == "tiff" or filelist[i][-3:] == "tif":
        tiflist.append(filelist[i])
filelist = natsorted(tiflist)

Commonpath = str(os.path.dirname(os.path.abspath(__file__)))

print(Commonpath)

DataPath = arg[1]

print(DataPath)

slashes = [i for i,ltr in enumerate(DataPath) if ltr == "\\"]
DataPathBase = DataPath[0:slashes[-2]]

print(DataPathBase)

#datarangelow = 1804 + 550

Backgroundimage = cv2.imread(DataPathBase + "/BG/raw_shortExp/File1.tif",cv2.IMREAD_UNCHANGED)
ListOfAllMatrices = []
ListOfAllSigmaMatrices = []
NumOfStripes = 1
NumOfSeries = 200

WidthOfStripe = 5#int((2304-datarangelow)/NumOfStripes)
print(WidthOfStripe)

MatrixOfSlit = []

maybeBadGamma = 0

targetfiles = ["Undupos.csv",
               "Amplitudes.csv",
               "LambdaRfits.csv",
               "dispersions.csv",
               "Ytrans.csv",
               "Offsets.csv",
               "deltads.csv",
               "deltadshisto.csv",
               "amplcenterCorrel.csv",
               "AmpldeltaCorrel.csv"]

badk = []

DoNotSave = 0

while k < int(len(filelist)/4):
    #Guessvalues = [4.4,-.0004,delta0+k*.001]
    if k == k0:
        checklist = os.listdir(DataPath + "\SharedP")
        maybeBadGamma = 0
        for ckid in range(len(targetfiles)):
            if targetfiles[ckid] in checklist:
                os.remove(DataPath + "\SharedP" + "//" + targetfiles[ckid])

    SigmaOneLineEachFromTenPicturesAtOnePosition=[]

    SeriesOneLineEachFromFourPicturesAtOnePosition=[]
    SeriesOneLineEachFromFourPicturesAtOnePositionForSigma=[]

    for ii in range(NumOfSeries):
        OneLineEachFromFourPicturesAtOnePosition = []
        OneLineEachFromFourPicturesAtOnePositionForSigma=[]
        for i in range(4):
            noti = i
            PictureMinusBG = cv2.imread(DataPath + "/" + filelist[i+ii*4+k*4],cv2.IMREAD_UNCHANGED).astype(int)-Backgroundimage
            PictureMinusBGForSigma = cv2.imread(DataPath + "/" + filelist[i+ii*4+k*4],cv2.IMREAD_UNCHANGED).astype(int)-Backgroundimage

            PictureMinusBG = np.float32(PictureMinusBG)
            PictureMinusBGForSigma = np.float32(PictureMinusBGForSigma)
            print(len(PictureMinusBG[0]))

            PictureMinusBG = nptp(nptp(PictureMinusBG)[(2304-WidthOfStripe):])
            print(len(PictureMinusBG[0]))
            PictureMinusBGForSigma = nptp(nptp(PictureMinusBGForSigma)[(2304-WidthOfStripe):])

            OneLineEachFromFourPicturesAtOnePosition.append(PictureMinusBG)
            OneLineEachFromFourPicturesAtOnePositionForSigma.append(PictureMinusBGForSigma)
        
        TenPixelFromTenPicturesTimesLineAtOnePosition = sum(OneLineEachFromFourPicturesAtOnePosition)/4.
        TenPixelFromTenPicturesTimesLineAtOnePositionForSigma = sum(OneLineEachFromFourPicturesAtOnePositionForSigma)/4.
        
        SeriesOneLineEachFromFourPicturesAtOnePosition.append(TenPixelFromTenPicturesTimesLineAtOnePosition)
        SeriesOneLineEachFromFourPicturesAtOnePositionForSigma.append(TenPixelFromTenPicturesTimesLineAtOnePositionForSigma)

    print(len(SeriesOneLineEachFromFourPicturesAtOnePosition))
    print(len(SeriesOneLineEachFromFourPicturesAtOnePosition[0]))
    print(len(SeriesOneLineEachFromFourPicturesAtOnePosition[0][0]))

    interpolatedpicture = []
    interpolatedsigmapicture = []
    
    for ii in range(NumOfSeries):
        stripe = SeriesOneLineEachFromFourPicturesAtOnePosition[ii]#nptp(nptp(SeriesOneLineEachFromFourPicturesAtOnePosition[ii])[WidthOfStripe*0:WidthOfStripe+WidthOfStripe*0])
        stripesigma = SeriesOneLineEachFromFourPicturesAtOnePositionForSigma[ii]#nptp(nptp(SeriesOneLineEachFromFourPicturesAtOnePositionForSigma[ii])[WidthOfStripe*0:WidthOfStripe+WidthOfStripe*0])
        
        #cv2.imshow("imahe",TenPixelFromTenPicturesTimesLineAtOnePosition)
        #cv2.waitKey(0)

        stripe = cv2.resize(stripe,(1,len(TenPixelFromTenPicturesTimesLineAtOnePosition)),interpolation=cv2.INTER_AREA) #1 zu 1000
        stripesigma = cv2.resize(stripesigma,(1,len(SeriesOneLineEachFromFourPicturesAtOnePositionForSigma)),interpolation=cv2.INTER_AREA) #1 zu 1000

        interpolatedpicture.append(nptp(stripe)[0])
        interpolatedsigmapicture.append(nptp(stripesigma))

        data2D = np.transpose(stripe)[0]/200.
        datasigma2D = np.transpose(stripesigma)[0]/200.
        
        print(stripe[500])

    data = np.array(interpolatedpicture)
    datasigma = np.array(interpolatedsigmapicture)


    dataLst = []
    for Strpidx in range(NumOfSeries):
        dataLst.append(data[Strpidx][500])

    fntsze = 13
    figA = plt.figure("Diffraction")
    axA = figA.add_subplot(111)
    axA.plot(range(len(dataLst)), dataLst, '-')
    axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
    axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
    axA.tick_params(axis='both', which='major', labelsize=fntsze)

    plt.show()


    print(k)
    print((NumOfSeries/5.))
    DoNotSave = 1
    
    Guessvalueslmfit = [4.4,-.00001,351.,.210+k*.001]

    fit_params = Parameters()
    for iy,y, in enumerate(data):
        fit_params.add('Amplfit_%i' % (iy+1),value=Guessvalueslmfit[0], min=3.5, max=5.)
        fit_params.add('Ytran_%i' % (iy+1),value=Guessvalueslmfit[1], min = Guessvalueslmfit[1]-.0001, max=Guessvalueslmfit[1]+.0001)
        fit_params.add('gamma_%i' % (iy+1),value=Guessvalueslmfit[2],min = Guessvalueslmfit[2]-1., max = Guessvalueslmfit[2]+1.)
        fit_params.add('delta0_%i' % (iy+1),value=Guessvalueslmfit[3], min = Guessvalueslmfit[3]-.013, max=Guessvalueslmfit[3]+.003)

    for iy in tuple(range(NumOfSeries+1)[2:]):
        fit_params['Amplfit_%i' % iy].expr='Amplfit_1'
        
    for iy in tuple(range(NumOfSeries+1)[2:]):
        fit_params['Ytran_%i' % iy].expr='Ytran_1'
        
    for iy in tuple(range(NumOfSeries+1)[2:]):
        fit_params['gamma_%i' % iy].expr='gamma_1'
        
    for iy in tuple(range(NumOfSeries+1)[2:]):
        fit_params['delta0_%i' % iy].expr='delta0_1'

    mini = Minimizer(objective_diffraction, fit_params,nan_policy = "raise",scale_covar=False, fcn_args=(zetaprime2,data,1./stripesigma),)

    print(data[0][500])
    print(data[1][500])

    print(len(zetaprime2))
    print(len(datasigma[0]))

    #"""

    result = mini.leastsq(ftol=1.e-3)

    report_fit(result.params)
    
    y_fitLst = []
    dataLst = []

    for Strpidx in range(NumOfSeries):

        WLofPixel = 408.184e-9 -.478e-9*1.#*(Strpidx-(int(NumOfStripes/2)))
        Sourcelst=[]

        Amplfit = result.params['Amplfit_1'].value
        Ytran = result.params['Ytran_1'].value
        gammafit = result.params['gamma_1'].value
        delta0fit = result.params['delta0_1'].value

        Sourcelst = UnDi.plainintensity(zetaprime2,Ytran,i,Amplfit,z0,408.3e-9,delta0fit+.001*Strpidx,gammafit)

        y_fitLst.append((UnDi.diffractionintensitySci(zetaprime2,Amplfit,Ytran,delta0fit+.001*Strpidx,gammafit))[500])
        dataLst.append(data[Strpidx][500])

    fntsze = 13
    figA = plt.figure("Diffraction")
    axA = figA.add_subplot(111)

    #axA.annotate(r'$\Delta d + d = %d $' % k, xy=(0, maxofplt), fontsize=fntsze)
    axA.plot(range(len(dataLst)), dataLst, '-', range(len(dataLst)), y_fitLst, '-')
    #axA.set_xlim(-7, 7)
    #axA.set_ylim(-10, maxofplt+10)
    axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
    axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
    axA.tick_params(axis='both', which='major', labelsize=fntsze)

    plt.show()
        

    for Strpidx in range(NumOfSeries):
        
        figA = plt.figure("Diffraction")
        axA = figA.add_subplot(111)
        figB = plt.figure("Sources")
        axB = figB.add_subplot(111)

        #y_fit = diffraction_dataset(result.params, Strpidx, zetaprime2)
        #y_fit = UnDi.diffractionintensitySci(zetaprime2,*popt)#
        #y_fit = UnDi.diffractionintensitySci_delta(zetaprime2,deltad)#diffraction_dataset(result.params, Strpidx, zetaprime2)

       
        WLofPixel = 408.184e-9 -.478e-9*1.#*(Strpidx-(int(NumOfStripes/2)))
        Sourcelst=[]
        
        Amplfit = result.params['Amplfit_1'].value
        Ytran = result.params['Ytran_1'].value
        gammafit = result.params['gamma_1'].value
        delta0fit = result.params['delta0_1'].value

        
        #deltad = delta0+Strpidx*.001

        Sourcelst = UnDi.plainintensity(zetaprime2,Ytran,i,Amplfit,z0,408.3e-9,delta0fit+.001*Strpidx,gammafit)

        SourcelstIntensity1D = 55000.*abs(np.array(Sourcelst))

        y_fit = UnDi.diffractionintensitySci(zetaprime2,Amplfit,Ytran,delta0fit+.001*Strpidx,gammafit)

        #residualbeforeplot = sum(((np.array(y_fit-data[Strpidx, :])[0:1000])**2.)/(1000.*datasigma[Strpidx][0:1000]**2.))
        #print(residualbeforeplot)

        #innerresidual = sum(((np.array(y_fit-data[Strpidx, :])[400:600])**2.)/(200.*datasigma[Strpidx][400:600]**2.))
        #print(innerresidual)

        #maxofdata = np.around(max(data[Strpidx, :]), decimals=-1)
        #maxofSourcelstIntensity1D = np.around(max(SourcelstIntensity1D), decimals=-1)
        #maxofplt = max([maxofdata,maxofSourcelstIntensity1D])

        fntsze = 13

        #axA.annotate(r'$\Delta d + d = %d $' % k, xy=(0, maxofplt), fontsize=fntsze)
        axA.plot(1000.*zetaprime2[0:1000], data[Strpidx, :], '-', 1000.*zetaprime2, y_fit, '-')
        axA.set_xlim(-7, 7)
        #axA.set_ylim(-10, maxofplt+10)
        axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
        axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
        axA.tick_params(axis='both', which='major', labelsize=fntsze)
        
        #axB.annotate(r'$\Delta d + d = %d $' % k, xy=(0, maxofplt), fontsize=fntsze)
        axB.plot(1000.*zetaprime2, SourcelstIntensity1D, '-')
        axB.set_xlim(-7, 7)
        #axB.set_ylim(-10, maxofplt+10)
        axB.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
        axB.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
        axB.tick_params(axis='both', which='major', labelsize=fntsze)

        plt.show()
    #"""
    k = k + 10

    if DoNotSave == 0:
        with open(DataPath + "\SharedP\\Amplfit.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([Amplfit])

        with open(DataPath + "\SharedP\Ytran.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([Ytran])

        #"""
        with open(DataPath + "\SharedP\deltad.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([deltad])
        #"""
        
        with open(DataPath + "\SharedP\\uncertAmpl.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([uncertAmpl])

        with open(DataPath + "\SharedP\\uncertYtra.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([uncertYtra])
        
        with open(DataPath + "\SharedP\\AmplfitYtran.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([AmplfitYtran])
            
        #"""
        with open(DataPath + "\SharedP\Amplfitdeltad.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([Amplfitdeltad])

        with open(DataPath + "\SharedP\Ytrandeltad.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([Ytrandeltad])
        #"""