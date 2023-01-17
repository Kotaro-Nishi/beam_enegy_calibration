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

import random

import imageio

ratio = 1000. / 606#606. #how many pixel are aperture
#UnDi.apy = .004
UnDi.CMOS = .0065/2.
UnDi.apy = .008/2.
#UnDi.apy = UnDi.CMOS

#0.00515
#0.00165

zetaprime2 = np.linspace(-ratio*UnDi.apy,ratio*UnDi.apy,1000)
CMOSprime2 = np.linspace(-UnDi.CMOS,UnDi.CMOS,1000)
fixedsize = .245


NumOfStripes = 200      #How many Wavelengths are scanned

NumOfSeries = 825       #How many of the datapoints of the interval are used for the fit.

LengthOfInterval = 825 #How many datapoints are included for each fit.

             #0.004808759252942803 Wert aus Auswertung
             #1.5629889792630894e-07
             #1315.3024745554976
             #1.0199041435307357
dispersion = 0.004808759252942803e-9 #Former 0.0048085643228767475e-9
pixelof404nm = 1315.3024745554976 #Former: 1243.8551272671261
CalibWl = 404.6565e-9

#(407.7837-404.6565)

global listOfboots
listOfboots = []

Yscreenfit = .3e-4#-.5e-4
LambdaRfit = 400e-9

def per_iteration(pars, iteration, resid, *args, **kws):
    print(" resid ", sum(resid))#sum(abs(np.array(resid)))/1000.)
    #print(" resid ", resid)
    print("Amplfit: " + str(pars['Amplfit_%i'% listOfboots[0]].value))
    print("bfit: " + str(pars['bfit_%i'% listOfboots[0]].value))
    #print("zpfit: " + str(pars['zpfit_%i'% listOfboots[0]].value))
    print('Ytran: ' + str(pars['Ytran_%i'% listOfboots[0]].value))
    print("gammafit: " + str(pars['gamma_%i'% listOfboots[0]].value))
    print("delta0fit: " + str(pars['delta0_%i'% listOfboots[0]].value))
    #print("Yscreenfit: " + str(pars['Yscreenfit_1'].value))
    #print("LambdaRfit: " + str(pars['LambdaRfit_1'].value))

    Amplfit = pars['Amplfit_%i'% listOfboots[0]].value
    Ytran = pars['Ytran_%i'% listOfboots[0]].value
    gammafit = pars['gamma_%i'% listOfboots[0]].value
    delta0fit = pars['delta0_%i'% listOfboots[0]].value
    bfit = pars['bfit_%i'% listOfboots[0]].value
    #zpfit = pars['zpfit_%i'% listOfboots[0]].value
    #Yscreenfit = pars['Yscreenfit_1'].value
    #LambdaRfit = pars['LambdaRfit_1'].value

    UnDi.SingleAmplitudeDone = 0
    UnDi.AmplXY = []
    

def diffraction_dataset(params,i,yy):
    Amplfit = params['Amplfit_%i' % (i)].value
    Ytran   = params['Ytran_%i' % (i)].value
    gammafit  = params['gamma_%i' % (i)].value
    delta0fit  = params['delta0_%i' % (i)].value
    bfit  = params['bfit_%i' % (i)].value
    #zpfit  = params['zpfit_%i' % (i)].value

    #Yscreenfit  = params['Yscreenfit_%i' % (i+1)].value
    #LambdaRfit  = params['LambdaRfit_%i' % (i+1)].value
    #print(delta0fit+i*.001)

    Diff_intensity = UnDi.diffractionintensity(yy,Ytran,i,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+i*.001,gammafit,Amplfit,bfit,5.0207,Yscreenfit)


    return Diff_intensity

def objective_diffraction(params, yy, data,weights):
    ndata,nx = data.shape
    resid = 0.0*data[:]
    
    for Strpidx in range(ndata):
        if Strpidx > 0:
            UnDi.SingleAmplitudeDone = 1
        resid[Strpidx] = (diffraction_dataset(params,listOfboots[Strpidx],yy) - data[Strpidx])*weights[Strpidx]#*(1.+50000.*yy**2.)#**2.
        
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
print(DataPathBase + "/BG/raw_shortExp/File1.tif")

WidthOfStripe = 5#int((2304-datarangelow)/NumOfStripes)
print(WidthOfStripe)

MatrixOfSlit = []

maybeBadGamma = 0

targetfiles = [#"residimg.csv",
               "Amplfit.csv",
               "Ytran.csv",
               "gammafit.csv",
               "delta0fit.csv",
               "bfit.csv",
               "Amplfitstderr.csv",
               "Ytranstderr.csv",
               "gammafitstderr.csv",
               "delta0fitstderr.csv",
               "bfitstderr.csv",
               "amplYtranCorrel.csv",
               "amplgammafitCorrel.csv",
               "ampldelta0fitCorrel.csv",
               "amplbfitCorrel.csv",
               #"amplzpfitCorrel.csv",
               "YtrangammafitCorrel.csv",
               "Ytrandelta0fitCorrel.csv",
               "YtranbfitCorrel.csv",
               "gammadelta0fitCorrel.csv",
               "gammabfitCorrel.csv",
               #"bfitzpfitCorrel.csv",
               #"delta0fitzpfitCorrel.csv",
               #"gammazpfitCorrel.csv",
               #"YtranzpfitCorrel.csv",
               #"zpfit.csv",
               #"zpfitstderr.csv",
               "delta0fitbfitCorrel.csv"]

badk = []

DoNotSave = 0

Amplfit = -100
Ytran = -100
gammafit = -100
delta0fit = -100
bfit = -100

kstep = 1

k0 = 000
k = k0

startpx = 950

SeriesListPictureMinusBGWidthOfStripe4Mean=[]
SeriesListPictureMinusBGWidthOfStripeSigma=[]

for ii in range(LengthOfInterval):

    chosenpx = startpx + int(NumOfStripes*WidthOfStripe/2)# + ii
    #print(chosenpx)
    #print((chosenpx-int(NumOfStripes*WidthOfStripe/2))-(chosenpx+int(NumOfStripes*WidthOfStripe/2)))
    #print(chosenpx-int(NumOfStripes*WidthOfStripe/2))
    #print(chosenpx+int(NumOfStripes*WidthOfStripe/2))

    ListPictureMinusBGWidthOfStripe = []
    ListPictureMinusBGWidthOfStripeForSigma=[]
    print(filelist[ii*4])
    for i in range(4):
        noti = i
        PictureMinusBG = cv2.imread(DataPath + "/" + filelist[i+ii*4],cv2.IMREAD_UNCHANGED).astype(int)-Backgroundimage
        PictureMinusBGForSigma = cv2.imread(DataPath + "/" + filelist[i+ii*4],cv2.IMREAD_UNCHANGED).astype(int)-Backgroundimage
        #PictureMinusBG = cv2.flip(PictureMinusBG, 0)
        #PictureMinusBGForSigma = cv2.flip(PictureMinusBGForSigma, 0)

        PictureMinusBG = np.float32(PictureMinusBG)
        PictureMinusBGForSigma = np.float32(PictureMinusBGForSigma)
        
        PictureMinusBGWidthOfStripe = nptp(nptp(PictureMinusBG)[(chosenpx-int(NumOfStripes*WidthOfStripe/2)):chosenpx+int(NumOfStripes*WidthOfStripe/2)])
        #print(len(PictureMinusBGWidthOfStripe[0]))
        #print(len(nptp(PictureMinusBG)))

        PictureMinusBGForSigmaWidthOfStripe = nptp(nptp(PictureMinusBGForSigma)[(chosenpx-int(NumOfStripes*WidthOfStripe/2)):chosenpx+int(NumOfStripes*WidthOfStripe/2)])
        
        ListPictureMinusBGWidthOfStripe.append(PictureMinusBGWidthOfStripe)
        ListPictureMinusBGWidthOfStripeForSigma.append(PictureMinusBGForSigmaWidthOfStripe)


    ListPictureMinusBGWidthOfStripeSigma = np.std(ListPictureMinusBGWidthOfStripeForSigma, axis=0)

    ListPictureMinusBGWidthOfStripe4Mean = sum(ListPictureMinusBGWidthOfStripe)/4.
    
    SeriesListPictureMinusBGWidthOfStripe4Mean.append(ListPictureMinusBGWidthOfStripe4Mean)
    SeriesListPictureMinusBGWidthOfStripeSigma.append(ListPictureMinusBGWidthOfStripeSigma)

checkstart = 1

random.seed(10000)
WL_lst=[]

for Saves in range(0,NumOfStripes):
    print("Saves " + str(Saves))

    chosenpx = startpx + WidthOfStripe/2. + Saves*WidthOfStripe
    WLofPixel = (chosenpx - pixelof404nm)*dispersion + CalibWl

    WL_lst.append(WLofPixel)

    zetaprime2Booted = []
    dataBooted = []
    datasigmaBooted = []
    listOfboots = []
    UnDi.SingleAmplitudeDone = 0
    for i in range(NumOfSeries):
        listOfboots.append(i)

    
    for SIt in range(len(listOfboots)):
        if listOfboots.count(listOfboots[SIt]) == 1:
            FirstSingularitem = SIt#listOfboots[SIt]
            break

    residLst = []

    interpolatedpicture = []
    interpolatedsigmapicture = []

    for ii in range(LengthOfInterval):
        stripe = SeriesListPictureMinusBGWidthOfStripe4Mean[ii]
        stripesigma = SeriesListPictureMinusBGWidthOfStripeSigma[ii]
        
        stripe = cv2.resize(stripe,(NumOfStripes,len(ListPictureMinusBGWidthOfStripe4Mean)),interpolation=cv2.INTER_AREA) #1 zu 1000

        ElectronSource = np.array(IFunc.ElectronSourcefunctionSet(zetaprime2,.00004))  #.00004

        stripe = IFunc.convolveplain(nptp(stripe)[Saves],ElectronSource).real

        stripejjj=[]
        for jjj in range(len(stripe)):
            stripejjj.append(stripe[jjj])
            
        stripe = nptp(np.array([stripejjj])*3000000.)

        stripesigma = cv2.resize(stripesigma,(NumOfStripes,len(ListPictureMinusBGWidthOfStripeSigma)),interpolation=cv2.INTER_AREA) #1 zu 1000

        interpolatedpicture.append(nptp(stripe)[0])
        interpolatedsigmapicture.append(nptp(stripesigma)[Saves])
        
    data = np.array(interpolatedpicture)
    datasigma = np.array(interpolatedsigmapicture)

    for Boot in range(NumOfSeries):
        dataBooted.append(data[listOfboots[Boot]])
        datasigmaBooted.append(datasigma[listOfboots[Boot]])
    dataBooted = np.array(dataBooted)
    datasigmaBooted = np.array(datasigmaBooted)

    dataLst = []
    datasigLst = []
    for Strpidx in range(NumOfSeries):
        dataLst.append(dataBooted[Strpidx][500])
        datasigLst.append(datasigmaBooted[Strpidx][500])
    
    """
    fntsze = 13
    figA = plt.figure("Diffraction")
    axA = figA.add_subplot(111)
    axA.plot(listOfboots, dataLst, 'b.')
    axA.plot(listOfboots, datasigLst, 'g.')
    axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
    axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
    axA.tick_params(axis='both', which='major', labelsize=fntsze)
    plt.show()
    #"""

    y_fitLst = []
    dataLst = []
    
    #"""

    for Strpidx in range(NumOfSeries):

        Sourcelst=[]

        Amplfit = 0.97864585
        Ytran = 1.4483e-04
        gammafit = 352.388261
        delta0fit = 0.20683398
        bfit = 0.51805446
        zpfit = 5.0207
        #Yscreenfit = result.params['Yscreenfit_1'].value
        #LambdaRfit = result.params['LambdaRfit_1'].value

        Sourcelst = UnDi.plainintensity(zetaprime2,Ytran,Strpidx,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+.001*Strpidx,gammafit,Amplfit,bfit,zpfit,Yscreenfit)
        
        #y_fitLst.append(sum(UnDi.diffractionintensity(zetaprime2,Ytran,Strpidx,AUndu.z0,404.3e-9,delta0fit+Strpidx*.001,gammafit,Amplfit,bfit)[475:485]))
        #dataLst.append(sum((dataBooted[Strpidx])[475:485]))
        y_fitLst.append(UnDi.diffractionintensity(zetaprime2,Ytran,Strpidx,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+Strpidx*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit)[480])
        dataLst.append((dataBooted[Strpidx])[480])

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
    
    axis_font = {'fontname':'Arial', 'size':'21.5'}
    numbers_font_Century = {'fontname':'Arial', 'size':'21.5'}
 
    y_fitLst = []
    yyy=[]

    for i in range(len(zetaprime2)):
        yyy.append(zetaprime2[i] - Yscreenfit)
    yyy = np.array(yyy)

    #for Strpidx in range(NumOfSeries):
    
    figA = plt.figure("Diffraction")
    axA = figA.add_subplot(111)
    #figB = plt.figure("Sources")
    #axB = figB.add_subplot(111)
    
    def plot_for_offset(Strpidx):

        Sourcelst=[]
        
        y_fit = UnDi.diffractionintensity(zetaprime2,Ytran,Strpidx,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+Strpidx*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit)

        fntsze = 13
        
        figA, axA = plt.subplots(figsize=(10,5))
        axA.set(title='Diffraction')

        axA.plot(1000.*yyy[0:1000], dataBooted[Strpidx, :], '-', 1000.*yyy, y_fit, '-')
        axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
        axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
        axA.tick_params(axis='both', which='major', labelsize=fntsze)

        figA.canvas.draw()       # draw the canvas, cache the renderer
        image = np.frombuffer(figA.canvas.tostring_rgb(), dtype='uint8')
        image  = image.reshape(figA.canvas.get_width_height()[::-1] + (3,))

        return image
        
    kwargs_write = {'fps':1.0, 'quantizer':'nq'}
    imageio.mimsave('./DiffractionVideo.gif', [plot_for_offset(i) for i in range(825)], fps=8)
    #"""

    
    checkstart = 0
print(WL_lst)
exit()