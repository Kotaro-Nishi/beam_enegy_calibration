import matplotlib.pyplot as plt
import UndulatorDiffractionOneDimensional_correct_XAmpl_InitialFit as UnDi
import InterferenceFunction_more_x as IFunc
import AmplitudeFormula_x as AUndu

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

ratio = 1000. / 606#606. #how many pixel are aperture
UnDi.apy = .0081/2.
zetaxx = np.array([0.])

zetaprime2 = np.linspace(-ratio*UnDi.apy,ratio*UnDi.apy,1000)

NumOfStripes = 1
NumOfSeries = 825
LengthOfInterval = 825

dispersion = 0.004808759252942803e-9 #Former 0.0048085643228767475e-9
pixelof404nm = 1315.3024745554976 #Former: 1243.8551272671261
CalibWl = 404.6565e-9

#(407.7837-404.6565)

chosenpx = 1300

global listOfboots
listOfboots = []

#Yscreenfit = (380)*10**-6
Yscreenfit = (460)*10**-6
LambdaRfit = 412e-9

def per_iteration(pars, iteration, resid, *args, **kws):
    
    print(" resid ", sum(resid))
    print(" residSqr ", sum(resid**2.))

    print("Amplfit: " + str(pars['Amplfit_%i'% listOfboots[0]].value))
    #print("bfit: " + str(pars['bfit_%i'% listOfboots[0]].value))
    print("zpfit: " + str(pars['zpfit_%i'% listOfboots[0]].value))
    print('Ytran: ' + str(pars['Ytran_%i'% listOfboots[0]].value))
    print("gammafit: " + str(pars['gamma_%i'% listOfboots[0]].value))
    print("delta0fit: " + str(pars['delta0_%i'% listOfboots[0]].value))
    print("elevationfit: " + str(pars['elevationfit_%i'% listOfboots[0]].value))
    print("cfit: " + str(pars['cfit_%i'% listOfboots[0]].value))

    Amplfit = pars['Amplfit_%i'% listOfboots[0]].value
    Ytran = pars['Ytran_%i'% listOfboots[0]].value
    gammafit = pars['gamma_%i'% listOfboots[0]].value
    delta0fit = pars['delta0_%i'% listOfboots[0]].value
    #bfit = pars['bfit_%i'% listOfboots[0]].value
    zpfit = pars['zpfit_%i'% listOfboots[0]].value
    elevationfit = pars['elevationfit_%i'% listOfboots[0]].value
    cfit = pars['cfit_%i'% listOfboots[0]].value

    UnDi.SingleAmplitudeDone = 0
    UnDi.AmplXY = []
    

def diffraction_dataset(params,i,yy):
    Amplfit = params['Amplfit_%i' % (i)].value
    Ytran   = params['Ytran_%i' % (i)].value
    gammafit  = params['gamma_%i' % (i)].value
    delta0fit  = params['delta0_%i' % (i)].value
    #bfit  = params['bfit_%i' % (i)].value
    zpfit  = params['zpfit_%i' % (i)].value
    elevationfit  = params['elevationfit_%i' % (i)].value
    cfit = params['cfit_%i' % (i)].value
    bfit = 0
    Diff_intensity = UnDi.diffractionintensityInitialFit_X(yy,zetaxx,Ytran,i,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+i*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit,elevationfit,whichdata,cfit)


    return Diff_intensity

def objective_diffraction(params, yy, data,weights):
    ndata,nx = data.shape
    
    resid = 0.0*data[:]

    for Strpidx in range(ndata):
        if Strpidx > 0:
            UnDi.SingleAmplitudeDone = 1
        diffractionresult = diffraction_dataset(params,listOfboots[Strpidx],yy)
        #1.5 is std of camera structure
        #1/.2 is conversion factor
        #1/2. is correction of convolutionsigma
        #.05 makes strict positive
        #1*.2 is conversion factor
        #1/.7 is QE
        resid[Strpidx] = (diffractionresult - data[Strpidx])*(1./(np.array([1.5]*len(weights[0]))+(1./.2)*(1./20.)*np.sqrt(0.05+abs(diffractionresult)*.2/.7)))
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

if DataPath[0:-33] == "D:\Doktor\BeamtimeFeb2020_CommonJuneYtoX\data\Beam20200703_1":
    if DataPath[0:-13] == "D:\Doktor\BeamtimeFeb2020_CommonJuneYtoX\data\Beam20200703_1\Energymeasurement_1":
        whichdata = 1
    if DataPath[0:-13] == "D:\Doktor\BeamtimeFeb2020_CommonJuneYtoX\data\Beam20200703_1\Energymeasurement_3":
        whichdata = 3
    if DataPath[0:-13] == "D:\Doktor\BeamtimeFeb2020_CommonJuneYtoX\data\Beam20200703_1\Energymeasurement_5":
        whichdata = 5
else:
    whichdata = 0

print(whichdata)

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

badk = []

DoNotSave = 1

Amplfit = -100
Ytran = -100
gammafit = -100
delta0fit = -100
#bfit = -100
bfit = 0
kstep = 1

FirstWL = (0 - pixelof404nm)*dispersion + CalibWl
WLofPixel = (chosenpx - pixelof404nm)*dispersion + CalibWl
LastWL = (2304 - pixelof404nm)*dispersion + CalibWl

print("First wavelength: ", FirstWL)
print("Wavelength at which the analysis is done: ", WLofPixel)
print("Last wavelength: ", LastWL)

k0 = 000
k = k0
SeriesListPictureMinusBGWidthOfStripe4Mean=[]
SeriesListPictureMinusBGWidthOfStripeSigma=[]

for ii in range(k0,k0+LengthOfInterval):
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
        
        PictureMinusBGWidthOfStripe = nptp(nptp(PictureMinusBG)[(chosenpx-int(WidthOfStripe/2)):chosenpx+int(WidthOfStripe/2)+1])
        PictureMinusBGForSigmaWidthOfStripe = nptp(nptp(PictureMinusBGForSigma)[(chosenpx-int(WidthOfStripe/2)):chosenpx+int(WidthOfStripe/2)+1])
        
        ListPictureMinusBGWidthOfStripe.append(PictureMinusBGWidthOfStripe)
        ListPictureMinusBGWidthOfStripeForSigma.append(PictureMinusBGForSigmaWidthOfStripe)

    ListPictureMinusBGWidthOfStripeSigma = np.std(ListPictureMinusBGWidthOfStripeForSigma, axis=0) + np.array(([[0]*5]*167)+([[0]*5]*250)+([[0]*5]*206)+([[0]*5]*250)+([[0]*5]*127))

    ListPictureMinusBGWidthOfStripe4Mean = sum(ListPictureMinusBGWidthOfStripe)/4.
    
    SeriesListPictureMinusBGWidthOfStripe4Mean.append(ListPictureMinusBGWidthOfStripe4Mean)
    SeriesListPictureMinusBGWidthOfStripeSigma.append(ListPictureMinusBGWidthOfStripeSigma)

SaveSharedP = "\SharedP_20210401_825_825_seed10000_new_zpfit_x"

print(NumOfSeries)
checkstart = 1

random.seed(10000)
#for i in range(NumOfSeries*6500):
#    random.randint(0,LengthOfInterval-1)

residsqrlst=[]

for Saves in range(0,1):

    #190 + 90
    #190 + 160 #17
    #190 + 270 #28
    print("Saves " + str(Saves))

    zetaprime2Booted = []
    dataBooted = []
    datasigmaBooted = []
    listOfboots = []
    UnDi.SingleAmplitudeDone = 0
    for i in range(NumOfSeries):
        listOfboots=list(range(NumOfSeries))#.append(random.randint(0,LengthOfInterval-1))
    
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

        stripe = cv2.resize(stripe,(1,len(ListPictureMinusBGWidthOfStripe4Mean)),interpolation=cv2.INTER_AREA) #1 zu 1000

        ElectronSource = np.array(IFunc.ElectronSourcefunctionSet(zetaprime2,.00004))  #.00004

        stripe = IFunc.convolveplain(nptp(stripe)[0],ElectronSource).real

        stripejjj=[]
        for jjj in range(len(stripe)):
            stripejjj.append(stripe[jjj])
            
        stripe = nptp(np.array([stripejjj])*3000000.)
        
        stripesigma = cv2.resize(stripesigma,(1,len(ListPictureMinusBGWidthOfStripeSigma)),interpolation=cv2.INTER_AREA) #1 zu 1000

        interpolatedpicture.append(nptp(stripe)[0])
        interpolatedsigmapicture.append(nptp(stripesigma))
        
    data = np.array(interpolatedpicture)
    datasigma = np.array(interpolatedsigmapicture)/np.sqrt(5.0)

    for Boot in range(NumOfSeries):
        dataBooted.append(data[listOfboots[Boot]])
        datasigmaBooted.append(datasigma[listOfboots[Boot]])
    dataBooted = np.array(dataBooted)
    datasigmaBooted = np.array(datasigmaBooted)

    dataLst = []
    for Strpidx in range(NumOfSeries):
        dataLst.append(dataBooted[Strpidx][500])
    
    """
    fntsze = 13
    figA = plt.figure("Diffraction")
    axA = figA.add_subplot(111)
    axA.plot(listOfboots, dataLst, 'b.')
    axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
    axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
    axA.tick_params(axis='both', which='major', labelsize=fntsze)
    plt.show()
    #"""

    #real gemessen 4.3m

    Guessvalueslmfit = [1.4,.0001,352.4,.190+.01+k*.001,.6,5.2,11.,1.]#,.35e-4]

    fit_params = Parameters()

    for iy,y, in enumerate(dataBooted):
        fit_params.add('Amplfit_%i' % (listOfboots[iy]),value=Guessvalueslmfit[0])#, max=2.5, min=1.23)#, min = Guessvalueslmfit[0]-1.1, max=Guessvalueslmfit[0]+1.1)
        fit_params.add('Ytran_%i' % (listOfboots[iy]),value=Guessvalueslmfit[1])#, min = Guessvalueslmfit[1]-.0001, max=Guessvalueslmfit[1]+.0001)
        fit_params.add('gamma_%i' % (listOfboots[iy]),value=Guessvalueslmfit[2])
        fit_params.add('delta0_%i' % (listOfboots[iy]),value=Guessvalueslmfit[3])#, min = Guessvalueslmfit[4]-.1, max=Guessvalueslmfit[4]+.1)
        #fit_params.add('bfit_%i' % (listOfboots[iy]),value=Guessvalueslmfit[4])#, min = Guessvalueslmfit[4]-.1, max=Guessvalueslmfit[4]+.1)
        fit_params.add('zpfit_%i' % (listOfboots[iy]),value=Guessvalueslmfit[5])#, min = Guessvalueslmfit[4]-.1, max=Guessvalueslmfit[4]+.1)
        fit_params.add('elevationfit_%i' % (listOfboots[iy]),value=Guessvalueslmfit[6])#, min = Guessvalueslmfit[4]-.1, max=Guessvalueslmfit[4]+.1)
        fit_params.add('cfit_%i' % (listOfboots[iy]),value=Guessvalueslmfit[7])#, min = Guessvalueslmfit[4]-.1, max=Guessvalueslmfit[4]+.1)

    listOfbootsnoFirst = [x for x in listOfboots if x != listOfboots[FirstSingularitem]]
    for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
        fit_params['Amplfit_%i' % listOfbootsnoFirst[iy]].expr='Amplfit_%i'% listOfboots[FirstSingularitem]
        
    for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
        fit_params['Ytran_%i' % listOfbootsnoFirst[iy]].expr='Ytran_%i'% listOfboots[FirstSingularitem]
        
    for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
        fit_params['gamma_%i' % listOfbootsnoFirst[iy]].expr='gamma_%i'% listOfboots[FirstSingularitem]
        
    for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
        fit_params['delta0_%i' % listOfbootsnoFirst[iy]].expr='delta0_%i'% listOfboots[FirstSingularitem]

    #for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
    #    fit_params['bfit_%i' % listOfbootsnoFirst[iy]].expr='bfit_%i'% listOfboots[FirstSingularitem]

    for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
        fit_params['zpfit_%i' % listOfbootsnoFirst[iy]].expr='zpfit_%i'% listOfboots[FirstSingularitem]

    for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
        fit_params['elevationfit_%i' % listOfbootsnoFirst[iy]].expr='elevationfit_%i'% listOfboots[FirstSingularitem]

    for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
        fit_params['cfit_%i' % listOfbootsnoFirst[iy]].expr='cfit_%i'% listOfboots[FirstSingularitem]

    mini = Minimizer(objective_diffraction, fit_params,nan_policy = "raise",scale_covar=False,iter_cb=per_iteration, fcn_args=(zetaprime2,dataBooted,1./datasigmaBooted),)
    
    print(len(zetaprime2))
    print(len(datasigmaBooted[0]))

    #try:
    result = mini.leastsq(ftol=1.e-2)
    #except RuntimeError:
    #    badSample+=1
    #    continue

    report_fit(result.params,show_correl=True, min_correl=0.01)
    
    #residLst.append(result.residual)
    residLst = result.residual
    
    residsqrlst.append(np.sum(residLst**2))

    #"""
    for nicestfit in range(1):
        y_fitLst = []
        dataLst = []
        everynth = 1
        for Strpidx in range(int(NumOfSeries/everynth)):

            Amplfit = result.params['Amplfit_%i'% listOfboots[FirstSingularitem]].value
            Ytran = result.params['Ytran_%i'% listOfboots[FirstSingularitem]].value
            gammafit = result.params['gamma_%i'% listOfboots[FirstSingularitem]].value
            delta0fit = result.params['delta0_%i'% listOfboots[FirstSingularitem]].value
            #bfit = result.params['bfit_%i'% listOfboots[FirstSingularitem]].value
            bfit = 0
            zpfit = result.params['zpfit_%i'% listOfboots[FirstSingularitem]].value
            elevationfit = result.params['elevationfit_%i'% listOfboots[FirstSingularitem]].value
            cfit = result.params['cfit_%i'% listOfboots[FirstSingularitem]].value
            #Yscreenfit = result.params['Yscreenfit_%i'% listOfboots[FirstSingularitem]].value
            #LambdaRfit = result.params['LambdaRfit_1'].value

            y_fitLst.append(UnDi.diffractionintensityInitialFit_X(zetaprime2,zetaxx,Ytran,listOfboots[Strpidx*everynth],AUndu.z0,WLofPixel,LambdaRfit,delta0fit+listOfboots[Strpidx*everynth]*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit,elevationfit,whichdata,cfit)[500+nicestfit*500])
            dataLst.append((dataBooted[Strpidx*everynth])[500+nicestfit*500])

        fntsze = 13
        figA = plt.figure("Diffraction")
        axA = figA.add_subplot(111)
        
        residLst = result.residual

        #axA.annotate(r'$\Delta d + d = %d $' % k, xy=(0, maxofplt), fontsize=fntsze)
        #axA.plot(range(len(dataLst)), dataLst, '-', range(len(dataLst)), y_fitLst, '-')
        axA.plot(range(len(dataLst)), dataLst, 'b.')
        axA.plot(range(len(y_fitLst)), y_fitLst,'r.')
        
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
        yyy.append(zetaprime2[i])# - Yscreenfit)
    yyy = np.array(yyy)

    for Strpidx in range(int(NumOfSeries/everynth)):
        if 1:#Strpidx == 100 or Strpidx == 225 or Strpidx == 350 or Strpidx == 475:#Strpidx > 600 and Strpidx > 610:

            #y_fit = diffraction_dataset(result.params, Strpidx, zetaprime2)
            #y_fit = UnDi.diffractionintensitySci(zetaprime2,*popt)#
            #y_fit = UnDi.diffractionintensitySci_delta(zetaprime2,deltad)#diffraction_dataset(result.params, Strpidx, zetaprime2)

            Sourcelst=[]
            
            Amplfit = result.params['Amplfit_%i'% listOfboots[FirstSingularitem]].value
            Ytran = result.params['Ytran_%i'% listOfboots[FirstSingularitem]].value
            gammafit = result.params['gamma_%i'% listOfboots[FirstSingularitem]].value
            delta0fit = result.params['delta0_%i'% listOfboots[FirstSingularitem]].value
            #bfit = result.params['bfit_%i'% listOfboots[FirstSingularitem]].value
            bfit = 0
            zpfit = result.params['zpfit_%i'% listOfboots[FirstSingularitem]].value
            elevationfit = result.params['elevationfit_%i'% listOfboots[FirstSingularitem]].value
            cfit = result.params['cfit_%i'% listOfboots[FirstSingularitem]].value
            #Yscreenfit = result.params['Yscreenfit_%i'% listOfboots[FirstSingularitem]].value
            #LambdaRfit = result.params['LambdaRfit_1'].value        

            #deltad = delta0+Strpidx*.001

            Sourcelst = UnDi.plainintensity(zetaprime2,zetaxx,Ytran,listOfboots[Strpidx*everynth],AUndu.z0,WLofPixel,LambdaRfit,delta0fit+.001*listOfboots[Strpidx*everynth],gammafit,Amplfit,bfit,zpfit,Yscreenfit,elevationfit)

            SourcelstIntensity1D = 55000.*abs(np.array(Sourcelst))

            y_fit = UnDi.diffractionintensityInitialFit_X(zetaprime2,zetaxx,Ytran,listOfboots[Strpidx*everynth],AUndu.z0,WLofPixel,LambdaRfit,delta0fit+listOfboots[Strpidx*everynth]*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit,elevationfit,whichdata,cfit)


            #residualbeforeplot = sum(((np.array(y_fit-dataBooted[Strpidx, :])[0:1000])**2.)/(1000.*datasigmaBooted[Strpidx][0:1000]**2.))
            #print(residualbeforeplot)

            #innerresidual = sum(((np.array(y_fit-dataBooted[Strpidx, :])[400:600])**2.)/(200.*datasigmaBooted[Strpidx][400:600]**2.))
            #print(innerresidual)

            #maxofdata = np.around(max(dataBooted[Strpidx, :]), decimals=-1)
            #maxofSourcelstIntensity1D = np.around(max(SourcelstIntensity1D), decimals=-1)
            #maxofplt = max([maxofdata,maxofSourcelstIntensity1D])
            
            y_fitLst.append(y_fit)

            fntsze = 13

            residLst = np.array(residLst)
            residimg = residLst.reshape((NumOfSeries, 1000))[Strpidx*everynth, :]#nptp(residLst.reshape((NumOfSeries, 1000)))[Strpidx, :]
            print(residLst)
            

            #print(datasigmaBootedAnalyt)
            #print(datasigmaBootedAnalyt[Strpidx])
            figA = plt.figure("Diffraction")
            axA = figA.add_subplot(111)
            figB = plt.figure("Sources")
            axB = figB.add_subplot(111)
            #axA.annotate(r'$\Delta d + d = %d $' % k, xy=(0, maxofplt), fontsize=fntsze)
            axA.set_title("Data and Fit and resid")
            #axA.plot(1E6*yyy[0:1000], dataBooted[Strpidx, :], 'b-', 1000.*yyy, y_fit, 'r-', 1000.*yyy,residimg, 1000.*yyy,datasigmaBootedAnalyt[Strpidx], 'k-')
            axA.plot(1E6*yyy[0:1000]/AUndu.z0, dataBooted[Strpidx*everynth, :], 'b-', 1E6*yyy/AUndu.z0, y_fit, 'r-')
            #axA.set_xlim(-7, 7)
            #axA.set_ylim(-10, maxofplt+10)
            #axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
            axA.set_xlabel(r'$\Theta_y$ [$\mu$rad]',fontsize=fntsze)
            axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
            axA.tick_params(axis='both', which='major', labelsize=fntsze)
            axA.legend(['Data','Fit'])
            


            axB.set_title("Intensity without diffraction")
            #axB.annotate(r'$\Delta d + d = %d $' % k, xy=(0, maxofplt), fontsize=fntsze)
            axB.plot(1E6*yyy/AUndu.z0, SourcelstIntensity1D, '-')
            #axB.set_xlim(-7, 7)
            #axB.set_ylim(-10, maxofplt+10)
            #axB.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
            axB.set_xlabel(r'$\Theta_y$ [$\mu$rad]',fontsize=fntsze)
            axB.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
            axB.tick_params(axis='both', which='major', labelsize=fntsze)
            #axB.legend(['Fit'])
            plt.show()
            #"""

    #exit()

    #"""

    UnDi.SingleAmplitudeDone = 0

    if 1:
    
        residLst = np.array(residLst)
        residimg = nptp(residLst.reshape((NumOfSeries, 1000)))

        residLstsigma = np.array(residLst)[0]*np.array(datasigmaBooted.flatten())
        residimgsigma = nptp(residLstsigma.reshape((NumOfSeries, 1000)))

        #residimg8bit = 254*(np.array(residLst[0])-min(residLst[0]))/(max(residLst[0])-min(residLst[0]))
        #residimg8bit = nptp(residimg8bit.reshape((NumOfSeries, 1000))).astype(np.uint8)

        #cv2.imshow("resid",residimg)
        #cv2.waitKey(0)
        
        #originalfilename = 'residimg%d.tif'%k
        filename = 'residimg490.tif'
        filenamesigma = 'residimg%dsigma.tif'%k
        # Using cv2.imwrite() method 
        # Saving the image 
        cv2.imwrite(filename, residimg) 
        cv2.imwrite(filenamesigma, residimgsigma.astype(np.float32))
        

        #residimgcolor = cv2.applyColorMap(residimg8bit, cv2.COLORMAP_JET)
        #cv2.imshow("residcolor",residimgcolor)
        #cv2.waitKey(0)

        filenamecolor = 'residimgcolor%d.tif'%k
        # Using cv2.imwrite() method 
        # Saving the image 
        #cv2.imwrite(filenamecolor, residimgcolor) 

print(residsqrlst)
fntsze = 13
fig = plt.figure("Sum of squared resuduals vs. Position on camera screen")
ax = fig.add_subplot(111)
ax.set_title("Sum of squared residuals vs. Position on camera screen")
ax.set_xlabel(r"Position on camera screen [$\mu$m]",fontsize=fntsze)
ax.set_ylabel(r"Sum of squared residuals [a.u.]",fontsize=fntsze)
ax.plot(np.array(range(len(residsqrlst)))*10.+190.,residsqrlst)
plt.show()

exit()
#[487560.6, 487007.5, 468904.1, 443914.88, 451031.88, 438164.12, 431601.34, 431989.78, 449037.34, 403003.22, 406521.25, 451916.78, 470088.06, 357902.25, 345925.72, 329069.53, 350599.38, 332844.4, 319250.5, 313760.56, 325451.53, 331350.6, 329683.5, 337734.8, 334034.03, 318151.28, 303200.8, 297692.97, 305630.03, 312934.9, 311024.1, 330209.16, 345175.2, 344390.4, 349874.16, 399000.38, 429087.12, 371948.12, 356589.5, 366830.06, 383465.16, 405741.7, 392860.25, 407189.12]