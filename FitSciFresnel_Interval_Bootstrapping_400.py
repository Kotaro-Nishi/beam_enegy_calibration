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

ratio = 1000. / 606#606. #how many pixel are aperture
UnDi.apy = .004


zetaprime2 = np.linspace(-ratio*UnDi.apy,ratio*UnDi.apy,1000)
fixedsize = .245


NumOfStripes = 1
NumOfSeries = 50#825
LengthOfInterval = 400#825

dispersion = 0.0048085643228767475e-9
pixelof404nm = 1243.8551272671261
CalibWl = 404.6565e-9

#(407.7837-404.6565)

chosenpx = 1920

global listOfboots
listOfboots = []

Yscreenfit = .2e-4#-.5e-4
LambdaRfit = 400e-9

def per_iteration(pars, iteration, resid, *args, **kws):
    print(" resid ", sum(resid))#sum(abs(np.array(resid)))/1000.)
    #print(" resid ", resid)
    print("Amplfit: " + str(pars['Amplfit_%i'% listOfboots[0]].value))
    print("bfit: " + str(pars['bfit_%i'% listOfboots[0]].value))
    print("zpfit: " + str(pars['zpfit_%i'% listOfboots[0]].value))
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
    zpfit = pars['zpfit_%i'% listOfboots[0]].value
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
               "amplzpfitCorrel.csv",
               "YtrangammafitCorrel.csv",
               "Ytrandelta0fitCorrel.csv",
               "YtranbfitCorrel.csv",
               "gammadelta0fitCorrel.csv",
               "gammabfitCorrel.csv",
               "bfitzpfitCorrel.csv",
               "delta0fitzpfitCorrel.csv",
               "gammazpfitCorrel.csv",
               "YtranzpfitCorrel.csv",
               "zpfit.csv",
               "zpfitstderr.csv",
               "delta0fitbfitCorrel.csv"]

badk = []

DoNotSave = 0

Amplfit = -100
Ytran = -100
gammafit = -100
delta0fit = -100
bfit = -100

kstep = 1

WLofPixel = (chosenpx - pixelof404nm)*dispersion + CalibWl

k0 = 000
k = k0

SeriesListPictureMinusBGWidthOfStripe4Mean=[]
SeriesListPictureMinusBGWidthOfStripeSigma=[]

for ii in range(400):#int(len(filelist)/4)):#range(k0,k0+LengthOfInterval):
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

    ListPictureMinusBGWidthOfStripeSigma = np.std(ListPictureMinusBGWidthOfStripeForSigma, axis=0)

    ListPictureMinusBGWidthOfStripe4Mean = sum(ListPictureMinusBGWidthOfStripe)/4.
    
    SeriesListPictureMinusBGWidthOfStripe4Mean.append(ListPictureMinusBGWidthOfStripe4Mean)
    SeriesListPictureMinusBGWidthOfStripeSigma.append(ListPictureMinusBGWidthOfStripeSigma)

SaveSharedP = "\SharedP_20210401_100_000400_825_seed10000_new_zpfit_after3160"


checkstart = 1

random.seed(10000)
for i in range(NumOfSeries*3160):
    random.randint(0,LengthOfInterval-1)

for Saves in range(3160,10000):
    print("Saves " + str(Saves))

    zetaprime2Booted = []
    dataBooted = []
    datasigmaBooted = []
    listOfboots = []
    UnDi.SingleAmplitudeDone = 0
    for i in range(NumOfSeries):
        listOfboots.append(random.randint(0,LengthOfInterval-1))
    
    for SIt in range(len(listOfboots)):
        if listOfboots.count(listOfboots[SIt]) == 1:
            FirstSingularitem = SIt#listOfboots[SIt]
            break

    residLst = []
    if checkstart == 1:
        checklist = os.listdir(DataPath + SaveSharedP)

        maybeBadGamma = 0
        for ckid in range(len(targetfiles)):
            if targetfiles[ckid] in checklist:
                os.remove(DataPath + SaveSharedP + "//" + targetfiles[ckid])

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
    datasigma = np.array(interpolatedsigmapicture)

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

    Guessvalueslmfit = [1.4,.0001,350.9,.190+.01+k*.001,.6,4.3]#,3.e-4]

    fit_params = Parameters()

    for iy,y, in enumerate(dataBooted):
        fit_params.add('Amplfit_%i' % (listOfboots[iy]),value=Guessvalueslmfit[0])#, max=2.5, min=1.23)#, min = Guessvalueslmfit[0]-1.1, max=Guessvalueslmfit[0]+1.1)
        fit_params.add('Ytran_%i' % (listOfboots[iy]),value=Guessvalueslmfit[1])#, min = Guessvalueslmfit[1]-.0001, max=Guessvalueslmfit[1]+.0001)
        fit_params.add('gamma_%i' % (listOfboots[iy]),value=Guessvalueslmfit[2])
        fit_params.add('delta0_%i' % (listOfboots[iy]),value=Guessvalueslmfit[3])#, min = Guessvalueslmfit[4]-.1, max=Guessvalueslmfit[4]+.1)

        fit_params.add('bfit_%i' % (listOfboots[iy]),value=Guessvalueslmfit[4])#, min = Guessvalueslmfit[4]-.1, max=Guessvalueslmfit[4]+.1)
        fit_params.add('zpfit_%i' % (listOfboots[iy]),value=Guessvalueslmfit[5])#, min = Guessvalueslmfit[4]-.1, max=Guessvalueslmfit[4]+.1)

    listOfbootsnoFirst = [x for x in listOfboots if x != listOfboots[FirstSingularitem]]
    for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
        fit_params['Amplfit_%i' % listOfbootsnoFirst[iy]].expr='Amplfit_%i'% listOfboots[FirstSingularitem]
        
    for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
        fit_params['Ytran_%i' % listOfbootsnoFirst[iy]].expr='Ytran_%i'% listOfboots[FirstSingularitem]
        
    for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
        fit_params['gamma_%i' % listOfbootsnoFirst[iy]].expr='gamma_%i'% listOfboots[FirstSingularitem]
        
    for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
        fit_params['delta0_%i' % listOfbootsnoFirst[iy]].expr='delta0_%i'% listOfboots[FirstSingularitem]

    for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
        fit_params['bfit_%i' % listOfbootsnoFirst[iy]].expr='bfit_%i'% listOfboots[FirstSingularitem]

    for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
        fit_params['zpfit_%i' % listOfbootsnoFirst[iy]].expr='zpfit_%i'% listOfboots[FirstSingularitem]

    mini = Minimizer(objective_diffraction, fit_params,nan_policy = "raise",scale_covar=False,iter_cb=per_iteration, fcn_args=(zetaprime2,dataBooted,[1.]*len(dataBooted)),)#,1./datasigmaBooted),)
    
    print(len(zetaprime2))
    print(len(datasigmaBooted[0]))

    #try:
    result = mini.leastsq(ftol=1.e-2)
    #except RuntimeError:
    #    badSample+=1
    #    continue

    report_fit(result.params,show_correl=True, min_correl=0.01)
    
    residLst.append(result.residual)

    y_fitLst = []
    dataLst = []
    
    """

    for Strpidx in range(NumOfSeries):

        Sourcelst=[]

        Amplfit = result.params['Amplfit_1'].value
        Ytran = result.params['Ytran_1'].value
        gammafit = result.params['gamma_1'].value
        delta0fit = result.params['delta0_1'].value
        bfit = result.params['bfit_1'].value
        zpfit = result.params['zpfit_1'].value
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

    for Strpidx in range(NumOfSeries):
        if 1:#Strpidx > 600 and Strpidx > 610:

            #y_fit = diffraction_dataset(result.params, Strpidx, zetaprime2)
            #y_fit = UnDi.diffractionintensitySci(zetaprime2,*popt)#
            #y_fit = UnDi.diffractionintensitySci_delta(zetaprime2,deltad)#diffraction_dataset(result.params, Strpidx, zetaprime2)

            Sourcelst=[]
            
            Amplfit = result.params['Amplfit_1'].value
            Ytran = result.params['Ytran_1'].value
            gammafit = result.params['gamma_1'].value
            delta0fit = result.params['delta0_1'].value
            bfit = result.params['bfit_1'].value
            zpfit = result.params['zpfit_1'].value
            #Yscreenfit = result.params['Yscreenfit_1'].value
            #LambdaRfit = result.params['LambdaRfit_1'].value        

            #deltad = delta0+Strpidx*.001

            Sourcelst = UnDi.plainintensity(zetaprime2,Ytran,Strpidx,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+.001*Strpidx,gammafit,Amplfit,bfit,zpfit,Yscreenfit)

            SourcelstIntensity1D = 55000.*abs(np.array(Sourcelst))

            y_fit = UnDi.diffractionintensity(zetaprime2,Ytran,Strpidx,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+Strpidx*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit)
                   #UnDi.diffractionintensitySci(zetaprime2,Ytran,delta0fit+.001*Strpidx,gammafit,Amplfit,bfit)

            #residualbeforeplot = sum(((np.array(y_fit-dataBooted[Strpidx, :])[0:1000])**2.)/(1000.*datasigmaBooted[Strpidx][0:1000]**2.))
            #print(residualbeforeplot)

            #innerresidual = sum(((np.array(y_fit-dataBooted[Strpidx, :])[400:600])**2.)/(200.*datasigmaBooted[Strpidx][400:600]**2.))
            #print(innerresidual)

            #maxofdata = np.around(max(dataBooted[Strpidx, :]), decimals=-1)
            #maxofSourcelstIntensity1D = np.around(max(SourcelstIntensity1D), decimals=-1)
            #maxofplt = max([maxofdata,maxofSourcelstIntensity1D])
            
            y_fitLst.append(y_fit)

            fntsze = 13


            figA = plt.figure("Diffraction")
            axA = figA.add_subplot(111)
            figB = plt.figure("Sources")
            axB = figB.add_subplot(111)

            #axA.annotate(r'$\Delta d + d = %d $' % k, xy=(0, maxofplt), fontsize=fntsze)
            axA.plot(1000.*yyy[0:1000], dataBooted[Strpidx, :], '-', 1000.*yyy, y_fit, '-')
            #axA.set_xlim(-7, 7)
            #axA.set_ylim(-10, maxofplt+10)
            axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
            axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
            axA.tick_params(axis='both', which='major', labelsize=fntsze)
            
            #axB.annotate(r'$\Delta d + d = %d $' % k, xy=(0, maxofplt), fontsize=fntsze)
            axB.plot(1000.*yyy, SourcelstIntensity1D, '-')
            #axB.set_xlim(-7, 7)
            #axB.set_ylim(-10, maxofplt+10)
            axB.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
            axB.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
            axB.tick_params(axis='both', which='major', labelsize=fntsze)

            plt.show()

    #"""

    UnDi.SingleAmplitudeDone = 0

    Amplfit = result.params['Amplfit_%i'% listOfboots[FirstSingularitem]].value
    Ytran = result.params['Ytran_%i'% listOfboots[FirstSingularitem]].value
    gammafit = result.params['gamma_%i'% listOfboots[FirstSingularitem]].value
    delta0fit = result.params['delta0_%i'% listOfboots[FirstSingularitem]].value
    bfit = result.params['bfit_%i'% listOfboots[FirstSingularitem]].value
    zpfit = result.params['zpfit_%i'% listOfboots[FirstSingularitem]].value

    Amplfitstderr = result.params['Amplfit_%i'% listOfboots[FirstSingularitem]].stderr
    Ytranstderr = result.params['Ytran_%i'% listOfboots[FirstSingularitem]].stderr
    gammafitstderr = result.params['gamma_%i'% listOfboots[FirstSingularitem]].stderr
    delta0fitstderr = result.params['delta0_%i'% listOfboots[FirstSingularitem]].stderr
    bfitstderr = result.params['bfit_%i'% listOfboots[FirstSingularitem]].stderr
    zpfitstderr = result.params['zpfit_%i'% listOfboots[FirstSingularitem]].stderr

    amplYtranCorrel = result.covar[0][1] / (result.params["Amplfit_%i"% listOfboots[FirstSingularitem]].stderr*result.params["Ytran_%i"% listOfboots[FirstSingularitem]].stderr)
    amplgammafitCorrel = result.covar[0][2] / (result.params["Amplfit_%i"% listOfboots[FirstSingularitem]].stderr*result.params["gamma_%i"% listOfboots[FirstSingularitem]].stderr)
    ampldelta0fitCorrel = result.covar[0][3] / (result.params["Amplfit_%i"% listOfboots[FirstSingularitem]].stderr*result.params["delta0_%i"% listOfboots[FirstSingularitem]].stderr)
    amplbfitCorrel = result.covar[0][4] / (result.params["Amplfit_%i"% listOfboots[FirstSingularitem]].stderr*result.params["bfit_%i"% listOfboots[FirstSingularitem]].stderr)
    amplzpfitCorrel = result.covar[0][5] / (result.params["Amplfit_%i"% listOfboots[FirstSingularitem]].stderr*result.params["zpfit_%i"% listOfboots[FirstSingularitem]].stderr)
    YtrangammafitCorrel = result.covar[1][2] / (result.params["Ytran_%i"% listOfboots[FirstSingularitem]].stderr*result.params["gamma_%i"% listOfboots[FirstSingularitem]].stderr)
    Ytrandelta0fitCorrel = result.covar[1][3] / (result.params["Ytran_%i"% listOfboots[FirstSingularitem]].stderr*result.params["delta0_%i"% listOfboots[FirstSingularitem]].stderr)
    YtranbfitCorrel = result.covar[1][4] / (result.params["Ytran_%i"% listOfboots[FirstSingularitem]].stderr*result.params["bfit_%i"% listOfboots[FirstSingularitem]].stderr)
    YtranzpfitCorrel = result.covar[1][5] / (result.params["Ytran_%i"% listOfboots[FirstSingularitem]].stderr*result.params["zpfit_%i"% listOfboots[FirstSingularitem]].stderr)
    gammadelta0fitCorrel = result.covar[2][3] / (result.params["gamma_%i"% listOfboots[FirstSingularitem]].stderr*result.params["delta0_%i"% listOfboots[FirstSingularitem]].stderr)
    gammabfitCorrel = result.covar[2][4] / (result.params["gamma_%i"% listOfboots[FirstSingularitem]].stderr*result.params["bfit_%i"% listOfboots[FirstSingularitem]].stderr)
    gammazpfitCorrel = result.covar[2][5] / (result.params["gamma_%i"% listOfboots[FirstSingularitem]].stderr*result.params["zpfit_%i"% listOfboots[FirstSingularitem]].stderr)
    delta0fitbfitCorrel = result.covar[3][4] / (result.params["delta0_%i"% listOfboots[FirstSingularitem]].stderr*result.params["bfit_%i"% listOfboots[FirstSingularitem]].stderr)
    delta0fitzpfitCorrel = result.covar[3][5] / (result.params["delta0_%i"% listOfboots[FirstSingularitem]].stderr*result.params["zpfit_%i"% listOfboots[FirstSingularitem]].stderr)
    bfitzpfitCorrel = result.covar[4][5] / (result.params["bfit_%i"% listOfboots[FirstSingularitem]].stderr*result.params["zpfit_%i"% listOfboots[FirstSingularitem]].stderr)

    if DoNotSave == 0:
    
        #residLst = np.array(residLst*(datasigmaBooted.reshape((1, NumOfSeries*1000))))
        #residimg = nptp(residLst.reshape((NumOfSeries, 1000)))
    
        #with open(DataPath + "\\SharedP\\residimg.csv", 'a+', newline='') as write_obj:
        #    # Create a writer object from csv module
        #    csv_writer = csv.writer(write_obj)
        #    # Add contents of list as last row in the csv file
        #    csv_writer.writerows(residimg)
        #exit()
            
        with open(DataPath + SaveSharedP + "\\Amplfit.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([Amplfit])

        with open(DataPath + SaveSharedP + "\\Ytran.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([Ytran])

        with open(DataPath + SaveSharedP + "\\gammafit.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([gammafit])

        with open(DataPath + SaveSharedP + "\\delta0fit.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([delta0fit])

        with open(DataPath + SaveSharedP + "\\bfit.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([bfit])
            
        with open(DataPath + SaveSharedP + "\\zpfit.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([zpfit])


        with open(DataPath + SaveSharedP + "\\Amplfitstderr.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([Amplfitstderr])

        with open(DataPath + SaveSharedP + "\\Ytranstderr.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([Ytranstderr])

        with open(DataPath + SaveSharedP + "\\gammafitstderr.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([gammafitstderr])

        with open(DataPath + SaveSharedP + "\\delta0fitstderr.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([delta0fitstderr])

        with open(DataPath + SaveSharedP + "\\bfitstderr.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([bfitstderr])

        with open(DataPath + SaveSharedP + "\\zpfitstderr.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([zpfitstderr])


        with open(DataPath + SaveSharedP + "\\amplYtranCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([amplYtranCorrel])

        with open(DataPath + SaveSharedP + "\\amplgammafitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([amplgammafitCorrel])
        
        with open(DataPath + SaveSharedP + "\\ampldelta0fitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([ampldelta0fitCorrel])
            
        with open(DataPath + SaveSharedP + "\\amplbfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([amplbfitCorrel])
            
        with open(DataPath + SaveSharedP + "\\amplzpfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([amplzpfitCorrel])

        with open(DataPath + SaveSharedP + "\\YtrangammafitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([YtrangammafitCorrel])

        with open(DataPath + SaveSharedP + "\\Ytrandelta0fitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([Ytrandelta0fitCorrel])

        with open(DataPath + SaveSharedP + "\\YtranbfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([YtranbfitCorrel])

        with open(DataPath + SaveSharedP + "\\YtranzpfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([YtranzpfitCorrel])
        
        with open(DataPath + SaveSharedP + "\\gammadelta0fitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([gammadelta0fitCorrel])
            
        with open(DataPath + SaveSharedP + "\\gammabfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([gammabfitCorrel])

        with open(DataPath + SaveSharedP + "\\gammazpfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([gammazpfitCorrel])

        with open(DataPath + SaveSharedP + "\\delta0fitbfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([delta0fitbfitCorrel])
            
        with open(DataPath + SaveSharedP + "\\delta0fitzpfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([delta0fitzpfitCorrel])

        with open(DataPath + SaveSharedP + "\\bfitzpfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([bfitzpfitCorrel])

    
    checkstart = 0
exit()