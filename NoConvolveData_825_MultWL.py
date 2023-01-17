import matplotlib.pyplot as plt
import UndulatorDiffractionOneDimensional_correct_XAmpl_InitialFit_NoConvolve as UnDi
import InterferenceFunction_more_x as IFunc
import AmplitudeFormula_x as AUndu

from itertools import chain

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

NumOfStripes = 65      #How many Wavelengths are scanned
NumOfSeries = 825       #How many of the datapoints of the interval are used for the fit.
LengthOfInterval = 825 #How many datapoints are included for each fit.

dispersion = 0.004808759252942803e-9 #Former 0.0048085643228767475e-9
pixelof404nm = 1315.3024745554976 #Former: 1243.8551272671261
CalibWl = 404.6565e-9

#(407.7837-404.6565)

global listOfboots
#global pltweights
listOfboots = []
#pltweights = []

#Yscreenfit = (380)*10**-6
#Yscreenfit = (460)*10**-6
Yscreenfit = (300)*10**-6
LambdaRfit = 412e-9

def per_iteration(pars, iteration, resid, *args, **kws):
    
    print(" resid ", sum(resid))
    print(" residSqr ", sum(resid**2.))

    print("Amplfit: " + str(pars['Amplfit_%i'% listOfboots[0]].value))
    print("zpfit: " + str(pars['zpfit_%i'% listOfboots[0]].value))
    print('Ytran: ' + str(pars['Ytran_%i'% listOfboots[0]].value))
    print("gammafit: " + str(pars['gamma_%i'% listOfboots[0]].value))
    print("delta0fit: " + str(pars['delta0_%i'% listOfboots[0]].value))
    print("elevationfit: " + str(pars['elevationfit_%i'% listOfboots[0]].value))
    print("cfit: " + str(pars['cfit_%i'% listOfboots[0]].value))
    #print("bfit: " + str(pars['bfit_%i'% listOfboots[0]].value))

    Amplfit = pars['Amplfit_%i'% listOfboots[0]].value
    Ytran = pars['Ytran_%i'% listOfboots[0]].value
    gammafit = pars['gamma_%i'% listOfboots[0]].value
    delta0fit = pars['delta0_%i'% listOfboots[0]].value
    zpfit = pars['zpfit_%i'% listOfboots[0]].value
    elevationfit = pars['elevationfit_%i'% listOfboots[0]].value
    cfit = pars['cfit_%i'% listOfboots[0]].value
    #bfit = pars['bfit_%i'% listOfboots[0]].value

    UnDi.SingleAmplitudeDone = 0
    UnDi.AmplXY = []
    

def diffraction_dataset(params,i,yy):
    Amplfit = params['Amplfit_%i' % (i)].value
    Ytran   = params['Ytran_%i' % (i)].value
    gammafit  = params['gamma_%i' % (i)].value
    delta0fit  = params['delta0_%i' % (i)].value
    zpfit  = params['zpfit_%i' % (i)].value
    elevationfit  = params['elevationfit_%i' % (i)].value
    cfit = params['cfit_%i' % (i)].value
    #bfit  = params['bfit_%i' % (i)].value
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
        
        #1/5 average over five pixel columns
        #1/4 average over four images
        
        #resid[Strpidx] = (diffractionresult - data[Strpidx])*(1./(np.array([1.5]*len(weights[0]))+(1./.2)*(1./20.)*np.sqrt(0.05+abs(diffractionresult)*.2/.7)))
        #resid[Strpidx] = (diffractionresult - data[Strpidx])*(1./(np.array([1.5/1.]*len(data[0]))+(1./.2)*(1/5.)*(1/4.)*np.sqrt(0.05+abs(diffractionresult)*.2/.7)))
        resid[Strpidx] = (diffractionresult - data[Strpidx])*(1./(np.array([1.5]*240+[1.5]*606+[1.5]*(1000-606-240))+np.sqrt(1/5.)*np.sqrt(1/4.)*(1./.2)*np.sqrt(0.05+abs(diffractionresult)*.2/.7)))
        #pltweights.append((np.array([1.5/1]*len(weights[0]))+(1./.2)*(1./1.)*np.sqrt(0.05+abs(diffractionresult)*.2/.7)))
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

targetfiles = [#"residimg.csv",
               "Amplfit.csv",
               "Ytran.csv",
               "gammafit.csv",
               "delta0fit.csv",
               "cfit.csv",
               "zpfit.csv",
               "elevationfit.csv",
               #"bfit.csv",

               "Amplfitstderr.csv",
               "Ytranstderr.csv",
               "gammafitstderr.csv",
               "delta0fitstderr.csv",
               "cfitstderr.csv",
               "zpfitstderr.csv",
               "elevationfitstderr.csv",
               #"bfitstderr.csv",

               "AmplfitYtranCorrel.csv",
               "AmplfitgammafitCorrel.csv",
               "Amplfitdelta0fitCorrel.csv",
               "AmplfitcfitCorrel.csv",
               "AmplfitzpfitCorrel.csv",
               "AmplfitelevationfitCorrel.csv",
               #"AmplfitbfitCorrel.csv",

               "YtrangammafitCorrel.csv",
               "Ytrandelta0fitCorrel.csv",
               "YtrancfitCorrel.csv",
               "YtranzpfitCorrel.csv",
               "YtranelevationfitCorrel.csv",
               #"YtranbfitCorrel.csv",
               
               "gammadelta0fitCorrel.csv",
               "gammacfitCorrel.csv",
               "gammazpfitCorrel.csv",
               "gammaelevationfitCorrel.csv",
               #"gammabfitCorrel.csv",

               "delta0fitcfitCorrel.csv"
               "delta0fitzpfitCorrel.csv"
               "delta0fitelevationfitCorrel.csv"
               #"delta0bfitCorrel.csv"
               
               "cfitzpfitCorrel.csv"
               "cfitelevationfitCorrel.csv"
               #"cfitbfitCorrel.csv"
               
               "zpfitelevationfitCorrel.csv"
               #"zpfitbfitCorrel.csv"
               ]

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

startpx = 0

SeriesListPictureMinusBGWidthOfStripe4Mean=[]

for ii in range(LengthOfInterval):


    #print(chosenpx)
    #print((chosenpx-int(NumOfStripes*WidthOfStripe/2))-(chosenpx+int(NumOfStripes*WidthOfStripe/2)))
    #print(chosenpx-int(NumOfStripes*WidthOfStripe/2))
    #print(chosenpx+int(NumOfStripes*WidthOfStripe/2))

    ListPictureMinusBGWidthOfStripe = []

    print(filelist[ii*4])
    for i in range(4):

        PictureMinusBG = cv2.imread(DataPath + "/" + filelist[i+ii*4],cv2.IMREAD_UNCHANGED).astype(int)-Backgroundimage

        PictureMinusBG = np.float32(PictureMinusBG)

        PictureMinusBGPatchwork = []
        for s in range(NumOfStripes):
            chosenpx = 3 + startpx + s*35
            #PictureMinusBGPatchwork.append(nptp(nptp(PictureMinusBG)[(chosenpx-int(WidthOfStripe/2)):chosenpx+int(WidthOfStripe/2)+1]))
            PictureMinusBGPatchwork.append(nptp(nptp(PictureMinusBG)[(chosenpx-int(WidthOfStripe/2)):chosenpx+int(WidthOfStripe/2)+1]))

        #merged_list=np.array([item for sublist in nptp(PictureMinusBGPatchwork) for item in sublist])
        merged_list=nptp(np.array([item for sublist in nptp(PictureMinusBGPatchwork) for item in sublist]))

        merged_list = nptp(np.reshape(merged_list,(WidthOfStripe*NumOfStripes,1000)))

        merged_list = np.reshape(merged_list,(1000,WidthOfStripe*NumOfStripes))        

        """
        cv2.imshow("fdsa",np.array(merged_list))
        cv2.waitKey(0)
        cv2.destroyAllWindows()
        
        print(len(merged_list))
        print(len(merged_list[0]))

        #print(len(PictureMinusBGPatchwork[0][0][0]))
        #exit()
        #"""
        
        
        ListPictureMinusBGWidthOfStripe.append(merged_list)

    ListPictureMinusBGWidthOfStripe = np.array(ListPictureMinusBGWidthOfStripe)

    ListPictureMinusBGWidthOfStripe4Mean = sum(ListPictureMinusBGWidthOfStripe)/4.

    """
    if 1:#ii > 20:
        cv2.imshow("fdsa",ListPictureMinusBGWidthOfStripe4Mean)
        cv2.waitKey(0)
        cv2.destroyAllWindows()
        fntsze = 13
        figA = plt.figure("Diffraction")
        axA = figA.add_subplot(111)
        #axA.plot(range(len(np.transpose(PictureMinusBG)[100])), np.transpose(PictureMinusBG)[100], 'b.')
        #axA.plot(range(len(np.transpose(ListPictureMinusBGWidthOfStripe4Mean)[100])), np.transpose(ListPictureMinusBGWidthOfStripe4Mean)[100], 'r.')
        axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
        axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
        axA.tick_params(axis='both', which='major', labelsize=fntsze)
        plt.show()
    #"""

    
    SeriesListPictureMinusBGWidthOfStripe4Mean.append(ListPictureMinusBGWidthOfStripe4Mean)

SaveSharedP = "\SharedP_20210401_825_seed10000_MultiWL_NoConv"


checkstart = 1

#random.seed(10000)
#for i in range(NumOfSeries*6500):
#    random.randint(0,LengthOfInterval-1)
random.seed(10000)
WL_lst=[]

for Saves in range(0,NumOfStripes):
    print("Saves " + str(Saves))

    chosenpx = 3 + startpx + Saves*35
    WLofPixel = (chosenpx - pixelof404nm)*dispersion + CalibWl

    WL_lst.append(WLofPixel)
    print(chosenpx)
    print(WL_lst)

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
    #"""
    if checkstart == 1:
        checklist = os.listdir(DataPath + SaveSharedP)

        maybeBadGamma = 0
        for ckid in range(len(targetfiles)):
            if targetfiles[ckid] in checklist:
                os.remove(DataPath + SaveSharedP + "//" + targetfiles[ckid])
    #"""

    interpolatedpicture = []
    interpolatedsigmapicture = []

    for ii in range(LengthOfInterval):
        stripe = SeriesListPictureMinusBGWidthOfStripe4Mean[ii]
        
        #print(len(stripe))
        #print(len(stripe[0]))

        stripe = cv2.resize(stripe,(NumOfStripes,len(ListPictureMinusBGWidthOfStripe4Mean)),interpolation=cv2.INTER_AREA) #1 zu 1000

        #print(len(stripe))
        #print(len(stripe[0]))
        #print(stripe[0])
        #exit()

        ElectronSource = np.array(IFunc.ElectronSourcefunctionSet(zetaprime2,.000004))  #.00004

        #stripe = IFunc.convolveplain(nptp(stripe)[0],ElectronSource).real
        stripe = nptp(stripe)[Saves]

        stripejjj=[]
        for jjj in range(len(stripe)):
            stripejjj.append(stripe[jjj])
            
        stripe = nptp(np.array([stripejjj])*1.) #3000000.
        interpolatedpicture.append(nptp(stripe)[0])
        
    data = np.array(interpolatedpicture)

    for Boot in range(NumOfSeries):
        dataBooted.append(data[listOfboots[Boot]])
    dataBooted = np.array(dataBooted)

    print(type(dataBooted))
    print(type(dataBooted[0]))

    
    dataLst = []
    for Strpidx in range(NumOfSeries):
        dataLst.append(dataBooted[Strpidx][500])

    if Saves == 0:
        #Guessvalueslmfit = [1.4,.0001,352.4,.190+.01+k*.001,.6,5.2,13.]
        Guessvalueslmfit = [8.4,.0001,352.4,0.1863,3.6,5.2,13.,10.06]
    else:
        Guessvalueslmfit = [Amplfit,.0001,352.4,delta0fit,cfit,zpfit,elevationfit]

    fit_params = Parameters()

    for iy,y, in enumerate(dataBooted):
        fit_params.add('Amplfit_%i' % (listOfboots[iy]),value=Guessvalueslmfit[0])
        fit_params.add('Ytran_%i' % (listOfboots[iy]),value=Guessvalueslmfit[1])
        fit_params.add('gamma_%i' % (listOfboots[iy]),value=Guessvalueslmfit[2])
        fit_params.add('delta0_%i' % (listOfboots[iy]),value=Guessvalueslmfit[3])
        fit_params.add('cfit_%i' % (listOfboots[iy]),value=Guessvalueslmfit[4])
        fit_params.add('zpfit_%i' % (listOfboots[iy]),value=Guessvalueslmfit[5])
        fit_params.add('elevationfit_%i' % (listOfboots[iy]),value=Guessvalueslmfit[6])

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
        fit_params['cfit_%i' % listOfbootsnoFirst[iy]].expr='cfit_%i'% listOfboots[FirstSingularitem]

    for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
        fit_params['zpfit_%i' % listOfbootsnoFirst[iy]].expr='zpfit_%i'% listOfboots[FirstSingularitem]

    for iy in tuple(range(len(listOfbootsnoFirst))):#tuple(range(NumOfSeries)[1:]):
        fit_params['elevationfit_%i' % listOfbootsnoFirst[iy]].expr='elevationfit_%i'% listOfboots[FirstSingularitem]
        

    mini = Minimizer(objective_diffraction, fit_params,nan_policy = "raise",scale_covar=False,iter_cb=per_iteration, fcn_args=(zetaprime2,dataBooted,1./1.),)
    
    print(len(zetaprime2))


    #try:
    result = mini.leastsq(ftol=1.e-2)
    #except RuntimeError:
    #    badSample+=1
    #    continue

    report_fit(result.params,show_correl=True, min_correl=0.01)
    
    residLst.append(result.residual)

    y_fitLst = []
    dataLst = []

    UnDi.SingleAmplitudeDone = 0

    Amplfit = result.params['Amplfit_%i'% listOfboots[FirstSingularitem]].value
    Ytran = result.params['Ytran_%i'% listOfboots[FirstSingularitem]].value
    gammafit = result.params['gamma_%i'% listOfboots[FirstSingularitem]].value
    delta0fit = result.params['delta0_%i'% listOfboots[FirstSingularitem]].value
    cfit = result.params['cfit_%i'% listOfboots[FirstSingularitem]].value
    zpfit = result.params['zpfit_%i'% listOfboots[FirstSingularitem]].value
    elevationfit = result.params['elevationfit_%i'% listOfboots[FirstSingularitem]].value    
    #bfit = result.params['bfit_%i'% listOfboots[FirstSingularitem]].value

    Amplfitstderr = result.params['Amplfit_%i'% listOfboots[FirstSingularitem]].stderr
    Ytranstderr = result.params['Ytran_%i'% listOfboots[FirstSingularitem]].stderr
    gammafitstderr = result.params['gamma_%i'% listOfboots[FirstSingularitem]].stderr
    delta0fitstderr = result.params['delta0_%i'% listOfboots[FirstSingularitem]].stderr
    cfitstderr = result.params['cfit_%i'% listOfboots[FirstSingularitem]].stderr
    zpfitstderr = result.params['zpfit_%i'% listOfboots[FirstSingularitem]].stderr
    elevationfitstderr = result.params['elevationfit_%i'% listOfboots[FirstSingularitem]].stderr
    #bfitstderr = result.params['bfit_%i'% listOfboots[FirstSingularitem]].stderr

    AmplfitYtranCorrel = result.covar[0][1] / (result.params["Amplfit_%i"% listOfboots[FirstSingularitem]].stderr*result.params["Ytran_%i"% listOfboots[FirstSingularitem]].stderr)
    AmplfitgammafitCorrel = result.covar[0][2] / (result.params["Amplfit_%i"% listOfboots[FirstSingularitem]].stderr*result.params["gamma_%i"% listOfboots[FirstSingularitem]].stderr)
    Amplfitdelta0fitCorrel = result.covar[0][3] / (result.params["Amplfit_%i"% listOfboots[FirstSingularitem]].stderr*result.params["delta0_%i"% listOfboots[FirstSingularitem]].stderr)
    AmplfitcfitCorrel = result.covar[0][4] / (result.params["Amplfit_%i"% listOfboots[FirstSingularitem]].stderr*result.params["cfit_%i"% listOfboots[FirstSingularitem]].stderr)
    AmplfitzpfitCorrel = result.covar[0][5] / (result.params["Amplfit_%i"% listOfboots[FirstSingularitem]].stderr*result.params["zpfit_%i"% listOfboots[FirstSingularitem]].stderr)
    AmplfitelevationfitCorrel = result.covar[0][6] / (result.params["Amplfit_%i"% listOfboots[FirstSingularitem]].stderr*result.params["elevationfit_%i"% listOfboots[FirstSingularitem]].stderr)
    #AmplfitbfitCorrel = result.covar[0][4] / (result.params["Amplfit_%i"% listOfboots[FirstSingularitem]].stderr*result.params["bfit_%i"% listOfboots[FirstSingularitem]].stderr)
    
    YtrangammafitCorrel = result.covar[1][2] / (result.params["Ytran_%i"% listOfboots[FirstSingularitem]].stderr*result.params["gamma_%i"% listOfboots[FirstSingularitem]].stderr)
    Ytrandelta0fitCorrel = result.covar[1][3] / (result.params["Ytran_%i"% listOfboots[FirstSingularitem]].stderr*result.params["delta0_%i"% listOfboots[FirstSingularitem]].stderr)
    YtrancfitCorrel = result.covar[1][4] / (result.params["Ytran_%i"% listOfboots[FirstSingularitem]].stderr*result.params["cfit_%i"% listOfboots[FirstSingularitem]].stderr)
    YtranzpfitCorrel = result.covar[1][5] / (result.params["Ytran_%i"% listOfboots[FirstSingularitem]].stderr*result.params["zpfit_%i"% listOfboots[FirstSingularitem]].stderr)
    YtranelevationfitCorrel = result.covar[1][6] / (result.params["Ytran_%i"% listOfboots[FirstSingularitem]].stderr*result.params["elevationfit_%i"% listOfboots[FirstSingularitem]].stderr)

    gammadelta0fitCorrel = result.covar[2][3] / (result.params["gamma_%i"% listOfboots[FirstSingularitem]].stderr*result.params["delta0_%i"% listOfboots[FirstSingularitem]].stderr)
    gammacfitCorrel = result.covar[2][4] / (result.params["gamma_%i"% listOfboots[FirstSingularitem]].stderr*result.params["cfit_%i"% listOfboots[FirstSingularitem]].stderr)
    gammazpfitCorrel = result.covar[2][5] / (result.params["gamma_%i"% listOfboots[FirstSingularitem]].stderr*result.params["zpfit_%i"% listOfboots[FirstSingularitem]].stderr)
    gammaelevationfitCorrel = result.covar[2][6] / (result.params["gamma_%i"% listOfboots[FirstSingularitem]].stderr*result.params["elevationfit_%i"% listOfboots[FirstSingularitem]].stderr)

    delta0fitcfitCorrel = result.covar[3][4] / (result.params["delta0_%i"% listOfboots[FirstSingularitem]].stderr*result.params["cfit_%i"% listOfboots[FirstSingularitem]].stderr)
    delta0fitzpfitCorrel = result.covar[3][5] / (result.params["delta0_%i"% listOfboots[FirstSingularitem]].stderr*result.params["zpfit_%i"% listOfboots[FirstSingularitem]].stderr)
    delta0fitelevationfitCorrel = result.covar[3][6] / (result.params["delta0_%i"% listOfboots[FirstSingularitem]].stderr*result.params["elevationfit_%i"% listOfboots[FirstSingularitem]].stderr)
    
    cfitzpfitCorrel = result.covar[4][5] / (result.params["cfit_%i"% listOfboots[FirstSingularitem]].stderr*result.params["zpfit_%i"% listOfboots[FirstSingularitem]].stderr)
    cfitelevationfitCorrel = result.covar[4][6] / (result.params["cfit_%i"% listOfboots[FirstSingularitem]].stderr*result.params["elevationfit_%i"% listOfboots[FirstSingularitem]].stderr)
    
    zpfitelevationfitCorrel = result.covar[5][6] / (result.params["zpfit_%i"% listOfboots[FirstSingularitem]].stderr*result.params["elevationfit_%i"% listOfboots[FirstSingularitem]].stderr)

    if DoNotSave == 0:

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

        with open(DataPath + SaveSharedP + "\\cfit.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([cfit])

        with open(DataPath + SaveSharedP + "\\zpfit.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([zpfit])

        with open(DataPath + SaveSharedP + "\\elevationfit.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([elevationfit])



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

        with open(DataPath + SaveSharedP + "\\cfitstderr.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([cfitstderr])

        with open(DataPath + SaveSharedP + "\\zpfitstderr.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([zpfitstderr])

        with open(DataPath + SaveSharedP + "\\elevationfitstderr.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([elevationfitstderr])








        with open(DataPath + SaveSharedP + "\\AmplfitYtranCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([AmplfitYtranCorrel])

        with open(DataPath + SaveSharedP + "\\AmplfitgammafitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([AmplfitgammafitCorrel])
        
        with open(DataPath + SaveSharedP + "\\Amplfitdelta0fitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([Amplfitdelta0fitCorrel])

        with open(DataPath + SaveSharedP + "\\AmplfitcfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([AmplfitcfitCorrel])

        with open(DataPath + SaveSharedP + "\\AmplfitzpfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([AmplfitzpfitCorrel])
            
        with open(DataPath + SaveSharedP + "\\AmplfitelevationfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([AmplfitelevationfitCorrel])
            
            
            
            
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

        with open(DataPath + SaveSharedP + "\\YtrancfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([YtrancfitCorrel])

        with open(DataPath + SaveSharedP + "\\YtranzpfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([YtranzpfitCorrel])

        with open(DataPath + SaveSharedP + "\\YtranelevationfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([YtranelevationfitCorrel])
            
            
            
            
        with open(DataPath + SaveSharedP + "\\gammadelta0fitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([gammadelta0fitCorrel])

        with open(DataPath + SaveSharedP + "\\gammacfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([gammacfitCorrel])

        with open(DataPath + SaveSharedP + "\\gammazpfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([gammazpfitCorrel])

        with open(DataPath + SaveSharedP + "\\gammaelevationfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([gammaelevationfitCorrel])


            
        with open(DataPath + SaveSharedP + "\\delta0fitcfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([delta0fitcfitCorrel])

        with open(DataPath + SaveSharedP + "\\delta0fitzpfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([delta0fitzpfitCorrel])
            
        with open(DataPath + SaveSharedP + "\\delta0fitelevationfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([delta0fitelevationfitCorrel])



        with open(DataPath + SaveSharedP + "\\cfitzpfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([cfitzpfitCorrel])

        with open(DataPath + SaveSharedP + "\\cfitelevationfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([cfitelevationfitCorrel])



        with open(DataPath + SaveSharedP + "\\zpfitelevationfitCorrel.csv", 'a+', newline='') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow([zpfitelevationfitCorrel])

    checkstart = 0
exit()