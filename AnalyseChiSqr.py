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
NumOfSeries = 100#825
LengthOfInterval = 425#825

dispersion = 0.004808759252942803e-9 #Former 0.0048085643228767475e-9
pixelof404nm = 1315.3024745554976 #Former: 1243.8551272671261
CalibWl = 404.6565e-9

#(407.7837-404.6565)

chosenpx = 1300

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
# [1.4,.0001,352.4,.190+.01+k*.001,.6,4.3]
def diffraction_datasetChiSqrAnalysis(params,i,yy):
    #Amplfit = params['Amplfit_%i' % (i)].value
    #Ytran   = params['Ytran_%i' % (i)].value
    #gammafit  = params['gamma_%i' % (i)].value
    #delta0fit  = params['delta0_%i' % (i)].value
    #bfit  = params['bfit_%i' % (i)].value
    #zpfit  = params['zpfit_%i' % (i)].value

    #Yscreenfit  = params['Yscreenfit_%i' % (i+1)].value
    #LambdaRfit  = params['LambdaRfit_%i' % (i+1)].value
    #print(delta0fit+i*.001)

    Diff_intensity = UnDi.diffractionintensity(yy,params[1],i,AUndu.z0,WLofPixel,LambdaRfit,params[3]+i*.001,params[2],params[0],params[4],params[5],Yscreenfit)
                          #diffractionintensity(yy,Ytran,i,AUndu.z0,WLofPixel,LambdaRfit,delta0fit+i*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit)


    return Diff_intensity

def objective_diffraction(params, yy, data,weights):
    ndata,nx = data.shape
    
    resid = 0.0*data[:]
    
    for Strpidx in range(ndata):
        if Strpidx > 0:
            UnDi.SingleAmplitudeDone = 1
        #resid[Strpidx] = (diffraction_dataset(params,listOfboots[Strpidx],yy) - data[Strpidx])*weights[Strpidx]#*(1.+50000.*yy**2.)#**2.
        resid[Strpidx] = (diffraction_datasetChiSqrAnalysis(params,listOfboots[Strpidx],yy) - data[Strpidx])*weights[Strpidx]#*(1.+50000.*yy**2.)#**2.
        
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

FolderList = [
"\SharedP_20210401_100_000425_825_seed10000_new_zpfit",
"\SharedP_20210401_100_025450_825_seed10000_new_zpfit",
"\SharedP_20210401_100_050475_825_seed10000_new_zpfit",
"\SharedP_20210401_100_075500_825_seed10000_new_zpfit",
"\SharedP_20210401_100_100525_825_seed10000_new_zpfit",
"\SharedP_20210401_100_125550_825_seed10000_new_zpfit",
"\SharedP_20210401_100_150575_825_seed10000_new_zpfit",
"\SharedP_20210401_100_175600_825_seed10000_new_zpfit",
"\SharedP_20210401_100_200625_825_seed10000_new_zpfit",
"\SharedP_20210401_100_225650_825_seed10000_new_zpfit",
"\SharedP_20210401_100_250675_825_seed10000_new_zpfit",
"\SharedP_20210401_100_275700_825_seed10000_new_zpfit",
"\SharedP_20210401_100_300725_825_seed10000_new_zpfit",
"\SharedP_20210401_100_325750_825_seed10000_new_zpfit",
"\SharedP_20210401_100_350775_825_seed10000_new_zpfit",
"\SharedP_20210401_100_375800_825_seed10000_new_zpfit",
"\SharedP_20210401_100_400825_825_seed10000_new_zpfit",
]

random.seed(10000)

for i_folder in range(0,1):#len(FolderList)):
    SeriesListPictureMinusBGWidthOfStripe4Mean=[]
    SeriesListPictureMinusBGWidthOfStripeSigma=[]
    for ii in range(i_folder*25,i_folder*25+LengthOfInterval):#int(len(filelist)/4)):#range(k0,k0+LengthOfInterval):
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



    SaveSharedP = FolderList[i_folder]#"\SharedP_20210401_100_300725_825_seed10000_new_zpfit"


    checkstart = 1

    for Saves in range(0,500):
        print("Saves " + str(Saves))

        zetaprime2Booted = []
        dataBooted = []
        datasigmaBooted = []
        listOfboots = []
        UnDi.SingleAmplitudeDone = 0
        #for i in range(NumOfSeries):
        #    listOfboots.append(i)
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
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
        datasigma = np.array(interpolatedsigmapicture)

        for Boot in range(NumOfSeries):
            dataBooted.append(data[listOfboots[Boot]-i_folder*25])
            datasigmaBooted.append(datasigma[listOfboots[Boot]-i_folder*25])
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
        """
        gamma_267:    352.445079 +/- 0.04167944 (0.01%) == 'gamma_295'
        delta0_267:   0.20615791 +/- 1.0607e-04 (0.05%) == 'delta0_295'
        bfit_267:     0.45970915 +/- 0.01294143 (2.82%) == 'bfit_295'
        zpfit_267:    5.38063293 +/- 0.00709303 (0.13%) == 'zpfit_295'
        Amplfit_424:  0.95072544 +/- 0.00187120 (0.20%) == 'Amplfit_295'
        Ytran_424:    1.4818e-04 +/- 1.0456e-06 (0.71%) == 'Ytran_295'

        Amplfit: 0.9507254433277184
        bfit: 0.45970915123822237
        zpfit: 5.380632929941211
        Ytran: 0.00014818383372392456
        gammafit: 352.445079279404
        delta0fit: 0.206157906592224

        Amplfit: 1.4
        bfit: 0.6
        zpfit: 4.3
        Ytran: 0.0001
        gammafit: 352.4
        delta0fit: 0.2
        
        141251.7

        """

        Guessvalueslmfit = [0.95072544,1.4818e-04,352.445079,0.20615791,0.45970915,5.3806329]
        #Guessvalueslmfit = [1.4,.0001,352.4,0.2,0.6,4.3]



        #print(objective_diffraction(Guessvalueslmfit, zetaprime2, dataBooted,datasigmaBooted))
        #print(len(objective_diffraction(Guessvalueslmfit, zetaprime2, dataBooted,datasigmaBooted)))
        print(sum(objective_diffraction(Guessvalueslmfit, zetaprime2, dataBooted,1/datasigmaBooted)))
        
        gammaChiLst = np.linspace(352.3,352.6,100)
        zpfitChiLst = np.linspace(4.95,5.45,100)
        
        fntsze = 13
        #figA = plt.figure("Diffraction")
        #axA = figA.add_subplot(111)
        
        
        AmpCent = 0.9507254433277184
        YtranCent = 0.00014818383372392456
        delta0Cent = 0.206157906592224
        bCent = 0.45970915123822237
        
        listOfAmp = []
        listOfYtran = []
        listOfdelta0 = []
        listOfb = []

        for i in range(len(zpfitChiLst)*len(gammaChiLst)):
            listOfAmp.append(random.uniform(AmpCent*.998, AmpCent/.998))
            listOfYtran.append(random.uniform(YtranCent*.9929, YtranCent/.9929))
            listOfdelta0.append(random.uniform(delta0Cent*.9995, delta0Cent/.9995))
            listOfb.append(random.uniform(bCent*.9718, bCent/.9718))

        ChiSQ_Lst2D = []
        for ii_Chi in range(len(zpfitChiLst)):
            ChiSQ_Lst = []
            for i_Chi in range(len(gammaChiLst)):
                #Guessvalueslmfit = [0.9507254433277184,0.00014818383372392456,gammaChiLst[i_Chi],0.206157906592224,0.45970915123822237,5.380632929941211]
                Guessvalueslmfit = [listOfAmp[i_Chi+len(gammaChiLst)*ii_Chi],listOfYtran[i_Chi+len(gammaChiLst)*ii_Chi],gammaChiLst[i_Chi],listOfdelta0[i_Chi+len(gammaChiLst)*ii_Chi],listOfb[i_Chi+len(gammaChiLst)*ii_Chi],zpfitChiLst[ii_Chi]]
                ChiSQ_Lst.append(sum((objective_diffraction(Guessvalueslmfit, zetaprime2, dataBooted,1/datasigmaBooted))**2.))
                print(i_Chi)
            ChiSQ_Lst2D.append(ChiSQ_Lst)

        plt.imshow(np.array(ChiSQ_Lst2D), cmap='jet', interpolation='nearest')
        plt.show()
        #axA.plot(gammaChiLst, ChiSQ_Lst, 'b.')
        #axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
        #axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
        #axA.tick_params(axis='both', which='major', labelsize=fntsze)
        
        #plt.show()
        
        exit()

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