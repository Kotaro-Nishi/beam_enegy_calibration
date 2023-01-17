
for nicestfit in range(1):
    y_fitLst = []
    dataLst = []
    everynth = 1
    whichdata = 0
    for Strpidx in range(int(NumOfSeries/everynth)):

        Amplfit = result.params['Amplfit_%i'% listOfboots[FirstSingularitem]].value
        Ytran = result.params['Ytran_%i'% listOfboots[FirstSingularitem]].value
        gammafit = result.params['gamma_%i'% listOfboots[FirstSingularitem]].value
        delta0fit = result.params['delta0_%i'% listOfboots[FirstSingularitem]].value
        bfit = 0
        zpfit = result.params['zpfit_%i'% listOfboots[FirstSingularitem]].value
        elevationfit = result.params['elevationfit_%i'% listOfboots[FirstSingularitem]].value
        cfit = result.params['cfit_%i'% listOfboots[FirstSingularitem]].value
        #Yscreenfit = result.params['Yscreenfit_%i'% listOfboots[FirstSingularitem]].value

        #y_fitLst.append(UnDi.diffractionintensityInitialFit_X(zetaprime2,zetaxx,Ytran,listOfboots[Strpidx*everynth],AUndu.z0,WLofPixel,LambdaRfit,delta0fit+listOfboots[Strpidx*everynth]*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit,elevationfit,whichdata,cfit)[500+nicestfit*500])
        y_fitLst.append(UnDi.diffractionintensityInitialFit_X(zetaprime2,zetaxx,Ytran,listOfboots[Strpidx*everynth],AUndu.z0,WLofPixel,LambdaRfit,delta0fit+posdata[listOfboots[Strpidx]]*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit,elevationfit,whichdata,cfit)[500])
        dataLst.append((dataBooted[Strpidx*everynth])[500+nicestfit*500])
    
    statuncert = np.array([3.5]*len(y_fitLst))+(np.sqrt(1/5.)*np.sqrt(1/4.)*(.7/.2)*np.sqrt(0.05+abs(np.array(y_fitLst))*.2/.7))

    fntsze = 2*13
    figA = plt.figure("Diffraction")
    axA = figA.add_subplot(111)

    residLst = result.residual

    #axA.annotate(r'$\Delta d + d = %d $' % k, xy=(0, maxofplt), fontsize=fntsze)
    #axA.plot(range(len(dataLst)), dataLst, '-', range(len(dataLst)), y_fitLst, '-')
    #axA.plot(range(len(dataLst)), dataLst, 'b.')
    #axA.plot(range(len(y_fitLst)), y_fitLst,'r.')
    
    posboot = [posdata[p] for p in listOfboots]

    axA.errorbar(np.array(posboot)/10., dataLst, yerr=statuncert,fmt='.',ecolor='green')#, capsize=4)
    #axA.plot(np.array(listOfboots)/10., dataLst, 'b.')
    axA.plot(np.array(posboot)/10., y_fitLst,'r.')
    
    #axA.plot(np.array(listOfboots)/10., statuncert,'g')

    #axA.plot(listOfboots, np.array(dataLst)-np.array(y_fitLst),'g.')

    #axA.set_xlim(-7, 7)
    #axA.set_ylim(-10, maxofplt+10)
    axA.set_xlabel(r'Undulator position $\Delta$ d [cm]',fontsize=fntsze)
    axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
    axA.tick_params(axis='both', which='major', labelsize=fntsze)
    axA.legend(['Fit','Data'],fontsize=fntsze)
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

        Sourcelst = UnDi.plainintensity(zetaprime2,zetaxx,Ytran,listOfboots[Strpidx*everynth],AUndu.z0,WLofPixel,LambdaRfit,delta0fit+posdata[listOfboots[Strpidx]]*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit,elevationfit)

        SourcelstIntensity1D = 55000.*abs(np.array(Sourcelst))[0]

        #y_fit = UnDi.diffractionintensityInitialFit_X(zetaprime2,zetaxx,Ytran,listOfboots[Strpidx*everynth],AUndu.z0,WLofPixel,LambdaRfit,delta0fit+listOfboots[Strpidx*everynth]*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit,elevationfit,whichdata,cfit)
        
        y_fit = UnDi.diffractionintensityInitialFit_X(zetaprime2,zetaxx,Ytran,listOfboots[Strpidx],AUndu.z0,WLofPixel,LambdaRfit,delta0fit+posdata[listOfboots[Strpidx]]*.001,gammafit,Amplfit,bfit,zpfit,Yscreenfit,elevationfit,whichdata,cfit)


        #print(residualbeforeplot)

        #print(innerresidual)

        #maxofdata = np.around(max(dataBooted[Strpidx, :]), decimals=-1)
        #maxofSourcelstIntensity1D = np.around(max(SourcelstIntensity1D), decimals=-1)
        #maxofplt = max([maxofdata,maxofSourcelstIntensity1D])
        
        y_fitLst.append(y_fit)

        fntsze = 13

        residLst = np.array(residLst)
        residimg = residLst.reshape((NumOfSeries, 1000))[Strpidx*everynth, :]#nptp(residLst.reshape((NumOfSeries, 1000)))[Strpidx, :]
        print(residLst)
        

        figA = plt.figure("Diffraction")
        axA = figA.add_subplot(111)
        figB = plt.figure("Sources")
        axB = figB.add_subplot(111)
        #axA.annotate(r'$\Delta d + d = %d $' % k, xy=(0, maxofplt), fontsize=fntsze)
        axA.set_title("Data and Fit and resid")
        axA.plot(1E6*yyy[0:1000]/AUndu.z0, dataBooted[Strpidx*everynth, :], 'b-', 1E6*yyy/AUndu.z0, y_fit, 'r-')
        
        #axA.plot(1E6*yyy[0:1000]/AUndu.z0, (np.array(pltweights[-(NumOfSeries+1):-1]).reshape((NumOfSeries, 1000)))[Strpidx],'k')
        #axA.set_xlim(-7, 7)
        #axA.set_ylim(-10, maxofplt+10)
        #axA.set_xlabel('Position on Aperture [mm]',fontsize=fntsze)
        axA.set_xlabel(r'$\Theta_y$ [$\mu$rad]',fontsize=fntsze)
        axA.set_ylabel('Intensity [a.u.]',fontsize=fntsze)
        axA.tick_params(axis='both', which='major', labelsize=fntsze)
        axA.legend(['Data','Fit'])



        axB.set_title("Intensity without diffraction")
        #axB.annotate(r'$\Delta d + d = %d $' % k, xy=(0, maxofplt), fontsize=fntsze)
        axB.plot(1E6*yyy[0:1000]/AUndu.z0, SourcelstIntensity1D, '-')
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