print("saving things")

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

print(DataPath + SaveSharedP)

if DoNotSave == 0:

    with open(DataPath + SaveSharedP + "/Amplfit.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([Amplfit])


    with open(DataPath + SaveSharedP + "/Ytran.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([Ytran])

    with open(DataPath + SaveSharedP + "/gammafit.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([gammafit])

    with open(DataPath + SaveSharedP + "/delta0fit.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([delta0fit])

    with open(DataPath + SaveSharedP + "/cfit.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([cfit])

    with open(DataPath + SaveSharedP + "/zpfit.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([zpfit])

    with open(DataPath + SaveSharedP + "/elevationfit.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([elevationfit])



    with open(DataPath + SaveSharedP + "/Amplfitstderr.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([Amplfitstderr])

    with open(DataPath + SaveSharedP + "/Ytranstderr.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([Ytranstderr])

    with open(DataPath + SaveSharedP + "/gammafitstderr.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([gammafitstderr])

    with open(DataPath + SaveSharedP + "/delta0fitstderr.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([delta0fitstderr])

    with open(DataPath + SaveSharedP + "/cfitstderr.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([cfitstderr])

    with open(DataPath + SaveSharedP + "/zpfitstderr.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([zpfitstderr])

    with open(DataPath + SaveSharedP + "/elevationfitstderr.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([elevationfitstderr])








    with open(DataPath + SaveSharedP + "/AmplfitYtranCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([AmplfitYtranCorrel])

    with open(DataPath + SaveSharedP + "/AmplfitgammafitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([AmplfitgammafitCorrel])
    
    with open(DataPath + SaveSharedP + "/Amplfitdelta0fitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([Amplfitdelta0fitCorrel])

    with open(DataPath + SaveSharedP + "/AmplfitcfitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([AmplfitcfitCorrel])

    with open(DataPath + SaveSharedP + "/AmplfitzpfitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([AmplfitzpfitCorrel])
        
    with open(DataPath + SaveSharedP + "/AmplfitelevationfitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([AmplfitelevationfitCorrel])
        
        
        
        
    with open(DataPath + SaveSharedP + "/YtrangammafitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([YtrangammafitCorrel])

    with open(DataPath + SaveSharedP + "/Ytrandelta0fitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([Ytrandelta0fitCorrel])

    with open(DataPath + SaveSharedP + "/YtrancfitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([YtrancfitCorrel])

    with open(DataPath + SaveSharedP + "/YtranzpfitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([YtranzpfitCorrel])

    with open(DataPath + SaveSharedP + "/YtranelevationfitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([YtranelevationfitCorrel])
        
        
        
        
    with open(DataPath + SaveSharedP + "/gammadelta0fitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([gammadelta0fitCorrel])

    with open(DataPath + SaveSharedP + "/gammacfitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([gammacfitCorrel])

    with open(DataPath + SaveSharedP + "/gammazpfitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([gammazpfitCorrel])

    with open(DataPath + SaveSharedP + "/gammaelevationfitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([gammaelevationfitCorrel])


        
    with open(DataPath + SaveSharedP + "/delta0fitcfitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([delta0fitcfitCorrel])

    with open(DataPath + SaveSharedP + "/delta0fitzpfitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([delta0fitzpfitCorrel])
        
    with open(DataPath + SaveSharedP + "/delta0fitelevationfitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([delta0fitelevationfitCorrel])



    with open(DataPath + SaveSharedP + "/cfitzpfitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([cfitzpfitCorrel])

    with open(DataPath + SaveSharedP + "/cfitelevationfitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([cfitelevationfitCorrel])



    with open(DataPath + SaveSharedP + "/zpfitelevationfitCorrel.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([zpfitelevationfitCorrel])


    residsqrlstmean = np.mean(residsqrlst)

    with open(DataPath + SaveSharedP + "/Residualmean.csv", 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow([residsqrlstmean])