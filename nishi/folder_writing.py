###
# making folder. It mey not be worked after directry changed. It need not work otherwise.
###
import os

path ="C:/Users/24nis/Documents/data/Beam20200320_1/Energymeasurement_2/raw_shortExp/All_Modular/"

chosenfolders = [
                "SharedP_20210401_825_825_seed10000_new_zpfit_boot000425",
                "SharedP_20210401_825_825_seed10000_new_zpfit_boot000825",
                "SharedP_20210401_825_825_seed10000_new_zpfit_boot400825",
                "SharedP_20210401_825_825_seed10000_new_zpfit_delta0fit",
                "SharedP_20210401_825_825_seed10000_new_zpfit_MultiWL",
                "SharedP_20210401_825_825_seed10000_new_zpfit_Yscreen"
                ]

targetfiles = ["Amplfit.csv",
               "Ytran.csv",
               "gammafit.csv",
               "delta0fit.csv",
               "cfit.csv",
               "Amplfitstderr.csv",
               "Ytranstderr.csv",
               "gammafitstderr.csv",
               "delta0fitstderr.csv",
               "cfitstderr.csv",
               "amplYtranCorrel.csv",
               "amplgammafitCorrel.csv",
               "ampldelta0fitCorrel.csv",
               "amplcfitCorrel.csv",
               "YtrangammafitCorrel.csv",
               "Ytrandelta0fitCorrel.csv",
               "YtrancfitCorrel.csv",
               "gammadelta0fitCorrel.csv",
               "gammacfitCorrel.csv",
               "delta0fitcfitCorrel.csv",
               "zpfit.csv",
               "Residualmean.csv"]

for Numdata in range(len(chosenfolders)):
    for i in range(len(targetfiles)):
        f = open(path + chosenfolders[Numdata] + "/" + targetfiles[i],"w")
        f.close()
