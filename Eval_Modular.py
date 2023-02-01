import tkinter as tk
from tkinter import ttk
from tkinter import Menu
from tkinter import Label
from tkinter.messagebox import showerror

import ModularFunctions as MoFu

from tkinter import filedialog as fd

####################################

import matplotlib.pyplot as plt
import numpy as np
from numpy import transpose as nptp
import csv
import os, sys
from natsort import natsorted, ns
from scipy.optimize import curve_fit
import pandas as pd

from scipy.spatial import distance

import random
import statistics

global dp
global Targetfolder
global mode
global posdata
Targetfolder = "default"
mode = "default"

import Analyse_functions as Ana

dp = 1 # do plot instant

ConfigLst = ["//data/Mainz/BeamEnergyCalib/data/Beam20200320_1/Energymeasurement_2/raw_shortExp",825,825,[0,825],[1.4,.0001,352.1,.1870,5.,11.,1.],0,330E-6,412E-9,0.0048088E-9,1315.3024,0,10000]

####################################

#root window
root = tk.Tk()
root.title("Choose Data")
root.geometry("1020x570")
root.resizable(False,False)

def fahrenheit_to_celsius(f):
    """Convert fahrenheit to celsius"""
    return (f-32)*5/9

frame  = ttk.Frame(root)
frame2 = ttk.Frame(root)

options = {"padx": 5, "pady": 5}

# temperature label
#temperature_label = ttk.Label(frame, text="Folder")
#temperature_label.grid(column=0, row=0, sticky="W", **options)

exec(open("tkFITentries.py").read())

# Loading window


def load_button_clicked():
    """ Button click event"""

    try:
        Folder = loadstring.get()
        print(Folder)

        #os.system("python ./Analyse_Gamma_Modular.py " + Folder)

        #exec(open("Analyse_Gamma_Modular.py ").read())

        Ana.ld("*gammafit*",Folder)
        Ana.ld("*delta0fit*",Folder)
        gWL = Ana.cl("gammafitSharedP_20210401_825_825_seed10000_new_zpfit_MultiWL*")
        gg0 = Ana.cl("gammafitSharedP_20210401_825_825_seed10000_new_zpfit_boot000425*")
        gg4 = Ana.cl("gammafitSharedP_20210401_825_825_seed10000_new_zpfit_boot400825*")
        gg8 = Ana.cl("gammafitSharedP_20210401_825_825_seed10000_new_zpfit_boot000825*")

        dd0 = Ana.cl("delta0fitSharedP_20210401_825_825_seed10000_new_zpfit_boot000425*")
        dd4 = Ana.cl("delta0fitSharedP_20210401_825_825_seed10000_new_zpfit_boot400825*")
        dd8 = Ana.cl("delta0fitSharedP_20210401_825_825_seed10000_new_zpfit_boot000825*")

        Ana.ht(gg0,20,"g")
        Ana.ht(gg4,20,"g")
        Ana.ht(gg8,20,"g")

        Ana.cp(gg0,dd0)
        Ana.cp(gg4,dd4)
        Ana.cp(gg8,dd8)

        Ana.cp([gg0[0],gg4[0],gg8[0]],[dd0[0],dd4[0],dd8[0]])

        Ana.pl(gWL)

        print("Systematic error of fit: ",.511/np.sqrt(2)*np.sqrt((np.mean(gg0)-np.mean(gg8))**2+(np.mean(gg4)-np.mean(gg8))**2))
        print("mean Energy 000425: ",.511 * np.mean(gg0))
        print("mean gamma 000425: ",np.mean(gg0))
        print("mean Energy 400825: ",.511 * np.mean(gg4))
        print("mean gamma 400825: ",np.mean(gg4))
        print("mean Energy 000825: ",.511 * np.mean(gg8))
        print("mean gamma 000825: ",np.mean(gg8))
        
        print("mean Energy WL: ",.511 * np.mean(gWL))
        print("mean gamma WL: ",np.mean(gWL))
        print("std gamma WL: ",np.std(gWL))

        #Ana.cp(gg4,20,"g")
        #Ana.cp(gg8,20,"g")

        #loadstring_entry.config(state="disabled") # disable text field


        #c = fahrenheit_to_celsius(f)
        #result = f'{f} Fahrenheit = {c:.2f} Celsius'
        #result_label.config(text=result)
    except ValueError as error:
        showerror(title="Error", message=error)

def plot_button_clicked():
    """ Button click event """
    
    print("now plot")

def fit_button_clicked():
    Folder = loadstring.get()
    NumOfSeries = int(NumOfSeries_string.get())
    LengthOfData = int(LengthOfData_string.get())

    Interval = [int(Interval_string1.get()), int(Interval_string2.get())]

    Guessvalueslmfit = [
    float(Guess_string1.get()),
    float(Guess_string2.get()),
    float(Guess_string3.get()),
    float(Guess_string4.get()),
    float(Guess_string5.get()),
    float(Guess_string6.get()),
    float(Guess_string7.get())]

    zetax = float(zetax_string.get())
    
    MoFu.zetaxx = np.array([zetax])
    
    Yscreen = float(Yscreen_string.get())
    MoFu.Yscreenfit = Yscreen
    
    lambdaR = float(lambdaR_string.get())
    MoFu.LambdaRfit = lambdaR
    MoFu.dispersion = float(dispersion_string.get())
    MoFu.pixelof404nm = float(pixelof404nm_string.get())
    MoFu.DoNotSave = int(DoNotSave_string.get())
    MoFu.randomseed = int(randomseed_string.get())

    MoFu.DataPath = Folder
    MoFu.SaveSharedP = Targetfolder

    if Targetfolder == "default":
        print("specify Targetfolder")
    else:
        print("Targetfolder: ", Targetfolder)
        exec(open("ModularEvaluation.py ").read())


"""
def load_button_2_clicked():
    
    try:
        f = float(loadstring.get())
        c = fahrenheit_to_celsius(f)
        result = f'{f} Fahrenheit = {c:.2f} Celsius'
        result_label.config(text=result)
    except ValueError as error:
        showerror(title="Error", message=error)
"""

# load button
load_button = ttk.Button(frame, text="Load")
load_button.grid(column=2, row=0, sticky="W", **options)
load_button.configure(command=load_button_clicked)

# plot button
plot_button = ttk.Button(frame, text="plot")
plot_button.grid(column=2, row=1, sticky="W", **options)
plot_button.configure(command=plot_button_clicked)

# fit button
fit_button = ttk.Button(frame,text="fit")
fit_button.grid(column=3, row=0, sticky="W", **options)
fit_button.configure(command=fit_button_clicked)

# result label
result_label = ttk.Label(frame)
result_label.grid(row=1, columnspan=3, **options)


#SaveSharedP = "\All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_MultiWL_sing"

def donothing():
    Entry4 = Label(frame, text = "Defaultfolder").grid(column=1, row=1, **options)

def Lambda_L000825():
    Entry4 = Label(frame, text = "Target folder: \All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_MultiWL").grid(column=1, row=1, **options)
    global Targetfolder
    global mode
    Targetfolder = "\All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_MultiWL"
    mode = "lamb"

def Interval000825():
    Entry4 = Label(frame, text = "Target folder: \All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_boot000825").grid(column=1, row=1, **options)
    global Targetfolder
    global mode
    Targetfolder = "\All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_boot000825"
    mode = "posi"

def Interval000425():
    Entry4 = Label(frame, text = "Target folder: \All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_boot000425").grid(column=1, row=1, **options)
    global Targetfolder
    global mode
    Targetfolder = "\All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_boot000425"
    mode = "posi"

def Interval400825():
    Entry4 = Label(frame, text = "Target folder: \All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_boot400825").grid(column=1, row=1, **options)
    global Targetfolder
    global mode
    Targetfolder = "\All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_boot400825"
    mode = "posi"

def YscreenCMD():
    Entry4 = Label(frame, text = "Target folder: \All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_Yscreen").grid(column=1, row=1, **options)
    global Targetfolder
    global mode
    Targetfolder = "\All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_Yscreen"
    mode = "Yscr"

def phaseCMD():
    Entry4 = Label(frame, text = "Target folder: \All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_delta0fit").grid(column=1, row=1, **options)
    global Targetfolder
    global mode
    Targetfolder = "\All_Modular\SharedP_20210401_825_825_seed10000_new_zpfit_delta0fit"
    mode = "delt"

def select_file():
    filetypes = (
        ('text files', '*.csv'),
        ('All files', '*.*')
    )

    filename = fd.askopenfilename(
        title='Open a file',
        initialdir='/',
        filetypes=filetypes)
    print(filename)

    with open(filename, "r") as csv_file:
        csvdata = csv.reader(csv_file, delimiter=',')

        ConfigLst = list(csvdata)
        ConfigLst = ConfigLst[0]

        print("ConfigLst: ",ConfigLst)

        loadstring.set(ConfigLst[0])
        NumOfSeries_string.set(int(ConfigLst[1]))
        LengthOfData_string.set(int(ConfigLst[2]))
        
        Interval_string1.set(int(ConfigLst[3]))
        Interval_string2.set(int(ConfigLst[4]))

        Guess_string1.set(float(ConfigLst[5]))
        Guess_string2.set(float(ConfigLst[6]))
        Guess_string3.set(float(ConfigLst[7]))
        Guess_string4.set(float(ConfigLst[8]))
        Guess_string5.set(float(ConfigLst[9]))
        Guess_string6.set(float(ConfigLst[10]))
        Guess_string7.set(float(ConfigLst[11]))

        zetax_string.set(float(ConfigLst[12]))
        Yscreen_string.set(float(ConfigLst[13]))
        lambdaR_string.set(float(ConfigLst[14]))
        dispersion_string.set(float(ConfigLst[15]))
        pixelof404nm_string.set(float(ConfigLst[16]))
        DoNotSave_string.set(int(ConfigLst[17]))
        randomseed_string.set(int(ConfigLst[18]))

    return ConfigLst

menubar = Menu(root)
filemenu = Menu(menubar, tearoff=0)
filemenu.add_command(label="New", command=donothing)

ConfigLst = filemenu.add_command(label="Open", command=select_file)

filemenu.add_command(label="Save", command=donothing)
filemenu.add_command(label="Save as...", command=donothing)

filemenu.add_separator()

filemenu.add_command(label="Exit", command=root.quit)
menubar.add_cascade(label="File", menu=filemenu)
fitmenu = Menu(menubar, tearoff=0)

fitmenu.add_command(label="Lambda_L000825", command=Lambda_L000825)
fitmenu.add_command(label="Interval000825", command=Interval000825)
fitmenu.add_command(label="Interval000425", command=Interval000425)
fitmenu.add_command(label="Interval400825", command=Interval400825)

fitmenu.add_separator()
fitmenu.add_command(label="Yscreen", command=YscreenCMD)
fitmenu.add_command(label="phase", command=phaseCMD)

menubar.add_cascade(label="Fit options", menu=fitmenu)
helpmenu = Menu(menubar, tearoff=0)
helpmenu.add_command(label="Help Index", command=donothing)
helpmenu.add_command(label="About...", command=donothing)
menubar.add_cascade(label="Help", menu=helpmenu)

root.config(menu=menubar)
#root.mainloop()







frame.grid(padx=10, pady=10)
frame2.grid(padx=10, pady=10)
root.mainloop()