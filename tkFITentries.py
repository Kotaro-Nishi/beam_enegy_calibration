# path to data set
loadstring = tk.StringVar()
loadstring.set(ConfigLst[0])
loadstring_entry = ttk.Entry(frame, width=100, textvariable=loadstring)
loadstring_entry.grid(column=1, row=0, **options)
loadstring_entry.focus()

ConstantsColumn = 1
ConstantsRow = 0

Entry1 = Label(frame2, text = "NumOfSeries").grid(column=ConstantsColumn-1, row=ConstantsRow, **options)

NumOfSeries_string = tk.StringVar()
NumOfSeries_string.set(ConfigLst[1])
NumOfSeries_string_entry = ttk.Entry(frame2, width=10, textvariable=NumOfSeries_string)
NumOfSeries_string_entry.grid(column=ConstantsColumn, row=ConstantsRow, **options)
NumOfSeries_string_entry.focus()

Entry2 = Label(frame2, text = "LengthOfData").grid(column=ConstantsColumn-1, row=ConstantsRow+1, **options)

LengthOfData_string = tk.StringVar()
LengthOfData_string.set(ConfigLst[2])
LengthOfData_string_entry = ttk.Entry(frame2, width=10, textvariable=LengthOfData_string)
LengthOfData_string_entry.grid(column=ConstantsColumn, row=ConstantsRow+1, **options)
LengthOfData_string_entry.focus()
LengthOfData_string_entry.config(state="disabled") # disable text field, even shorter intervals must load everything

Entry3 = Label(frame2, text = "Interval").grid(column=ConstantsColumn-1, row=ConstantsRow+2, **options)

Interval_string1 = tk.StringVar()
Interval_string1.set(ConfigLst[3][0])
Interval_string1_entry = ttk.Entry(frame2, width=10, textvariable=Interval_string1)
Interval_string1_entry.grid(column=ConstantsColumn, row=ConstantsRow+2, **options)
Interval_string1_entry.focus()

Interval_string2 = tk.StringVar()
Interval_string2.set(ConfigLst[3][1])
Interval_string2_entry = ttk.Entry(frame2, width=10, textvariable=Interval_string2)
Interval_string2_entry.grid(column=ConstantsColumn+1, row=ConstantsRow+2, **options)
Interval_string2_entry.focus()

Entry4 = Label(frame2, text = "Amplitude").grid(column=ConstantsColumn-1, row=ConstantsRow+3, **options)
Entry5 = Label(frame2, text = "Ytran").grid(column=ConstantsColumn+1, row=ConstantsRow+3, **options)
Entry6 = Label(frame2, text = "gamma").grid(column=ConstantsColumn+3, row=ConstantsRow+3, **options)
Entry7 = Label(frame2, text = "delta0fit").grid(column=ConstantsColumn+5, row=ConstantsRow+3, **options)
Entry8 = Label(frame2, text = "zpfit").grid(column=ConstantsColumn+7, row=ConstantsRow+3, **options)
Entry9 = Label(frame2, text = "elevation").grid(column=ConstantsColumn+9, row=ConstantsRow+3, **options)
Entry10 = Label(frame2, text = "cfit").grid(column=ConstantsColumn+11, row=ConstantsRow+3, **options)

Guess_string1 = tk.StringVar()
Guess_string1.set(ConfigLst[4][0])
Guess_string1_entry = ttk.Entry(frame2, width=10, textvariable=Guess_string1)
Guess_string1_entry.grid(column=ConstantsColumn, row=ConstantsRow+3, **options)
Guess_string1_entry.focus()

Guess_string2 = tk.StringVar()
Guess_string2.set(ConfigLst[4][1])
Guess_string2_entry = ttk.Entry(frame2, width=10, textvariable=Guess_string2)
Guess_string2_entry.grid(column=ConstantsColumn+2, row=ConstantsRow+3, **options)
Guess_string2_entry.focus()

Guess_string3 = tk.StringVar()
Guess_string3.set(ConfigLst[4][2])
Guess_string3_entry = ttk.Entry(frame2, width=10, textvariable=Guess_string3)
Guess_string3_entry.grid(column=ConstantsColumn+4, row=ConstantsRow+3, **options)
Guess_string3_entry.focus()

Guess_string4 = tk.StringVar()
Guess_string4.set(ConfigLst[4][3])
Guess_string4_entry = ttk.Entry(frame2, width=10, textvariable=Guess_string4)
Guess_string4_entry.grid(column=ConstantsColumn+6, row=ConstantsRow+3, **options)
Guess_string4_entry.focus()

Guess_string5 = tk.StringVar()
Guess_string5.set(ConfigLst[4][4])
Guess_string5_entry = ttk.Entry(frame2, width=10, textvariable=Guess_string5)
Guess_string5_entry.grid(column=ConstantsColumn+8, row=ConstantsRow+3, **options)
Guess_string5_entry.focus()

Guess_string6 = tk.StringVar()
Guess_string6.set(ConfigLst[4][5])
Guess_string6_entry = ttk.Entry(frame2, width=10, textvariable=Guess_string6)
Guess_string6_entry.grid(column=ConstantsColumn+10, row=ConstantsRow+3, **options)
Guess_string6_entry.focus()

Guess_string7 = tk.StringVar()
Guess_string7.set(ConfigLst[4][6])
Guess_string7_entry = ttk.Entry(frame2, width=10, textvariable=Guess_string7)
Guess_string7_entry.grid(column=ConstantsColumn+12, row=ConstantsRow+3, **options)
Guess_string7_entry.focus()

Entry11 = Label(frame2, text = "zetax").grid(column=ConstantsColumn-1, row=ConstantsRow+4, **options)

zetax_string = tk.StringVar()
zetax_string.set(ConfigLst[5])
zetax_string_entry = ttk.Entry(frame2, width=10, textvariable=zetax_string)
zetax_string_entry.grid(column=ConstantsColumn, row=ConstantsRow+4, **options)
zetax_string_entry.focus()

Entry12 = Label(frame2, text = "Yscreen").grid(column=ConstantsColumn-1, row=ConstantsRow+5, **options)

Yscreen_string = tk.StringVar()
Yscreen_string.set(ConfigLst[6])
Yscreen_string_entry = ttk.Entry(frame2, width=10, textvariable=Yscreen_string)
Yscreen_string_entry.grid(column=ConstantsColumn, row=ConstantsRow+5, **options)
Yscreen_string_entry.focus()

Entry13 = Label(frame2, text = "lambdaR").grid(column=ConstantsColumn-1, row=ConstantsRow+6, **options)

lambdaR_string = tk.StringVar()
lambdaR_string.set(ConfigLst[7])
lambdaR_string_entry = ttk.Entry(frame2, width=10, textvariable=lambdaR_string)
lambdaR_string_entry.grid(column=ConstantsColumn, row=ConstantsRow+6, **options)
lambdaR_string_entry.focus()

Entry14 = Label(frame2, text = "dispersion").grid(column=ConstantsColumn-1, row=ConstantsRow+7, **options)

dispersion_string = tk.StringVar()
dispersion_string.set(ConfigLst[8])
dispersion_string_entry = ttk.Entry(frame2, width=10, textvariable=dispersion_string)
dispersion_string_entry.grid(column=ConstantsColumn, row=ConstantsRow+7, **options)
dispersion_string_entry.focus()

Entry15 = Label(frame2, text = "pixelof404nm").grid(column=ConstantsColumn-1, row=ConstantsRow+8, **options)

pixelof404nm_string = tk.StringVar()
pixelof404nm_string.set(ConfigLst[9])
pixelof404nm_string_entry = ttk.Entry(frame2, width=10, textvariable=pixelof404nm_string)
pixelof404nm_string_entry.grid(column=ConstantsColumn, row=ConstantsRow+8, **options)
pixelof404nm_string_entry.focus()

Entry16 = Label(frame2, text = "DoNotSave").grid(column=ConstantsColumn-1, row=ConstantsRow+9, **options)

DoNotSave_string = tk.StringVar()
DoNotSave_string.set(ConfigLst[10])
DoNotSave_string_entry = ttk.Entry(frame2, width=10, textvariable=DoNotSave_string)
DoNotSave_string_entry.grid(column=ConstantsColumn, row=ConstantsRow+9, **options)
DoNotSave_string_entry.focus()

Entry17 = Label(frame2, text = "randomseed").grid(column=ConstantsColumn-1, row=ConstantsRow+10, **options)

randomseed_string = tk.StringVar()
randomseed_string.set(ConfigLst[11])
randomseed_string_entry = ttk.Entry(frame2, width=10, textvariable=randomseed_string)
randomseed_string_entry.grid(column=ConstantsColumn, row=ConstantsRow+10, **options)
randomseed_string_entry.focus()
