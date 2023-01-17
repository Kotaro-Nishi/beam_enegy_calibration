import numpy as np
import matplotlib.pyplot as plt

#D:\Doktor\BeamtimeFeb2020_Common\data\Beam20200319_1\Energymeasurement_1\HG\Mercury202003191818\Spectrum\Pos0
disp0 = 0.004807824803547447
disp0std = 3.9062962007447854e-07
loca0 = 1315.0000243860447
loca0std = 0.7341846176268909

#D:\Doktor\BeamtimeFeb2020_Common\data\Beam20200319_1\Energymeasurement_3\HG\Mercury202003192333\Spectrum\Pos0
disp1 = 0.004808888252396024
disp1std = 1.9062655655885548e-07
loca1 = 1315.7957379048294
loca1std = 1.0193786775716174

#D:\Doktor\BeamtimeFeb2020_Common\data\Beam20200320_1\Energymeasurement_1\HG\Mercury202003201438\Spectrum\Pos0
disp2 = 0.004808784982620716
disp2std = 2.514020840928887e-07
loca2 = 1315.33792502413
loca2std = 1.0167541956632948

#D:\Doktor\BeamtimeFeb2020_Common\data\Beam20200320_1\Energymeasurement_2\HG\Mercury202003202154\Spectrum\Pos0
disp3 = 0.004808759252942803
disp3std = 1.5629889792630894e-07
loca3 = 1315.3024745554976
loca3std = 1.0199041435307357

disp02lst = [disp0,disp1,disp2,disp3]
disp02stdlst = [disp0std,disp1std,disp2std,disp3std]
loca02lst = [loca0,loca1,loca2,loca3]
locastd02lst = [loca0std,loca1std,loca2std,loca3std]


#D:\Doktor\BeamtimeFeb2020_CommonJune\data\Beam20200702_1\Energymeasurement_7\HG
disp4 = 0.0047418919185504135
disp4std = 0.00012768303403048868
loca4 = 1244.8712801693136
loca4std = 0.8339000168254843

#D:\Doktor\BeamtimeFeb2020_CommonJune\data\Beam20200702_1\Hg-spektrumOhneLaserBlende
disp5 = 0.004784440666783993
disp5std = 3.753362973095752e-07
loca5 = 1244.4881130982933
loca5std = 0.5904886038498656

#D:\Doktor\BeamtimeFeb2020_CommonJune\data\Beam20200703_1\Energymeasurement_1\HG
disp6 = 0.00484885086891845
disp6std = 0.00019250962236681726
loca6 = 1244.7785595636628
loca6std = 0.8339075005626029

#D:\Doktor\BeamtimeFeb2020_CommonJune\data\Beam20200703_1\Energymeasurement_3\HG
disp7 = 0.004655975789964535
disp7std = 0.0005692242626492688
loca7 = 1244.7680267618518
loca7std = 0.8357238165201283

#D:\Doktor\BeamtimeFeb2020_CommonJune\data\Beam20200703_1\Energymeasurement_5\HG
disp8 = 0.004784707315424703
disp8std = 3.4242793529409127e-07
loca8 = 1240.3696567425088
loca8std = 0.8205751669134076

disp07lst = [disp4,disp5,disp6,disp7,disp8]
disp07stdlst = [disp4std,disp5std,disp6std,disp7std,disp8std]
loca07lst = [loca4,loca5,loca6,loca7,loca8]
locastd07lst = [loca4std,loca5std,loca6std,loca7std,loca8std]

print(sum(disp02lst)/len(disp02lst))
print(sum(loca02lst)/len(loca02lst))
print(np.std(loca02lst,axis=0, ddof=1))

print(sum(disp07lst)/len(disp07lst))
print(sum(loca07lst)/len(loca07lst))
print(np.std(loca07lst,axis=0, ddof=1))