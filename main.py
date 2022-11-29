# coding=utf-8
#! /usr/bin/python
# -*- coding: utf-8 -*-
import cmath
import math
#import scipy
import time
import sys
import datetime
import os
#import matplotlib.pyplot as plt
from scipy.integrate import odeint,quad
#from scipy import integrate
from datetime import datetime
from datetime import timedelta
from sys import argv as s
from decimal import *
import time
import re
import glob
import numpy as np
from matplotlib.widgets import SpanSelector
#from astropy.io import ascii
from scipy import interpolate
from scipy.integrate import simps
from scipy.signal import savgol_filter



import matplotlib.pyplot as plt
import matplotlib.colors

Filename="Zusammenfassung.plt"
s=[0,Filename]

#Polare = open(../Polaren_MoProMa_D2.plt)
#Polare = Polare.readlines()
PolarMatrix=np.genfromtxt("../Polaren_MoProMa_D2.plt",skip_header=58,delimiter="   ")
PolarMatrixTurb=np.genfromtxt("../Polaren_MoProMa_D2_Turb.plt",skip_header=12,delimiter="   ")
PolarMatrixvollTurb=np.genfromtxt("../Polaren_MoProMa_D2_vollTurb.plt",skip_header=12,delimiter="   ")
#print(PolarMatrix)

Messdatei = open(Filename, "r");
Messdaten= Messdatei.readlines()
Anzahl_Zeilen=len(Messdaten)
Anzahl_Spalten=len(Messdaten[102].strip(" ").rstrip().split(" "))
print("Anzahl_Zeilen: "+str(Anzahl_Zeilen))
print("Anzahl_Spalten: "+str(Anzahl_Spalten))
Matrix = np.zeros((Anzahl_Zeilen+1, Anzahl_Spalten+1), dtype=float);
#Matrix2 = np.zeros((25, 3), dtype=float);
print(Anzahl_Spalten)
i=0

Nue=0.0000153 # dyn. Viskosität Luft


for line in Messdaten[102:]:
    #print(line)
    i = i + 1
    Zeile = line.strip(" ").rstrip();

    Daten = (re.split(' ', Zeile));
    x = 0
    # print len(Daten)
    t = 0

    for t in range(0, len(Daten), 1):
        # print t
        Matrix[i][t] = float(Daten[t])
    #print(Matrix[i][96])

print("Schleife1)")
Event1=2693-102 #Flug1
#Event1=293-102 #Flug2
DauerEvent1=270
corr_Stat=np.zeros((5), dtype=float);
i=0
for i in range(0,5,1):
    print(i)
    #corr_Stat[i] = np.average(Matrix[Event1:Event1+DauerEvent1, 59+i]) - (np.average(Matrix[Event1:Event1+DauerEvent1, 59:96]))# )-(np.average(Matrix[Event1:Event1 + DauerEvent1, 30])*100)
    corr_Stat[i] = np.average(Matrix[Event1:Event1 + DauerEvent1, 59 + i]) - (np.average(Matrix[Event1:Event1 + DauerEvent1, 30])*100)
#corr_Stat = [190.17403658,165.74055795,179.0052416,157.09651233,167.07247849]
print(corr_Stat)

corr_Druck=np.zeros((32), dtype=float);
i=0

for i in range(0,32,1):
    corr_Druck[i] = np.average(Matrix[Event1:Event1 + DauerEvent1, 64 + i]) - (np.average(Matrix[Event1:Event1 + DauerEvent1, 59+(round(i/8))]))+ corr_Stat[(round(i/8))]
    #corr_Druck[i] = np.average(Matrix[Event1:Event1+DauerEvent1, 64+i]) - (np.average(Matrix[Event1:Event1+DauerEvent1, 30])*100)
    #corr_Druck[i] = np.average(Matrix[Event1:Event1 + DauerEvent1, 64 + i]) -  (np.average((Matrix[Event1:Event1+DauerEvent1, 59:63])))-(np.average(Matrix[Event1:Event1 + DauerEvent1, 29])*100)
#corr_Druck=[ 95.73202667,-34.10747975,12.16482024,-20.42145691,19.10439334,17.80408065,14.0040659,59.25360661,29.73998967,21.74893957,71.49066965,11.2840679,35.11962241,56.33791815,-1.24539648,55.81432354,1.29391858,95.39758493,29.02471263,-2.57815286,18.00788563,-2.38528471,28.49974716,-12.7860598,-27.59317983,-2.71555622,-75.29582768,-13.80548915,10.80049292,-3.31985468, -13.69906763,  0.        ]
print(corr_Druck)
time.sleep(3)
#print(corr_Druck)
Drucksensoren=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32)
Statikdrucksensoren= (1,2,3,4,5)
Statikdrucksensoren2= (1,8,16,24,32)
#Mittelwert über 100Werte in einem Zeitabschnitt
Delta=50
#Anfang=56000 #Flug2
#Ende2=80000 #Flug2
Anfang=107#000
Anfang2=107
Ende2=130000
Ende=Anfang+Delta
h=0
Anzahl_Cws=int((Ende2-Anfang)/Delta)
Beiwerte = np.zeros((Anzahl_Cws, 2), dtype=float);
Abstand_Sensoren = np.arange(0,32,1, dtype=float)* 0.003;
#print(Abstand_Sensoren)
Event=0
Position=0
CW_Nachlauf=0.0
CA_ges=0.0
Abstand_Messpunkt = 0
Re_Number=0.0
while Ende<=Ende2 and Anfang>=Anfang2:

    Anfang=Anfang+Delta
    Ende=Anfang+Delta
    #print(Anfang)
    #Druckbeiwerte:
    z=0
    Druckwerte=np.zeros((32), dtype=float);



    EndDruckwerte=np.zeros((32),dtype=float);

    i=0
    Statikdruckwerte = np.zeros((5), dtype=float);
    Statikdruckwerte_Interpol = np.zeros((32), dtype=float);

    for i in range(0,5,1):
        Statikdruckwerte[i]= (np.average(Matrix[Anfang:Ende, 59 + i]) )- corr_Stat[i]
    #print(Statikdruckwerte)
    s = interpolate.InterpolatedUnivariateSpline(Statikdrucksensoren2, Statikdruckwerte)
    Statiknew = s(Drucksensoren)
    for i in range(0,32,1):
        #print(i/8," ",z," ",i)
        Druckwerte[i] = (np.average(Matrix[Anfang:Ende, 64 + i] - corr_Druck[i]) )-((Statiknew[i]))
        #Druckwerte[i]=(np.average(Matrix[Anfang:Ende, 64+i]-corr_Druck[i])-np.average(Matrix[Anfang:Ende, 59+z]-corr_Stat[z]))-(np.average(Matrix[Anfang:Ende, 64+31]-corr_Druck[31])-np.average(Matrix[Anfang:Ende, 59+4]-corr_Stat[4]))
        #Druckwerte[i] = (np.average(Matrix[Anfang:Ende, 64 + i] - corr_Druck[i])) #- np.average(Matrix[Anfang:Ende, 95])) #- (np.average(Matrix[Anfang:Ende, 29]*100) )
        if i==1:
            Druckwerte[0]=Druckwerte[1]
        if i>29:
            Druckwerte[i]=Druckwerte[29]
            #print(i," ist gleich 27")
        #Statikdruckwerte= (np.average(Matrix[Anfang:Ende, 59 + i] - corr_Stat[i]))
        #print((np.average(Matrix[Anfang:Ende, 64+31]-corr_Druck[31])-np.average(Matrix[Anfang:Ende, 59+5]-corr_Stat[4])))
        #print(np.average(Matrix[Anfang:Ende, 64+i]))
        #print(np.average(Matrix[Anfang:Ende, 59+z]))


    #s2 = interpolate.InterpolatedUnivariateSpline(Drucksensoren, Druckwerte)
    #xnew = Drucksensoren

    #Drucknew=s2(Drucksensoren)
    #yhat = savgol_filter(Druckwerte, 5, 3)  # window size 51, polynomial order 3
    yhat = Druckwerte
    #plt.plot(Matrix[21400,1], data[:,10],'k+',label="1");
    Bildname="Nachlaufdelle%d.jpg"%h
    IAS = math.sqrt(abs((np.average(Matrix[Anfang:Ende, 29])*100-10.1))*2/1.225)*3.6
    Altitude = np.average(Matrix[Anfang:Ende, 7])
    Alpha = np.average(Matrix[Anfang:Ende,33])
    Sekunden=np.average(Matrix[Anfang:Ende,58])
    Position=np.average(Matrix[Anfang:Ende,96])

    Event=np.average(Matrix[Anfang:Ende,97])
    #print(Position)
    print(Event)
    a=math.sqrt(math.pow(np.average(Matrix[Anfang:Ende,26]),2)+math.pow(np.average(Matrix[Anfang:Ende,27]),2)+math.pow(np.average(Matrix[Anfang:Ende,28]),2))
    m=560.0

    #print(CA_ges)
    #print("Beschleunigung"+str(a))
    #plt.plot(Druckwerte, Drucksensoren,'o-',label="Delle");
    #plt.plot(Statikdruckwerte[4]-Statikdruckwerte, Statikdrucksensoren, 'o-', label="Statikdruckverlauf");
    #plt.plot(Statikdruckwerte-Statikdruckwerte[4], Statikdrucksensoren2, 'o-', label="Statikdruckverlauf");
    #plt.plot((yhat-np.average(yhat))-((Statiknew)-np.average(Statiknew))-50, Drucksensoren, 'o-', label="Druckwerte_Interpol");
    #EndDruckwerte = yhat - (yhat[0])

    EndDruckwerte = yhat - ((yhat[25]+yhat[25])/2)#((np.average(Matrix[Anfang:Ende, 29])-0.0935) * 100)

    #print("Druckwert vom Stau: ",((np.average(Matrix[Anfang:Ende, 29])-0.0935) * 100))
    #print("Druckwert vom Rechen 0: ", (np.average(yhat[0]) ))
    IAS2=np.sqrt(abs(((Druckwerte[0]+Druckwerte[31]+130)/2)*2/1.225))*3.6
    CA_ges = (m * a) / (1.225 / 2 * math.pow(IAS / 3.6, 2) * 11.36)
    TrapezFlaeche=simps(-EndDruckwerte,dx=0.003)
    #TrapezFlaeche1 = np.trapz(-EndDruckwerte,dx=0.003)/(1.225/2*math.pow(IAS/3.6,2)*0.806)
    #TrapezSimps=simps(-EndDruckwerte,dx=0.003)/(1.225/2*math.pow(IAS/3.6,2)*0.806)
    if Position ==100:
        CW_Nachlauf=TrapezFlaeche/((1.225 / 2) * pow(IAS / 3.6, 2)*0.814)
        Abstand_Messpunkt = 1.291  # Abstand des Messpunktes vom SP in Metern
        Re_Number=IAS2/3.6*0.814/Nue

    if Position ==250:
        CW_Nachlauf=TrapezFlaeche/((1.225 / 2) * pow(IAS / 3.6, 2)*0.806)
        Abstand_Messpunkt = 1.441
        Re_Number = IAS2/3.6 * 0.806 / Nue
    if Position == 500:
        CW_Nachlauf = TrapezFlaeche / ((1.225 / 2) * pow(IAS/ 3.6, 2) * 0.801)
        Abstand_Messpunkt = 1.691
        Re_Number = IAS2/3.6 * 0.801 / Nue
    #print(CW_Nachlauf)
    V_Zusatz = np.average(Matrix[Anfang:Ende, 23]) * (2 * Abstand_Messpunkt * np.pi / 360)
    Alfa_Zusatz = np.degrees(np.arctan(V_Zusatz / IAS2 / 3.6))
    #print(CA_ges/CW_Nachlauf)

    if np.isnan(Event)==True:
        Event=0
    #if True:
    if Event>=3: #and Event<14:
    #if int(Event)==3 or int(Event)==5 or int(Event)==7 or int(Event)==9 or int(Event)==11 or int(Event)==13: #Flug 1
    #if int(Event) == 1 or int(Event) == 3 or int(Event) == 5 or int(Event) == 7 : #Flug 2
        h += 1
        #if Position != -1:
        if True:
            Beiwerte[h - 1][0] = CA_ges
            Beiwerte[h - 1][1] = CW_Nachlauf
        #print(Alfa_Zusatz)
        #print(Re_Number)
        #print(CA_ges)
        #print(TrapezFlaeche)
        print("IAS=",IAS)
        print("IAS2=",IAS2)
        print(((Druckwerte[0] + Druckwerte[31])/2)-(np.average(Matrix[Anfang:Ende, 29])*100-10.1)+130)
        #print(CW_Nachlauf)
        #print(TrapezSimps)
        print("\n")
        Alfa_Zusatz=math.atan((V_Zusatz/IAS/3.6))
        plt.plot((Druckwerte) - (Druckwerte[0]+Druckwerte[31])/2  , Drucksensoren, 'o-',label="Druckwerte_UR");

        plt.plot(EndDruckwerte, Drucksensoren, 'o-', label="Druckwerte_Interpol");

        plt.plot(Statiknew-((Statiknew[0]+Statiknew[31])/2), Drucksensoren, 'o-', label="Statikdruckverlauf_Interpol");
        plt.plot(Statikdruckwerte-((Statikdruckwerte[0]+Statikdruckwerte[4])/2), Statikdrucksensoren2, 'o-', label="Statikdruck_UR");
        plt.axis([200,-500,0,33])
        #plt.axis([20, -20, 0, 33])
    #plt.axis([80000, 102000, 0, 6])
        plt.title("IAS: %3.d km/h; CA: %1.2f, CW %1.5f, t=%4.1f s" %(IAS,CA_ges,CW_Nachlauf,Sekunden));
        plt.xlabel("Abw. vom Stau-/Statikdruck in Pa");
        plt.ylabel("Sensornummer [-]");
    #plt.grid(color='k', linestyle='-', linewidth=1)
        plt.grid(b=True, which='major', color='k', linestyle='-',linewidth=1)
        plt.minorticks_on()
        plt.grid(b=True, which='minor', color='k', linestyle='--',linewidth=0.25)
        plt.legend(bbox_to_anchor=(0.6, 0.9), loc=2, borderaxespad=0.)
        #plt.savefig(Bildname);
        plt.show(block=False)
        plt.pause(.2)
        plt.close()
#print(Beiwerte)
plt.figure(figsize=(10, 8), dpi=80)
plt.plot(Beiwerte[:,1], Beiwerte[:,0], 'x',color='k', label="Beiwerte");
plt.plot(PolarMatrix[:,2],PolarMatrix[:,1],linestyle='solid',color='k',label="XFoil Xtr-Bot=78%")
plt.plot(PolarMatrixTurb[:,2],PolarMatrixTurb[:,1],linestyle='dashed',color='g',label="XFoil Xtr-Bot=60%")
plt.plot(PolarMatrixvollTurb[:,2],PolarMatrixvollTurb[:,1],linestyle='dotted',color='k',label="XFoil vollTurb")
#plt.grid(b=True,color='k', linestyle='-', linewidth=1)
plt.grid(b=True, which='major', color='k', linestyle='-', linewidth=1)
plt.minorticks_on()
plt.legend(bbox_to_anchor=(0.6, 0.9), loc=2, borderaxespad=0.)
plt.grid(b=True, which='minor', color='k', linestyle='--', linewidth=0.25)
plt.axis([0,0.025,0,1.3])
#plt.show(block=False)
#plt.pause(2.0)
plt.savefig("CA-CW-Plot.jpg");
plt.close()