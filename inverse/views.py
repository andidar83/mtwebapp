from django.shortcuts import render
from . import views

def index(request):
    context = {
        'greetings':'Welcome to forward modelling tool'
        }
   
   
    return render(request,'inverse/index.html', context)
import math
import cmath
import time
import numpy as np
import matplotlib.pyplot as plt
import base64
from io import BytesIO
import matplotlib
matplotlib.use('Agg')
import numpy as np


import matplotlib.pyplot as plt

from matplotlib.figure import Figure
import matplotlib.lines as line
from scipy.integrate import simps




def index(request):
    context = {
        'greetings':'Welcome to forward modelling tool'
        }
    return render(request,'inverse/index.html', context)


def forward (res,thic, frekuensi):


    mu = 4*math.pi*1E-7
    i = 1j
    apparentresistivity = []
    phase = []
    for f in frekuensi:
        w = 2*math.pi*f
        Zn = cmath.sqrt(i*w*mu*res[-1])
      
        for r,h in zip(reversed(res[:-1]),reversed(thic[:-1])):
            y = cmath.sqrt((i*w*mu)/r)
            kj = y*r
            Ej = np.exp(-2*y*h)
            Rj = (kj-Zn)/(kj+Zn)
            Zj = kj*((1-Rj*Ej)/(1+Rj*Ej))
            
            Zn = Zj
        apparentresis = (abs(Zn)**2)/(w*mu)
        fase = math.atan(Zn.imag/ Zn.real)
        phase.append((np.rad2deg(fase)))
        apparentresistivity.append(apparentresis)
    appresphase = apparentresistivity+phase
    appresphasemisah = [apparentresistivity,phase]
    return { "Apparentsesistivity" : apparentresistivity,
            "Phase" : phase,
            "appresphase":appresphase,
            "appresphasemisah":appresphasemisah
        }


def depthplot(thicknesses, Resistivities):

    depth = []
    Resistivitiess = Resistivities
    elev = 0
    for i in thicknesses:
        depth.append(elev)
        elev = elev + i/1000
    Resistivitiess.append(Resistivities[-1])
    depth.append(elev)
    out = [Resistivitiess,depth]
    return out


def add(request):
    resis = str(request.GET['var1'])
    thic = str(request.GET['var2'])
    thiclow = int(request.GET['var3'])
    thichi = int(request.GET['var4'])
    pangkatfrek = np.linspace(thiclow,thichi,60)
    frekuensi = 10**pangkatfrek

    resis = list(resis.split(","))
    thic = list(thic.split(","))

    resis = list(map(int, resis))
    thic =     list(map(int, thic))


    calculated = forward (resis,thic,frekuensi)["appresphasemisah"]
    Calappres= calculated[0]
    Calphase = calculated[1]
    

    fig = plt.figure(figsize=(12,8))
    ax1 = plt.subplot(221)
    ax1.plot(frekuensi,Calappres,'b''d')

    ax1.set_xscale('log')
    # ax1.set_yscale('log')
    ax1.invert_xaxis()
    ax1.set_xlabel("Frrekuensi(Hz)")
    ax1.set_ylabel("Resistivitas Apparent (Ohmm)")
    
       
    ax2 = plt.subplot(223)
    ax2.plot(frekuensi,Calphase,'g''d')


    ax2.set_xscale('log')
    # ax2.set_yscale('log')
    ax2.invert_xaxis()
    ax2.set_xlabel("Frrekuensi(Hz)")
    ax2.set_ylabel("Fase (derajat)")
    

    sounding = depthplot(thic, resis)
    
    

       
    ax3 = plt.subplot(122)
    ax3.step(sounding[0],sounding[1],'r')
   
    
    ax3.set_xscale('log')
    ax3.set_xlim(0.1,10000)
    ax3.set_ylim(0,7)
    ax3.set_xlabel('Resistvity (Ohmm)')
    ax3.set_ylabel('Depth (km)')
    plt.gca().invert_yaxis()
    plt.suptitle("Forward Modelling Result",
              size=20)
    plt.grid()

    plt.show()

    fig=plt.gcf()
    plt.close()
    buffer = BytesIO()
    fig.savefig(buffer, format='png')
    buffer.seek(0)
    image_png = buffer.getvalue()
    graph = base64.b64encode(image_png)
    graph = graph.decode('utf-8')
    buffer.close()

    out1={'graph':graph}
    


    

    return render(request,'inverse/output.html', out1)


