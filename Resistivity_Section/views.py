from django.shortcuts import render
from . import views
import matplotlib.pyplot as plt
import base64
from io import BytesIO
import numpy as np
import matplotlib


from django.http import HttpResponse

from django.core.files.storage import FileSystemStorage

import scipy.interpolate as inter
from scipy.signal import convolve2d


import glob
import math
import cmath
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as line
import scipy.stats as stats
import pandas as pd
import random
import statistics
from scipy import optimize
import mtpy.utils.calculator as MTcc
import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging
import gstools as gs
from scipy.signal import lfilter
from scipy import signal as sg
import matplotlib

import itertools
import matplotlib as mpl


def forward (res,thic, frekuensi):
    import math
    import cmath
    import numpy as np
    

    

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

def Jacobian( res, thic, frekuensi):
    import numpy as np


    dh = 1.0e-7
    JR=[]
    for i in range(len(res)):
        res[i]= res[i]+dh
        resphasetemp=forward(res,thic,frekuensi)["appresphase"]
        res[i]= res[i]-dh
        resphasetemp0=forward(res,thic,frekuensi)["appresphase"]
        finitdif = (np.array(resphasetemp) - np.array(resphasetemp0))/dh
        JR.append(finitdif)
         
    
    
    JT=[]
    for i in range(len(thic)-1):
        thic[i]= thic[i]+dh
        resphasetemp=forward(res,thic,frekuensi)["appresphase"]
        thic[i]= thic[i]-dh
        resphasetemp0=forward(res,thic,frekuensi)["appresphase"]
        finitdif = (np.array(resphasetemp) - np.array(resphasetemp0))/dh
        JT.append(finitdif)  
    
    
    Jacobian=np.hstack((np.array(JR).transpose(),np.array(JT).transpose()))
    return Jacobian


def RMS(calculated,obs):

    import numpy as np
    jumlah=0
    for h,g in zip(calculated,obs):
        err = (h - g)**2#(((g-statistics.mean(D))**2)/len(D))**0.5
        jumlah = jumlah + err
    Jumlah = jumlah
    ERROR = np.sqrt(np.mean(jumlah))
    return {"ERROR":ERROR, "JUMLAH":Jumlah}








def inversion( iteration, listresis, frek,Dr, pmarlow, pmarhi, depthincremental):


    import numpy as np
    from scipy import optimize
    Resistivities = listresis
    thictemp =  list(range(1,5000,depthincremental))
    thicknesses =  thictemp[:len(Resistivities)]
    ERROR = 200
    iterasi = 0
   
    while iterasi <iteration: 
        
        J = Jacobian(Resistivities,thicknesses,frek)
        Jt = J.transpose()
        
        I = np.identity(len(Resistivities+thicknesses)-1)
    
        
        Gm = forward(Resistivities,thicknesses,frek)["appresphase"]
        d_kurang_g =np.array(Dr).reshape(len(frek)*2,1)-np.array(Gm).reshape(len(frek)*2,1)
        
        
        def Inversi(mardquart): 
          
            
            DM = np.dot(np.dot(np.linalg.inv(np.dot(Jt,J)+np.dot(mardquart,I)),Jt),d_kurang_g)
        
            M = np.array(Resistivities + thicknesses[:-1]).reshape(len(Resistivities+thicknesses)-1,1)
            
            Iter1 = M+DM
            Resistivities1 = [item for sublist in Iter1[:int((len(Iter1)+1)/2)].tolist() for item in sublist]
            thicknesses1 = [item for sublist in Iter1[int((len(Iter1)+1)/2):].tolist() for item in sublist]+[7000]
            
            ERROR = RMS(forward (Resistivities1,thicknesses1,frek)["appresphase"],Dr)["ERROR"]
        
            return ERROR    
            
            
            
            
            
            
            
        mardquart= optimize.golden(Inversi,brack=(pmarlow,pmarhi))
        
        DM = np.dot(np.dot(np.linalg.inv(np.dot(Jt,J)+np.dot(mardquart,I)),Jt),d_kurang_g)
        
        M = np.array(Resistivities + thicknesses[:-1]).reshape(len(Resistivities+thicknesses)-1,1)
        
        Iter1 = M+DM
        
        
        
        Resistivities = [item for sublist in Iter1[:int((len(Iter1)+1)/2)].tolist() for item in sublist]
        thicknesses = [item for sublist in Iter1[int((len(Iter1)+1)/2):].tolist() for item in sublist]+[7000]
        
        ERROR = RMS(forward (Resistivities,thicknesses,frek)["appresphase"],Dr)["ERROR"]
        jumlahrms = RMS(forward (Resistivities,thicknesses,frek)["appresphase"],Dr)["JUMLAH"]
        iterasi = iterasi+1
    returndata = [Resistivities,thicknesses,ERROR,iterasi,jumlahrms]
    return { "returndata" : returndata}



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


def readedi(i,komponen):

    import math
    import cmath
    import numpy as np
    frek=[]
    ZYXR=[]
    ZYXI=[]
    with open(i) as input_data:
        for line in input_data:
            if line.startswith('>FREQ'):
                break
        for line in input_data:
            if line.startswith('>'):
                break
            frek.append(line.split())
        for line in input_data:
            if line.startswith(komponen):
                break
        for line in input_data:
            if line.startswith('>'):
                break 
            ZYXR.append(line.split())
        
        for line in input_data:
            if line.startswith('>'):
                break 
            ZYXI.append(line.split())
    ZYXR = [item for sublist in ZYXR for item in sublist]
    ZYXI = [item for sublist in ZYXI for item in sublist]
    frek = [item for sublist in frek for item in sublist]
    ZYXR = list(map(float, ZYXR))
    ZYXI = list(map(float, ZYXI))
    frek = list(map(float, frek))
    

    
    ZYX = []
    for r,i in zip(ZYXR,ZYXI):
        n = complex(r,i)
        ZYX.append(n)
    appres=[]
    phase = []
    for z,f in zip(ZYX,frek):
        mu = 4*math.pi*1E-7
        w = 2*math.pi*f
        apparentresis = (abs(z)**2)/(w)
        fase = math.atan(z.imag/ z.real)
        fase = (np.rad2deg(fase))
        appres.append(apparentresis)
        phase.append(fase)
    

    # def running_mean(x, N):
    #   """ x == an array of data. N == number of samples per average """
    #   cumsum = numpy.cumsum(numpy.insert(x, 0, 0)) 
    #   return (cumsum[N:] - cumsum[:-N]) / float(N)
    
    # raw_data = appres
    # raw_data2= phase

    # filtappres = running_mean(raw_data, 3)
    # filtphase = running_mean(raw_data2, 3)
    # filtfrek  = frek[:-(len(frek)-len(filtappres))]

    # raw_data = [frek,appres]
    # raw_data2= [frek,phase]
    # n = 4  # the larger n is, the smoother the curve will be
    # n2 = 4
    # b = [1.0 / n] * n
    # b2 = [1.0 / n2] * n2
    # a = 1
    # filtappres = sg.filtfilt(b, a, raw_data)
    # filtphase = sg.filtfilt(b2, a, raw_data2)
    
    # Stasiun = [list(filtappres[1]),list(filtphase[1]),list(filtappres[0]),list(filtappres[0]),list(filtappres[1]),list(filtphase[1])]
    # Stasiun = [filtappres,filtphase,filtfrek,filtfrek,filtappres,filtphase]

    Stasiun = [ appres,phase,frek,frek,appres,phase]
    return Stasiun










filename=[]
def index(request):
    context = {
        'greetings':'Welcome to  Resistivity Sectio tool'
        }
    if request.method == 'POST':
        uploaded_file = request.FILES['document']
       
        fs = FileSystemStorage()
        fs.save(uploaded_file.name, uploaded_file)
        filename.append(uploaded_file.name)
    return render(request,'Resistivity_Section/index.html', context)




def add(request):
    
    komponen = str(request.GET['var1'])
    resis = str(request.GET['var2'])
    resis = list(resis.split(","))
    resis = list(map(int, resis))
    thicinc = int(request.GET['var3'])
    iteration = int(request.GET['var4'])
    pmarlow = int(request.GET['var5'])
    pmarup = int(request.GET['var6'])
    targetdepth= int(request.GET['var7'])

    Data = []


    for files in glob.glob(r"C:\dev\Project\django_projectv2\forward\EDI2\*.*"):

        Stasiun = readedi(files,komponen)
        Data.append(Stasiun)










    ########################################################################################
    ############################################################################################
    def rhophi2rhodepth(rho, phase, frekuensi):
        resis = []
        depth0 = []
        for r,p,f in zip(rho,phase,frekuensi):
            depth = np.sqrt(r*(1/f)/2/np.pi/MTcc.mu0)
            rhonb = r * (np.pi/2/np.deg2rad(p%90) - 1)
            resis.append(rhonb)
            depth0.append(depth)
        depth1 = []
        depth1.append(depth0[0])
        for i in range(len(depth0)-1):
            d = depth0[i+1]-depth0[i]
            depth1.append(d)
        return {"rho_nb":resis, 
              "depth":depth1,
              "depth1":depth0,
              }












    dh = 1.0e-7




    value = []
    depthp = []
    displot = []
    displott=[]
    erracc = []
    jumlahrms =[]

    frekuensipseudo=[]
    apprescalculated=[]

    for i in Data:
    
        Dr = i[0]+i[1]
    
    
    
    

        inversout = inversion(iteration, resis, i[2], Dr,pmarlow,pmarup,thicinc)["returndata"]
    
        Resistivities = inversout[0]
        thicknesses = inversout[1]
        ERROR = round(inversout[2],2)
        iterasi = inversout[3]
        jumlahrms.append(inversout[4])
    
        calculated = forward (Resistivities,thicknesses,i[2])["appresphasemisah"]
        Calappres= calculated[0]
        Calphase = calculated[1]
    
        fig = plt.figure(figsize=(12,8))
        ax1 = plt.subplot(221)
        l1, = ax1.plot(i[2],Calappres,'r')
        l2, = ax1.plot(i[2],i[4],'b'"d")
        ax1.legend((l1,l2),('Calculated', 'Obs'))
        ax1.set_xscale('log')
        # ax1.set_yscale('log')
        ax1.invert_xaxis()
        ax1.set_xlabel("Frrekuensi(Hz)")
        ax1.set_ylabel("Resistivitas Apparent (Ohmm)")
    
       
        ax2 = plt.subplot(223)
        l3, = ax2.plot(i[2],Calphase,'r')
        l4, = ax2.plot(i[3],i[5],'g'"d")
        ax2.legend((l3,l4),('Calculated', 'Obs'))
        ax2.set_xscale('log')
        # ax2.set_yscale('log')
        ax2.invert_xaxis()
        ax2.set_xlabel("Frrekuensi(Hz)")
        ax2.set_ylabel("Fase (derajat)")
    

        sounding = depthplot(thicknesses, Resistivities)
    
    

       
        ax3 = plt.subplot(122)
        ax3.step(sounding[0],sounding[1],'r')
   
    
        ax3.set_xscale('log')
        ax3.set_xlim(0.1,10000)
        ax3.set_ylim(0,targetdepth)
        ax3.set_xlabel('Resistvity (Ohmm)')
        ax3.set_ylabel('Depth (km)')
        plt.gca().invert_yaxis()
        plt.suptitle('Iterasi {} | RMS = {}'.format(iteration,ERROR),
                  size=20)
        plt.grid()
        plt.show()
        value.append(sounding[0])
        depthp.append(sounding[1])
        erracc.append(ERROR)
        frekuensipseudo.append(i[2])
        apprescalculated.append(Calappres)

    RMStot = round( np.sqrt(np.mean(jumlahrms)),2)
  
    depthtemp =[]
    valuetemp=[]


    for o,g in zip(depthp,value) :

        xnt = []
        vnt = []
        for i in range(len(o)-1):
            u = list(np.arange(o[i]+0.00001,o[i+1],1))
            r = list(itertools.repeat(g[i], len(u)))
            if not u:
                xnt.append(np.NAN)
                vnt.append(np.NAN)
            else:
                xnt.append(u)
                vnt.append(r)
        depthtemp.append( [x for x in xnt if str(x) != 'nan'])
        valuetemp.append( [x for x in vnt if str(x) != 'nan']) 
 

    w = 0
    for u in range(len(valuetemp)):
        displot1=[]
        for i in range(len(valuetemp[u])):
            displot2=[]
            for x in range(len(valuetemp[u][i])):
                displot2.append(1+w)
            displot1.append(displot2)
        w = w+1
        displot.append(displot1)
    


    jarak = [1,1,1,1,1,1,1,1,1]
    w = 0
    for u,j in zip(range(len(value)),jarak):
        displot1=[]
        for i in range(len(value[u])):
            displot1.append(w)
        w = w+j
        displott.append(displot1)
    
    valuet = [item for sublist in valuetemp for item in sublist]
    depthpt = [item for sublist in depthtemp for item in sublist]
    displotft = [item for sublist in displot for item in sublist]

    value = [item for sublist in valuet for item in sublist]
    depthp = [item for sublist in depthpt for item in sublist]
    displotf = [item for sublist in displotft for item in sublist]






    stasiunx = []
    for i in displott :
        stasiunx.append(i[0])


    stasiunz= [0]*len(stasiunx)

    average_rms= round(np.mean(erracc),2)

    import scipy.interpolate as interp

    X = np.linspace(0, np.max(displotf)+1, 400)
    Y = np.linspace(0, targetdepth, 400)
    x, y = np.meshgrid(X,Y)





    cov_model = gs.Gaussian(dim=2, len_scale=4, anis=0.2, var=0.5, nugget=0.1)
    OK1 = OrdinaryKriging(displotf, depthp, value, cov_model)
    z1, ss1 = OK1.execute("grid", X, Y)



    color_map = plt.cm.get_cmap('jet')
    reversed_color_map = color_map.reversed()


    cmap = mpl.colors.ListedColormap(["red","orangered","yellow","greenyellow","lime","limegreen","forestgren","royalblue","blue","mediumblue"])

    # ['blue','blue', 'aquamarine','springgreen','lime','greenyellow', 'yellow','orange','orangered','deeppink']
    cmap.set_over('mediumblue')
    cmap.set_under('maroon')


    bounds = np.linspace(0, 150,11)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    fig1 = plt.figure(0)
    plt.scatter(stasiunx,stasiunz,s = 200, marker="v", color = "black")
    plt.contour( z1,extent=(0,np.max(displott),0,targetdepth),colors= "black",levels = 7, alpha=0.2)
    plt.imshow(z1,extent=[0,np.max(displott),0,targetdepth],origin="lower",cmap=reversed_color_map , aspect = np.max(displott)/10)

    plt.gca().invert_yaxis()

    plt.xlabel("Distence (Km)")
    plt.ylabel("Depth (Km)")
    plt.title("Resistivity Section")
    plt.suptitle('Iterasi {} | meanRMS = {} | RMS = {}'.format(iteration,average_rms,RMStot),size=8,x=0.4,y=0.15)
    plt.clim(0,150)
    plt.colorbar(cmap=reversed_color_map , norm=norm,
         ticks=bounds,boundaries=bounds,format='%1i').set_label("Ohmm")
    plt.show()





    

























    buffer = BytesIO()
    fig1.savefig(buffer, format='png')
    buffer.seek(0)
    image_png = buffer.getvalue()
    graph = base64.b64encode(image_png)
    graph = graph.decode('utf-8')
    buffer.close()

    out1={'graph':graph}
    
    return render(request,'Resistivity_Section/output.html', out1)


