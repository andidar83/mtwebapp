import numpy
cimport numpy
def forward (res,thic, frekuensi):
    import math
    import cmath
    import numpy as np
    
    cdef float mu, w
    cdef complex i
    cdef list apparentresistivity, phase
    cdef double r, h

    

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
    cdef float dh
    cdef list JR,JT
    cdef int i

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
    cdef float jumlah, ERROR
    cdef float h, g
    import numpy as np
    jumlah=0
    for h,g in zip(calculated,obs):
        err = (h - g)**2#(((g-statistics.mean(D))**2)/len(D))**0.5
        jumlah = jumlah + err
    Jumlah = jumlah
    ERROR = np.sqrt(np.mean(jumlah))
    return {"ERROR":ERROR, "JUMLAH":Jumlah}








def inversion( iteration, listresis, frek,Dr, pmarlow, pmarhi, depthincremental):
    # Resistivities = rhophi2rhodepth(i[0], i[1], i[2])["rho_nb"]
    # thicknesses = rhophi2rhodepth(i[0], i[1], i[2])["depth"]
    cdef float ERROR
    cdef int iterasi
    cdef list Resistivities,thicknesses,Gm, Resistivities1, thicknesses1
    cdef numpy.ndarray J , Jt, I, DM, M, Iter1, d_kurang_g
 

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
            cdef numpy.ndarray DM, M, Iter1
            
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
    cdef list depth
    cdef list Resistivitiess
    cdef double elev
    cdef double i
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


def readedi(i):
    cdef list frek, ZYXR, ZYXI, appres, phase, Stasiun
    cdef double mu, w
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
            if line.startswith('>ZYXR'):
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
    Stasiun = [ appres,phase,frek,frek,appres,phase]
    return Stasiun
