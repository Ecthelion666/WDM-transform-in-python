import scipy.special
import scipy
from scipy.fftpack import fft,ifft
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
from pycbc.waveform import get_td_waveform

dt=1/1024 #delta_t

apx = 'IMRPhenomD'
hp, hc = get_td_waveform(approximant=apx,
                        mass1=10,
                        mass2=10,
                        delta_t=dt,
                        f_lower=10,distance=10)
Nt=len(hp)

nx=4# filter steepness in frequency
Nf=128#frequency layers,also M in the paper
#V Necula et al 2012 J. Phys.: Conf. Ser. 363 012032
t_step=Nf*dt
mult=16#oversampling

OM=math.pi/dt#Omega
DOM=OM/Nf
insDOM= 1/math.sqrt(DOM)
B=DOM/2
A=(DOM-B)/2
K=mult*Nf*2#samples of phitilde
dom=2*math.pi/(dt*mult*2*Nf)


def phitilde(om,insDOM,A,B):
    if abs(om)<A:
        return insDOM
    elif A<=abs(om)<A+B:
        return insDOM*math.cos(math.pi/2*scipy.special.betainc(nx, nx, (abs(om)-A)/B))
    else:
        return 0

def Cmn(m,n):
    if (n+m)%2==0:
        return 1
    else:
        return i

#construct phitilde and transform to time domain
phif=[]
#zero frequency
phif.append(insDOM)
#postive frequencies
for i in range(1,int(K/2+1)):
    phif.append(phitilde(dom*i,insDOM,A,B))
#negative frequencies
for i in range(1,int(K/2+1)):
    phif.append(phitilde(dom*(K/2-i+1),insDOM,A,B))

phit=ifft(phif)#time step is dt/mult

#fastWDMtransform
X=[]
for i in range(int(Nt/Nf)-1):
    Xn=[]
    for j in range(2*Nf):
        Xn.append(hp[i*Nf+j-Nf]*phit[(j-Nf)])
    fftXn=fft(Xn)
    wdm=[]
    count=0
    for j in fftXn:
        if count>0 and count<Nf:
            wnm=j*Cmn(count,i)
            wdm.append(abs(wnm))
        count+=1
    X.append(wdm)

lenx=len(X)
leny=len(X[0])
X=np.matrix(X)
X=X.T

#plot
x=np.linspace(0,hp.duration,lenx)
y=np.linspace(0,OM/(2*math.pi),leny)
plt.figure()
norm=matplotlib.colors.LogNorm(vmin=10**(-19))
plt.pcolormesh(x,y,X,norm=norm)
plt.xlabel('time(s)')
plt.ylabel('freq(Hz)')
plt.colorbar()
plt.show()



