import numpy as np


class fftsig:

    def __init__(self):
        self.f=np.array
        self.FT=np.array
        self.amp=np.array
        self.phase=np.array

def myconv(a,b):
# return the convolution of the signals a and b

    na=np.size(a)
    nb=np.size(b)
#    c=
    nc=np.size(c)
    
    for ic in range(nc):
        for isum in range(max(na,nb)):
            if np.abs(isum) <na and np.abs(ic-isum)<nb and ic-isum>=0:
#                c[ic]=

                 
    return c

def mycorr(f,g):
# return the correlation of the signals f and g
    nf=np.size(f)
    ng=np.size(g)
    h=np.zeros([nf+ng-1])
    nh=np.size(h)
    index = 0    # index is 0 to len(h) to calculate h
    for ih in range(-(nf-1),ng): # ih is lag index
        for ind in range(max(nf,ng)):
            #You need to do 3 things in here:
            #  1. Check that all inidices will be within the associated array length (note that there are 3 checks to do)
            #  2. Calculate the cross-correlation (using the right indices) and add it into your array
            #  3. Increment the 'index'
            #It is probably easiest to do these steps in the order 2, 3 then 1
#            if ???
                #h[index]=
         index=index+1
        
    return h



def mysin(f,fs,tmax=0.1,tmin=0.0):
    #mysin(f,fs,tmax=0.1):
    # f is the frequency of the signal
    # tmax is the final time
    # fs is the sampling frequency, fs=1/dt
    ##HINT: You need to make a time vector with the appropriate dt, so first compute that then dt, then t, then y
    
    dt=1.0/fs
    t=np.arange(tmin,tmax,dt)
    y=np.sin(2*np.pi*f*t)

    return y,t


    
def bases_proj(a,b):
    #bases_proj(a,b) projects vector a onto vector b (or vice-versa provided all signals are real)

    return np.dot(a,b)

def myfft(t,s1):
#mfft(t,s1) computes the Fourier Transform, 
#the phase spectrum and the amplitude spectrum as well as the frequency sampling vector of an input signal
#outputs: 
# f -- the frequency vectory, shifted so that f=0 is in the centre of the array
# s1f -- the FFT of the signal shifted so that f=0 is at the centre of the array
# s1f_amp -- the amplitude spectrum 
# s1f_phase -- the phase spectrum
    nt=np.size(t)
    dt=t[2]-t[1]
    mysig=fftsig()
    s1f=np.fft.fftshift(np.fft.fft(s1))
    mysig.FT=s1f
    s1f_amp=np.abs(s1f)
    mysig.amp=s1f_amp
    s1f_phase=np.arctan2(np.imag(s1f),np.real(s1f))
    mysig.phase=s1f_phase
    f=np.fft.fftshift(np.fft.fftfreq(nt,dt))
    mysig.f=f
    
    
    return mysig

#%%Convolution
import numpy as np

a=np.array([1,-2,3,0])
na=np.size(a)
b=np.array([0,2,4])
nb=np.size(b)
c=np.zeros(na+nb-1)  


for i in np.arange(np.size(c)):
    c[i]=0
    for index in np.arange(na):
        if 0 <= i - index < nb:
            sum=a[index] * b[i - index]
            c[i]=c[i]+sum
    
       
