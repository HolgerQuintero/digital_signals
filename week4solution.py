#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

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
    c=np.zeros([na+nb-1])
    nc=np.size(c)
    
    for ic in range(nc):
        for isum in range(max(na,nb)):
            if np.abs(isum) <na and np.abs(ic-isum)<nb and ic-isum>=0:
               suma=a[isum] * b[ic - isum]
               c[ic]=c[ic]+suma

    return c




def mycorr(f,g):
# return the correlation of the signals f and g
    nf=np.size(f)
    ng=np.size(g)
    h=np.zeros([nf+ng-1])
    #nh=np.size(h)
    index = 0    # index is 0 to len(h) to calculate h
    for ih in range(-(nf-1),ng): # ih is lag index
        for ind in range(max(nf,ng)):
            if np.abs(ih+ind) <ng and ih+ind>-1 and ind <nf  : # if 1) g index is less than len(g)                                                              # if 2) g index is positive                                                   # if 3) f index is less than len(f)
                val = f[ind]*g[ih+ind]
                h[index]=h[index]+val
        index = index + 1
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

def ricker_equivalente_scipy(points, a):
    """Genera exactamente el mismo vector que scipy.signal.ricker"""
    A = (points - 1) / 2.0
    x = np.arange(-A, points - A)
    
    vec = 1 - (x / a) ** 2
    ts = np.exp(-0.5 * (x / a) ** 2)
    
    # Factor de normalización estricto de SciPy
    return 2 / (np.sqrt(3 * a) * np.pi**0.25) * vec * ts

#%%Convolution

a=np.array([1,-2,3,0])
na=np.size(a)
b=np.array([0,2,4])
nb=np.size(b)
c=np.zeros(na+nb-1)  

c[0]=a[0]*b[0]
c[1]=a[0]*b[1]+a[1]*b[0]
c[2]=a[0]*b[2]+a[1]*b[1]+a[2]*b[0]
c[3]=          a[1]*b[2]+a[2]*b[1]+a[3]*b[0]
c[4]=                    a[2]*b[2]+a[3]*b[1]
c[5]=                              a[3]*b[2]

print(c)


c=myconv(a,b)
print(c,np.convolve(a,b))

#%% sine 30Hz
(y,t)=mysin(30,600)


#flitered sine
filter=np.array([0,-1,-1,0])
y_filt=myconv(y,filter)
nb=np.size(filter)
dt=t[2]-t[1]
nt=np.size(t)
t2=np.arange(-int(nb/2)+1,nt+int(nb/2))*dt

plt.figure(1)
plt.plot(t2,y_filt,label='filtered signal')
plt.plot(t,y,label='30Hz sine')
plt.xlabel('time')
plt.ylabel('amplitude')
plt.legend()

example=np.arange(5)
myconv(example,filter)

# %% Set of signals
data = np.loadtxt(r"C:\Users\geral\OneDrive\Documentos\mun\repository\digital_signals\Sigs_for_shift.txt")

t = data[0]
s1 = data[1]
s2 = data[2]

plt.figure(2)
plt.plot(t, s1, label="Signal 1")
plt.plot(t, s2, label="Signal 2")
plt.xlabel("Time ($\mu$s)")
plt.ylabel("Amplitude")
plt.title("Signals")
plt.legend()

conv_my = myconv(s1, s2)
conv_full = np.convolve(s1, s2, mode='full')
conv_same = np.convolve(s1, s2, mode='same')
conv_valid = np.convolve(s1, s2, mode='valid')

plt.figure(3)
plt.subplot(2, 2, 1)
plt.plot(conv_my)
plt.title("My convolution")
plt.xlabel("Time")
plt.ylabel("Amplitude")

plt.subplot(2, 2, 2)
plt.plot(conv_full)
plt.title("Numpy convolution (full)")
plt.xlabel("Time")
plt.ylabel("Amplitude")

plt.subplot(2, 2, 3)
plt.plot(conv_same)
plt.title("Numpy convolution (same)")
plt.xlabel("Time")
plt.ylabel("Amplitude")

plt.subplot(2, 2, 4)
plt.plot(conv_valid)
plt.title("Numpy convolution (valid)")
plt.xlabel("Time")
plt.ylabel("Amplitude")

plt.tight_layout()
plt.show()



#%%Cross-correlation

print(mycorr(a,b),mycorr(b,a))#Same vector backwards

corr_s1_s2=mycorr(s1, s2)

nt=np.size(t)
dt=t[2]-t[1]
t_corr=np.arange(-nt+1,nt)*dt

plt.figure(4)
plt.plot(t_corr, corr_s1_s2)
plt.xlabel("Time ($\mu$s)")
plt.ylabel("Coefficient")
plt.title("Cross-correlation")


maxloc=np.argmax(np.abs(corr_s1_s2))
shift=t_corr[maxloc]
print(shift, "difference in arrival time")

#If we switch signals
t_corr2=np.arange(-nt+1,nt)*dt
shift=t_corr[np.argmax(np.abs(mycorr(s2, s1)))]
print(shift, "difference in arrival time (switch signals)")


#%%synthetic seismogram

r=np.zeros([1000])
r[300]=1
r[800]=.2
r[550]=-.4
wavelet = ricker_equivalente_scipy(30, 3)
#plt.figure(5)
#plt.plot(r)
#plt.plot(wavelet)

c1=myconv(r,wavelet)
c2=mycorr(r,wavelet)

plt.figure(6)
#plt.plot(c1,label=r'$r\ast wl$')
plt.plot(c2,label=r'$r\star wl$')
plt.legend()
plt.xlabel('index number')
plt.ylabel('amplitude')
plt.title('synthetic seismogram')

plt.figure(7)
plt.figure(figsize=(9,4))
plt.plot(c1)
plt.xlabel('Index number')
plt.ylabel('Amplitude')
plt.title('Convolution: Reflectivity * Wavelet')
plt.grid()
plt.show()


# Cross-correlation
plt.figure(figsize=(9,4))
plt.plot(c2)
plt.xlabel('Index number')
plt.ylabel('Amplitude')
plt.title('Cross-correlation: Reflectivity and Wavelet')
plt.grid()
plt.show()