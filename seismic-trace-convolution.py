#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


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

def ricker_wavelet(f, dt, tmax=0.05):
    #input: central frequency f, time step dt, maximum time tmax
    #output: ricker wavelet
    t = np.arange(-tmax, tmax + dt, dt)
    return (1 - 2 * (np.pi * f * t)**2) * np.exp(-(np.pi * f * t)**2)

wavelet = ricker_wavelet(f=30, dt=dt, tmax=0.05)

#%%synthetic seismogram

#Parametros del modelo y geometria
V = 2000            # Velocidad de la capa (2000 m/s)
x1 = 5            # Distancia a la fuente del Receptor 1 (100 m)
#x2 = 400            # Distancia a la fuente del Receptor 2 (400 m)

#reflectividad
r=np.zeros([1000])
depth_reflector=300
#r[800]=.2
#r[550]=-.4


# Calcular tiempos de llegada (t0, t1, t2)
nt = np.size(r)          # numero de muestras
t_max = 1.5         # Tiempo máximo de registro
dt = t_max / (nt - 1)
t = np.linspace(0, t_max, nt)

t0 = (2 * depth_reflector) / V          # tiempo de llegada vertical al reflector
t1 = np.sqrt(t0**2 + (x1**2 / V**2))    # tiempo de llegada al receptor 1

# Encontrar el índice de tiempo de la llegada más cercano a t1
idx1 = np.argmin(np.abs(t - t1))
r[idx1] = 1 #Representando la reflectividad del reflector en el tiempo de llegada al receptor 1

#Fuente de Ricker(30 Hz)
wavelet = ricker_wavelet(30,dt)

c1=np.convolve(r,wavelet,mode='same')


# Visualización de trazas simple: x = offset, y = tiempo
plt.figure(1, figsize=(10, 5))

source_x = 0
receiver_x = x1
trace_x = receiver_x +   c1

plt.scatter(source_x, -0.02, marker='*', s=200, color='red', zorder=5, label='Source')
plt.plot(trace_x, t, color='black', linewidth=1.2)
plt.scatter(receiver_x, -0.02, marker='v', s=100, color='blue', zorder=4, label='Receiver')

plt.xlim(-20, receiver_x + 10)
plt.ylim(t[-1], -0.05)
plt.xlabel('Offset (m)')
plt.ylabel('Time (s)')
plt.title('Sythetic traces')
plt.grid(True, alpha=0.3)
plt.legend(loc='upper right')
plt.show()




# %%
