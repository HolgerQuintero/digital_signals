# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 10:51:07 2026

@author: Holger Quintero
"""
import numpy as np
import matplotlib.pyplot as plt


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
    proj=np.dot(a,b)
    return proj
fs=600
sig1,t=mysin(30,600)
sig2,t2=mysin(60,600)
sig3=sig1+0.1*sig2

plt.figure(1)
plt.plot(t,sig1,label='30 Hz')
plt.plot(t2,sig2,label='60 Hz')
plt.plot(t,sig3,label='S3=S1+0.1S2')
plt.legend()
plt.xlabel('time (s)')
plt.ylabel('amplitude')
plt.title('Input Sine Waves')


projS3_S1=bases_proj(sig1, sig3)
projS3_S2=bases_proj(sig2, sig3)

sig4=2*sig1-0.5*sig2


projS4_S1=bases_proj(sig1, sig4)
projS4_S2=bases_proj(sig2, sig4)


def myfft(t,s1):
#mfft(t,s1) computes the Fourier Transform, 
#the phase spectrum and the amplitude spectrum as well as the frequency sampling vector of an input signal
#outputs: 
# mysig.f -- the frequency vectory, shifted so that f=0 is in the centre of the array
# mysig.FT -- the FFT of the signal shifted so that f=0 is at the centre of the array
# mysig.amp -- the amplitude spectrum 
# mysig.phase -- the phase spectrum
    nt=np.size(t)
    dt=t[2]-t[1]
    class mysign:
        def __init__(self,s1):
            self.FT = np.fft.fftshift(np.fft.fft(s1))
            self.amp = np.abs(self.FT)
            self.phase = np.angle(self.FT)
            self.f = np.fft.fftshift(np.fft.fftfreq(nt,dt))
            
    x=mysign(s1)
    
    return x

mysign=myfft(t,sig1)




dt=t[2]-t[1]
df=mysign.f[2]-mysign.f[1]
nt=np.size(t)
nf=np.size(mysign.f)
tmax=np.max(t)
fmax=np.max(mysign.f)



sig_cos=np.cos(2*np.pi*30*t)
mysign2=myfft(t,sig_cos)

plt.figure(2)
plt.plot(mysign.f, mysign.amp, label='sin')
plt.plot(mysign2.f, mysign2.amp, '--', label='cos')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.legend()
plt.title('Amplitude spectra')

plt.figure(3)
plt.plot(mysign.f, mysign.phase, label='sin')
plt.plot(mysign2.f, mysign2.phase, label='cos')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase (rad)')
plt.legend()
plt.title('Phase spectra')

#Ploteamos ahora cada parte de la transformada
plt.figure(4)
plt.plot(mysign.f, np.real(mysign.FT), label='Real part')
plt.plot(mysign.f, np.imag(mysign.FT), label='Imaginary part')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('Real and Imaginary Parts of the Fourier Transform')
plt.legend()
plt.grid()


nyquist_f=1/(2*dt)
f=np.arange(-nyquist_f,nyquist_f,df)

def myslowdft(f,samplerate,mysig):
#myslowdft(f,samplerate,mysig)
#f frequency vector at which to compute the coefficients
# samplerate, the sample rate of the signal mysig
# mysig the signal to take the dft of
# returns the imaginary part of the Fourier transform of mysig (unscaled)
    ii=0
    coeff=np.zeros_like(f)
    t=np.arange(np.size(mysig))
    for freq in f:
        coeff[ii]=1*np.sum(mysig*np.sin((2*np.pi*freq*t)/samplerate))
        ii=ii+1

        
    return coeff

slowft=myslowdft(f,fs,sig1)
plt.figure(5)
plt.plot(f, slowft)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Imaginary amplitude')
plt.title('Manual imaginary part of the DFT')
plt.grid()
plt.show()


#%%Speed comparison

import time

fs_test = 500
sig_test, t_test = mysin(30, fs_test, 1.0)

N = np.size(sig_test)
df_test = fs_test/N
f_test = np.arange(-fs_test/2, fs_test/2, df_test)

# Manual DFT
start = time.time()
slowft_test = myslowdft(f_test, fs_test, sig_test)
end = time.time()

time_slow = end - start


# NumPy FFT
start = time.time()
fastft_test = np.fft.fft(sig_test)
end = time.time()

time_fast = end - start


print('Number of samples:', N)
print('Manual DFT:', time_slow, 'seconds')
print('NumPy FFT:', time_fast, 'seconds')
print('FFT is approximately', time_slow/time_fast, 'times faster')