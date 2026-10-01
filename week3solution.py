import numpy as np

class fftsig:

    def __init__(self):
        self.f=np.array
        self.FT=np.array
        self.amp=np.array
        self.phase=np.array
        

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
# mysig.f -- the frequency vectory, shifted so that f=0 is in the centre of the array
# mysig.FT -- the FFT of the signal shifted so that f=0 is at the centre of the array
# mysig.amp -- the amplitude spectrum 
# mysig.phase -- the phase spectrum
    nt=np.size(t)
    dt=t[2]-t[1]
    class mysign:
        def __init__(self,s1):
            self.FT = np.fft.fftshift(np.fft.fft(s1))
            self.amp = np.fft.fftshift(np.abs(np.fft.fft(s1)))
            self.phase = np.angle(np.fft.fftshift(np.fft.fft(s1)))
            self.f = np.fft.fftshift(np.fft.fftfreq(nt,dt))
            
    x=mysign(s1)
    
    return x
    


    
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
