def function_name(parameters):
   """
   Brief description of the function.
   Parameters:
   parameters (type): Description of the parameters.
   Returns:
   return_type: Description of the return value.
   """
   # TODO: Implement the function logic here
   return

import numpy as np
from ows import fouriertransform as mathft
from matplotlib import pyplot as plt


def psd(dimmat, dxp, r0, wl, L0 = -1):
   """
   Generates a Phase power spectrum (PSD).
   Parameters:
   dimmat (int): Size of the phase screen matrix (should must be even).
   r0 (float): Fried's parameter [m]
   dxp (float): Pupil plane pixel size [m/rad].
   L0 (float): Outer scale of turbulence [m] (set to -1 for infinite scale).
   Returns:
   PSD (ndarray): 2D array representing the generated phase screen.
   """    
   # 2025.03.14 - Original version based on paola.pro (22.12.2022) by Laurent Jolissaint (HES-SO), lines 3659+
    
   # ------------ Input checks ------------
   # ------------ Input checks ------------#
   # ------------ Input checks ------------
   if type(dimmat) != int:
       raise TypeError(f"Expected {int} for dimmat, got {type(dimmat)}")
   if dimmat % 2 != 0:
       raise ValueError("dimmat must be an even integer")
   r0l = r0*(2*wl)**1.2
   Df = 1/(2*dxp)
   #fx = np.linspace(-dimmat//2,dimmat//2,dimmat)* Df
   fx = np.linspace(-1,1,dimmat)*Df*dimmat//2
   freqX, freqY = np.meshgrid(fx,fx)
   freqR= np.sqrt(freqX**2 + freqY**2)  # pupil plane spatial frequency radius
   w = freqR.nonzero()
   PSD = np.zeros_like(freqR, dtype=np.float64)  # Ensure PSD is the same shape and dtype as freqRtype
   PSD[w] = 0.0229*r0l**(-5.0/3) *(freqR[w]**2 + (L0 != -1) / L0**2)**(-11.0/6) # Spatial power spectrum with outer scale
   return PSD

def phase_screen(PSD, dxp, SEED = None, PSF=True, PUPIL = True, PD = [0,0]):
   """
   Generates a phase screen base on a Phase Spectrum Density (PSD)
   Parameters:
      dimmat (int): Size of the phase screen matrix. must be even?
      PSD (ndarray): phase spectrum density
      dxp (float): Pupil plane pixel size dxp = WL/(NA*Npx)
      SEED (int): Seed for the random number generator. (default None)
      PSF (bool): option to return the PSF associated with the phase screen
      PUPIL (bool): True will crop it to a pupil of diameter dimmat if PD = [0,0]
      PD (ndarray of int): Pupil diameters in pixels. PD[outer diameter, inner diameter)

   Returns (list):
      phase_screen (ndarray): real 2D array representing the generated phase screen.
      psf (ndarray): real 2D array representing the PSF associated with the phase screen

   TODO: Add a dimmat argument to allow more pixels on the phase screen and PSF than PSD.

   Based one WaveSeeingLimited.pro, Laurent Jolissaint, March 24, 2025 
   """
   # ------------ Input checks ------------
   # ------------ Input checks ------------
   # ------------ Input checks ------------
   if PSD.shape[0] != PSD.shape[1]:
      raise ValueError("PSD is not square")
   if dxp <= 0:
      raise ValueError("dxp must be positive")
   if type(PD[0]) != int:
      raise TypeError("Pupil diameters must be integers")
   if type(PD[1]) != int:
      raise TypeError("Pupil diameters must be integers")
   dimmat = PSD.shape[0]
   if PD[0]!= 0:
      if PD[0] >= dimmat:
        raise ValueError("Outer diameter must equal or smaller than the number of width of the phase screen.")
      if PD[0] < 0:
        raise ValueError(f"Outer diameter must be positive, got: {PD[0]}")
      if PD[1] < 0:
        raise ValueError(f"Outer diameter must be positive, got: {PD[1]}")
      if PD[0] <= PD[1]:
        raise ValueError(f"Outer diameter must be larger than inner diameter, got: OD={PD[0]}, ID={PD[1]}")

   # ------------ Phase screen ------------
   # ------------ Phase screen ------------
   # ------------ Phase screen ------------
   returns = [[None],[None]]
   
   rng = np.random.default_rng(seed=SEED)

   pxR = np.linspace(-1 ,1,dimmat)*dimmat//2
   xx, yy = np.meshgrid(pxR,pxR)
   R= np.sqrt(xx**2 + yy**2)  # pupil plane spatial frequency radius
   PP = np.zeros((dimmat+1, dimmat+1)) # Phase power [dimmat+1, dimmat+1] in order to have a pixcele in the middle [rad^2/m^-2]
   PP[0:dimmat,0:dimmat] = PSD*(dimmat*dxp)**2
   PP[0 : dimmat, dimmat] = PP[0 : dimmat, 0] # last row = 1st row
   PP[dimmat,0 : dimmat] = PP[0, 0 : dimmat] # last column = 1st column
   PP[dimmat, dimmat] = PP[0, 0] # bottom right corner = top left corner

   PSA = np.sqrt(PP)* rng.normal(0, 1, size=(dimmat+1, dimmat+1)) # random draw to create a random phase spectrum [rad/m^(-1)] AMPLITUDE
   PSA = np.sqrt(2)*(PSA+np.rot90(np.rot90(PSA,2)))/2 # forced even amplitude of wf spectrum
   PSop = rng.normal(0, 1, size=(dimmat+1, dimmat+1))*2*np.pi # random draw for the phase of the spectrum of the optical phase
   # in order to make sure that the phase is real,
   # (1) the real part of the phase spectrum is forced to be even (as in "f(x) is even")
   # (2) the imaginary part of phase spectrum is forced to be odd (as in "f(x) is odd")
   phaseft =  (PSA * np.cos(0.5 * (PSop - np.rot90(np.rot90(PSop)))))[0:dimmat,0:dimmat] + 1j*(PSA * np.sin(0.5 * (PSop - np.rot90(np.rot90(PSop)))))[0:dimmat,0:dimmat]
   
   # ------------ PUPIL OPTION ------------
   # ------------ PUPIL OPTION ------------
   # ------------ PUPIL OPTION ------------
   if PUPIL is True:
      MASK = np.ones((dimmat,dimmat))
      if PD[0] == 0:
        MASK[R >= dimmat//2 +1] = 0
      else:
         MASK[R >= PD[0]//2 +1] = 0
         MASK[R < PD[1]//2] = 0
   phase_screen = MASK*mathft.ift2(phaseft,dxp).real
   phase_screen[MASK == 1] -= np.mean(phase_screen[MASK == 1])
   returns[0] = phase_screen

   # ------------ PSF OPTION ------------
   # ------------ PSF OPTION ------------
   # ------------ PSF OPTION ------------
   if PSF == True:
      if dimmat <= 2**5: 
         N_pad = (2**6)
      elif (2**5 < dimmat <= 2**6): 
         N_pad = (2**8)
      elif (2**6 < dimmat <= 2**7): 
         N_pad = (2**10)
      elif (2**7 < dimmat <= 2**8): 
         N_pad = (2**12)
      elif (2**8 < dimmat <= 2**9): 
         N_pad = (2**14)
      padded_pupil = np.zeros((N_pad,N_pad))
      padded_pupil[:dimmat,:dimmat] = MASK*phase_screen
      Sp = np.sum(R)*dxp**2
      psf = mathft.ft2(padded_pupil,delta=1./dimmat)/Sp
      psf = psf[N_pad//2-dimmat//2:N_pad//2+dimmat//2,N_pad//2-dimmat//2:N_pad//2+dimmat//2]
      returns[1] = psf
    
   return returns