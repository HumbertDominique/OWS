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
    
    # Input checks
    if not isinstance(dimmat, int):
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

def phase_screen(dimmat, PSD, dxp, SEED = None, PUPIL = True):
    """
    Generates a phase screen base on a Phase Spectrum Density (PSD)
    Parameters:
        N (int): Size of the phase screen matrix. must be even?
        PSD (ndarray): phase spectrum density
        dxp (float): Pupil plane pixel size dxp = WL/(NA*Npx
        SEED (int): Seed for the random number generator. (default None)

    Returns:
        PSD (ndarray): complex 2D array representing the generated phase screen.
    """

    df=1.0/(dimmat - dimmat % 2)/dxp

    PP = np.zeros((dimmat+1, dimmat+1)) # Phase power [dimmat+1, dimmat+1] in order to have a pixcele in the middle [rad^2/m^-2]
   
    PP[0:dimmat,0:dimmat] = PSD*(dimmat*dxp)**2
    #TODO: understand what this is here
    PP[0 : dimmat, dimmat] = PP[0 : dimmat, 0] # last row = 1st row
    PP[dimmat,0 : dimmat] = PP[0, 0 : dimmat] # last column = 1st column
    PP[dimmat, dimmat] = PP[0, 0] # bottom right corner = top left corner

    if SEED is not None:
      np.random.seed(SEED)
    PSA = np.sqrt(PP)* np.random.rand(dimmat+1,dimmat+1) # random draw to create a random phase spectrum [rad/m^(-1)] AMPLITUDE
    PSA = np.sqrt(2)*(PSA+np.rot90(np.rot90(PSA,2)))/2 # forced even amplitude of wf spectrum
    PSop = np.random.rand(dimmat+1,dimmat+1)*2*np.pi # random draw for the phase of the spectrum of the optical phase
    # in order to make sure that the phase is real,
    # (1) the real part of the phase spectrum is forced to be even (as in "f(x) is even")
    # (2) the imaginary part of phase spectrum is forced to be odd (as in "f(x) is odd")
    phaseft =  (PSA * np.cos(0.5 * (PSop - np.rot90(np.rot90(PSop)))))[0:dimmat,0:dimmat] + 1j*(PSA * np.sin(0.5 * (PSop - np.rot90(np.rot90(PSop)))))[0:dimmat,0:dimmat]

    MASK = np.ones((dimmat,dimmat))
    if PUPIL is True:
      pxR = np.linspace(-1 ,1,dimmat)*dimmat//2
      xx, yy = np.meshgrid(pxR,pxR)
      R= np.sqrt(xx**2 + yy**2)  # pupil plane spatial frequency radius
      # Créer le masque pour exclure les pixels à l'intérieur du cercle
      MASK[R > dimmat/2] = 0
  
    phase_screen = ((np.fft.ifftshift(np.fft.ifft2(phaseft))*(dimmat*df)**2))*MASK
    phase_screen[MASK == 1] -= np.mean(phase_screen[MASK == 1])

    Sp = phase_screen.size*dxp**2
    apsf = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(MASK*np.exp(-1j*phase_screen))))/Sp

    psf = np.abs(apsf)**2
    return phase_screen, psf
