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

def phase_screen(N, PSD, dxp, SEED = None):
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
    # 2025.03.14 - Original version based on paola.pro (22.12.2022) by Laurent Jolissaint (HES-SO), lines 3659+
    # Not sure where all the coefficients come from
    #TODO: Change inherited variable names to something more meaningful

    tmp2 = PSD*(N*dxp)**2

    #TODO: understand what this is here
    tmp2[0 : N, N-1] = tmp2[0 : N, 0] # last row = 1st row
    tmp2[N-1, 0 : N] = tmp2[0, 0 : N] # last column = 1st column
    tmp2[N-1, N-1] = tmp2[0, 0] # bottom right corner = top left corner

    if SEED is not None:
      np.random.seed(SEED)
    tmp2 = np.sqrt(tmp2)* np.random.rand(N, N)
    tmp1 = np.sqrt(2)*(tmp2+np.rot90(tmp2,2))/2
    tmp2 = np.random.rand(N,N)*2*np.pi
    # in order to make sure that the phase is real,
    # (1) the real part of the phase spectrum is forced to be even
    # (2) the imaginary part of phase spectrum is forced to be odd
    tmp2 =  (tmp1 * np.cos(0.5 * (tmp2 - np.rot90(tmp2)))) +1j*(tmp1 * np.sin(0.5 * (tmp2 - np.rot90(tmp2))))
    phaseft = tmp2

    #phaseft -= np.mean(phaseft) # remove piston?
    return np.fft.ifft2(np.fft.ifftshift(phaseft))