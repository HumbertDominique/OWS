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

def psd(N, r0, dx, l0 = 0, L0 = -1):
    """
    Generates a Phase power spectrum (PSD).
    Parameters:
    N (int): Size of the phase screen matrix (must be even?).
    r0 (float): Fried's parameter
    dx (float): Pixel size [m].
    lo (float): Inner scale of turbulence.
    L0 (float): Outer scale of turbulence (set to -1 for infinite scale).

    Returns:
    PSD (ndarray): 2D array representing the generated phase screen.
    """    
    # 2025.03.14 - Original version based on paola.pro (22.12.2022) by Laurent Jolissaint (HES-SO), lines 3659+
    # Not sure where all the coefficients come from
    #TODO: Change inherited variable names to something more meaningful
    #TODO: decide weather to implement l0
    Df = 1/(N*dx)
    
    fx = np.linspace(-N/2,N/2-1,N)* Df
    fpcoohfx, fpcoohfy = np.meshgrid(fx,fx) # Since it will be computed anyways, fpcoohfy will be used instead of rotating the xmatrix.
    # Prepare the radius matrix
    xpcoohf = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            xpcoohf[i][j]= np.sqrt(fpcoohfx[i][j]**2 + fpcoohfy[i][j]**2)  # pupil plane coordianate radius

    # Von Karman corrected PSD (L0 != -1) PAOLA:3666
    if L0 == -1:
        # TODO: Evaluate wheather it is necessary to check weather (px.val <= 0)
        # dphi = 6.883877*(xpcoohf[w]/r0)**(5/3)
        Watm = np.zeros((N,N), dtype=np.double)
        for i in range(N):
            for j in range(N):
                if fpcoohfx[i][j]**2 + fpcoohfy[i][j]**2 != 0:
                    Watm[i][j] = 0.022896/(r0**(5/3))*(fpcoohfx[i][j]**2 + fpcoohfy[i][j]**2)**(-11/6)
    else:
        #w = np.where(fpcoohfx != 0) #TODO: implement actual logic if required
        # TODO: Evaluate wheather it is necessary to check weather (px.val <= 0)
        #dphi = 0.171661*(L0/r0)**(5/3)*((1.005635)-(2*np.pi*xpcoohf[w]/L0))**(5/6)
                               #beselk(5/6,2*np.pi*xpcoohf/L0)
        Watm = 0.022896/(r0**(5/3))*(fpcoohfx**2+fpcoohfy**2 +1/(L0**2))**(-11/3)

    PSD = Watm#*dphi
    PSD[int(N/2), int(N/2)] = 0
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