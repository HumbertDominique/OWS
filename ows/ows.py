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

def PSD(dim, r0, L0 = -1):
    """
    Generates a Phase power spectrum (PSD).

    Parameters:
    dim (int): Size of the phase screen matrix (must be even).
    dx (float): Pixel size in meters.
    r0 (float): Fried's parameter at 500 nm.
    L0 (float): Outer scale of turbulence (set to -1 for infinite scale).

    Returns:
    phase_screen (ndarray): 2D array representing the generated phase screen.
    """
    # 2025.03.16 - Original version based on paola.pro (22.12.2022) by Laurent Jolissaint (HES-SO), lines 3659+

    #TODO: Make it work with L0 != -1

    # Prepare the x direction coordinate matrix
    x = np.linspace(-1, 1, dim).astype(np.double)#*dim/(2*dim)
    y = np.linspace(-1, 1, dim).astype(np.double)#*dim/(2*dim)
    fpcoohfx, fpcoohfy = np.meshgrid(x,y)

    # Prepare the radius matrix
    xpcoohf = np.zeros((dim,dim))
    for i in range(dim):
        for j in range(dim):
            xpcoohf[i][j]= np.sqrt(fpcoohfx[i][j]**2 + fpcoohfy[i][j]**2)  # pupil plane coordianate radius

    fpcoohfx *=dim # Scale the grid to the size of the pupil matrix
    fpcoohfy *=dim
    xpcoohf *=dim # Scale the grid to the size of the pupil matrix

    # Von Karmàn corrected PSD (L0 != -1) PAOLA:3666
    #Watm = np.zeros((dim,dim), dtype=np.double)
    if L0 == -1:
        #w = np.where(fpcoohf[0]**2+fpcoohf.T[0]**2 > 0)
        dphi = 6.883877*(xpcoohf/r0)**(5/3)
        Watm = 0.022896/(r0**(5/3))*np.sqrt(fpcoohfx**2+fpcoohfy**2)**(-11/3)
    else:
        w = np.where(xpcoohf != 0)
        dphi = 0.171661*(L0/r0)**(5/3)*((1.005635)-(2*np.pi*xpcoohf/L0))**(5/6)
                               #beselk(5/6,2*np.pi*xpcoohf/L0)
        Watm = 0.022896/(r0**(5/3))*(fpcoohfx**2+fpcoohfy**2 +1/(L0**2))**(-11/3)

    PSD = Watm

    return PSD