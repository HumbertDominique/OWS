import numpy as np
from ows import fouriertransform as mathft
from scipy.special import kv
from matplotlib import pyplot as plt # for developpment only


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
   fx = np.linspace(-1,1,dimmat)*(1./(2*dxp))
   freqX, freqY = np.meshgrid(fx,fx)
   freqR = np.sqrt(freqX**2 + freqY**2)  # pupil plane spatial frequency radius
   w = freqR.nonzero()
   PSD = np.zeros_like(freqR, dtype=np.float64)  # Ensure PSD is the same shape and dtype as freqRtype
   PSD[w] = 0.0229*r0l**(-5.0/3) *(freqR[w]**2 + (L0 != -1)/L0**2)**(-11.0/6) # Spatial power spectrum with outer scale

   return PSD

def phase_screen(PSD, dxp, SEED = None, PSF=True, PUPIL = True, PD = [0,0]):
   """
   Generates a phase screen base on a Phase Spectrum Density (PSD). The padding should be taken into account with PD[0] <= PSD.shape[0]//2. Default is PD[0] = PSD.shape[0]/2
   Parameters:
      dimmat (int): Size of the phase screen matrix. must be even?
      PSD (ndarray): phase spectrum density
      dxp (float): Pupil plane pixel size dxp = WL/(NA*Npx)
      SEED (int): Seed for the random number generator. (default None)
      PSF (bool): option to return the PSF associated with the phase screen
      PUPIL (bool): True will crop it to a pupil of diameter dimmat if PD = [0,0]
      PD (ndarray of int): Pupil diameters in meters. PD[outer diameter-1, inner diameter)

   Returns (list):
      phase_screen (ndarray): real 2D array representing the generated phase screen.
      psf (ndarray): real 2D array representing the PSF associated with the phase screen
      MASK (ndarray): 0, 1 integer array representing the pupil.
      R (ndarray): matrix radius
   TODO: Add a dimmat argument to allow more pixels on the phase screen and PSF than PSD.
   TODO: Make it so the pupil diameters can be introduced either as pixel values or discances D [m] = D[px]/dxp (Needs to be checked)

   Based on WaveSeeingLimited.pro, Laurent Jolissaint, March 24, 2025 
   """
   # ------------ Input checks ------------
   # ------------ Input checks ------------
   # ------------ Input checks ------------
   if PSD.shape[0] != PSD.shape[1]:
      raise ValueError("PSD is not square")
   if dxp <= 0:
      raise ValueError("dxp must be positive")
   if type(PD[0]) != int:
      raise TypeError("Pupil diameters must be integers (for now)")
   if type(PD[1]) != int:
      raise TypeError("Pupil diameters must be integers (for now)")
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
   rng = np.random.default_rng(seed=SEED)

   pxR = np.linspace(-1 ,1,dimmat)*dxp*dimmat/2  
   xx, yy = np.meshgrid(pxR,pxR)
   R = np.sqrt(xx**2 + yy**2)  # pupil plane spatial frequency radius

   ## grid for pupil mask
   pxR = np.linspace(-dimmat//2 ,dimmat//2,dimmat)
   xx, yy = np.meshgrid(pxR,pxR)
   Rpx = np.sqrt(xx**2 + yy**2)  # Pupil radius in pixel

   PP = np.zeros((dimmat+1, dimmat+1)) # Phase power [dimmat+1, dimmat+1] in order to have a pixel in the middle [rad^2/m^-2]
   PP[0:dimmat,0:dimmat] = PSD*(dimmat*dxp)**2
   PP[0 : dimmat, dimmat] = PP[0 : dimmat, 0] # last row = 1st row
   PP[dimmat,0 : dimmat] = PP[0, 0 : dimmat] # last column = 1st column
   PP[dimmat, dimmat] = PP[0, 0] # bottom right corner = top left corner

   PSA = np.sqrt(PP)* rng.normal(size=(dimmat+1, dimmat+1)) # random draw to create a random phase spectrum [rad/m^(-1)] AMPLITUDE
   PSA = np.sqrt(2)*(PSA+np.rot90(PSA,2))/2 # forced even (in frequency) amplitude of wf spectrum
   PSop = rng.uniform(0, 1, size=(dimmat+1, dimmat+1))*2*np.pi # random draw for the phase of the spectrum of the optical phase
   
   # in order to make sure that the phase is real,
   # (1) the real part of the phase spectrum is forced to be even (as in "f(x) is even")
   # (2) the imaginary part of phase spectrum is forced to be odd (as in "f(x) is odd")
   phaseft =  (PSA * np.cos(0.5 * (PSop - np.rot90(PSop,2))))[0:dimmat,0:dimmat] + 1j*(PSA * np.sin(0.5 * (PSop - np.rot90(PSop,2))))[0:dimmat,0:dimmat]
   
   # ------------ PUPIL OPTION ------------
   # ------------ PUPIL OPTION ------------
   # ------------ PUPIL OPTION ------------
   MASK = np.ones((dimmat,dimmat))
   if PUPIL is True:
      if PD[0] == 0:
        PD[0] = dimmat//2
        MASK[Rpx >= PD[0]/2 +1] = 0
      else:
         MASK[Rpx >= PD[0]//2 +1] = 0
         MASK[Rpx < PD[1]//2] = 0

   phase_screen = mathft.ift2(phaseft,dxp).real
   phase_screen -= np.mean(phase_screen)
   #phase_screen[MASK == 1] -= np.mean(phase_screen[MASK == 1]) # if padding with the phsae screen istelf, cannot remove piston only un the pupil

   # ------------ PSF OPTION ------------
   # ------------ PSF OPTION ------------
   # ------------ PSF OPTION ------------
   if PSF == True:
      # padded_pupil = np.zeros((N_pad,N_pad))
      # padded_pupil[:dimmat,:dimmat] = MASK*phase_screen
      # padded_R = np.zeros((N_pad,N_pad))
      padded_pupil = MASK*phase_screen
      # padded_R = np.zeros((N_pad,N_pad))
      # padded_R[:dimmat,:dimmat] = MASK*R
      padded_R = MASK*R
      Sp = np.sum(R)*dxp**2
      apsf = mathft.ft2(padded_R*np.exp(-1j *padded_R*padded_pupil),delta=1./dxp)/Sp
      apsf = apsf[dimmat//2 - PD[0]//2:dimmat//2 + PD[0]//2,dimmat//2 - PD[0]//2:dimmat//2 + PD[0]//2]
      psf = np.abs(apsf)**2
   else :
      psf = None
   #return padded_pupil[dimmat//2 - PD[0]//2:dimmat//2 + PD[0]//2,dimmat//2 - PD[0]//2:dimmat//2 + PD[0]//2], psf, MASK, R
   return phase_screen, psf, MASK, R


def SHWFS(dimmat, PD,  N, F, wl, pupil_mask, R, phase_screen, dxp, lenslet_pad = 2):
    """
    This function calculates the SH wavefront sensor response.
    Parameters:
      dimmat: size of the sub-aperture in pixels (integer)
       N (int): number of lenslets along one dimension
       F (float): focal length of the wavefront sensor
       pupil_mask (ndarray): mask for the pupil, 1 inside and 0 outside
       phase_screen (ndarray): phase screen to be illuminated by the wavefront sensor
       lensled_pad (int): lenslet padding for the computation
    Returns:
       Lightfield (ndarray): Focal plane image of the phase screen through the lenslet array
       DwDx (ndarray): Wavefront slope along x
       DwDy (ndarray): Wavefront slope along y
    """

    #M = N # to implement on rectangular sub-apertures
    n = int(np.ceil(PD[0]/N)) # By definition, the number of lenslet must be an integer
    if N%2 == 0:
        padded_dim = n*(N) # N+1 in order to have a lenslet centered on the optical axis (not implemented here)
    else:
        padded_dim = n*(N+1)
    # Pad the pupil mask for illumination computation
    padded_mask = np.zeros((padded_dim,padded_dim),int)
    padded_mask = pupil_mask[dimmat//2-padded_dim//2:dimmat//2+padded_dim//2,dimmat//2-padded_dim//2:dimmat//2+padded_dim//2] # for use with "0-padding"
    # padded_mask = pupil_mask

    padded_phase_screen = np.zeros_like(padded_mask,np.float64)
    padded_phase_screen = phase_screen[dimmat//2-padded_dim//2:dimmat//2+padded_dim//2,dimmat//2-padded_dim//2:dimmat//2+padded_dim//2] # for use with "0-padding"
    # padded_phase_screen = phase_screen

    padded_R = np.zeros_like(padded_phase_screen,np.float64)
    padded_R = R[dimmat//2-padded_dim//2:dimmat//2+padded_dim//2,dimmat//2-padded_dim//2:dimmat//2+padded_dim//2] # for use with "0-padding"
    # padded_R = R
    padded_lightfield = np.zeros_like(padded_phase_screen,np.float64)

    lenslet_mask = np.zeros_like(padded_mask)
    active_lenslet_matrix = np.zeros((N,N))

	# Compute the "active" lenslets
    ii=0
    for i in range(0, padded_dim, n):
        jj=0
        for j in range(0, padded_dim, n):
            lenslet = padded_mask[i:i+n, j:j+n]
            if (lenslet.sum() >= (n*n)/2):
                lenslet_mask[i:i+2, j:j+n] = 1
                active_lenslet_matrix[ii,jj] = 1
            jj += 1
        ii +=1

    k_max = N*N

    sub_ps = np.zeros((lenslet_pad*n,lenslet_pad*n))
    sub_R = np.zeros_like(sub_ps)

    active_lenslet_matrix_1D = active_lenslet_matrix.reshape(N*N)
    Xc = np.zeros((active_lenslet_matrix_1D.shape[0]))
    Yc = np.zeros_like(Xc)
    k = 0
    for j in range(0, padded_dim, n):
        for i in range(0, padded_dim, n):
            if int(active_lenslet_matrix_1D[k]) == 1:
                Xc[k] = i+n/2
                Yc[k] = j+n/2
                sub_ps[0:n,0:n] = padded_phase_screen[i:i+n,j:j+n]
                sub_R[0:n,0:n] = padded_R[i:i+n,j:j+n]  
                sub_sp = np.sum(sub_R)*dxp**2
                apsf = (mathft.ft2(sub_R*np.exp(-1j *sub_R*sub_ps),delta=1./dxp)/sub_sp)/(wl*F) # focal plane through lenslet U = F{R*ps}/(lambda*f)
                padded_lightfield[i:i+n,j:j+n] = np.abs(apsf[lenslet_pad*n//2-n//2:lenslet_pad*n//2+n//2,lenslet_pad*n//2-n//2:lenslet_pad*n//2+n//2])**2
            k +=1

    ### 2. Computing the spot centres
    Xr = np.zeros_like(Xc)
    Yr = np.zeros_like(Xc)
    Xs = np.linspace(0 ,padded_dim-1,padded_dim)
    xx, yy = np.meshgrid(Xs,Xs)

    XI = xx*padded_lightfield
    YI = yy*padded_lightfield
    k = 0
    for i in range(0, padded_dim, n):
        for j in range(0, padded_dim, n):
            if int(active_lenslet_matrix_1D[k]) == 1:
                Xr[k] = XI[i:i+n,j:j+n].sum()/padded_lightfield[i:i+n,j:j+n].sum()
                Yr[k] = YI[i:i+n,j:j+n].sum()/padded_lightfield[i:i+n,j:j+n].sum()
            k +=1

    ### 3. Extract the WF slopes
    DwDx = ((Xr-Xc)/F).reshape(N,N)
    DwDy = ((Yr-Yc)/F).reshape(N,N)

    return padded_lightfield, DwDx, DwDy


def phase_structure_function(rho, r0, L0):
   """
   Generates a phase screen base on a Phase Spectrum Density (PSD)
   Parameters:
       rho (int): Size of the phase screen matrix. must be even?
       r0 (ndarray): Fried's parameter
       L0 (float): Turbulence outer scale
   Returns (list):
       D_phi (ndarray): 2D array representing the phase structure
   TODO: Error check
   Based on "STEP-BY-STEP PROCEDURE FOR COMPUTING NUMERICALLY THE SEEING LIMITED POINT SPREAD FUNCTION FROM THE OPTICAL TURBULENCE PHASE STRUCTURE FUNCTION EQUATION" Laurent Jolissaint
   """
   # Constants for the von Karman model
   a = 0.17166136 * (L0/r0)**(5/3)
   b = 1.0056349
   # Handle the case where rho is zero
   D_phi = np.zeros_like(rho)
   mask = rho > 0
   # Apply the formula for non-zero rho
   x = 2*np.pi*rho[mask]/L0
   D_phi[mask] = a * (b - x**(5/6) * kv(5/6, x))
   return D_phi

def telescope_otf(nu_n, epsilon=None):
   """
   Computes a telescope's OTF
   Parameters:
       nu_n (ndarray): Normalized spatial frequencies
       epsilon (float): obtruation ratio
   Returns (list):
       otf_tsc (ndarray): real 2D array representing the telescope OTF.
   TODO: Error check
   TODO: Implement anular pupils
   Based on "STEP-BY-STEP PROCEDURE FOR COMPUTING NUMERICALLY THE SEEING LIMITED POINT SPREAD FUNCTION FROM THE OPTICAL TURBULENCE PHASE STRUCTURE FUNCTION EQUATION" Laurent Jolissaint
   """
   # Plain pupil case (no central obtruation)
   otf_tsc = np.zeros_like(nu_n)
   mask = nu_n <= 1
   otf_tsc[mask] = (2/np.pi) * (np.arccos(nu_n[mask]) - nu_n[mask] * np.sqrt(1 - nu_n[mask]**2))
   return otf_tsc
    
def atmospheric_otf(nu, r0, L0, wavelength):
   """
   Computes the atmoshpere's OTF
   Parameters:
       nu_n (ndarray): Normalized spatial frequencies
   Returns (list):
       otf_tsc (ndarray): real 2D array representing the telescope OTF.
   TODO: Error check
   Based on "STEP-BY-STEP PROCEDURE FOR COMPUTING NUMERICALLY THE SEEING LIMITED POINT SPREAD FUNCTION FROM THE OPTICAL TURBULENCE PHASE STRUCTURE FUNCTION EQUATION" Laurent Jolissaint
   """
   rho = wavelength * nu
   D_phi = phase_structure_function(rho, r0, L0)
   otf_atm = np.exp(-0.5 * D_phi)
   return otf_atm

def compute_diffLim_psf(PD, dimmat, dxp):
    MASK = np.ones((dimmat,dimmat))
    pxR = np.linspace(-dimmat//2 ,dimmat//2,dimmat)
    xx, yy = np.meshgrid(pxR,pxR)
    Rpx = np.sqrt(xx**2 + yy**2)  # Pupil radius in pixel for mas

    pxR = np.linspace(-1 ,1,dimmat)*dxp*dimmat/2  
    xx, yy = np.meshgrid(pxR,pxR)
    R = np.sqrt(xx**2 + yy**2)  # pupil plane spatial frequency radius

    if PD[0] == 0:
        PD[0] = dimmat//2
        MASK[Rpx >= PD[0]/2 +1] = 0
    else:
        MASK[Rpx >= PD[0]//4 +1] = 0
        MASK[Rpx < PD[1]//4] = 0

    padded_R = MASK*R
    Sp = np.sum(R)*dxp**2

    apsf = mathft.ft2(padded_R*np.exp(-1j *padded_R),delta=1./dxp)/Sp
    diff_psf = np.abs(apsf)**2
    return diff_psf

def normalize(data):
    min_val = np.min(data)
    max_val = np.max(data)
    return (data - min_val) / (max_val - min_val)


def pixel_adder(data, scale_factor = [1,1], final_shape = None):
    """
    Generates a phase screen base on a Phase Spectrum Density (PSD)
    Parameters:
        data (ndarray): source data
        scale_factor (floar): scalling factor. Can be > or < o
        final_shape (ndarray): Option to input directly the desired shape. final_shape/data.shape must be an integer or the final shape will be data.shape*(final_shape//data.shape).
    Returns (list):
        D_phi (ndarray): 2D array representing the phase structure
    TODO: Error check
    """

    datashape = data.shape

    if final_shape == None:
        iscale = scale_factor[0]
        jscale = scale_factor[1]
    else:
        iscale = datashape[0]*(final_shape//datashape[0])
        jscale = datashape[1]*(final_shape//datashape[1])
    
    image_scaled = np.zeros((int(scale_factor[0]*datashape[0]), int(scale_factor[1]*datashape[1])))
    for i in range(0, datashape[0]):
        for j in range(0, datashape[1]):
            iscaled = i*iscale
            jscaled = j*jscale
            image_scaled[int(iscaled):int(iscaled+iscale), int(jscaled):int(jscaled+jscale)] = data[i,j]

    return image_scaled

