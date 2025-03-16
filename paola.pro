;+
;FUNCTION PAOLA - December 22, 2022
;
;    Laurent Jolissaint
;    University of Applied Sciences Western Switzerland, HES-SO
;    Haute Ecole d'Ingenierie et de Gestion du Canton de Vaud, HEIG-VD
;    Yverdon-les-Bains, Switzerland
;
;       Seeing limited mode: see section SEEING LIMITED MODE
;               NGS AO mode: see section NATURAL GUIDE STAR MODE
;                 GLAO mode: see section GROUND LAYER AO MODE
;
;    To compute PSF metrics: see section PSF METRICS CALCULATION
;         Saving the result: see section SAVING THE RESULT
;
;CALLING SEQUENCE
;
;    ====================================SEEING LIMITED=============================
;
;    res=PAOLA('seli',DIM,TSC,W0,L0,ZA,[options])
;    --------------------------------------------
;    res=PAOLA('seli',DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,WIND,[options])
;
;    =========================================NGS===================================
;
;    res=PAOLA('ngs',DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,WIND,DM_PARAMS,WFS_PARAMS,$
;              NGS_ANG,NGS_ORI,WFS_INT,LAG,LOOP_MODE,LOOP_FREQ,LOOP_GAIN,WFS_NEA,$
;              [options])
;    -------------------------------------------------------------------------------
;    res=PAOLA('ngs',DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,WIND,DM_PARAMS,WFS_PARAMS,$
;              NGS_ANG,NGS_ORI,WFS_INT,LAG,LOOP_MODE,LOOP_FREQ,LOOP_GAIN,NGS_MAG,$
;              FILTER,NGS_TEM,[options])
;
;    ========================================GLAO===================================
;
;    "full" and "edge" modes:
;    res=PAOLA('glao',DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,DM_PARAMS,WFS_PARAMS,$
;              SO_ANG,SO_ORI,GLAO_WFS,[options])
;    -------------------------------------------------------------------------------
;    "star" mode, where the WFS NEA is given:
;    res=PAOLA('glao',DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,WIND,DM_PARAMS,WFS_PARAMS,$
;              SO_ANG,SO_ORI,WFS_INT,LAG,LOOP_MODE,LOOP_GAIN,GLAO_WFS,GS_WEIGHT,$
;              WFS_NEA,[options])
;    -------------------------------------------------------------------------------
;    "star" mode, where NGS mag & temp & WFS RON & WFS_TAU are given:
;    res=PAOLA('glao',DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,WIND,DM_PARAMS,WFS_PARAMS,$
;              SO_ANG,SO_ORI,WFS_INT,LAG,LOOP_MODE,LOOP_GAIN,GLAO_WFS,GS_WEIGHT,$
;              NGS_MAG,FILTER,NGS_TEM,[options])
;
;    =======================================OPTIONS=================================
;
;                 Set...   To get...
;            /ANTI_ALIAS   WFS aliasing spatial filtering
;            DISPERSION=   Refractive index dispersion parameter. Default is 1.
;      /EE_ANALYSIS_DISC   Encircled energy width (for 50% and 80% EE)
;      /EE_ANALYSIS_SQUA   Ensquared energy width (for 50% and 80% EE)
;      /EE_ANALYSIS_SLIT   Enslited energy width  (for 50% and 80% EE)
;              EEW_DISC=   Fraction of energy within disc diameter EEW_DISC
;              EEW_SQUA=   Fraction of energy within square width EEW_SQUA
;              EEW_SLIT=   Fraction of energy within slit width EEW_SLIT
;              FITSCODE=   Create PSF/OTF/SF/PSD FITS files, XXX<FITSCODE>.fits
;                 /FRFFT   Use the FRACTIONAL_FFT & FFTW algorithms, if installed.
;         /FWHM_ANALYSIS   Min/max FWHM elliptic section on PSF
;                  /INFO   Main parameters and results displayed on screen.
;               LOGCODE=   Main parameters and results written in 'paola<LOGCODE>.log'
;         MAX_LOOP_FREQ=   Maximum acceptable loop frequency during optimization
;    /OPTIMIZE_LOOP_FREQ   Loop frequency optimization
;    /OPTIMIZE_LOOP_GAIN   Loop gain optimization
;      /OPTIMIZE_WFS_INT   WFS integration time optimization, GLAO mode
;          /OPTIMIZE_ALL   Loop frequency AND loop gain optimization
;                   /OTF   OTF in the output structure variable
;                   /PSD   Phase PSDs in the output structure variable
;              /ONLY_PSD   ONLY the phase PSD are computed, not the SF, OTF, PSFs
;                   /PSF   PSF in the output structure variable
;                 /X_PSF   PSF section along x-axis in the output structure variable
;                 /Y_PSF   PSF section along y-axis in the output structure variable
;          /POST_TIPTILT   Post AO G-tilt correction, as with an OT-CCD
;       TILT_ANGLE_STEP=   G-tilt correction by integer step of width TILT_ANGLE_STEP
;               /RADDMTF   To get a radial DM transfer function at the AO radius 0.5/pitch
;         /SCINTILLATION   To include scintillation effects in the OTF
;                    /SF   Phase structure function in output structure variable
;             VIBRATION=   Jitter of the PSF to simulate instrument vibration (all modes)
;                  /WAVE   Prepare a residual phase PSD array for use
;                          in wave.pro. You get the PSD in a fits file 
;                          PSD_WAVE***.fits and in the output structure variable.
;
;NOTE ----- TO SPEED UP THE COMPUTATION
;
;    FRACTIONAL FFT is optionally implemented in this version of the
;    code, as well as the FFTW algorithm. These tools makes the code
;    significantly faster for VERY large off-axis angles (close to
;    half a degree or more). IN NORMAL CONDITIONS, THERE IS NO NEED TO 
;    USE THESE OPTIONS.
;
;    Now if you need this, you will have to
;    (1) install the fractional FFT and FFTW (C code) libraries which
;    are not part of the basic PAOLA distribution, so you can ask me
;    at laurent.jolissaint@heig-vd.ch;
;    (2) set the keyword /FRFFT in the call to paola.
;
;    By default, the fractional FFT and FFTW are not used in the code.
;    Note that we are running an implementation of the fractional FFT
;    developed by Visa Korkiakoski and Christophe Verinaud.
;
;NOTE ----- DEFAULT VALUES FOR THE INPUTS : SET THE INPUT TO -1
;
;    Many inputs may have default values, for instance representing a
;    standard in the AO business, or the most reasonable first guess
;    of what must be the optimal value of a given parameter. In this
;    case, of your are willing to use the default value for a given
;    input, maybe because you are not familiar enough with AO, or any
;    other reason, set the said input to the value -1 and PAOLA will
;    do the rest. You can then check the screen output (/info) or the
;    output variable of PAOLA to see what is the default value.
;
;###################################################################################
;###################################################################################
;#############################                               #######################
;#############################      SEEING LIMITED MODE      #######################
;#############################                               #######################
;###################################################################################
;###################################################################################
;
;CALLING SEQUENCE
;
;    Usage 1: as the seeing limited PSF does not depend on the Cn2 distribution
;    neither on the wind velocity, but only on the seeing angle and the outer
;    scale, the following sequence is enough to get the seeing limited PSF:
;
;    res=PAOLA(MODE,DIM,TSC,W0,L0,ZA,[options])
;
;    where MODE is the calculation mode (seeing limited here), DIM indicates the
;    dimensions of the calculation, TSC gives the telescope OTF, W0 is the
;    seeing angle, L0 the outer scale of optical turbulence, and ZA the zenith
;    angle. 
;
;    Usage 2: if you have the Cn2 and wind velocity profiles and want to know
;    the associated mean layers altitude and isoplanatic angle, the average 
;    wind velocity and the life time of the speckles, use the following calling
;    sequence: 
;
;    res=PAOLA(MODE,DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,WIND,[options])
;
;    where HEIGHT is the layers altitudes w.r.t. the telescope pupil, DISTCN2 is
;    the vertical distribution of the optical turbulence strength, and WIND
;    indicates the WIND vertical profile in both directions (x & y).
;
;    See batch.pro file for examples - included in PAOLA toolbox distribution.
;
;INPUTS name | type | unit
;
;    MODE | STRING | -
;    Enter 'seli' for seeing limited.
;
;    DIM | structure variable | -
;    DIM is a structure variable where the pixels and matrices sizes are
;    kept. See output of function PIXMATSIZE.PRO
;
;    TSC | structure variable | -
;    TSC is a structure variable including the telescope optical parameters and
;    OTF/PSF/PHASE/PUPIL. See output of function PSFOTFTSC.PRO
;
;    W0 | REAL SCALAR | ARCSEC
;    Median SEEING angle @ 500 nm, at zenith. This value is used to compute the
;    Fried parameter r0 for the specified wavelength. Zenith value must be given.
;
;    L0 | REAL SCALAR | METRE
;    Outer scale of optical turbulence. To get an infinite outer scale, set this 
;    parameter to -1. Note that an outer scale larger than 10 km is considered 
;    infinite.
;
;    ZA | REAL SCALAR | DEGREE
;    Zenith angle of the observation. The apparent seeing angle and turbulent
;    layers altitude change with the zenith angle and are modified
;    accordingly within the code. ZA must be in the range [0,90[ degrees,
;    i.e. the 90 deg value is excluded (it gives an infinite seeing angle).
;
;OPTIONAL INPUTS name | type | unit | default if any
;
;    In seeing limited mode, in principle the only thing you need to know to compute 
;    the long exposure OTF/PSF is the Fried parameter r0 and the outer scale L0,
;    described above. Therefore, the optional inputs below are useful only if
;    you want PAOLA to compute
;    (1) the mean turbulent layers altitude dispersion relative to the telescope
;    entrance pupil and the isoplanatic angle, and
;    (2) the mean turbulent layers velocities and the optical phase aberration
;    life-time.
;    Scalar values can be provided, if for instance you assume the turbulence to
;    be concentrated into a single layer (in the telescope pupil or at any
;    altitude). 
;
;    HEIGHT | REAL ARRAY(number of layers) OR SCALAR | METRE
;    altitude of each of the NL considered turbulent layers, and this is
;    >>>>> RELATIVE TO THE TELESCOPE PUPIL <<<<<<<. Note that it can be a scalar 
;    (for a simplified 1-layer profile, for instance ground layer).
;
;    DISTCN2 | REAL ARRAY(number of layers) OR SCALAR | 1
;    relative distribution of the Cn^2*dh profile over the NL turbulent layers (dh 
;    is the layer thickness). It does not have to be normalized, i.e total(DISTCN2) 
;    can be different that 1, and the Cn2*dh profile itself can be given here. The 
;    actual amplitude of the Cn2*dh profile will be set inside PAOLA according to 
;    the provided seeing angle W0. Can be a scalar (1-layer profile).
;
;    WIND | REAL ARRAY(number of layers,2) | METRE/SEC
;    WIND[*,0]: x-axis component (east->west dir/sky). 
;    WIND[*,1]: y-axis component (south->north dir/sky).
;    Must be introduced that way: [[x components],[y components]]
;    example for 3 layers : [[10,0,5],[0,-1,-5]]. Can be given as a
;    scalar, too. Then we will consider that WIND corresponds to the
;    x-coordinate of the wind velocity of a single layer (1-layer
;    profile case).
;    NOTE: WIND VELOCITY VECTOR IS NOT ADAPTED TO THE ZENITH DISTANCE
;    BECAUSE THIS ADAPTATION WOULD DEPENDS ON THE POINTING POSITION ON THE
;    SKY, AND THAT WOULD MAKE EVERYTHING UNNECESSARY COMPLICATED. IF YOU
;    WANT TO STUDY THIS EFFECT, JUST ADAPT THE WIND VECTOR COEFFICIENTS TO
;    YOUR PARTICULAR POINTING POSITION IN THE SKY, WHEN YOU DECLARE YOUR
;    WIND VECTOR (OUTSIDE OF PAOLA, IN YOUR BATCH).
;
;    VIBRATION= | A 3 COMPONENTS STRUCTURE VARIABLE | - | -
;    |  Vibrations induced on the beam, as seen in the image plane. Compatible
;    |  with the vibration of the telescope beam. This input is useful to
;    |  simulate a jitter occuring in the science arm, after the telescope,
;    |  independant from the telescope jitter. So you may have separate telescope
;    |  and instrument vibrations.
;    |
;    |__.SIGMA_TTX | POSITIVE REAL SCALAR | ASEC
;    |  Tip jitter RMS, as seen in the focal plane (sky).
;    |
;    |__.SIGMA_TTY | POSITIVE REAL SCALAR | ASEC
;    |  Tilt jitter RMS, as seen in the focal plane (sky). 
;    |
;    |__.ORIENTATION | REAL SCALAR | DEGREE
;    |  Orientation angle of the tip (TTx) perturbation. Counterclockwise.
;    |  +90 corresponds to the positive y-axis, +180 to the negative
;    |  x-axis. The tilt jitter is perpendicular to this direction. 
;    |
;    |--SYNTAX : VIBRATION={SIGMA_TTX:value,SIGMA_TTY:value,ORIENTATION:angle}
;
;KEYWORDS
;
;    /SCINTILLATION to take into account the impact of the electromagnetic field
;    amplitude fluctuation on the PSF. This is a 2nd order effect that grows
;    with the layers altitudes. Works only if a Cn2 profile is provided. See
;    HEIGHT and DISTCN2 optional inputs above. 
;
;    /INFO print the parameters and results of the calculation.
;
;    /PSF to get the PSF at output.
;
;    /X_PSF If set, only the section of the PSF along the x-axis (horizontal) is
;    calculated and given in the output structure variable. Saves a lot of
;    computation time.
;
;    /Y_PSF same as above, for y-axis (vertical).
;
;    /OTF to get the OTF at output.
;
;    /SF  to get the phase Structure Function at output.
;
;    /PSD to get the Phase Power Spectrum at output.
;
;    /ONLY_PSD same as above, but no structure function nor PSF or OTF are 
;    calculated, just the PSDs. Use this option if you just want PAOLA to compute 
;    the PSD, and nothing else. Incompatible with the keywords FWHM_ANALYSIS, EE 
;    ANALYSIS & co., PSF, OTF, SF.
;
;    /WAVE to get the phase power spectral density inside the output
;    structure variable. This is useful to run the function wave.pro
;    to generate instantaneous corrected phase screens.
;
;OUTPUTS name | type | unit
;
;    The output is a structure variable, whose components are defined
;    as follow:
;
;    .lam | REAL SCALAR | MICROMETERS
;    Wavelength used in the calculation.
;
;    .inst | STRING | -
;    Type of focal instrumentation. 'IMAGER' or 'INTENSITY'. If 'IMAGER', we
;    assume that the focal plane pixel size is the one set in DIM input (dxf),
;    and the PSF is averaged within the detector pixels (it's actually a
;    convolution). If 'INTENSITY', there is no averaging, the PSF is
;    given as a sampling of the intensity field in the image plane.
;
;    .mode | STRING | -
;    Calculation mode, copy of input MODE.
;
;    .w0ZA | REAL SCALAR | METER
;    Seeing angle @ 500 nm for zenith angle ZA.
;
;    .r05 | REAL SCALAR | METER
;    Fried's parameter @ 500 nm for zenith angle ZA.
;
;    .r0l | REAL SCALAR | METER
;    Fried's parameter at the imaging wavelength for zenith angle ZA.
;
;    .L0  | REAL SCALAR | METER
;    Outer scale of turbulence. -1 if infinite.
;
;    .strehl | REAL SCALAR | -
;    >>> NOT IF /ONLY_PSD IS SET <<<
;    Strehl ratio.
;
;    .cfr | REAL SCALAR | -
;    >>> NOT IF /ONLY_PSD IS SET <<<
;    Ratio between the practical cut-off angular frequency in the focal plane, and 
;    the cut-off frequency of the telescope without aberrations. Between 0 and 1.
;
;OPTIONAL OUTPUTS name | type | unit
;
;    If the turbulence profile is specified:
;
;      .alt | REAL SCALAR | METER
;      mean altitude of the turbulent layers for zenith angle ZA.
;
;      .ani | REAL SCALAR | ASEC
;      anisoplanatism angle for zenith angle ZA.
;
;      .win | REAL SCALAR | METER/SECOND
;      mean wind speed over the layers (NOT adapted to zenith angle).
;
;      .lft | REAL SCALAR | MILLISECOND
;      life time of the speckle field for zenith angle ZA.
;
;    If you have requested the calculation of several metrics of the AO PSF,
;    these metrics will be added to the output structure variable as well. See
;    section << PSF METRICS CALCULATION >> for details.
;
;    IF /SCINTILLATION IS SET
;
;      .sci_index | REAL SCALAR | 1
;      Scintillation index - variance of the relative irradiance fluctuations, 4 
;      times the variance of the relative amplitude fluctuations). See Ref [9].
;
;    IF /PSF OR /X_PSF OR /Y_PSF ARE SET
;
;      .dxf_usr | REAL SCALAR | ASEC/PX
;      PSF angular pixel scale.
;
;    IF /PSF is set
;
;      .psf | REAL ARRAY(N,N) | STREHL RATIO
;      PSF for the selected mode.
;
;    IF /X_PSF OR /Y_PSF are set
;
;      .psfx | REAL ARRAY(N) | STREHL RATIO
;      x-cut section of the PSF for the selected mode.
;
;      .psfy | REAL ARRAY(N) | STREHL RATIO
;      y-cut section of the PSF for the selected mode.
;
;    IF /OTF IS SET
;
;      .dff | REAL SCALAR | 1/RAD/PX | ALL
;      Spatial frequency pixel scale in the OTF matrix.
;
;      .otf | REAL ARRAY(N,N) | MAX OTF = 1
;      OTF for the selected mode.
;
;    IF /SF IS SET
;
;      .dxp | REAL SCALAR | METERS/PIXEL
;      Pixel scale in the pupil plane associated with the structure functions 
;      matrices.
;
;      .sf | REAL ARRAY(N,N) | RAD^2
;      seeing limited SF.
;
;    IF /PSD OR /ONLY_PSD ARE SET
;
;      .psd_atm | REAL ARRAY(N,N) | RAD^2*M^2
;      Turbulent phase spatial power spectrum.
;
;      .dfp_hf | REAL SCALAR | 1/METERS/PX
;      Spatial frequency Pixel scale in the pupil, associated with the
;      PSD matrices. The tag _hf means high frequency and is here for
;      compatibility with the NGS and GLAO modes.
;
;    IF /PSD OR /ONLY_PSD AND /SCINTILLATION ARE SET
;
;      .psd_amp | REAL 2D ARRAY | RAD^2*M^2
;      Power spectrum of the relative fluctuation of the wave amplitude.
;
;    IF /WAVE IS SET
;
;      .psd_wave | REAL ARRAY(N,N) | RAD^2*M^2
;      Turbulent phase spatial power spectrum.
;
;      .pup_wave | REAL ARRAY(N,N) | RAD^2*M^2
;      >>> ONLY IF /PUPIL WAS SET IN THE CALL TO psfotftsc.pro <<<
;      Pupil mask.
;
;      .dfp_wave | REAL SCALAR | 1/METERS/PX
;      Spatial frequency pixel scale in the pupil, associated with the
;      psd_wave matrix.
;
;      .dxp_wave | REAL SCALAR | METERS/PX
;      Spatial pixel scale in the pupil.
;
;    IF VIBRATION IS SET
;
;      .VIBRATION : a copy of optional input VIBRATION.
;
;###################################################################################
;###################################################################################
;############################                               ########################
;############################    NATURAL GUIDE STAR MODE    ########################
;############################                               ########################
;###################################################################################
;###################################################################################
;
;CALLING SEQUENCE
;
;  USAGE 1:
;
;    Here a WFS noise equivalent angle (NEA) is provided, from which
;    the WFS noise phase variance is calculated and used to scale the WFS noise
;    spatial power spectrum.
;
;    res=PAOLA('ngs',DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,WIND,$
;              DM_PARAMS,WFS_PARAMS,DM_HEIGHT,NGS_ANG,NGS_ORI,WFS_INT,LAG,$
;              1000,0.5,'closed',WFS_NEA,[options])
;
;    Note that in the example above, loop is closed, loop frequency is
;    1 kHz, gain is 0.5. Open loop option is also available.
;
;  USAGE 2:
;
;    Here, instead of the WFS NEA, the magnitude in a given filter and
;    the NGS spectrum equivalent black body temperature are provided, plus the WFS
;    detector read noise and the total optical throughput from the star to the
;    WFS detector (see components of the WFS structure variable WFS_PARAMS). The
;    WFS noise phase variance is calculated from a physical model of the noise
;    propagation in the WFS.
;
;    res=PAOLA('ngs',DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,WIND,$
;              DM_PARAMS,WFS_PARAMS,DM_HEIGHT,NGS_ANG,NGS_ORI,WFS_INT,LAG,$
;              1000,0.5,'closed',NGS_MAG,FILTER,NGS_TEM,[options])
;
;  USAGE 3 & 4:
;
;    Same as above, but loop is OPEN, and gain is necessary 1.
;    LOOP_GAIN input is actually ignored, but is kept as an input ... to make
;    the code simpler.
;
;    res=PAOLA('ngs',DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,WIND,$
;              DM_PARAMS,WFS_PARAMS,DM_HEIGHT,NGS_ANG,NGS_ORI,WFS_INT,LAG,1000,1,$
;              'open',WFS_NEA,[options])
;
;    res=PAOLA('ngs',DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,WIND,$
;              DM_PARAMS,WFS_PARAMS,DM_HEIGHT,NGS_ANG,NGS_ORI,WFS_INT,LAG,1000,1,$
;              'open',NGS_MAG,FILTER,NGS_TEM,[options])
;
;    See the file batch.pro for examples (included in PAOLA toolbox distribution).
;
;INPUTS name | type | unit
;
;    MODE | STRING | -
;    AO mode. Here must be set to 'ngs'.
;
;    DIM | structure variable | -
;    DIM is a structure variable where the pixels and matrices sizes are
;    kept. See output of function PIXMATSIZE.PRO
;
;    TSC | structure variable | -
;    TSC is a structure variable including the telescope optical parameters and
;    OTF/PSF/PHASE/PUPIL. See output of function PSFOTFTSC.PRO
;
;    W0 | REAL SCALAR | ARCSEC
;    Median SEEING angle @ 500 nm, at zenith. This value is used to compute the
;    Fried parameter r0 for the specified wavelength. Zenith value must be given.
;
;    L0 | REAL SCALAR | METRE
;    Outer scale of optical turbulence. To get an infinite outer scale, set L0
;    to -1. Note that an outer scale larger than 10 km is considered infinite.
;
;    ZA | REAL SCALAR | DEGREE
;    Zenith angle of the observation. The apparent seeing angle and turbulent
;    layers altitude change with the zenith angle and are modified
;    accordingly within the code. ZA must be in the range [0,90[ degrees,
;    i.e. the 90 deg value is excluded (it gives an infinite seeing angle).
;
;    HEIGHT | REAL ARRAY(number of layers) OR SCALAR | METRE
;    altitude of each of the NL considered turbulent layers, and this is
;    >>>>> RELATIVE TO THE TELESCOPE PUPIL <<<<<<<. Note that it can be a scalar
;    (for a simplified 1-layer profile, for instance ground layer).
;
;    DISTCN2 | REAL ARRAY(number of layers) OR SCALAR | 1
;    relative distribution of the Cn^2*dh profile over the NL turbulent layers (dh 
;    is the layer thickness). It does not have to be normalized, i.e total(DISTCN2) 
;    can be different than 1, and the Cn2*dh profile itself can be given here. The 
;    actual amplitude of the Cn2*dh profile will be set inside PAOLA according to 
;    the provided seeing angle W0. Can be a scalar (1-layer profile).
;
;    WIND | REAL ARRAY(number of layers,2) | METRE/SEC
;    WIND[*,0]: x-axis component (east->west dir/sky). 
;    WIND[*,1]: y-axis component (south->north dir/sky).
;    Must be introduced this way: [[x components],[y components]]
;    example for 3 layers : [[10,0,5],[0,-1,-5]]. Can be given as a
;    scalar, too. Then we will consider that WIND corresponds to the
;    x-coordinate of the wind velocity of a single layer.
;    NOTE: WIND VELOCITY VECTOR IS NOT ADAPTED TO THE ZENITH DISTANCE
;    BECAUSE THIS ADAPTATION WOULD DEPENDS ON THE POINTING POSITION ON THE
;    SKY, AND THAT WOULD MAKE EVERYTHING UNNECESSARY COMPLICATED. IF YOU
;    WANT TO STUDY THIS EFFECT, JUST ADAPT THE WIND VECTOR COEFFICIENTS TO
;    YOUR PARTICULAR POINTING POSITION IN THE SKY, WHEN YOU DECLARE YOUR
;    WIND VECTOR (OUTSIDE OF PAOLA, IN YOUR BATCH).
;
;    DM_PARAMS | a 3 components structure variable | -
;    DM parameters structure variable. Check IDL user manual to see how to create a 
;    structure variable. See also example in batch.pro file. The components are:
;    |
;    |__DM_PARAMS.DM_HEIGHT | REAL SCALAR | METERS
;    |  Conjugated altitude of the DM with respect to the ground level.
;    |
;    |__DM_PARAMS.DMTF | REAL ARRAY(N,N) OR -1 | 1
;    |  DM **AMPLITUDE** spatial transfer function. To get a perfect DM filter, set 
;    |  this input to -1.
;    |  See function DM_MODEL.PRO user manual for instructions to build a DM
;    |  transfer function using a given set of influence functions.
;    |
;    |__DM_PARAMS.ACTPITCH | REAL SCALAR OR -1 | M
;    |  DM actuator pitch as seen in the pupil plane. To get the default value (r0 
;    |  at the imaging wavelength and at zenith), set this input to -1.
;    |
;    |  An important limitation: the transition between the diffraction limited core 
;    |  of the PSF and the residual atmospheric halo occurs at an off-axis angle 
;    |  lambda/(2*actpitch). Therefore, it is important that you choose
;    |  a PSF field-of-view (see PIXMATSIZE inputs) large enough to see this
;    |  transition. Otherwise, the Strehl performance will be overestimated. In my
;    |  opinion a FoV radius at least twice the transition angle (radius), i.e. a 
;    |  total FoV width of 4*lambda/(2*actpitch) is really a minimum. Mama PAOLA 
;    |  will warn you if this is not the case, and the calculation will not start.
;    |
;    |  There are three possible options for input DM_PARAMS:
;    |
;    |  1st option: DM_PARAMS.DMTF = -1
;    |              DM_PARAMS.ACTPITCH = -1
;    |
;    |  Here, we use a perfect DM filtering of the phase aberration below the DM 
;    |  cut-off frequency 1/(2*ACTPITCH), i.e. the DM transfer function is = 1 for
;    |  |fx|<f_c and |fy|<f_c and 0 outside of the square domain. Such a model is 
;    |  fairly accurate and works for most applications. Actuator pitch is set to 
;    |  the default value (see above).
;    |
;    |  2nd option: DM_PARAMS.DMTF = -1
;    |              DM_PARAMS.ACTPITCH = a positive value.
;    |
;    |  The DM is still an ideal one, and the actuator pitch is set to a specific 
;    |  value.
;    |
;    |  3rd option: DM_PARAMS.DMTF = a 2D DM transfer function model
;    |              DM_PARAMS.ACTPITCH = a positive value.
;    |
;    |  Here, a 2 dimensional realistic model of a DM transfer function is given, 
;    |  along 
;    |  with the corresponding actuator pitch. See usage of function DM_MODEL.PRO, 
;    |  that you should use to build the DM transfer function. DM_PARAMS.DMTF is 
;    |  dimensionless (it's a transfer function).
;    |
;    |  Note that for most AO correction simulation cases, the perfect DM model is 
;    |  good enough and does not produce PSF that are significantly different than 
;    |  PSF produced using a more sophisticated DM filter. ExAO modeling might be 
;    |  the only case where using sophisticated DM modeling matters (see usage of 
;    |  DM_MODEL.PRO function for more details on how to build such a DM transfer 
;    |__function).
;
;    WFS_PARAMS | a structure variable | -
;    WFS parameters structure variable. Check IDL user manual to see how to create a 
;    structure variable. See also example in batch.pro file. The components are:
;    |
;    |__WFS_PARAMS.WFS_PITCH | REAL SCALAR | METERS
;    |  WFS lenslet array pitch size, as seen from the entrance pupil plane. If set 
;    |  to -1, the WFS lenslet pitch is assumed to be equal to the DM actuator 
;    |  pitch, i.e. WFS_PITCH is set to DM_PARAMS.ACTPITCH within the code.
;    |
;    |__IF THE WFS NOISE IS GIVEN AS A NOISE EQUIVALENT ANGLE WITH THE INPUT
;       WFS_NEA, ONLY SET WFS_PARAMS.WFS_PITCH AND IGNORE THE OTHER
;       COMPONENTS BELOW, AS THEY ARE USELESS IN THIS CASE.
;       |
;       |__WFS_PARAMS.MIRCOATING | VECTOR OF STRING | [-]
;       |  This input set the number and type of mirrors from (and including) the
;       |  telescope primary mirror until the WFS lenslet array. It is used in the
;       |  calculation of the number of photons received by the WFS detector /
;       |  lenslet / frame.
;       |  There are three possible mirror coating: Aluminum, Silver, and Gold.
;       |  Set MIRCOATING to the suite of coatings along the optical train: 
;       |  For instance if your telescope is Al-Al-Ag then Al-Al-Au inside,
;       |  WFS_PARAMS.MIRCOATING=['Al','Al','Ag','Al','Al','Au']
;       |  but if you want to use the default case, which is 3 Aluminum mirrors & 3 Silver
;       |  mirrors with normal incidence, simply set MIRCOATING to -1.
;       |
;       |__WFS_PARAMS.MIRANGLE | VECTOR OF REAL | DEGREE
;       |  Beam incidence angle of the mirrors introduced above. For instance
;       |  WFS_PARAMS.MIRANGLE=[0,0,45,30,30,45]
;       |  and if you use the MIRCOATING default, set this input to -1 as well.
;       |
;       |__WFS_PARAMS.NBLENSES | INTEGER SCALAR | [0,N]
;       |  Number of REFRACTIVE surfaces from the entrance pupil (included) to the
;       |  detector. Set it to 0 if you want to ignore this input for now.
;       |
;       |__WFS_PARAMS.WFS_RON | REAL SCALAR | ELECTRON/PIXEL
;       |  Pixel read noise of the WFS CCD detector. If WIND has been set to 0, 
;       |  WFS_RON is ignored.
;       |
;       |__WFS_PARAMS.SKY_RADIANCE | REAL SCALAR | W/M^2/SR
;       |  Sky background radiance in watt / square meter / steradian
;       |
;       |__WFS_PARAMS.ALGORITHM | CHAR STRING | -
;       |  Type of algorithm for SH spot centroid calculation:
;       |  Enter '4Q' for a 4-quadrant mode;
;       |  Enter 'CG' for a center-of-gravity mode and diffraction limited spot.
;       |  Enter 'GS' for a center-of-gravity mode and Gaussian spot
;       |             whose FWHM is set by residual turbulence.
;       |
;       |__WFS_PARAMS.EXTRAFILTER - a structure variable
;       |  If, on top of the transmission losses due to atmosphere and optical
;       |  surfaces, there is in your system an other transmission loss that cannot 
;       |  be modeled with the inputs above, you can use the EXTRAFILTER input.
;       |  EXTRAFILTER is a structure variable, with two components:
;       |  .TAUEXT | REAL VECTOR | 1
;       |  The filter transmission, in the range [0,1] (i.e NOT in %)
;       |  .LAMBDA | REAL VECTOR | MICROMETERS
;       |  The wavelength associated with the filter data. Should not be
;       |  outside of the range 300 - 5000 nm.
;       |  The number of data points in TAUEXT and LAMBDA does not matter:
;       |  the code will either fit a spline to the data or rebin the data
;       |  to adjust to the internal wavelength sampling, which is 10 nm.
;       |  Syntax inside the setting of WFS_PARAMS:
;       |  WFS_PARAMS.EXTRAFILTER:{TAUEXT:some vector,LAMBDA:some vector}
;       |
;       |  -----> IF YOU DO NOT NEED THIS INPUT, SET IT TO 'no'
;       |         WFS_PARAMS.EXTRAFILTER:'no'
;       |
;       |__THE TWO INPUTS BELOW ARE ONLY TO BE SET IF YOU HAVE CHOSEN
;          A CENTER-OF-GRAVITY ALGORITHM (CG)
;          OR A GAUSSIAN SPOT (GS)FOR THE CENTROID MEASUREMENT. SKIP THESE IF
;          YOU HAVE A 4-QUADRANT ALGORITHM.
;          |
;          |__WFS_PARAMS.WFS_PXFOV | SCALAR INTEGER | PIXELS
;          |  Field of view for each lenslets, in pixels. Must be >= 6, because with 
;          |  less than that, non-linearities due to pixelization become an issue.
;          |
;          |__WFS_PARAMS.WFS_PXSIZE | SCALAR OR -1 | ASEC/PIXEL
;             WFS detector pixel size (projected on sky). Must be >= Nyquist limit
;             lambda/(2*WFS_PITCH), expressed in asec. Set WFS_PXSIZE to -1 to
;             get the Nyquist limit, which is the best choice anyway.
;
;    *********
;    IMPORTANT - WFS NOISE ERROR WITH CoG ALGORITHM WILL NEVER BE 0 EVEN
;                WITH AN INFINITE NUMBER OF PHOTONS : THIS IS BECAUSE
;                THE MODEL TAKES INTO ACCOUNT THE EFFECT OF SPOT
;                TRUNCATION DUE TO THE LIMITED WFS FIELD-OF-VIEW
;                SEE Sandrine Thomas et al. MNRAS 371, Issue 1, pp. 323-336
;    *********
;
;    NGS_ANG | REAL SCALAR | ASEC
;    science object off-axis angle.
;
;    NGS_ORI | REAL SCALAR | DEGREES
;    orientation of the science object position vector in the field of view, where 0 
;    correspond to an alignment with the +x axis, and +90 with the +y axis.
;
;    WFS_INT | REAL SCALAR OR -1 | MILLISECOND
;    integration time of the wave-front sensor. MUST be > 0.
;    NOTE If WFS_INT optimization is required, to balance the servo-lag and the noise
;    errors, set the keyword /OPTIMIZE_LOOP_FREQ. If you already have an initial guess
;    for the optimal WFS integration time, set WFS_INT to this initial value, and the
;    optimization will start with this. In this case, the user must provide the NGS
;    magnitude, black body temperature, and WFS CCD read noise.
;    WARNING : cannot be > 1 / loop frequency (see LOOP_FREQ).
;
;    LAG | REAL SCALAR | MILLISECOND
;    Technical servo system time lag  = read-out time of the WFS + computation
;    time. In other words, this is the time spent between the end of the
;    WFS exposure and the application of the updated command to the DM. 
;
;      CAUTION: the TOTAL DELAY between the WFS measurement and the DM
;      correction is equal to the WFS exposure time, or the loop
;      period, plus the technical LAG defined above. So, when entering
;      the LAG parameter, make sure that it is the technical delay,
;      and NOT the total delay. Inside PAOLA, the total delay is
;      computed by adding together the WFS integration time and LAG.
;
;    LOOP_MODE | CHAR STRING | -
;    Loop mode. 'OPEN' or 'CLOSED'.
;
;    LOOP_FREQ | scalar integer | Hz | DEFAULT: set it to -1
;    Loop frequency. It can be set to a different value than the
;    inverse of the WFS integration time in order to explore the
;    effect of the WFS phase averaging. But in principle you want this
;    frequency to be 1/(WFS integration time). If this is the case,
;    then set LOOP_FREQ to -1, so the default value
;    1000/(WFS integration time [ms]) will be set.
;    WARNING: cannot be > 1 / WFS integration time (see WFS_INT)
;
;    LOOP_GAIN | REAL SCALAR | 1
;    Loop gain. Must be within ]0,4.9348] because above this range
;    there are instabilities.
;    This input is ignored if LOOP_MODE is 'open', i.e. in
;    open loop the gain is necessarily 1. See /OPTIMIZE_LOOP_GAIN if you want to
;    optimize this value. In the later case, LOOP_GAIN is used as the starting
;    value for the optimization procedure.
;
;    The WFS noise error can be calculated either from the NGS properties, or
;    from a pre-defined Noise Equivalent Angle input. In the first case, use the
;    following three inputs:
;    |
;    |__NGS_MAG | REAL SCALAR | 1
;    |  Magnitude of the NGS, for a given filter (see FILTER input below). We use
;    |  here the Johnson filter system, and the Vega zero point magnitude
;    |  reference values.
;    |
;    |__FILTER | CHARACTER STRING | -
;    |  Tag name of the filter associated to the magnitude given in
;    |  NGS_MAG. Possible values are (Johnson system): U B V R I J H K L M N Q.
;    |  Syntax: do not forget to enter the tag name as a string, for instance 'V'.
;    |
;    |__NGS_TEM | REAL SCALAR | KELVIN
;       Black-body temperature associated to the NGS spectrum.
;    
;    In the second case, use the following input:
;    |
;    |__WFS_NEA | REAL SCALAR | ASEC
;       WFS Noise Equivalent Angle (same for each lenslet).
;
;OPTIONAL INPUTS name | type | unit | default
;
;    VIBRATION= | A 3 COMPONENTS STRUCTURE VARIABLE | - | -
;    |  Vibrations induced on the beam, as seen in the image plane. Compatible
;    |  with the vibration of the telescope beam. This input is useful to
;    |  simulate a jitter occuring in the science arm, after AO correction,
;    |  independant from the telescope jitter. So you may have separate telescope
;    |  and post-AO instrument vibrations. This is useful for coronagraphs analysis.
;    |
;    |__.SIGMA_TTX | POSITIVE REAL SCALAR | ASEC
;    |  Tip jitter RMS, as seen in the focal plane (sky).
;    |
;    |__.SIGMA_TTY | POSITIVE REAL SCALAR | ASEC
;    |  Tilt jitter RMS, as seen in the focal plane (sky). 
;    |
;    |__.ORIENTATION | REAL SCALAR | DEGREE
;    |  Orientation angle of the tip (TTx) perturbation. Counterclockwise.
;    |  +90 corresponds to the positive y-axis, +180 to the negative
;    |  x-axis. The tilt jitter is perpendicular to this direction. 
;    |
;    |--SYNTAX : VIBRATION={SIGMA_TTX:value,SIGMA_TTY:value,ORIENTATION:angle}
;
;    DISPERSION= | SCALAR REAL > 0 | 1 | 1
;    Dispersion factor = N(LAM_WFS) / N(imaging wavelength)
;                      = [n(LAM_WFS)-1] / [n(imaging wavelength)-1]
;    Note that users must use their own refractivity models N(lambda) to get the
;    dispersion factor. Syntax: DISPERSION=value > 0, but can be < 1.
;
;KEYWORDS
;
;    /RADDMTF To get a radial only DM transfer function at the AO radius 1/2/pitch
;
;    /ANTI_ALIAS if set, simulate spatial filtering of WFS aliasing by not adding 
;    the aliased component to the total power spectrum of the corrected phase.
;
;    OPTIMIZATION OF WFS INTEGRATION TIME OR LOOP GAIN
;    =================================================
;
;      IF THE WFS NOISE IS ALREADY FIXED WITH THE INPUT WFS_NEA, OR 
;      IF THE WIND SPEED IS NULL, OPTIMIZATION MAKES NO SENSE AND PAOLA
;      RETURNS AN ERROR MESSAGE.
;
;      /OPTIMIZE_LOOP_FREQ Loop frequency optimization: PAOLA seeks the optimal 
;      loop frequency to minimize the combined effect of WFS noise (NGS photon 
;      noise and CCD read noise), which increases with frequency and the servo-lag 
;      error, which decreases with the loop frequency.
;      **** NOTE ****
;      When using this option, IT IS RECOMMENDED to set the keyword
;      MAX_LOOP_FREQ: this will define a top value limiting the
;      exploration of the AMOEBA algorithm, to avoid unrealistically
;      high frequency values.
;
;      /OPTIMIZE_LOOP_GAIN optimization of loop gain to minimize the combined
;      effect of WFS noise and servo-lag. The loop gain input value, in the
;      iterative optimization procedure, is used as the starting value.
;
;      /OPTIMIZE_ALL optimization of both loop frequency and loop gain.
;      **** NOTE ****
;      When using this option, IT IS RECOMMENDED to set the keyword
;      MAX_LOOP_FREQ: this will define a top value limiting the
;      exploration of the AMOEBA algorithm, to avoid unrealistically
;      high frequency values.
;
;    /SCINTILLATION to take into account wave amplitude fluctuation. Works only if a 
;    Cn2 profile is provided. See HEIGHT and DISTCN2 optional inputs above.
;
;    /INFO print the parameters and results of the calculation.
;
;    /PSF to get the PSF at output.
;
;    /X_PSF If set, only the section of the PSF along the x-axis (horizontal) is
;    calculated and given in the output structure variable. Saves a lot of 
;    computation time.
;
;    /Y_PSF same as above, for y-axis (vertical).
;
;    /OTF to get the OTF at output.
;
;    /SF  to get the phase Structure Function at output.
;
;    /PSD to get the Phase Power Spectrum at output.
;
;    /ONLY_PSD same as above, but no structure function nor PSF or OTF are 
;    calculated,
;    just the PSDs. Use this option if you just want PAOLA to compute the PSD, and
;    nothing else. Incompatible with the keywords FWHM_ANALYSIS, EE ANALYSIS & co.,  
;    PSF, OTF, SF and PSD.
;
;    /WAVE to get the phase power spectral density inside the output
;    structure variable. This is useful to run the function wave.pro
;    to generate instantaneous corrected phase screens.
;
;OUTPUTS name | type | unit
;
;    The output is a structure variable, whose components are defined
;    as follow:
;
;    .lam | REAL SCALAR | MICROMETERS
;    Wavelength used in the calculation. A copy of input DIM.LAMBDA.
;
;    .inst | STRING | -
;    Type of focal instrumentation. 'IMAGER' or 'INTENSITY'. A copy of input TSC.INST.
;
;    .mode | STRING | -
;    Calculation mode.
;
;    .w0ZA | REAL SCALAR | METER
;    Seeing angle @ 500 nm for zenith angle ZA.
;
;    .r05 | REAL SCALAR | METER
;    Fried's parameter @ 500 nm for zenith angle ZA.
;
;    .r0l | REAL SCALAR | METER
;    Fried's parameter at the imaging wavelength for zenith angle ZA.
;
;    .L0  | REAL SCALAR | METER
;    Outer scale of turbulence. -1 if infinite.
;
;    .strehl | REAL SCALAR | -
;    >>> NOT IF /ONLY_PSD IS SET <<<
;    Strehl ratio.
;
;    .cfr | REAL SCALAR | -
;    >>> NOT IF /ONLY_PSD IS SET <<<
;    Ratio between the practical cut-off angular frequency in the focal plane, and
;    the cut-off frequency of the telescope without aberrations. Between 0 and 1.
;
;    .alt | REAL SCALAR | METER
;    ONLY IF HEIGHT IS SPECIFIED.
;    mean altitude of the turbulent layers for zenith angle ZA.
;
;    .ani | REAL SCALAR | ASEC
;    ONLY IF HEIGHT IS SPECIFIED.
;    anisoplanatism angle for zenith angle ZA.
;
;    .win | REAL SCALAR | METER/SECOND
;    ONLY IF WIND IS SPECIFIED.
;    mean wind speed over the layers (NOT adapted to zenith angle).
;
;    .lft | REAL SCALAR | MILLISECOND
;    ONLY IF WIND IS SPECIFIED.
;    life time of the speckle field for zenith angle ZA.
;
;    .c2h | REAL SCALAR | ASEC
;    On the PSF, radius where the diffraction limited core starts to be burden by
;    the residual non-corrected halo.
;
;    .var | REAL ARRAY | RAD^2
;    Variances associated with the different components of the residual phase power 
;    spectrum. These numbers are useful to identify the main sources of disturbances 
;    that cause the residual features on the final PSF.
;
;    |  .var[0] | Fitting error
;    |  .var[1] | Aliasing
;    |  .var[2] | Total aniso-servo and/or dispersion error.
;    |  .var[3] | WFS noise
;
;    .rms | REAL ARRAY | NM
;    Standard deviation associated with the different components of the residual wavefront 
;    spectrum. These numbers are useful to identify the main sources of disturbances 
;    that cause the residual features on the final PSF.
;
;    |  .rms[0] | Fitting error
;    |  .rms[1] | Aliasing
;    |  .rms[2] | Total aniso-servo and/or dispersion error.
;    |  .rms[3] | WFS noise
;
;    .dm_height | REAL SCALAR | METER
;    conjugation altitude of the deformable mirror.
;
;    .nna | REAL SCALAR | 1
;    nominal number of actuators across the deformable mirror diameter 
;    (corresponding to one act per r0-cell in the primary mirror).
;
;    .ana | REAL SCALAR | 1
;    actual number of actuators across the deformable mirror diameter. Equal to .nna 
;    if DM_PARAMS.ACTPITCH has been set to -1.
;
;    .tna | REAL SCALAR | 1
;    total number of actuators inside pupil = (M1 surface)/(act pitch in M1)^2. 
;
;    .act | REAL SCALAR | METER
;    Actuator pitch projected into the primary mirror plane.
;
;    .anl | REAL SCALAR | 1
;    actual number of WFS lenslets across the pupil diameter. Equal to .ana if 
;    WFS_PITCH has been set to -1.
;
;    .tnl | REAL SCALAR | 1
;    total number of WFS lenslets inside pupil = (M1 surface)/(WFS pitch in M1)^2. 
;
;    .wlp | REAL SCALAR | METER
;    WFS lenslet pitch projected into the primary mirror plane.
;
;    .ang | REAL SCALAR | ASEC
;    NGS off-axis angle.
;
;    .ori | REAL SCALAR | DEGREE
;    NGS off-axis orientation relative to focal plane x-axis.
;
;    .int | REAL SCALAR | MILLISEC
;    WFS integration time.
;
;    .lag | REAL SCALAR | MILLISEC
;    System's loop technical time lag.
;    
;    .servolag | REAL SCALAR | MILLISEC
;    servolag delay = technical lag (see LAG) + integration time
;    
;    .lfr | REAL SCALAR | 1
;    Loop frequency.
;    
;    .lpg | REAL SCALAR | 1
;    Loop gain. Only in closed loop mode.
;    
;    .stability | STRING | -
;    Indicates loop stability according to Nyquist criterion. Answer
;    is 'yes' or 'no'.
;
;    .mag | REAL SCALAR | 1
;    NGS magnitude for the specified filter.
;
;    .filter | STRING | -
;    Tag name of the filter for which the NGS magnitude is given.
;
;    .tem | REAL SCALAR | K
;    BB temperature associated with the NGS spectrum.
;
;    .ron | REAL SCALAR | e/px
;    WFS CCD read noise
;
;    .sbg | REAL SCALAR | ph/px/dt
;    Photon sky background per WFS pixel per integration time.
;
;    .wfslam | REAL SCALAR | MICRONS
;    Average wavelength of the guide star light distribution after atm+tsc+WFS 
;    optics filtering.
;
;    .wfsbdw | REAL SCALAR | MICRONS
;    RMS width of the guide star light distribution after atm+tsc+WFS optics
;    filtering.
;
;    .wfstau | REAL SCALAR | ]0,1]
;    Average optical transmission from the star to the WFS CCD inside the WFS
;    optical bandwidth (wfsbdw).
;
;    .nph | REAL SCALAR | ]0,infinite[
;    Number of photons/lenslet/frame.
;
;    .nea | REAL SCALAR | ASEC
;    WFS noise equivalent angle.
;
;OPTIONAL OUTPUTS name | type | unit
;
;    If you have requested the calculation of several metrics of the AO PSF,
;    these metrics will be added to the output structure variable as well. See
;    section << PSF METRICS CALCULATION >> for details.
;
;    If DISPERSION modeling has been used:
;
;      .dispersion | scalar real > 0 | 1 | 1
;      Dispersion factor, see input.
;
;    IF /SCINTILLATION IS SET
;
;      .sci_index | REAL SCALAR | 1
;      Scintillation index - variance of the relative irradiance fluctuations (4
;      times the variance of the relative amplitude fluctuations). See Ref [9].
;
;    IF /PSF OR /X_PSF OR /Y_PSF ARE SET
;
;      .dxf_usr | REAL SCALAR | ASEC/PX
;      PSF angular pixel scale.
;
;    IF /PSF is set
;
;      .psf | REAL ARRAY(N,N) | STREHL RATIO
;      PSF for the selected mode.
;
;    IF /X_PSF OR /Y_PSF are set
;
;      .psfx | REAL ARRAY(N) | STREHL RATIO
;      x-cut section of the PSF for the selected mode.
;
;      .psfy | REAL ARRAY(N) | STREHL RATIO
;      y-cut section of the PSF for the selected mode.
;
;    IF /OTF IS SET
;
;      .dff | REAL SCALAR | 1/RAD/PX | ALL
;      Spatial frequency pixel scale in the OTF matrix.
;
;      .otf | REAL ARRAY(N,N) | MAX OTF = 1
;      OTF of the selected mode.
;
;    IF /SF IS SET
;
;      .sf | REAL ARRAY(N,N) | RAD^2
;      Total low frequency plus high frequency structure function.
;
;      .lfsf | REAL ARRAY(N,N) | RAD^2
;      Low frequency phase structure function.
;
;      .hfsf | REAL ARRAY(N,N) | RAD^2
;      High frequency phase structure function.
;
;      .dxp | REAL SCALAR | METERS/PIXEL
;      Pixel scale in the pupil plane associated with the structure functions 
;      matrices.
;
;      .dfp | REAL SCALAR | 1/METERS/PX
;      Spatial frequency Pixel scale in the pupil, associated with the PSD matrices.
;
;    IF /PSD OR /ONLY_PSD AND /SCINTILLATION ARE SET
;
;      .psd_amp | REAL ARRAY(N,N) | RAD^2*M^2
;      Power spectrum of the relative fluctuation of the wave amplitude.
;
;    IF /PSD OR /ONLY_PSD ARE SET
;
;      .psd_fe | REAL 2D ARRAY | RAD^2*M^2
;      Fitting error power spectrum.
;    
;      .psd_al | REAL 2D ARRAY | RAD^2*M^2
;      WFS aliasing error power spectrum.
;      >>> NOT GIVEN IF /ANTI_ALIAS HAS BEEN SET <<<
;    
;      .psd_as | REAL 2D ARRAY | RAD^2*M^2
;      Aniso + servo error power spectrum.
;    
;      .psd_ns | REAL 2D ARRAY | RAD^2*M^2
;      WFS noise error power spectrum.
;
;      .psd_disp | REAL 2D ARRAY | RAD^2*M^2
;      Dispersion effect power spectrum.
;      >>> ONLY IF DISPERSION HAS BEEN SET <<<
;    
;      .dfp_lf | REAL SCALAR | M^(-1)
;      low spatial frequency domain (below 1/(2*pitch)) pixel size.
;
;      .dfp_hf | REAL SCALAR | M^(-1)
;      high spatial frequency domain (above 1/(2*pitch)) pixel size.
;
;    IF /PSD OR /ONLY_PSD AND /SCINTILLATION ARE SET
;
;      .psd_amp | REAL 2D ARRAY | RAD^2*M^2
;      Power spectrum of the relative fluctuation of the wave
;      amplitude.
;
;    IF /WAVE IS SET
;
;      .psd_wave | REAL ARRAY(N,N) | RAD^2*M^2
;      Residual phase spatial power spectrum, including in a single
;      matrix the low order and high order PSD.
;
;      .pup_wave | REAL ARRAY(N,N) | RAD^2*M^2
;      >>> ONLY IF /PUPIL WAS SET IN THE CALL TO psfotftsc.pro <<<
;      Pupil mask.
;
;      .dfp_wave | REAL SCALAR | 1/METERS/PX
;      Spatial frequency pixel scale in the pupil, associated with the
;      psd_wave matrix.
;
;      .dxp_wave | REAL SCALAR | METERS/PX
;      Spatial pixel scale in the pupil.
;
;    IF VIBRATION IS SET
;
;      .VIBRATION : a copy of optional input VIBRATION.
;
;###################################################################################
;###################################################################################
;############################                               ########################
;############################      GROUND LAYER AO MODE     ########################
;############################                               ########################
;###################################################################################
;###################################################################################
;
;CALLING SEQUENCE
;
;    USAGE 1:
;
;      We assume that we have a perfect knowledge of the turbulence above 
;      the telescope in all directions, and that the DM command is defined by the 
;      average of the phase over all possible directions on the WFS FoV. It is 
;      equivalent to assume that we have an infinite number of natural guide stars 
;      spread over the system FoV. Is it realistic ? not really, but at least it 
;      gives an upper bound on what we can get in GLAO correction. Historically, 
;      this ideal wavefront sensing scheme was the first studied in the early 
;      analysis of GLAO systems performance (Rigaut, 2002).
;
;      GLAO_WFS.TYPE='full'
;
;      res=PAOLA('glao',DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,
;                DM_PARAMS,WFS_PARAMS,SO_ANG,SO_ORI,GLAO_WFS,[options])
;
;    USAGE 2:
;
;      We assume that we know the phase in all the directions along a circle, as
;      if there was an infinite number of NGS aligned on a circle (A. Tokovinin
;      proposed this mode in 2003). The DM command is 
;      the average of the phase measured from each direction along the circle. In 
;      such a mode, the correction is made more homogeneous across the FoV, but less 
;      effective than in "full" mode in the center of the FoV. This mode is a bit 
;      more realistic though, as we can imagine to implement an annulus of N>>1 LGS, 
;      or rotate a single LGS very fast in the sky (M. Britton idea). Note that LGS 
;      cone effect is not taken into account, so we are considering only NGS for 
;      now.
;
;      GLAO_WFS.TYPE='edge'
;
;      res=PAOLA('glao',DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,$
;                DM_PARAMS,WFS_PARAMS,SO_ANG,SO_ORI,GLAO_WFS,[options])
;
;    USAGE 3: REALISTIC
;
;      Here we have the realistic case, where the phase is measured only on a couple 
;      of individual NGS spread across the field. GS_WEIGHT input is used to give a
;      weight to the different GS, depending for instance on the distance to the 
;      science object, or the GS magnitude. We take into account the
;      correlation of the aniso and servo errors (aniso-servo), so the WIND
;      input is needed, the WFS aliasing, and the WFS noise. So this is the most
;      realistic case of NGS-based GLAO. The NGS magnitude, for a specified
;      filter and BB temperature, the WFS light transmission and detection
;      characteristics must be given, or instead, the WFS noise equivalent angle
;      NEA can be given.
;
;      GLAO_WFS.TYPE='star'
;
;      res=PAOLA('glao',DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,WIND,$
;                DM_PARAMS,WFS_PARAMS,SO_ANG,SO_ORI,WFS_INT,LAG,LOOP_MODE,$
;                LOOP_GAIN,GLAO_WFS,GS_WEIGHT,NGS_MAG,FILTER,NGS_TEM,[options])
;
;      Alternatively, the WFS noise equivalent angle (lenslet slope noise) can be
;      given instead of the NGS and WFS parameters:
;
;      res=PAOLA('glao',DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,WIND,$
;                DM_PARAMS,WFS_PARAMS,SO_ANG,SO_ORI,WFS_INT,LAG,LOOP_MODE,$
;                LOOP_GAIN,GLAO_WFS,GS_WEIGHT,WFS_NEA,[options])
;
;    For examples, see the section GLAO in batch.pro
;
;INPUTS name | type | unit
;
;    MODE | STRING | -
;    AO mode. Here it must be set to 'glao'.
;
;    DIM | structure variable | -
;    DIM is a structure variable where the pixels and matrices sizes are
;    kept. See output of function PIXMATSIZE.PRO
;
;    TSC | structure variable | -
;    TSC is a structure variable including the telescope optical parameters and
;    OTF/PSF/PHASE/PUPIL. See output of function PSFOTFTSC.PRO
;
;    W0 | REAL SCALAR | ARCSEC
;    Median SEEING angle @ 500 nm, at zenith. This value is used to compute the
;    Fried parameter r0 for the specified wavelength. Zenith value must be given.
;
;    L0 | REAL SCALAR | METRE
;    Outer scale of optical turbulence. To get an infinite outer scale, set L0 
;    to -1. Note that an outer scale larger than 10 km is considered infinite.
;
;    ZA | REAL SCALAR | DEGREE
;    Zenith angle of the observation. The apparent seeing angle and turbulent
;    layers altitude change with the zenith angle and are modified
;    accordingly within the code. ZA must be in the range [0,90[ degrees,
;    i.e. the 90 deg value is excluded (it gives an infinite seeing angle).
;
;    HEIGHT | REAL ARRAY(number of layers) OR SCALAR | METRE
;    altitude of each of the NL considered turbulent layers, and this is
;    >>>>> RELATIVE TO THE TELESCOPE PUPIL <<<<<<<. Note that it can be a scalar 
;    (for a simplified 1-layer profile, for instance ground layer).
;
;    DISTCN2 | REAL ARRAY(number of layers) OR SCALAR | 1
;    relative distribution of the Cn^2*dh profile over the NL turbulent layers (dh 
;    is the 
;    layer thickness). It does not have to be normalized, i.e total(DISTCN2) can be 
;    different that 1, and the Cn2*dh profile itself can be given here. The actual 
;    amplitude of the Cn2*dh profile will be set inside PAOLA according to the 
;    provided seeing angle W0. Can be a scalar (1-layer profile).
;
;    WIND | REAL ARRAY(number of layers,2) | METRE/SEC
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    WIND[*,0]: x-axis component (east->west dir/sky). 
;    WIND[*,1]: y-axis component (south->north dir/sky).
;    Must be introduced that way: [[x components],[y components]]
;    example for 3 layers : [[10,0,5],[0,-1,-5]]. Can be given as a
;    scalar, too. Then we will consider that WIND corresponds to the
;    x-coordinate of the wind velocity of a single layer (1-layer
;    profile case).
;    NOTE: WIND VELOCITY VECTOR IS NOT ADAPTED TO THE ZENITH DISTANCE
;    BECAUSE THIS ADAPTATION WOULD DEPENDS ON THE POINTING POSITION ON THE
;    SKY, AND THAT WOULD MAKE EVERYTHING UNNECESSARY COMPLICATED. IF YOU
;    WANT TO STUDY THIS EFFECT, JUST ADAPT THE WIND VECTOR COEFFICIENTS TO
;    YOUR PARTICULAR POINTING POSITION IN THE SKY, WHEN YOU DECLARE YOUR
;    WIND VECTOR (OUTSIDE OF PAOLA, IN YOUR BATCH).
;
;    DM_PARAMS | a 3 components structure variable | -
;    DM parameters structure variable. Check IDL user manual to see how to create a 
;    structure variable. See also example in batch.pro file. The components are:
;    |
;    |__DM_PARAMS.DM_HEIGHT | REAL SCALAR | METERS
;    |  Conjugated altitude of the DM with respect to the ground level.
;    |
;    |__DM_PARAMS.DMTF | REAL ARRAY(N,N) OR -1 | 1
;    |  DM **AMPLITUDE** spatial transfer function. To get a perfect DM filter, set 
;    |  this input to -1.
;    |  See function DM_MODEL.PRO user manual for instructions to build a DM
;    |  transfer function using a given set of influence functions.
;    |
;    |__DM_PARAMS.ACTPITCH | REAL SCALAR OR -1 | M
;    |  DM actuator pitch as seen in the pupil plane. To get the default value (r0 
;    |  at the imaging wavelength and at zenith), set this input to -1.
;    |
;    |  An important limitation: the transition between the diffraction limited core 
;    |  of the PSF and the residual atmospheric halo occurs at an off-axis angle 
;    |  lambda/(2*actpitch). Therefore, it is important that you choose
;    |  a PSF field-of-view (see PIXMATSIZE inputs) large enough to see this
;    |  transition. Otherwise, the Strehl performance will be overestimated. In my 
;    |  opinion a FoV radius at least twice the transition angle (radius), i.e. a 
;    |  total FoV width of 4*lambda/(2*actpitch) is really a minimum. Mama PAOLA 
;    |  will warn you if this is not the case, and the calculation will not start.
;    |
;    |  There are three possible options for input DM_PARAMS:
;    |
;    |  1st option: DM_PARAMS.DMTF = -1
;    |              DM_PARAMS.ACTPITCH = -1
;    |
;    |  Here, we use a perfect DM filtering of the phase aberration below the DM 
;    |  cut-off frequency 1/(2*ACTPITCH), i.e. the DM transfer function is = 1 for
;    |  |fx|<f_c and |fy|<f_c and 0 outside of the square domain. Such a model is 
;    |  fairly accurate and works for most applications. Actuator pitch is set to 
;    |  the default value (see above).
;    |
;    |  2nd option: DM_PARAMS.DMTF = -1
;    |              DM_PARAMS.ACTPITCH = a positive value.
;    |
;    |  The DM is still an ideal one, and the actuator pitch is set to a specific 
;    |  value.
;    |
;    |  3rd option: DM_PARAMS.DMTF = a 2D DM transfer function model
;    |              DM_PARAMS.ACTPITCH = a positive value.
;    |
;    |  Here, a 2 dimensional realistic model of a DM transfer function is given, 
;    |  along 
;    |  with the corresponding actuator pitch. See usage of function DM_MODEL.PRO, 
;    |  that you should use to build the DM transfer function. DM_PARAMS.DMTF is 
;    |  dimensionless (it's a transfer function).
;    |
;    |  Note that for most AO correction simulation cases, the perfect DM model is 
;    |  good enough and does not produce PSF that are significantly different than 
;    |  PSF produced using a more sophisticated DM filter. ExAO modeling might be 
;    |  the only case where using sophisticated DM modeling matters (see usage of 
;    |  DM_MODEL.PRO function for more details on how to build such a DM transfer 
;    |__function).
;
;    WFS_PARAMS | a structure variable | -
;    WFS parameters structure variable. Check IDL user manual to see how to create a 
;    structure variable. See also example in batch.pro file. The components are:
;    |
;    |__WFS_PARAMS.WFS_PITCH | REAL SCALAR | METERS
;    |  WFS lenslet array pitch size, as seen from the entrance pupil plane. If set 
;    |  to -1, the WFS lenslet pitch is assumed to be equal to the DM actuator 
;    |  pitch, i.e. WFS_PITCH is set to DM_PARAMS.ACTPITCH within the code.
;    |
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    |__IF THE WFS NOISE IS GIVEN AS A NOISE EQUIVALENT ANGLE WITH THE INPUT
;       WFS_NEA, ONLY SET WFS_PARAMS.WFS_PITCH AND IGNORE THE OTHER
;       COMPONENTS BELOW, AS THEY ARE USELESS IN THIS CASE.
;       |
;       |__WFS_PARAMS.MIRCOATING | VECTOR OF STRING | [-]
;       |  This input set the number and type of mirrors from (and including) the
;       |  telescope primary mirror until the WFS lenslet array. It is used in the
;       |  calculation of the number of photons received by the WFS detector /
;       |  lenslet / frame.
;       |  There are three possible mirror coating: Aluminum, Silver, and Gold.
;       |  Set MIRCOATING to the suite of coatings along the optical train: 
;       |  For instance if your telescope is Al-Al-Ag then Al-Al-Au inside,
;       |  WFS_PARAMS.MIRCOATING=['Al','Al','Ag','Al','Al','Au']
;       |  but if you want to use the default case, which is 3 Aluminum mirrors & 3 Silver
;       |  mirrors with normal incidence, simply set MIRCOATING to -1.
;       |
;       |__WFS_PARAMS.MIRANGLE | VECTOR OF REAL | DEGREE
;       |  Beam incidence angle of the mirrors introduced above. For instance
;       |  WFS_PARAMS.MIRANGLE=[0,0,45,30,30,45]
;       |
;       |__WFS_PARAMS.NBLENSES | INTEGER SCALAR | [0,N]
;       |  Number of REFRActive surfaces from the entrance pupil (included) to the
;       |  detector.
;       |
;       |__WFS_PARAMS.WFS_RON | REAL SCALAR | ELECTRON/PIXEL
;       |  Pixel read noise of the WFS CCD detector. If WIND has been set to 0,  
;       |  WFS_RON is ignored.
;       |
;       |__WFS_PARAMS.SKY_RADIANCE | REAL SCALAR | W/M^2/SR
;       |  Sky background radiance in watt / square meter / steradian
;       |
;       |__WFS_PARAMS.ALGORITHM | CHAR STRING | -
;       |  Type of algorithm for SH spot centroid calculation:
;       |  Enter '4Q' for a 4-quadrant mode;
;       |  Enter 'CG' for a center-of-gravity mode.
;       |  Enter 'GS' for a center-of-gravity mode and Gaussian spot
;       |             whose FWHM is set by residual turbulence.
;       |
;       |__WFS_PARAMS.EXTRAFILTER - a structure variable
;       |  If, on top of the transmission losses due to atmosphere and optical
;       |  surfaces, there is in your system an other transmission loss that cannot 
;       |  be modeled with the inputs above, you can use the EXTRAFILTER input.
;       |  EXTRAFILTER is a structure variable, with two components:
;       |  .TAUEXT | REAL VECTOR | 1
;       |  The filter transmission, in the range [0,1] (i.e NOT in %)
;       |  .LAMBDA | REAL VECTOR | MICROMETERS
;       |  The wavelength associated with the filter data. Should not be
;       |  outside of the range 300 - 5000 nm.
;       |  The number of data points in TAUEXT and LAMBDA does not matter:
;       |  the code will either fit a spline to the data or rebin the data
;       |  to adjust to the internal wavelength sampling, which is 10 nm.
;       |  Syntax inside the setting of WFS_PARAMS:
;       |  WFS_PARAMS.EXTRAFILTER:{TAUEXT:some vector,LAMBDA:some vector}
;       |
;       |  -----> IF YOU DO NOT NEED THIS INPUT, SET IT TO 'no'
;       |         WFS_PARAMS.EXTRAFILTER:'no'
;       |
;       |__THE TWO INPUTS BELOW ARE ONLY TO BE SET IF YOU HAVE CHOSEN A
;          CENTER-OF-GRAVITY ALGORITHM FOR THE CENTROID MEASUREMENT. SKIP THESE IF
;          YOU HAVE A 4-QUADRANT ALGORITHM.
;          |
;          |__WFS_PARAMS.WFS_PXFOV | SCALAR INTEGER | PIXELS
;          |  Field of view for each lenslets, in pixels. Must be >= 6, because with 
;          |  less than that, non-linearities due to pixelization become an issue.
;          |
;          |__WFS_PARAMS.WFS_PXSIZE | SCALAR OR -1 | ASEC/PIXEL
;             WFS detector pixel size (projected on sky). Must be >= Nyquist limit
;             lambda/(2*WFS_PITCH), expressed in asec. Set WFS_PXSIZE to -1 to
;             get the Nyquist limit, which is the best choice anyway.
;
;    SO_ANG | REAL SCALAR | ASEC
;    science object off-axis angle. If HEIGHT has been set to 0, SO_ANG is ignored.
;
;    SO_ORI | REAL SCALAR | DEGREES
;    orientation of the science object position vector in the field of view, where 0 
;    correspond to an alignment with the +x axis, and +90 with the +y axis.
;    If HEIGHT has been set to 0, SO_ORI is ignored.
;
;    WFS_INT | REAL SCALAR | MILLISECOND
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    integration time of the wave-front sensor. MUST be > 0.
;    NOTE If WFS_INT optimization is required, to balance the servo-lag and the noise
;    errors, set the keyword /OPTIMIZE_WFS_INT. If you already have an initial guess
;    for the optimal WFS integration time, set WFS_INT to this initial value, and the
;    optimization will start with this. In this case, the user must provide the NGS
;    magnitude, black body temperature, and WFS CCD read noise.
;
;    LAG | REAL SCALAR | MILLISECOND
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    servo system technical time lag  = read-out time of the WFS + computation
;    time. In other words, this is the time lapse between the end of the
;    WFS exposure and the application of the updated command to the DM. 
;
;      CAUTION: the TOTAL DELAY between the WFS measurement and the DM
;      correction is equal to the WFS exposure time, or the loop
;      period, plus the technical LAG defined above. So, when entering
;      the LAG parameter, make sure that it is the technical delay,
;      and not the total delay. Inside PAOLA, the total delay is
;      computed by adding together the WFS integration time and LAG.
;
;    LOOP_MODE | CHAR STRING | -
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    Loop mode. 'OPEN' for now.
;
;    LOOP_GAIN | REAL SCALAR | 1
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    Loop gain. Ignored for now. Set it to 1.
;
;    GLAO_WFS: a structure variable, with components:
;
;    .TYPE | CHAR STRING | -
;    type of WFS measurement:
;    'full': average over the full WFS system FoV,
;    'edge': average over the edge of the WFS FoV,
;    'star': average over several guide stars in different directions
;            in the WFS FoV.
;
;    .ANG | REAL SCALAR OR ARRAY(2,number of guide star) | ARCSEC
;    depends on the type of WFS measurement:
;    'full': .ANG is the WFS FoV
;    'edge': .ANG is the WFS FoV
;    'star': .ANG[0,*] = guide stars angular x-coordinate (East->West)
;            .ANG[1,*] = guide stars angular y-coordinate (South->North)
;
;    GS_WEIGHT | REAL ARRAY (Number of NGS) | 1
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    Weight to be associated to each guide star, according to some
;    user defined optimization scheme. Does not have to be normalized
;    to 1. This is done inside the code by default. If you want to
;    give the same weight to each star (a good start), set this input
;    to -1.
;
;    The WFS noise error can be calculated either from the NGS properties, or
;    from a pre-defined Noise Equivalent Angle input. In the first case, use the
;    following four inputs:
;
;      NGS_MAG | REAL VECTOR (NB of NGS) | 1
;      ##################################
;      **** ONLY IN GLAO 'star' MODE ****
;      ##################################
;      >>> Incompatible with WFS_NEA <<<
;      Magnitude of the NGS, for a given filter (see FILTER input below). We use
;      here the Johnson filter system, and the Vega zero point magnitude
;      reference values.
;
;      FILTER | CHARACTER STRING VECTOR (NB of NGS) | -
;      ##################################
;      **** ONLY IN GLAO 'star' MODE ****
;      ##################################
;      >>> Incompatible with WFS_NEA <<<
;      Tag name of the filter associated to the magnitude given in
;      NGS_MAG. Possible values are (Johnson system): U B V R I J H K L M N Q.
;      Syntax: do not forget to enter the tag name as a string, for instance 'V'.
;
;      NGS_TEM | REAL VECTOR (NB of NGS) | KELVIN
;      ##################################
;      **** ONLY IN GLAO 'star' MODE ****
;      ##################################
;      >>> Incompatible with WFS_NEA <<<
;      Black-body temperature associated to the NGS spectrum.
;
;    and in the second case, use the following input:
;
;      WFS_NEA | REAL VECTOR (NB NGS) | ASEC
;      ##################################
;      **** ONLY IN GLAO 'star' MODE ****
;      ##################################
;      >>> Incompatible with NGS_MAG,FILTER,NGS_TEM <<<
;      >>> Incompatible with WFS_INT OPTIMIZATION <<<
;      WFS Noise Equivalent Angle (same for each lenslet).
;
;OPTIONAL INPUT
;
;    VIBRATION= | A 3 COMPONENTS STRUCTURE VARIABLE | - | -
;    |  Vibrations induced on the beam, as seen in the image plane. Compatible
;    |  with the vibration of the telescope beam. This input is useful to
;    |  simulate a jitter occuring in the science arm, after AO correction,
;    |  independant from the telescope jitter. So you may have separate telescope
;    |  and post-AO instrument vibrations. This is useful for coronagraphs analysis.
;    |
;    |__.SIGMA_TTX | POSITIVE REAL SCALAR | ASEC
;    |  Tip jitter RMS, as seen in the focal plane (sky).
;    |
;    |__.SIGMA_TTY | POSITIVE REAL SCALAR | ASEC
;    |  Tilt jitter RMS, as seen in the focal plane (sky). 
;    |
;    |__.ORIENTATION | REAL SCALAR | DEGREE
;    |  Orientation angle of the tip (TTx) perturbation. Counterclockwise.
;    |  +90 corresponds to the positive y-axis, +180 to the negative
;    |  x-axis. The tilt jitter is perpendicular to this direction. 
;    |
;    |--SYNTAX : VIBRATION={SIGMA_TTX:value,SIGMA_TTY:value,ORIENTATION:angle}
;
;KEYWORDS
;
;    /RADDMTF To get a radial only DM transfer function at the AO radius 1/2/pitch
;
;    /ANTI_ALIAS if set, simulate spatial filtering of WFS aliasing by not adding 
;    the aliased component to the total power spectrum of the corrected phase.
;
;    OPTIMIZATION OF WFS INTEGRATION TIME
;    ====================================
;
;      IF THE WFS NOISE IS ALREADY FIXED WITH THE INPUT WFS_NEA, OR 
;      IF THE WIND SPEED IS NULL, OPTIMIZATION MAKES NO SENSE AND PAOLA
;      RETURNS AN ERROR MESSAGE.
;
;      /OPTIMIZE_WFS_INT WFS integration time optimization: PAOLA seeks the optimal 
;      WFS integration time to minimize the combined effect of WFS noise (NGS photon 
;      noise and CCD read noise), which decreases with WFS_INT and the servo-lag 
;      error, which increases with WFS_INT.
;
;    /SCINTILLATION to take into account wave amplitude fluctuation. Works only if a 
;    Cn2 profile is provided. See HEIGHT and DISTCN2 optional inputs above.
;
;    /INFO print the parameters and results of the calculation.
;
;    /PSF to get the PSF at output.
;
;    /X_PSF If set, only the section of the PSF along the x-axis (horizontal) is
;    calculated and given in the output structure variable. Saves a lot of 
;    computation time.
;
;    /Y_PSF same as above, for y-axis (vertical).
;
;    /OTF to get the OTF at output.
;
;    /SF  to get the phase Structure Function at output.
;
;    /PSD to get the Phase Power Spectrum at output.
;
;    /ONLY_PSD same as above, but no structure function nor PSF or OTF are 
;    calculated,
;    just the PSDs. Use this option if you just want PAOLA to compute the PSD, and
;    nothing else. Incompatible with the keywords FWHM_ANALYSIS, EE ANALYSIS & co.,  
;    PSF, OTF, SF and PSD.
;
;    /WAVE to get the phase power spectral density inside the output
;    structure variable. This is useful to run the function wave.pro
;    to generate instantaneous corrected phase screens.
;
;OUTPUTS name | type | unit
;
;    The output is a structure variable, whose components are defined as follow:
;
;    .lam | REAL SCALAR | MICROMETERS
;    Wavelength used in the calculation. A copy of input DIM.LAMBDA.
;
;    .inst | STRING | -
;    Type of focal instrumentation. 'IMAGER' or 'INTENSITY'. A copy of input TSC.INST.
;
;    .mode | STRING | -
;    Calculation mode.
;
;    .w0ZA | REAL SCALAR | METER
;    Seeing angle @ 500 nm for zenith angle ZA.
;
;    .r05 | REAL SCALAR | METER
;    Fried's parameter @ 500 nm for zenith angle ZA.
;
;    .r0l | REAL SCALAR | METER
;    Fried's parameter at the imaging wavelength for zenith angle ZA.
;
;    .L0  | REAL SCALAR | METER
;    Outer scale of turbulence. -1 if infinite.
;
;    .strehl | REAL SCALAR | -
;    >>> NOT IF /ONLY_PSD IS SET <<<
;    Strehl ratio.
;
;    .cfr | REAL SCALAR | -
;    >>> NOT IF /ONLY_PSD IS SET <<<
;    Ratio between the practical cut-off angular frequency in the focal plane, and 
;    the cut-off frequency of the telescope without aberrations. Between 0 and 1.
;
;    .alt | REAL SCALAR | METER
;    ONLY IF HEIGHT IS SPECIFIED.
;    mean altitude of the turbulent layers for zenith angle ZA.
;
;    .ani | REAL SCALAR | ASEC
;    ONLY IF HEIGHT IS SPECIFIED.
;    anisoplanatism angle for zenith angle ZA.
;
;    .win | REAL SCALAR | METER/SECOND
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    mean wind speed over the layers (NOT adapted to zenith angle ZA).
;
;    .lft | REAL SCALAR | MILLISECOND
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    life time of the speckle field for zenith angle ZA.
;
;    .c2h | REAL SCALAR | ASEC
;    On the PSF, radius where the diffraction limited core starts to be burden by
;    the residual non-corrected halo.
;
;    .var | REAL ARRAY | RAD^2
;    Variances associated with the different components of the residual phase power 
;    spectrum. These numbers are useful to identify the main sources of disturbances 
;    that generate residual features on the PSF.
;
;       .var[0] | Fitting error
;       .var[1] | Aliasing
;       .var[2] | Total aniso-servo
;       .var[3] | WFS noise
;
;    .rms | REAL ARRAY | NM
;    Standard deviation associated with the different components of the residual wavefront 
;    spectrum. These numbers are useful to identify the main sources of disturbances 
;    that cause the residual features on the final PSF.
;
;    |  .rms[0] | Fitting error
;    |  .rms[1] | Aliasing
;    |  .rms[2] | Total aniso-servo and/or dispersion error.
;    |  .rms[3] | WFS noise
;
;    .dm_height | REAL SCALAR | METER
;    conjugation altitude of the deformable mirror.
;
;    .nna | REAL SCALAR | 1
;    nominal number of actuators across the deformable mirror diameter 
;    (corresponding to one act per r0-cell in the primary mirror).
;
;    .ana | REAL SCALAR | 1
;    actual number of actuators across the deformable mirror diameter. Equal to .nna 
;    if DM_PARAMS.ACTPITCH has been set to -1.
;
;    .tna | REAL SCALAR | 1
;    total number of actuators inside pupil = (M1 surface)/(act pitch in M1)^2. 
;
;    .act | REAL SCALAR | METER
;    Actuator pitch projected into the primary mirror plane.
;
;    .anl | REAL SCALAR | 1
;    actual number of WFS lenslets across the pupil diameter. Equal to .ana if 
;    WFS_PITCH has been set to -1.
;
;    .tnl | REAL SCALAR | 1
;    total number of WFS lenslets inside pupil = (M1 surface)/(WFS pitch in M1)^2. 
;
;    .wlp | REAL SCALAR | METER
;    WFS lenslet pitch projected into the primary mirror plane.
;
;    .ang | REAL SCALAR | ASEC
;    Science object off-axis angle.
;
;    .ori | REAL SCALAR | DEGREE
;    Science object off-axis orientation relative to focal plane x-axis.
;
;    .int | REAL SCALAR | MILLISECOND
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    GS WFS integration time.
;
;    .lag | REAL SCALAR | MILLISECOND
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    System's technical time lag.
;
;    .servolag | REAL SCALAR | MILLISEC
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    servolag delay = technical lag (see LAG) + integration time
;    
;    .gs_weight | REAL ARRAY (Number of NGS) | 1
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    Weight associated to each guide star.
;
;    .glao_wfs: a structure variable, with components:
;
;      .glao_wfs.type | CHAR STRING | -
;      type of WFS measurement:
;      'full': average over the full WFS system FoV,
;      'edge': average over the edge of the WFS FoV,
;      'star': average over several guide stars in different directions
;              in the WFS FoV.
;
;      .glao_wfs.ang | REAL SCALAR OR ARRAY(2,number of guide star) | ARCSEC
;      depends on the type of WFS measurement:
;      'full': .ang is the WFS FoV
;      'edge': .ang is the WFS FoV
;      'star': .ang[1,*] = guide stars angular x-coordinate (East-West)
;              .ang(2,*) = guide stars angular y-coordinate (North-South)
;
;    .mag | VECTOR REAL | 1
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    NGS magnitudes in the given filter.
;
;    .filter | VECTOR STRING | -
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    Filters tags associated to the magnitudes of the NGS.
;
;    .tem | VECTOR REAL | K
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    BB temperatures associated to the NGS.
;
;    .ron | REAL SCALAR | e/px
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    WFS CCD read noise.
;
;    .sbg | REAL SCALAR | ph/px/dt
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    Photon sky background per WFS pixel per integration time.
;
;    .wfslam | VECTOR REAL | MICRONS
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    Average wavelengths of the NGS light distributions after atm+tsc+WFS optics
;    filtering.
;
;    .wfsbdw | VECTOR REAL | MICRONS
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    RMS width of the guide star light distribution after atm+tsc+WFS optics
;    filtering.
;
;    .wfstau | VECTOR REAL | ]0,1]
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    Average optical transmission from the star to the WFS CCD inside the WFS
;    optical bandwidth (wfsbdw).
;
;    .nph | VECTOR REAL | ]0,infinite[
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    Number of photons/lenslet/frame.
;
;    .nea | VECTOR REAL | ASEC
;    ##################################
;    **** ONLY IN GLAO 'star' MODE ****
;    ##################################
;    WFS noise equivalent angle.
;
;OPTIONAL OUTPUTS name | type | unit
;
;    If you have requested the calculation of several metrics of the AO PSF,
;    these metrics will be added to the output structure variable as well. See
;    section << PSF METRICS CALCULATION >> for details.
;
;    IF /SCINTILLATION IS SET
;
;    .sci_index | REAL SCALAR | 1
;    Scintillation index - variance of the relative irradiance fluctuations (4 times 
;    the variance of the relative amplitude fluctuations). See Ref [9].
;
;    IF /PSF OR /X_PSF OR /Y_PSF ARE SET
;
;    .dxf_usr | REAL SCALAR | ASEC/PX
;    PSF angular pixel scale.
;
;    IF /PSF is set
;
;    .psf | REAL ARRAY(N,N) | STREHL RATIO
;    PSF for the selected mode.
;
;    IF /X_PSF OR /Y_PSF are set
;
;    .psfx | REAL ARRAY(N) | STREHL RATIO
;    x-cut section of the PSF for the selected mode.
;
;    .psfy | REAL ARRAY(N) | STREHL RATIO
;    y-cut section of the PSF for the selected mode.
;
;    IF /OTF IS SET
;
;    .dff | REAL SCALAR | 1/RAD/PX | ALL
;    Spatial frequency pixel scale in the OTF matrix.
;
;    .otf | REAL ARRAY(N,N) | MAX OTF = 1
;    OTF of the selected mode.
;
;    IF /SF IS SET
;
;    .sf | REAL ARRAY(N,N) | RAD^2
;    Total low frequency plus high frequency structure function.
;
;    .lfsf | REAL ARRAY(N,N) | RAD^2
;    Low frequency phase structure function.
;
;    .hfsf | REAL ARRAY(N,N) | RAD^2
;    High frequency phase structure function.
;
;    .dxp | REAL SCALAR | METERS/PIXEL
;    Pixel scale in the pupil plane associated with the structure functions 
;    matrices.
;
;    .dfp | REAL SCALAR | 1/METERS/PX
;    Spatial frequency Pixel scale in the pupil, associated with
;    the PSD matrices.
;
;    IF /PSD OR /ONLY_PSD ARE SET
;
;    .psd_fe | REAL 2D ARRAY | RAD^2*M^2
;    Fitting error power spectrum.
;    
;    .psd_al | REAL 2D ARRAY | RAD^2*M^2
;    WFS aliasing error power spectrum.
;    
;    .psd_as | REAL 2D ARRAY | RAD^2*M^2
;    Aniso + servo error power spectrum.
;    
;    .psd_ns | REAL 2D ARRAY | RAD^2*M^2
;    WFS noise error power spectrum.
;
;    .dfp_lf | REAL SCALAR | M^(-1)
;    low spatial frequency domain (below 1/(2*pitch)) pixel size.
;
;    .dfp_hf | REAL SCALAR | M^(-1)
;    high spatial frequency domain (above 1/(2*pitch)) pixel size.
;
;    IF /PSD OR /ONLY_PSD AND /SCINTILLATION ARE SET
;
;    .psd_amp | REAL 2D ARRAY | RAD^2*M^2
;    Power spectrum of the relative fluctuation of the wave
;    amplitude.
;
;    IF /WAVE IS SET
;
;      .psd_wave | REAL ARRAY(N,N) | RAD^2*M^2
;      Residual phase spatial power spectrum, including in a single
;      matrix the low order and high order PSD.
;
;      .pup_wave | REAL ARRAY(N,N) | RAD^2*M^2
;      >>> ONLY IF /PUPIL WAS SET IN THE CALL TO psfotftsc.pro <<<
;      Pupil mask.
;
;      .dfp_wave | REAL SCALAR | 1/METERS/PX
;      Spatial frequency pixel scale in the pupil, associated with the
;      psd_wave matrix.
;
;      .dxp_wave | REAL SCALAR | METERS/PX
;      Spatial pixel scale in the pupil.
;
;    IF VIBRATION IS SET
;
;      .VIBRATION : a copy of optional input VIBRATION.
;
;###################################################################################
;###################################################################################
;############################                               ########################
;############################    PSF METRICS CALCULATION    ########################
;############################                               ########################
;###################################################################################
;###################################################################################
;
;    Several options are available to measure the characteristic of the PSF:
;
;    ==========================
;    FULL WIDTH at HALF MAXIMUM
;    ==========================
;
;        An ellipse is fitted to a section of the PSF at half maximum, and the 
;        maximum and minimum width of the ellipse are given, with the orientation of 
;        the PSF elongation w.r.t. the horizontal axis.
;
;      INPUTS
;
;        Set the keyword /FWHM_ANALYSIS. Not compatible with keywords /X_PSF & 
;        /Y_PSF.
;
;      OUTPUTS
;
;        .fmi | REAL SCALAR | ASEC
;        Minor axis of an ellipse fitted to the PSF at half maximum
;        level. -1 is returned if the PSF is so large that the FWHM does not
;        fit into the dedicated matrix (FWHM estimate failure).
;
;        .fma | REAL SCALAR | ASEC
;        Major axis of an ellipse fitted to the PSF at half maximum
;        level. -1 is returned if the PSF is so large that the FWHM does not
;        fit into the dedicated matrix (FWHM estimate failure).
;
;        .dir | REAL SCALAR | DEG
;        Major axis orientation of an ellipse fitted to the PSF at half maximum 
;        level. Between -90deg and +90deg. If >0 the PSF ellipse points towards the 
;        1st quadrant. -1 is returned if the PSF is so large that the FWHM does not 
;        fit into the dedicated matrix (FWHM estimate failure).
;
;    ============================================
;    50% and 80% INTEGRATED ENERGY APERTURE WIDTH
;    ============================================
;
;        The aperture can be circular, square, or a slit. We give the disc diameter, 
;        or the square width, or the slit width width corresponding to 50% and 80% 
;        of the total intensity in the PSF.
;
;        Not compatible with keywords /X_PSF & /Y_PSF.
;
;      INPUTS
;
;        Set the keywords /EE_ANALYSIS_DISC or
;                         /EE_ANALYSIS_SQUA or
;                         /EE_ANALYSIS_SLIT
;
;      OUTPUTS - depending on the keyword chosen,
;
;        .e50disc | REAL SCALAR | ASEC
;        Angular diameter on the PSF including 50% of the total energy.
;
;        .e80disc | REAL SCALAR | ASEC
;        Angular diameter on the PSF including 80% of the total energy.
;
;        .e50squa | REAL SCALAR | ASEC
;        Square width angle on the PSF including 50% of the total energy.
;
;        .e80squa | REAL SCALAR | ASEC
;        Square width angle on the PSF including 80% of the total energy.
;
;        .e50slit | REAL SCALAR | ASEC
;        Angular slit width on the PSF including 50% of the total energy.
;    
;        .e80slit | REAL SCALAR | ASEC
;        Angular slit width on the PSF including 80% of the total energy.
;
;    ======================================================================
;    APERTURE WIDTH ASSOCIATED TO A CERTAIN PROPORTION OF INTEGRATED ENERGY
;    ======================================================================
;
;        Here the proportion of energy within a given disc diameter, square
;        width or slit width is calculated.
;
;        Not compatible with keywords /X_PSF & /Y_PSF.
;
;      INPUTS
;
;        EEW_DISC= | REAL | ARCSEC | -
;        EEW_SQUA= | REAL | ARCSEC | -
;        EEW_SLIT= | REAL | ARCSEC | -
;        Computes the proportion of energy within a disc of diameter EEW_DISC, or a 
;        square of width EEW_SQUA, or a slit of width EEW_SLIT.
;
;      OUTPUTS - depending on the chosen input chosen,
;
;        .eew | REAL SCALAR | ASEC
;        Aperture width (copy of input EEW_****)
;
;        .epdisc | REAL SCALAR | 1
;        Proportion of energy within an aperture of diameter EEW_DISC.
;
;        .epsqua | REAL SCALAR | 1
;        Proportion of energy within a square aperture of width EEW_SQUA.
;
;        .epslit | REAL SCALAR | 1
;        Proportion of energy within a slit width EEW_SLIT.
;
;###################################################################################
;###################################################################################
;############################                               ########################
;############################       SAVING THE RESULT       ########################
;############################                               ########################
;###################################################################################
;###################################################################################
;
;    LOGCODE= | STRING
;    If given, PAOLA writes in a log file 'paola'+LOGCODE+'.log' the main parameters 
;    value and results of the calculation, including the metrics above.
;
;    FITSCODE= | STRING
;    If a code is given, PAOLA saves the PSF, OTF, SF and all PSD matrices in FITS 
;    format, depending on which of the optional keywords /PSF,/OTF,/SF,/PSD have 
;    been set. The FITS file names are constructed according to 'PSF_'+FITSCODE
;    +'.fits', 'OTF_'+FITSCODE+'.fits', etc...
;
;FITS FILE HEADER CONTENT
;
;    MODE      Calculation mode                          [-]
;    INST      Includes imager pixel averaging           YES/NO
;    DXF       Focal plane pixel size                    [ASEC/PX]
;    DFF       Focal plane angular frequency pixel size  [1/RAD/PX]
;    DXP       Pupil plane pixel size                    [M/PX]
;    DFP       Pupil plane spatial frequency pixel size  [1/M/PX]
;    LAMBDA    Wavelength                                [MICRONS]
;    TSC_DIAM  Telescope primary mirror diameter         [METERS]
;    TSC_SURF  Telescope primary mirror surface          [METERS^2]
;    SEEING    Seeing angle @ 500 nm                     [ASEC]
;    r0_500nm  Fried's parameter @ 500 nm                [METERS]
;    L0        Outer scale of turbulence                 [METERS]
;    ALTITUDE  Mean altitude of the turbulent layers     [METERS]
;    ISOPLANA  Isoplanatic angle  @ <LAMBDA>             [ASEC]
;    VELOCITY  Mean velocity of the turbulent layers     [METERS/S]
;    TIME_SCL  Life time of the speckle field @ <LAMBDA> [MS]
;    >>>>>>>>>>ONLY FOR AO MODES-------------------------
;    NOMINAL   Nominal number of actuators over DM dia.  [1]
;    DM_NACT   Number of actuators over DM diameter      [1]
;    TOT_NACT  Total number of actuators under DM surf.  [1] 
;    ACTPITCH  Actuator pitch in the M1 plane            [METERS]
;    NLS_DIA   Number of WFS lenslets across pupil       [1]
;    NLS_PUP   Total number of WFS lenslets within pupil [1] 
;    WFS_PITC  WFS lenslet pitch in the M1 plane         [METERS]
;    DM_HEIGH  Conjugation altitude of the DM            [METERS]
;    ANG       NGS off-axis angle                        [ASEC]
;    ORI       NGS off-axis orientation                  [DEG]
;    INT       Wave-front sensor integration time        [MS]
;    LAG       Control system technical time-lag         [MS]
;    SERVOLAG  Total time lag of the system LAG+WFSINT   [MS]
;    LOOPFREQ  Loop frequency                            [Hz]
;    LOOPGAIN  Loop gain. Only in closed loop mode.      [1]
;    STABLE    Loop stability - 'yes' or 'no'            [-]
;    MAG       Bolometric magnitude of the NGS           [1]
;    TEM       Black-body temperature of the NGS         [K]
;    RON       Read-noise of the GS WFS CCD              [e/px]
;    TAU       WFS optics overall transmission           [1]
;    NPH       Number of photons/frame/lenslet           [1]
;    NEA       NGS WFS noise equivalent angle RMS        [ASEC]
;    DISP      Dispersion factor                         [1]
;    >>>>>>>>>>ONLY IN GLAO MODE-------------------------
;    GLAOTY    GLAO WFS measurement type                 [-]
;    WFSFOV    Average GLAO WFS Field-of-view            [ASEC]
;
;###################################################################################
;###################################################################################
;############################                               ########################
;############################          OTHER THINGS         ########################
;############################                               ########################
;###################################################################################
;###################################################################################
;
;ASSOCIATED FUNCTIONS, PROCEDURES AND FILES NEEDED TO RUN THE TOOLBOX PAOLA
;
;    aob_readpar.pro
;    attos_fwhm.pro
;    coogrid.pro
;    coomax.pro
;    discft.pro
;    dm_model.pro
;    err_exit.pro
;    fwhm_psf.pro
;    genshift.pro
;    hardcopy.pro
;    hexaft.pro
;    indzer.pro
;    valid_input.pro
;    intene.pro
;    intx.pro
;    mathft.pro
;    nbphotons.pro
;    nbphotons.sav
;    paola.pro
;    pixmatsize.pro
;    psd_alias_ngs.pro
;    psd_alias_ngs_glao.pro
;    psd_aniso_glao_edge.pro
;    psd_aniso_glao_full.pro
;    psd_aniso_servo.pro
;    psd_aniso_servo_glao_star.pro
;    psd_hferr_ngs.pro
;    psd_lferr_ngs.pro
;    psd_noise_ngs.pro
;    psd_noise_ngs_glao.pro
;    psfotftsc.pro
;    pupilarchitecture.pro
;    rectft.pro
;    remseg.pro
;    seemir.pro
;    sh_nea.pro
;    sinc.pro
;    spider.pro
;    triaft.pro
;    wave.pro
;
;MODIFICATIONS HISTORY
;
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    VERSION 1.0 March 20, 2001 Laurent Jolissaint HIA/NRC
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    
;    Apr 19, 2001 LJ real AO correction : aliasing, servo-lag,
;                 anisoplanatism and WFS noise implemented
;    Jan 29, 2003 LJ bug fixed : the deformable mirror spatial
;                 cutting frequency (1/(2d)) was limited to be < or =
;                 to the spatial sampling frequency of the phase in
;                 the matrices used along the computation. This is not
;                 the case anymore. 
;    Jan 30, 2003 LJ better error management of a fitting process
;                 when calculating the alias component of the phase
;                 power spectrum. Improved comments.
;    Jan 31, 2003 LJ bug corrected in polychromatic PSF calculation.
;    Feb  7, 2003 LJ replaced function unitgridfft() by COOGRID()
;    Feb  7, 2003 LJ visible wavelength @ 500 nm instead of 550 nm
;    Feb  8, 2003 LJ all in double precision
;
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    VERSION 2.0 March 13, 2003 Laurent Jolissaint HIA/NRC
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;    Mar 18, 2003 LJ included 'GLAO' mode & new FWHM.PRO function
;    Mar 28, 2003 LJ moved the LAMBDA input in a structure
;                 variable component of the input TSC.
;    Mar 28, 2003 LJ added the component .lam at the output
;                 structure variable.
;    Apr 04, 2003 LJ Added off-axis effect for 'GLAO' mode. The
;                 off-axis angle is given by the parameters OFFAXIS
;                 and OFFARG.
;    Apr 07, 2003 LJ Added the DHEIGHT input, thickness of the
;                 turbulent layers.
;    Apr 07, 2003 LJ Replaced the whites ' ' in the FITS file header
;                 keywords by underscores '_'. 
;    Apr 09, 2003 LJ layers are defined by their r0 rather than
;                 the relative Cn2 distribution.
;    Apr 09, 2003 LJ removed the NYQUIST warning, because FWHM.PRO
;                 gives good results even in the NYQUIST case.
;    Apr 10, 2003 LJ updated header for 'GLAO' mode
;    Apr 10, 2003 LJ implemented the Single Conjugate AO mode,
;             i.e. the DM can be conjugated to any
;                 altitude. See input HDM.
;    Apr 10, 2003 LJ changed the modes names to 'SELI' & 'SCAO'
;                 instead 
;                 of 'OFF' & 'ON' - because now, with GLAO, there
;                 are 2 AO modes.
;    Apr 10, 2003 LJ set the WIND input as an optional option for
;                 the seeing limited case. So, the layers mean
;                 velocity and the phase life time can be given on
;                 output.
;    Apr 10, 2003 LJ corrected a bug on the minimum number of
;                 actuators, when the default case is chosen
;                 (-1).
;    Apr 11, 2003 LJ changed the /LOG keyword for LOGCODE=, so a
;                 code name can be set for the log file.
;    Apr 11, 2003 LJ added on output the actual number of actuators
;                 over the 
;                 DM diameter and the total number of actuators
;                 under the DM surface (.ana & .tna).
;
;    VERSION 2.1 April 11, 2003 Laurent Jolissaint HIA/NRC
;
;    May 06, 2003 LJ added the residual phase total variance info @
;                 output 
;    Jul 31, 2003 LJ included the outer scale as an input parameter
;
;    VERSION 2.2 July  31, 2003 Laurent Jolissaint HIA/NRC
;
;    Oct 31, 2003 LJ corrected a bug in the logical unit management
;
;    VERSION 2.3 Nov 4, 2003 Laurent Jolissaint HIA/NRC
;
;    VERSION 2.3.1 Nov 7, 2003 Laurent Jolissaint HIA/NRC
;
;    Nov 12, 2003 LJ improved 50% and 80% encircled energy radius
;                 calculation, normalizing by lambda^2 / mirror
;                 surface rather than max of encircled
;                 energy.
;    Nov 13, 2003 LJ Introduced the TSC.INST input variable
;                 for focal plane instrumentation type.
;    Nov 13, 2003 LJ Removed the 50% and 80% energy radius
;                 calculation in the case of an imager (CCD or IR cam)
;                 instrumentation focal plane.
;    Nov 13, 2003 LJ added the .inst and .mode variables at the
;                 output structure variable.
;    Nov 21, 2003 LJ corrected a minor bug in GLAO mode
;
;    VERSION 2.3.2 Nov 21, 2003 Laurent Jolissaint HIA/NRC release of
;
;    Nov 26, 2003 LJ removed keyword STATUS in the POLY_FIT function
;    because it is not recognized by IDL version 5.6 and following
;    versions. 
;    Nov 28, 2003 LJ optimization of memory use.
;
;    VERSION 2.3.3 Nov 28, 2003 Laurent Jolissaint HIA/NRC
;
;    Dec 02, 2003 LJ & Rodolphe Conan (LAOG/ONERA) : implemented the
;                 exact bi-dimensional expression for the WFS aliasing
;                 power spectrum and the WFS noise power
;                 spectrum. Changed the spatial frequency domain
;                 corrected by the DM from the circle f<1/2d to the
;                 square |fx|<1/2d & |fy|<1/2d, where d is the
;                 actuator pitch as seen from the entrance pupil.
;                
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    VERSION 3.0 December 2, 2003 Laurent Jolissaint HIA/NRC
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;    Dec 15, 2003 LJ added keyword /DOUBLE to save the FITS files in
;                 DOUBLE rather than FLOAT precision. 
;    Apr 29, 2004 LJ removed constraint on power of 2 for matrix size.
;    May  3, 2004 JL (Jeff LeDue) corrected a bug in noise
;                 calculation. 
;    May  5, 2004 LJ changed function FWHM.PRO name to FWHM_PSF.PRO to
;                 avoid conflict with SIGMA_FILTER.PRO function.
;    May  7, 2004 LJ introduced ENCENE.PRO function to compute 50% and
;                 80% encircled energy radius. 
;    May  7, 2004 LJ removed constraint on integer number for actuator
;                 number across the pupil.
;    May 23, 2004 LJ TSC.INST can also be 'ccd' or 'spectro', not only
;                 in upper case letters.
;    Oct 13, 2004 LJ implemented option /ANTI_ALIAS, to simulate WFS
;                 aliasing filtering.
;    Feb 03, 2005 LJ corrected a bug in aliasing power spectrum
;    Mar 04, 2005 LJ improved time-lag and WFS integration time inputs
;                 management.
;    Mar 07, 2005 LJ defined a default value for entry DHEIGHT (-1)
;    Mar 13, 2005 LJ added keyword ANALYSIS for PSF analysis.
;    Mar 16, 2005 LJ improvement of the IDL on-line user manual, and
;                 re-organization of the code architecture.
;    Mar 18, 2005 LJ added ensquared energy and enslited energy
;                 metrics. 
;    Apr 06, 2005 LJ implemented LGS and TT-NGS.
;    May 02, 2005 LJ implemented GS_NEA & TT_NEA keywords.
;    May 04, 2005 LJ allowed EE calculation when INST is CCD type.
;    May 26, 2005 LJ introduced amplitude effects, via \SCINTILLATION
;                 keyword.
;    Jun 10, 2005 LJ introduced cross-coupling between anisoplanatism
;                 and servo-lag in NGS mode.
;    Jun 10, 2005 LJ made individual power spectrum available.
;    Jun 13, 2005 LJ force PAOLA to use an old release of writefits
;                 (Feb 22, 1992), because the current one (June 13,
;                 2005) simply does not work. 
;    Jun 14, 2005 LJ default actuator pitch set to r0(<lambda>)
;                 instead of r0(500nm)
;    Jun 16, 2005 LJ added aniso-servo tip-tilt correlation power 
;                 spectrum.
;
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    VERSION 4.0 June 24, 2005 Laurent Jolissaint HIA/NRC
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;    Jun 24, 2005 LJ corrected a buglet in a writefits routine.
;    Jun 25, 2005 VERSION 4.0.1
;
;    Jul 15, 2005 LJ replaced writefits920222 with fits_write.
;    Jul 15, 2005 VERSION 4.0.2
;
;    Aug 10, 2005 LJ added pupil plane spatial frequency pixel size to
;                 the output structure variable and in the fits files
;                 headers (DFP).
;    Aug 11, 2005 LJ added filtering of piston component in the WFS 
;                 noise power spectrum, because piston cannot be seen, 
;                 then not build, from the SH WFS.
;    Aug 11, 2005 LJ added NEA value to the output structure variable
;                 when GS_NEA directly specified. NEA written onto
;                 screen and log file too, in both cases (NGS_MAG and
;                 GS_NEA).
;    Aug 11, 2005 LJ noise was calculated twice. Bug corrected.
;    Aug 11, 2005 VERSION 4.0.3 released.
;
;    Sep 12, 2005 LJ added warning when matrix size is too small to
;                 get a good representation of the high frequency
;                 phase power spectrum wings.
;    Sep 19, 2005 VERSION 4.0.4 released.
;
;    Sep 23, 2005 LJ improved Scintillation calculation by adding
;                 correct amplitude fluctuation filtering within
;                 telescope pupil.
;    Sep 23, 2005 VERSION 4.0.5 released.
;
;    Feb 01, 2006 LJ corrected a bug in the scintillation routine
;                 (thanks Remi !).
;    Feb 01, 2006 VERSION 4.0.6 released.
;    Aug 30, 2006 LJ updated for new version of MATHFT.PRO
;    Sep 14, 2006 LJ added test for undefined inputs.
;
;    Sep 14, 2006 VERSION 4.1 released.
;
;    Sep 19, 2006 LJ improved NGS aliasing spectrum calculation, and 
;                 added forced Cn2 distribution normalization.
;
;    Sep 19, 2006 VERSION 4.1.1 released.
;
;    Sep 27, 2006 LJ corrected stupid buglet in 4.1.1 release.
;
;    Sep 27, 2006 VERSION 4.1.2 released.
;
;    Oct 05, 2006 LJ turbulence profiles - height, cn2 etc - can be 
;                 single scalar, now. This allows the study of
;                 servo-lag effects simply in terms of a given V/D
;                 ratio, or the case with a single layer at altitude
;                 "height".
;    Oct 06, 2006 LJ implemented WFS integration time optimization for
;                 a given NGS magnitude. See keyword /OPTIMIZE_WFS_INT.
;
;    Oct 10, 2006 VERSION 4.1.3 released.
;
;    Feb 20, 2007 LJ MAJOR reorganization of the code.
;                 Inputs cannot be given as if they were optional
;                 inputs anymore.
;                 Replaced NUMAOACT with ACTPITCH.
;                 Phase power spectrums removed for the main code, and
;                 set as independent external functions.
;                 Added 'edge' and 'star' GLAO modes.
;
;    February 2007 release of an unofficial beta version of the 5.0
;    to some people.
;
;    Apr 04, 2007 LJ squeezed the OTF and the structure functions
;                 arrays to the minimal matrix size, i.e. null rows
;                 and columns above the telescope cutoff frequency
;                 D/lambda are removed. This significantly speeds up
;                 the calculation time.
;    Apr 10, 2007 LJ added check of compatibility between psfotftsc's
;                 value of W0 and PAOLA input W0.
;    Apr 11, 2007 LJ computation of final PSF set as an option: active 
;                 only if the /PSF keyword is set, or /FWHM_ANALYSIS
;                 or /EE_ANALYSIS are set. This reduces execution time.
;
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    VERSION 5.0 May 2007 Laurent Jolissaint, Leiden Observatory, NL.
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;    May 25, 2007 LJ replaced /DOUBLE with /SINGLE. Double is the default for output.
;    May 25, 2007 LJ introduced aniso-servo power spectrum in GLAO star mode.
;    May 25, 2007 LJ introduced GS weighting in GLAO star mode.
;    May 29, 2007 LJ rewritten paola.pro header.
;    Jun 06, 2007 LJ removed loop over wavelength to simplify the code.
;    Jun 06, 2007 LJ introduced DM_PARAMS (replaces ACTPITCH).
;    Jun 07, 2007 LJ added test for undefined inputs.
;    Jun 07, 2007 LJ introduced DM filtering for all AO modes.
;    Jun 22, 2007 LJ added WFS_PITCH input.
;    Oct 01, 2007 LJ HEIGHT is now relative to pupil level, not sea level.
;    Nov 28, 2007 LJ improved routine's header, implemented new 
;                 versions of PIXMATSIZE and PSFOTFTSC functions.
;    Nov 28, 2007 LJ removed piston filtering and TT LGS mode in all modes.
;    Nov 28, 2007 LJ removed the prefix "gs_" in the outputs.
;
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    Dec 10, 2007 released an unofficial version with LGS mode set off
;                 to some people.
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;    Jan 08, 2008 LJ fixed buglet with var name definition conflict:
;                 PSF elongation orientation in output .dir, instead
;                 of .ori (conflict with .ori for GS position orientation).
;    Jan 19, 2008 LJ changed test of seeing angle compatibility with
;                 tsc.w0 by w0 < or = to tsc.w0, because if the seeing
;                 angle at PAOLA input is smaller than the seeing
;                 angle used when running PIXMATSIZE, it's not an issue.
;    Jan 21, 2008 LJ adapted to new version of FWHM_PSF.PRO function.
;
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    Jan 31, 2008 released an unofficial version with LGS mode set off
;                 to some people. See paola_LGS_OFF_v6.0.tar
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;    Feb 11, 2008 LJ high frequency structure function was missing for
;                 very high Strehl ratio.
;    Feb 12, 2008 LJ added check for actpitch and wfs_pitch: must be >= 2*lambda/FoV
;                 to keep the core to halo transition radius at least
;                 smaller or equal to half the field radius.
;    Feb 12, 2008 LJ WFS_PITCH input check wasn't implemented in GLAO mode.
;    Feb 27, 2008 LJ improved FITS file saving. PSD components saved
;                 as FITS files only if different from 0.
;    Mar 14, 2008 LJ changed SAVECODE for FITSCODE (more meaningful).
;    Mar 14, 2008 LJ Aniso PSD was saved twice (FITS) in GLAO mode. Corrected.
;    Mar 14, 2008 LJ Buglet corrected when using FITSCODE with /PSD set.
;    Mar 28, 2008 LJ split keyword EE_ANALYSIS into the three cases
;                 EE_ANALYSIS_DISC, ..._SQUA, ..._SLIT, and did the same
;                 for EEW -> EEW_DISC, EEW_SQUA and EEW_SLIT.
;    Apr 21, 2008 LJ PUT LGS MODE ON HOLD until bug fixed.
;
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    VERSION 6.0 April 2008 Laurent Jolissaint, Leiden Observatory, NL.
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;    Apr 29, 2008 LJ fixed bug when low frequency structure function
;                 is null.
;    May 02, 2008 LJ check that pupil plane spatial frequency
;                 resolution is < DM or WFS cutoff frequency. 
;    May 02, 2008 LJ added PARDISP optional input, to take into
;                 account effect of dispersion from WFS wavelength to
;                 correction (science) wavelength.
;
;  May 02, 2008 RELEASE OF VERSION 6.1
;
;    Jun 14, 2008 LJ added Richard Mathar's model of
;                 refractivity - see input PARDISP.
;    Jul 12, 2008 LJ improved GLAO "star" mode: introduced WFS noise
;                 power spectrum.
;    Jul 12, 2008 LJ added GLAO WFS integration time optimization.
;    Jul 12, 2008 LJ removed high spatial frequency warning message:
;                 wasn't useful anymore.
;    Jul 13, 2008 LJ added a check of the WFS integration time, loop
;                 time lag and GS to science object angular separation
;                 in order to detect bad sampling of the aniso-servo
;                 PSD terms.
;    Jul 13, 2008 LJ added WFS optics overall transmission factor WFS_TAU
;
;  Jul 13, 2008 RELEASE OF VERSION 6.2
;
;    Aug 30, 2008 LJ added GLAO WFS spatial aliasing. Now GLAO in star
;                 mode is complete.
;    Sep 10, 2008 LJ corrected a bug in x-coo matrix definition for seeing limited
;                 case, which was given Strehl of 1 for L0 not infinite.
;    Sep 30, 2008 LJ decreased limit for aniso-servo undersampling
;                 error. It has to be 8 times smaller than initially set.
;    Oct 11, 2008 LJ corrected a bug in GLAO star mode when setting WFS_NEA to 0.
;    Oct 11, 2008 LJ GS_WEIGHT couldn't be used when WFS noise and
;                 servo-lag error were neglected in GLAO star mode. Corrected now.  
;    Nov 15, 2008 LJ added keyword ONLY_PSD, improved GLAO user manual.
;    Nov 16, 2008 LJ changed routine AMOEBA for MINF_PARABOLIC in determination 
;                 of optimum integration time (a bit faster).
;    Dec 16, 2008 LJ improved servo-lag undersampling warning estimation.
;    Dec 17, 2008 LJ corrected a bug: WFS wavelength is now the one associated 
;                    with the NGS black body spectrum, not the science wavelength.
;    Dec 18, 2008 LJ added DXF, DFF, DXP, DFP in fits file header.
;    Dec 23, 2008 LJ improved off-axis maximal angle check in both NGS & GLAO.
;    Dec 24, 2008 LJ Now giving a wfs integration time <> -1 when
;                    optimize_wfs_int is set is possible. The value given
;                    will be used as the starting value for optimization.
;    Jan 07, 2009 LJ corrected a bug with GLAO full & edge modes: GS_WEIGHT not
;                    needed for these cases.
;    Jan 08, 2009 LJ removed SPECTRO mode request for integrated energy
;                    calculation.
;    Jan 08, 2009 LJ rewritten function's header.
;    Jan 08, 2009 LJ moved input GS_WEIGHT after GLAO_WFS input.
;    Feb 26, 2009 LJ added inputs DIM and ZA.
;    May 01, 2009 LJ ignoring layers for which Cn2 is null.
;    May 17, 2009 LJ removed OPTIMIZEWFSINTGLAO_BLOCK common block.
;    May 17, 2009 LJ taking into account NGS angle for WFS integration time optimization in
;                    NGS mode. I had forgotten to do it, but it was done in GLAO.
;    May 29, 2009 LJ introduced CLOSED LOOP in NGS mode, with loop gain
;                    optimization.
;    Jun 19, 2009 LJ introduced dispersion effect in NGS closed loop
;                    mode. Removed 'Ciddor' & 'Mathar' dispersion modes.
;    Jul 22, 2009 LJ buglet corrected (check of positivity of GS weight).
;    Jul 29, 2009 LJ improved computation of number of photo-electrons detected
;                    by lenslets/integration time.
;    Aug 09, 2009 LJ introduced VALID_INPUT.PRO to check inputs more easily.
;    Aug 09, 2009 LJ introduced WFS_PARAMS input, with mirror coating type, and
;                    more sophisticated WFS light sensitivity modeling.
;    Aug 17, 2009 LJ introduced NGS magnitude per filter instead of bolometric.
;    Aug 22, 2009 LJ checked 4-quadrant modeling of SH-WFS and added a
;                    complete model of centre-of-gravity based SH.
;    Aug 25, 2009 LJ moved DM_HEIGHT into DM_PARAMS structure variable input.
;    Sep 21, 2009 LJ LOOP_GAIN used as the starting value for the optimization.
;    Oct 20, 2009 LJ introduced turbulent layers height stretch with zenith
;                    angle.
;    Nov 24, 2009 LJ introduced cosmetic changes on the way the
;                     results are printed.
;    Mar 23, 2010 LJ aquilAOptics, introduced fractional Fourier transform.
;    Mar 26, 2010 LJ improved a lot the PSF FWHM computation thanks to
;                    fractional_fft algorithm usage.
;    Mar 30, 2010 LJ improved FWHM code once again for very large FWHM PSF.
;    Apr 03, 2010 LJ implemented G-tilt correction.
;    Jun 09, 2010 LJ improved computation of HF structure function for
;                    large WFS pitch values.
;    Aug 28, 2010 LJ added /FRFFT keyword for people who want to use
;                    the FRACTIONAL FFT stuff. By default fractional
;                    FFT is not used.
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    VERSION 7.0 Sept 2010 Laurent Jolissaint, aquilAOptics, Switzerland
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;    Sep 04, 2010 LJ corrected a bug with FRFFT call (for FWHM computation).
;    Sep 05, 2010 LJ released VERSION 7.0.1
;    Sep 20, 2010 LJ added test of DM transfer function array size
;    Oct 13, 2010 LJ added dispersion to the noise PSD model in closed loop.
;    Nov 19, 2010 LJ improved user manual.
;    Mar 02, 2011 LJ improved G-tilt correction model.
;    Mar 04, 2011 LJ there was some incompatible keywords in the FITS header.
;    May 26, 2011 LJ added DFP_LF & DFP_HF outputs.
;    Sep 05, 2011 LJ released VERSION 7.0.2
;    Sep 19, 2011 LJ there was a bug when wind=0. Corrected.
;    Dec 01, 2011 LJ improved significantly the WFS integration time
;                    and loop gain optimization part, by taking into
;                    account the validity domain for the WFS
;                    integration time and loop gain.
;    Dec 01, 2011 LJ included loop stability verification according to
;                    the Nyquist ]-1,0[ criterion, in NGS mode.
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    VERSION 7.1 Dec 2011 Laurent Jolissaint, aquilAOptics, Switzerland
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;    Feb 11, 2012 LJ buglet corrected in seeing limited mode when
;                    using the FITSCODE tag to save fits files (F. Patru)
;    Apr 05, 2012 LJ added w0ZA on output.
;    Jul 13, 2012 LJ re-introduced piston filtering in noise PSD.
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    VERSION 7.2 Jul 2012 Laurent Jolissaint, aquilAOptics, Switzerland
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;    Feb 15, 2014 LJ comment for loop stability removed in open loop
;    May 21, 2014 LJ added keyword /WAVE to prepare a residual phase
;                 PSD for use with the function wave.pro to create
;                 instantaneous AO-corrected phase screens.
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    VERSION 7.3 May 2014 Laurent Jolissaint, HEIG-VD, Switzerland
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;    Nov 10, 2014 LJ added /RADDMTF for radial DM transfer function.
;    May 12, 2015 LJ changed PSF option for OTF option in the
;                 computation of the /E50 & /E80.
;    Apr 12, 2015 LJ changed definition of WFS FoV in GLAO star mode.
;    Apr 19, 2017 LJ improved SH WFS noise model, removed JITTER, now
;                 automatically computed using the Fried parameter r0.
;                 And added a Gaussian spot model (Rousset et al. 1999).
;    Oct 10, 2017 LJ corrected a bug at line 3548 in /POST_TIPTILT option.
;    Jun 04, 2018 LJ total time delay due to servo lag was not
;                 computed correctly, and this made the performance
;                 too good at high magnitude. The time delay was too
;                 short by a factor 2 ! Now it is corrected.
;    Jun 30, 2018 LJ made sure that /PSF is set when /FWHM_AMALYSIS is
;                 required.
;    Oct 15, 2018 LJ added background noise option for the SH WFS.
;    Feb 27, 2020 LJ re-introduced piston filtering.
;    Feb 27, 2020 LJ ()->[] for arrays 
;    Feb 28, 2020 LJ removed PRECISION and SINGLE optional input, makes no sense
;                 anymore.
;    Feb 29, 2020 LJ I forgot to add the GS option for loop optimisation.
;    Mar 01, 2020 LJ set off sky background noise option until model is debuged.
;    Mar 01, 2020 LJ added wavefront RMS into the output structure variable.
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    VERSION 7.4 Mar 2020 Laurent Jolissaint, HEIG-VD, Switzerland
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;    May 01, 2020 LJ corrected a bug in /WAVE mode when using the tel pupil mask     
;    May 01, 2020 LJ added turbulent phase coherence time tau_0 on output
;    May 11, 2020 LJ rewritten the user manual, implemented vibrations option,
;                 added an input validity check to the aniso-servo
;                 error PSD, allowing now an overall optimization of
;                 loop gain and loop frequency at the same time.
;    May 11, 2020 LJ included many new options in the calculation of the telescope
;                 PSF (static post-AO, pre-AO, vibrations).
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    VERSION 7.5 May 2020 Laurent Jolissaint, HEIG-VD, Switzerland
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;    Dec 08, 2020 LJ removed a condition in the function
;                 PSD_ANISO_SERVO that was too restrictive and
;                 generated crashes in optimisation mode.
;    Jan 22, 2021 LJ limited the WFS integration time optimisation. It
;                 was not realistic to do an optimisation beyond a few 10's of
;                 ms, because Taylor frozen flow assumption does not hold too
;                 long. After some time, boiling dominates, and PAOLA servo-lag model
;                 does not hold anymore. So I have limited the WFS integration time
;                 optimisation to the phase time scale at which the observation takes
;                 place, and this might be already a bit too optimistic. Issue was the
;                 same for GLAO. Also the loop gain optimization has been limited to
;                 values <= 1 as overcorrecting has no real sense here.
;    Jan 27, 2021 LJ WFS integration time cannot be 0, because in this case
;                 the WFS noise is not defined. An input check has been added.
;    Jan 27, 2021 LJ changed PSF mode SPECTRO to INTENSITY because it was
;                 misleading.
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    VERSION 7.6 January 2021 Laurent Jolissaint, HEIG-VD, Switzerland
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;    Jan 29, 2021 LJ an integrator gain greater than 1 is not forbiden
;                 from control theory, so this limit has been lifted.
;    Jan 29, 2021 LJ removed warning relative to outer scale. It had
;                 no interest.
;    Aug 30, 2021 LJ massively improved the calculation of the number
;                 of photons received by the wavefront sensor (better
;                 atmosphere transmission model, better mirrors
;                 reflectivity models, polarization effects).
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    VERSION 7.6.2 September 2021 Laurent Jolissaint, HEIG-VD, Switzerland
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;    Apr 01, 2022 LJ WFS NEA keyword is not compatible with loop
;                 optimization, nor computation of the number of
;                 photons, so if one ask this and at the same time set
;                 the WFS noise equivalent angle, an error message is
;                 generated. 
;    Oct 27, 2022 LJ added the loop frequency input, to make it
;                 independant from the WFS integration time.
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;    VERSION 7.6.3 November 2022 Laurent Jolissaint, HEIG-VD, Switzerland
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;    Dec 21, 2022 LJ Set option for loop frequency to -1
;                 to get the default value which is 1000/wfsint.
;    Dec 22, 2022 LJ added an optional limit for the loop frequency
;                 in optimization mode, in order to avoid the
;                 code to produce unrealistic values in loop
;                 frequency and loop gain optimization.
;
;BUGS & QUESTIONS
;
;    Write to laurent.jolissaint@heig-vd.ch
;
;ACKNOWLEDGMENTS
;
;    Special thanks to Jean-Pierre Veran, François Rigaut, Jeff LeDue, Rodolphe
;    Conan, Olivier Lardiere, Jean-Francois Lavigne, Christophe Verinaud, Miska
;    Le Louarn, Remi Soummer, Bruno Femenia, Visa Korskiakoski,
;    Jeff Stoesz and many others, for helping me developing the
;    toolbox by providing ideas, bugs reports and great software packages
;    to implement into the toolbox. And teaching me adaptive optics !
;
;-
FUNCTION OPTIMIZE_NGS,P

  COMMON OPT_BLOCK,B

  if B.optmode eq 1 then begin ; optimization of loop frequency
    dTmax = B.timescale
    if strlowcase(B.lm) eq 'closed' then begin ; CLOSED loop
      dTmin = B.lag*tan(0.7692056d*B.lg-0.054163105d*B.lg^3+0.0030530390d*B.lg^5-9.3010157e-05*B.lg^7+1.1401983d-06*B.lg^9) ; required for loop stability
      dTmin = dTmin > 0.09 ; required to avoid unrealistically short exposure time
      if where(tag_names(B) eq 'MAX_LOOP_FREQ') ne -1 then dTmin = max([dTmin,1d3/B.max_loop_freq])
      if where(tag_names(B) eq 'MAX_LOOP_FREQ') ne -1 then P[0] = max([P[0],alog10(1d3/B.max_loop_freq)])
      dt = (10.d)^P[0]
      t0 = min([dTmin,dTmax]) ; this is because dTmin might be > dTmax depending on loop gain and lag
      t1 = max([dTmin,dTmax]) ; same
      if dt lt t0 then return,1d99
      if dt gt t1 then return,1d99
    endif else begin ; OPEN loop
      if where(tag_names(B) eq 'MAX_LOOP_FREQ') ne -1 then P[0] = alog10(1d3/B.max_loop_freq)
      dt = (10.d)^P[0]
    endelse
    if dt ge dTmax then return,1d99 ; to avoid going too far, which is unrealistic because boiling occurs
    var = total(PSD_ANISO_SERVO(B.lm,1000/dt,B.lg,B.n_psd,B.fpcoo,B.lambda,B.height,B.wind,B.r0500_i,B.l0,0,B.wfs_pitch,B.dmtf,B.dmh,B.ang,B.ori,dt,B.lag,B.ddisp,B.dextmax))
    tmp = NBPHOTONS(B.wfs_pitch^2,B.mirsv,B.nblenses,B.za,B.ngs_mag,B.filter,B.ngs_tem,dt,-1,EXTRAFILTER = B.extrafilter)
    if B.lm eq   'open' and strlowcase(B.algor) eq '4q' then $
       nea = (SH_NEA(tmp.nph,B.r0500,tmp.lam,B.wfs_pitch,B.wfs_ron,!dpi*tmp.lam^3*B.sky_radiance*dt/3.1783133d*1e3,B.algor)).nea
    if B.lm eq 'closed' and strlowcase(B.algor) eq '4q' then $
       nea = (SH_NEA(tmp.nph,B.r0500,tmp.lam,B.wfs_pitch,B.wfs_ron,!dpi*tmp.lam^3*B.sky_radiance*dt/3.1783133d*1e3,B.algor,dt,B.lag,B.lg)).nea
    if B.lm eq   'open' and (strlowcase(B.algor) eq 'cg' or strlowcase(B.algor) eq 'gs') then $
       nea = (SH_NEA(tmp.nph,B.r0500,tmp.lam,B.wfs_pitch,B.wfs_ron,!dpi*tmp.lam^3*B.sky_radiance*dt/3.1783133d*1e3,B.algor,B.wfs_fov,B.wfs_pxs)).nea
    if B.lm eq 'closed' and (strlowcase(B.algor) eq 'cg' or strlowcase(B.algor) eq 'gs') then $
       nea = (SH_NEA(tmp.nph,B.r0500,tmp.lam,B.wfs_pitch,B.wfs_ron,!dpi*tmp.lam^3*B.sky_radiance*dt/3.1783133d*1e3,B.algor,B.wfs_fov,B.wfs_pxs,dt,B.lag,B.lg)).nea
    var = var+total(PSD_NOISE_NGS(B.n_psd,B.fpcoo,B.wfs_pitch,B.dmtf,nea,B.lambda,B.ddisp,B.dextmax))
    if B.alias then var = var+total(PSD_ALIAS_NGS(B.n_psd,B.fpcoo,B.lambda,B.r0500_i,B.L0,0,B.wind,B.wfs_pitch,dt,B.lag,B.lg,B.lm,B.dmtf,B.ddisp,B.dextmax))
    if B.prt then printf,-1,format = '(2e20.12)',[1d3/dt,sqrt(var*B.dfpLAM^2)*B.lambda*1e3/2/!dpi]
  endif

  if B.optmode eq 2 then begin ; optimization of loop gain
    atdt = atan(B.wfs_int/B.lag)
    g_max = 1.3890312d*atdt-0.017740598d*atdt^3-0.20128461d*atdt^5+0.58676804d*atdt^7-0.34550100d*atdt^9+0.076581176d*atdt^11
    gain = (atan(P[0])+!dpi/2)/!dpi*g_max
    if gain le 0 or gain gt 4.9348d then return,1d99  ; this is a brick wall to avoid forbiden values of loop gain
    var = total(PSD_ANISO_SERVO(B.lm,B.lf,gain,B.n_psd,B.fpcoo,B.lambda,B.height,B.wind,B.r0500_i,B.l0,0,B.wfs_pitch,B.dmtf,B.dmh,B.ang,B.ori,B.wfs_int,B.lag,B.ddisp,B.dextmax))
    tmp = NBPHOTONS(B.wfs_pitch^2,B.mirsv,B.nblenses,B.za,B.ngs_mag,B.filter,B.ngs_tem,B.wfs_int,-1,EXTRAFILTER = B.extrafilter)
    ;noise equivalent angle
    if strlowcase(B.algor) eq '4q' then nea = (SH_NEA(tmp.nph,B.r0500,tmp.lam,B.wfs_pitch,B.wfs_ron,!dpi*tmp.lam^3*B.sky_radiance*B.wfs_int/3.1783133d*1e3,B.algor,B.wfs_int,B.lag,gain)).nea
    if strlowcase(B.algor) eq 'cg' or strlowcase(B.algor) eq 'gs' then $
       nea = (SH_NEA(tmp.nph,B.r0500,tmp.lam,B.wfs_pitch,B.wfs_ron,!dpi*tmp.lam^3*B.sky_radiance*B.wfs_int/3.1783133d*1e3,B.algor,B.wfs_fov,B.wfs_pxs,B.wfs_int,B.lag,gain)).nea
    var = var+total(PSD_NOISE_NGS(B.n_psd,B.fpcoo,B.wfs_pitch,B.dmtf,nea,B.lambda,B.ddisp,B.dextmax))
    if B.alias then var = var+total(PSD_ALIAS_NGS(B.n_psd,B.fpcoo,B.lambda,B.r0500_i,B.L0,0,B.wind,B.wfs_pitch,B.wfs_int,B.lag,gain,B.lm,B.dmtf,B.ddisp,B.dextmax))
    if B.prt then printf,-1,format = '(2e20.12)',[gain,sqrt(var*B.dfpLAM^2)*B.lambda*1e3/2/!dpi]
  endif

  if B.optmode eq 3 then begin ; optimization of WFS integration time & loop gain
    dTmax = B.timescale
    dTmin = B.lag*tan(0.7692056d*B.lg-0.054163105d*B.lg^3+0.0030530390d*B.lg^5-9.3010157e-05*B.lg^7+1.1401983d-06*B.lg^9) ; required for loop stability
    dTmin = dTmin > 0.09 ; required to avoid unrealistically short exposure time
    if where(tag_names(B) eq 'MAX_LOOP_FREQ') ne -1 then dTmin = max([dTmin,1d3/B.max_loop_freq])
    if where(tag_names(B) eq 'MAX_LOOP_FREQ') ne -1 then P[0] = max([P[0],1d3/B.max_loop_freq])
    dt = P[0]
    t0 = min([dTmin,dTmax]) ; this is because dTmin might be > dTmax depending on loop gain and lag
    t1 = max([dTmin,dTmax]) ; same
    if dt lt t0 then return,1d99
    if dt gt t1 then return,1d99
    gain = P[1]
    if gain le 0 or gain gt 4.9348d then return,1d99  ; this is a brick wall to avoid forbiden values of loop gain
    var = total(PSD_ANISO_SERVO(B.lm,1000/dt,gain,B.n_psd,B.fpcoo,B.lambda,B.height,B.wind,B.r0500_i,B.l0,0,B.wfs_pitch,B.dmtf,B.dmh,B.ang,B.ori,dt,B.lag,B.ddisp,B.dextmax))
    tmp = NBPHOTONS(B.wfs_pitch^2,B.mirsv,B.nblenses,B.za,B.ngs_mag,B.filter,B.ngs_tem,dt,-1,EXTRAFILTER = B.extrafilter)
    ;noise equivalent angle
    if strlowcase(B.algor) eq '4q' then nea = (SH_NEA(tmp.nph,B.r0500,tmp.lam,B.wfs_pitch,B.wfs_ron,!dpi*tmp.lam^3*B.sky_radiance*dt/3.1783133d*1e3,B.algor,dt,B.lag,gain)).nea
    if strlowcase(B.algor) eq 'cg' or strlowcase(B.algor) eq 'gs' then $
       nea = (SH_NEA(tmp.nph,B.r0500,tmp.lam,B.wfs_pitch,B.wfs_ron,!dpi*tmp.lam^3*B.sky_radiance*dt/3.1783133d*1e3,B.algor,B.wfs_fov,B.wfs_pxs,dt,B.lag,gain)).nea
    var = var+total(PSD_NOISE_NGS(B.n_psd,B.fpcoo,B.wfs_pitch,B.dmtf,nea,B.lambda,B.ddisp,B.dextmax))
    if B.alias then var = var+total(PSD_ALIAS_NGS(B.n_psd,B.fpcoo,B.lambda,B.r0500_i,B.L0,0,B.wind,B.wfs_pitch,dt,B.lag,gain,B.lm,B.dmtf,B.ddisp,B.dextmax))
    if B.prt then printf,-1,format = '(4e20.12)',[1d3/dt,gain,sqrt(var*B.dfpLAM^2)*B.lambda*1e3/2/!dpi]
  endif

  return,var*B.dfpLAM^2

end
; GLAO WFS integration time optimization
FUNCTION OPTIMIZE_GLAO,P

  COMMON OPT_BLOCK,B

  ; limiting dT to phase time scale at science lambda time
  dTmax = B.timescale
  dt = P[0]
  if dTmax le 0.01 then message,'IMPOSSIBLE OPTIMISATION OF WFS EXPOSURE TIME: THE MINIMUM ACCEPTABLE EXPOSURE TIME IS LARGER THAN THE MAXIMUM ACCEPTABLE EXPOSURE TIME'
  dt = dt > 0.01  ; required to avoid unrealistically short exposure time
  if dt ge dTmax then return,1d99 ; to avoid going too far, which is unrealistic because boiling occurs

  var = total(PSD_ANISO_SERVO_GLAO_STAR(B.n_psd,B.fpcoo,B.lambda,B.height,B.wind,B.r0500_i,B.l0,B.wfs_pitch,B.dmtf,B.dmh,B.ang,B.ori,dt,B.lag,B.gs_weight,B.glao_wfs,B.dextmax))
  tmp = NBPHOTONS(B.wfs_pitch^2,B.mirsv,B.nblenses,B.za,B.ngs_mag,B.filter,B.ngs_tem,abs(double(dt)),-1,EXTRAFILTER = B.extrafilter)
  if strlowcase(B.algor) eq '4q' then nea = (SH_NEA(tmp.nph,B.r0500,tmp.lam,B.wfs_pitch,B.wfs_ron,!dpi*tmp.lam^3*B.sky_radiance[0]*dt/3.1783133d*1e3,B.algor)).nea
  if strlowcase(B.algor) eq 'cg' or strlowcase(B.algor) eq 'gs' then nea = (SH_NEA(tmp.nph,B.r0500,tmp.lam,B.wfs_pitch,B.wfs_ron,!dpi*tmp.lam^3*B.sky_radiance[0]*dt/3.1783133d*1e3,B.algor,B.wfs_fov,B.wfs_pxs)).nea
  var = var+total(PSD_NOISE_NGS_GLAO(B.n_psd,B.fpcoo,B.wfs_pitch,B.dmtf,B.gs_weight,nea,B.lambda,B.dextmax))
  if B.alias then var = var+total(PSD_ALIAS_NGS_GLAO(B.n_psd,B.fpcoo,B.lambda,B.height,B.wind,B.r0500_i,B.L0,0,B.wfs_pitch,B.dmtf,B.dmh,dt,B.gs_weight,B.glao_wfs,B.dextmax))
  if B.prt then printf,-1,format = '(2e20.12)',[dt,sqrt(var*B.dfpLAM^2)*B.lambda*1e3/2/!dpi]

  return,var*B.dfpLAM^2

end
;;;;;;;;;;;;;;;
; MAIN FUNCTION
;;;;;;;;;;;;;;;
FUNCTION PAOLA,MODE,DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,P09,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,$
               ANTI_ALIAS = ANTI_ALIAS,$
               DISPERSION = DISPERSION,$
               EE_ANALYSIS_DISC = EE_ANALYSIS_DISC,$
               EE_ANALYSIS_SQUA = EE_ANALYSIS_SQUA,$
               EE_ANALYSIS_SLIT = EE_ANALYSIS_SLIT,$
               EEW_DISC = EEW_DISC,$
               EEW_SQUA = EEW_SQUA,$
               EEW_SLIT = EEW_SLIT,$
               FITSCODE = FITSCODE,$
               FRFFT = FRFFT,$
               FWHM_ANALYSIS = FWHM_ANALYSIS,$
               INFO = INFO,$
               LOGCODE = LOGCODE,$
               MAX_LOOP_FREQ = MAX_LOOP_FREQ,$
               ONLY_PSD = ONLY_PSD,$
               OPTIMIZE_WFS_INT = OPTIMIZE_WFS_INT,$
               OPTIMIZE_LOOP_FREQ = OPTIMIZE_LOOP_FREQ,$
               OPTIMIZE_LOOP_GAIN = OPTIMIZE_LOOP_GAIN,$
               OPTIMIZE_ALL = OPTIMIZE_ALL,$
               OTF = OTF,$
               VIBRATION = VIBRATION,$
               POST_TIPTILT = POST_TIPTILT,$
               PSD = PSD,$
               PSF = PSF,$
               RADDMTF = RADDMTF,$
               SCINTILLATION = SCINTILLATION,$
               SF = SF,$
               TILT_ANGLE_STEP = TILT_ANGLE_STEP,$
               X_PSF = X_PSF,$
               Y_PSF = Y_PSF,$
               WAVE = WAVE

; seeing limited inputs
;               6 params: MODE,DIM,TSC,W0,L0,ZA
;               9 params: MODE,DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,WIND
;
; NGS AO inputs
;               19 params: MODE,DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,P09,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19
;               21 params: MODE,DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,P09,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21
;
; GLAO inputs
;               13 params: MODE,DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,P09,P10,P11,P12,P13
;               20 params: MODE,DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,P09,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20
;               22 params: MODE,DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,P09,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22

  COMMON OPT_BLOCK,block

  if strlowcase(strtrim(MODE,2)) eq 'lgs' then message,'THIS IS A VERSION OF THE PAOLA TOOLBOX WHERE LGS IS NOT DEBUGGED, SO CANNOT BE USED YET.'

  ;==================================================================================== 
  ;================================= ARGUMENTS CHECK ==================================== 
  ;==================================================================================== 

  ;------------------------------
  ;CHECKING IF INPUTS ARE DEFINED
  ;------------------------------
  if where(n_params() eq [6,9,13,19,20,21,22]) eq -1 then message,'UNCORRECT NUMBER OF INPUTS: MUST BE 6 OR 9 IN SEEING LIMITED MODE, 19 OR 21 IN NGS MODE, 13, 20 OR 22 IN GLAO MODE'
  if float(size(MODE,/type))*float(size(DIM,/type))*float(size(TSC,/type))*float(size(W0,/type))*float(size(L0,/type))*float(size(ZA,/type)) eq 0 then $
    message,'ONE OR SEVERAL INPUTS ARE UNDEFINED: MODE, DIM, TSC, W0, L0 OR ZA'
  if n_params() ge  9 then if float(size(HEIGHT,/type))*float(size(DISTCN2,/type))*float(size(P09,/type)) eq 0 then $
    message,'ONE OR SEVERAL INPUTS ARE UNDEFINED: CHECK "HEIGHT", "DISTCN2", AND 9th INPUT'
  if n_params() ge 13 then if float(size(P10,/type))*float(size(P11,/type))*float(size(P12,/type))*float(size(P13,/type)) eq 0 then $
    message,'INPUTS # 10, 11, 12, OR 13 ARE UNDEFINED'
  if n_params() ge 19 then if float(size(P14,/type))*float(size(P15,/type))*float(size(P16,/type))*float(size(P17,/type))*float(size(P18,/type))*float(size(P19,/type)) eq 0 then $
    message,'INPUTS # 14, 15, 16, 17, 18 OR 19 ARE UNDEFINED'
  if n_params() ge 20 then if float(size(P20,/type)) eq 0 then message,'INPUT # 20 UNDEFINED'
  if n_params() ge 21 then if float(size(P21,/type)) eq 0 then message,'INPUT # 21 UNDEFINED'
  if n_params() eq 22 then if float(size(P22,/type)) eq 0 then message,'INPUT # 22 UNDEFINED'

  ;---------------------
  ;CHECKING IMAGING MODE
  ;---------------------
  VALID_INPUT,'PAOLA.PRO','MODE',MODE,'string',0,'no',['seli','ngs','glao','SELI','NGS','GLAO'],'free'

  ;-------------------------
  ;CHECKING NUMBER OF INPUTS
  ;-------------------------
  if strlowcase(strtrim(MODE,2)) eq 'seli' and n_params() ne 6 and n_params() ne 9 then $
    message,'NEED 6 OR 9 INPUTS IN SEEING LIMITED MODE'
  if strlowcase(strtrim(MODE,2)) eq 'ngs'  and n_params() ne 19 and n_params() ne 21 then $
    message,'NEED 19 OR 21 INPUTS IN NGS MODE'
  if strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() ne 13 and n_params() ne 20 and n_params() ne 22 then $
    message,'NEED 13, 20 OR 22 INPUTS IN GLAO MODE'

  ;--------------------------------------
  ;CHECKING VARIABLES COMMON TO ALL MODES
  ;--------------------------------------
  VALID_INPUT,'PAOLA.PRO','DIM',DIM,'structure',[1,-1],'no','free','free'
  VALID_INPUT,'PAOLA.PRO','TSC',TSC,'structure',[1,-1],'no','free','free'
  VALID_INPUT,'PAOLA.PRO','W0',float(W0[0]),'real',0,'no','++','free'
  VALID_INPUT,'PAOLA.PRO','L0',float(L0[0]),'real',0,'yes','++','free'
  VALID_INPUT,'PAOLA.PRO','ZA',float(ZA[0]),'real',0,'no',[0,90],'free'
  if n_params() gt 6 then begin
    VALID_INPUT,'PAOLA.PRO','HEIGHT',float(HEIGHT),'real','free','no','free','free'
    if n_elements(HEIGHT) gt 1 then VALID_INPUT,'PAOLA.PRO','HEIGHT',float(HEIGHT),'real',[1,-1],'no','free','free'
    VALID_INPUT,'PAOLA.PRO','DISTCN2',DISTCN2,'real','free','no','++','free'
    if n_elements(DISTCN2) gt 1 then VALID_INPUT,'PAOLA.PRO','DISTCN2',DISTCN2,'real',1,'no','++','free'
    if n_elements(HEIGHT) ne n_elements(DISTCN2) then message,'DISTCN2 AND ALTITUDE INPUTS MUST HAVE THE SAME NUMBER OF LAYERS'
  endif
  if n_params() ge 9 and not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then begin
    VALID_INPUT,'PAOLA.PRO','WIND',float(P09),'real','free','no','free','free'
    if (size(P09))[0] eq 0 then P09 = [[P09],[0]]
    VALID_INPUT,'PAOLA.PRO','WIND',float(P09),'real',[2,-1,2],'no','free','free'
    if (size(P09))[1] ne n_elements(HEIGHT) then message,'WIND AND ALTITUDE INPUTS MUST HAVE THE SAME NUMBER OF LAYERS'
  endif
  if (where(tag_names(DIM) eq 'W0'))[0] eq -1 then message,$
    'SEEING ANGLE HAS NOT BEEN TAKEN INTO ACCOUNT IN THE DETERMINATION OF MATRIX SIZES AND PIXEL SIZES. RE-RUN PIXMATSIZE FUNCTION WITH W0 KEYWORD SET.'
  if DIM.W0 lt W0 then message,$
    'THE SEEING ANGLE YOU HAVE SPECIFIED HERE IS LARGER THAN THE ONE USED TO DEFINE PIXEL AND MATRIX SIZES, '+$
    'THEREFORE THE PSF WOULD BE ALIASED: RE-RUN PIXMATSIZE WITH THE CORRECT SEEING ANGLE, OR DECREASE THE SEEING ANGLE AT PAOLA INPUT'
  if size(VIBRATION,/type) ne 0 then begin
    VALID_INPUT,'PAOLA.PRO','VIBRATION',VIBRATION,'structure',{dim:[1,3],tags:['SIGMA_TTX','SIGMA_TTY','ORIENTATION']},'no','free','free'
    VALID_INPUT,'PAOLA.PRO','VIBRATION.SIGMA_TTX',VIBRATION.SIGMA_TTX,'real',0,'no','0+','free'
    VALID_INPUT,'PAOLA.PRO','VIBRATION.SIGMA_TTY',VIBRATION.SIGMA_TTY,'real',0,'no','0+','free'
    VALID_INPUT,'PAOLA.PRO','VIBRATION.ORIENTATION',VIBRATION.ORIENTATION,'real',0,'no','free','free'
  endif

  ;------------------------------
  ;CHECKING SEEING LIMITED MODE
  ; is included in previous tests
  ;------------------------------

  ;-----------------
  ;CHECKING NGS MODE
  ;-----------------
  if strlowcase(strtrim(MODE,2)) eq 'ngs' then begin
    VALID_INPUT,'PAOLA.PRO','DM_PARAMS',P10,'structure',{dim:[1,3],tags:['dm_height','dmtf','actpitch']},'no','free','free'
    VALID_INPUT,'PAOLA.PRO','DM_PARAMS.DM_HEIGHT',P10.DM_HEIGHT,'real',0,'no','free','free'
    VALID_INPUT,'PAOLA.PRO','DM_PARAMS.DMTF',P10.DMTF,'real',[2,DIM.N_LF,DIM.N_LF],'yes','free','free'
    VALID_INPUT,'PAOLA.PRO','DM_PARAMS.ACTPITCH',P10.ACTPITCH,'real',0,'yes','++','free'
    VALID_INPUT,'PAOLA.PRO','WFS_PARAMS',P11,'structure',[1,-1],'no','free','free'
    if (where([1,8,10] eq n_tags(P11)))[0] eq -1 then message,'WFS_PARAMS STRUCTURE VARIABLE INPUT MUST HAVE 1, 8, OR 10 COMPONENTS'
    if n_tags(P11) eq 1 then VALID_INPUT,'PAOLA.PRO','WFS_PARAMS',P11,'structure',{dim:[1,1],tags:['wfs_pitch']},'no','free','free'
    if n_tags(P11) eq 8 then VALID_INPUT,'PAOLA.PRO','WFS_PARAMS',P11,'structure',{dim:[1,8],tags:['mircoating','mirangle','nblenses','extrafilter','wfs_pitch','wfs_ron','sky_radiance','algorithm']},'no','free','free'
    if n_tags(P11) eq 10 then VALID_INPUT,'PAOLA.PRO','WFS_PARAMS',P11,'structure',$
                                          {dim:[1,10],tags:['mircoating','mirangle','nblenses','extrafilter','wfs_pitch','wfs_ron','sky_radiance','algorithm','wfs_pxfov','wfs_pxsize']},'no','free','free'
    VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.WFS_PITCH',P11.WFS_PITCH,'real',0,'yes','++','free'
    VALID_INPUT,'PAOLA.PRO','NGS_ANG',P12,'real',0,'no','0+','free'
    VALID_INPUT,'PAOLA.PRO','NGS_ORI',P13,'real',0,'no','free','free'
    VALID_INPUT,'PAOLA.PRO','WFS_INT',P14,'real',0,'yes','++','free'
    VALID_INPUT,'PAOLA.PRO','LAG',P15,'real',0,'no','0+','free'
    VALID_INPUT,'PAOLA.PRO','LOOP_MODE',P16,'string',0,'no',['open','closed','OPEN','CLOSED'],'free'
    VALID_INPUT,'PAOLA.PRO','LOOP_FREQ',P17,'real',0,'yes','++',0
    VALID_INPUT,'PAOLA.PRO','LOOP_GAIN',P18,'real',0,'no',[0,4.9348d],0
    if P17 ne -1 then if P14-1000.d/P17 gt 0.001*P14 then message,'WFS INTEGRATION TIME IS > THAN THE LOOP SAMPLING TIME: CHANGE EITHER WFS_INT OR THE LOOP FREQUENCY'
    if n_params() eq 19 then VALID_INPUT,'PAOLA.PRO','WFS_NEA',P19,'real',0,'no','0+','free'
    if n_params() eq 21 then begin
      if n_tags(P11) eq 1 then message,'IF YOU WANT TO EVALUATE THE EFFECT OF THE NGS MAGNITIUDE, YOU NEED TO SET THE WFS PARAMETERS TOO: SEE INPUT WFS_PARAMS IN THE USER MANUAL (FILE HEADER) '
      VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.MIRCOATING',P11.MIRCOATING,'string',1,'yes',['al','ag','au','Al','Ag','Au','AL','AG','AU','aL','aG','aU'],'free'
      VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.MIRANGLE',P11.MIRANGLE,'real',1,'yes',[0,90],90
      VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.NBLENSES',P11.NBLENSES,'integer',0,'yes','0+','free'
      if size(P11.EXTRAFILTER,/type) ne 7 and size(P11.EXTRAFILTER,/type) ne 8 then message,'WFS_PARAMS.EXTRAFILTER SHALL BE EITHER A STRUCTURE VARIABLE OR ''NO'' '
      if size(P11.EXTRAFILTER,/type) eq 7 then VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.EXTRAFILTER',P11.EXTRAFILTER,'string',0,'no','no','free'
      if size(P11.EXTRAFILTER,/type) eq 8 then begin
        VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.EXTRAFILTER',P11.EXTRAFILTER,'structure',{dim:[1,2],tags:['lambda','tauext']},'no','free','free'
        VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.EXTRAFILTER.LAMBDA',P11.EXTRAFILTER.LAMBDA,'real',1,'yes',[0.3d,5.d],'free'
        VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.EXTRAFILTER.TAUEXT',P11.EXTRAFILTER.TAUEXT,'real',1,'yes',[0,1],'free'
      endif
      VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.WFS_RON',P11.WFS_RON,'real',0,'no','0+','free'
      VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.SKY_RADIANCE',P11.SKY_RADIANCE,'real',0,'no','0+','free'
      VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.ALGORITHM',P11.ALGORITHM,'string',0,'no',['4q','4Q','cg','Cg','cG','CG','gs','Gs','gS','GS'],'free'
      if strlowcase(P11.ALGORITHM) eq 'cg' then begin
        VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.WFS_PXFOV',P11.WFS_PXFOV,'integer',0,'no','++','free'
        if P11.WFS_PXFOV lt 6 then message,'NUMBER OF PIXEL ACROSS THE WFS LENSLET FOV MUST BE > OR  =  TO 6'
        if P11.WFS_PXFOV mod 2 eq 1 then message,'NUMBER OF PIXEL ACROSS THE WFS LENSLET FOV MUST BE EVEN'
        VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.WFS_PXSIZE',P11.WFS_PXSIZE,'real',0,'yes','++','free'
      endif
      VALID_INPUT,'PAOLA.PRO','NGS_MAG',P19,'real',0,'no','free','free'
      VALID_INPUT,'PAOLA.PRO','FILTER',P20,'string',0,'no',['u','b','v','r','i','j','h','k','l','m','n','q','U','B','V','R','I','J','H','K','L','M','N','Q'],'free'
      VALID_INPUT,'PAOLA.PRO','NGS_TEM',P21,'real',0,'no','++','free'
    endif
  endif

  ;------------------
  ;CHECKING GLAO MODE
  ;------------------
  if strlowcase(strtrim(MODE,2)) eq 'glao' then begin
    if (size(HEIGHT))[0] ne 1 then message,'LAYERS ALTITUDES (HEIGHT) MUST BE A VECTOR IN GLAO MODE'
    if (size(DISTCN2))[0] ne 1 then message,'LAYERS CN2DH DISTRIBUTION (DISTCN2) MUST BE A VECTOR IN GLAO MODE'
    if n_params() eq 13 then begin
      VALID_INPUT,'PAOLA.PRO','DM_PARAMS',P09,'structure',{dim:[1,3],tags:['dm_height','dmtf','actpitch']},'no','free','free'
      VALID_INPUT,'PAOLA.PRO','DM_PARAMS.DM_HEIGHT',P09.DM_HEIGHT,'real',0,'no','free','free'
      VALID_INPUT,'PAOLA.PRO','DM_PARAMS.ACTPITCH',P09.ACTPITCH,'real',0,'yes','++','free'
      VALID_INPUT,'PAOLA.PRO','DM_PARAMS.DMTF',P09.DMTF,'real',2,'yes','free','free'
      VALID_INPUT,'PAOLA.PRO','WFS_PARAMS',P10,'structure',{dim:[1,1],tags:['wfs_pitch']},'no','free','free'
      VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.WFS_PITCH',P10.WFS_PITCH,'real',0,'yes','++','free'
      VALID_INPUT,'PAOLA.PRO','SO_ANG',P11,'real',0,'no','0+','free'
      VALID_INPUT,'PAOLA.PRO','SO_ORI',P12,'real',0,'no','free','free'
      VALID_INPUT,'PAOLA.PRO','GLAO_WFS',P13,'structure',{dim:[1,2],tags:['type','ang']},'no','free','free'
      VALID_INPUT,'PAOLA.PRO','GLAO_WFS.TYPE',P13.TYPE,'string',0,'no',['full','edge','FULL','EDGE'],'free'
      VALID_INPUT,'PAOLA.PRO','GLAO_WFS.ANG',P13.ANG,'real',0,'no','0+','free'
    endif else begin
      VALID_INPUT,'PAOLA.PRO','DM_PARAMS',P10,'structure',{dim:[1,3],tags:['dm_height','dmtf','actpitch']},'no','free','free'
      VALID_INPUT,'PAOLA.PRO','DM_PARAMS.DM_HEIGHT',P10.DM_HEIGHT,'real',0,'no','free','free'
      VALID_INPUT,'PAOLA.PRO','DM_PARAMS.DMTF',P10.DMTF,'real',2,'yes','free','free'
      VALID_INPUT,'PAOLA.PRO','DM_PARAMS.ACTPITCH',P10.ACTPITCH,'real',0,'yes','++','free'
      VALID_INPUT,'PAOLA.PRO','WFS_PARAMS',P11,'structure',[1,-1],'no','free','free'
      if n_params() eq 20 and n_tags(P11) ne 1 then message,'WFS_PARAMS STRUCTURE VARIABLE INPUT MUST HAVE 1 COMPONENTS'
      if n_params() eq 22 then if (where([1,8,10] eq n_tags(P11)))[0] eq -1 then message,'WFS_PARAMS STRUCTURE VARIABLE INPUT MUST HAVE 1, 8, OR 10 COMPONENTS'
      if n_tags(P11) eq 1 then VALID_INPUT,'PAOLA.PRO','WFS_PARAMS',P11,'structure',{dim:[1,1],tags:['wfs_pitch']},'no','free','free'
      if n_tags(P11) eq 8 then VALID_INPUT,'PAOLA.PRO','WFS_PARAMS',P11,'structure',{dim:[1,8],tags:['mircoating','mirangle','nblenses','extrafilter','wfs_pitch','wfs_ron','sky_radiance','algorithm']},'no','free','free'
      if n_tags(P11) eq 10 then VALID_INPUT,'PAOLA.PRO','WFS_PARAMS',P11,'structure',$
        {dim:[1,10],tags:['mircoating','mirangle','nblenses','extrafilter','wfs_pitch','wfs_ron','sky_radiance','algorithm','wfs_pxfov','wfs_pxsize']},'no','free','free'
      VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.WFS_PITCH',P11.WFS_PITCH,'real',0,'yes','++','free'
      VALID_INPUT,'PAOLA.PRO','SO_ANG',P12,'real',0,'no','0+','free'
      VALID_INPUT,'PAOLA.PRO','SO_ORI',P13,'real',0,'no','free','free'
      VALID_INPUT,'PAOLA.PRO','WFS_INT',P14,'real',0,'yes','++','free'
      VALID_INPUT,'PAOLA.PRO','LAG',P15,'real',0,'no','0+','free'
      VALID_INPUT,'PAOLA.PRO','LOOP_MODE',P16,'string',0,'no',['open','closed','OPEN','CLOSED'],'free'
      VALID_INPUT,'PAOLA.PRO','LOOP_GAIN',P17,'real',0,'no',[0,4.9348],0
      VALID_INPUT,'PAOLA.PRO','GLAO_WFS',P18,'structure',{dim:[1,2],tags:['type','ang']},'no','free','free'
      VALID_INPUT,'PAOLA.PRO','GLAO_WFS.TYPE',P18.TYPE,'string',0,'no',['star','STAR'],'free'
      VALID_INPUT,'PAOLA.PRO','GLAO_WFS.ANG',P18.ANG,'real',[2,2,-1],'no','free','free'
      VALID_INPUT,'PAOLA.PRO','GS_WEIGHT',P19,'real',1,'yes','++','free'
      if n_elements(P19) gt 1 then if n_elements(P19) ne n_elements(P18.ANG[0,*]) then message,'THE NUMBER OF GUIDE STARS IN GS_WEIGHT AND GLAO_WFS.ANG ARE NOT THE SAME'
      if n_params() eq 20 then begin
        VALID_INPUT,'PAOLA.PRO','WFS_NEA',P20,'real',1,'no','0+','free'
        if n_elements(P20) ne n_elements(P18.ANG[0,*]) then $
          message,'THE NUMBER OF GUIDE STARS IN GLAO_WFS.ANG AND IN THE WFS_NEA ARRAY ARE NOT THE SAME - CHECK THE NUMBER OF ELEMENTS IN WFS_NEA'
      endif
      if n_params() eq 22 then begin
        if n_tags(P11) eq 1 then message,'IF YOU WANT TO EVALUATE THE EFFECT OF THE NGS MAGNITIUDE, YOU NEED TO SET THE WFS PARAMETERS TOO: SEE INPUT WFS_PARAMS IN THE USER MANUAL (FILE HEADER) '
        VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.MIRCOATING',P11.MIRCOATING,'string',1,'yes',['al','ag','au','Al','Ag','Au','AL','AG','AU','aL','aG','aU'],'free'
        VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.MIRANGLE',P11.MIRANGLE,'real',1,'yes',[0,90],90
        VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.NBLENSES',P11.NBLENSES,'integer',0,'yes','0+','free'
        if size(P11.EXTRAFILTER,/type) ne 7 and size(P11.EXTRAFILTER,/type) ne 8 then message,'WFS_PARAMS.EXTRAFILTER SHALL BE EITHER A STRUCTURE VARIABLE OR ''NO'' '
        if size(P11.EXTRAFILTER,/type) eq 7 then VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.EXTRAFILTER',P11.EXTRAFILTER,'string',0,'no','no','free'
        if size(P11.EXTRAFILTER,/type) eq 8 then begin
          VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.EXTRAFILTER',P11.EXTRAFILTER,'structure',{dim:[1,2],tags:['lambda','tauext']},'no','free','free'
          VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.EXTRAFILTER.LAMBDA',P11.EXTRAFILTER.LAMBDA,'real',1,'yes',[0.3d,5.d],'free'
          VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.EXTRAFILTER.TAUEXT',P11.EXTRAFILTER.TAUEXT,'real',1,'yes',[0,1],'free'
        endif
        VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.WFS_RON',P11.WFS_RON,'real',0,'no','0+','free'
        VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.SKY_RADIANCE',P11.SKY_RADIANCE,'real',0,'no','0+','free'
        if P11.SKY_RADIANCE ne 0 then begin
           P11.SKY_RADIANCE = 0
           print,'WFS_PARAMS.SKY_RADIANCE OPTION REMOVED FOR NOW, NEED BETTER MODELING'
        endif
        VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.ALGORITHM',P11.ALGORITHM,'string',0,'no',['4q','4Q','cg','Cg','cG','CG','gs','Gs','gS','GS'],'free'
        if strlowcase(P11.ALGORITHM) ne '4q' then begin
          VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.WFS_PXFOV',P11.WFS_PXFOV,'integer',0,'no','++','free'
          if P11.WFS_PXFOV lt 6 then message,'NUMBER OF PIXEL ACROSS THE WFS LENSLET FOV MUST BE > OR  =  TO 6'
          if P11.WFS_PXFOV mod 2 eq 1 then message,'NUMBER OF PIXEL ACROSS THE WFS LENSLET FOV MUST BE EVEN'
          VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.WFS_PXSIZE',P11.WFS_PXSIZE,'real',0,'yes','++','free'
        endif
        VALID_INPUT,'PAOLA.PRO','NGS_MAG',P20,'real',1,'no','free','free'
        VALID_INPUT,'PAOLA.PRO','FILTER',P21,'string',[1,-1],'no',['u','b','v','r','i','j','h','k','l','m','n','q','U','B','V','R','I','J','H','K','L','M','N','Q'],'free'
        VALID_INPUT,'PAOLA.PRO','NGS_TEM',P22,'real',1,'no','++','free'
        tmp = [n_elements(P20),n_elements(P21),n_elements(P22)]
        if min(tmp) ne max(tmp) then message,'INPUTS NGS_MAG, FILTER AND NGS_TEM DO NOT HAVE THE SAME NUMBER OF ELEMENTS'
        if min(tmp) ne n_elements(P18.ANG[0,*]) then $
          message,'THE NUMBER OF GUIDE STARS IN GLAO_WFS.ANG AND IN THE NGS PARAMETERS ARE NOT THE SAME - CHECK THE LENGHT OF INPUTS NGS_MAG, FILTER AND NGS_TEM'
      endif
    endelse
  endif

  ;----------------
  ;CHECKING OPTIONS
  ;----------------
  if size(EEW_DISC,/type) ne 0 then VALID_INPUT,'PAOLA.PRO','EEW_DISC',EEW_DISC,'real',0,'no','++','free'
  if size(EEW_SQUA,/type) ne 0 then VALID_INPUT,'PAOLA.PRO','EEW_SQUA',EEW_SQUA,'real',0,'no','++','free'
  if size(EEW_SLIT,/type) ne 0 then VALID_INPUT,'PAOLA.PRO','EEW_SLIT',EEW_SLIT,'real',0,'no','++','free'
  if keyword_set(DISPERSION) and strlowcase(strtrim(MODE,2)) ne 'ngs' then message,'DISPERSION EFFECTS ONLY AVAILABLE IN NGS MODE FOR NOW'
  if keyword_set(DISPERSION) then VALID_INPUT,'PAOLA.PRO','DISPERSION',DISPERSION,'real',0,'no','++','free'
  if keyword_set(OPTIMIZE_WFS_INT) and strlowcase(strtrim(MODE,2)) eq 'seli' then message,'OPTIMIZE_WFS_INT KEYWORD AND SEEING LIMITED MODE ARE NOT COMPATIBLE'
  if keyword_set(OPTIMIZE_WFS_INT) and strlowcase(strtrim(MODE,2)) eq 'ngs' then message,'OPTIMIZE_WFS_INT KEYWORD AND NGS MODE ARE NOT COMPATIBLE'
  if keyword_set(OPTIMIZE_LOOP_FREQ) and strlowcase(strtrim(MODE,2)) eq 'seli' then message,'OPTIMIZE_LOOP_FREQ KEYWORD AND SEEING LIMITED MODE ARE NOT COMPATIBLE'
  if keyword_set(OPTIMIZE_LOOP_FREQ) and strlowcase(strtrim(MODE,2)) eq 'glao' then message,'OPTIMIZE_LOOP_FREQ KEYWORD AND GLAO MODE ARE NOT COMPATIBLE'
  if keyword_set(OPTIMIZE_LOOP_GAIN) and strlowcase(strtrim(MODE,2)) eq 'seli' then message,'OPTIMIZE_LOOP_GAIN KEYWORD AND SEEING LIMITED MODE ARE NOT COMPATIBLE'
  if keyword_set(OPTIMIZE_ALL) and strlowcase(strtrim(MODE,2)) ne 'ngs' then message,'OPTIMIZE_ALL KEYWORD AND SEEING LIMITED OR GLAO MODES ARE NOT COMPATIBLE'
  if keyword_set(OPTIMIZE_LOOP_FREQ) and n_params() eq 18 then message,'OPTIMIZE_LOOP_FREQ MAKES NO SENSE WHEN THE WFS NOISE EQUIVALENT ANGLE IS DEFINED (SEE INPUT WFS_NEA)'
  if keyword_set(OPTIMIZE_LOOP_GAIN) and n_params() eq 18 then message,'OPTIMIZE_LOOP_FREQ MAKES NO SENSE WHEN THE WFS NOISE EQUIVALENT ANGLE IS DEFINED (SEE INPUT WFS_NEA)'
  if keyword_set(OPTIMIZE_ALL) and n_params() eq 18 then message,'OPTIMIZE_LOOP_FREQ MAKES NO SENSE WHEN THE WFS NOISE EQUIVALENT ANGLE IS DEFINED (SEE INPUT WFS_NEA)'
  if size(MAX_LOOP_FREQ,/type) ne 0 and not (keyword_set(OPTIMIZE_ALL) or keyword_set(OPTIMIZE_LOOP_FREQ)) then message,'MAX_LOOP_FREQ OPTIONAL INPUT ONLY REQUIRED IN OPTIMZE ALL OR OPTIMIZE FREQ MODES'
  if size(MAX_LOOP_FREQ,/type) ne 0 then VALID_INPUT,'PAOLA.PRO','MAX_LOOP_FREQ',MAX_LOOP_FREQ,'real',0,'no','++','free'
  if (strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) and (keyword_set(OPTIMIZE_LOOP_FREQ)+keyword_set(OPTIMIZE_LOOP_GAIN)+keyword_set(OPTIMIZE_ALL) gt 1b) then $
    message,'OPTIMIZATION OF LOOP GAIN AND/OR WFS INTEGRATION TIME NOT POSSIBLE IN GLAO <FULL> AND <EDGE> MODES'
  if strlowcase(strtrim(MODE,2)) ne 'seli' then begin
    if keyword_set(OPTIMIZE_LOOP_FREQ) or keyword_set(OPTIMIZE_LOOP_GAIN) or keyword_set(OPTIMIZE_ALL) then if mean(P09^2) eq 0 then $
       message,'YOU CANNOT OPTIMIZE THE LOOP WHEN THE WIND SPEED IS ZERO: THE SERVO-LAG IS NATURALLY NULL, AND THE WFS INTEGRATION TIME CAN BE SET TO ANYTHING ARBITRARILY LARGE '+$
               'TO CANCEL THE NOISE, THEREFORE THERE IS NO OPTIMAL VALUE OF THE INTEGRATION TIME, AND FOR THE LOOP GAIN NEITHER. DO AT RETALL AND CHANGE THE PARAMETERS.'
    if strlowcase(strtrim(MODE,2)) ne 'glao' then begin
      if keyword_set(OPTIMIZE_LOOP_GAIN) and strlowcase(P18) eq 'open' then message,'OPTIMIZE_LOOP_GAIN KEYWORD NOT COMPATIBLE WITH OPEN LOOP MODE'
      if keyword_set(OPTIMIZE_ALL) and strlowcase(P18) eq 'open' then message,'OPTIMIZE_ALL KEYWORD NOT COMPATIBLE WITH OPEN LOOP MODE'
    endif
    if keyword_set(OPTIMIZE_LOOP_FREQ)+keyword_set(OPTIMIZE_LOOP_GAIN)+keyword_set(OPTIMIZE_ALL) gt 1b then $
    message,'AMBIGUOUS KEYWORDS: PLEASE MAKE A CHOICE BETWEEN /OPTIMIZE_LOOP_FREQ - /OPTIMIZE_LOOP_GAIN - /OPTIMIZE_ALL (ONLY ONE CHOICE)'
    if (keyword_set(OPTIMIZE_LOOP_FREQ) or keyword_set(OPTIMIZE_LOOP_GAIN) or keyword_set(OPTIMIZE_ALL)) and strlowcase(strtrim(MODE,2)) eq 'ngs' then if n_params() ne 21 then $
      message,'NUMBER OF INPUTS IS INCOMPATIBLE WITH THE KEYWORDS /OPTIMIZE_LOOP_FREQ /OPTIMIZE_LOOP_GAIN AND /OPTIMIZE_ALL'
    if (keyword_set(OPTIMIZE_LOOP_FREQ) or keyword_set(OPTIMIZE_LOOP_GAIN) or keyword_set(OPTIMIZE_ALL)) and strlowcase(strtrim(MODE,2)) eq 'glao' then if n_params() ne 22 then $
      message,'NUMBER OF INPUTS IS INCOMPATIBLE WITH THE KEYWORDS /OPTIMIZE_LOOP_FREQ /OPTIMIZE_LOOP_GAIN AND /OPTIMIZE_ALL'
    if not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then if size(P15,/type) ne 0 then if P15 eq -1 and not (keyword_set(OPTIMIZE_LOOP_FREQ) or keyword_set(OPTIMIZE_ALL)) then $
      message,'DO YOU WANT WFS INTEGRATION TIME OPTIMIZATION ? IN THIS CASE, SET THE KEYWORD /OPTIMIZE_LOOP_FREQ, OTHERWISE, GIVE A POSITIVE VALUE TO WFS_INT'
  endif
  if keyword_set(ANTI_ALIAS) and strlowcase(strtrim(MODE,2)) eq 'seli' then message,'WFS ANTI-ALIASING KEYWORD /ANTI_ALIAS AND SEEING LIMITED MODE ARE NOT COMPATIBLE'
  if size(FITSCODE,/type) ne 0 then begin
    VALID_INPUT,'PAOLA.PRO','FITSCODE',FITSCODE,'string',0,'no','free','free'
    if not keyword_set(PSF) and not keyword_set(OTF) and not keyword_set(SF) and not keyword_set(PSD) then $
      message,'WHEN FITSCODE IS SET, YOU MUST CHOOSE AT LEAST ONE OF THE OPTIONS /PSF - /OTF - /SF - /PSD'
  endif
  if size(LOGCODE,/type) ne 0 then VALID_INPUT,'PAOLA.PRO','LOGCODE',LOGCODE,'string',0,'no','free','free'
  if n_params() eq 6 and keyword_set(SCINTILLATION) then message,'IF SCINTILLATION KEYWORD IS SET, LAYERS ALTITUDE AND CN2 PROFILE INPUTS ARE NEEDED'
  if keyword_set(ONLY_PSD) and (keyword_set(OTF) or keyword_set(PSF) or keyword_set(FWHM_ANALYSIS) or keyword_set(EE_ANALYSIS_DISC) or keyword_set(EEW_DISC) or $
                                keyword_set(EEW_SQUA) or keyword_set(EE_ANALYSIS_SQUA) or keyword_set(EE_ANALYSIS_SLIT) or keyword_set(EEW_SLIT) or keyword_set(SF)) then $
    message,'/ONLY_PSD OPTION IS NOT COMPATIBLE WITH OPTIONS PSF, OTF, STRUCTURE FUNCTION, ENCIRCLED/ENSQUARED/ENSLITED ENERGY ANALYSIS AND FWHM ANALYSIS'
  if keyword_set(PSF) and (keyword_set(X_PSF) or keyword_set(Y_PSF)) then message,'/PSF NOT COMPATIBLE WITH /X_PSF OR /Y_PSF'
  if keyword_set(FWHM_ANALYSIS) and (keyword_set(X_PSF) or keyword_set(Y_PSF)) then message,'/FWHM_ANALYSIS KEYWORD NOT COMPATIBLE WITH /X_PSF OR /Y_PSF'
  if keyword_set(FWHM_ANALYSIS) and not keyword_set(PSF) then message,'TO COMPUTE THE PSF FWHM, PAOLA NEEDS TO BUILD THE PSF, SO IF YOU SET /FWHM_ANALYSIS YOU ALSO NEED TO SET THE KEYWORD /PSF'
  if (keyword_set(X_PSF) or keyword_set(Y_PSF)) and (keyword_set(EE_ANALYSIS_DISC) or keyword_set(EEW_DISC) or keyword_set(EEW_SQUA) or keyword_set(EE_ANALYSIS_SQUA) or $
      keyword_set(EE_ANALYSIS_SLIT) or keyword_set(EEW_SLIT) or keyword_set(SF)) then $
    message,'/X_PSF & /Y_PSF OPTIONS ARE NOT COMPATIBLE WITH OPTIONS ENCIRCLED/ENSQUARED/ENSLITED ENERGY ANALYSIS'
  if strlowcase(strtrim(MODE,2)) eq 'seli' and keyword_set(FRFFT) then message,'FRACTIONAL FFT MODE (/FRFFT KEYWORD) IS NOT COMPATIBLE WITH SEEING LIMITED MODE (SELI)'
  if strlowcase(strtrim(MODE,2)) ne 'seli' then if DIM.N_LF_PADDED ge 6000 and not keyword_set(FRFFT) then begin
    print,'################################################'
    print,'WARNING WARNING WARNING WARNING WARNING WARNING'
    print,'################################################'
    print,''
    print,'The matrix size of the low spatial frequency PSD is very large (>6000): to speed up the computation,'
    print,'you should consider using the optional fractional FFT algorithm. Note that this requires the installation'
    print,'of the fractional FFT and FFTW C-code libraries on your computer. The IDL code which is calling the C codes'
    print,'is available in the PAOLA distribution - see fractional_fft.pro. To get the C codes and the installation'
    print,'procedure, send me a message at laurent.jolissaint@heig-vd.ch.com'
    print,''
    print,'################################################'
    print,'WARNING WARNING WARNING WARNING WARNING WARNING'
    print,'################################################'
  endif
  if keyword_set(FRFFT) then if file_which('fractional_fft.pro') eq '' then begin ; test if the fractional FFT program is available
    FRFFT = 0b
    print,'################################################'
    print,'WARNING WARNING WARNING WARNING WARNING WARNING'
    print,'################################################'
    print,''
    print,'FRACTIONAL FFT IDL CODE (fractional_fft.pro) NOT FOUND. REVERTING TO THE CLASSICAL IDL-FFT CODE.'
    print,''
    print,'################################################'
    print,'WARNING WARNING WARNING WARNING WARNING WARNING'
    print,'################################################'
  endif

  ;============================================================================== 
  ;========================= END OF ARGUMENTS CHECK =============================== 
  ;============================================================================== 

  ;-----------------------------------
  ;GENERAL SETTINGS
  ;-----------------------------------
  rad2asec = 3600.d*180.d/!dpi
  asec2rad = 1.d/rad2asec

  ;-----------------------------------
  ;SETTING VARIABLES ACCORDING TO MODE
  ;-----------------------------------
  if strlowcase(strtrim(MODE,2)) eq 'seli' and n_params() eq 9 then WIND = double(P09)
  if strlowcase(strtrim(MODE,2)) eq 'ngs' then begin
    WIND = double(P09)
    DM_HEIGHT = double(P10.DM_HEIGHT)
    ACTPITCH = double(P10.ACTPITCH)
    DMTF = double(P10.DMTF)
    WFS_PITCH = double(P11.WFS_PITCH)
    GS_ANG = double(P12)
    GS_ORI = double(P13)
    WFS_INT = double(P14)
    LAG = double(P15)
    LOOP_MODE = P16
    if P17 ne -1 then LOOP_FREQ = double(P17)
    if P17 eq -1 then LOOP_FREQ = 1000.d/WFS_INT
    LOOP_GAIN = double(P18)
    if n_params() eq 21 then begin
      MIRCOATING = P11.MIRCOATING
      MIRANGLE = P11.MIRANGLE
      if n_elements(MIRANGLE) eq 1 then if MIRANGLE eq -1 then begin
        MIRCOATING = ['Al','Al','Al','Ag','Ag','Ag']
        MIRANGLE = [0,0,0,0,0,0] 
      endif
      NBLENSES = P11.NBLENSES
      EXTRAFILTER = P11.EXTRAFILTER
      WFS_RON = double(P11.WFS_RON)
      SKY_RADIANCE = double(P11.SKY_RADIANCE)
      ALGORITHM = P11.ALGORITHM
      if strlowcase(ALGORITHM) eq 'cg' or strlowcase(ALGORITHM) eq 'gs' then begin
        WFS_PXFOV = long(P11.WFS_PXFOV)
        WFS_PXSIZE = double(P11.WFS_PXSIZE)
      endif
      NGS_MAG = double(P19)
      FILTER = strlowcase(strtrim(P20,2))
      NGS_TEM = double(P21)
    endif else WFS_NEA = double(P19)
  endif
  if strlowcase(strtrim(MODE,2)) eq 'glao' then begin
    if n_params() eq 13 then begin
      DM_HEIGHT = double(P09.DM_HEIGHT)
      ACTPITCH = double(P09.ACTPITCH)
      DMTF = double(P09.DMTF)
      WFS_PITCH = double(P10.WFS_PITCH)
      GS_ANG = double(P11)
      GS_ORI = double(P12)
      GLAO_WFS = P13
    endif
    if n_params() eq 20 or n_params() eq 22 then begin
      WIND = double(P09)
      DM_HEIGHT = double(P10.DM_HEIGHT)
      ACTPITCH = double(P10.ACTPITCH)
      DMTF = double(P10.DMTF)
      WFS_PITCH = double(P11.WFS_PITCH)
      GS_ANG = double(P12)
      GS_ORI = double(P13)
      WFS_INT = double(P14)
      LAG = double(P15)
      LOOP_MODE = P16
      if strlowcase(LOOP_MODE) eq 'closed' then message,'CLOSED LOOP NOT IMPLEMENTED IN GLAO MODE, YET. SORRY. SET LOOP MODE TO ''OPEN'' FOR NOW.'
      LOOP_GAIN = double(P17)
      GLAO_WFS = P18
      GS_WEIGHT = double(P19)
    endif
    if n_params() eq 20 then WFS_NEA = double(P20)
    if n_params() eq 22 then begin
      MIRCOATING = P11.MIRCOATING
      MIRANGLE = P11.MIRANGLE 
      if n_elements(MIRANGLE) eq 1 then if MIRANGLE eq -1 then begin
        MIRCOATING = ['Al','Al','Al','Ag','Ag','Ag']
        MIRANGLE = [0,0,0,0,0,0] 
      endif
      NBLENSES = P11.NBLENSES
      EXTRAFILTER = P11.EXTRAFILTER
      WFS_RON = double(P11.WFS_RON)
      SKY_RADIANCE = double(P11.SKY_RADIANCE)
      ALGORITHM = P11.ALGORITHM
      if strlowcase(ALGORITHM) eq 'cg' or strlowcase(ALGORITHM) eq 'gs' then begin
        WFS_PXFOV = long(P11.WFS_PXFOV)
        WFS_PXSIZE = double(P11.WFS_PXSIZE)
      endif
      NGS_MAG = double(P20)
      FILTER = strlowcase(strtrim(P21,2))
      NGS_TEM = double(P22)
    endif
  endif

  ;---------------------------------------------------
  ;SEPARE WFS_INT LOOP_FREQ LOOP_GAIN INPUTS FROM OUTSIDE WORLD,
  ;BECAUSE OPTIMIZATION CHANGES INPUT VALUES
  ;---------------------------------------------------
  if size(WFS_INT,/type) ne 0 then gsint = double(WFS_INT)
  if size(LOOP_FREQ,/type) ne 0 then loopfreq = double(LOOP_FREQ)
  if size(LOOP_GAIN,/type) ne 0 then loopgain = double(LOOP_GAIN)
  if strlowcase(strtrim(MODE,2)) eq 'glao' then if strlowcase(GLAO_WFS.TYPE) ne 'star' then gsint = 0

  ;---------------------------------
  ;FRIED PARAMETER VERSUS WAVELENGTH
  ;---------------------------------
  W0rad = W0*asec2rad
  r0500 = 0.98d*0.5d-6/W0rad*cos(double(ZA)/180*!dpi)^(3.d/5) ; Fried's r0 @ 500 nm @ z = ZA
  r0LAM = r0500*(DIM.LAMBDA/0.5)^(1.2d)

  ;-------------------------------------------
  ;SEEING ANGLE ATTENUATION DUE TO OUTER SCALE
  ;-------------------------------------------
  if L0 eq -1 then aos = 1
  if L0 ne -1 then aos = (ATTOS_FWHM(double(TSC.DEXTMAX)/r0LAM,double(TSC.DEXTMAX)/L0))[0]

  ;----------------------------------------------------------------------------
  ;REMOVING LAYERS FOR WHICH CN2  =  0 AND ADAPTING LAYERS HEIGHT TO ZENITH ANGLE
  ;----------------------------------------------------------------------------
  if n_params() gt 6 then begin
    tmp = DISTCN2/total(DISTCN2) ; normalization of cn2 profile distribution
    w = where(tmp gt 0)
    if w[0] ne -1 then begin
      newdistcn2 = tmp[w]
      newheight = HEIGHT[w]/cos(double(ZA)/180*!dpi)
      if not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then begin
        newwindx = (reform(WIND[*,0]))[w]
        newwindy = (reform(WIND[*,1]))[w]
        newwind = [[newwindx],[newwindy]]
      endif
    endif else begin
      newdistcn2 = tmp
      newheight = HEIGHT/cos(double(ZA)/180*!dpi)
      if not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then newwind = WIND
    endelse
  endif

  ;-----------------------------------------
  ;FRIED PARAMETER AND Cn2*dh FOR EACH LAYER
  ;-----------------------------------------
  if n_params() gt 6 then begin
    cn2dh_i = (500d-9/2/!dpi)^2/0.423d*r0500^(-5.d/3)*newdistcn2
    r0500_i = r0500*newdistcn2^(-3.d/5)
  endif

  ;------------------------------------------
  ;LAYERS MEAN ALTITUDE AND ISOPLANATIC ANGLE
  ;------------------------------------------
  if n_params() gt 6 then begin
    meanalti = total(newdistcn2*abs(newheight)^(5.d/3))^(3.d/5)
    if meanalti ne 0 then anisoang = 0.314d/meanalti*r0500*(DIM.LAMBDA/0.5d)^(1.2d)*rad2asec
    if meanalti eq 0 then anisoang = 0.d
  endif

  ;--------------------------------------------------
  ;LAYERS MEAN VELOCITY, TURBULENT PHASE LIFE TIME
  ;--------------------------------------------------
  if n_params() gt 6 then if not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then begin
    if (size(newwind))[0] ne 0 then begin
      meanwind = total(newdistcn2*sqrt(newwind[*,0]^2+newwind[*,1]^2)^(5.d/3))^(3.d/5)
      if meanwind ne 0 then timescale = 0.314d*r0500*(DIM.LAMBDA/0.5d)^(1.2d)/meanwind*1d3
      if meanwind eq 0 then timescale = 1d99
    endif else begin
      meanwind = newwind
      if meanwind ne 0 then timescale = 0.314d*r0500*(DIM.LAMBDA/0.5d)^(1.2d)/meanwind*1d3
      if meanwind eq 0 then timescale = 1d99
    endelse
    tau0 = timescale*(DIM.LAMBDA/0.5d)^(-1.2d)
  endif

  ;----------------------
  ;ACTUATOR AND WFS PITCH
  ;----------------------
  if strlowcase(strtrim(MODE,2)) ne 'seli' then begin
    ACTPITCHnom = r0LAM ; NOMINAL PITCH  =  r0 @ <lambda>
    if ACTPITCH eq -1 then ACTPITCH = ACTPITCHnom ; THE DEFAULT VALUE (-1) IS SET TO THE NOMINAL VALUE ACTPITCHnom
    dm_nactnom = TSC.DEXTMAX/ACTPITCHnom
    dm_nact = double(TSC.DEXTMAX)/ACTPITCH
    if WFS_PITCH eq -1 then WFS_PITCH = ACTPITCH
  endif

  ;---------------------------------
  ;PUPIL PLANE COORDINATE RADIUS [m]
  ;---------------------------------
  xpcoohf = (COOGRID(DIM.N_OTF,DIM.N_OTF,SCALE = DIM.N_OTF/2*DIM.DXP,/FT,/RADIUS)).r

  ;-----------------------------------
  ;PUPIL PLANE SPATIAL FREQUENCY [1/m]
  ;-----------------------------------
  if strlowcase(strtrim(MODE,2)) ne 'seli' then fpcoolf = (COOGRID(DIM.N_LF,DIM.N_LF,SCALE = (DIM.N_LF-1)/2*DIM.DFP_LF,/COO_X)).x
  fpcoohf = (COOGRID(DIM.N_OTF,DIM.N_OTF,SCALE = DIM.N_OTF/2*DIM.DFP,/FT,/COO_X)).x

  ;----------------------------
  ;DM SPATIAL TRANSFER FUNCTION
  ;----------------------------
  if strlowcase(strtrim(MODE,2)) ne 'seli' then begin
    if (size(DMTF))[0] eq 0 then FTMD = double(abs(fpcoolf) le (1+1e-3)*0.5d/ACTPITCH and abs(rotate(fpcoolf,1)) le (1+1e-3)*0.5d/ACTPITCH)
    if (size(DMTF))[0] eq 0 then if keyword_set(RADDMTF) then FTMD = sqrt(fpcoolf^2+rotate(fpcoolf,1)^2) le (1+1e-3)*0.5d/ACTPITCH
    if (size(DMTF))[0] eq 2 then FTMD = DMTF
  endif

  ;---------------------------------
  ;INIT VARIANCE OF PHASE COMPONENTS
  ;---------------------------------
  var = dblarr(4)

  ;------------------------
  ;DISPERSION DEFAULT VALUE
  ;------------------------
  if not keyword_set(DISPERSION) then DISPERSION = 1

  ;========================================================================== 
  ;====================== END OF SETTINGS SECTION ============================= 
  ;========================================================================== 

  ;========================================================================== 
  ;===================== PHASE POWER SPECTRUM, OTF & PSF ====================== 
  ;========================================================================== 

  ;---------------
  ;RUNNING MESSAGE
  ;---------------
  if keyword_set(INFO) then begin
    if strlowcase(strtrim(MODE,2)) eq 'seli' then print,'...RUNNING SEEING LIMITED MODE......'
    if strlowcase(strtrim(MODE,2)) eq  'ngs' then print,'...RUNNING NATURAL GUIDE STAR MODE......'
    if strlowcase(strtrim(MODE,2)) eq 'glao' then print,'...RUNNING GROUND LAYER ADAPTIVE OPTICS MODE......'
    print,''
  endif

  ;----------------------------------------
  ;SINGLE CONJUGATE NATURAL GUIDE STAR MODE
  ;----------------------------------------
  if strlowcase(strtrim(MODE,2)) eq 'ngs' then begin
    ; get WFS central wavelength
    if size(WFS_NEA,/type) eq 0 then begin
       tmp = NBPHOTONS(WFS_PITCH^2,{mircoating:MIRCOATING,mirangle:MIRANGLE},NBLENSES,ZA,NGS_MAG,FILTER,NGS_TEM,gsint,-1,EXTRAFILTER = EXTRAFILTER)
       wfslam = tmp.lam
    endif
    ; LOOP FREQUENCY AND/OR LOOP GAIN OPTIMIZATION [OPTION]
    if keyword_set(OPTIMIZE_LOOP_FREQ) or keyword_set(OPTIMIZE_LOOP_GAIN) or keyword_set(OPTIMIZE_ALL) then begin
      if keyword_set(INFO) then begin
        print,' ================================================ '
        if keyword_set(OPTIMIZE_LOOP_FREQ) then print,'LOOP FREQUENCY OPTIMIZATION...'
        if keyword_set(OPTIMIZE_LOOP_GAIN) then print,'LOOP GAIN OPTIMIZATION...'
        if keyword_set(OPTIMIZE_ALL)       then print,'LOOP FREQUENCY & LOOP GAIN OPTIMIZATION...'
        if keyword_set(OPTIMIZE_LOOP_FREQ) then print,' LOOP FREQ [Hz] OPD RMS [nm]'
        if keyword_set(OPTIMIZE_LOOP_GAIN) then print,' LOOP GAIN    OPD RMS [nm]'
        if keyword_set(OPTIMIZE_ALL)       then print,' LOOP FREQ [Hz] LOOP GAIN    OPD RMS [NM]'
      endif
      block = {optmode:0,$
               lm:strlowcase(strtrim(LOOP_MODE,2)),lf:loopfreq,lg:loopgain,n_psd:DIM.N_LF,fpcoo:fpcoolf,dfplam:DIM.DFP_LF,$
               lambda:DIM.LAMBDA,height:newheight,wind:newwind,r0500:r0500,r0500_i:r0500_i,l0:double(L0),za:ZA,$
               mirsv:{mircoating:MIRCOATING,mirangle:MIRANGLE},nblenses:NBLENSES,extrafilter:EXTRAFILTER,wfs_pitch:WFS_PITCH,wfs_ron:WFS_RON,sky_radiance:SKY_RADIANCE,$
               dmtf:FTMD,dmh:DM_HEIGHT,ang:GS_ANG,ori:GS_ORI,wfs_int:gsint,lag:LAG,algor:ALGORITHM,ngs_mag:NGS_MAG,filter:FILTER,ngs_tem:NGS_TEM,$
               ddisp:DISPERSION[0],prt:keyword_set(INFO),alias:1-keyword_set(ANTI_ALIAS),dextmax:TSC.DEXTMAX,timescale:timescale}
      if size(MAX_LOOP_FREQ,/type) ne 0 then block = create_struct(block,'max_loop_freq',MAX_LOOP_FREQ)
      if strlowcase(ALGORITHM) eq 'cg' or strlowcase(ALGORITHM) eq 'gs' then block = create_struct(block,'wfs_fov',WFS_PXFOV,'wfs_pxs',WFS_PXSIZE)
      ; LOOP FREQUENCY OPTIMIZATION [OPTION]
      if keyword_set(OPTIMIZE_LOOP_FREQ) then begin
        block.optmode = 1
        tmp = dblarr(12) ; we first explore the valid dT space to find a rough estimate of the minima location
        if strlowcase(LOOP_MODE) eq 'closed' then begin ; we start in the range 0.1 to 10 ms
          dTmin = LAG*tan(0.7692056d*loopgain-0.054163105d*loopgain^3+0.0030530390d*loopgain^5-9.3010157e-05*loopgain^7+1.1401983d-06*loopgain^9) ; required to avoid loop instabilities
          dTmin = dTmin > 0.09                                                                                                                    ; required to avoid non sense exposure time
          if size(MAX_LOOP_FREQ,/type) ne 0 then dTmin = dTmin > 1d3/MAX_LOOP_FREQ
          dTmax = timescale
          if dTmin ge dTmax then message,'IMPOSSIBLE OPTIMISATION OF WFS EXPOSURE TIME: THE MINIMUM ACCEPTABLE EXPOSURE TIME IS LARGER THAN THE MAXIMUM ACCEPTABLE EXPOSURE TIME'
          dTint = 10^INTX(12,alog10(dTmin),alog10(dTmax))
          for i_tmp = 0,11 do tmp[i_tmp] = OPTIMIZE_NGS(alog10(dTint[i_tmp]))
          COOMIN,tmp,pos,/silent
          pp = AMOEBA(1e-7,P0 = alog10(dTint[pos.i]),scale = alog10(dTint[pos.i])-1,function_name = 'OPTIMIZE_NGS',function_value = fval)
          gsint = (10.d)^pp[0]
          loopfreq = 1000.d/gsint
        endif
        if strlowcase(LOOP_MODE) eq 'open'   then begin
          dTint = 10^INTX(12,-2,1)
          if size(MAX_LOOP_FREQ,/type) ne 0 then dTint = 10^INTX(12,-2,alog10(1d3/MAX_LOOP_FREQ))
          for i_tmp = 0,11 do tmp[i_tmp] = OPTIMIZE_NGS(alog10(dTint[i_tmp]))
          COOMIN,tmp,pos,/silent
          pp = AMOEBA(1e-7,P0 = alog10(dTint[pos.i]),scale = alog10(dTint[pos.i])-1,function_name = 'OPTIMIZE_NGS',function_value = fval)
          gsint = (10.d)^pp[0]
          loopfreq = 1000.d/gsint
        endif
        if keyword_set(INFO) then print,' ======================================= '
        if keyword_set(INFO) then printf,-1,format = '("OPTIMIZED LOOP FREQUENCY [Hz]",d13.6)',loopfreq
        if keyword_set(INFO) then printf,-1,format = '("OPTIMIZED VARIANCE    [RAD^2]",d13.6)',fval[0]
      endif
      ; LOOP GAIN OPTIMIZATION [OPTION]
      if keyword_set(OPTIMIZE_LOOP_GAIN) then begin
        block.optmode = 2
        tmp = dblarr(12) ; we first explore the loop gain space [0.01,g_max] to find a rough estimate of the minima location
        atdt = atan(gsint/LAG)
        g_max = 1.3890312d*atdt-0.017740598d*atdt^3-0.20128461d*atdt^5+0.58676804d*atdt^7-0.34550100d*atdt^9+0.076581176d*atdt^11
        if g_max ge 4.9348022d then g_max = 4.9348d
        g_tmp = INTX(12,0.01,g_max)
        for i_tmp = 0,11 do tmp[i_tmp] = OPTIMIZE_NGS(tan(!dpi*g_tmp[i_tmp]/g_max-!dpi/2))
        COOMIN,tmp,pos,/silent
        pp = AMOEBA(1e-7,P0 = tan(!dpi*g_tmp[pos.i]/g_max-!dpi/2),scale = 0.1*abs(tan(!dpi*g_tmp[pos.i]/g_max-!dpi/2)),function_name = 'OPTIMIZE_NGS',function_value = fval)
        loopgain = (atan(pp[0])+!dpi/2)/!dpi*g_max
        if keyword_set(INFO) then print,' ======================================= '
        if keyword_set(INFO) then printf,-1,format = '("OPTIMIZED LOOP GAIN    [1]",d13.6)',loopgain
        if keyword_set(INFO) then printf,-1,format = '("OPTIMIZED VARIANCE [RAD^2]",d13.6)',fval[0]
      endif    
      ; WFS INT TIME AND LOOP GAIN OPTIMIZATION [OPTION]
      if keyword_set(OPTIMIZE_ALL) then begin
        block.optmode = 3
        tmp = dblarr(8,8) ; we first explore the [dt,loopgain] space [{1e-3:1e+3},{0.01,4.9}] to find a rough estimate of the minima location
        dT_tmp = 10^INTX(8,-3,3)
        g_tmp = 10^INTX(8,-2,alog10(4.9)) ; note that 4.9 is just below the absolute limit for the gain
        for i_tmp = 0,7 do for j_tmp = 0,7 do tmp[i_tmp,j_tmp] = OPTIMIZE_NGS([dT_tmp[i_tmp],g_tmp[j_tmp]])
        COOMIN,tmp,pos,/silent ; here we look for the best initial dT and loop gain values, then we start the real AMOEBA optimization
        pp = AMOEBA(1e-7,P0 = [dT_tmp[pos.i],g_tmp[pos.j]],scale = 0.1*[dT_tmp[pos.i],g_tmp[pos.j]],function_name = 'OPTIMIZE_NGS',function_value = fval)
        gsint = pp[0]
        loopfreq = 1000.d/gsint
        loopgain = pp[1]
        if keyword_set(INFO) then print,' ======================================= '
        if keyword_set(INFO) then printf,-1,format = '("OPTIMIZED LOOP FREQUENCY [Hz]",d13.6)',loopfreq
        if keyword_set(INFO) then printf,-1,format = '("OPTIMIZED LOOP GAIN       [1]",d13.6)',loopgain
        if keyword_set(INFO) then printf,-1,format = '("OPTIMIZED VARIANCE    [RAD^2]",d13.6)',fval[0]
      endif
      if keyword_set(INFO) then print,' ================================================ '
      if keyword_set(INFO) then print,''
      if DIM.DFP_LF gt 0.1/((gsint+LAG)*1e-3*meanwind) then begin
        print,'################################################'
        print,'WARNING WARNING WARNING WARNING WARNING WARNING'
        print,'################################################'
        print,''
        print,'Loop frequency is now too low for the'
        print,'current servo-lag power spectrum: numerical'
        print,'undersampling will occur in the calculation'
        print,'of the servo-lag PSD, and the servo-lag'
        print,'variance will be strongly underestimated.'
        print,'To avoid this:'
        print,'Re-run the function PIXMATSIZE.PRO with the'
        print,'optimized loop frequency value.'
        print,''
        print,'COMPUTATION IS CONTINUING ANYWAY'
        print,''
        print,'################################################'
        print,''
        warning = 'yes'
      endif
    endif

    ;checking for loop stability - single pole & pure integrator
    ;based on Nyquist criterion of Re(OL transfer function) withing ]-1,0[
    nu = 10^INTX(10000,-8,8) ; loop frequency up to sampling frequency
    jw = dcomplex(0,1)*2*!dpi*nu ; s = j*w Laplace transform argument
    nolTF = loopgain*exp(-jw*LAG*1e-3)*(1-exp(-jw*gsint*1e-3))
    dolTF = (jw*gsint*1e-3)^2
    olTF = nolTF/dolTF
    signflip = (shift(imaginary(olTF) gt 0,1)-(imaginary(olTF) gt 0))[1:*]
    signflip = [0,signflip]
    w1 = where(signflip gt 0); indicates where the olTF crosses the Im = 0 line
    if w1[0] eq -1 and (where(imaginary(olTF) lt 0))[0] eq -1 then stability = 'no' ; if the OL TF does not cross the Im = 0 axis and is always > 0 -> instability
    if w1[0] eq -1 and (where(imaginary(olTF) gt 0))[0] eq -1 then stability = 'yes' ; if the OL TF does not cross the Im = 0 axis and is always < 0 -> stability
    if w1[0] ne -1 then begin
      olTF = olTF(w1)
      if min(double(olTF)) le -1  or min(double(olTF)) ge 0 then stability = 'no' ; if the leftmost sign change is outside ]-1,0[ -> instability
      if min(double(olTF)) gt -1 and min(double(olTF)) lt 0 then stability = 'yes' ; if the leftmost sign change is  inside ]-1,0[ -> stability
    endif
    ;
    ; initialisation of low order phase spatial power spectrum
    Wlf = dblarr(DIM.N_LF,DIM.N_LF)
    ;
    ; Fitting error PSD
    tmp = PSD_HFERR_NGS(DIM.N_HF,fpcoohf,r0LAM,double(L0),0,WFS_PITCH,TSC.DEXTMAX)
    Whf = tmp
    var[0] = 0.2313d*(WFS_PITCH/r0LAM)^(5.d/3)
    ;
    ; WFS aliasing PSD
    if not keyword_set(ANTI_ALIAS) then begin
      if gsint eq 0 then tmp = PSD_ALIAS_NGS(DIM.N_LF,fpcoolf,DIM.LAMBDA,r0500_i,double(L0),0,[[0],[0]],WFS_PITCH,    0,  0,loopgain,strtrim(LOOP_MODE,2),FTMD,DISPERSION[0],TSC.DEXTMAX)
      if gsint gt 0 then tmp = PSD_ALIAS_NGS(DIM.N_LF,fpcoolf,DIM.LAMBDA,r0500_i,double(L0),0,  newwind,WFS_PITCH,gsint,LAG,loopgain,strtrim(LOOP_MODE,2),FTMD,DISPERSION[0],TSC.DEXTMAX)
      if keyword_set(PSD) or keyword_set(ONLY_PSD) then Walias = tmp
      Wlf = Wlf+tmp
      var[1] = total(tmp)*DIM.DFP_LF^2
    endif
    ;
    ; aniso error
    if GS_ANG ne 0 and gsint eq 0 and LAG eq 0 then begin
      tmp = PSD_ANISO_SERVO(strtrim(LOOP_MODE,2),loopfreq,loopgain,DIM.N_LF,fpcoolf,DIM.LAMBDA,newheight,dblarr(n_elements(newheight),2),r0500_i,double(L0),0,WFS_PITCH,FTMD,$
          DM_HEIGHT,GS_ANG,GS_ORI,0,0,DISPERSION[0],TSC.DEXTMAX)
      if keyword_set(PSD) or keyword_set(ONLY_PSD) then Wanisoservo = tmp
      Wlf = Wlf+tmp
      var[2] = total(tmp)*DIM.DFP_LF^2
    endif
    ; servo error
    if GS_ANG eq 0 and (gsint ne 0 or LAG ne 0) then begin
      tmp = PSD_ANISO_SERVO(strtrim(LOOP_MODE,2),loopfreq,loopgain,DIM.N_LF,fpcoolf,DIM.LAMBDA,dblarr(n_elements(newheight)),$
                          newwind,r0500_i,double(L0),0,WFS_PITCH,FTMD,0,0,0,gsint,LAG,DISPERSION[0],TSC.DEXTMAX)
      if keyword_set(PSD) or keyword_set(ONLY_PSD) then Wanisoservo = tmp
      Wlf = Wlf+tmp
      var[2] = total(tmp)*DIM.DFP_LF^2
    endif
    ; aniso-servo error
    if GS_ANG ne 0 and (gsint ne 0 or LAG ne 0) then begin
      tmp = PSD_ANISO_SERVO(strtrim(LOOP_MODE,2),loopfreq,loopgain,DIM.N_LF,fpcoolf,DIM.LAMBDA,newheight,newwind,$
                          r0500_i,double(L0),0,WFS_PITCH,FTMD,DM_HEIGHT,GS_ANG,GS_ORI,gsint,LAG,DISPERSION[0],TSC.DEXTMAX)
      if keyword_set(PSD) or keyword_set(ONLY_PSD) then Wanisoservo = tmp
      Wlf = Wlf+tmp
      var[2] = total(tmp)*DIM.DFP_LF^2
    endif
    ; dispersion error
    if GS_ANG eq 0 and gsint eq 0 and LAG eq 0 and DISPERSION ne 1 then begin
      tmp = (1-DISPERSION[0])^2*PSD_LFERR_NGS(DIM.N_LF,fpcoolf,r0LAM,double(L0),0,WFS_PITCH,TSC.DEXTMAX)
      if keyword_set(PSD) or keyword_set(ONLY_PSD) then Wdisp = tmp
      Wlf = Wlf+tmp
      var[2] = total(tmp)*DIM.DFP_LF^2
    endif
    ; WFS noise PSD
    if size(NGS_MAG,/type) ne 0 then begin
      tmp = NBPHOTONS(WFS_PITCH^2,{mircoating:MIRCOATING,mirangle:MIRANGLE},NBLENSES,ZA,NGS_MAG,FILTER,NGS_TEM,gsint,-1,EXTRAFILTER = EXTRAFILTER)
      wfsnph = tmp.nph
      wfslam = tmp.lam
      wfsbdw = tmp.bdw
      wfstau = tmp.aot
      if strlowcase(strtrim(LOOP_MODE,2)) eq   'open' and strlowcase(ALGORITHM) eq '4q' then nea = (SH_NEA(wfsnph,r0500,wfslam,WFS_PITCH,WFS_RON,!dpi*wfslam^3*SKY_RADIANCE*gsint/3.1783133d*1e3,ALGORITHM)).nea
      if strlowcase(strtrim(LOOP_MODE,2)) eq 'closed' and strlowcase(ALGORITHM) eq '4q' then nea = (SH_NEA(wfsnph,r0500,wfslam,WFS_PITCH,WFS_RON,!dpi*wfslam^3*SKY_RADIANCE*gsint/3.1783133d*1e3,ALGORITHM,gsint,LAG,loopgain)).nea
      if strlowcase(strtrim(LOOP_MODE,2)) eq   'open' and strlowcase(ALGORITHM) ne '4q' then nea = (SH_NEA(wfsnph,r0500,wfslam,WFS_PITCH,WFS_RON,!dpi*wfslam^3*SKY_RADIANCE*gsint/3.1783133d*1e3,ALGORITHM,WFS_PXFOV,WFS_PXSIZE)).nea
      if strlowcase(strtrim(LOOP_MODE,2)) eq 'closed' and strlowcase(ALGORITHM) ne '4q' then nea = (SH_NEA(wfsnph,r0500,wfslam,WFS_PITCH,WFS_RON,!dpi*wfslam^3*SKY_RADIANCE*gsint/3.1783133d*1e3,ALGORITHM,WFS_PXFOV,WFS_PXSIZE,gsint,LAG,loopgain)).nea
      tmp = PSD_NOISE_NGS(DIM.N_LF,fpcoolf,WFS_PITCH,FTMD,nea,DIM.LAMBDA,DISPERSION[0],TSC.DEXTMAX)
      if keyword_set(PSD) or keyword_set(ONLY_PSD) then Wnoise = tmp
      Wlf = Wlf+tmp
      var[3] = total(tmp)*DIM.DFP_LF^2
    endif
    if size(WFS_NEA,/type) ne 0 then if WFS_NEA gt 0 then begin
      tmp = PSD_NOISE_NGS(DIM.N_LF,fpcoolf,WFS_PITCH,FTMD,WFS_NEA,DIM.LAMBDA,DISPERSION[0],TSC.DEXTMAX)
      if keyword_set(PSD) or keyword_set(ONLY_PSD) then Wnoise = tmp
      Wlf = Wlf+tmp
      var[3] = total(tmp)*DIM.DFP_LF^2
    endif
  endif

  ;---------
  ;GLAO MODE
  ;---------
  if strlowcase(strtrim(MODE,2)) eq 'glao' then begin
    if keyword_set(OPTIMIZE_LOOP_FREQ) then begin ; WFS integration time optimization
      if keyword_set(INFO) then print,'...GLAO WFS INTEGRATION TIME OPTIMIZATION...'
      if keyword_set(INFO) then print,''
      if keyword_set(INFO) then print,'    WFS INT [MS]     WFE SERVO+NOISE [NM] '
      nstar = n_elements(reform(GLAO_WFS.ANG[0,*]))
      if n_elements(GS_WEIGHT) eq 1 then gsw = dblarr(nstar)+1.d/nstar
      if n_elements(GS_WEIGHT) ne 1 then gsw = GS_WEIGHT
      block = {n_psd:DIM.N_LF,fpcoo:fpcoolf,dfplam:DIM.DFP_LF,lambda:DIM.LAMBDA,height:newheight,$
             wind:newwind,r0500:r0500,r0500_i:r0500_i,l0:double(L0),za:ZA,mirsv:{mircoating:MIRCOATING,mirangle:MIRANGLE},nblenses:NBLENSES,$
             extrafilter:EXTRAFILTER,wfs_pitch:WFS_PITCH,wfs_ron:WFS_RON,sky_radiance:SKY_RADIANCE,$
             dmtf:FTMD,dmh:DM_HEIGHT,ang:GS_ANG,ori:GS_ORI,wfs_int:gsint,lag:LAG,$
             algor:ALGORITHM,ngs_mag:NGS_MAG,filter:FILTER,ngs_tem:NGS_TEM,gs_weight:gsw,glao_wfs:GLAO_WFS,$
             prt:keyword_set(INFO),alias:1-keyword_set(ANTI_ALIAS),dextmax:TSC.DEXTMAX,timescale:timescale}
      if strlowcase(ALGORITHM) eq 'cg' or strlowcase(ALGORITHM) eq 'gs' then block = create_struct(block,'wfs_fov',WFS_PXFOV,'wfs_pxs',WFS_PXSIZE)
      if gsint gt 0 then minF_parabolic,1e-4,gsint,1e4,gsint,fmin,FUNC_NAME = 'OPTIMIZE_GLAO',TOLERANCE = 1e-5
      if keyword_set(INFO) then print,' ======================================= '
      if keyword_set(INFO) then printf,-1,format = '("OPTIMIZED WFS_INT     [MS]",d13.6)',gsint
      if keyword_set(INFO) then printf,-1,format = '("OPTIMIZED VARIANCE [RAD^2]",d13.6)',fmin[0]
      if keyword_set(INFO) then print,' ================================================ '
      if keyword_set(INFO) then print,''
      if DIM.DFP_LF gt 0.1/((gsint+LAG)*1e-3*meanwind) then begin
        print,'################################################'
        print,'WARNING WARNING WARNING WARNING WARNING WARNING'
        print,'################################################'
        print,''
        print,'WFS integration time is now too large for the'
        print,'current servo-lag power spectrum: numerical'
        print,'undersampling will occur in the calculation'
        print,'of the servo-lag PSD, and the servo-lag'
        print,'variance will be strongly underestimated.'
        print,'To avoid this:'
        print,'Re-run the function PIXMATSIZE.PRO with the'
        print,'optimized WFS integration time value.'
        print,''
        print,'COMPUTATION IS CONTINUING ANYWAY'
        print,''
        print,'################################################'
        print,''
        warning = 'yes'
      endif
    endif
    ; Fitting error PSD
    tmp = PSD_HFERR_NGS(DIM.N_HF,fpcoohf,r0LAM,double(L0),0,WFS_PITCH,TSC.DEXTMAX) ; fitting error PSD
    Whf = tmp
    var[0] = 0.2313d*(WFS_PITCH/r0LAM)^(5.d/3)
    Wlf = dblarr(DIM.N_LF,DIM.N_LF)
    ; WFS aliasing PSD
    if not keyword_set(ANTI_ALIAS) and strlowcase(GLAO_WFS.TYPE) eq 'star' then begin
      if size(gsint,/type) eq 0 then begin
        tmp = DISPERSION[0]^2*PSD_ALIAS_NGS_GLAO(DIM.N_LF,fpcoolf,DIM.LAMBDA,newheight,[[0],[0]],r0500,double(L0),0,WFS_PITCH,FTMD,DM_HEIGHT,0,GS_WEIGHT,GLAO_WFS,TSC.DEXTMAX)
      endif else tmp = DISPERSION[0]^2*PSD_ALIAS_NGS_GLAO(DIM.N_LF,fpcoolf,DIM.LAMBDA,newheight,newwind,r0500_i,double(L0),0,WFS_PITCH,FTMD,DM_HEIGHT,gsint,GS_WEIGHT,GLAO_WFS,TSC.DEXTMAX)
      if keyword_set(PSD) or keyword_set(ONLY_PSD) then Walias = tmp
      Wlf = Wlf+tmp
      var[1] = total(tmp)*DIM.DFP_LF^2
    endif
    ; aniso-servo PSD
    if strlowcase(GLAO_WFS.TYPE) eq 'full' then $
      tmp = PSD_ANISO_GLAO_FULL(DIM.N_LF,fpcoolf,DIM.LAMBDA,newheight,r0500_i,double(L0),WFS_PITCH,FTMD,DM_HEIGHT,GS_ANG,GS_ORI,GLAO_WFS,TSC.DEXTMAX)
    if strlowcase(GLAO_WFS.TYPE) eq 'edge' then $
      tmp = PSD_ANISO_GLAO_EDGE(DIM.N_LF,fpcoolf,DIM.LAMBDA,newheight,r0500_i,double(L0),WFS_PITCH,FTMD,DM_HEIGHT,GS_ANG,GS_ORI,GLAO_WFS,TSC.DEXTMAX)
    if strlowcase(GLAO_WFS.TYPE) eq 'star' then begin
      if n_params() eq 13 then tmp = PSD_ANISO_SERVO_GLAO_STAR(DIM.N_LF,fpcoolf,DIM.LAMBDA,newheight,dblarr(n_elements(newheight),2),r0500_i,double(L0),WFS_PITCH,FTMD,$
                                                             DM_HEIGHT,GS_ANG,GS_ORI,0,0,GS_WEIGHT,GLAO_WFS,TSC.DEXTMAX)
      if n_params() ne 13 then tmp = PSD_ANISO_SERVO_GLAO_STAR(DIM.N_LF,fpcoolf,DIM.LAMBDA,newheight,newwind,r0500_i,double(L0),WFS_PITCH,FTMD,$
                                                             DM_HEIGHT,GS_ANG,GS_ORI,gsint,LAG,GS_WEIGHT,GLAO_WFS,TSC.DEXTMAX)
    endif
    if keyword_set(PSD) or keyword_set(ONLY_PSD) then Wanisoservo = tmp
    Wlf = Wlf+tmp
    var[2] = total(tmp)*DIM.DFP_LF^2
    ; WFS noise PSD
    if strlowcase(GLAO_WFS.TYPE) eq 'star' then begin
      tmp = dblarr(DIM.N_LF,DIM.N_LF)
      if size(NGS_MAG,/type) ne 0 then nstar = n_elements(NGS_MAG)
      if size(WFS_NEA,/type) ne 0 then nstar = n_elements(WFS_NEA)
      if n_elements(GS_WEIGHT) eq 1 then gsw = dblarr(nstar)+1.d/nstar
      if n_elements(GS_WEIGHT) ne 1 then gsw = GS_WEIGHT
      ; WFS noise PSD from the GS magnitudes and spectrums and RON
      if size(NGS_MAG,/type) ne 0 then begin
        tmp = NBPHOTONS(WFS_PITCH^2,{mircoating:MIRCOATING,mirangle:MIRANGLE},NBLENSES,ZA,NGS_MAG,FILTER,NGS_TEM,gsint,-1,EXTRAFILTER = EXTRAFILTER)
        wfsnph = tmp.nph
        wfslam = tmp.lam
        wfsbdw = tmp.bdw
        wfstau = tmp.aot
        if strlowcase(ALGORITHM) eq '4q' then nea = (SH_NEA(wfsnph,r0500,wfslam,WFS_PITCH,WFS_RON,!dpi*wfslam^3*SKY_RADIANCE*gsint/3.1783133d*1e3,ALGORITHM)).nea
        if strlowcase(ALGORITHM) eq 'cg' or strlowcase(ALGORITHM) eq 'gs' then nea = (SH_NEA(wfsnph,r0500,wfslam,WFS_PITCH,WFS_RON,!dpi*wfslam^3*SKY_RADIANCE*gsint/3.1783133d*1e3,ALGORITHM,WFS_PXFOV,WFS_PXSIZE)).nea
        tmp = PSD_NOISE_NGS_GLAO(DIM.N_LF,fpcoolf,WFS_PITCH,FTMD,gsw,nea,DIM.LAMBDA,TSC.DEXTMAX)
        if keyword_set(PSD) or keyword_set(ONLY_PSD) then Wnoise = tmp
        Wlf = Wlf+tmp
        var[3] = total(tmp)*DIM.DFP_LF^2
      endif
      ; WFS noise PSD from GS Noise Equivalent Angle
      if size(WFS_NEA,/type) ne 0 then begin
        tmp = PSD_NOISE_NGS_GLAO(DIM.N_LF,fpcoolf,WFS_PITCH,FTMD,gsw,WFS_NEA,DIM.LAMBDA,TSC.DEXTMAX)
        if keyword_set(PSD) or keyword_set(ONLY_PSD) then Wnoise = tmp
        Wlf = Wlf+tmp
        var[3] = total(tmp)*DIM.DFP_LF^2
      endif
    endif
  endif

  ;------------------------------
  ;SCINTILLATION EFFECTS [OPTION]
  ;------------------------------
  if keyword_set(SCINTILLATION) then begin
    w = where(fpcoohf^2+rotate(fpcoohf,1)^2 gt 0)
    tmp = dblarr(DIM.N_OTF,DIM.N_OTF)
    tpf = dblarr(DIM.N_OTF,DIM.N_OTF)
    for indtmp = 0,n_elements(cn2dh_i)-1 do tmp[w] = tmp[w]+(sin(!dpi*1d-6*DIM.LAMBDA*abs(newheight(indtmp))*sqrt(fpcoohf[w]^2+(rotate(fpcoohf,1))[w]^2)^2))^2*cn2dh_i(indtmp)
    tmp[w] = 0.38*(1d-6*DIM.LAMBDA)^(-2)*sqrt(fpcoohf[w]^2+(rotate(fpcoohf,1))[w]^2)^(-11.d/3)*tmp[w] ; amplitude fluctuation spatial power spectrum (relative to 1)
    tpf[w] = (2*beselj(!dpi*TSC.DEXTMAX*sqrt(fpcoohf[w]^2+(rotate(fpcoohf,1))[w]^2),1)/(!dpi*TSC.DEXTMAX*sqrt(fpcoohf[w]^2+(rotate(fpcoohf,1))[w]^2)))^2 ; pupil filt. for scint. index comp.
    sci_index = 4*total(tmp*tpf)*DIM.DFP^2
    if keyword_set(PSD) or keyword_set(ONLY_PSD) then Wamp = tmp
    if not keyword_set(ONLY_PSD) then begin
      B0 = total(tmp)*DIM.DFP^2
      B = double(MATHFT(tmp,DX = DIM.DXP,IC = DIM.N_OTF/2,JC = DIM.N_OTF/2,/INVERSE))
      dphiamp = 2*(B0-B)
      dphiamp(DIM.N_OTF/2,DIM.N_OTF/2) = 0.d
    endif
  endif

  ;----------------------------------
  ;ADAPTIVE OPTICS STRUCTURE FUNCTION
  ;----------------------------------
  ;it is computed from the Fourier transform of the phase PSD
  if strlowcase(strtrim(MODE,2)) ne 'seli' then if not keyword_set(ONLY_PSD) then begin
    ;low frequency structure function
    varLF = total(Wlf)*DIM.DFP_LF^2
    if varLF gt 0 then begin
      if not keyword_set(FRFFT) then begin
        dimPSD = DIM.N_LF_PADDED>DIM.N_HF ; we take the largest dimension
        tmp = dblarr(dimPSD,dimPSD)
        tmp[dimPSD/2-(DIM.N_LF-1)/2:dimPSD/2+(DIM.N_LF-1)/2,dimPSD/2-(DIM.N_LF-1)/2:dimPSD/2+(DIM.N_LF-1)/2] = Wlf ; we insert the LF PSD into the matrix tmp
        tmp = double(MATHFT(tmp,ic = dimPSD/2,jc = dimPSD/2)) ; the phase covariance function B
      endif else begin
        dimPSD = 2*DIM.N_LF>DIM.N_HF ; we take the largest dimension
        tmp = dblarr(dimPSD,dimPSD)
        tmp[dimPSD/2-(DIM.N_LF-1)/2:dimPSD/2+(DIM.N_LF-1)/2,dimPSD/2-(DIM.N_LF-1)/2:dimPSD/2+(DIM.N_LF-1)/2] = Wlf ; we insert the LF PSD into the matrix tmp
        tmp = double(FRACTIONAL_FFT(tmp,dimPSD,DIM.N_LF_PADDED,dimPSD/2,dimPSD/2,dimPSD/2,dimPSD/2)) ; the phase covariance function B
      endelse
      tmp = tmp[dimPSD/2-DIM.N_HF/2:dimPSD/2+DIM.N_HF/2-1,dimPSD/2-DIM.N_HF/2:dimPSD/2+DIM.N_HF/2-1]/tmp[dimPSD/2,dimPSD/2]*varLF ; we keep the DIM.N_HF part and renormalize with variance
      dphiLF = 2*(varLF-tmp)
      dphiLF = dphiLF-min(dphiLF)
    endif else dphiLF = 0
    ;high frequency structure function
    tmp = double(MATHFT(Whf,dx = DIM.DXP,ic = DIM.N_HF/2,jc = DIM.N_HF/2,/inverse))
    tmp = tmp/tmp[DIM.N_OTF/2,DIM.N_OTF/2]*var[0]
    dphiHF = 2*(var[0]-tmp)
    dphi = dphiLF+dphiHF
    if keyword_set(SCINTILLATION) then dphi = dphi+dphiamp
  endif

  ;-------------------
  ;SEEING LIMITED MODE
  ;-------------------
  if strlowcase(strtrim(MODE,2)) eq 'seli' then begin
    if L0 eq -1 then begin
      w = where(fpcoohf^2+rotate(fpcoohf,1)^2 gt 0)
      if not keyword_set(ONLY_PSD) then dphi = 6.883877d*(xpcoohf/r0LAM)^(5.d/3)
      Watm = dblarr(DIM.N_OTF,DIM.N_OTF)
      if w[0] ne -1 then Watm[w] = 0.022896d/r0LAM^(5.d/3)*sqrt(fpcoohf[w]^2+(rotate(fpcoohf,1))[w]^2)^(-11.d/3)
    endif
    if L0 ne -1 then begin
      if not keyword_set(ONLY_PSD) then begin
        w = where(xpcoohf ne 0)
        dphi = dblarr(DIM.N_OTF,DIM.N_OTF)
        if w[0] ne -1 then dphi[w] = 0.171661d*(double(L0)/r0LAM)^(5.d/3)*((1.005635d)-(2*!dpi*xpcoohf[w]/double(L0))^(5.d/6)*$
                                 beselk(2*!dpi*xpcoohf[w]/double(L0),5.d/6))
        if keyword_set(SCINTILLATION) then dphi = dphi+dphiamp
      endif
      Watm = 0.022896d/r0LAM^(5.d/3)*(fpcoohf^2+rotate(fpcoohf,1)^2+1.d/double(L0)^2)^(-11.d/6)
    endif
  endif

  ;---------------------------------------------
  ;INDEPENDANT POST TIP-TILT CORRECTION [OPTION]
  ;---------------------------------------------
  if keyword_set(POST_TIPTILT) then begin

    ;first task is to compute the G-tilt covariance Zernike matrix
    ;note that we compute both tilt components because
    ;the covariance in x and y are not necessarily the same
    ;for instance <a2*a8> can be different from <a3*a7>
    covmat = dblarr(8,8)
    ;we start with the low spatial frequency part
    if strlowcase(strtrim(MODE,2)) ne 'seli' then begin
      Q02 = ZERFT( 2,fpcoolf*TSC.DEXTMAX*0.5,transpose(fpcoolf)*TSC.DEXTMAX*0.5) ; and Q03 = transpose(Q02)
      Q08 = ZERFT( 8,fpcoolf*TSC.DEXTMAX*0.5,transpose(fpcoolf)*TSC.DEXTMAX*0.5) ;     Q07 = transpose(Q08)
      Q16 = ZERFT(16,fpcoolf*TSC.DEXTMAX*0.5,transpose(fpcoolf)*TSC.DEXTMAX*0.5) ;     Q17 = transpose(Q16)
      Q30 = ZERFT(30,fpcoolf*TSC.DEXTMAX*0.5,transpose(fpcoolf)*TSC.DEXTMAX*0.5) ;     Q29 = transpose(Q30)
      ;covariance matrix elements
      covmat[0,0] = total(Q02*Wlf*conj(Q02))*DIM.DFP_LF^2
      covmat[1,1] = total(transpose(Q02)*Wlf*conj(transpose(Q02)))*DIM.DFP_LF^2
      covmat[0,3] = total(Q02*Wlf*conj(Q08))*DIM.DFP_LF^2
      covmat[3,3] = total(Q08*Wlf*conj(Q08))*DIM.DFP_LF^2
      covmat[1,2] = total(transpose(Q02)*Wlf*conj(transpose(Q08)))*DIM.DFP_LF^2
      covmat[2,2] = total(transpose(Q08)*Wlf*conj(transpose(Q08)))*DIM.DFP_LF^2
      covmat[0,4] = total(Q02*Wlf*conj(Q16))*DIM.DFP_LF^2
      covmat[3,4] = total(Q08*Wlf*conj(Q16))*DIM.DFP_LF^2
      covmat[4,4] = total(Q16*Wlf*conj(Q16))*DIM.DFP_LF^2
      covmat[1,5] = total(transpose(Q02)*Wlf*conj(transpose(Q16)))*DIM.DFP_LF^2
      covmat[2,5] = total(transpose(Q08)*Wlf*conj(transpose(Q16)))*DIM.DFP_LF^2
      covmat[5,5] = total(transpose(Q16)*Wlf*conj(transpose(Q16)))*DIM.DFP_LF^2
      covmat[0,7] = total(Q02*Wlf*conj(Q30))*DIM.DFP_LF^2
      covmat[3,7] = total(Q08*Wlf*conj(Q30))*DIM.DFP_LF^2
      covmat[4,7] = total(Q16*Wlf*conj(Q30))*DIM.DFP_LF^2
      covmat[7,7] = total(Q30*Wlf*conj(Q30))*DIM.DFP_LF^2
      covmat[1,6] = total(transpose(Q02)*Wlf*conj(transpose(Q30)))*DIM.DFP_LF^2
      covmat[2,6] = total(transpose(Q08)*Wlf*conj(transpose(Q30)))*DIM.DFP_LF^2
      covmat[5,6] = total(transpose(Q16)*Wlf*conj(transpose(Q30)))*DIM.DFP_LF^2
      covmat[6,6] = total(transpose(Q30)*Wlf*conj(transpose(Q30)))*DIM.DFP_LF^2
    endif
    ;and now the high spatial frequency part
    Q02 = ZERFT( 2,fpcoohf*TSC.DEXTMAX*0.5,transpose(fpcoohf)*TSC.DEXTMAX*0.5) ; and Q03 = transpose(Q02)
    Q08 = ZERFT( 8,fpcoohf*TSC.DEXTMAX*0.5,transpose(fpcoohf)*TSC.DEXTMAX*0.5) ;     Q07 = transpose(Q08)
    Q16 = ZERFT(16,fpcoohf*TSC.DEXTMAX*0.5,transpose(fpcoohf)*TSC.DEXTMAX*0.5) ;     Q17 = transpose(Q16)
    Q30 = ZERFT(30,fpcoohf*TSC.DEXTMAX*0.5,transpose(fpcoohf)*TSC.DEXTMAX*0.5) ;     Q29 = transpose(Q30)
    ;covariance matrix elements HF added to LF ones
    if strlowcase(strtrim(MODE,2)) eq 'seli' then Whf = Watm
    covmat[0,0] = covmat[0,0]+total(Q02*Whf*conj(Q02))*DIM.DFP^2
    covmat[1,1] = covmat[1,1]+total(transpose(Q02)*Whf*conj(transpose(Q02)))*DIM.DFP^2
    covmat[0,3] = covmat[0,3]+total(Q02*Whf*conj(Q08))*DIM.DFP^2
    covmat[3,3] = covmat[3,3]+total(Q08*Whf*conj(Q08))*DIM.DFP^2
    covmat[1,2] = covmat[1,2]+total(transpose(Q02)*Whf*conj(transpose(Q08)))*DIM.DFP^2
    covmat[2,2] = covmat[2,2]+total(transpose(Q08)*Whf*conj(transpose(Q08)))*DIM.DFP^2
    covmat[0,4] = covmat[0,4]+total(Q02*Whf*conj(Q16))*DIM.DFP^2
    covmat[3,4] = covmat[3,4]+total(Q08*Whf*conj(Q16))*DIM.DFP^2
    covmat[4,4] = covmat[4,4]+total(Q16*Whf*conj(Q16))*DIM.DFP^2
    covmat[1,5] = covmat[1,5]+total(transpose(Q02)*Whf*conj(transpose(Q16)))*DIM.DFP^2
    covmat[2,5] = covmat[2,5]+total(transpose(Q08)*Whf*conj(transpose(Q16)))*DIM.DFP^2
    covmat[5,5] = covmat[5,5]+total(transpose(Q16)*Whf*conj(transpose(Q16)))*DIM.DFP^2
    covmat[0,7] = covmat[0,7]+total(Q02*Whf*conj(Q30))*DIM.DFP^2
    covmat[3,7] = covmat[3,7]+total(Q08*Whf*conj(Q30))*DIM.DFP^2
    covmat[4,7] = covmat[4,7]+total(Q16*Whf*conj(Q30))*DIM.DFP^2
    covmat[7,7] = covmat[7,7]+total(Q30*Whf*conj(Q30))*DIM.DFP^2
    covmat[1,6] = covmat[1,6]+total(transpose(Q02)*Whf*conj(transpose(Q30)))*DIM.DFP^2
    covmat[2,6] = covmat[2,6]+total(transpose(Q08)*Whf*conj(transpose(Q30)))*DIM.DFP^2
    covmat[5,6] = covmat[5,6]+total(transpose(Q16)*Whf*conj(transpose(Q30)))*DIM.DFP^2
    covmat[6,6] = covmat[6,6]+total(transpose(Q30)*Whf*conj(transpose(Q30)))*DIM.DFP^2
    ;building the other triangle of the covariance matrix
    tmp = covmat+transpose(covmat) ; the covariance matrix is symetric / diagonal
    for i = 0,7 do tmp[i,i] = 0.5*tmp[i,i]
    covmat = tmp

    ;here we build the Fourier transform of the Zernike products that
    ;appears in the G-tilt structure function formula
    a0202 = ZPZQTOZJ(02,02) & Q0202 = 0
    for i = 0,n_elements(a0202.j)-1 do Q0202 = Q0202+a0202.aj[i]*ZERFT(a0202.j[i],fpcoohf*TSC.DEXTMAX*0.5,transpose(fpcoohf)*TSC.DEXTMAX*0.5)
    a0208 = ZPZQTOZJ(02,08) & Q0208 = 0
    for i = 0,n_elements(a0208.j)-1 do Q0208 = Q0208+a0208.aj[i]*ZERFT(a0208.j[i],fpcoohf*TSC.DEXTMAX*0.5,transpose(fpcoohf)*TSC.DEXTMAX*0.5)
    a0216 = ZPZQTOZJ(02,16) & Q0216 = 0
    for i = 0,n_elements(a0216.j)-1 do Q0216 = Q0216+a0216.aj[i]*ZERFT(a0216.j[i],fpcoohf*TSC.DEXTMAX*0.5,transpose(fpcoohf)*TSC.DEXTMAX*0.5)
    a0230 = ZPZQTOZJ(02,30) & Q0230 = 0
    for i = 0,n_elements(a0230.j)-1 do Q0230 = Q0230+a0230.aj[i]*ZERFT(a0230.j[i],fpcoohf*TSC.DEXTMAX*0.5,transpose(fpcoohf)*TSC.DEXTMAX*0.5)

    ;the alpha coefficients - see theory
    if keyword_set(TILT_ANGLE_STEP) then begin ; this is with quantification of the tilt correction, by TILT_ANGLE_STEP steps
      covxGxG = covmat[0,0]+2*sqrt(2)*covmat[0,3]+2*covmat[3,3]
      covyGyG = covmat[1,1]+2*sqrt(2)*covmat[1,2]+2*covmat[2,2]
      facX = QUANTCOVARIANCEGTILT(TILT_ANGLE_STEP,TSC.DEXTMAX,DIM.LAMBDA,TSC.DEXTMAX/r0LAM)
      facY = QUANTCOVARIANCEGTILT(TILT_ANGLE_STEP,TSC.DEXTMAX,DIM.LAMBDA,TSC.DEXTMAX/r0LAM)
      alpha0202 = 2*facX.f1[0]*(covmat[0,0]+sqrt(2)*covmat[0,3])-facX.f2[0]*covxGxG ; alpha[2,2]
      alpha0303 = 2*facY.f1[0]*(covmat[1,1]+sqrt(2)*covmat[1,2])-facY.f2[0]*covyGyG ; alpha[3,3]
      alpha0208 = facX.f3[0]*(covmat[0,3]+sqrt(2)*covmat[3,3]) ; alpha[2,8]
      alpha0307 = facY.f3[0]*(covmat[1,2]+sqrt(2)*covmat[2,2]) ; alpha[3,7]
    endif else begin
      alpha0202 = covmat[0,0]-2*covmat[3,3]-3*covmat[4,4]-4*covmat[7,7]-2*sqrt(6)*covmat[3,4]-2*sqrt(8)*covmat[3,7]-2*sqrt(12)*covmat[4,7] ; alpha[2,2]
      alpha0208 = covmat[0,3]+sqrt(2)*covmat[3,3]+sqrt(3)*covmat[3,4]+2*covmat[3,7] ; alpha[2,8]
      alpha0303 = covmat[1,1]-2*covmat[2,2]-3*covmat[5,5]-4*covmat[6,6]-2*sqrt(6)*covmat[2,5]-2*sqrt(8)*covmat[2,6]-2*sqrt(12)*covmat[5,6] ; alpha[3,3]
      alpha0307 = covmat[1,2]+sqrt(2)*covmat[2,2]+sqrt(3)*covmat[2,5]+2*covmat[2,6] ; alpha[3,7]
    endelse
    alpha0216 = covmat[0,4]+sqrt(2)*covmat[3,4]+sqrt(3)*covmat[4,4]+2*covmat[4,7] ; alpha[2,16]
    alpha0230 = covmat[0,7]+sqrt(2)*covmat[3,7]+sqrt(3)*covmat[4,7]+2*covmat[7,7] ; alpha[2,30]
    alpha0317 = covmat[1,5]+sqrt(2)*covmat[2,5]+sqrt(3)*covmat[5,5]+2*covmat[5,6] ; alpha[3,17]
    alpha0329 = covmat[1,6]+sqrt(2)*covmat[2,6]+sqrt(3)*covmat[5,6]+2*covmat[6,6] ; alpha[3,29]

    ;sum of the Fourier transform of the products Uij*Ap
    ;where we multiply the Qs by !dpi*(DIAMETER*0.5)^2
    ;to go back in the regular space (not-normalized)
    pupilFT = DISCFT(fpcoohf,transpose(fpcoohf),TSC.DEXTMAX,0) ; pupil Fourier transform
    tmp = alpha0202*2*          double(pupilFT*conj(Q0202)*!dpi*(TSC.DEXTMAX*0.5)^2-Q02*conj(Q02)*(!dpi*(TSC.DEXTMAX*0.5)^2)^2) +$
        alpha0303*2*transpose(double(pupilFT*conj(Q0202)*!dpi*(TSC.DEXTMAX*0.5)^2-Q02*conj(Q02)*(!dpi*(TSC.DEXTMAX*0.5)^2)^2))+$
        alpha0208*2*          double(pupilFT*conj(Q0208)*!dpi*(TSC.DEXTMAX*0.5)^2-Q02*conj(Q08)*(!dpi*(TSC.DEXTMAX*0.5)^2)^2) +$
        alpha0307*2*transpose(double(pupilFT*conj(Q0208)*!dpi*(TSC.DEXTMAX*0.5)^2-Q02*conj(Q08)*(!dpi*(TSC.DEXTMAX*0.5)^2)^2))+$
        alpha0216*2*          double(pupilFT*conj(Q0216)*!dpi*(TSC.DEXTMAX*0.5)^2-Q02*conj(Q16)*(!dpi*(TSC.DEXTMAX*0.5)^2)^2) +$
        alpha0317*2*transpose(double(pupilFT*conj(Q0216)*!dpi*(TSC.DEXTMAX*0.5)^2-Q02*conj(Q16)*(!dpi*(TSC.DEXTMAX*0.5)^2)^2))+$
        alpha0230*2*          double(pupilFT*conj(Q0230)*!dpi*(TSC.DEXTMAX*0.5)^2-Q02*conj(Q30)*(!dpi*(TSC.DEXTMAX*0.5)^2)^2) +$
        alpha0329*2*transpose(double(pupilFT*conj(Q0230)*!dpi*(TSC.DEXTMAX*0.5)^2-Q02*conj(Q30)*(!dpi*(TSC.DEXTMAX*0.5)^2)^2))

    ;inverse FT to get the average structure function
    tmp = double(MATHFT(tmp-mean(tmp),dx = DIM.DXP,ic = DIM.N_OTF/2,jc = DIM.N_OTF/2,/INVERSE))
    sfGtilt = tmp-tmp[DIM.N_OTF/2,DIM.N_OTF/2]

    ;now, we divide with the pupil auto-correlation to follow the
    ;definition of the pupil averaged structure function
    pupilAC = TSCOTF(xpcoohf/TSC.DEXTMAX,0)*!dpi*(TSC.DEXTMAX*0.5)^2 ; pupil autocorrelation (OTF times pupil surface, and we do not need to take into account the actual pupil shape here)
    mask = xpcoohf lt TSC.DEXTMAX-DIM.DXP
    sfGtilt[where(mask eq 1)] = sfGtilt[where(mask eq 1)]/pupilAC[where(mask eq 1)]
    sfGtilt[where(mask eq 0)] = 0

    ;and now we remove the G-tilt SF from the overall SF
    ;to model the G-tilt correction
    dphi = (dphi-sfGtilt)*mask
    tmp = exp(-0.5*dphi)
    tmp[where(1-mask)] = max(tmp)
    COOMIN,abs(tmp),pos,/silent
    kk = 5.d/(6*xpcoohf[pos.i[0],pos.j[0]]^6)
    filterOTF = exp(-xpcoohf^8*kk)

  endif

  ;---------------------------------------------------------------
  ;BUILDING A RESIDUAL PHASE PSD TO BE USED WITH FUNCTION WAVE.PRO
  ;TO GENERATE INSTANTANEOUS CORRECTED PHASE SCREENS
  ;---------------------------------------------------------------
  if keyword_set(WAVE) then begin
    if size(VIBRATION,/type) ne 0 then vibe_wave = VIBRATION
    if strlowcase(strtrim(MODE,2)) eq 'seli' then begin
      dim_wave = DIM.N_OTF
      dfp_wave = DIM.DFP
      dxp_wave = DIM.DXP
      psd_wave = Watm
      if (where(tag_names(TSC) eq 'PUPIL'))[0] ne -1 then pup_wave = TSC.PUPIL gt 0.5
      if (where(tag_names(TSC) eq 'PHASE'))[0] ne -1 then sta_wave = TSC.PHASE
    endif
    if strlowcase(strtrim(MODE,2)) ne 'seli' then begin
      dim_wave = DIM.N_LF_PADDED
      dfp_wave = DIM.DFP_LF
      dxp_wave = DIM.DXP
      psd_wave = dblarr(dim_wave,dim_wave)
      psd_wave[dim_wave/2-(DIM.N_LF-1)/2:dim_wave/2+(DIM.N_LF-1)/2,dim_wave/2-(DIM.N_LF-1)/2:dim_wave/2+(DIM.N_LF-1)/2] = Wlf ; we insert the LF PSD into the matrix
      tmp = congrid(Whf,dim_wave,dim_wave,cubic = -0.5,/inter) ; adapting HF PSD to LF resolution
      fpup4wave = COOGRID(dim_wave,dim_wave,/ft,scale = dim_wave/2*DIM.DFP_LF)
      wtmp = where(abs(fpup4wave.x) le 1.d/(2*WFS_PITCH) and abs(fpup4wave.y) le 1.d/(2*WFS_PITCH)) ; selecting pixels inside LF domain
      tmp[wtmp] = 0 ; removing HF PSD array possibly leaking into LF domain
      psd_wave = psd_wave+tmp ; the full LF+HF PSD
      if (where(tag_names(TSC) eq 'PUPIL'))[0] ne -1 then begin
        pup_wave = bytarr(dim_wave,dim_wave)
        pup_wave[dim_wave/2-DIM.N_OTF/2:dim_wave/2+DIM.N_OTF/2-1,dim_wave/2-DIM.N_OTF/2:dim_wave/2+DIM.N_OTF/2-1] = TSC.PUPIL gt 0.5
      endif
      if (where(tag_names(TSC) eq 'PHASE'))[0] ne -1 then begin
        sta_wave = dblarr(dim_wave,dim_wave)
        sta_wave[dim_wave/2-DIM.N_OTF/2:dim_wave/2+DIM.N_OTF/2-1,dim_wave/2-DIM.N_OTF/2:dim_wave/2+DIM.N_OTF/2-1] = TSC.PHASE
      endif
    endif
  endif

  ;[OPTION] INSTRUMENT VIBRATION EFFECT ON THE OTF
  vibeOTF = 1
  if size(VIBRATION,/type) ne 0 then begin
     vibeSF = STRUCTURE_FUNCTION_JITTER(DIM.N_OTF,DIM.DXP,VIBRATION.SIGMA_TTX,VIBRATION.SIGMA_TTY,VIBRATION.ORIENTATION)
     vibeOTF = exp(-0.5*vibeSF*(2*!dpi/(DIM.LAMBDA*1e-6))^2)
     instvib_asec = sqrt(VIBRATION.SIGMA_TTX^2+VIBRATION.SIGMA_TTY^2)
     instvib_nm = instvib_asec*asec2rad*TSC.DEXTMAX/2*1e9
  endif
  if where(tag_names(TSC) eq 'VIBRATION') ne -1 then tscvib_asec = sqrt(TSC.VIBRATION.SIGMA_TTX^2+TSC.VIBRATION.SIGMA_TTY^2)
  if where(tag_names(TSC) eq 'VIBRATION') ne -1 then tscvib_nm = tscvib_asec*asec2rad*TSC.DEXTMAX/2*1e9

  ;------------------------------------------
  ;FINAL MONOCHROMATIC OTF & PSF (ATM+TSC+AO)
  ;------------------------------------------
  if not keyword_set(ONLY_PSD) then begin
    tmp = dphi
    tmp = alog(TSC.OTF*vibeOTF)/alog(2.0)-0.5*tmp/alog(2.0) ; LOG2(OTF_TSC*OTF_AO) and including vibration if any
    w = where(double(tmp) lt -1022)
    if w[0] ne -1 then tmp[w] = -1022 ; removing values of OTF that are too close to zero (to avoid underflow message)
    tototf = 2.d^tmp ; final OTF
    if size(filterOTF,/type) ne 0 then tototf = tototf*filterOTF ; optional OT-CCD correction
    strehl = total(abs(tototf))/total(abs(TSC.OTF0)) ; overall Strehl ratio
    cutfreq = sqrt(total(abs(tototf))*DIM.DFF^2/!dpi)*2*DIM.LAMBDA*1d-6/TSC.DEXTMAX<1 ; relative angular cutoff frequency
    if strlowcase(strtrim(MODE,2)) ne 'seli' then co2halo = DIM.LAMBDA*1d-6/2/max([ACTPITCH,WFS_PITCH])*rad2asec
    if keyword_set(PSF) or keyword_set(X_PSF) or keyword_set(Y_PSF) or keyword_set(EE_ANALYSIS_DISC) or $
       keyword_set(EE_ANALYSIS_SQUA) or keyword_set(EE_ANALYSIS_SLIT) then begin
      tmp = dcomplexarr(DIM.N_PSF_USR,DIM.N_PSF_USR)
      tmp[DIM.N_PSF_USR/2-DIM.N_OTF/2:DIM.N_PSF_USR/2+DIM.N_OTF/2-1,DIM.N_PSF_USR/2-DIM.N_OTF/2:DIM.N_PSF_USR/2+DIM.N_OTF/2-1] = tototf
      if keyword_set(X_PSF) then begin
        totpsfx = strehl*NORMAL(abs(MATHFT(reform(total(tmp,2))*DIM.DFF,DX = DIM.DXF_USR,IC = DIM.N_PSF_USR/2,/INVERSE)))
        totpsfx = congrid(totpsfx,DIM.N_PSF_USR,cubic = -0.4,/center)
      endif
      if keyword_set(Y_PSF) then begin
        totpsfy = strehl*NORMAL(abs(MATHFT(reform(total(tmp,1))*DIM.DFF,DX = DIM.DXF_USR,IC = DIM.N_PSF_USR/2,/INVERSE)))
        totpsfy = congrid(totpsfy,DIM.N_PSF_USR,cubic = -0.4,/center)
      endif
      if not keyword_set(X_PSF) and not keyword_set(Y_PSF) then begin
        totpsf = strehl*NORMAL(abs(MATHFT(tmp,DX = DIM.DXF_USR,IC = DIM.N_PSF_USR/2,JC = DIM.N_PSF_USR/2,/INVERSE)))
      endif
    endif
  endif

  ;========================================================================== 
  ;=============== END OF PSD, OTF & PSF CALCULATION SECTION ================== 
  ;========================================================================== 

  ;========================================================================== 
  ;================== PSF & OTF POST-PROCESSING ANALYSIS ====================== 
  ;========================================================================== 

  ;-----------------------
  ; FWHM ANALYSIS [OPTION]
  ;-----------------------
  if keyword_set(FWHM_ANALYSIS) then begin
    if keyword_set(FRFFT) then begin ; if fractional FFT is implemented, we can use the code below
      ; we are using a high resolution small FoV PSF matrix just for this
      tmp = dblarr(2*DIM.N_OTF,2*DIM.N_OTF)
      tmp[DIM.N_OTF-DIM.N_OTF/2:DIM.N_OTF+DIM.N_OTF/2-1,DIM.N_OTF-DIM.N_OTF/2:DIM.N_OTF+DIM.N_OTF/2-1] = tototf
      tmp = NORMAL(abs(FRACTIONAL_FFT(tmp,2*DIM.N_OTF,8*DIM.N_OTF,DIM.N_OTF,DIM.N_OTF,DIM.N_OTF,DIM.N_OTF)))*strehl
      dxf = rad2asec/(DIM.N_OTF*DIM.DFF*8)
      mask = tmp gt 0.5*strehl
      maskprojX = total(mask,1) & maskprojY = total(mask,2)
      wx = where(maskprojX gt 0) & wy = where(maskprojY gt 0)
      ; if the FoV of this small PSF is not large enough, we double the FoV
      if max(wx)-min(wx)+1 gt 2*DIM.N_OTF-6 or max(wy)-min(wy)+1 gt 2*DIM.N_OTF-6 then begin
        tmp = dblarr(4*DIM.N_OTF,4*DIM.N_OTF)
        tmp[2*DIM.N_OTF-DIM.N_OTF/2:2*DIM.N_OTF+DIM.N_OTF/2-1,2*DIM.N_OTF-DIM.N_OTF/2:2*DIM.N_OTF+DIM.N_OTF/2-1] = tototf
        tmp = NORMAL(abs(FRACTIONAL_FFT(tmp,4*DIM.N_OTF,8*DIM.N_OTF,2*DIM.N_OTF,2*DIM.N_OTF,2*DIM.N_OTF,2*DIM.N_OTF)))*strehl
        dxf = rad2asec/(DIM.N_OTF*DIM.DFF*8)
        mask = tmp gt 0.5*strehl
        maskprojX = total(mask,1) & maskprojY = total(mask,2)
        wx = where(maskprojX gt 0) & wy = where(maskprojY gt 0)
        ;if this is still not enough we use the whole FoV
        if max(wx)-min(wx)+1 gt 4*DIM.N_OTF-6 or max(wy)-min(wy)+1 gt 4*DIM.N_OTF-6 then begin
          tmp = dblarr(DIM.N_PSF_USR,DIM.N_PSF_USR)
          tmp[DIM.N_PSF_USR/2-DIM.N_OTF/2:DIM.N_PSF_USR/2+DIM.N_OTF/2-1,DIM.N_PSF_USR/2-DIM.N_OTF/2:DIM.N_PSF_USR/2+DIM.N_OTF/2-1] = tototf
          tmp = NORMAL(abs(MATHFT(tmp)))*strehl
          dxf = DIM.DXF_USR
          mask = tmp gt 0.5*strehl
          maskprojX = total(mask,1) & maskprojY = total(mask,2)
          wx = where(maskprojX gt 0) & wy = where(maskprojY gt 0)
          if max(wx)-min(wx)+1 gt DIM.N_PSF_USR-6 or max(wy)-min(wy)+1 gt DIM.N_PSF_USR-6 then begin
            print,'################################################'
            print,'WARNING WARNING WARNING WARNING WARNING WARNING'
            print,'################################################'
            print,''
            PRINT,'PSF FULL WIDTH AT HALF MAXIMUM IS TOO LARGE WITH RESPECT TO THE MATRIX SIZE: THE FWHM CANNOT BE COMPUTED, RETURNING -1.'
            print,''
            print,'################################################'
            print,'WARNING WARNING WARNING WARNING WARNING WARNING'
            print,'################################################'
            fwhmmax = -1 ; and if even this is not a sufficient FoV, then we give up,
            fwhmmin = -1 ; because the PSF FWHM is larger than the maximal possible FoV
            asymori = 0
            goto,postfwhm
          endif
        endif
      endif
      tmp = FWHM_PSF(tmp,DXF = dxf)
      fwhmmax = tmp.max
      fwhmmin = tmp.min
      asymori = tmp.ang
      if abs(fwhmmax-fwhmmin)/(fwhmmax+fwhmmin) le 0.01 then begin
        fwhmmax = (fwhmmax+fwhmmin)*0.5d
        fwhmmin = fwhmmax
        asymori = 0.d
      endif
    endif else begin ; this is the case where the fractional FFT algorithm is not available
      tmp = FWHM_PSF(totpsf,DXF = DIM.DXF_USR)
      fwhmmax = tmp.max
      fwhmmin = tmp.min
      if abs(fwhmmax-fwhmmin)/(fwhmmax+fwhmmin) le 0.01 then begin
        fwhmmax = (fwhmmax+fwhmmin)*0.5d
        fwhmmin = fwhmmax
        asymori = 0.d
      endif else asymori = tmp.ang
    endelse
  endif
  postfwhm:

  ;-------------------------------------------
  ; ENCIRCLED ENERGY METRICS ANALYSIS [OPTION]
  ;-------------------------------------------
  if keyword_set(EE_ANALYSIS_DISC) then begin
    ede = INTENE(tototf,DIM.DFF,'disc',/E50,/E80)
    eedia50 = ede.e50 ; diameter of 50% total energy
    eedia80 = ede.e80 ; diameter of 80% total energy
  endif
  if keyword_set(EE_ANALYSIS_SQUA) then begin
    ese = INTENE(tototf,DIM.DFF,'squa',/E50,/E80)
    eesqu50 = ese.e50 ; square witdh at 50% total energy
    eesqu80 = ese.e80 ; square width at 80% total energy
  endif
  if keyword_set(EE_ANALYSIS_SLIT) then begin
    efe = INTENE(tototf,DIM.DFF,'slit',/E50,/E80)
    eeslt50 = efe.e50 ; slit witdh at 50% total energy
    eeslt80 = efe.e80 ; slit width at 80% total energy
  endif

  ;------------------------------------------------------------
  ; PROPORTION OF ENERGY WITHIN A GIVEN APERTURE WIDTH [OPTION]
  ;------------------------------------------------------------
  if keyword_set(EEW_DISC) then epdisc = (INTENE(tototf,DIM.DFF,'disc',WIDTH = EEW_DISC)).ief
  if keyword_set(EEW_SQUA) then epsqua = (INTENE(tototf,DIM.DFF,'squa',WIDTH = EEW_SQUA)).ief
  if keyword_set(EEW_SLIT) then epslit = (INTENE(tototf,DIM.DFF,'slit',WIDTH = EEW_SLIT)).ief

  ;=================================== 
  ;===== END OF POST-PROCESSING SECTION
  ;=================================== 

  ;=================================== 
  ;SAVING/PRINT RESULTS OF CALCULATION
  ;=================================== 

  ;-------------------------
  ;WRITE FITS FILES [OPTION]
  ;-------------------------
  if keyword_set(FITSCODE) then begin
    ;building fits header
    mkhdr,hdr,dblarr(DIM.N_PSF_USR,DIM.N_PSF_USR)
    if strlowcase(strtrim(MODE,2)) eq 'seli' then sxaddpar,hdr,'MODE','SELI'
    if strlowcase(strtrim(MODE,2)) eq 'ngs'  then sxaddpar,hdr,'MODE','NGS'
    if strlowcase(strtrim(MODE,2)) eq 'glao' then sxaddpar,hdr,'MODE','GLAO'
    if (where(tag_names(TSC) eq 'VIBRATION'))[0] ne -1 then sxaddpar,hdr,'TEL_VIB',tscvib_asec
    if size(VIBRATION,/type) ne 0 then sxaddpar,hdr,'INST_VIB',instvib_asec
    sxaddpar,hdr,'INST',TSC.INST
    if keyword_set(PSF) then sxaddpar,hdr,'DX_FOCUS',DIM.DXF_USR[0]
    if keyword_set(OTF) then sxaddpar,hdr,'DF_FOCUS',DIM.DFF[0]
    if keyword_set(SF) then sxaddpar,hdr, 'DX_PUPIL',DIM.DXP[0]
    if keyword_set(PSD) or keyword_set(ONLY_PSD) then sxaddpar,hdr,'DFP_HIGH',DIM.DFP[0]
    if keyword_set(PSD) or keyword_set(ONLY_PSD) then sxaddpar,hdr,'DFP_LOW',DIM.DFP_LF[0]
    sxaddpar,hdr,'LAMBDA',DIM.LAMBDA
    sxaddpar,hdr,'TSC_DIAM',TSC.DEXTMAX
    sxaddpar,hdr,'TSC_SURF',TSC.SURF
    sxaddpar,hdr,'SEEING  ',double(W0)
    sxaddpar,hdr,'r0_500nm',r0500
    if keyword_set(DISPERSION) then sxaddpar,hdr,'DISP',DISPERSION[0]
    if L0 eq -1 then sxaddpar,hdr,'L0','INF'
    if L0 ne -1 then sxaddpar,hdr,'L0',double(L0)
    if n_params() gt 6 then sxaddpar,hdr,'ALTITUDE',meanalti
    if n_params() gt 6 then sxaddpar,hdr,'ISOPLANA',anisoang
    if n_params() gt 6 and not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then begin
      sxaddpar,hdr,'VELOCITY',meanwind
      sxaddpar,hdr,'TIME_SCL',timescale
      sxaddpar,hdr,'TIME_COE',tau0
    endif
    if strlowcase(strtrim(MODE,2)) ne 'seli' then begin
      sxaddpar,hdr,'NOMINAL',dm_nactnom
      sxaddpar,hdr,'DM_NACT',dm_nact
      sxaddpar,hdr,'TOT_NACT',TSC.SURF/ACTPITCH^2
      sxaddpar,hdr,'ACTPITCH',ACTPITCH
      sxaddpar,hdr,'NLS_DIA',TSC.DEXTMAX/WFS_PITCH
      sxaddpar,hdr,'NLS_PUP',TSC.SURF/WFS_PITCH^2
      sxaddpar,hdr,'WFS_PITC',WFS_PITCH
      sxaddpar,hdr,'DM_HEIGH',DM_HEIGHT
      sxaddpar,hdr,'GS_ANG',GS_ANG
      sxaddpar,hdr,'GS_ORI',GS_ORI
      if not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then begin
        sxaddpar,hdr,'WFS_INT',gsint
        sxaddpar,hdr,'LAG',LAG
        if strlowcase(strtrim(MODE,2)) eq 'ngs' then begin
          sxaddpar,hdr,'LOOPFREQ',loopfreq
          if strlowcase(strtrim(LOOP_MODE,2)) eq 'closed' then sxaddpar,hdr,'LOOPGAIN',loopgain
          if strlowcase(strtrim(LOOP_MODE,2)) eq 'closed' then sxaddpar,hdr,'STABLE  ',stability
          if n_params() eq 21 then begin
            sxaddpar,hdr,'NGS_MAG',NGS_MAG
            sxaddpar,hdr,'NGS_TEM',NGS_TEM
            sxaddpar,hdr,'WFS_RON',WFS_RON
            sxaddpar,hdr,'SKY_RADIANCE',SKY_RADIANCE
            sxaddpar,hdr,'WFS_TAU',wfstau
            sxaddpar,hdr,'NPHOTON',wfsnph
            sxaddpar,hdr,'WFS_LAM',wfslam
            sxaddpar,hdr,'WFS_NEA',nea
          endif
          if n_params() eq 19 then sxaddpar,hdr,'WFS_NEA',WFS_NEA
        endif
        if strlowcase(strtrim(MODE,2)) eq 'glao' then begin
          if n_params() eq 22 then begin
            sxaddpar,hdr,'NGS_MAG',mean(NGS_MAG)
            sxaddpar,hdr,'NGS_TEM',mean(NGS_TEM)
            sxaddpar,hdr,'WFS_RON',WFS_RON
            sxaddpar,hdr,'SKY_RAD',SKY_RADIANCE
            sxaddpar,hdr,'WFS_TAU',mean(wfstau)
            sxaddpar,hdr,'NPHOTON',mean(wfsnph)
            sxaddpar,hdr,'WFS_LAM',mean(wfslam)
            sxaddpar,hdr,'WFS_NEA',mean(nea)
          endif
          if n_params() eq 20 then sxaddpar,hdr,'WFS_NEA',mean(WFS_NEA)
        endif
      endif
    endif
    if strlowcase(strtrim(MODE,2)) eq 'glao' then sxaddpar,hdr,'WFSTYP',GLAO_WFS.TYPE
    ;writing arrays
    if keyword_set(WAVE) then fits_write,'PSD_WAVE_'+FITSCODE+'.fits',psd_wave,hdr
    if keyword_set(WAVE) then fits_write,'PUP_WAVE_'+FITSCODE+'.fits',pup_wave,hdr
    if keyword_set(PSF)  then fits_write,'PSF_'+FITSCODE+'.fits',totpsf,hdr
    if keyword_set(OTF)  then fits_write,'OTF_'+FITSCODE+'.fits',double(tototf),hdr
    if keyword_set(SF)   then if strlowcase(strtrim(MODE,2)) ne 'seli' then if varLF gt 0 then fits_write,'LFSF_'+FITSCODE+'.fits',dphiLF,hdr
    if keyword_set(SF)   then fits_write,'HFSF_'+FITSCODE+'.fits',dphiHF,hdr
    if keyword_set(SF)   then if keyword_set(POST_TIPTILT) then fits_write,'GTILTSF_'+FITSCODE+'.fits',sfGtilt,hdr
    if keyword_set(SF)   then fits_write,'SF_'+FITSCODE+'.fits',dphi,hdr
    if keyword_set(PSD) or keyword_set(ONLY_PSD) then begin
      if strlowcase(strtrim(MODE,2)) eq 'seli' then fits_write,'PSD_ATM_'+FITSCODE+'.fits',Watm,hdr
      if strlowcase(strtrim(MODE,2)) ne 'seli' then fits_write,'PSD_FE_'+FITSCODE+'.fits',Whf,hdr
      if size(Wamp,/type) ne 0        then fits_write,'PSD_AMP_'+FITSCODE+'.fits',Wamp,hdr
      if size(Walias,/type) ne 0      then fits_write,'PSD_AL_'+FITSCODE+'.fits',Walias,hdr
      if size(Wanisoservo,/type) ne 0 then fits_write,'PSD_AS_'+FITSCODE+'.fits',Wanisoservo,hdr
      if size(Wnoise,/type) ne 0      then fits_write,'PSD_NS_'+FITSCODE+'.fits',Wnoise,hdr
      if size(Wdisp,/type) ne 0       then fits_write,'PSD_DISP_'+FITSCODE+'.fits',Wdisp,hdr
    endif
  endif

  ;------------------------------------------------------------
  ;PRINT THE RESULTS ON THE SCREEN OR INTO A LOG FILE [OPTIONS]
  ;------------------------------------------------------------
  if keyword_set(INFO) or keyword_set(LOGCODE) then begin
    if keyword_set(INFO) then unit = -1
    if keyword_set(LOGCODE) then openw,unit,'paola_'+LOGCODE+'.log',/get_lun
    aff:

    printf,unit,'ATMOSPHERIC PARAMETERS'
    printf,unit,'----------------------'
    printf,unit,format = '("   SEEING AT ZENITH @ 500 NM [asec]",d13.6)',W0
    if L0 eq -1 then printf,unit,         '                 OUTER SCALE            INFINITE'
    if L0 ne -1 then printf,unit,format = '("                 OUTER SCALE L0 [m]",d13.6)',double(L0)
    printf,unit,format = '("                ZENITH ANGLE  [deg]",d13.6)',ZA
    printf,unit,format = '("APPROX. SEEING(LAMBDA,L0,ZA) [asec]",d13.6)',0.98*DIM.LAMBDA*1e-6/r0LAM*rad2asec*aos
    printf,unit,format = '("            r0 @ 500 nm @ ZA    [m]",d13.6)',r0500
    printf,unit,format = '("            r0 @ LAMBDA @ ZA    [m]",d13.6)',r0LAM
    if n_params() ge 9 then printf,unit,format = '("           <LAYERS ALTITUDE>    [m]",d13.6)',meanalti
    if n_params() ge 9 then printf,unit,format = '("  ISOPLANATIC ANGLE @ LAMBDA [asec]",d13.6)',anisoang
    if n_params() ge 9 and not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then begin
      printf,unit,format = '("           <LAYERS VELOCITY>  [m/s]",d13.6)',meanwind
      if meanwind gt 0 then printf,unit,format = '("   PHASE TIME SCALE @ LAMBDA   [ms]",d13.6)',timescale
      if meanwind gt 0 then printf,unit,format = '("   PHASE TIME SCALE @ 500 nm   [ms]",d13.6)',tau0
      if meanwind eq 0 then          printf,unit,'   PHASE TIME SCALE @ LAMBDA            INFINITE'
      if meanwind eq 0 then          printf,unit,'   PHASE TIME SCALE @ 500 nm            INFINITE'
    endif
    if keyword_set(SCINTILLATION) then printf,unit,format = '("            SCINTILLATION INDEX [-]",e13.6)',sci_index
    printf,unit,''

    printf,unit,'TELESCOPE & FOCAL PLANE PARAMETERS'
    printf,unit,'----------------------------------'
    printf,unit,format = '("            TELESCOPE STREHL    [1]",d13.6)',TSC.STREHL
    printf,unit,format = '("             MIRROR DIAMETER    [m]",d13.6)',TSC.DEXTMAX
    printf,unit,format = '("             MIRROR  SURFACE  [m^2]",d13.6)',TSC.SURF
    printf,unit,format = '("       ANGULAR FIELD OF VIEW [asec]",d13.6)',DIM.FOV_PSF
    printf,unit,format = '("             PSF MATRIX SIZE   [px]",i13)',DIM.N_PSF_USR
    printf,unit,format = '("             OTF MATRIX SIZE   [px]",i13)',DIM.N_OTF
    printf,unit,format = '("           PUPIL MATRIX SIZE   [px]",i13)',DIM.N_OTF
    printf,unit,format = '("          IMAGING WAVELENGTH [mu-m]",d13.6)',DIM.LAMBDA
    printf,unit,format = '("         ANGULAR PIXEL SCALE  [mas]",d13.6)',DIM.DXF_USR*1d3
    printf,unit,format = '("        PSF THEORETICAL FWHM  [mas]",d13.6)',DIM.LAMBDA*1e-6/TSC.DEXTMAX*rad2asec*1e3
    if where(tag_names(TSC) eq 'VIBRATION') ne -1 then printf,unit,format = '("          TELESCOPE VIBRATION [mas]",d13.6)',tscvib_asec*1e3
    if where(tag_names(TSC) eq 'VIBRATION') ne -1 then printf,unit,format = '("          TELESCOPE VIBRATION  [nm]",d13.6)',tscvib_nm
    if where(tag_names(TSC) eq 'PHASE') ne -1 then printf,unit,format = '("  TELESCOPE SURFACE ERROR      [nm]",d13.6)',sqrt(total(TSC.PHASE^2)/total(TSC.PUPIL gt 0.5))*DIM.LAMBDA*1e3/2/!dpi*0.5
    if where(tag_names(TSC) eq 'PHASE') ne -1 then printf,unit,format = '("  TELESCOPE STATIC ABERRATION  [nm]",d13.6)',sqrt(total(TSC.PHASE^2)/total(TSC.PUPIL gt 0.5))*DIM.LAMBDA*1e3/2/!dpi
    printf,unit,''

    printf,unit,'PSF PIXEL VALUES'
    printf,unit,'-----------------------'
    if strtrim(TSC.INST,2) eq 'IMAGER'   then  printf,unit,'                 AVERAGE INTENSITY WITHIN PIXELS'
    if strtrim(TSC.INST,2) eq 'INTENSITY' then printf,unit,'            SAMPLED INTENSITY IN THE IMAGE PLANE'
    printf,unit,''

    if size(VIBRATION,/type) ne 0 then begin
      printf,unit,'INSTRUMENT VIBRATIONS'
      printf,unit,'---------------------'
      printf,unit,format = '("                    VIBRATION [mas]",d13.6)',instvib_asec*1e3
      printf,unit,''
    endif
      
    printf,unit,'AO CORRECTION MODE'
    printf,unit,'------------------'
    if strlowcase(strtrim(MODE,2)) eq 'seli' then printf,unit,'                           SEEING LIMITED, NO AO'
    if strlowcase(strtrim(MODE,2)) eq 'ngs' then printf,unit,'                       SINGLE NATURAL GUIDE STAR'
    if strlowcase(strtrim(MODE,2)) eq 'glao' then begin
      printf,unit,'       GROUND LAYER ADAPTIVE OPTICS'
      if strlowcase(GLAO_WFS.TYPE) eq 'full' then printf,unit,'             FULL WFS FOV AVERAGING'
      if strlowcase(GLAO_WFS.TYPE) eq 'edge' then printf,unit,'          EDGE OF WFS FOV AVERAGING'
      if strlowcase(GLAO_WFS.TYPE) eq 'star' then printf,unit,' GUIDE STAR CONSTELLATION AVERAGING'
    endif
    if keyword_set(POST_TIPTILT) and strlowcase(strtrim(MODE,2)) ne 'seli' then printf,unit,'                  WITH POST AO G-TILT CORRECTION'
    if keyword_set(POST_TIPTILT) and strlowcase(strtrim(MODE,2)) eq 'seli' then printf,unit,'                          WITH G-TILT CORRECTION'
    if keyword_set(TILT_ANGLE_STEP) then printf,unit,format = '("   G-TILT STEP CORRECTION ANGLE [-]",e13.6)',TILT_ANGLE_STEP
    printf,unit,''

    if strlowcase(strtrim(MODE,2)) eq 'ngs' then begin
      printf,unit,'NATURAL GUIDE STAR PARAMETERS'
      printf,unit,'-----------------------------'
      if n_params() eq 19 then if WFS_NEA eq 0 then printf,unit,'                  BRIGHT STAR MODE, NO WFS NOISE'
      if n_params() eq 19 then if WFS_NEA gt 0 then printf,unit,'      NGS MAGNITUDE NOT GIVEN, SEE WFS NEA INPUT'
      if n_params() eq 21 then printf,unit,format = '("                      MAGNITUDE [1]",d13.6)',NGS_MAG
      if n_params() eq 21 then printf,unit,         '                         FILTER           '+strupcase(FILTER)+'-BAND'
      if n_params() eq 21 then printf,unit,format = '("           BLK-BODY TEMPERATURE [K]",d13.6)',NGS_TEM
      printf,unit,''
    endif

    if strlowcase(strtrim(MODE,2)) eq 'glao' then if strlowcase(GLAO_WFS.TYPE) eq 'star' and n_params() eq 22 then begin
      printf,unit,'NATURAL GUIDE STARS PARAMETERS'
      printf,unit,'------------------------------'
      printf,unit,format = '("                 NGS MAGNITUDES [1]",'+strtrim(string(n_elements(NGS_MAG)),2)+'d13.6)',NGS_MAG
      printf,unit,format = '("FILTER ASSOCIATED TO MAGNITUDES    ",'+strtrim(string(n_elements(FILTER)),2)+'a13)',strupcase(strtrim(FILTER,2))
      printf,unit,format = '("          BLK-BODY TEMPERATURES [K]",'+strtrim(string(n_elements(NGS_TEM)),2)+'d13.6)',NGS_TEM
      if max(wfsnph) le 99999.0 then printf,unit,format = '("   NUMBER PHOTONS/FRAME/LENSLET [1]",'+strtrim(string(n_elements(wfsnph)),2)+'d13.6)',wfsnph
      if max(wfsnph) gt 99999.0 then printf,unit,format = '("   NUMBER PHOTONS/FRAME/LENSLET [1]",'+strtrim(string(n_elements(wfsnph)),2)+'e13.6)',wfsnph
      printf,unit,''
    endif

    if strlowcase(strtrim(MODE,2)) ne 'seli' then begin
      printf,unit,'SCIENCE TARGET POSITION / OPTICAL AXIS'
      printf,unit,'--------------------------------------'
      printf,unit,format = '("              OFF-AXIS ANGLE [asec]",d13.6)',GS_ANG
      printf,unit,format = '("              OFF-AXIS ANGLE [amin]",d13.6)',GS_ANG/60.d
      if GS_ANG gt 0 then printf,unit,format = '(" ORIENTATION/X-AXIS (->WEST)  [deg]",d13.6)',GS_ORI
      printf,unit,''
    endif

    if strlowcase(strtrim(MODE,2)) ne 'seli' then begin
      printf,unit,'LOOP PARAMETERS'
      printf,unit,'---------------'
      if not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then begin
        if strlowcase(strtrim(MODE,2)) ne 'glao' then printf,unit,format = '("                LOOP FREQUENCY [Hz]",d13.6)',loopfreq
        if strlowcase(strtrim(MODE,2)) ne 'glao' then printf,unit,format = '("               SAMPLING PERIOD [ms]",d13.6)',1e3/loopfreq
        printf,unit,format = '("              INTEGRATION TIME [ms]",d13.6)',gsint
        printf,unit,format = '("            TECHNICAL TIME LAG [ms]",d13.6)',LAG
        if strlowcase(strtrim(MODE,2)) ne 'glao' then printf,unit,format = '("           TOTAL LOOP TIME LAG [ms]",d13.6)',LAG+max([gsint,1e3/loopfreq])
        if strlowcase(strtrim(MODE,2)) eq 'glao' then printf,unit,format = '("           TOTAL LOOP TIME LAG [ms]",d13.6)',LAG+gsint
        printf,unit,         '                          LOOP MODE     '+strupcase(LOOP_MODE)
        if strlowcase(strtrim(LOOP_MODE,2)) eq 'closed' then begin
          printf,unit,format = '("                     LOOP GAIN  [1]",d13.6)',loopgain
          printf,unit,         '               IS THE LOOP STABLE ?     '+strupcase(stability)
        endif
      endif
      printf,unit,''
    endif

    if strlowcase(strtrim(MODE,2)) ne 'seli' then begin
      printf,unit,'WAVEFRONT SENSOR PARAMETERS'
      printf,unit,'---------------------------'
      printf,unit,                                                   '              WAVEFRONT SENSOR TYPE               SH'
      if keyword_set(ANTI_ALIAS) then printf,unit,                   '   WFS ANTI-ALIASING SPATIAL FILTER               ON'
      if not keyword_set(ANTI_ALIAS) then printf,unit,               '   WFS ANTI-ALIASING SPATIAL FILTER              OFF'
      printf,unit,                                          format = '("   LENSLET PITCH IN M1 PLANE        [m]",d13.6)',WFS_PITCH
      printf,unit,                                          format = '(" # OF LENSLETS / M1 DIAMETER        [1]",d13.6)',TSC.DEXTMAX/WFS_PITCH
      printf,unit,                                          format = '("  # OF LENSLETS / M1 SURFACE        [1]",d13.6)',TSC.SURF/WFS_PITCH^2
      if keyword_set(DISPERSION) then printf,unit,          format = '(" WFS TO DM DISPERSION FACTOR        [1]",d13.6)',DISPERSION[0]
      if not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then begin
        if size(WFS_RON,/type) ne 0 then printf,unit,       format = '("          DETECTOR READNOISE     [e/px]",d13.6)',WFS_RON
        if size(SKY_RADIANCE,/type) ne 0 then printf,unit,  format = '("     SKY BACKGROUND RADIANCE [W/m^2/sr]",d13.6)',SKY_RADIANCE
        if size(SKY_RADIANCE,/type) ne 0 then printf,unit,  format = '("     SKY BACKGROUND             [ph/px]",d13.6)',!dpi*SKY_RADIANCE*(wfslam*1e-6)^3*gsint*1e-3/16/6.62607004d-34/299792458.d
        if size(wfslam,/type) ne 0 then printf,unit,        format = '("          CENTRAL WAVELENGTH       [mu]",'+strtrim(string(n_elements(wfslam)),2)+'d13.6)',wfslam
        if size(wfstau,/type) ne 0 then printf,unit,        format = '("           OPTICAL BANDWIDTH       [mu]",'+strtrim(string(n_elements(wfsbdw)),2)+'d13.6)',wfsbdw
        if size(wfstau,/type) ne 0 then printf,unit,                               format = '("AVERAGE OPTICAL TRANSMISSION        [1]",'+strtrim(string(n_elements(wfstau)),2)+'d13.6)',wfstau
        if size(wfsnph,/type) ne 0 then if max(wfsnph) le 99999.0 then printf,unit,format = '("NUMBER PHOTONS/FRAME/LENSLET        [1]",'+strtrim(string(n_elements(wfsnph)),2)+'d13.6)',wfsnph
        if size(wfsnph,/type) ne 0 then if max(wfsnph) gt 99999.0 then printf,unit,format = '("NUMBER PHOTONS/FRAME/LENSLET        [1]",'+strtrim(string(n_elements(wfsnph)),2)+'e13.6)',wfsnph
        if size(nea,/type) ne 0 then printf,unit,                                  format = '("      NOISE EQUIVALENT ANGLE     [asec]",'+strtrim(string(n_elements(nea)),2)+'d13.6)',nea
        if size(WFS_NEA,/type) ne 0 then printf,unit,                              format = '("      NOISE EQUIVALENT ANGLE     [asec]",d13.6)',WFS_NEA
      endif
      if strlowcase(strtrim(MODE,2)) eq 'glao' then begin
        if strlowcase(GLAO_WFS.TYPE) ne 'star' then printf,unit,format = '("   GLAO WAVEFRONT SENSOR FoV [asec]",d13.6)',GLAO_WFS.ANG
        if strlowcase(GLAO_WFS.TYPE) ne 'star' then printf,unit,format = '("   GLAO WAVEFRONT SENSOR FoV [amin]",d13.6)',GLAO_WFS.ANG/60.d
        if strlowcase(GLAO_WFS.TYPE) eq 'star' then printf,unit,format = '("   GLAO WAVEFRONT SENSOR FoV [asec]",d13.6)',2*max(sqrt(GLAO_WFS.ANG[0,*]^2+GLAO_WFS.ANG[1,*]^2))
        if strlowcase(GLAO_WFS.TYPE) eq 'star' then printf,unit,format = '("   GLAO WAVEFRONT SENSOR FoV [amin]",d13.6)',2*max(sqrt(GLAO_WFS.ANG[0,*]^2+GLAO_WFS.ANG[1,*]^2))/60.d
        if strlowcase(GLAO_WFS.TYPE) eq 'star' then printf,unit,format = '("  NUMBER OF GLAO GUIDE STARS    [1]",i13)',n_elements(GLAO_WFS.ANG[0,*])
        if size(WFS_NEA,/type) ne 0 then printf,unit,format = '("  WFS NOISE EQUIVALENT ANGLE [asec]",d13.6)',WFS_NEA
      endif
      printf,unit,''

      printf,unit,'DEFORMABLE MIRROR PARAMETERS'
      printf,unit,'----------------------------'
      if (size(DMTF))[0] eq 0 then printf,unit,'          DM TRANSFER FUNCTION TYPE      PERFECT'
      if (size(DMTF))[0] eq 2 then printf,unit,'          DM TRANSFER FUNCTION TYPE    REALISTIC'
      printf,unit,format = '("           CONJUGATION ALTITUDE [m]",d13.6)',DM_HEIGHT
      printf,unit,format = '("     ACTUATOR PITCH IN M1 PLANE [m]",d13.6)',ACTPITCH
      printf,unit,format = '("   # OF ACTUATORS / M1 DIAMETER [1]",d13.6)',DM_NACT
      printf,unit,format = '("    # OF ACTUATORS / M1 SURFACE [1]",d13.6)',TSC.SURF/ACTPITCH^2
      printf,unit,''
      printf,unit,'VARIANCE OF RESIDUAL AO PHASE COMPONENTS'
      printf,unit,'----------------------------------------'
      printf,unit,format = '("   HIGH FREQUENCY (FITTING) [rad^2]",e13.3)',var[0]
      printf,unit,format = '("                   ALIASING [rad^2]",e13.3)',var[1]
      printf,unit,format = '("  ANISO-SERVO OR DISPERSION [rad^2]",e13.3)',var[2]
      printf,unit,format = '("                  WFS NOISE [rad^2]",e13.3)',var[3]
      printf,unit,'........................................'
      printf,unit,format = '("    TOTAL RESIDUAL VARIANCE [rad^2]",e13.3)',total(var)
      printf,unit,''
      printf,unit,'RMS OF RESIDUAL AO WAVEFRONT COMPONENTS'
      printf,unit,'---------------------------------------'
      printf,unit,format = '("      HIGH FREQUENCY (FITTING) [nm]",d13.3)',sqrt(var[0])*DIM.LAMBDA*1e3/2/!dpi
      printf,unit,format = '("                      ALIASING [nm]",d13.3)',sqrt(var[1])*DIM.LAMBDA*1e3/2/!dpi
      printf,unit,format = '("     ANISO-SERVO OR DISPERSION [nm]",d13.3)',sqrt(var[2])*DIM.LAMBDA*1e3/2/!dpi
      printf,unit,format = '("                     WFS NOISE [nm]",d13.3)',sqrt(var[3])*DIM.LAMBDA*1e3/2/!dpi
      printf,unit,'.......................................'
      printf,unit,format = '("TOTAL RESIDUAL WAVEFRONT ERROR [nm]",d13.3)',sqrt(total(var))*DIM.LAMBDA*1e3/2/!dpi

      printf,unit,''
      if where(tag_names(TSC) eq 'VIBRATION') ne -1 or size(VIBRATION,/type) ne 0 then begin
        printf,unit,'VIBRATIONS'
        printf,unit,'------------------------------------'
        if where(tag_names(TSC) eq 'VIBRATION') ne -1 then printf,unit,format = '("           TELESCOPE VIBRATION [nm]",d13.6)',tscvib_nm
        if size(VIBRATION,/type) ne 0 then     printf,unit,format = '("          INSTRUMENT VIBRATION [nm]",d13.6)',instvib_nm
        if where(tag_names(TSC) eq 'VIBRATION') ne -1 and size(VIBRATION,/type) eq 0 then $
           printf,unit,format = '("     TOTAL RESIDUAL AO+VIB WFE [nm]",d13.3)',sqrt((sqrt(total(var))*DIM.LAMBDA*1e3/2/!dpi)^2+tscvib_nm^2)
        if where(tag_names(TSC) eq 'VIBRATION') eq -1 and size(VIBRATION,/type) ne 0 then $
           printf,unit,format = '("     TOTAL RESIDUAL AO+VIB WFE [nm]",d13.3)',sqrt((sqrt(total(var))*DIM.LAMBDA*1e3/2/!dpi)^2+instvib_nm^2)
        if where(tag_names(TSC) eq 'VIBRATION') ne -1 and size(VIBRATION,/type) ne 0 then $
           printf,unit,format = '("     TOTAL RESIDUAL AO+VIB WFE [nm]",d13.3)',sqrt((sqrt(total(var))*DIM.LAMBDA*1e3/2/!dpi)^2+tscvib_nm^2+instvib_nm^2)
        printf,unit,''
      endif
    endif

    if not keyword_set(ONLY_PSD) then begin
      printf,unit,'FINAL PSF CHARACTERISTICS'
      printf,unit,'-------------------------'
      printf,unit,format = '("                STREHL RATIO    [1]",d13.6)',strehl
      printf,unit,format = '("  RELATIVE CUTTING FREQUENCY    [1]",d13.6)',cutfreq
      if strlowcase(strtrim(MODE,2)) ne 'seli' then printf,unit,format = '("     TRANSITION CORE -> HALO [asec]",d13.6)',co2halo
    endif
    if keyword_set(FWHM_ANALYSIS) then begin
      printf,unit,format = '("   FWHM ELLIPSOID SMALL AXIS  [mas]",d13.6)',fwhmmin*1000
      printf,unit,format = '("   FWHM ELLIPSOID LARGE AXIS  [mas]",d13.6)',fwhmmax*1000
      printf,unit,format = '("  ELLIPSE ORIENTATION/X-AXIS  [deg]",d13.6)',asymori
    endif
    if keyword_set(EE_ANALYSIS_DISC) then begin
        printf,unit,format = '("    50% ENERGY DISC DIAMETER  [mas]",d13.6)',eedia50*1000
        printf,unit,format = '("    80% ENERGY DISC DIAMETER  [mas]",d13.6)',eedia80*1000
    endif
    if keyword_set(EE_ANALYSIS_SQUA) then begin
        printf,unit,format = '("     50% ENERGY SQUARE WIDTH  [mas]",d13.6)',eesqu50*1000
        printf,unit,format = '("     80% ENERGY SQUARE WIDTH  [mas]",d13.6)',eesqu80*1000
    endif
    if keyword_set(EE_ANALYSIS_SLIT) then begin
        printf,unit,format = '("       50% ENERGY SLIT WIDTH  [mas]",d13.6)',eeslt50*1000
        printf,unit,format = '("       80% ENERGY SLIT WIDTH  [mas]",d13.6)',eeslt80*1000
    endif
    if keyword_set(EEW_DISC) then begin
        printf,unit,format = '(" FOR AN APERTURE OF DIAMETER [asec]",d13.6)',EEW_DISC
        printf,unit,format = '("  PROPORTION OF TOTAL ENERGY IS [1]",d13.6)',epdisc
    endif
    if keyword_set(EEW_SQUA) then begin
        printf,unit,format = '("FOR A SQUARE SPAXEL OF WIDTH [asec]",d13.6)',EEW_SQUA
        printf,unit,format = '("  PROPORTION OF TOTAL ENERGY IS [1]",d13.6)',epsqua
    endif
    if keyword_set(EEW_SLIT) then begin
        printf,unit,format = '("            FOR A SLIT WIDTH [asec]",d13.6)',EEW_SLIT
        printf,unit,format = '("  PROPORTION OF TOTAL ENERGY IS [1]",d13.6)',epslit
    endif
    printf,unit,' ================================================ '
    printf,unit,''
    if keyword_set(INFO) and keyword_set(LOGCODE) then begin
      if unit ne -1 then begin
        close,unit
        free_lun,unit
      endif
      INFO = 0
      unit = -1
      goto,aff
    endif
    if unit ne -1 then free_lun,unit
  endif

  ;----------------------------------------------------
  ;SAVING THE RESULTS IN A STRUCTURE VARIABLE -> result
  ;----------------------------------------------------
  if not keyword_set(ONLY_PSD) then begin
    result = create_struct('lam',DIM.LAMBDA,'inst',TSC.INST,'mode',MODE,'w0ZA',W0,'r05',r0500,'r0l',r0LAM,'L0',double(L0),'strehl',strehl,'cfr',cutfreq)
    if strlowcase(strtrim(MODE,2)) ne 'seli' then result = create_struct(result,'c2h',co2halo,'var',var,'rms',sqrt(var)*DIM.LAMBDA*1e3/2/!dpi)
  endif else begin
    result = create_struct('lam',DIM.LAMBDA,'inst',TSC.INST,'mode',MODE,'w0ZA',W0,'r05',r0500,'r0l',r0LAM,'L0',double(L0))
    if strlowcase(strtrim(MODE,2)) ne 'seli' then result = create_struct(result,'var',var,'rms',sqrt(var)*DIM.LAMBDA*1e3/2/!dpi)
  endelse
  if n_params() gt 6 then result = create_struct(result,'alt',meanalti,'ani',anisoang)
  if n_params() ge 9 and not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then result = create_struct(result,'win',meanwind,'lft',timescale,'tau0',tau0)
  if strlowcase(strtrim(MODE,2)) ne 'seli' then begin
    result = create_struct(result,'nna',DM_NACTnom,'ana',DM_NACT,'tna',TSC.SURF/ACTPITCH^2,'act',ACTPITCH,$
                         'anl',TSC.DEXTMAX/WFS_PITCH,'tnl',TSC.SURF/WFS_PITCH^2,'wlp',WFS_PITCH,'dm_height',DM_HEIGHT,'ang',GS_ANG,'ori',GS_ORI)
    if size(WFS_INT,/type) ne 0 then result = create_struct(result,'int',gsint,'lag',LAG)
    if strlowcase(strtrim(MODE,2)) eq 'ngs' then result = create_struct(result,'lfr',loopfreq,'lpg',loopgain)
    if strlowcase(strtrim(MODE,2)) eq 'ngs' then if strlowcase(strtrim(LOOP_MODE,2)) eq 'closed' then result = create_struct(result,'stability',stability)
    if size(NGS_MAG,/type) ne 0 then result = create_struct(result,'mag',NGS_MAG,'tem',NGS_TEM,'ron',WFS_RON,'sbg',!dpi*SKY_RADIANCE*(wfslam*1e-6)^3*gsint*1e-3/16/6.62607004d-34/299792458.d,$
                                                            'wfstau',wfstau,'nph',wfsnph,'nea',nea,'wfslam',wfslam,'wfsbdw',wfsbdw)
    if size(WFS_NEA,/type) ne 0 then result = create_struct(result,'nea',WFS_NEA)
    if strlowcase(strtrim(MODE,2)) eq 'glao' then result = create_struct(result,'glao_wfs',GLAO_WFS)
    if size(GS_WEIGHT,/type) ne 0 then result = create_struct(result,'gs_weight',GS_WEIGHT)
  endif
  if keyword_set(FWHM_ANALYSIS) then result = create_struct(result,'fmi',fwhmmin,'fma',fwhmmax,'dir',asymori)
  if keyword_set(EE_ANALYSIS_DISC) then result = create_struct(result,'e50disc',eedia50,'e80disc',eedia80)
  if keyword_set(EE_ANALYSIS_SQUA) then result = create_struct(result,'e50squa',eesqu50,'e80squa',eesqu80)
  if keyword_set(EE_ANALYSIS_SLIT) then result = create_struct(result,'e50slit',eeslt50,'e80slit',eeslt80)
  if keyword_set(EEW_DISC) then result = create_struct(result,'eew_disc',double(EEW_DISC),'epdisc',epdisc)
  if keyword_set(EEW_SQUA) then result = create_struct(result,'eew_squa',double(EEW_SQUA),'epsqua',epsqua)
  if keyword_set(EEW_SLIT) then result = create_struct(result,'eew_slit',double(EEW_SLIT),'epslit',epslit)
  if keyword_set(SCINTILLATION) then result = create_struct(result,'sci_index',sci_index)
  if keyword_set(DISPERSION) then result = create_struct(result,'dispersion',DISPERSION[0])
  if keyword_set(TILT_ANGLE_STEP) then result = create_struct(result,'tilt_angle_step',TILT_ANGLE_STEP)

  ; ADDITION OF THE OPTIONAL OUTPUTS /PSF, /OTF, /SF, AND /PSD
  if keyword_set(PSF) then result = create_struct(result,'psf',totpsf,'dxf',DIM.DXF_USR)
  if keyword_set(OTF) then result = create_struct(result,'otf',tototf,'dff',DIM.DFF)
  if keyword_set(SF) then begin
    result = create_struct(result,'sf',dphi,'dxp',DIM.DXP)
    if strlowcase(strtrim(MODE,2)) ne 'seli' then result = create_struct(result,'lfsf',dphiLF,'hfsf',dphiHF)
    if keyword_set(POST_TIPTILT) then result = create_struct(result,'sfGtilt',sfGtilt)
    if keyword_set(SCINTILLATION) then result = create_struct(result,'amsf',dphiamp)
  endif
  if keyword_set(PSD) or keyword_set(ONLY_PSD) then begin
    if strlowcase(strtrim(MODE,2)) eq 'seli' then result = create_struct(result,'PSD_ATM',Watm,'dfp_hf',DIM.DFP[0])
    if strlowcase(strtrim(MODE,2)) ne 'seli' then result = create_struct(result,'PSD_FE',Whf,'dfp_hf',DIM.DFP[0],'dfp_lf',DIM.DFP_LF[0])
    if size(Walias,/type) ne 0 then result = create_struct(result,'PSD_AL',Walias)
    if size(Wanisoservo,/type) ne 0 then result = create_struct(result,'PSD_AS',Wanisoservo)
    if size(Wnoise,/type) ne 0 then result = create_struct(result,'PSD_NS',Wnoise)
    if keyword_set(SCINTILLATION) then result = create_struct(result,'PSD_AMP',Wamp)
    if size(Wdisp,/type) ne 0 then result = create_struct(result,'PSD_DISP',Wdisp)
  endif
  if keyword_set(WAVE) then begin
    result = create_struct(result,'PSD_WAVE',psd_wave,'DFP_WAVE',dfp_wave,'DXP_WAVE',dxp_wave,'DIM_WAVE',dim_wave)
    if (where(tag_names(TSC) eq 'PUPIL'))[0] ne -1 then result = create_struct(result,'PUP_WAVE',pup_wave)
    if (where(tag_names(TSC) eq 'PHASE'))[0] ne -1 then result = create_struct(result,'STA_WAVE',sta_wave)
  endif
  if size(VIBRATION,/type) ne 0 then result = create_struct(result,'INST_VIBRATION',VIBRATION)
  if (where(tag_names(TSC) eq 'VIBRATION'))[0] ne -1 then result = create_struct(result,'TSC_VIBRATION',TSC.VIBRATION)
  if keyword_set(LOGCODE) and keyword_set(INFO) then print,'***********MESSAGE : A LOG FILE HAS BEEN CREATED : paola_'+LOGCODE+'.log'

  ;----------------------
  ;LOOP STABILITY WARNING
  ;----------------------
  if strlowcase(strtrim(MODE,2)) eq 'ngs' then if strlowcase(strtrim(LOOP_MODE,2)) eq 'closed' then if stability eq 'no' then begin
    print,'################################################'
    print,'WARNING WARNING WARNING WARNING WARNING WARNING'
    print,'################################################'
    print,''
    PRINT,'THE CURRENT LOOP FREQUENCY, LOOP LAG AND GAIN MAKE'
    PRINT,'THE LOOP INSTABLE ACCORDING TO THE NYQUIST CRITERION.'
    PRINT,''
    PRINT,'********* THIS AO SYSTEM WILL NOT BE STABLE *********'
    PRINT,''
    PRINT,'SOLUTION: DECREASE THE LOOP GAIN OR/AND DECREASE THE'
    PRINT,'LOOP FREQUENCY OR LOOP TIME LAG'
    print,''
    print,'################################################'
    print,'WARNING WARNING WARNING WARNING WARNING WARNING'
    print,'################################################'
    print,''
  endif
  
  ;----------------------
  ;LOW ORDER PSD SAMPLING WARNING
  ;----------------------
  if size(warning,/type) ne 0 then if warning eq 'yes' then begin
    print,'################################################'
    print,'WARNING WARNING WARNING WARNING WARNING WARNING'
    print,'################################################'
    print,''
    print,'YOU MUST BE AWARE THAT BECAUSE EITHER THE WFS'
    print,'INTEGRATION TIME, OR THE LOOP TIME LAG, OR THE'
    print,'ANGULAR SEPARATION BETWEEN THE SCIENCE OBJECT'
    print,'AND THE GUIDE STARS WAS TOO LARGE, AN'
    print,'UNDERSAMPLING OF THE SERVO-LAG AND/OR THE'
    print,'ANISOPLANATISM POWER SPECTRUM HAS OCCURED.'
    print,'PLEASE CAREFULLY CHECK THE RESULT OF YOUR'
    print,'CALCULATION: VARIANCES ARE CERTAINLY'
    print,'UNDERESTIMATED.'
    print,''
    print,'################################################'
    print,'WARNING WARNING WARNING WARNING WARNING WARNING'
    print,'################################################'
    print,''
  endif
  
  ;----------------------
  ;WFS EXPOSURE TIME LIMIT WARNING
  ;----------------------
  if strlowcase(strtrim(MODE,2)) ne 'seli' and not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then if gsint/timescale gt 1.01 then begin
    print,'################################################'
    print,'WARNING WARNING WARNING WARNING WARNING WARNING'
    print,'################################################'
    print,''
    print,'THE WAVEFRONT SENSOR EXPOSURE TIME THAT YOU HAVE'
    print,'SET IS LONGER THAN THE EXPECTED PHASE TIME SCALE' 
    print,'AT THE SCIENCE SENSOR WAVELENGHT WHERE WE CAN'
    print,'CONSIDER THAT THE TAYLOR ASSUMPTION STILL HOLDS,'
    print,'I.E. BEFORE BOILING STARTS TO DOMINATE.'
    print,'TO SECURE THE PREDICTION, IT IS THEREFORE HIGHLY'
    print,'RECOMMENDED TO LIMIT THE WAVEFRONT SENSOR EXPOSURE'
    print,'TIME TO '+trim(timescale)+' ms'
    print,''
    print,'################################################'
    print,'WARNING WARNING WARNING WARNING WARNING WARNING'
    print,'################################################'
    print,''
  endif

  return,result

end
