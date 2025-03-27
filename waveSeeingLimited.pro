FUNCTION WAVE,DIMMAT,DXP,R0,OSCALE,D1=D1,D2=D2,LAMBDA=LAMBDA

  ;SOME SETTINGS
  rad2asec = 3600.d*180.d/!dpi
  asec2rad = 1.d/rad2asec

    DFP = 1.d/DXP/DIMMAT
    r0l = double(R0*(LAMBDA/0.5d)^1.2)
    freq = (COOGRID(DIMMAT,DIMMAT,/ft,scale = 1.d/(2*DXP),/radius)).r ; spatial frequency pupil plane radius
    tmp1 = dblarr(DIMMAT,DIMMAT)
    w = where(freq ne 0)
    tmp1[w] = 0.0229d*r0l^(-5.d/3)*(freq[w]^2+double(OSCALE ne -1)/OSCALE^2)^(-11.d/6) ; spatial power spectrum with outer scale

  ;TELESCOPE PUPIL MASK CHECK [OPTIONAL INPUT]
  if keyword_set(TSCPUP) then begin
    radius = (COOGRID(DIMMAT,DIMMAT,/ft,scale = DXP*DIMMAT/2,/radius)).r
    puprad = max(radius[where(TSCPUP gt 0.5)]) ; get pupil maximum radius
  endif

  ;INSTANTANEOUS PHASE FOURIER TRANSFORM
  tmp2 = dblarr(DIMMAT+1,DIMMAT+1) ; puissance de phase rad^2/m^(-2)
  tmp2[0:DIMMAT-1,0:DIMMAT-1] = tmp1*(DIMMAT*DXP)^2
  tmp2[0:DIMMAT-1,DIMMAT] = tmp2[0:DIMMAT-1,0]
  tmp2[DIMMAT,0:DIMMAT-1] = tmp2[0,0:DIMMAT-1]
  tmp2[DIMMAT,DIMMAT] = tmp2[0,0]
  tmp2 = sqrt(tmp2)*randomn(S,DIMMAT+1,DIMMAT+1) ; random draw to create a random phase spectrum rad/m^(-1) AMPLITUDE
  tmp1 = sqrt(2)*(tmp2+rotate(tmp2,2))/2 ; forced even amplitude of wf spectrum (we put it in variable tmp1 to save some memory space)
  tmp2 = randomu(S,DIMMAT+1,DIMMAT+1)*2*!dpi ; random draw for the phase of the spectrum of the optical phase
  ; in order to make sure that the phase is real, 
  ; (1) the real part of the phase spectrum is forced to be even
  ; (2) the imaginary part of phase spectrum is forced to be odd
  tmp2 = dcomplex((tmp1*cos(0.5*(tmp2-rotate(tmp2,2))))[0:DIMMAT-1,0:DIMMAT-1],$
                  (tmp1*sin(0.5*(tmp2-rotate(tmp2,2))))[0:DIMMAT-1,0:DIMMAT-1]) ; phase FT
  if keyword_set(SPECTRUM) then phaseft = tmp2

  ;INSTANTANEOUS PHASE
  tmp2 = double(MATHFT(tmp2,dx = dxp,/inverse,ic = DIMMAT/2,jc = DIMMAT/2))

  ;PUPIL MASK [OPTION]
  xyrtpupil = COOGRID(DIMMAT,DIMMAT,/ft,scale = DXP*DIMMAT/2)
  if size(D1,/type) ne 0 then pupil = xyrtpupil.r le double(D1)/2 and xyrtpupil.r ge double(D2)/2
  if size(TSCPUP,/type) ne 0 then pupil = TSCPUP

  ;PHASE AVERAGE SUBTRACTION
  tmp2 = tmp2-mean(tmp2)
  if size(D1,/type) ne 0 or size(TSCPUP,/type) ne 0 then tmp2[where(pupil gt 0.5)] = tmp2[where(pupil gt 0.5)]-mean(tmp2[where(pupil gt 0.5)])
  if size(D1,/type) ne 0 or size(TSCPUP,/type) ne 0 then tmp2[where(pupil lt 0.5)] = 0

  ;INSTANTANEOUS PSF [OPTION]
  if size(pupil,/type) ne 0 then begin
    tmp2 = pupil*tmp2
    tmp2[where(pupil ne 0)] = tmp2[where(pupil ne 0)]-mean(tmp2[where(pupil ne 0)])
    Sp = total(pupil)*DXP^2
    apsf = MATHFT(pupil*exp(-dcomplex(0,1)*tmp2),dx = DXP,ic = DIMMAT/2,jc = DIMMAT/2)/Sp
    psf = abs(apsf)^2
    if n_params() eq 5 then dxf = DFP*1e-6*LAMBDA/!dpi*180*3600
  endif

  ;RETURN
  if not keyword_set(SPECTRUM) then begin
    if size(pupil,/type) eq 0 then return,{phase:tmp2,dxp:DXP}
    if size(pupil,/type) ne 0 then if n_params() eq 2 then return,{phase:tmp2,dxp:DXP,apsf:apsf,psf:psf,pupil:pupil}
    if size(pupil,/type) ne 0 then if n_params() eq 5 then return,{phase:tmp2,dxp:DXP,apsf:apsf,psf:psf,pupil:pupil,dxf:dxf}
  endif else begin
    if size(pupil,/type) eq 0 then return,{phase:tmp2,dxp:DXP,phaseft:phaseft,psd:psd,dfp:DFP}
    if size(pupil,/type) ne 0 then if n_params() eq 2 then return,{phase:tmp2,dxp:DXP,apsf:apsf,psf:psf,pupil:pupil,phaseft:phaseft,psd:psd,dfp:DFP}
    if size(pupil,/type) ne 0 then if n_params() eq 5 then return,{phase:tmp2,dxp:DXP,apsf:apsf,psf:psf,pupil:pupil,dxf:dxf,phaseft:phaseft,psd:psd,dfp:DFP}
  endelse

end
