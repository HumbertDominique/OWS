;+
;FUNCTION MATHFT(OBJECT[,DX=,IC=,JC=,/INVERSE,/COO,/SINGLE])
;
;  MATHFT computes the discrete Fourier transform (DFT) (or 
;  inverse DFT) of a 1D or 2D object. The difference with the 
;  IDL FFT function is that the amplitude of the result is such 
;  that the DFT gives exactly the theoretical value of the 
;  continuous Fourier transform at each frequency discrete value.
;  In other words, the OBJECT energy is conserved, and Parseval's
;  theorem is directly applicable:
;  1 dimension : total(abs(object)^2)*DX=total(abs(MATHFT(OBJECT,...))^2)*df
;  2 dimensions: total(abs(object)^2)*DX^2=total(abs(MATHFT(OBJECT,...))^2)*df^2
;
;  Variance and Fourier transform:
;  1 dimension : total(abs(MATHFT(OBJECT,...))^2)*df = variance(OBJECT) * T
;  where T is the length of the signal, for instance for a time signal
;  that would be the duration of the signal.
;  2 dimensions : total(abs(MATHFT(OBJECT,...))^2)*df^2 = variance(OBJECT) * S
;  where S is the surface of the support of the signal.
;
;  MATHFT works on any matrix size, even not a power of 2. But at
;  least, matrix size must be even. If the input object matrix as 
;  an odd size, then the top line and the right column will be 
;  truncated.
;
;  Note that output is DOUBLE precision by default, unless SINGLE
;  is set.
;
;INPUT name | type | unit | default
;
;  OBJECT | array 1D or 2D | - | -
;  array from which we compute the Fourier or inverse Fourier
;  transform. Size of the array must be even (not necessarily a
;  power of 2). OBJECT type must be a numerical IDL type
;  i.e. from 1 to 6, or 9.
;
;OPTIONS
;
;  DX | SCALAR | LENGTH | 1
;  spatial sampling in the object space
;
;  IC,JC | SCALAR REAL | PX | 0,0
;  pixel position associated with the (0,0) coordinate in the
;  linear space (the conjugate of the frequency space). JC is not
;  needed if OBJECT is 1D.
;  *******IMPORTANT*******
;  IC and JC can be real, in a sense that the 0,0 DOES NOT HAVE to be
;  located on a pixel but can be anywhere in between pixels.
;  *******IMPORTANT******* 
;
;  /INVERSE for computing the inverse FT of OBJECT.
;
;  /COO if set, the output is added a coordinate array, spatial 
;  frequency or direct (x,y) space.
;
;  /SINGLE if set, calculations are done in single precision, and 
;  the result is given in single precision too.
;
;NOTE
;
;  I have checked if the amplitude of the result is correct
;  following the mathematical definition of the Fourier
;  transform, and this is OK. Beside, I have checked if the sign
;  of the result is OK too. Concerning the energy
;  conservation, I have checked that
;  total(abs(object)^2)*DX^2=total(abs(objet FT)^2)*df^2.
;  All the check have been done in 1D and 2D. Last point : the
;  aliasing. This is usual concern with Discrete Fourier
;  Transform. The user has to be aware of that problem, and
;  should use MATHFT knowing that the high frequencies are going
;  to be folded in the (-Fn,Fn) domain. Nothing can be done
;  against this, except choosing a sufficiently small sampling of
;  the object.
;
;OUTPUT name | type | unit
;
;  res | single/double complex array | -
;  Inpt object Fourier transform
;
;OPTIONAL OUTPUT name | type | unit
; 
;  If /COO is set, the output becomes a structure variable:
;
;  *.y is the Fourier transform, and
;
;  *.coo | double complex array | unit(DX) or 1/unit(DX)
;  x abscissa in the FT space (frequency if direct FT, spatial
;  position if reverse FT)
;
;EXTERNAL PROCEDURES/FUNCTION NEEDED
;
;  INTX.PRO
;  COOGRID.PRO
;
;HISTORY
;
;  Jul  8, 2002 Laurent Jolissaint HIA/NRC, written.
;  Jul 17, 2002 LJ modified to accept non power of 2 matrix sizes.
;  Sep 16, 2002 LJ added checking of the object type.
;  Oct 01, 2002 LJ improved error management.
;  Oct 24, 2002 LJ improved documentation.
;  Feb 06, 2003 LJ changed UNITGRID & UNITGRIDFFT for COOGRID.
;  Feb 07, 2003 LJ changed the output structure variable x to coo.
;  May 16, 2003 LJ introduced the ERR_EXIT procedure for error handling.
;  May 16, 2003 LJ modified to accept column vectors.
;  Jun 11, 2003 LJ corrected documentation.
;  Nov 27, 2003 LJ removed output .coo.y, .coo.r and .coo.t
;                  and optimized memory space use.
;  Aug 07, 2005 LJ corrected a bug when OBJECT size is odd.
;  Jul 10, 2006 LJ set /COO as an option.
;  Jul 11, 2006 LJ optimization of memory space use.
;  Aug 30, 2006 LJ Some minor changes (internal syntax).
;  Sep 14, 2006 LJ added test for undefined input array.
;  Mar 20, 2008 LJ introduced /SINGLE, improved doc.
;  Mar 27, 2008 LJ code reorganization to minimize memory use.
;  Mar 27, 2008 LJ removed ECLAT function.
;  Sep 14, 2010 LJ changed OBJECT to reform(OBJECT).
;  Sep 14, 2010 LJ changed ERR_EXIT procedure for MESSAGE.
;  Nov 14, 2012 LJ improved user manual.
;  Mar 07, 2014 LJ ()->[] for arrays
;  Mar 11, 2015 LJ improved user manual.
;
;BUGS report to laurent.jolissaint@heig-vd.ch
;
;-
FUNCTION MATHFT,OBJECT,DX=DX,IC=IC,JC=JC,INVERSE=INVERSE,COO=COO,SINGLE=SINGLE

  ; argument check
  if n_params() ne 1 then message,'WRONG NUMBER OF ARGUMENTS. NEEDS ONLY ONE.'
  if size(OBJECT,/type) eq 0 then message,'INPUT ARRAY UNDEFINED'
  sz=size(OBJECT)
  if sz[0] eq 0 then message,'OBJECT MUST BE 1D OR 2D'
  if n_elements(where(size(OBJECT,/dim) gt 1)) gt 2 then message,'OBJECT MUST BE 1D OR 2D'
  sz=size(reform(OBJECT))
  if sz[0] eq 1 and sz[1] eq 1 then message,'OBJECT MUST BE 1D OR 2D'
  if sz[0] eq 2 then if sz[1] ne sz[2] then message,'OBJECT MATRIX MUST BE SQUARE'
  type=size(reform(OBJECT),/type)
  if not (type ge 1 and type le 6 or type eq 9) then message,'WRONG OBJECT TYPE. MUST BE EITHER INT/REAL/COMPLEX'
  if size(DX,/type) ne 0 then if DX le 0 then message,'DX MUST BE > 0'
  if sz[0] eq 1 and size(JC,/type) ne 0 then message,'JC NOT NEEDED HERE: ARRAY IS 1D.'
  if size(IC,/type) eq 0 and size(JC,/type) ne 0 then message,'IF JC IS SET, IC MUST BE SET TOO.'

  ;settings
  if not keyword_set(DX) then DX=1.0
  if not keyword_set(IC) then IC=0.0
  if not keyword_set(JC) then JC=0.0
  fNyqu=0.5/float(DX)
  df=1.0/(sz[1]-(sz[1] mod 2))/float(DX)

  ;Fourier transform: case w/o shift element (IC=0 and JC=0)
  if ((sz[0] eq 1 and IC eq 0) or (sz[0] eq 2 and IC eq 0 and JC eq 0)) and not keyword_set(COO) then begin
    if not keyword_set(INVERSE) then begin
      if sz[0] eq 1 then return,shift(fft((reform(OBJECT))[0:sz[1]-1-(sz[1] mod 2)]),(sz[1]-(sz[1] mod 2))/2)/df
      if sz[0] eq 2 then return,shift(fft((reform(OBJECT))[0:sz[1]-1-(sz[1] mod 2),0:sz[1]-1-(sz[1] mod 2)]),(sz[1]-(sz[1] mod 2))/2,(sz[1]-(sz[1] mod 2))/2)/df^2
    endif else begin
      if sz[0] eq 1 then return,fft(  df*shift((reform(OBJECT))[0:sz[1]-1-(sz[1] mod 2)],(sz[1]-(sz[1] mod 2))/2),/inverse)
      if sz[0] eq 2 then return,fft(df^2*shift((reform(OBJECT))[0:sz[1]-1-(sz[1] mod 2),0:sz[1]-1-(sz[1] mod 2)],(sz[1]-(sz[1] mod 2))/2,(sz[1]-(sz[1] mod 2))/2),/inverse)
    endelse
  endif

  ;abscissa construction OPTIONAL
  if keyword_set(COO) then begin
    if sz[0] eq 1 then begin
      if not keyword_set(INVERSE) then abscissa=INTX(sz[1]-(sz[1] mod 2),-fNyqu,fNyqu-df)
      if     keyword_set(INVERSE) then abscissa=INTX(sz[1]-(sz[1] mod 2),-IC*DX,(sz[1]-(sz[1] mod 2)-1-IC)*DX)
    endif
    if sz[0] eq 2 then  begin
      if keyword_set(SINGLE) then begin
        if not keyword_set(INVERSE) then abscissa=(COOGRID(sz[1]-(sz[1] mod 2),sz[1]-(sz[1] mod 2),SCALE=fNyqu,/FT,/SINGLE,/COO_X)).x
        if     keyword_set(INVERSE) then abscissa=(COOGRID(sz[1]-(sz[1] mod 2),sz[1]-(sz[1] mod 2),XMIN=-IC*float(DX),XMAX=(sz[1]-(sz[1] mod 2)-1-IC)*float(DX),$
                                                                                                   YMIN=-JC*float(DX),YMAX=(sz[1]-(sz[1] mod 2)-1-JC)*float(DX),/SINGLE,/COO_X)).x
      endif else begin
        if not keyword_set(INVERSE) then abscissa=(COOGRID(sz[1]-(sz[1] mod 2),sz[1]-(sz[1] mod 2),SCALE=fNyqu,/FT,/COO_X)).x
        if     keyword_set(INVERSE) then abscissa=(COOGRID(sz[1]-(sz[1] mod 2),sz[1]-(sz[1] mod 2),XMIN=-IC*DX,XMAX=(sz[1]-(sz[1] mod 2)-1-IC)*DX,$
                                                                                                   YMIN=-JC*DX,YMAX=(sz[1]-(sz[1] mod 2)-1-JC)*DX,/COO_X)).x
      endelse
    endif
  endif

  ;Fourier transform: case w/o shift element (IC=0 and JC=0) and /COO set
  if ((sz[0] eq 1 and IC eq 0) or (sz[0] eq 2 and IC eq 0 and JC eq 0)) and keyword_set(COO) then begin
    if not keyword_set(INVERSE) then begin
      if sz[0] eq 1 then return,{y:shift(fft((reform(OBJECT))[0:sz[1]-1-(sz[1] mod 2)]),(sz[1]-(sz[1] mod 2))/2)/df,coo:abscissa}
      if sz[0] eq 2 then return,{y:shift(fft((reform(OBJECT))[0:sz[1]-1-(sz[1] mod 2),0:sz[1]-1-(sz[1] mod 2)]),(sz[1]-(sz[1] mod 2))/2,(sz[1]-(sz[1] mod 2))/2)/df^2,coo:abscissa}
    endif else begin
      if sz[0] eq 1 then return,{y:fft(  df*shift((reform(OBJECT))[0:sz[1]-1-(sz[1] mod 2)],(sz[1]-(sz[1] mod 2))/2),/inverse),coo:abscissa}
      if sz[0] eq 2 then return,{y:fft(df^2*shift((reform(OBJECT))[0:sz[1]-1-(sz[1] mod 2),0:sz[1]-1-(sz[1] mod 2)],(sz[1]-(sz[1] mod 2))/2,(sz[1]-(sz[1] mod 2))/2),/inverse),coo:abscissa}
    endelse
  endif

  ;shift element -> tmp=exp(i2pi*f*r)
  if sz[0] eq 1 then begin
    tmp=intx(sz[1]-(sz[1] mod 2),-fNyqu,fNyqu-df)
    xc=IC*DX
    if keyword_set(SINGLE) then begin
      tmp=2*!pi*xc*tmp
      tmp=complex(cos(tmp),sin(tmp))
    endif else begin
      tmp=2*!dpi*xc*tmp
      tmp=dcomplex(cos(tmp),sin(tmp))
    endelse
  endif
  if sz[0] eq 2 then begin
    if keyword_set(SINGLE) then begin
      tmp=(COOGRID(sz[1]-(sz[1] mod 2),sz[1]-(sz[1] mod 2),SCALE=fNyqu,/FT,/SINGLE,/COO_X)).x
      xc=IC*DX
      yc=JC*DX
      tmp=2*!pi*(xc*tmp+yc*rotate(tmp,1))
      tmp=complex(cos(tmp),sin(tmp))
    endif else begin
      tmp=(COOGRID(sz[1]-(sz[1] mod 2),sz[1]-(sz[1] mod 2),SCALE=fNyqu,/FT,/COO_X)).x
      xc=IC*DX
      yc=JC*DX
      tmp=2*!dpi*(xc*tmp+yc*rotate(tmp,1))
      tmp=dcomplex(cos(tmp),sin(tmp))
    endelse
  endif

  ;Fourier transform
  if not keyword_set(INVERSE) then begin
     if sz[0] eq 1 then tmp=tmp*shift(fft((reform(OBJECT))[0:sz[1]-1-(sz[1] mod 2)]),(sz[1]-(sz[1] mod 2))/2)/df
     if sz[0] eq 2 then tmp=tmp*shift(fft((reform(OBJECT))[0:sz[1]-1-(sz[1] mod 2),0:sz[1]-1-(sz[1] mod 2)]),(sz[1]-(sz[1] mod 2))/2,(sz[1]-(sz[1] mod 2))/2)/df^2
  endif else begin
     if sz[0] eq 1 then tmp=fft(  df*shift((reform(OBJECT))[0:sz[1]-1-(sz[1] mod 2)]/tmp,(sz[1]-(sz[1] mod 2))/2),/inverse)
     if sz[0] eq 2 then tmp=fft(df^2*shift((reform(OBJECT))[0:sz[1]-1-(sz[1] mod 2),0:sz[1]-1-(sz[1] mod 2)]/tmp,(sz[1]-(sz[1] mod 2))/2,(sz[1]-(sz[1] mod 2))/2),/inverse)
  endelse

  ;return
  if not keyword_set(COO) then return,tmp
  return,{y:tmp,coo:abscissa}

end
