;+
; NAME:
;	MEASURE_EW
; PURPOSE:
;	Calculate equivalenth widths (EWs).
; EXPLANATION:
;	Measure the EW given the feature and continuum windows. Return EW in Angstroms.
;
; CALLING SEQUENCE:
;      ew = MEASURE_EW(lambda, flux, continuum, feature,
;		[RVOFFSET=, /NAN, /ZERO, 
;		/SHOWPLOT, /QUIET])
;
; INPUTS:
;	lambda = Wavelength array
;	flux = Flux array
;	continuum = Array with continuum window defintion
;	feature = Array with feature window definition
;	
; OPTIONAL INPUTS:
;	None
;
; OPTIONAL KEYWORD INPUTS:
;	rvoffset = An RV shift to apply [None]
;	nan = Flag to remove NaNs from the data
;	zero = Flag to remove zeros from the data
;	showplot = Flag to show plot
;	quiet = Flag to limit messages
;
; OUTPUTS:
;	ew = measure EW in Angstroms
;
; EXAMPLE:
;	data = MRDFITS('spec/J0455+0440W_tc.fits')
;  	continuum = [[blue1, blue2]],[red1, red2]]
;  	feature = [f1,f2]
;	ew = measure_ew(data[*,0,0], data[*,1,0], $
;		continuum, feature, /showplot)
;
; METHOD:
;	Fits a straight line to the continuum regions to use as a
;	pseudo-continuum. Uses the trapezoidal rule to calculate the
;	equivalenth width. IT IS RECOMMENDED THAT YOU OVERSAMPLE THE
;	SPECTRUM TO MITIGATE EDGE EFFECTS. See examples.
;
; PROCEDURES USED:
;	nirew
;
;-


FUNCTION MEASURE_EW, $
  lambda, flux, $
  continuum, $
  feature, $
  rvoffset=rvoffset, $ ; apply an RV offset
  nan=nan, $ ; bad measurements are NaNs
  zero=zero, $ ; bad measurements are zeros
  showplot=showplot, $
  quiet=quiet

  slambda = lambda ; don't modify the original value
  sflux = flux
  IF KEYWORD_SET(zero) THEN BEGIN
    slambda = slambda[WHERE(sflux GT 0)]
    sflux = sflux[WHERE(sflux GT 0)]  
    IF ~KEYWORD_SET(quiet) THEN print, "MEASURE_EW: using non-zero fluxes."
  IF KEYWORD_SET(nan) THEN BEGIN
    slambda = slambda[WHERE(FINITE(sflux))]
    sflux = sflux[WHERE(FINITE(sflux))]
    IF ~KEYWORD_SET(quiet) THEN print, "MEASURE_EW: using finite fluxes."
  IF ~KEYWORD_SET(nan) AND ~KEYWORD_SET(zero) THEN IF ~KEYWORD_SET(quiet) THEN print, "MEASURE_EW: using all provided data."


  ; do rv offset
  IF KEYWORD_SET(rvoffset) THEN BEGIN
    slambda = slambda - rvoffset/(3.e5)*slambda
    print, "MEASURE_EW: applying providing RV offset."
  ENDIF

  ; find pseudocontinuum 
  pseudo=ew_pseudo(slambda, sflux, continuum)

  ; feature region of interest
  roi=WHERE(slambda GE feature[0] AND slambda LE feature[1], count)
  IF pseudo[0] EQ -1 OR count LE 3 THEN $
    out=!values.f_nan $
  ELSE BEGIN	
    ; do trapezoidal sum
    width=TSUM(slambda[roi],1-sflux[roi]/pseudo[roi])

    ; output
    out=width*10^4.
  ENDELSE

  IF KEYWORD_SET(showplot) AND FINITE(out) THEN BEGIN
  
    plot, lambda, sflux, xrange=[continuum[0,0]-.012,continuum[1,1]+.012], thick=2, /ynozero, /xsty
    oplot,[feature[0],feature[0]],[0,100],color=2
    oplot,[feature[1],feature[1]],[0,100],color=2
    oplot,[continuum[0,0],continuum[0,0]],[0,10],color=3,linestyle=2
    oplot, [continuum[0,1],continuum[0,1]],[0,100],color=3,linestyle=2
    oplot, [continuum[1,1],continuum[1,1]],[0,100],color=3,linestyle=2
    oplot, [continuum[1,0],continuum[1,0]],[0,100],color=3,linestyle=2
    oplot, lambda, pseudo, color=3

    print, 'MEASURE_EW: Measured EW in Angstroms: ', out
  ENDIF
  
  RETURN, out
END