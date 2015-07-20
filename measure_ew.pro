; MEASURE_EW
; Measure the EW given the feature and continuum windows. Return EW in Angstroms.


FUNCTION MEASURE_EW, $
  lambda, flux, $
  continuum, $
  feature, $
  showplot=showplot, $
  rvoffset=rvoffset, $ ; apply an RV offset
  nan=nan, $ ; bad measurements are NaNs
  zero=zero ; bad measurements are zeros


  IF KEYWORD_SET(zero) THEN BEGIN
    slambda = lambda[WHERE(flux GT 0)]
    sflux = flux[WHERE(flux GT 0)]  
  ENDIF ELSE IF KEYWORD_SET(nan) THEN BEGIN
    slambda = lambda[WHERE(FINITE(flux))]
    sflux = flux[WHERE(FINITE(flux))]
  ENDIF ELSE BEGIN
    slambda = lambda ; don't modify the original value
    sflux = flux
  ENDELSE

  ; do rv offset
  IF KEYWORD_SET(rvoffset) THEN $
    slambda = slambda - rvoffset/(3.e5)*slambda

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
  
    mkct
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