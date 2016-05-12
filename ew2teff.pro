FUNCTION EW2TEFF, mg1, ala, alb, emg1, eala, ealb, err=err, erand=erand

;   res = [271.4, 392.7, 2427]
  res = [261.9, 413.8, 2402]
  teff = res[0] * ala $
    + res[1] * mg1/alb $
    + res[2]
    
  sigt = 73
  syst = 0.

  IF KEYWORD_SET(emg1) AND KEYWORD_SET(eala) AND KEYWORD_SET(ealb) THEN BEGIN
    erand = SQRT( (res[0]*eala)^2 + (res[1]*mg1/alb^2*ealb)^2 + (res[1]*emg1/alb)^2 )
    err = SQRT( erand^2 + sigt^2 )
  ENDIF
  
RETURN, teff
END

; difference: sigma = 9 K, median = 4 K