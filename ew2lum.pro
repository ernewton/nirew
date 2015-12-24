FUNCTION EW2LUM, mg1, mg3, emg1, emg3, err=err, erand=erand

  res = [0.832, -0.176, 0.266, -3.491]
  loglum =  res[0] * mg3 $
    + res[1] * mg3^2 $
    + res[2] * mg1 $
    + res[3]
    
  sigl = 0.049
  sysl = 0.02

  IF KEYWORD_SET(emg1) AND KEYWORD_SET(emg3) THEN BEGIN
    erand = SQRT( (res[0]*emg3 + 2.*res[1]*mg3*emg3)^2 + (res[2]*emg1)^2)
    err = SQRT( erand^2 + sigl^2 )
  ENDIF
  
  RETURN, loglum
END