FUNCTION EW2RAD, mg2, ala, emg2, eala, err=err, erand=erand

  res = [-0.0489, 0.275,  0.201, - 0.216]
  rad = res[0] * mg2 $
    + res[1] * ala $
    + res[2] * mg2/ala $
    + res[3]
    
  sigr = 0.027
  sysr = 0.006

  IF KEYWORD_SET(emg2) AND KEYWORD_SET(eala) THEN BEGIN
    erand = SQRT( (res[0]*emg2 + res[2]*emg2/ala)^2 + (res[1]*eala - res[2]*eala*mg2/ala^2)^2 )
    err = SQRT( erand^2 + sigr^2 )
  ENDIF
  
  RETURN, rad
  
END