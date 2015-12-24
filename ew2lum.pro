FUNCTION EW2LUM, mg1, mg3
  loglum = 0.832 * mg3 $
    - 0.176 * mg3^2 $
    + 0.266 * mg1 $
    - 3.491
  sigl = 0.049
  sysl = 0.02
  
  RETURN, loglum
END