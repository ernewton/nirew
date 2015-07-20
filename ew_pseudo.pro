; find the pseudocontinuum
FUNCTION EW_PSEUDO, $
  lambda, flux, $
  continuum

  c1 = WHERE(lambda GE continuum[0,0] AND lambda LE continuum[1,0], $
    count1)
  c2 = WHERE(lambda GE continuum[0,1] AND lambda LE continuum[1,1], $
    count2)

  IF count1 LE 3 OR count2 LE 3 THEN pseudo=-1 $
  ELSE BEGIN
    x1 = median(lambda[c1])
    y1 = median(flux[c1])
    x2 = median(lambda[c2])
    y2 = median(flux[c2])
    pseudo = (y2-y1)/(x2-x1)*(lambda-x1)+y1
  ENDELSE

  RETURN, pseudo
END
