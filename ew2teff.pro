FUNCTION EW2TEFF, mg1, ala, alb

  teff = 271.4 * ala $
    + 392.7 * mg1/alb $
    + 2427
sigt = 73
syst = 0.
RETURN, teff
END