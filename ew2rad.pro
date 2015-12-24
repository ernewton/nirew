FUNCTION EW2RAD, mg2, ala
rad = -0.0489 * mg2 $
  + 0.275 * ala $
  + 0.201 * mg2/ala $
  - 0.216
sigr = 0.027
sysr = 0.006
RETURN, rad
END