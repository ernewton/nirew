;============================================
; Estimate spectral type from water index

FUNCTION HIND2SP, hind

  RETURN, 25.4 - 24.2*hind

END