;============================================
; Estimate spectral type from water index

FUNCTION SP_HIND, hind

  RETURN, 25.4 - 24.2*hind

END