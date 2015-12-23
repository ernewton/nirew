;============================================
; return 1-sigma confidence interval from a sample

FUNCTION confidence_interval, x
  n = N_ELEMENTS(x)
  left = 0.157*n
  middle = 0.5*n
  right = (1.0 - 0.157)*n
  y = x[SORT(x)]
  RETURN, y[[left, middle, right]]
END
