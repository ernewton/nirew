;============================================
; find the pseudocontinuum

FUNCTION EW_PSEUDO, $
  lambda, flux, $
  continuum, $
  verbose=verbose, $
  flat=flat, $
  mean=mean
  
  IF ~KEYWORD_SET(verbose) THEN verbose = 0

  c1 = WHERE(lambda GE continuum[0,0] AND lambda LE continuum[1,0], $
    count1)
  c2 = WHERE(lambda GE continuum[0,1] AND lambda LE continuum[1,1], $
    count2)

  IF count1 LE 3 OR count2 LE 3 THEN BEGIN
    pseudo=-1 
    IF verbose GT 0 THEN print, "EW_PSEUDO: <= 3 elements in array. Check your lambda array and continuum interval."
  ENDIF ELSE BEGIN
  
    ; remove zero-order correction only (mean or median)
    IF KEYWORD_SET(flat) THEN BEGIN
      IF KEYWORD_SET(verbose) THEN print, "EW_PSEUDO: removing zero-order term only."
      IF KEYWORD_SET(mean) THEN $ 
	offset = mean(flux[[c1[*],c2[*]]]) $
      ELSE $
	offset = median(flux[[c1[*],c2[*]]])
      pseudo = FLTARR(N_ELEMENTS(lambda)) + offset
      
    ; remove a first-order fit (line)
    ENDIF ELSE BEGIN
      IF KEYWORD_SET(mean) THEN BEGIN
	x1 = mean(lambda[c1])
	y1 = mean(flux[c1])
	x2 = mean(lambda[c2])
	y2 = mean(flux[c2])
      ENDIF ELSE BEGIN
	x1 = median(lambda[c1])
	y1 = median(flux[c1])
	x2 = median(lambda[c2])
	y2 = median(flux[c2])
      ENDELSE
      pseudo = (y2-y1)/(x2-x1)*(lambda-x1)+y1
    ENDELSE
    
  ENDELSE

  RETURN, pseudo
END
