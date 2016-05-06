FUNCTION water_index, lambda, flux, sp=sp, verbose=verbose

    IF ~KEYWORD_SET(verbose) THEN verbose = 0

    IF KEYWORD_SET(sp) THEN BEGIN
	IF (size(sp))[0] EQ 2 AND (size(sp))[1] EQ 2 AND (size(sp))[2] EQ 3 THEN BEGIN
	    IF verbose GT 0 THEN print, "WATER_INDEX: Using provided index"
	ENDIF ELSE BEGIN
	    IF verbose GT 0 THEN print, "WATER_INDEX: Provided value failed. Defaulting to H20-K2. Input array: ", sp
	    sp = [[2.070,2.090],[2.235,2.255],[2.360,2.380]]
	ENDELSE
    ENDIF ELSE BEGIN
	    IF verbose GT 0 THEN print, "WATER_INDEX: Defaulting to H20-K2"
	    sp = [[2.070,2.090],[2.235,2.255],[2.360,2.380]]
    ENDELSE 

    out = fltarr(3)
    FOR i=0,2 DO BEGIN
	roi = WHERE(lambda GE sp[0,i] AND lambda LE sp[1,i], count)
	IF count GT 0 THEN $
	  out[i] = MEDIAN(flux[roi]) $
	ELSE $
	  RETURN, !values.f_nan
    ENDFOR
    
    index = out[0] * out[2] / (out[1]^2.)

    RETURN, index
    
END
