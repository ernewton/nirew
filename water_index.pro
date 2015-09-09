FUNCTION water_index, lambda, spectra, sp=sp, verbose=verbose

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

    FOR i=0,2 DO BEGIN
        sp[*,i] = MEDIAN(spectra[WHERE(lambda GE sp[0,i] AND lambda LE sp[1,i])])
    ENDFOR
    sp = sp[0,*]
    index = sp[0] * sp[2] / (sp[1]^2.)

    RETURN, index
    
END
