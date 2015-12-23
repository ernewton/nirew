;============================================
; Get metallicity from Na line

FUNCTION NA_FEH, na, na_low, na_high, feh_low=feh_low, feh_high=feh_high, feh_err=feh_err

res00=0.596111
res01=-0.0392013
res02=-1.96390
feh=res00*na+res01*na^2+res02

IF KEYWORD_SET(na_low) AND KEYWORD_SET(na_high) THEN BEGIN

	READCOL,'na_boot.res', res0, res1, const, /silent
	
	fehna_low=res0[0]*na_low+res1[0]*na_low^2+const[0]
	fehna_high=res0[0]*na_high+res1[0]*na_high^2+const[0]
	fehna_err=(fehna_high-fehna_low)/2.

	boot_low=FLTARR(N_ELEMENTS(na))
	boot_high=FLTARR(N_ELEMENTS(na))
	feh_boot=res0#na + res1#(na^2) + CMREPLICATE(const,N_ELEMENTS(na))
	FOR i=0,N_ELEMENTS(na)-1 DO BEGIN
		conf=confidence_interval(feh_boot[*,i])
		boot_low[i]=conf[0]
		boot_high[i]=conf[2]
	ENDFOR
	boot_err=(boot_high-boot_low)/2.

	feh_err=SQRT(fehna_err^2+boot_err^2)

END

RETURN, feh

END