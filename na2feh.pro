;+
; NAME:
;	NA2FEH
; PURPOSE:
;	Return an M dwarfs' metallicity from its Na EW.
; EXPLANATION:
;	Apply Newton et al. (2014) relationship to estimate metallicity. 
;
; CALLING SEQUENCE:
;      feh = NA2FEH(na, [e_na_low, e_na_high, 
;		EFEH=, BFEH=, /BOOTSTRAP])
;
; INPUTS:
;	na = Equivalent Width of Na line at 2.2 microns in Angstroms
;	
; OPTIONAL INPUTS:
;	e_na_low = negative error on Na EW, error from EW calculated
;		only if provided
;	e_na_high = positive error on Na EW [e_na_low]
;
; OPTIONAL KEYWORD INPUTS:
;	bootstrap = estimate error using the best-fit bootstraps
;
; OUTPUTS:
;	feh = Metallicity estimated from Na EW
;	efeh = Error on metallicity
;	bfeh = All points from Monte Carlo simulation
;		See measure_na for example of use
;
; EXAMPLE:
;	na = 5.0
;	ena = 0.1
;	feh = NA2FEH(na, ena, efeh=efeh)
;
; METHOD:
;	Calculates metallicity from Na EW from Newton et al. (2014) relation.
;	Optionally calculates error due to Na EW error, and error due to
;	uncertainty in the fit (with /BOOTSTRAP). Combines the two in quadrature.
;
; PROCEDURES USED:
;	nirew
;	tellrv
;
;-
FUNCTION NA2FEH, na, e_na_low, e_na_high, efeh=feh_err, bfeh=feh_boot, bootstrap=bootstrap

  ; calculate metallicity
  res00=0.596111
  res01=-0.0392013
  res02=-1.96390
  feh=res00*na+res01*na^2+res02

  ; uncertainty from Na EW uncertainty
  IF KEYWORD_SET(e_na_low)  THEN BEGIN
    IF ~KEYWORD_SET(e_na_high) THEN BEGIN
      e_na_high = e_na_low
      if KEYWORD_SET(verbose) THEN print, "NA_FEH: Assuming symmetric errors."
    ENDIF
    fehna_low=res00*(na-e_na_low)+res01*(na-e_na_low)^2+res02
    fehna_high=res00*(na+e_na_high)+res01*(na+e_na_high)^2+res02
    fehna_err=(fehna_high-fehna_low)/2.
    feh_err = fehna_err
  ENDIF

  ; uncertainty from calibration
  IF KEYWORD_SET(bootstrap) THEN BEGIN
    READCOL,'na_boot.res', res0, res1, const, /silent
    boot_low=FLTARR(N_ELEMENTS(na))
    boot_high=FLTARR(N_ELEMENTS(na))
    feh_boot=res0#na + res1#(na^2) + CMREPLICATE(const,N_ELEMENTS(na))    
    FOR i=0,N_ELEMENTS(na)-1 DO BEGIN
      conf=confidence_interval(feh_boot[*,i])
      boot_low[i]=conf[0]
      boot_high[i]=conf[2]
    ENDFOR
    boot_err=(boot_high-boot_low)/2.
    feh_err = boot_err
  END

  ; combine the errors
  IF KEYWORD_SET(fehna_err) AND KEYWORD_SET(boot_err) THEN $
    feh_err=SQRT(fehna_err^2+boot_err^2)

  RETURN, feh

END