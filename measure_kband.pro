;+
; NAME:
;	MEASURE_KBAND
; PURPOSE:
;	Calculate an M dwarfs' metallicity
; EXPLANATION:
;	Calculates Na EW and uses Newton et al. (2014) relationship to
;	estimate metallicity. Also calculates H2O-K2 index (Rojas-Ayala et al.
;	2012), and NIR spectral (Rojas-Ayala et al. 2012, Newton et al. 2014)
;	type.
;
; CALLING SEQUENCE:
;      MEASURE_KBAND, data, [FEH=, EFEH=,
;		NA=, ENA=, 
;		HINDK2=, SPTYPE=,
;		NITERS=, /ERROR,
;		CCORR=, /CONTF, /ATREST,
;		/SHOWPLOT, /QUIET]
;
; INPUTS:
;	data = Array containing [[wave],[flux],[eflux]]
;		eflux is required if ERROR is set
;	
; OPTIONAL INPUTS:
;	None
;
; OPTIONAL KEYWORD INPUTS:
;	niters = number of iterations for Monte Carlo error simulation [100]
;	error = flag to do Monte Carlo error simulation [0]
;	ccorr = cross-correlation routine to use ['c_correlate']
;	contf = flag to use contf continuum function
;	atrest = flag to skip shift to rest velocity
;	showplot = flag show plots for debugging
;	quiet = flag to suppress messages
;
; OUTPUTS:
;	feh = Metallicity estimate from Na EW in dex
;	efeh = Error on metallicity
;	na = Na EW in Angstroms
;	ena = Error on Na EW
;	hindk2 = H2O-K2 index
;	sptype = Spectral type estimated from hindk2
;	
;
; EXAMPLE:
;	data = MRDFITS('spec/J0200+1303_tc.fits')
;	measure_kband, data, na=naew, ena=enaew, feh=nafeh, efeh=enafeh, /error
;	print, naew, '+-', enaew
;	print, nafeh, '+-', enafeh
;
; METHOD:
;	Shifts system to rest velocity by cross-correlating with a standard.
;	Measures the Na EW and calculates metallicity using the calibration
;	from Newton et al. (2014). Optionally perform a Monte Carlo analysis
;	to estimate errors in Na EW and in metallicity.
;
; PROCEDURES USED:
;	nirew
;	tellrv
;
;-

PRO measure_kband, data, $
  feh=nafeh, efeh=efeh, $
  na=ew, ena=eew, $
  hindk2=hind, sptype=sp_hind, $
  niters=ni, error=doerrors, $
  ccorr=ccorr, contf=contf, atrest=atrest, $
  showplot=showplot, quiet=quiet, $
  std=std

  IF KEYWORD_SET(doerrors) THEN BEGIN
    IF ~KEYWORD_SET(ni) THEN ni=100
    ewvec = FLTARR(ni)
;     hvec = FLTARR(ni)
  ENDIF ELSE ni = 0
  
  ; standard RV file
  IF ~KEYWORD_SET(std) THEN $
    std = MRDFITS('$NIREW/spec/J0727+0513_rest.fits',0, /silent)

  ; line definitions
  READCOL, '$NIREW/linedefs.txt', lineall, f1all,f2all, c1all, c2all, c3all, c4all,  format='A,F,F,F,F,F,F'
  k = 7 ; Na at 2.2 microns
  continuum = [[c1all[k], c2all[k]],[c3all[k],c4all[k]]]
  feature = [f1all[k],f2all[k]]

  ; settings
  sorder = 0
  IF SIZE(data, /N_DIMEN) EQ 2 THEN order = 0 ELSE order = sorder
  wrange = [2.18, 2.3]

  FOR i=0, ni DO BEGIN
  
    ; add errors
    mydata = data[*,*,order]
    IF KEYWORD_SET(doerrors) THEN IF i GT 0 THEN BEGIN
      rand=RANDOMN(seed,n_elements(data[*,0,order]))
      mydata[*,1] = mydata[*,1]+rand*mydata[*,2]
    ENDIF
    
    ; shift to rest
    IF KEYWORD_SET(quiet) THEN iquiet = quiet ELSE IF i GT 0 THEN iquiet=1
    IF KEYWORD_SET(atrest) THEN rv0 = 0 ELSE $
      ERN_RV, mydata, std[*,*,sorder], wrange=wrange, rv0=rv0, ccorr=ccorr, contf=contf, quiet=iquiet

    ; oversample flux
    inc0 = N_ELEMENTS(mydata[*,0])*10.
    lambda0 = REBIN(mydata[*,0]*(1. - rv0/(3.*10.^5)),inc0)
    flux0 = REBIN(mydata[*,1],inc0)

    ; calculate parameters of interest
    ewi = measure_ew(lambda0,flux0,continuum,feature, quiet=iquiet, showplot=showplot)
;     hi = water_index(lambda0, flux0)
    
    IF i EQ 0 THEN BEGIN ; measured value
      ew = ewi
      nafeh = NA2FEH(ewi)
      hind = water_index(lambda0, flux0)
      sp_hind = HIND2SP(hind)
    ENDIF ELSE BEGIN ; to calcualte errors
      ewvec[i-1] = ewi
;       hvec[i-1] = hi
    ENDELSE


  ENDFOR
 
  IF KEYWORD_SET(doerrors) THEN BEGIN
    
    ; error in Na EW
    conf_na = confidence_interval(ewvec)
    eew = (conf_na[2]-conf_na[0])/2.
    
    ; this takes the sample of EW measurements and runs each through the bootstrapped sample of fits
    feh = NA2FEH(ewvec, bfeh=bfeh, bootstrap=doerrors)
    conf_feh = confidence_interval(bfeh)
    rand_efeh = (conf_feh[2]-conf_feh[0])/2. ; random error
    intr_efeh = 0.12 ; intrinsic scatter in relation
    
    efeh = SQRT(rand_efeh^2 + intr_efeh^2)
    
;     ; the following combines the error from the EWs with the error from the bootstrapping as if they are independent, equivalent to NA2FEH(ew, e_ew, e_feh=e_feh, /bootstrap)
;     temp = confidence_interval(fehvec) ; EW Na error
;     print, sqrt((median(e_feh))^2 + ((temp[2]-temp[1])/2.)^2) ; e_feh is bootstrap error

  ENDIF
 
END
