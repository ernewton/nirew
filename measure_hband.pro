;+
; NAME:
;	MEASURE_HBAND
; PURPOSE:
;	Calculate an M dwarfs' stellar properties
; EXPLANATION:
;	Calculates H-band EWs and uses Newton et al. (2015) relationship to
;	estimate temperature, radius, and luminosity. 
;
; CALLING SEQUENCE:
;      MEASURE_HBAND, data, [TEFF=, RAD=, LUM=, 
;		ETEFF=, ERAD=, ELUM=,
;		EW=, EEW=,
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
;	teff = Estimated effective temperature, in K
;	rad = Estimated radius, in Rsun
;	lum = Estimated log luminosity, in Lsun
;	eteff = Error on teff
;	erad = Error on rad
;	elum = Error on luminosity;
;	ew = Array of EWs in Angstroms, including:
;		mg1 = EW of Mg at 1.50 microns
;		mg2 = EW of Mg at 1.57 microns
;		mg3 = EW of Mg at 1.71 microns
;		ala = EW of Al at 1.67 microns (blue line)
;		alb = EW of Mg at 1.67 microns (red line)
;	eew = Array of errors on EWs
;	
; EXAMPLE:
; 	data = MRDFITS('spec/J0200+1303_tc.fits')
; 	d0 = MRDFITS('spec/J0200+1303.fits')
; 	data[*,0,*] = d0[*,0,*]
; 	measure_hband, data, teff=teff, rad=rad, lum=lum, eteff=eteff, erad=erad, elum=elum, ew=ew
;
; METHOD:
;	Shifts system to rest velocity by cross-correlating with a standard.
;	Measures the EW and calculates stellar properties using the calibration
;	from Newton et al. (2015). Optionally perform a Monte Carlo analysis
;	to estimate errors in Na EW and in metallicity.
;
; PROCEDURES USED:
;	nirew
;	tellrv
; 
; NOTE: EW defintions as originally published in Table 1 of Newton et al. (2015) are 
;	from an earlier version of the code. The correct definitions, which reproduce
;	the EWs used in that paper, are used here.
;-

PRO MEASURE_HBAND, data, $
  teff=teff, rad=rad, lum=lum, $
  eteff=eteff, erad=erad, elum=elum, $
  ew=ew, eew=eew, $
  niters=ni, error=doerrors, $
  ccorr=ccorr, contf=contf, atrest=atrest, $
  showplot=showplot, quiet=quiet, $
  std=std
  
  ; standard RV file
  IF ~KEYWORD_SET(std) THEN $
    std = MRDFITS('$NIREW/spec/J0727+0513_rest.fits',0, /silent)

  ; line definitions
  READCOL, '$NIREW/linedefs.txt', lineall, f1all,f2all, c1all, c2all, c3all, c4all,  format='A,F,F,F,F,F,F'
  lines = [1,2,4,5,6]
  ; 1 Mg(1.50) mg1
  ; 2 Mg(1.58) mg2
  ; 4 Al-a(1.67) ala
  ; 5 Al-b(1.67) alb
  ; 6 Mg(1.71) mg3

  IF KEYWORD_SET(doerrors) THEN BEGIN
    IF ~KEYWORD_SET(ni) THEN ni=100
    ewvec = FLTARR(ni, N_ELEMENTS(lines))
  ENDIF ELSE ni = 0


  ; settings
  sorder = 1
  IF SIZE(data, /N_DIMEN) EQ 2 THEN order = 0 ELSE order = sorder
  wrange = [1.49, 1.73]

  FOR i=0, ni DO BEGIN
  
    ; add errors
    mydata = data[*,*,order]
    IF KEYWORD_SET(doerrors) THEN BEGIN
      rand=RANDOMN(seed,n_elements(data[*,0,order]))
      mydata[*,1] = mydata[*,1]+rand*mydata[*,2]
    ENDIF
    
    ; shift to rest
    IF KEYWORD_SET(quiet) THEN iquiet = quiet ELSE IF i GT 0 THEN iquiet=1
    IF KEYWORD_SET(atrest) THEN rv0 = 0 ELSE $
      ERN_RV, mydata, std[*,*,sorder], wrange=wrange, rv0=rv0, ccorr=ccorr, contf=contf, quiet=iquiet

    ; oversample flux
    inc0 = N_ELEMENTS(mydata[*,0])*10.
    lambda0 = REBIN(mydata[*,0]*(1. - rv0/(3.e5)),inc0)
    flux0 = REBIN(mydata[*,1],inc0)
  
    ; calculate parameters of interest
    j = 0
    ewveci = FLTARR(N_ELEMENTS(lines))
    FOREACH k, lines DO BEGIN
      continuum = [[c1all[k], c2all[k]],[c3all[k],c4all[k]]]
      feature = [f1all[k],f2all[k]]
      ewi = measure_ew(lambda0,flux0,continuum,feature, showplot=showplot, quiet=iquiet)
      IF KEYWORD_SET(showplot) THEN wait, 1
      ewveci[j] = ewi
      j++
    ENDFOREACH
    
    IF i EQ 0 THEN BEGIN ; measured value
      ew = ewveci
    ENDIF ELSE BEGIN ; to calculate errors
      ewvec[i-1,*] = ewveci
    ENDELSE

  ENDFOR

  ; calculate stellar parameters
  mg1 = ew[0]
  mg2 = ew[1]
  ala = ew[2]
  alb = ew[3]
  mg3 = ew[4]
  IF ~KEYWORD_SET(quiet) THEN BEGIN
    print, "MG1 EW = ", mg1
    print, "MG2 EW = ", mg2
    print, "ALA EW = ", ala
    print, "ALB EW = ", alb
    print, "MG3 EW = ", mg3
  ENDIF
  
  teff = EW2TEFF(mg1, ala, alb)
  rad = EW2RAD(mg2, ala)
  lum = EW2LUM(mg1, mg3)

  ; errors on things
  IF KEYWORD_SET(doerrors) THEN BEGIN
    
    ; error in EW
    eew = FLTARR(N_ELEMENTS(ew))
    FOR i=0, N_ELEMENTS(ew)-1 DO BEGIN
      conf = confidence_interval(ewvec[*,i])
      eew[i] = (conf[2]-conf[0])/2.
    ENDFOR
    
    ; this takes the sample of EW measurements and runs each through the best fit
    arr = EW2TEFF(ewvec[*,0], ewvec[*,2], ewvec[*,3])
    conf = confidence_interval(arr)
    rand = (conf[2]-conf[0])/2. ; EW error
    intr = 73. ; intrinsic scatter in relation
    eteff = SQRT(rand^2 + intr^2)
    ; Al-b and Al-b EWs are somewhat correlated so error is smaller this way
    ;	than if supplying EW errors to program

    arr = EW2RAD(ewvec[*,1], ewvec[*,2])
    conf = confidence_interval(arr)
    rand = (conf[2]-conf[0])/2. ; EW error
    intr = 0.027 ; intrinsic scatter in relation
    erad = SQRT(rand^2 + intr^2)

    arr = EW2LUM(ewvec[*,0], ewvec[*,4])
    conf = confidence_interval(arr)
    rand = (conf[2]-conf[0])/2. ; EW error
    intr = 0.049 ; intrinsic scatter in relation
    elum = SQRT(rand^2 + intr^2)

    IF ~KEYWORD_SET(quiet) THEN BEGIN
      print, "teff =   ", teff, " +-", eteff
      print, "rad =    ", rad, " +-", erad
      print, "loglum = ", lum, " +-", elum
      print, "Random errors from EWs combined in quadrature with intrinsic scatter."
    ENDIF

  ENDIF ELSE BEGIN
  
    IF ~KEYWORD_SET(quiet) THEN BEGIN
      print, "teff =   ", teff
      print, "rad  =   ", rad
      print, "loglum = ", lum
    ENDIF

  ENDELSE

END