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
;      MEASURE_HBAND, data, [MG1=, MG2=, MG3=,
;		ALA=, ALB=, CCORR=, /CONTF,
;		TEFF=, RAD=, LUM=, 
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
;	ccorr = cross-correlation routine to use ['c_correlate']
;	contf = flag to use contf continuum function
;	atrest = flag to skip shift to rest velocity
;
; OUTPUTS:
;	mg1 = EW of Mg at 1.50 microns, in Angstroms
;	mg2 = EW of Mg at 1.57 microns, in Angstroms
;	mg3 = EW of Mg at 1.71 microns, in Angstroms
;	ala = EW of Al at 1.67 microns (blue line), in Angstroms
;	alb = EW of Mg at 1.67 microns (red line), in Angstroms
;	teff = Estimated effective temperature, in K
;	rad = Estimated radius, in Rsun
;	lum = Estimated log luminosity, in Lsun
;
; EXAMPLE:
; 	data = MRDFITS('spec/J0200+1303_tc.fits')
; 	d0 = MRDFITS('spec/J0200+1303.fits')
; 	data[*,0,*] = d0[*,0,*]
; 	measure_hband, data, teff=teff, rad=rad, lum=lum, ew=ew
;	print, teff
;	print, rad
;	print, logl
;
; METHOD:
;	Shifts system to rest velocity by cross-correlating with a standard.
;	Measures the EW and calculates stellar properties using the calibration
;	from Newton et al. (2015).
;
; PROCEDURES USED:
;	nirew
;	tellrv
; 
; NOTE: EW defintions as originally published in Newton et al. (2015) are from
; 	an earlier version of the code. The correct definitions, which reproduce
;	the EWs used in that paper, are used here.
;-

PRO MEASURE_HBAND, data, $
  ew=ew, eew=eew, $
  teff=teff, rad=rad, lum=lum, $
  ccorr=ccorr, contf=contf, atrest=atrest, $
  showplot=showplot, quiet=quiet
  
  ; standard RV file
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
    IF KEYWORD_SET(atrest) THEN rv0 = 0 ELSE $
      ERN_RV, mydata, std[*,*,sorder], wrange=wrange, rv0=rv0, ccorr=ccorr, contf=contf, quiet=quiet

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
      ewi = measure_ew(lambda0,flux0,continuum,feature, showplot=showplot, quiet=quiet)
      IF KEYWORD_SET(showplot) THEN wait, 1
      ewveci[j] = ewi
      j++
    ENDFOREACH
    
    IF i EQ 0 THEN BEGIN ; measured value
      ew = ewveci
    ENDIF ELSE BEGIN ; to calculate errors
      ewvec[i-1] = ewveci
    ENDELSE


  ENDFOR

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

END