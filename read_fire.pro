; READ_FIRE
; Read a spectrum produced by FIRE

FUNCTION read_fire, star, dir=dir, suffix=suffix, noerr=noerr, irtf=irtf, struct=struct, verbose=verbose, hdr=hd

  IF KEYWORD_SET(verbose) THEN silent = 0 ELSE silent = 1

  IF ~KEYWORD_SET(suffix) THEN $
    suffix = '_F.fits'
  IF ~KEYWORD_SET(dir) THEN $
    dir = '' $
  ELSE $
    dir = dir+'/'
  f_file = dir+star+suffix

  ; read in the flux file, wavelengths in microns
  flx = readfits(f_file,hd, silent=silent)
  wave = 10.^(findgen(n_elements(flx))*sxpar(hd,'cdelt1')+sxpar(hd,'crval1'))/10000

  ; lower resolution of FIRE to IRTF
  IF KEYWORD_SET(irtf) THEN BEGIN
    resol = 4.5 ; approx
    flx[WHERE(flx GT 0)] = gauss_smooth(flx[WHERE(flx GT 0)], resol)
  ENDIF

  ; no error file
  IF KEYWORD_SET(noerr) THEN BEGIN
    IF KEYWORD_SET(struct) THEN $
      data = {wave:wave, flux:flx} $
    ELSE $
      data = [[wave],[flx]]
    
  ; read in the error file
  ENDIF ELSE BEGIN
    e_file = dir+sxpar(hd,'SIGFILE')
    err = readfits(e_file,ehd, silent=silent)
    IF KEYWORD_SET(irtf) THEN $
      err[WHERE(flx GT 0)] = gauss_smooth(err[WHERE(flx GT 0)], resol)

    IF KEYWORD_SET(struct) THEN $
      data = {wave:wave, flux:flx, err:err} $
    ELSE $
      data = [[wave],[flx],[err]]
      
  ENDELSE
    
  RETURN, data

END
