; READ_FIRE
; Read a spectrum produced by FIRE

FUNCTION read_fire, star, dir=dir, suffix=suffix, noerr=noerr, irtf=irtf, struct=struct, verbose=verbose, hdr=hd, merged=merged

  IF N_ELEMENTS(merged) EQ 0 OR KEYWORD_SET(merged) THEN comb = 1 ELSE comb = 0
  IF KEYWORD_SET(verbose) THEN silent = 0 ELSE silent = 1

  IF ~KEYWORD_SET(suffix) THEN $
    suffix = '_F.fits'
  IF ~KEYWORD_SET(dir) THEN $
    d = '' $
  ELSE $
    d = dir+'/'
    
  ; orders read in individually
  IF comb EQ 0 THEN BEGIN
    orders = 11+INDGEN(21)
    ; not sure how big biggest array will be! pad with NaNs...
    all_data = FLTARR(3000, 3, N_ELEMENTS(orders))*!values.f_nan
    
  ; read in merged data instead
  ENDIF ELSE BEGIN
    orders = -1
    all_data = !values.f_nan
  ENDELSE
  
  ; begin loop over orders
  FOR i = 0, N_ELEMENTS(orders)-1 DO BEGIN
  
    IF TOTAL(orders) EQ -1 THEN $
      f_file = STRCOMPRESS(d+star+suffix, /REMOVE_ALL) $
    ELSE $
      f_file = STRCOMPRESS(d+star+'_order'+string(orders[i])+suffix, /REMOVE_ALL)
      
    ; does the file exist??
    IF FILE_TEST(f_file) THEN BEGIN
      print, "FILE: ", f_file
   
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
      
	IF TOTAL(ORDERS) NE -1 THEN $
	  e_file = STRCOMPRESS(d+star+'_order'+string(orders[i])+'_E.fits', /REMOVE_ALL) $
	ELSE $
	  e_file = d+fxpar(hd,'SIGFILE')
	IF ~FILE_TEST(e_file) THEN $
	    e_file = STRCOMPRESS(d+star+'_E.fits', /REMOVE_ALL) 
	    
	print, "SIGFILE: ", e_file
	err = readfits(e_file,ehd, silent=silent)
	
	; smooth data to approx. IRTF resolution
	IF KEYWORD_SET(irtf) THEN $
	  err[WHERE(flx GT 0)] = gauss_smooth(err[WHERE(flx GT 0)], resol)

	IF KEYWORD_SET(struct) THEN $
	  data = {wave:wave, flux:flx, err:err} $
	ELSE $
	  data = [[wave],[flx],[err]]
	  
      ENDELSE

      ; return combined data
      IF KEYWORD_SET(comb) THEN $
	RETURN, data

      ; return array of orders (numbered like SpEX data)
      IF KEYWORD_SET(noerr) THEN $
	all_data[0:N_ELEMENTS(wave)-1,0:1,i] = data $
      ELSE $
	all_Data[0:N_ELEMENTS(wave)-1,*,i] = data

    ENDIF ELSE print, "Failed to find files."
    ; end file exist

  ENDFOR
  ; end loop over orders

  RETURN, all_data
  ; returns array of zeros if it failed, or just 0 in the case of merged data

END
