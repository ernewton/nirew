;============================================
; test code to demonstrate use in loops
; developed by J.A. Lewis

PRO JL_EXAMPLE, showplot=showplot, ps=ps

  IF N_ELEMENTS(ps) EQ 0 THEN ps = 0
  IF N_ELEMENTS(showplot) EQ 0 THEN showplot = 1
  ccorr='c_correlate'
  contf=0

  ; line definitions
  READCOL, 'linedefs.txt', lineall, f1all,f2all, c1all, c2all, c3all, c4all,format='A,F,F,F,F,F,F', /silent
  
  ; standard RV file (K-band only for now)
  std_tc = MRDFITS('spec/J0727+0513_rest.fits',0, shdr, /silent)
  wlcal=1
  atrest=1
  
  ; IRTF telluric corrected spectra (unmerged)
  star = 'J0200+1303'
  data_tc = MRDFITS('spec/'+star+'_tc.fits', 0, hdr, /silent)
  data = MRDFITS('spec/'+star+'.fits', 0, hdr, /silent)
  data_tc[*,0,*] = data[*,0,*]
  orders = STRSPLIT(SXPAR(hdr,'ORDERS'),',',/extract)
  
  ; Determine spectral type using N14 relation
  index = water_index(data_tc[*,0,0],data_tc[*,1,0])
  print, 'Water index: ',index

  ; shift to rest
  rv_arr = fltarr(6)
  rv0_arr = fltarr(6)
  mshft_arr = fltarr(6)
  FOR order=0,5 DO BEGIN
    CASE order OF
      0:wrange = [2.18, 2.41]
      1:wrange = [1.49, 1.73]
      2:wrange = [1.15,1.32]
      3:wrange = [1.0, 1.1]
      4:wrange = [0.82, 0.93]
      else: wrange = [MIN(data_tc[*,0,order]),MAX(data_tc[*,0,order])]
    ENDCASE
    order_variables, hdr, order, wrange, trange, pixscale, polydegree, instrument="spex"
    dtc = data_tc[*,*,order]
    d = data[*,*,order]
    stc = std_tc[*,*,order]
;     s = std[*,*,order]
    NIR_RV, dtc,hdr, d, $
	stc,shdr, $
	wlcal=wlcal, atrest=atrest, stdrv=stdrv, $ ; already wavelength calibrated?
	pixscale=pixscale, polydegree=polydegree, $
	spixscale=pixscale, spolydegree=polydegree, $ ; standard is from same set-up 
	wrange=wrange, trange=trange, $
	mshft=myshft, shftarr=shftarr, torest=myrv0, rv=myrv, $
	quiet=1, contf=contf, ccorr=ccorr
    rv_arr[order] = myrv ; measured absolute RV
    rv0_arr[order] = myrv0 ; shift to zero velocity
    mshft_arr[order] = myshft ; shift to abs wavelength
  ENDFOR
  print,''
  print,'Radial Velocity (km/s): ',rv_arr 
  print,''

  ; Run for all lines
  FOR k=0,N_ELEMENTS(lineall)-1 DO BEGIN
  
    continuum = [[c1all[k], c2all[k]],[c3all[k],c4all[k]]]
    feature = [f1all[k],f2all[k]]
    
    FOR order=0,N_ELEMENTS(orders)-1 DO BEGIN
      CASE order OF
        0:wrange = [2.18, 2.41]
        1:wrange = [1.49, 1.73]
        2:wrange = [1.15,1.32]
        3:wrange = [1.0, 1.1]
        4:wrange = [0.82, 0.93]
        else: wrange = [MIN(data_tc[*,0,order]),MAX(data_tc[*,0,order])]
      ENDCASE

      lambda = (data_tc[*,0,order]+mshft_arr[order])*(1-rv0_arr[order]/3.e5)
      IF (continuum[0] ge min(wrange)) and (continuum[3] le max(wrange)) THEN BEGIN
        print, lineall[k], ' in Order ',orders[order]
        flux = data_tc[*,1,order]
        eflux = data_tc[*,2,order]
        ew = measure_ew(lambda,flux,continuum,feature)
        print, "EW, no oversampling: ", ew
        
        ; oversample flux
        inc0 = N_ELEMENTS(data_tc[*,0,order])*10.
        lambda0 = REBIN(lambda,inc0)
        flux0 = REBIN(data_tc[*,1,order],inc0)
        ew = measure_ew(lambda0,flux0,continuum,feature)
        print, "EW, oversampling: ", ew
        print, ''
      
        ; plots for good measure
        IF KEYWORD_SET(showplot) THEN BEGIN
          pseudo=ew_pseudo(lambda, flux,continuum)
          xrange=[continuum[0] - 0.01,continuum[3] + 0.01]
          y = flux/pseudo
          y = y[where((lambda ge xrange[0]) and (lambda le xrange[1]))]
          yrange = [MIN(y)*.9,MAX(y)*1.1]
          IF KEYWORD_SET(ps) THEN BEGIN
            set_plot,'PS'
            device,filen = sxpar(hdr,'OBJECT')+'_Feature_'+lineall[k]+'.eps',/enca,/col
          ENDIF
          title = sxpar(hdr,'OBJECT')+' Feature: '+lineall[k]+ ' EW: '+STRING(ew,FORMAT='(F0.1)')+' A
          plot, lambda,flux/pseudo, xrange=xrange, yrange = yrange,/nodata, $
            title = title, xtitle = sxpar(hdr,'XTITLE'), ytitle = 'Normalized flux'
          oplot, lambda,flux/pseudo, psym = 1
          oplot, lambda,flux/pseudo
          
          oplot,[feature[0],feature[0]],[0,100]
          oplot,[feature[1],feature[1]],[0,100]
          oplot,[continuum[0,0],continuum[0,0]],[0,10], linestyle=2
          oplot, [continuum[0,1],continuum[0,1]],[0,100],linestyle=2
          oplot, [continuum[1,1],continuum[1,1]],[0,100],linestyle=2
          oplot, [continuum[1,0],continuum[1,0]],[0,100],linestyle=2
          oplot, lambda, pseudo, color=3
          IF PS THEN BEGIN
            device,/close
            set_plot,'X'
          ENDIF ELSE $
	    wait,1
        ENDIF
      ENDIF
    ENDFOR
  ENDFOR
  
END
