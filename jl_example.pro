PRO JL_EXAMPLE
  ;mkct ;load indexed color table
  ; line definitions
  READCOL, 'linedefs.txt', lineall, f1all,f2all, c1all, c2all, c3all, c4all,format='A,F,F,F,F,F,F'
  
  ; standard RV file (K-band only for now)
  std = MRDFITS('spec/J0727+0513_rest.fits',0)
  
  ; IRTF telluric corrected spectra (unmerged)
  star = 'spec/J0455+0440W_tc.fits'
  ;star = '~/EPIC/proc/corrected.fits'
  data_tc = MRDFITS(star, 0, hdr)
  orders = STRSPLIT(SXPAR(hdr,'ORDERS'),',',/extract)
  
  ; K-band specifics
  order = 0
  wrange = [2.18, 2.41] ;optional
  ;pixscale = 0.00053752800 ;order dependent
  pixscale = SXPAR(hdr, STRING(orders[order],FORMAT='("DISPO",I02)')) ; assume IRTF standard header
  
  ; Determine spectral type using N14 relation
  ;merged = mrdfits('~/Dropbox/Merged-New-JJ1629.fits',0)
  index = waterindex(data_tc[*,0,0],data_tc[*,1,0])
  print, 'Water index: ',index
  
  ; shift to rest
  ERN_RV, data_tc[*,*,order], std, wrange=wrange, pixscale=pixscale, rv0=rv0, ccorr=1
  print,''
  print,'Radial Velocity (km/s): ',rv0
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
        ; 4:wrange = [0.82, 1.01]
        ; 4:wrange = [0.95, 1.01]
        4:wrange = [0.82, 0.93]
        else: wrange = [MIN(data_tc[*,0,order]),MAX(data_tc[*,0,order])]
      ENDCASE
      lambda = data_tc[*,0,order] - rv0/(3.*10.^5)*data_tc[*,0,order]
      IF (continuum[0] ge min(wrange)) and (continuum[3] le max(wrange)) THEN BEGIN
        print, lineall[k], ' in Order ',orders[order]
        flux = data_tc[*,1,order]
        eflux = data_tc[*,2,order]
        ew = measure_ew(lambda,flux,continuum,feature)
        print, "EW, no oversampling: ", ew
        
        ; oversample flux
        inc0 = N_ELEMENTS(data_tc[*,0,order])*10.
        lambda0 = REBIN(data_tc[*,0,order],inc0)
        flux0 = REBIN(data_tc[*,1,order],inc0)
        ew = measure_ew(lambda0,flux0,continuum,feature)
        print, "EW, oversampling: ", ew
        print, ''
        
        
        
        ; plots for good measure
        PLT = 1
        PS = 0
        IF PLT THEN BEGIN
          pseudo=ew_pseudo(lambda, flux,continuum)
          xrange=[continuum[0] - 0.01,continuum[3] + 0.01]
          y = flux/pseudo
          y = y[where((lambda ge xrange[0]) and (lambda le xrange[1]))]
          yrange = [MIN(y)*.9,MAX(y)*1.1]
          IF PS THEN BEGIN
            set_plot,'PS'
            device,filen = sxpar(hdr,'OBJECT')+'_Feature_'+lineall[k]+'.eps',/enca,/col
          ENDIF
          title = sxpar(hdr,'OBJECT')+' Feature: '+lineall[k]+ ' EW: '+STRING(ew,FORMAT='(F0.1)')+' A
          plot, lambda,flux/pseudo, xrange=xrange, yrange = yrange,/nodata, $
            title = title, xtitle = sxpar(hdr,'XTITLE'), ytitle = 'Normalized flux'
          oplot, lambda,flux/pseudo, psym = 1, color = 5
          oplot, lambda,flux/pseudo
          ; pseudo2=ew_pseudo(data[*,0],data[*,1],continuum)
          ; oplot, data[*,0] - rv0/(3.*10.^5)*data[*,0],data[*,1]/pseudo2, co=5
          
          oplot,[feature[0],feature[0]],[0,100],color=2
          oplot,[feature[1],feature[1]],[0,100],color=2
          oplot,[continuum[0,0],continuum[0,0]],[0,10],color=3,linestyle=2
          oplot, [continuum[0,1],continuum[0,1]],[0,100],color=3,linestyle=2
          oplot, [continuum[1,1],continuum[1,1]],[0,100],color=3,linestyle=2
          oplot, [continuum[1,0],continuum[1,0]],[0,100],color=3,linestyle=2
          oplot, lambda, pseudo, color=3
          IF PS THEN BEGIN
            device,/close
            set_plot,'X'
          ENDIF
          ;wait,3
        ENDIF
      ENDIF
    ENDFOR
  ENDFOR
  
END
