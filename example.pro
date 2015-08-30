PRO EXAMPLE

  ; line definitions
  READCOL, 'linedefs.txt', lineall, f1all,f2all, c1all, c2all, c3all, c4all,  format='A,F,F,F,F,F,F'

  ; standard RV file (K-band only for now)
  std = MRDFITS('spec/J0727+0513_rest.fits',0)

  ; example -- Na at 2.2 microns
  k = 5
  continuum = [[c1all[k], c2all[k]],[c3all[k],c4all[k]]]
  feature = [f1all[k],f2all[k]]


  ; IRTF spectra (unmerged)
  star = 'spec/J0455+0440W_tc.fits'
  data_tc = MRDFITS(star, 0, hdr)

  ; K-band specifics
  order = 0
  wrange = [2.18, 2.41]
  pixscale = 0.00053752800

  ; shift to rest
  ERN_RV, data_tc[*,*,order], std, wrange=wrange, pixscale=pixscale, rv0=rv0, ccorr=1
  ; don't oversample
  lambda = data_tc[*,0,order] - rv0/(3.*10.^5)*data_tc[*,0,order]
  flux = data_tc[*,1,order]
  ew = measure_ew(lambda,flux,continuum,feature)
  print, "SpEX, no oversampling: ", ew

  ; oversample flux
  inc0 = N_ELEMENTS(data_tc[*,0,order])*10.
  lambda0 = REBIN(data_tc[*,0,order],inc0)
  flux0 = REBIN(data_tc[*,1,order],inc0)
  ew = measure_ew(lambda0,flux0,continuum,feature)
  print, "SpEX, oversampling: ", ew


  ; FIRE spectra (merged)
  data = read_fire('J04555445+0440164',dir='spec')
  ERN_RV, data, std, wrange=wrange, pixscale=pixscale, rv0=rv0, ccorr=1
  ew = measure_ew(data[*,0] - rv0/(3.*10.^5)*data[*,0],data[*,1],continuum,feature)
  print, "FIRE, no smoothing: ", ew


  ; resolution degraded
  data_lowres = read_fire('J04555445+0440164',dir='spec', /irtf)
  ew = measure_ew(data[*,0] - rv0/(3.*10.^5)*data[*,0],data_lowres[*,1],continuum,feature)
  print, "FIRE, degraded to SpEX: ", ew


  ; plots for good measure
  pseudo=ew_pseudo(lambda, flux,continuum)
  plot, lambda,flux/pseudo, xrange=[2.18,2.23], yrange=[0.6,1.05]
;   pseudo2=ew_pseudo(data[*,0],data[*,1],continuum)
;   oplot, data[*,0] - rv0/(3.*10.^5)*data[*,0],data[*,1]/pseudo2, co=5
  pseudo3=ew_pseudo(data_lowres[*,0],data_lowres[*,1],continuum)
  oplot, data_lowres[*,0] - rv0/(3.*10.^5)*data[*,0],data_lowres[*,1]/pseudo3, co=4

  oplot,[feature[0],feature[0]],[0,100],color=2
  oplot,[feature[1],feature[1]],[0,100],color=2
  oplot,[continuum[0,0],continuum[0,0]],[0,10],color=3,linestyle=2
  oplot, [continuum[0,1],continuum[0,1]],[0,100],color=3,linestyle=2
  oplot, [continuum[1,1],continuum[1,1]],[0,100],color=3,linestyle=2
  oplot, [continuum[1,0],continuum[1,0]],[0,100],color=3,linestyle=2
  oplot, lambda, pseudo, color=3

END

