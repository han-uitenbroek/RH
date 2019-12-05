FUNCTION getWavelengthTable, files, waveRange

  ;; Reads line transitions from files specified in files parameter
  ;; and returns lines only in specified wavenumber region.
  ;; Input parameters:

  ;;   files     -- String array of filenames
  ;;   waveRange -- Range of wavenumbers [cm^-1]

  molecular_rt = {wavenumber: 0.0, dipole: 0.0, Aji: 0.0, Ei: 0.0, gf: 0.0, $
                  strength: 0.0, vj: 0, vi: 0, type: 'x', Ji: 0, isotope: 0}
  line_format  = '(F9.4, 2E10.3, F11.4, 2E10.3, 2I3, A2, I4, I3)'

  Nfiles = n_elements(files)

  Nmrt = 0L
  FOR n=0, Nfiles-1 DO BEGIN
    wtdum = readmollines(files[n])
    IF (n EQ 0) THEN  wtable = wtdum  ELSE  wtable = [wtable, wtdum]
  ENDFOR

  indices = where(wtable.wavenumber GE waveRange[0]  AND $
                  wtable.wavenumber LE waveRange[1], Nmrt)
  IF (Nmrt GT 0) THEN wtable = wtable[indices]  ELSE  wtable = 0

  return, wtable
END

PRO cospectrum, spectrum, wavenumber, wtable, $
                XRANGE=xrange, MARK=mark, NOPLOT=noplot, IMAXTOP=imaxtop, $
                W_ATLAS=w_atlas, I_ATLAS=I_atlas, _EXTRA=extra

  C12_O16_FUNDAMENTAL = '/home/uitenbr/src/rh/Molecules/CO/goorvitch.dv1_26'
  C13_O16_FUNDAMENTAL = '/home/uitenbr/src/rh/Molecules/CO/goorvitch.dv1_36'


  IF (n_params(0) LT 1) THEN BEGIN
    print, "Usage: cospectrum, spectrum, wave, wtable, $"
    print, "                   XRANGE=xrange, MARK=mark, NOPLOT=noplot, $"
    print, "                   IMAXTOP=imaxtop"
    return
  ENDIF

  CLIGHT     = 2.99792458E+08
  HPLANCK    = 6.626176E-34
  KBOLTZMANN = 1.380662E-23

  NM_TO_M    = 1.0E-09

  lambda = 1.0E7 / wavenumber
  IF (NOT keyword_set(XRANGE)) THEN $
   range = [wavenumber(0), wavenumber(n_elements(lambda)-1)] $
  ELSE $
   range = xrange
  wrange = [min(range), max(range)]

  IF (NOT keyword_set(IMAXTOP)) THEN imaxtop = 1.15

  imax = max(spectrum, MIN=imin)
  irange = [0.925*imin, imaxtop*imax]

  lambda_avg = 0.5*(lambda[0] + lambda[n_elements(lambda)-1]) * NM_TO_M
  trange = (HPLANCK*CLIGHT) / ((lambda_avg*KBOLTZMANN) * $
             alog(1.0 + (2.0*HPLANCK*CLIGHT)/(lambda_avg^3 * irange)))

  IF (keyword_set(MARK)  AND  (n_tags(wtable) EQ 0)) THEN BEGIN
    wtable = getWavelengthTable([C12_O16_FUNDAMENTAL, $
                                 C13_O16_FUNDAMENTAL], wrange)
    order = sort(wtable.wavenumber)
    wtable = wtable[order]
  ENDIF
  IF (keyword_set(NOPLOT)) THEN return

  xtl = 'Wavelength [nm]'  &  xtw = 'Wave number [cm!U-1!N]'
  yt = 'Intensity [J s!U-1!N m!U-2!N Hz!U-1!N sr!U-1!N]'

  plot, wavenumber, spectrum, XTITLE=xtw, YTITLE=yt, $
   XRANGE=wrange, XSTYLE=9, YRANGE=irange, YSTYLE=9, $
   XMARGIN=[14, 10], YMARGIN=[5, 5], /YNOZERO, /NODATA

  oplot, wavenumber, spectrum, _EXTRA=extra

  IF (keyword_set(W_ATLAS) AND keyword_set(I_ATLAS)) THEN $
   oplot, w_atlas, I_atlas, COLOR=16B, _EXTRA=extra

  axis, XRANGE=1.0E7/wrange, /XAXIS, XSTYLE=1, XTITLE=xtl
  axis, yrange=trange, /YAXIS, YSTYLE=1, YTITLE='Brightness Temperature [K]'

  IF (keyword_set(MARK)  AND  n_tags(wtable) GT 0) THEN BEGIN
    ymark1 = 1.05*imax
    ymark2 = 1.03*imax
    ytext  = 1.06*imax

    Nmrt = n_elements(wtable)
    dx0 = 0.05
    dx  = (1.0 - 2*dx0) / (Nmrt-1)  

    FOR n=0, Nmrt-1 DO BEGIN
      x = wtable(n).wavenumber
      label = string(FORMAT='(I2, "-", I2, A2, I3)', $
                     wtable(n).vj, wtable(n).vi, $
                     wtable(n).type, wtable(n).Ji)
      IF (wtable(n).isotope EQ 36) THEN label = label + ' [13]'
      tabinv, wavenumber, x, weff
      ym = linear(spectrum, weff)
      oplot, [x, x], [ymark2, 0.4*ymark2 + 0.6*ym], COLOR=100B
      oplot, [xann(dx0+n*dx), x], [ymark1, ymark2], COLOR=100B
      xyouts, xann(dx0+n*dx), ytext, label, ORIENT=89.99, $
       COLOR=100B, CHARSIZE=0.7 ;;;0.833
    ENDFOR
  ENDIF

  return
END

PRO create_wavelength_input, wtable, $
                             DELTAL=deltal, NPOINT=Npoint, FILENAME=filename

  IF (n_params(0) LT 1) THEN BEGIN
    print, "Usage: create_wavelength_input, wtable, $"
    print, "                 DELTAL=deltal, NPOINT=Npoint, FILENAME=filename"
    return
  ENDIF

  IF (NOT keyword_set(DELTAL))   THEN deltal = 0.25
  IF (NOT keyword_set(NPOINT))   THEN Npoint = 30
  IF (NOT keyword_set(FILENAME)) THEN filename = 'wavelength.input'

  lambda_central = 1.0E+07 / wtable.wavenumber
  Nline = n_elements(wtable)

  lambda = fltarr(Npoint * Nline)
  FOR n=0, Nline-1 DO $
   lambda[n*Npoint:(n+1)*Npoint-1] = lambda_central[n] - deltal/2.0 + $
                                     deltal * findgen(Npoint) / (Npoint-1)

  lambda = lambda(sort(lambda))

  openw, unit, filename, /GET_LUN
  writeu, unit, n_elements(lambda), lambda
  free_lun, unit

END
