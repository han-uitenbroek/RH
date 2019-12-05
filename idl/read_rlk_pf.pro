FUNCTION fcn, x, c

  y = c[3] + x*(c[2] + x*(c[1] + x*c[0]))
  return, y
END

FUNCTION fcn0, x

  COMMON data_share, c, pf

  return, fcn(x, c) - pf[0]
END
  
FUNCTION dfcn, x, c, k

  CASE k OF
    0: return, x^3
    1: return, x^2
    2: return, x
    3: return, 1.0
  ENDCASE
END

PRO read_rlk_pf, FIT=fit, PLOT=plot

  COMMON data_share, c, pf

  pf_data_file = 'pf_all_Kurucz.dat'
  ion_pot_file = 'ionpot.input'

  Ntemp = 201L  &  Ncolumn = 6L
  N_periodic_table = 99

  IF (keyword_set(FIT)) THEN BEGIN 
    Nfit = 4
    epsilon = fltarr(Nfit) + 1.0E-3
    Niter = 25
  ENDIF
  openr, lun_in, pf_data_file, /GET_LUN
  openr, lun_ip, ion_pot_file, /GET_LUN, XDR

  pf_output_file = 'pf_Kurucz.input'
  openw, lun_out, pf_output_file, /GET_LUN, /XDR

  pf_str = {pt_index: 0L,  T: 0.0,  pf: dblarr(Ncolumn)}
  pf_data = replicate(pf_str, Ntemp)

  IF (keyword_set(PLOT)) THEN BEGIN
    psopen, XSIZE=7.0, YSIZE=8.0/6.0 * 7.0
    !P.MULTI = [0, 4, 6]
    pswindow, XSIZE=600, YSIZE=1000
    nplot = 0
  ENDIF

  pfi = 0L  &  Nip = 0L
  FOR i = 0, N_periodic_table-1 DO BEGIN
    readf, lun_in, pf_data 
    IF (i EQ 0) THEN writeu, lun_out, Ntemp, double(pf_data.T)

    readu, lun_ip, pfi, Nip
    ionpot = fltarr(Nip)
    readu, lun_ip, ionpot
    IF (pfi NE pf_data[0].pt_index) THEN print, "Mismatch for element:", pfi

    s_max = Ncolumn - 1 
    WHILE (pf_data[0].pf[s_max] EQ 0.0) DO s_max = s_max - 1
    writeu, lun_out, pf_data[0].pt_index, s_max + 1
    writeu, lun_out, transpose(pf_data.pf[0:s_max])
    writeu, lun_out, double(ionpot[0:s_max])

    IF (keyword_set(FIT)) THEN BEGIN 
      FOR s=0, s_max DO BEGIN 
        c = fltarr(Nfit)
        t = alog(5040.0/pf_data.T)
        pf = cmprss(alog(pf_data.pf[s]))

        index = where(pf GT pf[0], count)
        c[0] = pf[0]
        IF (count gt 8) THEN BEGIN
          result = dlsq_fit(t[index], pf[index], c, epsilon, Niter, $
                            FUNCTION_NAME='fcn', DERIVATIVE_NAME='dfcn')
          print, c
          plot, t, pf, /PSYM
          oplot, t, fcn(t, c), COLOR=150

          t0 = fx_root(t[index[0:2]-1], 'fcn0')
          oplot, [t0], fcn([t0], c), PSYM=5, col=125, SYMSIZE=1.5
        ENDIF
      ENDFOR
    ENDIF
    IF (keyword_set(PLOT)) THEN BEGIN
      FOR s=0, s_max DO BEGIN
        plot, pf_data.T, pf_data.pf[s], /YLOG, /XLOG        
        rhannotate, xann(0.05), yann(0.9), TEXT="element:" + $
         string(FORMAT='(I3)', pf_data.pt_index), CHARSIZE=0.7
        rhannotate, xann(0.05), yann(0.8), TEXT="ion num:" + $
         string(FORMAT='(I3)', s), CHARSIZE=0.7
        nplot = nplot + 1
        IF (nplot MOD 24 EQ 0) THEN stop
      ENDFOR
    ENDIF
  ENDFOR

  IF (keyword_set(PLOT)) THEN psclose
  free_lun, lun_in, lun_ip, lun_out
  !P.MULTI = 0
END
  
