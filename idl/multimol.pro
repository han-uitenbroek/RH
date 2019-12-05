FUNCTION eqcCO, temp

  EV = 1.60219E-19
  KBOLTZMANN = 1.380662E-23
  CM_TO_M = 1.0E-2

  Ediss = 11.091 * EV
  kT = KBOLTZMANN * temp
  a = reverse([-49.0414, 14.0306, -26.6341, 35.3827, -26.5424, 8.32385])

  t = temp*1.0E-4
  eqc = a(0)
  FOR n=1, n_elements(a)-1 DO eqc = eqc*t + a(n)

  return, exp(Ediss/kT + eqc - 1.5*alog(temp)) * CM_TO_M^3
END

FUNCTION eqcCH, temp

  EV = 1.60219E-19
  KBOLTZMANN = 1.380662E-23
  CM_TO_M = 1.0E-2

  Ediss = 3.470 * EV
  kT = KBOLTZMANN * temp
  a = reverse([-4.5506E+01, 1.7112E-03, -3.6319E-07, 5.0164E-11, -2.8716E-15])

  eqc = a(0)
  FOR n=1, n_elements(a)-1 DO eqc = eqc*temp + a(n)

  return, exp(Ediss/kT + eqc - 1.5*alog(temp)) * CM_TO_M^3
END

FUNCTION eqcCN, temp

  EV = 1.60219E-19
  KBOLTZMANN = 1.380662E-23
  CM_TO_M = 1.0E-2

  Ediss = 8.109 * EV
  kT = KBOLTZMANN * temp
  a = reverse([-4.7853E+01, 1.8656E-03, -4.6185E-07, 7.1497E-11, -4.3750E-15])

  eqc = a(0)
  FOR n=1, n_elements(a)-1 DO eqc = eqc*temp + a(n)

  return, exp(Ediss/kT + eqc - 1.5*alog(temp)) * CM_TO_M^3
END

FUNCTION fn, n

  COMMON EQC, eCO, eCH, eCN

  f = [n(0) + n(4) + n(5) + n(6), $
       n(1) + n(4), $
       n(2) + n(5), $
       n(3) + n(6), $
       n(4) - n(0)*n(1)*eCO, $
       n(5) - n(0)*n(2)*eCN, $
       n(6) - n(0)*n(3)*eCH]

  return, f
END

FUNCTION dfdn, n

  COMMON EQC, eCO, eCH, eCN

  dfdn = [[1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0], $
          [0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0], $
          [0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0], $
          [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0], $
          [-n(1)*eCO, -n(0)*eCO, 0.0, 0.0, 1.0, 0.0, 0.0], $
          [-n(2)*eCN, 0.0, -n(0)*eCN, 0.0, 0.0, 1.0, 0.0], $
          [-n(3)*eCH, 0.0, 0.0, -n(0)*eCH, 0.0, 0.0, 1.0]]

  return, dfdn
END

PRO multimol, temp, nH, nC, nO, nN, nHI, nCO, nCN, nCH

  COMMON EQC, eCO, eCH, eCN

  ;; Toy program to solve equilibrium C, O, N, HI, CO, CN, and CH at
  ;; temperature temp and total hydrogen density nH

  Nmax_iter = 8

  ;; Abundances of Carbon and Oxygen

  A_C = 3.958E-04
  A_O = 8.451E-04
  A_N = 9.930E-05

  eCO = eqcCO(temp)
  eCN = eqcCN(temp)
  eCH = eqcCH(temp)
  a = [A_C*nH, A_O*nH, A_N*nH, nH, 0.0, 0.0, 0.0]

  ;; Initial solution without any molecules.

  n = a  &  nold = n

  ;; Iterate with Newton-Raphson.

  FOR niter=1, Nmax_iter DO BEGIN
    df = dfdn(n)
    ludc, df, index
    dn = lusol(df, index, fn(n) - a)
    n = n - dn

    IF (niter GT 1) THEN BEGIN
      print, -dn/n, FORMAT='(7(E9.2, 1X))'
      IF (max(abs(dn/nold)) LT 1.0E-4) THEN GOTO, break
    ENDIF
    nold = n
  ENDFOR

break:
  nC  = n(0)
  nO  = n(1)
  nN  = n(2)
  nHI = n(3)
  nCO = n(4)
  nCN = n(5)
  nCH = n(6)
END

@ATMOS

  Nt = 51
  T0 = 2.0E3  &  TN = 6.0E3
  T = T0 + (TN - T0) * findgen(Nt)/(Nt-1)

  Ndens = 50
  nH0 = 16.0  &  nHN = 24.0
  nHtot = double(10.0^(nH0 + (nHN - nH0) * findgen(Ndens)/(Ndens-1)))

  CO   = fltarr(Nt, Ndens)  &  ch = co
  chco = fltarr(Nt, Ndens)  &  cnco = chco
  FOR j=0, Ndens-1 DO BEGIN
    FOR i=0, Nt-1 DO BEGIN
      multimol, T(i), nHtot(j), nC, nO, nN, nHI, nCO, nCN, nCH
      co(i, j) = nCO
      ch(i, j) = nCH
      chco(i, j) = nCH / nCO
      cnco(i, j) = nCN / nCO
    ENDFOR
  ENDFOR

  read_atmos, '~/src/rh_v2/Atmos/FALC_82.atmos'

  nHfal = fltarr(Ndep)
  FOR k=0, Ndep-1 DO $
    nHfal(k) = (total(nH(k, *)) + np(k)) * 1.0E+6

  index = where(temp GE 1.0E+3  AND  temp LE 6.0E+3) 
  Tfal = temp(index)
  nHfal = double(nHfal(index))

  COfal = nHfal  &  CHfal = COfal  &  CNfal = COfal
  FOR k=0, n_elements(index)-1 DO BEGIN
    multimol, Tfal(k), nHfal(k), nC, nO, nN, nHI, nCO, nCN, nCH
    COfal(k) = nCO
    CHfal(k) = nCH
    CNfal(k) = nCN
  ENDFOR

  MAKE_PS = 1

  IF (MAKE_PS) THEN $
   psopen, FILENAME='molecule.ps', /COLOR, /FONT

  loadct, 1
  !P.MULTI = [0, 2, 2, 0, 0]  &  !P.CHARSIZE=1.4

  shade_surf, co, T, nHtot, /YLOG, /ZLOG, $
   XTITLE='Temperature [K]', YTITLE='Hydrogen density [m^-3]', $
   ZTITLE='n!DCO!N [m^-3]', /SAVE
  plots, Tfal, nHfal, COfal, /T3D, THICK=2, COLOR=250B

  shade_surf, ch, T, nHtot, /YLOG, /ZLOG, $
   XTITLE='Temperature [K]', YTITLE='Hydrogen density [m^-3]', $
   ZTITLE='n!DCH!N [m^-3]', /SAVE
  plots, Tfal, nHfal, CHfal, /T3D, THICK=2, COLOR=250B

  shade_surf, chco, T, nHtot, /YLOG, /ZLOG, $
   XTITLE='Temperature [K]', YTITLE='Hydrogen density [m^-3]', $
   ZTITLE='n!DCH!N/n!DCO!N', /SAVE
  plots, Tfal, nHfal, CHfal/COfal, /T3D, THICK=2, COLOR=250B

  shade_surf, cnco, T, nHtot, /YLOG, /ZLOG, $
   XTITLE='Temperature [K]', YTITLE='Hydrogen density [m^-3]', $
   ZTITLE='n!DCN!N/n!DCO!N', /SAVE
  plots, Tfal, nHfal, CNfal/COfal, /T3D, THICK=2, COLOR=250B

  IF (MAKE_PS) THEN $
   psclose

  !P.MULTI=0  &  !P.CHARSIZE=1.0
END
