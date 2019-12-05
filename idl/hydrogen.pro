FUNCTION Gaunt_bf, E, quantumNo

  COMMON const, ERydb, hPlanck, cLight, qElect, mElect, epsil0

  x    = (ERydb - E) / ERydb
  x3   = x^0.33333333
  nsqx = 1.0 / (quantumNo^2 * x)

  return,  1.0 + 0.1728*x3 * (1.0 - 2.0*nsqx) - $
               0.0496*x3^2 * (1.0 - (1.0 - nsqx)*0.66666667*nsqx)
END

PRO hydrogen, fileName, NLEVEL=Nlevel, NLINE=Nline

  COMMON const, ERydb, hPlanck, cLight, qElect, mElect, epsil0

  IF (NOT keyword_set(NLEVEL)) THEN Nlevel = 16
  IF (NOT keyword_set(NLINE))  THEN Nline  = 0

  ERydb     = 1.0967875E+7
  ERydb_inf = 1.09737E+7
  hPlanck = 6.626176E-34
  cLight  = 2.99792458E+08
  mElect  = 9.109534E-31
  qElect  = 1.60219E-19
  epsil0  = 8.85419E-12 

  E = fltarr(Nlevel)  &  E(Nlevel-1) = Erydb
  E(0:Nlevel-2) = (1.0 - 1.0/(1.0+findgen(Nlevel-1))^2) * Erydb

  labels = ["H I 1S 2SE          ", $
            "H I 2P 2PO          ", $
            "H I 3D 2DE          ", $
            "H I 4F 2FO          ", $
            "H I 5G 2GE          ", $
            "H I 6               ", $
            "H I 7               ", $
            "H I 8               ", $
            "H I 9               ", $
            "H I 10              ", $
            "H I 11              ", $
            "H I 12              ", $
            "H I 13              ", $
            "H I 14              ", $
            "H I 15              ", $
            "H II GROUND TERM    " ]

  g = fltarr(Nlevel)  &  g(Nlevel-1) = 1.0
  g(0:Nlevel-2) = 2.0*(1.0 + findgen(Nlevel-1))^2

  stage = intarr(Nlevel)  &  stage(Nlevel-1) = 1

  sigma0 = 32.0/3.0^1.5 * (qElect/sqrt(4.0*!PI*epsil0))^2 / $
    (mElect * cLight) / (2.0*ERydb_inf*cLight)

  alpha = fltarr(Nlevel-1)
  FOR i=0,Nlevel-2 DO alpha(i) = $
   sigma0 * (i + 1) * Gaunt_bf(E(i), i + 1)

  openw, Hunit, fileName, /GET_LUN
  printf, Hunit, FORMAT='("* Hydrogen", /"  H")

  abundance = 12.0  &  atomWeight = 1.008
  printf, Hunit, abundance, atomWeight, $
   FORMAT='(/"* abundance[Dex]     weight[a.u.]", /4X, F5.2, 13X, F7.3)'
  
  printf, Hunit, Nlevel, Nline, Nlevel-1, $
   FORMAT='(/"* Nlevel  Nline   Ncont", /2X, I3, 5X, I3, 5X, I3)'

  printf, Hunit, $
   FORMAT='(/"*  E[cm^-1]    g           label[20]        stage")'
  printf, Hunit, FORMAT='("*                     |----|----|----|----")'

  FOR i=0, Nlevel-1 DO BEGIN
    printf, Hunit, E(i)*1.0E-2, g(i), labels(i), stage(i), $
     FORMAT='(F11.3, 1X, F6.2, 4X, A20, 3X, I2)'
  ENDFOR

  printf, Hunit, $
   FORMAT='(/"* j   i       f     type  Nlambda symmetr  qcore  qwing' + $
          '  vdWapprx           vdWaals         radiative  Starck", ' + $
          '/"*                                                      ' + $
          '                 H          He")'

  printf, Hunit, FORMAT=$
   '(/"* j   i   alpha [m^2]   Nlambda     Wavel. Dep.   lamb_min [nm]", /"*")'

  FOR i=0, Nlevel-2 DO BEGIN
    printf, Hunit, Nlevel-1, i, alpha(i), 1, "HYDROGENIC", 50.0, $
     FORMAT='(X, I2, 2X, I2, 4X, E10.4, 5X, I2, 8X, A, 7X, F5.1)'
  ENDFOR

  free_lun, Hunit
END
