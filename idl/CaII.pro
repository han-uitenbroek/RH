PRO getColl4s_4p, E, Q, PLOT=plot

  COMMON constant, cLight,  hPlanck, kBoltz, rBohr, ERydb, mElect, $
   one_eV, pia02

  ;; --- Collisional strength 4S1/2 - 4P3/2

  ;; --- From J. Mitroy, D.C. Griffin, D.W. Norcross 1988,
  ;;     Phys. Rev. 38, 3339-3350. First Fig. 1 (E in Rydb, Q in 10^-20 m^2)

  E1 = [0.24, 0.28, 0.31, 0.325, 0.34, 0.37, 0.38, 0.39, 0.40, 0.44, 0.48, $
       0.50, 0.58, 0.64, 0.68]
  Q1 = [17.8, 18.0, 17.9, 19.5,  17.7, 17.0, 16.3, 17.0, 19.8, 16.5, 17.7, $
       16.0, 14.7, 14.4, 14.2] * 1.0E-20

  ;; --- Ibidem Fig. 2

  E2 = [1.0,  1.2,  1.5,  2.0]
  Q2 = [13.2, 12.7, 12.0, 10.6] * 1.0E-20

  ;; --- Higher energy data from Alfred Z. Msezane 1988,
  ;;     J. Phys. B 21, L61-L66. Figure 3 (E in eV, Q in 10^-20 m^2)

  E3 = [40.0, 60.0, 100.0, 1000.0] * one_eV
  Q3 = [ 9.1,  7.2,   4.8,    1.0] * 1.0E-20

  E = [E1, E2, E3]
  Q = [Q1, Q2, Q3] / pia02
  s = sort(E)  &  E = E(s)  &  Q = Q(s)

  IF (keyword_set(PLOT)) THEN BEGIN
    plot_oi, E, Q, /NODATA, XTITLE='Energy [Rydb]', $
     YTITLE='Cross Section 4S!D1/2!N - 4P!D3/2!N [!7p!xa!U2!N]', CHARSIZE=1.4

    oplot, E1, Q1/pia02, PSYM=1
    oplot, E2, Q2/pia02, PSYM=2
    oplot, E3, Q3/pia02, PSYM=4
  ENDIF
END

PRO getColl4p_4p, E, Q, PLOT=plot

  COMMON constant, cLight,  hPlanck, kBoltz, rBohr, ERydb, mElect, $
   one_eV, pia02

  ;; --- Collisional strength fine structure 4P1/2 - 4P3/2
  ;;     From H. E. Saraph 1970, J. Phys. B 3, 952-958
  ;;     Energies in Rydberg, cross_section in units of pi * a_0^2

  ;; --- From table 2.

  gi = 2.0
  E1 = [0.25, 0.35, 0.50, 0.70]
  Q1 = [10.1, 11.1, 11.7, 11.5] / (gi * E1)

  ;; --- Next values are extrapolated, assuming that Ohmega as defined
  ;;     by Saraph remains constant which seems to be a good approximation
  ;;     for the fine structure collisional strengths.

  E2 = [1.0, 1.2, 1.6, 2.0, 2.5, 4.0]
  Q2 = fltarr(6) + 11.5 / (gi * E2)

  E = [E1, E2]
  Q = [Q1, Q2]

  IF (keyword_set(PLOT)) THEN BEGIN
    plot_oi, E, Q, /NODATA, XTITLE='Energy [Rydb]', $
     YTITLE='Cross Section 4P!D1/2!N - 4P!D3/2!N [!7p!xa!U2!N]', CHARSIZE=1.4

    oplot, E1, Q1, PSYM=1
    oplot, E2, Q2, PSYM=2

  ENDIF
END

PRO getColl3d_3d, E, Q, PLOT=plot

  COMMON constant, cLight,  hPlanck, kBoltz, rBohr, ERydb, mElect, $
   one_eV, pia02

  ;; --- Collisional strength fine structure 4D3/2 - 4D5/2
  ;;     From H. E. Saraph 1970, J. Phys. B 3, 952-958
  ;;     Energies in Rydberg, cross_section in units of pi * a_0^2

  ;; --- From table 2.

  gi = 4.0
  E1 = [0.25, 0.35, 0.50, 0.70]
  Q1 = [12.2, 13.0, 13.6, 13.8] / (gi * E1)

  ;; --- Next values are extrapolated, assuming that Ohmega as defined
  ;;     by Saraph remains constant which seems to be a good approximation
  ;;     for the fine structure collisional strengths.

  E2 = [1.0, 1.2, 1.6, 2.0, 2.5, 4.0]
  Q2 = fltarr(6) + 13.8 / (gi * E2)

  E = [E1, E2]
  Q = [Q1, Q2]

  IF (keyword_set(PLOT)) THEN BEGIN
    plot_oi, E, Q, /NODATA, XTITLE='Energy [Rydb]', $
     YTITLE='Cross Section 3D!D3/2!N - 3D!D5/2!N [!7p!xa!U2!N]', CHARSIZE=1.4

    oplot, E1, Q1, PSYM=1
    oplot, E2, Q2, PSYM=2

  ENDIF
END

PRO getColl4s_3d, E, Q, PLOT=plot

  COMMON constant, cLight,  hPlanck, kBoltz, rBohr, ERydb, mElect, $
   one_eV, pia02

  ;; --- Collisions 4s - 3d --

  ;; --- From J. Mitroy, D.C. Griffin, D.W. Norcross 1988,
  ;;     Phys. Rev. 38, 3339-3350. First Fig. 5 (E in Rydb, Q in 10^-20 m^2)

  E1 = [0.142, 0.164, 0.168, 0.182, 0.190, 0.20,  0.218, 0.222, 0.268]
  Q1 = [20.0,  20.0,  72.0,  10.0,  46.0,  14.0,  40.0,  11.0,  9.0] * 1.0E-20

  ;; --- Higher energy data from Alfred Z. Msezane 1988,
  ;;     J. Phys. B 21, L61-L66. Figure 1 (E in eV, Q in 10^-20 m^2)

  E2 = [20.0, 40.0, 60.0, 80.0, 100.0] * one_eV
  Q2 = [1.8,  0.8,  0.58, 0.4,  0.31] * pia02

  ;; See also: Chidichimo 1981, J. Phys. B 14, 4149-4164, Figure 2
  ;; Energy in Rydberg, cross section in 1.0E-20 m^2.

  E3 = [0.25, 0.40, 0.50, 0.70, 0.90, 1.30, 1.7,  2.1,  2.5]
  Q3 = [10.0, 7.0,  4.8,  3.7,  3.0,  2.0,  1.6,  1.3,  1.2] * 1.0E-20

  E = [E1, E2, E3]
  Q = [Q1, Q2, Q3] / pia02
  s = sort(E)  &  E = E(s)  &  Q = Q(s)

  IF (keyword_set(PLOT)) THEN BEGIN
    plot_oo, E, Q, /NODATA, XTITLE='Energy [Rydb]', $
     YTITLE='Cross Section 4s-3d [!7p!xa!U2!N]', CHARSIZE=1.4

    oplot, E1, Q1/pia02, PSYM=1
    oplot, E2, Q2/pia02, PSYM=2
    oplot, E3, Q3/pia02, PSYM=4
  ENDIF
END

PRO getColl3d_4p, E, Q, PLOT=plot

  COMMON constant, cLight,  hPlanck, kBoltz, rBohr, ERydb, mElect, $
   one_eV, pia02

  ;; --- Collisions 3d - 4p --

  ;; Chidichimo 1981, J. Phys. B 14, 4149-4164, Figure 3
  ;; Energy in Rydberg, cross section in 1.0E-20 m^2.

  E = [0.12, 0.2,  0.3,  0.5,  0.7,  1.0,  1.4,  1.8,  2.2]
  Q = [28.0, 21.0, 16.0, 12.2, 10.0, 8.0,  6.2,  5.0,  4.1] * 1.0E-20
  Q = Q / pia02

  ;; Fine structure collision strengths from H. E. Saraph,
  ;; J Phys. B 3, 952-958

  k2 = [0.25, 0.35, 0.5, 0.7]
  Q13 = [16.6, 22.0, 27.6, 32.2]/(4.0*k2)
  Q14 = [7.42, 8.07, 8.69, 9.17]/(4.0*k2)
  Q23 = [3.41, 3.05, 2.65, 2.30]/(6.0*k2)
  Q24 = [32.6, 42.0, 51.7, 59.7]/(6.0*k2)

  c14 = (alog(Q14(0))*alog(k2(3)) - alog(Q14(3))*alog(k2(0))) / $
   (alog(k2(3)) - alog(k2(0)))
  alph14 = (alog(Q14(3)) - alog(Q14(0))) / (alog(k2(3)) - alog(k2(0)))
  Qp14 = exp(c14) * E^alph14

  c23 = (alog(Q23(0))*alog(k2(3)) - alog(Q23(3))*alog(k2(0))) / $
   (alog(k2(3)) - alog(k2(0)))
  alph23 = (alog(Q23(3)) - alog(Q23(0))) / (alog(k2(3)) - alog(k2(0)))
  Qp23 = exp(c23) * E^alph23

  Qp13 = (Q - 0.6*Qp23)
  Qp24 = (1.5*Q - Qp14)
  IF (keyword_set(PLOT)) THEN BEGIN
    plot_oo, E, Q, /NODATA, XTITLE='Energy [Rydb]', $
     YTITLE='Cross Section 3d-4p [!7p!xa!U2!N]', CHARSIZE=1.4
    oplot, E, Q, PSYM=1

    oplot, k2, Q13  &  rhannotate, k2(1), Q13(1), text='Q13'
    oplot, k2, Q14  &  rhannotate, k2(1), Q14(1), text='Q14'
    oplot, k2, Q23  &  rhannotate, k2(1), Q23(1), text='Q23'
    oplot, k2, Q24  &  rhannotate, k2(1), Q24(1), text='Q24'

    oplot, E, Qp14, COLOR=175B
    oplot, E, Qp23, COLOR=175B
    oplot, E, Qp13, COLOR=175B
    oplot, E, Qp24, COLOR=175B
  ENDIF
  NEnergy = n_elements(E)  &  Q = fltarr(NEnergy, 4)
  Q(*,0) = Qp13  &  Q(*,1) = Qp14  &  Q(*,2) = Qp23  &  Q(*,3) = Qp24
END

PRO CaII

  COMMON constant, cLight,  hPlanck, kBoltz, rBohr, ERydb, mElect, $
   one_eV, pia02

  cLight  = 2.99792458E+08
  hPlanck = 6.626176E-34
  kBoltz  = 1.380662E-23
  rBohr   = 5.29177E-11
  ERydb   = 2.17991E-18
  mElect  = 9.109534E-31
  one_eV  = 1.0 / 13.605         ;; in Rydberg
  pia02   = 8.798E-21            ;; in m^2

  ;; Wavelengths from C. E. Moore, M. G. J. Minnaert, & J. Houtgast,
  ;; "The solar spectrum 2935 /AA/ to 8770 /AA/"

  ;;   3933.682  --  4S_1/2 - 4P_3/2
  ;;   3968.492  --  4S_1/2 - 4P_1/2
  ;;   8498.062  --  4P_3/2 - 3D_5/2
  ;;   8542.144  --  4P_3/2 - 3D_3/2
  ;;   8662.170  --  4P_1/2 - 3D_3/2

  Nlevel = 6  &  Nline = 5

  E = [0.0, 13654.062, 13714.813, 25198.488, 25421.477, 95785.470]
  g = [2.0, 4.0, 6.0, 2.0, 4.0, 1.0]

  labels = ['CA II 3P6 4S 2SE    ', $
            'CA II 3P6 3D 2DE 3  ', $
            'CA II 3P6 3D 2DE 5  ', $
            'CA II 3P6 4P 2PO 1  ', $
            'CA II 3P6 4P 2PO 3  ', $
            'CA III GROUND TERM  ']

  Z = 2
  EJoule = E*1.0E2 * (hPlanck * cLight)
  n_eff = sqrt(Z^2 * ERydb / (EJoule(Nlevel-1) - EJoule(0:Nlevel-2)))

  T = [3000.0, 5000.0, 7000.0, 15000.0, 50000.0, 100000.0]
  Ntemp = n_elements(T)

  ohm = fltarr(Ntemp)
  fmt = '(" OMEGA", 2X, I2, X, I2, 2X, 6E11.3)'

  E_edge = (g(1)*E(1) + g(2)*E(2))/(g(1) + g(2)) * 1.0E2 * $
   (hPlanck * cLight) / ERydb
  getcoll4s_3d, E_exc, Q, /PLOT

  i = 0  &  j = 1 
  FOR n=0,Ntemp-1 DO $
   ohm(n) = maxwell_avg(g(i), T(n), E_edge, 0.3, E_exc, Q)
  print, FORMAT=fmt, j, i, 0.4*ohm
  i = 0  &  j = 2
  print, FORMAT=fmt, j, i, 0.6*ohm

  E_edge = (E(2) - E(1))* 1.0E2 * (hPlanck * cLight) / ERydb
  i = 1  &  j = 2
  FOR n=0,Ntemp-1 DO $
   ohm(n) = maxwell_avg(g(i), T(n), E_edge, 0.0, E_exc, Q)
  print, FORMAT=fmt, j, i, ohm

  E_edge = (g(3)*E(3) + g(4)*E(4))/(g(3) + g(4)) * 1.0E2 * $
   (hPlanck * cLight) / ERydb
  getcoll4s_4p, E_exc, Q, /PLOT

  i = 0  &  j = 3 
  FOR n=0,Ntemp-1 DO $
   ohm(n) = maxwell_avg(g(i), T(n), E_edge, 0.5, E_exc, Q)
  print, FORMAT=fmt, j, i, 0.5*ohm
  i = 0  &  j = 4
  print, FORMAT=fmt, j, i, ohm

  E_edge = (E(4) - E(3))* 1.0E2 * (hPlanck * cLight) / ERydb
  i = 3  &  j = 4
  FOR n=0,Ntemp-1 DO $
   ohm(n) = maxwell_avg(g(i), T(n), E_edge, 0.0, E_exc, Q)
  print, FORMAT=fmt, j, i, ohm


  E_edge = ( (g(3)*E(3) + g(4)*E(4))/(g(3) + g(4)) - $
             (g(1)*E(1) + g(2)*E(2))/(g(1) + g(2)) ) * $
             1.0E2 * (hPlanck * cLight) / ERydb
  getcoll3d_4p, E_exc, Q, /PLOT

  i = 1  & j = 3
  FOR n=0,Ntemp-1 DO $
   ohm(n) = maxwell_avg(g(i), T(n), E_edge, 0.0, E_exc, Q(*,0))
  print, FORMAT=fmt, j, i, ohm
  i = 1  & j = 4
  FOR n=0,Ntemp-1 DO $
   ohm(n) = maxwell_avg(g(i), T(n), E_edge, 0.0, E_exc, Q(*,1))
  print, FORMAT=fmt, j, i, ohm
  i = 2  & j = 3
  FOR n=0,Ntemp-1 DO $
   ohm(n) = maxwell_avg(g(i), T(n), E_edge, 0.0, E_exc, Q(*,2))
  print, FORMAT=fmt, j, i, ohm
  i = 2  & j = 4
  FOR n=0,Ntemp-1 DO $
   ohm(n) = maxwell_avg(g(i), T(n), E_edge, 0.0, E_exc, Q(*,3))
  print, FORMAT=fmt, j, i, ohm

  fmt2 = '(" CI", 5X, I2, X, I2, 2X, 6E11.3)'
  alpha = [2.0363E-23, 6.1484E-22, 6.1484E-22, 2.3823E-22, 2.3823E-22]
  C2 = 1.55E+11     ; s^-1 m K^1/2
  j  = 5
  FOR i=0,4 DO BEGIN
    CI = 0.2 * C2 * alpha(i) * (kBoltz/(EJoule(j) - EJoule(i)))
    print, FORMAT=fmt2, i, j, CI + fltarr(Ntemp)
  ENDFOR
END
