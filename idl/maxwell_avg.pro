FUNCTION maxwell_avg, gi, T, E0, Elin, E, Q

  ;; --- Compute Maxwellian average of the collisional cross sections
  ;;     Q (given in units of pi*a_0^2) for impact energies E given
  ;;     in Rydberg. --

  COMMON constant, cLight,  hPlanck, kBoltz, rBohr, ERydb, mElect, $
   one_eV, pia02

  ;; --- Use 8-point Gauss-Laguerre quadrature --

  uLaguerre = [ 0.170279632305D+00,  0.903701776799D+00, $
                2.251086629866D+00,  4.266700170288D+00, $
                7.045905402393D+00, 10.758516010181D+00, $
               15.740678641278D+00, 22.863131736889D+00 ]

  wLaguerre = [3.69188589342D-01, 4.18786780814D-01, $
               1.75794986637D-01, 3.33434922612D-02, $
               2.79453623523D-03, 9.07650877336D-05, $
               8.48574671627D-07, 1.04800117487D-09 ]

  kT = kBoltz * T
  ulin = Elin*ERydb / kT  &  u = E*ERydb / kT  &  u0 = E0*ERydb / kT

  linearPart = where(u LE ulin, count)
  IF (count GT 1) THEN BEGIN
   funct = u*Q(linearPart)*exp(u0 - u(linearPart))
   OhmegaLin = int_tabulated([u0, u(linearPart)], [u0*Q(0), funct])
   uc = u(count)
  ENDIF ELSE BEGIN
   uc = u0
   OhmegaLin = 0.0
  ENDELSE

  huspline, u, Q, uc + uLaguerre, QLaguerre
  OhmegaLag = exp(u0-uc) * total((uc + uLaguerre) * $
                                       QLaguerre * wLaguerre)

  print, 'Ohm_lin = ', OhmegaLin, ',  Ohm_Lag = ', OhmegaLag
  return, gi * kT/ERydb * (OhmegaLin + OhmegaLag)
END

