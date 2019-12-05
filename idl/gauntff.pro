FUNCTION gauntff, lambda, charge, T
 
  HPLANCK    = 6.6260755E-34  
  CLIGHT     = 2.99792458E+08
  E_RYDBERG  = 2.1798741E-18
  KBOLTZMANN = 1.380658E-23
  NM_TO_M    = 1.0E-9 

  x  = ((HPLANCK * CLIGHT)/(lambda * NM_TO_M)) / (E_RYDBERG * charge^2)
  x3 = x^0.33333333
  y  = (2.0 * lambda * NM_TO_M * KBOLTZMANN*T) / (HPLANCK*CLIGHT)

  gIII = (1.0 + 0.1728*x3 * (1.0 + y) - $
   0.0496*x3^2 * (1.0 + (1.0 + y)*0.33333333*y)) > 1.0

  return, gIII
END
