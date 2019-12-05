FUNCTION PLANCK, T, lambda, HZ=hz

  ;; --- lambda in [nm] --                              -------------- ;

  CLIGHT     = 2.99792458E+08
  HPLANCK    = 6.626176E-34
  KBOLTZMANN = 1.380662E-23

  CM_TO_M      = 1.0E-02
  NM_TO_M      = 1.0E-09
  ERG_TO_JOULE = 1.0E-07

  ;; --- Convert lambda to [m] --                       -------------- ;

  lambda_m   = NM_TO_M * lambda

  hc_kl      = double((HPLANCK * CLIGHT) / (KBOLTZMANN * lambda_m))
  twohnu3_c2 = (2.0*HPLANCK*CLIGHT) / lambda_m^3

  IF (keyword_set(HZ)) THEN BEGIN

    ;; --- In [J m^-2 s^-1 Hz^-1 sr^-1] --             --------------- ;

    return, twohnu3_c2 / (exp(hc_kl/T) - 1.0)
  ENDIF ELSE BEGIN

    ;; --- In [erg cm^-2 s^-1 nm^-1 sr^-1] --          --------------- ;

    C = NM_TO_M * CM_TO_M^2 / ERG_TO_JOULE

    return, C * twohnu3_c2 * (CLIGHT/lambda_m^2) / (exp(hc_kl/T) - 1.0)
  ENDELSE
END
