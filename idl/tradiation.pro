FUNCTION tradiation, intensity, lambda0

;+
; NAME:
;	TRADIATION
;
; PURPOSE:
;	This function returns the radiation temperature of corresponding
;       to the intensities INTENSITY at central wavelength LAMBDA0.
;
; CATEGORY:
;	Radiative transfer
;
; CALLING SEQUENCE:
;	Result = TRADIATION(intensity, lambda0)
;
; INPUTS:
;	Intensity:  intensity [J m^-2 s^-1 Hz^-1 sr^-1]
;       Lambda0:    central wavelength [nm]
;
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	Returns radiation temperature in K.
;
; OPTIONAL OUTPUTS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; EXAMPLE:
;
; MODIFICATION HISTORY:
;
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Thu Apr  9 12:33:46 1998 --
;-

  IF (n_params(0) LT 2) THEN BEGIN
    print, "Usage: Trad = TRADIATION(intensity, lambda0)"
    return, 0.0
  ENDIF

  CLIGHT     = 2.99792458E+08
  HPLANCK    = 6.626176E-34
  KBOLTZMANN = 1.380662E-23
  NM_TO_M    = 1.0E-09

  IF (n_elements(lambda0) GT 1) THEN BEGIN
    lambda = avg(lambda0)
    print, FORMAT='("Using average wavelength: ", F8.2, " [nm])', lambda
    lambda = lambda * NM_TO_M
  ENDIF ELSE $ 
   lambda = lambda0 * NM_TO_M

  trad = (HPLANCK*CLIGHT) / ((lambda*KBOLTZMANN) * $
             alog(1.0 + (2.0*HPLANCK*CLIGHT)/(lambda^3 * intensity)))

  return, trad
END

