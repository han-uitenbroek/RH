FUNCTION airtovacuum, wave, TO_VACUUM_LIMIT=to_vacuum_limit                   
;+
; NAME:
;	AIRTOVACUUM
; PURPOSE:
;	Convert air wavelengths to vacuum wavelengths, i.e. correct 
;	for the index of refraction of air under standard conditions.  
;	Wavelength values below 200.0 nm will not be altered.  Uses the IAU 
;	standard for conversion given in Morton (1991 Ap.J. Suppl. 77, 119)
;
; CALLING SEQUENCE:
;	W_VAC = AIRTOVACUUM(WAVE [, /TO_VACUUM_LIMIT)
;
; INPUT/OUTPUT:
;	WAVE - Wavelength in Angstroms, scalar or vector
;		WAVE should be input as air wavelength(s), it will be
;		returned as vacuum wavelength(s).
;
; EXAMPLE:
;	If the air wavelength is W = 605.6125 (a Krypton line), 
;	then AIRTOVACUUM yields a vacuum wavelength of W = 605.78019
;
; METHOD:
;	See Morton (Ap. J. Suppl. 77, 119) for the formula used
;
; REVISION HISTORY
;	Written W. Landsman                November 1991
;       Revised H. Uitenbroek, Jul 23 1996
;-
  On_error, 2

  IF N_params() EQ 0 THEN BEGIN
    print,'Syntax - w_vac = AIRTOVACUUM(WAVE [, /TO_VACUUM_LIMIT)'
    print,'WAVE (Input) is the air wavelength in nm'
    print,'On output WAVE contains the vacuum wavelength in nm'
    return, 0.0
  ENDIF

  IF (keyword_set(TO_VACUUM_LIMIT)) THEN $
   limit = to_vacuum_limit $
  ELSE $
   limit = 199.9352

  sigma2 = (1.0D+7 / wave)^2
  fact = 1.0000834213D+0 + $
   2.406030D+6/(1.30D+10 - sigma2) + 1.5997D+4/(3.89D+9 - sigma2)
  fact = fact * (wave GE limit) + 1.0 * (wave LT limit)
  
  return, wave * fact
END
