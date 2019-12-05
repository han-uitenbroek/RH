FUNCTION vacuumtoair, wave, TO_AIR_LIMIT=to_air_limit
;+
; NAME:
;	VACUUMTOAIR
; PURPOSE:
;	Convert vacuum wavelengths to air wavelengths, i.e. correct
;	for the index of refraction of air under standard conditions.  
;	Wavelength values below 200 nm will not be altered.  Accurate to 
;	about 0.005 A 
;
; CALLING SEQUENCE:
;	w_air = VACUUMTOAIR(W [,/TO_AIR_LIMIT])
;
; INPUT/OUTPUT:
;	WAVE - Wavelength in Angstroms, scalar or vector
;		WAVE should be input as vacuum wavelength(s), it will be
;		returned as air wavelength(s).  WAVE is always converted to
;		double precision
;
; EXAMPLE:
;	If the vacuum wavelength is  W = 200.0, then 
;
;	IDL> w_air = VACUUMTOAIR(W [,/TO_AIR_LIMIT])
;
;	yields an air wavelength of W = 199.9353 nm
;
; METHOD:
;	An approximation to the 4th power of inverse wavenumber is used
;	See IUE Image Processing Manual   Page 6-15.
;
; REVISION HISTORY
;	Written, D. Lindler 1982 
;	Documentation W. Landsman  Feb. 1989
;       Modified, Han Uitenbroek, Jul 23 1996
;-
  On_error, 2
  IF (N_params() EQ 0) THEN BEGIN
    print,'Syntax - w_air = VACUUMTOAIR(Wave [, /TO_AIR_LIMIT])'
    return, 0.0
  ENDIF

  IF (keyword_set(TO_AIR_LIMIT)) THEN $
   limit = to_air_limit $
  ELSE $
   limit = 200.0

  wave2 = double(wave)*wave 
  fact = 1.0 + 2.735182D-4 + (1.314182D0 + 2.76249D+4/wave2) / wave2
  fact = fact * (wave GE limit) + 1.0 * (wave LT limit)

; Convert wavelengths

  return, float(wave/fact)
END                        
