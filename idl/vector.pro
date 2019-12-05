;-------------------------------------------------------------
;+
; NAME:
;       VECTOR
; PURPOSE:
;       Force given argument to be an vector.
; CATEGORY:
; CALLING SEQUENCE:
;       y = vector(x)
; INPUTS:
;       x = input which may be an vector or scalar.      in
; KEYWORD PARAMETERS:
; OUTPUTS:
;       y = out which is an vector.                      out
; COMMON BLOCKS:
; NOTES:
; MODIFICATION HISTORY:
;       Han Uitenbroek, Oct 92
;-
;-------------------------------------------------------------
 
	function VECTOR, x
 
	if (n_params(0) lt 1) THEN BEGIN
	  print,' Force given argument to be a vector.'
	  print,' y = vector(x)'
	  print,'   x = input which may be an vector or scalar.      in'
	  print,'   y = out which is an vector.                      out'
	  return, -1
	ENDIF
 
	IF n_elements(x) eq 0 THEN BEGIN
	  print,' Error in ARRAY: argument undefined.'
	  stop, ' Stopping in ARRAY.'
	  return, -1
	ENDIF

	s = size(x)
	if s(0) gt 0 then return, x	; already a vector.
	n = s(s(0)+2)			; number of elements.
	type = s(s(0)+1) 
 
	CASE type OF
          1: y = bytarr(n) + x
          2: y = intarr(n) + x
          3: y = lonarr(n) + x
          4: y = fltarr(n) + x
          5: y = dblarr(n) + x
          6: y = complexarr(n) + x
	ENDCASE
 
	return, y
	end
