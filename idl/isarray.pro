;-------------------------------------------------------------
;+
; NAME:
;       ISARRAY
; PURPOSE:
;       Tests if the argument is an array.
; CATEGORY:
; CALLING SEQUENCE:
;       flag = isarray(a)
; INPUTS:
;       a = variable to test.                                in
; KEYWORD PARAMETERS:
; OUTPUTS:
;       flag = test result: 0 if not array, else non-zero.   out
; COMMON BLOCKS:
; NOTES:
; MODIFICATION HISTORY:
;       R. Sterner  20 Mar, 1986.
;       Checked for undefined variables.  RES 25 Aug, 1989.
;       Johns Hopkins Applied Physics Lab.
;
; Copyright (C) 1986, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
 
	FUNCTION ISARRAY,A, help = h
 
	if (n_params(0) lt 1) or keyword_set(h) then begin
	  print,' Tests if the argument is an array.' 
	  print,' flag = isarray(a)' 
	  print,'   a = variable to test.                                in'
	  print,'   flag = test result: 0 if not array, else non-zero.   out' 
	  return, -1
	endif
 
	if n_elements(a) eq 0 then return, 0
	S = SIZE(A)
	RETURN, S(0) NE 0
	END
