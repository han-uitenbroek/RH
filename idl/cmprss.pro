function CMPRSS, array
;+
; NAME:
;       CMPRSS
; PURPOSE:
;       Remove trivial dimensions from input array
; CATEGORY:
; CALLING SEQUENCE:
;       cmparr = CMPRSS(array)
; INPUTS:
;       array  --  Array to be compressed
; KEYWORD PARAMETERS:
; OUTPUTS:
;       cmparr --  Array with trivial dimensions removed
; COMMON BLOCKS:
; NOTES:
;       
; MODIFICATION HISTORY:
;       Han Uitenbroek
;-
  dims = size(array)
  ndim = dims[0]
  dims = dims[1:ndim]
  dims = dims[where(dims gt 1)]

  return, reform(array, dims, /OVERWRITE)
END
