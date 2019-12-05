FUNCTION read1datmos, fileName, NHYDR=NHydr

;+
; NAME:
;	READ1DATMOS
;
; PURPOSE:
;	This routine reads 1-dimensional model atmospheres in MULTI format.
;       This is also the format used by RHF1D.
;
; CATEGORY:
;	Data reduction
;
; CALLING SEQUENCE:
;       atmos = READ1DATMOS(fileName)
;
; INPUTS:
;	fileName:  Name of the file to be read.
;
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;	atmos:  Model atmosphere in atmos structure.
;                (see routine READATMOS).
;
; OPTIONAL OUTPUTS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; EXAMPLE:
; MODIFICATION HISTORY:
;
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Fri Oct 11 16:18:01 2013 --
;-

  CM_TO_M = 1.0D-02
  KM_TO_M = 1.0D+03
  G_TO_KG = 1.0E-3
  

  openr, atmosUnit, fileName, /GET_LUN

  atmosID = ''
  inputLine = getline(atmosUnit, COMMENT_CHAR='*', /EXIT_ON_EOF)
  reads, inputLine, FORMAT='(A)', atmosID

  massScale = ''
  inputLine = getline(atmosUnit, COMMENT_CHAR='*', /EXIT_ON_EOF)
  reads, inputLine, FORMAT='(A)', massScale
  type = strupcase(strmid(massScale, 0, 1))
  CASE (type) OF
    'M': scaleType = 'MASS_SCALE'
    'T': scaleType = 'TAU500_SCALE'
    'H': scaleType = 'GEOMETRIC_SCALE'
    ELSE:
  ENDCASE

  gravitation = 0.0D+0
  inputLine = getline(atmosUnit, COMMENT_CHAR='*', /EXIT_ON_EOF)
  reads, inputLine, FORMAT='(E)', gravitation

  Ndep = 0L
  inputLine = getline(atmosUnit, COMMENT_CHAR='*', /EXIT_ON_EOF)
  reads, inputLine, FORMAT='(I)', Ndep

  depth = dblarr(Ndep)
  T = depth  &  n_elec = depth  &  v = depth  &  vturb = depth
  dummy = dblarr(5)
  FOR n=0, Ndep-1 DO BEGIN
    inputLine = getline(atmosUnit, COMMENT_CHAR='*', /EXIT_ON_EOF)
    reads, inputLine, dummy

    depth(n)  = dummy(0)
    T(n)      = dummy(1)
    n_elec(n) = dummy(2)
    v(n)      = dummy(3)
    vturb(n)  = dummy(4)
  ENDFOR

  inputLine = getline(atmosUnit, COMMENT_CHAR='*')
  s = size(inputLine)

  IF (s(1) NE 7) THEN BEGIN
    HLTE = 1
    nH = 0.0D+0
  ENDIF ELSE BEGIN
    HLTE = 0

    if (not keyword_set(Nhydr)) then Nhydr = 6

    nH = dblarr(Ndep, Nhydr)
    dummy = dblarr(Nhydr)
    reads, inputLine, dummy
    nH(0, *) = dummy
    FOR n=1, Ndep-1 DO BEGIN
      inputLine = getline(atmosUnit, COMMENT_CHAR='*', /EXIT_ON_EOF)
      reads, inputLine, dummy
      nH(n, *) = dummy
    ENDFOR
  ENDELSE
  free_lun, atmosUnit

  ;; Convert to proper units

  nH /= CM_TO_M^3
  n_elec /= CM_TO_M^3

  v     *= KM_TO_M
  vturb *= KM_TO_M

  CASE (scaleType) OF
    'MASS_SCALE': BEGIN
      depth = (10.0^depth) * G_TO_KG / CM_TO_M^2

      atmos = {atmosID: atmosID,  gravitation: 10.0^gravitation, $
               Ndep: long(Ndep),  scale: scaleType,  cmass: depth, $
               T: T,  n_elec: n_elec,  v: v,  vturb: vturb, $
               nH: nH}
    END
    'TAU500_SCALE': BEGIN
      atmos = {atmosID: atmosID,  gravitation: 10.0^gravitation, $
               Ndep: long(Ndep),  scale: scaleType,  tau500: 10.0^depth, $
               T: T,  n_elec: n_elec,  v: v,  vturb: vturb, $
               nH: nH}
    END
    'GEOMETRIC_SCALE': BEGIN
      atmos = {atmosID: atmosID,  gravitation: 10.0^gravitation, $
               Ndep: long(Ndep),  scale: scaleType,  height: depth, $
               T: T,  n_elec: n_elec,  v: v,  vturb: vturb, $
               nH: nH}
    END
  ENDCASE

  return, atmos
END
