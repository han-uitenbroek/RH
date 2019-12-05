PRO rhannotate, xc, yc, $
                TEXT=text, $
                MARKTYPE=marktype, $
                NOMOUSE=nomouse, $
                CHARSIZE=charsize, $
                MARKCOLOR=markcolor, $
                CHARCOLOR=charcolor, $
                THICK=thick, $
                MARKLENGTH=marklength, $
                _EXTRA=extra_keywords
;+
; NAME:
;       RHANNOTATE
; PURPOSE:
;       Annotate plot
; CATEGORY:
; CALLING SEQUENCE:
;       ANNOTATE, xc, yc, TEXT=text, MARK=mark, /NOMOUSE, ....
; INPUTS:
;       XC  --  X-position of annotation in data coordinates
;       YC  --  Y-position 
; KEYWORD PARAMETERS:
;       TEXT       --  Text of annotation
;       MARKTYPE   --  Mark annotation, number is used as linestyle or
;                      symbol style (marktype<0).
;       NOMOUSE    --  Ask for XC and YC through keybord, not mouse
;       CHARSIZE   --  Scaling factor for character size
;       CHARCOLOR  --  Color of annotation
;       MARKCOLOR  --  Color of marks
;       THICK      --  Thickness of marks
;       MARKLENGTH --  Scaling length for marks
; OUTPUTS:
; COMMON BLOCKS:
; NOTES:
;       If device = PS, /NOMOUSE is in effect
; MODIFICATION HISTORY:
;       Han Uitenbroek, May 1991, last revised Mar 22 1996
;-

  IF (NOT keyword_set(MARKCOLOR))  THEN markcolor  = !P.COLOR
  IF (NOT keyword_set(CHARCOLOR))  THEN charcolor  = !P.COLOR
  IF (NOT keyword_set(THICK))      THEN thick      = 1
  IF (NOT keyword_set(MARKLENGTH)) THEN marklength = 1.0
  IF (NOT keyword_set(CHARSIZE))   THEN charsize   = 1.0

  IF (!D.NAME EQ 'X') AND (NOT keyword_set(NOMOUSE)) THEN BEGIN
    IF (n_params(0) LT 2) THEN BEGIN
      print, 'Use mouse to mark text position'
      cursor, xc, yc, /DOWN
      print, xc, yc
    ENDIF
  ENDIF ELSE BEGIN
    IF (n_params(0) LT 2) THEN $
     read, 'Give x and y position label (data coordinates):', xc, yc
  ENDELSE

  IF (NOT keyword_set(TEXT)) THEN BEGIN
    text= ''
    read, 'Give text for label: ', text
  ENDIF

  xt = xc
  IF (keyword_set(MARKTYPE) OR n_elements(MARKTYPE) GT 0) THEN BEGIN

    IF (!X.TYPE EQ 1) THEN BEGIN
      xm = xc * 10^(0.05*marklength * (!X.CRANGE(1) - !X.CRANGE(0)))
      xt = xc * 10^((0.05*marklength + 0.02) * (!X.CRANGE(1) - !X.CRANGE(0)))
    ENDIF ELSE BEGIN
      xm = xc + 0.05*marklength * (!X.CRANGE(1) - !X.CRANGE(0))
      xt = xc + (0.05*marklength + 0.02) * (!X.CRANGE(1) - !X.CRANGE(0))
    ENDELSE

    IF (marktype LT 0) THEN $
     oplot, [xc, 0.5*(xc + xm),  xm], [yc, yc, yc], $
     PSYM=-marktype, COLOR=markcolor $
    ELSE $
     oplot, [xc, xm], [yc, yc], LINESTYLE=marktype, COLOR=markcolor, $
     THICK=thick
  ENDIF

  xyouts, xt, yc, text, /NOCLIP, CHARSIZE=charsize, COLOR=charcolor, $
   _EXTRA=extra_keywords
END
