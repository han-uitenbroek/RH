FUNCTION angstrom

  ;; Define the Angstrom symbol (for people who do NOT want to use nm).

  IF (!D.NAME EQ 'PS') THEN BEGIN
    help, /DEVICE, OUTPUT=help_text
    index = where(strpos(help_text, 'ISOLatin1', 0) ge 0, count)
    IF (count EQ 0) THEN BEGIN
      message, "Setting font encoding to ISOLatin1", /INFORMATIONAL
      device, /ISOLATIN1
    ENDIF
  ENDIF

  CASE (!P.FONT) OF
  -1: BEGIN                               ;; System fonts
    message, "Cannot print Angstrom in system font", /CONTINUE
    return, "\AA"
  END
  0: return, string("305B)                ;; PS fonts
  1: return, string("305B)                ;; Hershey fonts
  ENDCASE

END
