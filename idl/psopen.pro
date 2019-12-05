PRO psopen, FILENAME=filename, $
            LANDSCAPE=landscape, $
            COLOR=color, $
            FONT=font, $
            ENCAPSULATED=encapsulated, $
            BITS_PER_PIXEL=bits_per_pixel, $
            ADOBEDEFAULT=AdobeDefault, $
            XSIZE=xsize, YSIZE=ysize, XOFFSET=xoffset, YOFFSET=yoffset
  
@ps.common

  P_old = !P
  OldDevice = !D.NAME
  ColorPrint = 0

  set_plot, 'PS'
  IF (keyword_set(FONT)) THEN !P.FONT = font ELSE !P.FONT = 0
  IF (keyword_set(LANDSCAPE)) THEN $
   device, /LANDSCAPE $
  ELSE $
   device, /PORTRAIT

  IF (keyword_set(ENCAPSULATED)) THEN device, /ENCAPSULATED
  IF (NOT keyword_set(FILENAME)) THEN $
   filename = '/tmp/idl-' + getenv('USER') + $
   ((keyword_set(ENCAPSULATED)) ? '.eps' : '.ps')
  PS_FileName = filename

  IF (keyword_set(COLOR)) THEN BEGIN
    device, /COLOR
    ColorPrint = 1
  ENDIF
  IF (keyword_set(BITS_PER_PIXEL))       THEN $
   device, BITS_PER_PIXEL=bits_per_pixel $
  ELSE $
   device, BITS_PER_PIXEL=24

  IF (NOT keyword_set(ENCAPSULATED)) THEN BEGIN
    IF (NOT keyword_set(XSIZE)) THEN $
     IF (keyword_set(LANDSCAPE)) THEN  xsize=9.0  ELSE  xsize=7.0
    IF (NOT keyword_set(YSIZE)) THEN $
     IF (keyword_set(LANDSCAPE)) THEN  ysize=7.0  ELSE  ysize=5.0
    IF (NOT keyword_set(XOFFSET)) THEN  xoffset=0.75
    IF (NOT keyword_set(YOFFSET)) THEN $
     IF (keyword_set(LANDSCAPE)) THEN  yoffset=10.0  ELSE  yoffset=1.0

    device, XSIZE=xsize, YSIZE=ysize, XOFFSET=xoffset, YOFFSET=yoffset, /INCHES
  ENDIF

  device, FILENAME=FileName

  IF (!P.FONT EQ 0) THEN BEGIN
    device, /HELVETICA, FONT_INDEX=3, /NARROW
    device, /TIMES,     FONT_INDEX=5
    device, /SYMBOL,    FONT_INDEX=7
  ENDIF
  IF (NOT keyword_set(ADOBEDEFAULT)) THEN $
   device, /ISOLATIN1 $
  ELSE $
   device, ISOLATIN1=0

  return
END
