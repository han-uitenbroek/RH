PRO vv, u, v, x, y, LENGTH=length, V_COLORS=v_colors, _EXTRA=extra

;; Adaptation of IDL's velovect to get the velocity vectors at the
;; proper x and y location without extending the x- and y axes.
;;
;; Han Uitenbroek,  Feb 3, 1997

  IF (n_params(0) LT 4) THEN BEGIN
    print, "Usage: vv, u, v, x, y, LENGTH=length, _EXTRA=extra"
    return
  ENDIF
  IF (NOT keyword_set(LENGTH)) THEN length = 1.0

  s = size(u)  &  sx = size(x)
  t = size(v)  &  sy = size(y)
  IF (s(0) NE 2  OR  t(0) NE 2  OR  sx(0) NE 1  OR  sy(0) NE 1) THEN BEGIN
    print, "Error: u and v have to be 2-dimensionale arrays,"
    print, "       and x and y have to be 1-dimensionale arrays"
    help, u, v, x, y
    return
  ENDIF ELSE BEGIN
    IF (s(1) NE t(1)  OR  s(2) NE t(2)) THEN BEGIN
      print, "Error: u and v should have the same dimensions"
      help, u, v
      return
    ENDIF
    IF (s(1) NE n_elements(x)) THEN BEGIN
      print, $
      "Error: first dimension of u and v should be the same as that of x"
      help, u, v, x
      return
    ENDIF
    IF (s(2) NE n_elements(y)) THEN begin
      print, $
      "Error: second dimension of u and v should be the same as that of y"
      help, u, v, x
      return
    ENDIF
  ENDELSE

  magnitude = sqrt(u^2.+v^2.)

  x0 = min(x)  &  x1 = max(x)
  y0 = min(y)  &  y1 = max(y)

  x_step = float(x1-x0) / float(s(1)-1)
  y_step = float(y1-y0) / float(s(2)-1)

  max_mag = max([max(abs(u/x_step)), max(abs(v/y_step))])
  sina = length * (u/max_mag)
  cosa = length * (v/max_mag)

  plot, [x0, x1], [y0, y1], /NODATA, _EXTRA = extra

  r = 0.3
  angle = 22.5 * !dtor
  st = r * sin(angle)
  ct = r * cos(angle)

  thick = !P.THICK  &  thecolor = !P.COLOR
  clip  = !P.CLIP   &  t3d   = !P.T3D
  IF (n_elements(extra) GT 0) THEN BEGIN
    names = tag_names(extra)
    FOR n=0, n_tags(extra)-1 DO BEGIN
      IF (strpos("THICK", names(n)) EQ 0) THEN thick = extra.(n)
      IF (strpos("COLOR", names(n)) EQ 0) THEN thecolor = extra.(n)
      IF (strpos("CLIP", names(n)) EQ 0)  THEN clip  = extra.(n)
      IF (strpos("T3D", names(n)) EQ 0)   THEN t3d   = extra.(n)
    ENDFOR
  ENDIF
  FOR i=0, s(1)-1 DO BEGIN
    x0 = x(i)
    FOR j=0, s(2)-1 DO BEGIN
      dx = sina(i, j)
      x1 = x0 + dx
      y0 = y(j)
      dy = cosa(i, j)
      y1 = y0 + dy
      xd = x_step
      yd = y_step

      IF (keyword_set(V_COLORS)) THEN $
       color = v_colors(i, j) $
      ELSE $
       color = thecolor

      plots, [x0, x1, x1-(ct*dx/xd-st*dy/yd)*xd, $
              x1, x1-(ct*dx/xd+st*dy/yd)*xd], $
             [y0, y1, y1-(ct*dy/yd+st*dx/xd)*yd, $
              y1, y1-(ct*dy/yd-st*dx/xd)*yd], $
       THICK=thick, COLOR=color, CLIP=clip, T3D=t3d
    ENDFOR
  ENDFOR
END
