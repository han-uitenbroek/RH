FUNCTION monoint, x, y, xp

  i = where(x GT xp, count)
  IF (count LT 1) THEN BEGIN
    print, "Interpolation point outside domain"
    return, 0.0
  ENDIF

  Nx = n_elements(x)
  IF (i[0] GT 1  AND  i[0] LT Nx-1) THEN BEGIN
    i1 = i[0]
    i0 = i1 - 1
    h1 = x[i1] - x[i0]
    h0 = x[i0] - x[i0-1]
    s1 = (y[i1] - y[i0]) / h1
    s0 = (y[i0] - y[i0-1]) / h0

    yp = (s0*h1 + s1*h0) / (h0 + h1)
    IF (s0*s1 LE 0.0) THEN $
     yp = 0.0 $
    ELSE BEGIN
      IF (abs(yp) GT 2*abs(s0) OR abs(yp) GT 2*abs(s1)) THEN $
       yp = 2 * ((s0 GT 0.0) ? 1.0 : -1.0) * min(abs([s0, s1]))
    ENDELSE

    ai = (s1 - yp) / h1
    bi = yp
    ti = (xp - x[i0])

    ymono = y[i0] + (bi + ai*ti) * ti
  ENDIF ELSE BEGIN
    IF (i[0] EQ 1) THEN BEGIN
      ymono = y[0]
    ENDIF ELSE BEGIN
      ymono = y[Nx-1]
    ENDELSE
  ENDELSE

  return, ymono
END

x = [0.0, 1.0, 2.1, 2.9, 3.1, 4.0, 5.2, 6.0, 7.1, 9.0, 10.0]
y = [10.0, 9.5, 9.1, 8.3, 3.0, 2.0, 1.6, 1.2, 1.0, 1.05, 0.9]

xx=findgen(100)/9.9
yy=fltarr(100)
for i=0,99 do yy[i]=monoint(x, y, xx[i])
