FUNCTION cc_kernel, s
 
  s2   = s * s
  s3   = s2 * s
  u = [-0.5*(s3 + s) + s2,       1.5*s3 - 2.5*s2 + 1.0, $
       -1.5*s3 + 2.0*s2 + 0.5*s, 0.5*(s3 - s2) ]
 
  return, u
END
FUNCTION MOD_B_spline, x, y, xx, DY=dy, NY=ny, MONOTONIC=monotonic, $
                       CUBIC_CONVOLUTION=cubic_convolution

  Nx = n_elements(x)
  ysp1 = fltarr(5)  &  ysp1[2] = 2.0/3.0
  ysp2 = fltarr(3)  &  ysp2[1] = 1.0  

  spl_str = {a: fltarr(4),  b: fltarr(4), c: fltarr(4),  d: fltarr(4)}
  spl = replicate(spl_str, Nx)
  h = x - shift(x, 1)

  FOR n=2, Nx-3 DO BEGIN
    M = huspline(x[n-2:n+2], ysp1, YP0=0.0, YPN=0.0)
    spl[n].a[[0, 3]] = ysp1[[0, 3]]
    spl[n].b[[0, 3]] = -(2*M[[0, 3]] + M[[1, 4]]) * h[[n-1, n+2]] / 6
    spl[n].c[[0, 3]] = M[[0, 3]] / 2
    spl[n].d[[0, 3]] = (M[[1, 4]] - M[[0, 3]]) / (6 * h[[n-1, n+2]])

    xx1 = x[n-2] + (x[n+2] - x[n-2]) * findgen(100)/99
    plot, xx1, splint(x[n-2:n+2], ysp1, m, xx1)

    yp0 =  (M[1]/3 + M[0]/6) * h[n-1]
    ypN = -(M[3]/3 + M[4]/6) * h[n+2]
    M2 = huspline(x[n-1:n+1], ysp2, YP0=yp0, YPN=ypN)
    spl[n].a[1:2] = ysp2[0:1]
    spl[n].b[1:2] = (ysp2[1:2] - ysp2[0:1]) / h[n:n+1] - $
     (2*M2[0:1] + M2[1:2]) * h[n:n+1] / 6
    spl[n].c[1:2] = M2[0:1] / 2
    spl[n].d[1:2] = (M2[1:2] - M2[0:1]) / (6 * h[n:n+1])

    xx2 = x[n-1] + (x[n+1] - x[n-1]) * findgen(100)/99
    oplot, xx2, splint(x[n-1:n+1], ysp2, m2, xx2), color=15
    wait, 2
  ENDFOR

  Nxx = n_elements(xx)
  tabinv, x, xx, i_eff
  yy = fltarr(Nxx)
  dy = fltarr(Nxx)
  ny = fltarr(Nxx)
  FOR i=0, Nxx-1 DO BEGIN
    n = long(i_eff[i])
    IF (keyword_set(CUBIC_CONVOLUTION)) THEN BEGIN
      rx = i_eff[i] - n
      yy[i] = total(y[n-1:n+2] * cc_kernel(rx))
    ENDIF ELSE BEGIN
      rx = xx[i] - x[n]
      FOR j=-1, 2 DO BEGIN
        s = spl[n+j]
        u = s.a[2-j] + rx*(s.b[2-j] + rx*(s.c[2-j] + rx*s.d[2-j]))
        yy[i] = yy[i] + y[n+j] * u
        ny[i] = ny[i] + u
        u = s.b[2-j] + rx*(2*s.c[2-j] + rx*3*s.d[2-j])
        dy[i] = dy[i] + y[n+j] * u
      ENDFOR
    ENDELSE

;    IF (keyword_set(MONOTONIC) AND $
;        NOT keyword_set(CUBIC_CONVOLUTION)) THEN BEGIN
;      IF (D2[i] GT 0.0) THEN BEGIN
;        print, "MONOTONIC active for xx = ", xx[i], ", positive second"
;        IF (yy[i] LT y[n] AND yy[i] LT y[n+1]) THEN yy[i] = min(y[n:n+1])
;      ENDIF ELSE BEGIN
;        print, "MONOTONIC active for xx = ", xx[i], ", negative second"
;        IF (yy[i] GT y[n] AND yy[i] GT y[n+1]) THEN yy[i] = max(y[n:n+1])
;      ENDELSE
;    ENDIF   
  ENDFOR

  return, yy
END
