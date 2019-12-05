FUNCTION cc_kernel, s

  s2   = s * s
  s3   = s2 * s
  u = [-0.5*(s3 + s) + s2,       1.5*s3 - 2.5*s2 + 1.0, $
       -1.5*s3 + 2.0*s2 + 0.5*s, 0.5*(s3 - s2) ]

  return, u
END

FUNCTION modulo, n, Np

  m = n 
  index = where(n GE 0, count)
  IF (count GT 0) THEN m(index) = n(index) MOD Np

  index = where(n LT 0, count)
  IF (count GT 0) THEN m(index) = Np - (abs(n(index)) MOD Np)

  return, m
END

FUNCTION cubeshift_1D, f, x, INVERSE=inverse

  IF (keyword_set(INVERSE)) THEN x = -x

  N = n_elements(f)
  dn = fix(x)
  xoffset = x - fix(dn)
  IF (xoffset LT 0.0) THEN BEGIN
    dn = dn - 1
    xoffset = 1.0 + xoffset
  ENDIF
  kernel = cc_kernel(xoffset)

  g = fltarr(N)
  FOR i=0, N-1 DO BEGIN
    index = modulo(i + [-1, 0, 1, 2] + dn, N-1)
    fsub = f(index)
    fsubmin = min(fsub, MAX=fsubmax)
    fxx  = (fsub(3) - fsub(2)) - (fsub(1) - fsub(0))
    g(i) = total(f(index) * kernel)
    IF (fxx GT 0) THEN g(i) = max([g(i), fsubmin])
    IF (fxx LT 0) THEN g(i) = min([g(i), fsubmax])
  ENDFOR

  return, g
END

FUNCTION cubeshift_2D, A, x, y, INVERSE=inverse, MONOTONIC=monotonic

  IF (keyword_set(INVERSE)) THEN BEGIN
    x = -x
    y = -y
  ENDIF

  s = size(A)
  Nx = s(1)  &  Ny = s(2)

  dlx = fix(x)
  xoffset = x - fix(dlx)
  IF (xoffset LT 0.0) THEN BEGIN
    dlx = dlx - 1
    xoffset = 1.0 + xoffset
  ENDIF
  xkernel = cc_kernel(xoffset)

  dly = fix(y)
  yoffset = y - fix(dly)
  IF (yoffset LT 0.0) THEN BEGIN
    dly = dly - 1
    yoffset = 1.0 + yoffset
  ENDIF
  ykernel = cc_kernel(yoffset)

  g = fltarr(Nx, Ny)
  FOR ly=0, Ny-1 DO BEGIN
    yindex = modulo(ly + [-1, 0, 1, 2] + dly, Ny)
    FOR lx=0, Nx-1 DO BEGIN
      xindex = modulo(lx + [-1, 0, 1, 2] + dlx, Nx)

      subA = fltarr(4, 4)
      FOR m=0, 3 DO subA(*, m) = A(xindex, yindex(m))

      Amin = min(subA(1:2, 1:2), MAX=Amax)
      sum  = total((xkernel#ykernel) * subA)

      IF (keyword_set(MONOTONIC)) THEN BEGIN
        IF (xoffset GE 0.5  AND  NOT keyword_set(inverse)) THEN $
         IF (yoffset GE 0.5  AND  NOT keyword_set(inverse)) THEN $
          Apref = subA(2, 2) ELSE Apref = subA(2, 1) $
        ELSE $
         IF (yoffset GE 0.5  AND  NOT keyword_set(inverse)) THEN $
          Apref = subA(1, 2) ELSE Apref = subA(1, 1)

        IF (sum LE 0.0) THEN sum = Amin
      ENDIF
      g(lx, ly) = sum
    ENDFOR
  ENDFOR

  return, g
END
