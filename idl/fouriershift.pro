FUNCTION cwindow, relr

  c = [0.074, 0.302, 0.233, 0.390]

  IF (relr GT 1.5) THEN return, 0.0

  r  = 1.0 - relr^2
  r2 = r^2
  r3 = r2 * r
  cwind = ((c(0) + (r * c(1))) + (r2 * c(2))) + (r3 * c(3))

  return, (cwind > 0.0)
END

FUNCTION rolloff, Nx, Ny, cutoff

  Nx2 = Nx/2
  Ny2 = Ny/2
  r0  = sqrt(Nx^2 + Ny^2) * cutoff;

  filter = fltarr(Nx, Ny)

  FOR k=0, Ny2-1 DO BEGIN
    for l=0, Nx2-1 DO BEGIN
      d1 = sqrt(l^2 + k^2) / r0;
      d2 = sqrt((Nx2-l)^2 + k^2) / r0; 
      d3 = sqrt(l^2 + (Ny2-k)^2) / r0;
      d4 = sqrt((Nx2-l)^2 + (Ny2-k)^2) / r0;

      filter(l, k) = cwindow(d1)
      filter(Nx2+l, k) = cwindow(d2)
      filter(l, Ny2+k) = cwindow(d3)
      filter(Nx2+l, Ny2+k) = cwindow(d4)
    ENDFOR
  ENDFOR

  return, filter
end

FUNCTION fouriershift, f, dx, dy, LOWPASS_CUTOFF=lowpass_cutoff, $
  ATTENUATION=attenuation

  ;; Use Fourier transform to shift (periodic) one or two dimensional
  ;; Arrays. A shift in physical space corresponds to a multiplication
  ;; with a (frequency dependent) phase factor in Fourier space. Take
  ;; care that the phase shifted transform is Hermitian
  ;; (ie F_n = (F_N-n)^*) to ensure a real valued shifted function. The
  ;; lowest and highest frequncies are left unchanged.

  ;; Result can be either filtered with a low-pass filter to get rid of
  ;; high frequency noise introduced interpolation (set low-pass cutoff
  ;; with keyword LOWPASS_CUTOFF, or attenuation factors can be used
  ;; (see Stoer & Bulirsch, Introduction to Numerical Analysis, p. 72-93).
  ;; ATTENUATION can be set to either 1 or 2 for linear or cubic spline
  ;; behavior of the interpolating function.

  s = size(f)
  IF (s(0) LT 1  OR  s(0) GT 2  OR  n_params(0) NE s(0)+1) THEN BEGIN
    print, $
     "Usage: fs = FOURIERSHIFT(f, dx [, dy, LOWPASS_CUTOFF=lowpass_cutoff])"
    print, "  where f is either a one or two-dimensional array"
  ENDIF
  IF (NOT keyword_set(ATTENUATION)) THEN $
   attenuation = 0 $
  ELSE BEGIN
    IF (attenuation NE 1 AND attenuation NE 2) THEN BEGIN
      print, "Invalid value of keyword ATTENUATION: ", attenuation
      return, 0.0
    ENDIF
  ENDELSE

  Nx = s(1)
  IF (s(0) EQ 2) THEN Ny = s(2)

  ft = fft(f)

  zx = complexarr(Nx)
  zx(0) = complex(1.0, 0.0)
  FOR n=1, Nx/2 DO BEGIN
    alpha    = 2*!PI * (dx / Nx) * n
    zx(n)    = complex(cos(alpha), -sin(alpha))
    zx(Nx-n) = conj(zx(n))

    IF (attenuation GT 0) THEN BEGIN
      beta = (!PI*n) / Nx
      sinn = (sin(beta)/beta)^2
      CASE (attenuation) OF
        1: tau_n = sinn
        2: tau_n = (sinn)^2 * 3.0/(3.0 - 2.0*beta^2 * sinn)
      ENDCASE
      zx(n) = zx(n) * tau_n
      IF (n lt Nx/2) THEN zx(Nx-n) = zx(Nx-n) * tau_n
    END
  ENDFOR

  IF (s(0) EQ 2) THEN BEGIN
    zy = complexarr(Ny)
    zy(0) = complex(1.0, 0.0)
    FOR n=1, Ny/2 DO BEGIN
      alpha    = 2*!PI * (dy / Ny) * n
      zy(n)    = complex(cos(alpha), -sin(alpha))
      zy(Ny-n) = conj(zy(n))

      IF (attenuation GT 0) THEN BEGIN
        beta = (!PI*n) / Ny
        sinn = (sin(beta)/beta)^2
        CASE (attenuation) OF
          1: tau_n = sinn
          2: tau_n = (sinn)^2 * 3.0/(3.0 - 2.0*beta^2 * sinn)
        ENDCASE
        zy(n) = zy(n) * tau_n
        IF (n LT Ny/2) THEN zy(Ny-n) = zy(Ny-n) * tau_n
      END
    ENDFOR

    IF (keyword_set(LOWPASS_CUTOFF)) THEN BEGIN
     filter = rolloff(Nx, Ny, lowpass_cutoff)
     ftl = filter * ((zx#zy)*ft)
    ENDIF ELSE $
     ftl = (zx#zy)*ft

  ENDIF ELSE BEGIN
    IF (keyword_set(LOWPASS_CUTOFF)) THEN BEGIN
      filter = rolloff(Nx, lowpass_cutoff)
      ftl = filter * (zx*ft)
    ENDIF ELSE $
     ftl = zx*ft
  ENDELSE

  return, float(fft(ftl, /INVERSE))
END

