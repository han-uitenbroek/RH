FUNCTION tauone, z, chi, TAUVALUE=tauvalue

  IF (n_params() LT 2) THEN BEGIN
    print, "Usage: result = tauone(z, chi, TAUVALUE=tauvalue)"
    return, 0
  ENDIF
  IF (NOT keyword_set(TAUVALUE)) THEN tauvalue = 1.0

  Nz = n_elements(z)
  s  = size(chi)
  IF (s[0] EQ 2) THEN BEGIN
    Nx = s[1]

    IF (s[2] NE Nz) THEN BEGIN
      print, "Second dimension of chi should be the same as that of z"
      return, 0
    ENDIF
  ENDIF

  IF (s[0] EQ 1) THEN BEGIN
    tau = fltarr(Nz)
    dz = shift(z, 1) - z
    FOR k=1, Nz-1 DO $
      tau[k] = tau[k-1] + 0.5*(chi[k-1] + chi[k]) * dz[k]

    tabinv, tau, tauvalue, zeff
    zeff = zeff
  ENDIF ELSE BEGIN
    tau = fltarr(Nx, Nz)
    dz = shift(z, 1) - z

    FOR k=1, Nz-1 DO $
      tau[*, k] = tau[*, k-1] + 0.5*(chi[*, k-1] + chi[*, k]) * dz[k]

    tau  = transpose(tau)
    zeff = fltarr(Nx)

    FOR l=0, Nx-1 DO BEGIN
      tabinv, tau[*, l], tauvalue, z_eff
      zeff[l] = z_eff
    ENDFOR
  ENDELSE

  return, zeff
END

 
    
