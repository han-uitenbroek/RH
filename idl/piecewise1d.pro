FUNCTION w2, dtau

  w = dblarr(2)

  if (dtau < 5.0E-4) THEN BEGIN
    w[0] = dtau*(1.0 - 0.5*dtau)
    w[1] = dtau^2 * (0.5 - 0.33333333*dtau)
    return, w
  ENDIF

  IF (dtau > 50.0) THEN BEGIN
    w[1] = 1.0
    w[0] = 1.0
  ENDIF ELSE BEGIN
    expdt = exp(-dtau)
    w[0]  = 1.0 - expdt
    w[1]  = w[0] - dtau*expdt
  ENDELSE
  return, w
END

FUNCTION w3, dtau

  w = dblarr(3)

  IF (dtau < 5.0E-4) THEN BEGIN
    w[0]   = dtau*(1.0 - 0.5*dtau)
    delta  = dtau^2
    w[1]   = delta*(0.5 - 0.33333333*dtau)
    delta  = delta * dtau
    w[2]   = delta*(0.33333333 - 0.25*dtau)
    return, w
  ENDIF

  IF (dtau > 50.0) THEN begin
    w[1] = 1.0
    w[0] = 1.0
    w[2] = 2.0
  ENDIF ELSE BEGIN
    expdt = exp(-dtau)
    w[0]  = 1.0 - expdt
    w[1]  = w[0] - dtau*expdt
    w[2]  = 2.0*w[1] - SQ(dtau) * expdt
  ENDELSE

  return, w
END

PRO piecewise1d, height, reverse, xmu, chi, S, I, Psi

  ;; Piecewise quadratic integration of the transfer equation.

  Ndep = n_elements(height)
  I = dblarr(Ndep)
  Psi = I

  IF (reverse) THEN begin

    ;; --- From TOP to BOTTOM --                       -------------- ;;

    I[0] = 0.0
    dtau_uw = 0.5/xmu * (chi[0] + chi[1]) * (height[0] - height[1])
    dS_uw   = (S[0] - S[1]) / dtau_uw

    FOR k = 1, Ndep-2 DO BEGIN
      dtau_dw = 0.5/xmu * (chi[k] + chi[k+1]) *	(height[k] - height[k+1])
      dS_dw   = (S[k] - S[k+1]) / dtau_dw
      w = w3(dtau_uw)

      c1 = (dS_uw*dtau_dw + dS_dw*dtau_uw)
      c2 = (dS_uw - dS_dw)
      I[k] = w[0]*S[k] + (w[1]*c1 + w[2]*c2) / (dtau_uw + dtau_dw) + $
	I[k-1] * exp(-dtau_uw)

      c1 = dtau_uw - dtau_dw
      Psi[k] = w[0] + (w[1]*c1 - w[2]) / (dtau_uw * dtau_dw)

      dS_uw   = dS_dw
      dtau_uw = dtau_dw
    ENDFOR 
    w = w2(dtau_uw)
    c1 = (S[Ndep-2] - S[Ndep-1]) / dtau_uw
    I[Ndep-1] = w[0]*S[Ndep-1] + w[1]*c1 + I[Ndep-2] * exp(-dtau_uw)
    Psi[Ndep-1] = w[0] - w[1] / dtau_uw
  ENDIF ELSE BEGIN

    ;; --- From BOTTOM to TOP --                       -------------- ;;

    I[Ndep-1] = S[Ndep-1];
    dtau_uw = 0.5/xmu * (chi[Ndep-2] + chi[Ndep-1]) * $
      (height[Ndep-2] - height[Ndep-1])
    dS_uw   = (S[Ndep-1] - S[Ndep-2]) / dtau_uw

    for k = Ndep-2, 1, -1 DO BEGIN
      dtau_dw = 0.5/xmu * (chi[k] + chi[k-1]) *	(height[k-1] - height[k])
      dS_dw   = (S[k] - S[k-1]) / dtau_dw
      w = w3(dtau_uw)

      c1 = (dS_uw*dtau_dw + dS_dw*dtau_uw)
      c2 = (dS_uw - dS_dw)
      I[k] = w[0]*S[k] + (w[1]*c1 + w[2]*c2) / (dtau_uw + dtau_dw) + $
	I[k+1] * exp(-dtau_uw)

      c1 = dtau_uw - dtau_dw
      Psi[k] = w[0] + (w[1]*c1 - w[2]) / (dtau_uw * dtau_dw)

      dS_uw   = dS_dw
      dtau_uw = dtau_dw
    ENDFOR
    w = w2(dtau_uw)
    c1 = (S[1] - S[0]) / dtau_uw
    I[0] = w[0]*S[0] + w[1]*c1 + I[1] * exp(-dtau_uw)
    Psi[0] = w[0] - w[1] / dtau_uw
  ENDELSE

END
