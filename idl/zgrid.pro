FUNCTION zgrid, Nz, zbottom, ztop, alpha

  IF (alpha LT 1.0) THEN BEGIN
    print, "Scaling parameter alpha has to be larger than 1.0, not ", alpha
    return, 0
  ENDIF
  IF (alpha EQ 1.0) THEN $
    return, zbottom + dindgen(Nz) * (ztop - zbottom) / float(Nz - 1)

  z1 = (ztop - zbottom) / (alpha^(Nz-1) - 1.0)
  z0 = zbottom - z1

  return, z0 + z1 * alpha^dindgen(Nz)
END

  
