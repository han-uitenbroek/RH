function combine_ns_fal, ns_model, fal_model, OUTPUT_NAME=output_name, $
                         H0=h0, H1=h1

  IF (n_params() LT 2) THEN BEGIN
    print, "Usage:"
    print, "  output = COMBINE_NS_FAL(ns_model, fal_model, $"
    print, "                          OUTPUT_NAME=output_name, H0=h0, H1=h1)"
    return, 0
  ENDIF

  IF (NOT keyword_set(H0)) THEN h0 = 500.0D+0
  IF (NOT keyword_set(H1)) THEN h1 = 650.0D+0

  IF (h0 GE h1) THEN BEGIN
    print, " H1 should be greater than H0"
    return, 0
  ENDIF

@ATMOS

  CM_TO_M = 1.0E-2

  nsa = read2datmos(ns_model)
  IF (h1 GT nsa.z[0]) THEN BEGIN
    print, " H1 should be below top of Nordlund-Stein atmosphere"
    return, 0
  ENDIF

  readfal, fal_model
  height = -height
  nHtot = (total(nh, 2) + np) / CM_TO_M^3
  n_elec = n_elec / CM_TO_M^3

  z = double([height[where(height GT h1)], nsa.z[where(nsa.z LE h1)]])
  Nz = n_elements(z)  &  Nx = nsa.Nx

  ;; Define structure combined atmosphere

  a_ext = {Nx: nsa.Nx, Nz: Nz, NHydr: 1L, boundary: nsa.boundary, $
           dx: nsa.dx, z: z, T: dblarr(Nx, Nz), n_elec: dblarr(Nx, Nz), $
           vturb: dblarr(Nx, Nz), vx: dblarr(Nx, Nz), vz: dblarr(Nx, Nz), $
           nH: dblarr(Nx, Nz)}

  ;; Copy Nordlund-Stein values in bottom of atmosphere.

  ns_index = where(z LE h0)
  copy     = where(nsa.z LE h0)
  a_ext.T[*, ns_index]      = nsa.T[*, copy]
  a_ext.n_elec[*, ns_index] = nsa.n_elec[*, copy]
  a_ext.vturb[*, ns_index]  = nsa.vturb[*, copy]
  a_ext.vx[*, ns_index]     = nsa.vx[*, copy]
  a_ext.vz[*, ns_index]     = nsa.vz[*, copy]
  a_ext.nH[*, ns_index]     = nsa.nH[*, copy]

  ;; Replicate FAL model values horizontally in top of atmosphere.

  fal_index = where(z GT h1)
  copy = where(height GT h1)  &  Ncopy = n_elements(copy)
  a_ext.T[*, fal_index]      = $
   rebin(reform(temp[copy], 1, Ncopy), Nx, Ncopy) 
  a_ext.n_elec[*, fal_index] = $
   rebin(reform(n_elec[copy], 1, Ncopy), Nx, Ncopy) 
  a_ext.vturb[*, fal_index]  = $
   rebin(reform(vturb[copy], 1, Ncopy), Nx, Ncopy)
  a_ext.nH[*, fal_index]     = $
   rebin(reform(nHtot[copy], 1, Ncopy), Nx, Ncopy)

  ;; Interpolate between heights h0 and h1 and add linear combination
  ;; of the two models.

  mix_index = where(a_ext.z LE h1  AND  a_ext.z GT h0)
  Nmix = n_elements(mix_index)
  ns_index = where(nsa.z LE h1  AND  nsa.z GT h0)

  ;; Fraction of FAL model to be applied

  f = rebin(reform((z[mix_index] - h0)/(h1 - h0), 1, Nmix), Nx, Nmix)
  tabinv, height, a_ext.z[mix_index], z_eff

  fal = interpolate(temp, z_eff, CUBIC=0.5)
  fal = rebin(reform(fal, 1, Nmix), Nx, Nmix)
  a_ext.T[*, mix_index] = f * fal + (1.0 - f) * nsa.T[*, ns_index]

  fal = interpolate(alog(n_elec), z_eff, CUBIC=0.5)
  fal = rebin(reform(fal, 1, Nmix), Nx, Nmix)
  a_ext.n_elec[*, mix_index] = f * exp(fal) + $
   (1.0 - f) * nsa.n_elec[*, ns_index]

  fal = interpolate(alog(nHtot), z_eff, CUBIC=0.5)
  fal = rebin(reform(fal, 1, Nmix), Nx, Nmix)
  a_ext.nH[*, mix_index] = f * exp(fal) + (1.0 - f) * nsa.nH[*, ns_index]

  a_ext.vx[*, mix_index] = (1.0 - f) * nsa.vx[*, ns_index]
  a_ext.vz[*, mix_index] = (1.0 - f) * nsa.vz[*, ns_index]

  fal = interpolate(vturb, z_eff, CUBIC=0.5)
  fal = rebin(reform(fal, 1, Nmix), Nx, Nmix)
  a_ext.vturb[*, mix_index] = f * fal

  IF (keyword_set(OUTPUT_NAME)) THEN BEGIN
    openw, lun, OUTPUT_NAME, /GET_LUN, /XDR
    writeu, lun, a_ext
    free_lun, lun
  ENDIF

  return, a_ext
END
