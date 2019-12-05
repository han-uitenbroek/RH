PRO replic_atmos, multiAtmos, OUTPUT_BASE=output_base, NX=Nx, DX=dx, VX=vx, $
     NZ1=Nz1

@ATMOS
  read_atmos, multiAtmos

  CM_TO_M = 1.0E-02
  KM_TO_M = 1.0E+03
  G_TO_KG = 1.0E-03

  FIXED      = 0L  &  PERIODIC = 1L
  IRRADIATED = 0L  &  ZERO = 1L  &  THERMALIZED = 2L

  WGHT_PER_H = 1.4350
  AMU        = 1.6605402E-27

  s = size(nH)  &  NHydr = s(2) + 1
  nHtotal = (total(nH, 2) + np) / (CM_TO_M)^3
  rho     = WGHT_PER_H * nHtotal * AMU
  cmass   = 10^(cmass) * (G_TO_KG / CM_TO_M^2);

  height = fltarr(Ndep)
  FOR k=1, Ndep-1 DO $
   height(k) = height(k-1) - 2.0*(cmass(k) - cmass(k-1)) / (rho(k-1) + rho(k))
  height = height - height(Ndep-1)

  IF (NOT keyword_set(NX)) THEN BEGIN
    Nx = 3
    print, FORMAT='(/" Using default value Nx = ", I3)', Nx
  ENDIF
  IF (NOT keyword_set(DX)) THEN BEGIN
    dx = 100.0 + fltarr(Nx)
    print, FORMAT='(" Using default constant value dx = ", F6.1, " [km]")', $
     dx(0)
  ENDIF ELSE BEGIN
    IF (n_elements(dx) NE Nx) THEN BEGIN
      print, "Number of elements in dx has to be Nx"
      return
    ENDIF
  ENDELSE
  IF (NOT keyword_set(VX)) THEN BEGIN
    vx = fltarr(Ndep)
    print, FORMAT='(" Using default constant value vx = ", F6.1, " [km/s]")', $
     vx(0)
  ENDIF ELSE BEGIN
    IF (n_elements(vx) NE Ndep) THEN BEGIN
      print, "Number of elements in vx has to be Ndep"
      return
    ENDIF
  ENDELSE

  IF (keyword_set(NZ1)) THEN BEGIN
    Nz = Ndep - Nz1
    z  = height(Nz1:*) / KM_TO_M

    temp = temp(Nz1:*)  &  vturb = vturb(Nz1:*)

    n_elec = n_elec(Nz1:*) / (CM_TO_M)^3
    nh     = nh(Nz1:*, *) / (CM_TO_M)^3
    np     = np(Nz1:*) / (CM_TO_M)^3

    vx = vx(Nz1:*)
    vz = vel(Nz1:*)
  ENDIF ELSE BEGIN
    Nz = Ndep
    z  = height / KM_TO_M

    n_elec = n_elec / (CM_TO_M)^3
    nh     = nh / (CM_TO_M)^3
    np     = np / (CM_TO_M)^3

    vz = vel
  ENDELSE

  outputFile = $
   string(FORMAT='(A, "_", I2.2, "x", I2.2, ".atmos")', output_base, Nx, Nz)
  openw, unit, outputFile, /GET_LUN, /XDR

  lBoundCond = PERIODIC
  lBoundVal  = [ZERO,  THERMALIZED]

  writeu, unit, long([Nx, Nz, NHydr])
  writeu, unit, lBoundCond, lBoundVal
  writeu, unit, double(dx), double(z)

  writeu, unit, double(rebin(reform(temp, 1, Nz, /OVERWRITE), Nx, Nz))
  writeu, unit, double(rebin(reform(n_elec, 1, Nz, /OVERWRITE), Nx, Nz))
  writeu, unit, double(rebin(reform(vturb, 1, Nz, /OVERWRITE), Nx, Nz))
  writeu, unit, double(rebin(reform(vx, 1, Nz, /OVERWRITE), Nx, Nz))
  writeu, unit, double(rebin(reform(vz, 1, Nz, /OVERWRITE), Nx, Nz))

  writeu, unit, double(rebin(reform(nH, 1, Nz, NHydr-1, /OVERWRITE), $
                             Nx, Nz, NHydr-1))
  writeu, unit, double(rebin(reform(np, 1, Nz, /OVERWRITE), Nx, Nz))

  free_lun, unit
END
