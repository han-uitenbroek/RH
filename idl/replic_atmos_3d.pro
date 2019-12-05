PRO replic_atmos_3d, multiAtmos, $
                     OUTPUT_BASE=output_base, $
                     NX=Nx, DX=dx, VX=vx, $
                     Ny=Ny, DY=dy, VY=vy, $
                     NZ1=Nz1

  a1D = read1datmos(multiAtmos)

  IRRADIATED = 0L  &  ZERO = 1L  &  PLANCK = 2L

  WGHT_PER_H = 1.4350D0
  AMU        = 1.6605402D-27

  s = size(a1D.nH)  &  NHydr = s[2]
  nHtotal = total(a1D.nH, 2)
  rho = WGHT_PER_H * nHtotal * AMU

  IF (a1D.scale EQ 'MASS_SCALE') THEN BEGIN
    height = dblarr(a1D.Ndep)
    FOR k=1, a1D.Ndep-1 DO $
     height[k] = height[k-1] - 2.0*(a1D.cmass[k] - a1D.cmass[k-1]) / $
     (rho[k-1] + rho[k])
    height = height - height[a1D.Ndep-1]
  ENDIF ELSE IF (a1D.scale EQ 'GEOMETRIC_SCALE') THEN $
   height = a1D.height

  IF (NOT keyword_set(NX)) THEN BEGIN
    Nx = 4
    print, FORMAT='(/" Using default value Nx = ", I3)', Nx
  ENDIF
  IF (NOT keyword_set(DX)) THEN BEGIN
    dx = 100.0
    print, FORMAT='(" Using default constant value dx = ", F6.1, " [km]")', dx
  ENDIF ELSE BEGIN
    IF (n_elements(dx) NE 1) THEN BEGIN
      print, "dx has to be a scalar"
      return
    ENDIF
  ENDELSE
  IF (NOT keyword_set(VX)) THEN BEGIN
    vx = dblarr(a1D.Ndep)
    print, FORMAT='(" Using default constant value vx = ", F6.1, " [km/s]")', $
     vx(0)
  ENDIF ELSE BEGIN
    IF (n_elements(vx) NE a1D.Ndep) THEN BEGIN
      print, "Number of elements in vx has to be Ndep"
      return
    ENDIF
  ENDELSE

  IF (NOT keyword_set(NY)) THEN BEGIN
    Ny = 4
    print, FORMAT='(/" Using default value Ny = ", I3)', Ny
  ENDIF
  IF (NOT keyword_set(DY)) THEN BEGIN
    dy = 100.0
    print, FORMAT='(" Using default constant value dx = ", F6.1, " [km]")', dy
  ENDIF ELSE BEGIN
    IF (n_elements(dy) NE 1) THEN BEGIN
      print, "dx has to be a scalar"
      return
    ENDIF
  ENDELSE
  IF (NOT keyword_set(VY)) THEN BEGIN
    vy = dblarr(a1D.Ndep)
    print, FORMAT='(" Using default constant value vy = ", F6.1, " [km/s]")', $
     vy(0)
  ENDIF ELSE BEGIN
    IF (n_elements(vy) NE a1D.Ndep) THEN BEGIN
      print, "Number of elements in vy has to be Ndep"
      return
    ENDIF
  ENDELSE

  IF (keyword_set(NZ1)) THEN BEGIN
    Nz = a1D.Ndep - Nz1
    z  = height[Nz1:*]
    T = a1D.T[Nz1:*]
    vturb = a1D.vturb[Nz1:*]
    n_elec = a1D.n_elec[Nz1:*]
    nh     = a1D.nh[Nz1:*, *]
    vx = vx[Nz1:*]  &  vy = vy[Nz1:*]
    vz = a1D.v[Nz1:*]
  ENDIF ELSE BEGIN
    Nz = a1D.Ndep
    z = height
    T = a1D.T
    vturb = a1D.vturb
    n_elec = a1D.n_elec
    nh = a1D.nH
    vz = a1D.v
  ENDELSE

  outputFile = $
   string(FORMAT='(A, "_", I3.3, "x", I3.3, "x", I3.3, ".atmos")', $
          output_base, Nx, Ny, Nz)
  openw, unit, outputFile, /GET_LUN, /XDR

  lBoundVal  = [ZERO,  PLANCK]

  writeu, unit, long([Nx, Ny, Nz, NHydr])
  writeu, unit, lBoundVal
  writeu, unit, double([dx, dy]), z

  writeu, unit, rebin(reform(T, 1, 1, Nz), Nx, Ny, Nz)
  writeu, unit, rebin(reform(n_elec, 1, 1, Nz), Nx, Ny, Nz)
  writeu, unit, rebin(reform(vturb, 1, 1, Nz), Nx, Ny, Nz)
  writeu, unit, rebin(reform(vx, 1, 1, Nz), Nx, Ny, Nz)
  writeu, unit, rebin(reform(vy, 1, 1, Nz), Nx, Ny, Nz)
  writeu, unit, rebin(reform(vz, 1, 1, Nz), Nx, Ny, Nz)

  writeu, unit, rebin(reform(nH, 1, 1, Nz, NHydr), $
                      Nx, Ny, Nz, NHydr)

  free_lun, unit
END
