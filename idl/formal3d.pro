
  angleSet = 0L  &  mu = 0L
  Ninclination = 4L  &  Nazimuth = 2L
  Nx = 25L  &  Ny = 10L  &  Nz =100L
  Nplane = Nx * Ny
  dx = 2.5*12.5D0
  dy = 2.5*12.5D0
  zmin = 0.0  &  zmax = 600.0


;;  z = reverse(zgrid(Nz, zmin, zmax, 1.00))
  z = zmax + (zmin - zmax)*dindgen(Nz)/double(Nz - 1)

  chi = dblarr(Nx, Ny, Nz) + 1.0D-5
  tau = z * chi[0, 0, *] * 1.0D3
  S   = rebin(reform(1.0 + tau*(0.1 + tau*0.01), 1, 1, Nz), Nx, Ny, Nz)

  Itop_x = dblarr(Nx)
;;  Itop_x[20+dindgen(6)-3] = abs(sin(!PI*dindgen(6)/5.0))
  Itop_y = dblarr(Ny)
;;  Itop_y[20+dindgen(6)-3] = abs(sin(!PI*dindgen(6)/5.0))
  Itop = (Itop_x # Itop_y)

  Ibottom = Itop

;;openw, pipe, 'formal.in', /GET_LUN
  spawn, './formal3d', UNIT=pipe, /NOSHELL

  writeu, pipe, angleSet, mu
  IF (angleSet EQ 1) THEN writeu, pipe, Ninclination, Nazimuth
  writeu, pipe, Nx, Ny, Nz, dx, dy, z
  writeu, pipe, Itop, Ibottom
  writeu, pipe, chi, S
;;stop
  I_d   = dblarr(Nx, Ny, Nz)  &  I_u   = I_d
  Psi_d = dblarr(Nx, Ny, Nz)  &  Psi_u = Psi_d

;;free_lun, pipe
;;openr, pipe, 'formal.out', /GET_LUN

  readu, pipe, I_d, Psi_d, I_u, Psi_u
  close, pipe
  free_lun, pipe

END
