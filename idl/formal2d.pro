
  angleSet = 4L  &  mu = 6L
  Ninclination = 4L  &  Nazimuth = 3L
  Nx = 250L  &  Nz = 100L
  dx = 12.5 + dblarr(Nx)
  zmin = 0.0  &  zmax = 600.0

  dS_x = 0.0  &  dS_z = 0.0

  FIXED = 0L  &  PERIODIC = 1L
  hboundCond = PERIODIC

  x = dblarr(Nx)  &  FOR l=0, Nx-2 DO x(l+1) = x(l) + dx(l)
  z = zmax + (zmin - zmax)*dindgen(Nz)/double(Nz - 1)

  chi = dblarr(Nx, Nz) + 1.0E-6
;;  chi = dblarr(Nx, Nz) + (1.0 + fltarr(Nx)) # $
;;        (1.0E-7 + z*(1.0E-9 + z*(1.0E-12 + z*1.0E-15)))
;;  S   = (1.0 + dS_x*cos(4*!pi*dindgen(Nx)/(Nx-1))) # $
;;   (1.0 + dS_z*sin(!pi*dindgen(Nz)/(Nz-1)))

  S = dblarr(Nx, Nz)
  S[*, 40:49] = 1.0

;;  tau = z * chi[0, *] * 1.0E3
;;  S = (1.0 + fltarr(Nx)) # (1.0 + tau*(0.1 + tau*0.01))

  Itop  = dblarr(Nx)  &  Ibottom = dblarr(Nx)
  Ileft = dblarr(Nz)  &  Iright  = dblarr(Nz)
  index = 20 + (dindgen(6)-3)

;  Itop(index) = (sin(!PI*dindgen(6)/5.0) > 0.0) 
;  Ibottom(index) = (sin(!PI*dindgen(6)/5.0) > 0.0)
  Ibottom[index] = 1.0
;  Ileft(index) = (sin(!PI*dindgen(6)/5.0) > 0.0)

;;  openw, pipe, 'formal.in', /GET_LUN
  spawn, './formal2d', UNIT=pipe, /NOSHELL

  writeu, pipe, angleSet, mu
  IF (angleSet EQ 1) THEN writeu, pipe, Ninclination, Nazimuth

  writeu, pipe, hboundCond
  writeu, pipe, Nx, Nz, dx, z
  writeu, pipe, Itop, Ibottom;;, Ileft, Iright
  writeu, pipe, chi, S

  I_d   = dblarr(Nx, Nz)  &  I_u   = I_d
  Psi_d = dblarr(Nx, Nz)  &  Psi_u = Psi_d

;;  free_lun, pipe
;;  openr, pipe, 'formal.out', /GET_LUN

  readu, pipe, I_d, Psi_d, I_u, Psi_u
  close, pipe
  free_lun, pipe

  pswindow, XSIZE=400 + 150, YSIZE=400 + 100 
  panel, scaleimg_idl(I_u, 400, 400), x, z, $
   XPOS=50, YPOS=50, XTITLE='x [km]', YTITLE='z [km]', /ORDER, TITLE='Top'

;  FOR k=0,Nz-1 DO oplot, [x(0), x(Nx-1)], [z(k), z(k)], $
;;   THICK=1+(k MOD 10 EQ 0), COLOR=240B
;  FOR l=0,Nx-1 DO oplot, [x(l), x(l)], [z(0), z(Nz-1)], $
;   THICK=1+(l MOD 10 EQ 0), COLOR=240B
END
