
  angleSet = 4L
  Nx = 100L  &  Nz = 100L
  zmin = 0.0  &  zmax = 600  &  H = 60.0

  Nlambda  = 51L
  NmaxIter = 50L  &  iterLimit = 1.0D-4
  Ngdelay  =  6L  &  Ngorder   = 3L   &  Ngperiod  = 5L

  FIXED = 0L  &  PERIODIC = 1L
  boundCond = PERIODIC

  IRRADIATED = 0L  &  ZERO = 1L  &  THERMALIZED = 2L  &  REFLECTIVE = 3L
  boundVal = [ZERO, IRRADIATED]

  dx = 12.2 + dblarr(Nx)
  x  = dblarr(Nx)
  FOR l=0, Nx-2 DO x(l+1) = x(l) + dx(l)
  z  = double(zmax + (zmin - zmax)*findgen(Nz)/float(Nz - 1))

  dBp_x = 0.15  &  dBp_z = 0.05
;;  dBp_x = 0.0  &  dBp_z = 0.0
  
  eps0  = 1.0D-04

  rho0  = 2.0D+3
  rho = double(rho0 * 10.0^(-z/H))
  chi = dblarr(Nx, Nz)
  FOR n=0, Nz-1 DO chi(*, n) = rho(n)

  Bp = (1.0 + dBp_x*cos(4*!pi*dindgen(Nx)/(Nx-1))) # $
   (1.0 + dBp_z*sin(!pi*dindgen(Nz)/(Nz-1)))
  epsilon = dblarr(Nx, Nz) + eps0

  Adamp  = 1.0D-03
  lambda = 15.0 * dindgen(Nlambda)/double(Nlambda)

  Ileft  = dblarr(Nz, Nlambda)
;;  Ileft[0:9, *] = Ileft[0:9, *] * $
;;   ((reverse(cos(!PI *dindgen(10)/9.0)) + 1)/2. # (dblarr(Nlambda) + 1.0))
  Iright = dblarr(Nz, Nlambda)
;;  FOR la=0,Nlambda-1 DO $
;;   Ileft(*, la) = Ileft(*, la) + 0.2 * abs(sin(2*!PI * findgen(Nz)/(Nz-1)))
  Ibottom = rebin(reform(Bp[*, Nz-1], Nx, 1), Nx, Nlambda)

  openw, fpin, /GET_LUN, '2dinput.dat', /XDR
  writeu, fpin, angleSet, Nx, Nz, NmaxIter, Ngdelay, Ngorder, Ngperiod, Nlambda
  writeu, fpin, iterLimit, Adamp
  writeu, fpin, boundCond, boundVal
  writeu, fpin, dx, z
  writeu, fpin, chi, Bp, epsilon
  writeu, fpin, Ibottom
;;  writeu, fpin, Ileft, Ileft
  writeu, fpin, lambda
  free_lun, fpin

  spawn, './solve2d', /NOSHELL, UNIT=fpout

  Nrays_in_set = [1, 2, 6, 6, 12, 12, 20]
  Nmu = Nrays_in_set(angleSet)

  S = dblarr(Nx, Nz)
  readu, fpout, S
  free_lun, fpout

  pswindow, 0, xsize = 500, ysize = 600
  loadct, 1, /SILENT
  
  panel, scaleimg_idl(Bp, 350, 250), [x(0), x(Nx-1)], [z(Nz-1), z(0)], $
   xpos=50, ypos=325, /ORDER, ytitle='height [km]', title='B'
  panel, scaleimg_idl(S, 350, 250), $
   [x(0), x(Nx-1)], [z(Nz-1), z(0)], $
   xpos=50, ypos=50, /ORDER, xtitle='x [km]', ytitle='height [km]', $
   title='S'

  window, 2, xs=1000, ys=1200
  shade_surf, S
END
