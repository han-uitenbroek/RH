PRO make2datmos, file, XCUT=xcut, YCUT=ycut, HGRID=hgrid, MICRO=micro, $
     DEBUG=debug
;;
;;  Create 2D cut through granulation cube of Nordlund 
;;
;;  For the microturbulence there are two possibilities.
;;  1. Standard option: zero microturbulence. Note that Nordlund uses this.
;;     It is very easy to introduce another CONSTANT value for microturbulence
;;     by changing variable vmicro.
;;  2. Read microturbulence as function of height from file. Use e.g. the
;;     VAL values, which can easily be extracted from some out file of Multi.
;;
;;  Since Nordlund has all ionization equilibria in LTE, use rho to 
;;  generate input LTE hydrogen populations.
;;
  IF (n_params() EQ 0) THEN BEGIN
    print,'Usage: make2datmos,file,xcut=xcut,ycut=ycut,hgrid=hgrid,micro=micro'
    print,'Parameters:'
    print,' file  : name of the atmosphere file in RH format'
    print,'Keywords:'
    print,' xcut  : Value of the x-index along which to take cut'
    print,' ycut  : Value of the y-index along which to take cut'
    print,'   Note that either xcut OR ycut may be specified, but not both'
    print,' hgrid : The new height grid (in km) to be used in the model'
    print,'         A good choice is: 650.0 - (650.0+200.0) * findgen(64)/63.0'
    print,' micro : The microturbulent velocity to be used'
    print,'         This may be a scalar, in which case it is '
    print,'         assumed depth-independent'
    print,'         Or it may be an array with the same length as h or hgrid'
    print,'         resp. nz or nz2. In case nz and nz2 are equal, it is '
    print,'         assumed that micro is specified on the final grid'
    print,' debug  : If set, stops the execution before exiting the routine'
    return
  ENDIF
;; 
;;----------- First get the full data cube -------------------
;; All arrays (te,rho,nne,abso,ux,uy,uz) have sizes (nx,ny,nz)
;; with nx=ny=63 and nz=47
;;
@readaake
;;
;;----------- Check whether only XCUT or YCUT is specified
;;
  IF (keyword_set(XCUT) AND keyword_set(YCUT)) THEN BEGIN
    print,'INPUT ERROR: specify either XCUT or YCUT, not both simultaneously'
    return
  ENDIF
  IF (NOT(keyword_set(XCUT)) AND NOT(keyword_set(YCUT))) THEN BEGIN 
    print,'INPUT ERROR: specify either XCUT or YCUT'
    return
  ENDIF
;;
;;----------- Prepare the correct cut from the data cube
;; Careful with the macroscopic velocities; for XCUT specified, 
;; horizontal velocity is UY component, and the other way around.
;;
  IF keyword_set(XCUT) THEN BEGIN
    x_p         = y
    nx_p        = n_elements(x_p)
    temperature = reform(te(xcut, *, *), nx_p,nz)
    n_elec      = reform(nne(xcut, *, *), nx_p,nz)
    velx        = reform(uy(xcut, *, *), nx_p,nz)
    velz        = reform(uz(xcut, *, *), nx_p,nz)
    dens        = reform(rho(xcut, *, *), nx_p,nz)
  ENDIF ELSE BEGIN
    x_p         = x
    nx_p        = n_elements(x_p)
    temperature = reform(te(*,ycut,*),nx_p,nz)
    n_elec      = reform(nne(*,ycut,*),nx_p,nz)
    velx        = reform(ux(*,ycut,*),nx_p,nz)
    velz        = reform(uz(*,ycut,*),nx_p,nz)
    dens        = reform(rho(*,ycut,*),nx_p,nz)
  ENDELSE

;; note the sign change in velz that is necessary to get + = upflow

  velz = -velz

;;------------ If new height grid then interpolate data in height

  z_p  = h
  nz_p = n_elements(z_p)
  IF keyword_set(hgrid) THEN BEGIN
    z_p    = hgrid 
    nz_p   = n_elements(z_p)
    te_p   = fltarr(nx_p, nz_p)
    ne_p   = te_p
    velx_p = te_p
    velz_p = te_p
    rho_p  = te_p
    FOR i=0,nx_p-1 DO BEGIN
      te_p(i,0:nz_p-1)   = interpol(temperature(i,0:nz-1),h,z_p)
      ne_p(i,0:nz_p-1)   = exp(interpol(alog(n_elec(i,0:nz-1)),h,z_p))
      rho_p(i,0:nz_p-1)  = exp(interpol(alog(dens(i,0:nz-1)),h,z_p))
      velx_p(i,0:nz_p-1) = interpol(velx(i,0:nz-1),h,z_p)
      velz_p(i,0:nz_p-1) = interpol(velz(i,0:nz-1),h,z_p)
    ENDFOR
  ENDIF ELSE BEGIN
    te_p = temperature
    ne_p = n_elec
    rho_p = dens
    velx_p = velx
    velz_p = velz
  ENDELSE
  IF keyword_set(MICRO) THEN BEGIN
    nmicro = n_elements(micro)
    IF (nmicro EQ 1) THEN micro=fltarr(nz_p)+micro
    IF (nmicro EQ nz_p) THEN BEGIN
      vturb_p = micro 
    ENDIF ELSE BEGIN
      vturb_p = interpol(micro,h,z_p)
    ENDELSE
  ENDIF ELSE BEGIN
    vturb_p = fltarr(nz_p)
  ENDELSE
  vturb_p = (fltarr(nx_p)+1)#vturb_p

;;-------------- Take care of hydrogen pops, run code with hydrogen set
;;               to LTE

  NHydr = 1L
  XMY   = 1.433896
  AMU   = 1.660531E-24
  WH    = 1.008

  nH = rho_p / (XMY * WH * AMU)

;;-------------- Conversions of units

  CM_TO_M = 1.0E-02
  ne_p = ne_p / (CM_TO_M)^3
  nH   = nH   / (CM_TO_M)^3

;;-------------- dx rather than x is needed

  dx = shift(x, -1) - x
  dx(Nx-1) = dx(0)

  openw, unit, /GET_LUN, file, /XDR

;;-------------- Boundary conditions

  FIXED = 0L  &  PERIODIC = 1L
  IRRADIATED =  0L  &  ZERO = 1L  &  THERMALIZED = 2L
  lBoundCond = PERIODIC
  lBoundVal  = [ZERO,  THERMALIZED]

  IF keyword_set(debug) THEN stop

;;-------------- Now write to file

  writeu, unit, long([Nx, Nz_p, NHydr])
  writeu, unit, lBoundCond, lBoundVal
  writeu, unit, double(dx), double(z_p)
  writeu, unit, double(Te_p)
  writeu, unit, double(ne_p)
  writeu, unit, double(vturb_p)
  writeu, unit, double(velx_p)
  writeu, unit, double(velz_p)
  writeu, unit, double(nH)

  free_lun, unit
END
