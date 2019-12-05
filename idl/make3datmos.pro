PRO make3datmos, file, HGRID=hgrid, MICRO=micro, DEBUG=debug
;;
;;  Create 3D atmospheric input file from granulation cube of Nordlund 
;;
;;  For the microturbulence there are two possibilities.
;;  1. Standard option: zero microturbulence. Note that Nordlund uses this.
;;     It is very easy to introduce another CONSTANT value for microturbulence
;;     by changing variable vmicro.
;;  2. Read microturbulence as function of height from file. Use e.g. the
;;     VAL values, which can easily be extracted from some out file of Multi.
;;
;;  Nordlund has all ionization equilibria in LTE.


  IF (n_params() EQ 0) THEN BEGIN
    print,'Usage: make3datmos, file, HGRID=hgrid, MICRO=micro, /DEBUG'
    print,'Parameters:'
    print,' file  : name of the atmosphere file in RH format'

    print,'Keywords:'
    print,' HGRID : The new height grid (in km) to be used in the model'

    print,' MICRO : The microturbulent velocity to be used'
    print,'         This may be a scalar, in which case it is '
    print,'         assumed depth-independent'
    print,'         Or it may be an array with the same length as h or hgrid'
    print,'         resp. nz or nz2. In case nz and nz2 are equal, it is '
    print,'         assumed that micro is specified on the final grid'

    print,' DEBUG : If set, stops the execution before exiting the routine'
    return
  ENDIF

;;----------- First get the full data cube -------------------

;; All arrays (Te, rho, nne, abso, ux, uy, uz) have sizes (Nx, Ny, Nz)
;; with Nx = Ny = 63 and Nz = 47

@readaake

;; note the sign change in uz that is necessary to get + = upflow

  uz = -uz

;;------------ If new height grid then interpolate data in height

  IF keyword_set(HGRID) THEN BEGIN
    z  = hgrid 
    Nh = Nz
    Nz = n_elements(z)
  ENDIF ELSE $
   z = h

  Te_p  = fltarr(Nz, Nx*Ny)
  nne_p = Te_p
  rho_p = Te_p
  ux_p  = Te_p
  uy_p  = Te_p
  uz_p  = Te_p
  
  IF (keyword_set(HGRID)) THEN BEGIN
    Te  = transpose(reform(Te, Nx*Ny, Nh, /OVERWRITE))
    nne = transpose(reform(nne, Nx*Ny, Nh, /OVERWRITE))
    rho = transpose(reform(rho, Nx*Ny, Nh, /OVERWRITE))
    ux  = transpose(reform(ux, Nx*Ny, Nh, /OVERWRITE))
    uy  = transpose(reform(uy, Nx*Ny, Nh, /OVERWRITE))
    uz  = transpose(reform(uz, Nx*Ny, Nh, /OVERWRITE))

    FOR i=0, Nx*Ny-1 DO BEGIN
      Te_p(*, i)  = interpol(Te(*, i), h, z)
      nne_p(*, i) = exp(interpol(alog(nne(*, i)), h, z))
      rho_p(*, i) = exp(interpol(alog(rho(*, i)), h, z))
      ux_p(*, i)  = interpol(ux(*, i), h, z)
      uy_p(*, i)  = interpol(uy(*, i), h, z)
      uz_p(*, i)  = interpol(uz(*, i), h, z)
    ENDFOR
  ENDIF

  IF keyword_set(MICRO) THEN BEGIN
    nmicro = n_elements(micro)
    IF (nmicro EQ 1) THEN micro = fltarr(Nz) + micro
    IF (nmicro EQ Nz) THEN BEGIN
      vturb = micro 
    ENDIF ELSE BEGIN
      vturb = interpol(micro, h, z)
    ENDELSE
  ENDIF ELSE BEGIN
    vturb = fltarr(Nz)
  ENDELSE
  vturb = reform(fltarr(Nx*Ny) # vturb, Nx, Ny, Nz, /OVERWRITE)


;;-------------- Take care of hydrogen pops, run code with hydrogen set
;;               to LTE

  CM_TO_M = 1.0E-02

  NHydr = 1L
  XMY   = 1.433896
  AMU   = 1.660531E-24
  WH    = 1.008

  IF (keyword_set(HGRID)) THEN BEGIN
    nne_p = nne_p / (CM_TO_M)^3
    nH_p  = rho_p / (XMY * WH * AMU) / (CM_TO_M)^3
  ENDIF ELSE BEGIN
    nne = nne / (CM_TO_M)^3
    nH  = rho / (XMY * WH * AMU) / (CM_TO_M)^3
  ENDELSE

;;-------------- dx rather than x is needed

  dx = x(1) - x(0)
  dy = y(1) - y(0)

;;-------------- Boundary conditions

  IRRADIATED =  0L  &  ZERO = 1L  &  THERMALIZED = 2L
  lBoundVal  = [ZERO,  THERMALIZED]

  IF (keyword_set(DEBUG)) THEN stop

;;-------------- Now write to file

  openw, unit, /GET_LUN, file, /XDR

  writeu, unit, long([Nx, Ny, Nz, NHydr])
  writeu, unit, lBoundVal
  writeu, unit, double(dx), double(dy), double(z)

  IF (keyword_set(HGRID)) THEN BEGIN
    Te  = transpose(Te_p)
    nne = transpose(nne_p)
    nH  = transpose(nH_p)
    ux  = transpose(ux_p)
    uy  = transpose(uy_p)
    uz  = transpose(uz_p)
  ENDIF

  writeu, unit, double(Te)
  writeu, unit, double(nne)
  writeu, unit, double(vturb)
  writeu, unit, double(ux)
  writeu, unit, double(uy)
  writeu, unit, double(uz)
  writeu, unit, double(nH)

  free_lun, unit
END
