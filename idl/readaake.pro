;=======================================================

;
; Read the XDR formatted granulation cube (Nordlund 1989)
;
datadir = '~/src/rh/Atmos'
;
openr, unit, datadir + '/aake.xdr', /XDR, /GET_LUN
nx = 0L
ny = 0L
nz = 0L
readu, unit, nx, ny, nz
x = fltarr(nx)
y = fltarr(ny)
h = fltarr(nz)
readu, unit, x, y, h
;
te   = fltarr(nx,ny,nz)
rho  = te
nne  = te
abso = te
ux   = te
uy   = te
uz   = te
;
readu, unit, te, rho, nne, abso, ux, uy, uz
free_lun,unit
