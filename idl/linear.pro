function DSUP, r,i,j
;
ni= n_elements(i)
nj= n_elements(j)
sr= size(r)  &  type= sr(sr(0)+1)
;
IF (nj eq 1) THEN return, r(i,scalar(j))
rij= make_array(dimension=[ni,nj],type=type)
FOR n=0,nj-1 DO rij(*,n)= r(i,j(n))
;
return, rij
end
function LINEAR, rr,ix,jy
;+
; NAME:
;       LINEAR
; PURPOSE:
;       (Bi)linear interpolation of 1- or 2-D array on effective 
;        subscripts IX and JY
; CATEGORY:
;       Interpolation
; CALLING SEQUENCE:
;       Z= BILINEAR(r,ix,iy)
; INPUTS:
;       R   --  1- or 2-D array
;       IX  --  Effective indices of first dimension
;       JY  --  Effective indices of second dimension
; KEYWORD PARAMETERS:
; OUTPUTS:
;       Z   --  Interpolated array
; COMMON BLOCKS:
; NOTES:
;       Uses matrix multiplication #
; MODIFICATION HISTORY:
;       Han Uitenbroek, Sept. 1991
;-
srr= size(rr)
ndim= srr(0)
;
nx= srr(1)
i= fix(ix)  &  i= i < (nx-1)
ip= i+ 1  &  ip= ip < (nx-1)
dx= vector(ix- i)  &  dx1= 1.0- dx
;
IF (ndim eq 1) THEN $
  return, rr(i)*dx1+ rr(ip)*dx $
ELSE BEGIN
  ny= srr(2)
  j= fix(jy)  &  j= j < (ny-1)
  jp= j+ 1  &  jp= jp < (ny-1)
  dy= vector(jy- j)  &  dy1= 1.0- dy
  return, dsup(rr,i,j)*(dx1#dy1)+ dsup(rr,i,jp)*(dx1#dy)+ $
   dsup(rr,ip,j)*(dx#dy1)+ dsup(rr,ip,jp)*(dx#dy)
ENDELSE
;
end
