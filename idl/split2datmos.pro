PRO split2datmos, twod_in, oned_out, x_index

@ATMOS

  CM_TO_M = 1.0E-2

  IF (n_params(0) LT 3) THEN BEGIN
    print, "Usage: split2datmos, twod_in, oned_out, x_index"
    return
  ENDIF

  atmos = read2datmos(twod_in)
  IF (x_index LT 0  OR  x_index GE atmos.Nx) THEN BEGIN
    print, FORMAT='("Invalid index for x! Should be between 0 and ", I3)', $
     atmos.Nx-1
    return
  ENDIF

  openw, unit, oned_out, /GET_LUN

  printf, unit, $
   FORMAT='(2X, "Column #", I2, ", extracted from: ", A)', $
   x_index, twod_in
  printf, unit, $
   FORMAT='(2X,"Height scale")'
  printf, unit, $
   FORMAT='("*",/"* lg g",/F6.2,/"*",/"* Ndep",/I4,/"*")', $
   4.44, atmos.Nz
  printf, unit, $
   FORMAT='("*height [km]   ", 5X,"Temperature",8X,"Ne",9X,"V",14X,"Vturb")'

  FOR k=0, atmos.Nz-1 DO $
   printf, unit, FORMAT='(E17.8, 4E15.6)', $
   atmos.z(k), atmos.T(x_index, k), atmos.n_elec(x_index, k) * CM_TO_M^3, $
   atmos.vz(x_index, k), atmos.vturb(x_index, k)

  IF (atmos.NHydr EQ 1) THEN $
   printf, unit, FORMAT='("*",/"* Hydrogen populations (LTE)")' $
  ELSE $
   printf, unit, FORMAT='("*",/"* Hydrogen populations")' 

  printf, unit, FORMAT=$
   '("*",4X,"nh(1)",7X,"nh(2)",7X,"nh(3)",7X,"nh(4)",7X,"nh(5)",7X,"np")'

  nH = dblarr(atmos.Nz, 6)
  FOR i=0,atmos.NHydr-1 DO nH(*, i) = atmos.nH(x_index, *, i)

  FOR k=0, atmos.Nz-1 DO $
   printf, unit, FORMAT='(6E12.4)', nH(k, *) * CM_TO_M^3

  free_lun, unit
END

