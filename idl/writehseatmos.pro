PRO writehseatmos, atmos, geometry,FILENAME=filename

  CM_TO_M = 1.0E-2
  KM_TO_M = 1.0E+3

  KG_TO_G = 1.0E+3
  
  openw, unit, filename, /GET_LUN

  printf, unit, atmos.ID

  printf, unit, $
   FORMAT='(1X,"MASS SCALE")'
  printf, unit, $
   FORMAT='("*",/"* lg g [cm s^-2]",/F6.2,/"*",/"* Ndep",/I4,/"*")', $
   4.44, n_elements(atmos.T)
  printf, unit, $
   FORMAT='("*lg column Mass   ", 5X,"Temperature",8X,"Ne",9X,"V",14X,"Vturb")'

  FOR k=0, n_elements(atmos.T)-1 DO BEGIN
    printf, unit, FORMAT='(E17.8, 4E15.6)', $
            ALOG10(geometry.cmass[k] * CM_TO_M^2 * KG_TO_G), $
            atmos.T[k], atmos.n_elec[k] * CM_TO_M^3 , $
            geometry.vz[k]/KM_TO_M, atmos.vturb[k]/KM_TO_M
  ENDFOR

  printf, unit, FORMAT='("*",/"* Hydrogen populations (LTE)")'
 
  printf, unit, FORMAT=$
   '("*",4X,"nh(1)",7X,"nh(2)",7X,"nh(3)",7X,"nh(4)",7X,"nh(5)",7X,"np")'

  FOR k=0, n_elements(atmos.T)-1 DO $
   printf, unit, FORMAT='(6E12.4)', atmos.nH[k, 0:5] * CM_TO_M^3

  free_lun, unit
END
