PRO dump1datmos, geometry, atmos, filename, $
                 MODELNAME=modelname, COMMENT=comment, GRAVITY=gravity

  CM_TO_M = 1.0E-02
  G_TO_KG = 1.0E-03

  IF (n_params() LT 3) THEN BEGIN
    message, "Usage: DUMP1DATMOS, geometry, atmos, filename", /INFORMATIONAL
    return
  ENDIF

  openw, atmos_unit, /GET_LUN, filename

  IF (NOT keyword_set(COMMENT)) THEN comment = 'Model saved from data'

  IF (NOT keyword_set(MODELNAME)) THEN modelname = 'XXX'
  Ndep_string = strtrim(string(geometry.Ndep), 2)
  modelname = modelname + '_' + Ndep_string

  IF (NOT keyword_set(GRAVITY)) THEN gravity = 10.0^4.44


  format1 = '("* ", A, /"*", /2X, A, /"  Mass scale", /"*")'
  format2 = '("* log g [cm s^-2]", /F6.2, /"*")'
  format3 = '("* Ndep", /I4, /"*")'
  format4 = $
     '("*lg column Mass", 5X, "Temperature", 8X, "Ne", 9X, "V", 14X, "Vturb")'

  printf, atmos_unit, FORMAT=format1, comment, modelname
  printf, atmos_unit, FORMAT=format2, alog10(gravity)
  printf, atmos_unit, FORMAT=format3, geometry.Ndep

  printf, atmos_unit, FORMAT=format4
  FOR k=0, geometry.Ndep-1 DO $
   printf, atmos_unit, FORMAT='(E17.8, 4E15.6)', $
   alog10(geometry.cmass[k] * (CM_TO_M^2 / G_TO_KG)), $
   atmos.T[k], atmos.n_elec[k] * CM_TO_M^3, geometry.vz[k], atmos.vturb[k]

  IF (atmos.NHydr EQ 1) THEN $
   printf, atmos_unit, FORMAT='("*", /"* Hydrogen populations (LTE)")' $
  ELSE $
   printf, atmos_unit, FORMAT='("*", /"* Hydrogen populations")' 
  printf, atmos_unit, FORMAT=$
   '("*", 4X, "nh(1)", 7X, "nh(2)", 7X, "nh(3)", 7X, "nh(4)"' + $
   ', 7X, "nh(5)", 7X, "np")'
  FOR k=0, geometry.Ndep-1 DO $
   printf, atmos_unit, FORMAT='(6E12.4)', atmos.nH[k, *] * CM_TO_M^3

  free_lun, atmos_unit
END
