PRO labsneckatlas, Fatlas, latlas

  ;; Returns selected wavelength range from the Labs & Neckel
  ;; absolute flux atlas.
  ;;
  ;; Parameters:
  ;;   Input:    lambdaMin  -- Minimum wavelength [nm]
  ;;             lambdaMax  -- Maximum wavelength [nm]

  ;;   Output:   Fatlas     -- Float array of flux values [J s^-1 m^-2 Hz^-1]
  ;;             latlas     -- Float array of corresponding wavelengths [nm]

  MICROWATT = 1.0E-06
  CM_TO_M   = 1.0E-02
  NM_TO_M   = 1.0E-09
  CLIGHT    = 2.99792458E+08

  if (n_params(0) lt 1) then begin
    print, "Usage: labsneckatlas, Fatlas, [latlas]"
    print, "NOTE:  Wavelengths are in nanometers!!"
    return
  ENDIF
 
  atlasfile = getenv('RH_ATLAS_PATH') + '/labsneck.atlas'
  openr, unit, atlasfile, /GET_LUN

  Nln = 0L
  readf, unit, Nln
  lnFlux = fltarr(3, Nln)
  readf, unit, lnFlux
  free_lun, unit

  lnFlux = transpose(lnFlux)
  latlas = lnFlux(*, 0)
  Fatlas = MICROWATT/CLIGHT * (latlas*NM_TO_M) * (latlas*lnFlux(*, 1))

  ;; Transparency correction

  Fatlas = Fatlas * (1.0 - lnFlux(*, 2)) * 1.0E+08
END
