PRO opacandsource, WAVENO=waveNo, RAY=ray, chi, S, JNU=Jnu

@files.common
@geometry.common
@atmos.common
@opacity.common
@spectrum.common

  IF (n_params() LT 2) THEN BEGIN
    print, "Usage: OPACANDSOURCE, WAVENO=waveNo, RAY=ray, chi, S, JNU=Jnu"
    return
  ENDIF


  metalFile      = "metals.out"
  moleculeFile   = "molecules.out"
  JFile          = "J.dat"
  opacFile       = "opacity.out"
  backgroundFile = "background.dat"

  result = readInput('input.out')
  result = readgeometry('geometry.out')
  result = readatmos('atmos.out')
  result = readspectrum('spectrum.out')

  result = openJ(JFile)
  result = openOpacity(opacFile)

  readJ, waveNo
  readOpacity, waveNo, ray

  chi = chi_as + chi_c
  S   = (eta_as + eta_c + J*scatt) / chi

  Jnu = J

  free_lun, Junit
  free_lun, opacunit
  free_lun, backgroundunit
END
