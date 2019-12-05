PRO KPIRatlas, lambdamin, lambdamax, Iatlas, latlas

  atlasdir = getenv('RH_ATLAS_PATH') + '/KPIR/'
  file = 'KP_IR.atlas'
  openr, lun, atlasdir + file, /GET_LUN, /XDR

  Nwave = 0L
  readu, lun, Nwave
  Iatlas = fltarr(Nwave, /NOZERO)
  latlas = fltarr(Nwave, /NOZERO)

  readu, lun, latlas, Iatlas
  free_lun, lun

  tabinv, latlas, [lambdamin, lambdamax], lambda_eff

  latlas = latlas[long(lambda_eff[0]):long(lambda_eff[1])+1]
  Iatlas = Iatlas[long(lambda_eff[0]):long(lambda_eff[1])+1]
END
