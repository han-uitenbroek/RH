FUNCTION read_rlk_pf_table, pf_input, TPF=Tpf

  pf_struct = {pti: 0L,  Nstage: 0L,  pf: ptr_new(),  ionpot: ptr_new()}

  Nelem = 99
  pf_elem = replicate(pf_struct, Nelem)

  Npf = 0L
  openr, lun, /XDR, /GET_LUN, pf_input
  readu, lun, Npf
  Tpf = dblarr(Npf)
  readu, lun, Tpf

  pti    = 0L
  Nstage = 0L
  FOR n=0, Nelem-1 DO BEGIN
    readu, lun, pti, Nstage
    pf     = dblarr(Npf, Nstage)
    ionpot = dblarr(Nstage)
    readu, lun, pf, ionpot

    pf_elem[n].pti    = pti
    pf_elem[n].Nstage = Nstage
    pf_elem[n].pf     = ptr_new(pf)
    pf_elem[n].ionpot = ptr_new(ionpot)
  ENDFOR

  free_lun, lun

  return, pf_elem
END

    
