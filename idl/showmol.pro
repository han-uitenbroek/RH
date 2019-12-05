
  nHtot = total(atmos.nh, 2)

  CASE geometryType OF
    "ONE_D_PLANE": BEGIN
      hkm = geometry.height / 1.0E+03
      nrel = dblarr(geometry.Ndep, n_elements(molecules))
    END
    "SPHERICAL_SYMMETRIC": BEGIN
      hkm = geometry.r / 1.0E+03
      nrel = dblarr(geometry.Nradius, n_elements(molecules))
    END
  ENDCASE

  

  FOR nm=0, n_elements(molecules)-1 DO $
   nrel[*, nm] = *(molecules[nm].n_ptr) / nHtot
  ymin = min(nrel(where(nrel GT 0.0)), MAX=ymax)

  plot, hkm, nrel(*, 0), /YLOG, $
   YRANGE=[ymin, ymax], YSTYLE=1, XTITLE='Height [km]', $
   YTITLE='n!Dmolecule!N / n!S!DH!N!R!Utot!N', $
   TITLE=atmos.ID
  rhannotate, xann(0.02), yann(0.95), text=molecules[0].ID, MARKTYPE=0

  dcolor = 12
  dstyle = 4
  Nmolecule = n_elements(molecules)
  FOR nm=1, Nmolecule-1 DO BEGIN
    color = 200-dcolor*nm
    style = nm MOD dstyle
    oplot, hkm, nrel[*, nm], COLOR=color, LINE=style
    rhannotate, xann(0.02), yann(0.95-0.05*nm), TEXT=molecules[nm].ID, $
     MARKCOL=color, MARKTYPE=style
  ENDFOR

  color  = 200-dcolor*Nmolecule
  srtyle = Nmolecule MOD dstyle
  oplot, hkm, nHmin/nHtot, COLOR=color, LINE=style
  rhannotate, xann(0.02), yann(0.95-0.05*Nmolecule), TEXT='H!U-!N', $
     MARKCOL=color, MARKTYPE=style
END

