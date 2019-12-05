FUNCTION raytrace, geometry, mu, x_index, SHOW=show, XRAY=xray, ZRAY=zray

  mux = geometry.xmu[mu]
  muz = sqrt(1.0 - (geometry.xmu[mu]^2 + geometry.ymu[mu]^2))

  dztotal = geometry.z[0] - geometry.z[geometry.Nz-1]
  dxtotal = geometry.x[geometry.Nx-1] - geometry.x[0]
  xbottom = geometry.x[x_index] - (dztotal/muz) * mux

  x = geometry.x
  index = x_index
  Nxgrid = geometry.Nx
  IF (mux GT 0.0) THEN BEGIN
    WHILE (xbottom LT 0.0) DO BEGIN
      xbottom += dxtotal
      x = [geometry.x, x + dxtotal]
      index += geometry.Nx
      Nxgrid += geometry.Nx
    ENDWHILE
  ENDIF ELSE BEGIN
    WHILE (xbottom GT x[Nxgrid-1]) DO BEGIN
      x = [x, x[Nxgrid-1] + geometry.x]
      Nxgrid += geometry.Nx
    ENDWHILE
  ENDELSE

  ;; Define the structure of intersections

  intersects = {s: 0.0, type: "horizontal", x: 0.0,  z: 0.0, k: 0, l: 0}
  is = replicate(intersects, Nxgrid + geometry.Nz)

  ;; The first intersection is the point of origin

  is[0].type = "horizontal"
  is[0].s = 0.0
  is[0].x = x[index]
  is[0].z = geometry.z[0]
  is[0].k = 0
  is[0].l = index
  Nis = 1

  ;; Find all the intersections with horizontal planes

  FOR k=1, geometry.Nz-1 DO BEGIN
    dz = geometry.z[0] - geometry.z[k]

    is[Nis].type = "horizontal"
    is[Nis].s = dz / muz
    is[Nis].x = x[index] - mux * is[Nis].s
    is[Nis].z = geometry.z[k]
    is[Nis].k = k
    Nis++
  ENDFOR

  ;; Find all the intersections with vertical planes

  tabinv, x, xbottom, xeff
  IF (mux GT 0.0) THEN BEGIN
    lbottom = fix(xeff[0]+1)
    FOR l=index-1, lbottom, -1 DO BEGIN
      is[Nis].type = "vertical"

      is[Nis].s = (x[index] - x[l]) / mux
      is[Nis].x = x[l]
      is[Nis].l = l
      is[Nis].z = geometry.z[0] - muz * is[Nis].s
      Nis++
    ENDFOR
  ENDIF ELSE BEGIN
    lbottom = fix(xeff[0])
    FOR l=index+1, lbottom DO BEGIN
      is[Nis].type = "vertical"

      is[Nis].s = (x[index] - x[l]) / mux
      is[Nis].x = x[l]
      is[Nis].l = l
      is[Nis].z = geometry.z[0] - muz * is[Nis].s
      Nis++    
    ENDFOR
  ENDELSE

  ;; Shorten to to include only used intersections and sort on s

  is = is[0:Nis-1]
  is = is[sort(is.s)]

  ;; Eliminate intersections that coincide (i.e., have the same s)

  valid = intarr(Nis)
  valid[0] = 1
  FOR n=1, Nis-1 DO $
   IF (is[n].s GT is[n-1].s) THEN valid[n] = 1
  is = is[where(valid)]
  Nvalid = n_elements(is)

  vis = replicate({type: "vertical", k: [0, 0], l: [0, 0], $
                   coeff: [0.0, 0.0]}, Nvalid)

  FOR n=0, Nvalid-1 DO BEGIN
    vis[n].type = is[n].type
    CASE (vis[n].type) OF
      "vertical": BEGIN
        vis[n].l[0] = is[n].l MOD geometry.Nx
        tabinv, geometry.z, is[n].z, z_eff
        k0 = fix(z_eff) + 1
        vis[n].k = [k0, k0-1]
        frac = $
         (is[n].z - geometry.z[k0]) / (geometry.z[k0-1] - geometry.z[k0])
      END
      "horizontal": BEGIN
        vis[n].k = is[n].k
        tabinv, x, is[n].x, x_eff
        l0 = fix(x_eff)
        frac = (l0 EQ Nxgrid-1) ? 1.0 : (is[n].x - x[l0]) / (x[l0+1] - x[l0])
        vis[n].l = ([l0, l0+1] MOD geometry.Nx)
      END
    ENDCASE

    vis[n].coeff = [1.0 - frac, frac]
  ENDFOR

  IF (keyword_set(SHOW)) THEN BEGIN
    KM_TO_M = 1.0E3

    plot, [geometry.x[0], geometry.x[geometry.Nx-1]] / KM_TO_M, $
          [geometry.z[geometry.Nz-1], geometry.z[0]] / KM_TO_M, /NODATA, $
     XTITLE='x [km]', YTITLE='z [km]'

    FOR k=0, geometry.Nz-1 DO $
     oplot, [geometry.x[0], geometry.x[geometry.Nx-1]] / KM_TO_M, $
     [geometry.z[k], geometry.z[k]] / KM_TO_M, LINESTYLE=1

    FOR l=0, geometry.Nx-1 DO $
     oplot, [geometry.x[l], geometry.x[l]] / KM_TO_M, $
     [geometry.z[geometry.Nz-1], geometry.z[0]] / KM_TO_M, LINESTYLE=1
  ENDIF

  xray = fltarr(Nvalid)
  zray = xray
  FOR n=0, Nvalid-1 DO BEGIN
    CASE (vis[n].type) OF
      "vertical": BEGIN
        xray[n] = geometry.x[vis[n].l[0]]
        zray[n] = vis[n].coeff[0]*geometry.z[vis[n].k[0]] + $
         vis[n].coeff[1]*geometry.z[vis[n].k[1]]
      END
      "horizontal": BEGIN
        xray[n] = vis[n].coeff[0]*geometry.x[vis[n].l[0]] + $
         vis[n].coeff[1]*geometry.x[vis[n].l[1]]
        zray[n] = geometry.z[vis[n].k[0]]
      END
    ENDCASE
  ENDFOR

  IF (keyword_set(SHOW)) THEN oplot, xray/KM_TO_M, zray/KM_TO_M, THICK=3

  return, vis
END
