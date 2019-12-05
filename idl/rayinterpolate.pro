FUNCTION rayinterpolate, f, vis

  s = size(f)
  IF (s[0] NE 2) THEN BEGIN
    print, "RAYINTERPOLATE: f should be a 2-D function"
    return, 0
  ENDIF

  Nvalid = n_elements(vis)
  fray = dblarr(Nvalid)

  FOR n=0, Nvalid-1 DO BEGIN
    CASE (vis[n].type) OF
      "vertical": BEGIN
        fray[n] = vis[n].coeff[0]*f[vis[n].l[0], vis[n].k[0]] + $
         vis[n].coeff[1]*f[vis[n].l[0], vis[n].k[1]]
      END
      "horizontal": BEGIN
        fray[n] = vis[n].coeff[0]*f[vis[n].l[0], vis[n].k[0]] + $
         vis[n].coeff[1]*f[vis[n].l[1], vis[n].k[0]]
      END
    ENDCASE
  ENDFOR

  return, fray
END
