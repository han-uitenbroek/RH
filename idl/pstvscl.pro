pro PSTVSCL, image, x, y

  common SCREEN_COMMON, ScreenSize, ScaleFactor

  if (n_elements(x)) eq 0 then xpos = 0
  if (n_elements(y)) eq 0 then ypos = 0

  if (!D.NAME eq 'PS') then begin
    xps = float(x)/ScreenSize(0)*!D.X_SIZE
    yps = float(y)/ScreenSize(1)*!D.Y_SIZE
    sz  = size(image)
    xsz = float(sz(1))/ScreenSize(0)*!D.X_SIZE
    ysz = float(sz(2))/ScreenSize(1)*!D.Y_SIZE
    tvscl, image, xps, yps, xsize=xsz, ysize=ysz, /DEVICE
  endif else $
    tvscl, image, x, y

  return
end
