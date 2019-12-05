pro PSTV, image, x, y, order=order

  common SCREEN_COMMON, ScreenSize, ScaleFactor

  if (n_elements(x)) eq 0 then xpos = 0
  if (n_elements(y)) eq 0 then ypos = 0

  if(keyword_set(order)) then order = 1  else order = 0

  if (!D.NAME eq 'PS') then begin
    xps = float(x)/ScreenSize[0]*!D.X_SIZE
    yps = float(y)/ScreenSize[1]*!D.Y_SIZE
    sz  = size(image)
    xsz = float(sz[1])/ScreenSize[0]*!D.X_SIZE
    ysz = float(sz[2])/ScreenSize[1]*!D.Y_SIZE
    tv, image, xps, yps, XSIZE=xsz, YSIZE=ysz, /DEVICE, ORDER=order
  endif else $
   tv, image, x, y, ORDER=order

  return
end
