pro PSWINDOW, WindowNumber, xsize=xsize, ysize=ysize, nopop=nopop, $
               title=title, free=free

  common SCREEN_COMMON, ScreenSize, ScaleFactor

  if (n_params(0) lt 1) then WindowNumber = 0

  if (not keyword_set(xsize)) then xsize = 640
  if (not keyword_set(ysize)) then ysize = 512
  if (not keyword_set(title)) then $
   title = 'PSWindow' + string(WindowNumber, FORMAT='(I2)')

  if (!D.NAME eq 'X' and (not keyword_set(nopop))) then begin
    if (keyword_set(free)) then begin
      window, xsize=xsize, ysize=ysize, title=title, /FREE
    endif else begin
      window, WindowNumber, xsize=xsize, ysize=ysize, title=title
      WindowNumber = !D.WINDOW
    endelse
  endif

  ScreenSize  = [xsize, ysize]
  ScaleFactor = [!D.X_SIZE, !D.Y_SIZE]/float(ScreenSize)

  return
end
