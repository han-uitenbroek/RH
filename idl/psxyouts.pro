pro PSXYOUTS, xpos, ypos, string, alignment=alignment, charsize=charsize

  common SCREEN_COMMON, ScreenSize, ScaleFactor

  if (not keyword_set(alignment)) then alignment = 0.0
  if (not keyword_set(charsize))  then charsize  = 1.0

  xd = xpos*ScaleFactor(0)
  yd = ypos*ScaleFactor(1)

  xyouts, xd, yd, /DEVICE, string, alignment=alignment, charsize=charsize

end
