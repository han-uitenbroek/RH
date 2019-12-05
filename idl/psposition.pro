function PSPOSITION, position

  common SCREEN_COMMON, ScreenSize, ScaleFactor

  if n_elements(ScaleFactor) eq 0 then $
    message, 'scaling factor not established, call swindow first'

  return, position*[ScaleFactor, ScaleFactor]

end
