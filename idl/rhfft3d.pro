
;;  Nx = 256L  &  Ny = 256L
  Nx = 16L  &  Ny = 16L
  xoffset = double(56.7639)
  yoffset = double(27.5679)

  x = findgen(Nx)/(Nx-1)
  y = findgen(Ny)/(Ny-1)
  data   = double((1.0 - 4*(x-0.5)^2) # (1.0 - 4*(y-0.5)^2))
  data_s = dblarr(Nx, Ny)

  spawn, /NOSHELL, UNIT=pipe, './rhfft3d_test'
  writeu, pipe, Nx, Ny, xoffset, yoffset, data

  readu, pipe, data_s

  pswindow
  panel, scaleimg(data_s, 256, 256)
END
