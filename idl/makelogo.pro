; -------- begin -------------------------- makeLogo_1D.pro ---------- ;

PRO makeLogo_1D, fileName

  IF (n_elements(fileName) LE 0) THEN fileName = '~/idl/lib/rh/logo_1D.gif'

  window, 2, XSIZE=450, YSIZE=100
  plot, [0,1], XMARGIN=[2,2], YMARGIN=[2,2], CHARSIZE=0.3
  backgroundColor = 225B  &  textColor = 75B  &  rayColor = 208B
  device, FONT='-linotype-*-*-o-*-*-40-200-*-*-*-*-*-*'
  !P.BACKGROUND = backgroundColor

  tv, bytarr(450, 100) + backgroundColor
  plots, [0.3, 0.6], [0.0, 0.0], color=textColor
  plots, [0.3, 0.6], [0.85, 0.85], color=textColor
  plots, [0.35, 0.35], [0.0, 1.0], color=textColor
  xmu = [0.99, 0.98, 0.95]
  FOR mu=0, n_elements(xmu)-1 DO BEGIN
    sinTheta = sqrt(1.0 - xmu(mu)^2)
    tanTheta = sinTheta / xmu(mu)
    arrow, 0.35, 0.0, 0.35 + tanTheta, 1.0, COLOR=rayColor, /DATA, THICK=2
  ENDFOR

  xyouts, /DEVICE, 25, 40, 'Analyze', FONT=1, COLOR=textColor
  xyouts, /DEVICE, 300, 20, '1-D', FONT=1, COLOR=textColor
  logo = tvrd()
  write_gif, fileName, logo
  !P.BACKGROUND = 0B

  tv, logo
END
; -------- end ---------------------------- makeLogo_1D.pro ---------- ;

; -------- begin -------------------------- makeLogo_sphere.pro ------ ;

PRO makeLogo_sphere, fileName

  IF (n_elements(fileName) LE 0) THEN $
   fileName = '~/idl/lib/rh/logo_sphere.gif'

  window, 2, XSIZE=450, YSIZE=100
  plot, [0,1], XMARGIN=[2,2], YMARGIN=[2,2], CHARSIZE=0.3
  backgroundColor = 225B  &  textColor = 75B  &  rayColor = 208B
  device, FONT='-linotype-*-*-o-*-*-40-200-*-*-*-*-*-*'
  !P.BACKGROUND = backgroundColor

  tv, bytarr(450, 100) + backgroundColor
  plots, [0.3, 0.6], [0.0, 0.0], color=textColor
  plots, [0.35, 0.35], [0.0, 1.2], color=textColor

  x = 0.3*sin(0.5*!PI*findgen(50)/49.0)
  y = cos(0.5*!PI*findgen(50)/49.0)
  radius = [0.2, 0.4, 0.7, 0.8]
  FOR k= 0, n_elements(radius)-1 DO $
    plots, 0.35 + radius(k)*x, radius(k)*y, COLOR=textColor
  FOR k= 0, n_elements(radius)-1 DO $
   arrow, 0.35, radius(k), 0.65, radius(k), COLOR=rayColor, /DATA, THICK=2

  xyouts, /DEVICE, 25, 40, 'Analyze', FONT=1, COLOR=textColor
  xyouts, /DEVICE, 300, 20, 'SPHERE', FONT=1, COLOR=textColor
  logo = tvrd()
  write_gif, fileName, logo
  !P.BACKGROUND = 0B

  tv, logo
END
; -------- end ---------------------------- makeLogo_sphere.pro ------ ;

