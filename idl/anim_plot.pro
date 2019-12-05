PRO EHANDLER, ev
  widget_control, /DESTROY, ev.TOP
END

PRO anim_plot, x, movie, $
               TRACK=track, CYCLE=cycle, DISPLAY_MEAN=display_mean, $
               PNG_NAME=png_name, _EXTRA=extra_keywords

  COMMON colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  movieSize = size(movie)
  IF (movieSize[0] NE 2) THEN BEGIN
    print, 'Not a two-dimensional array'
    return
  ENDIF

  ;; --- Some cosmetic parameters --

  xwindowSize = 500
  ywindowSize = 380
  meanColor   = 16B


  moviemin = min(movie, MAX=moviemax)
  mrange   = moviemax - moviemin 
  yrange   = [moviemin, moviemax]

  Ntime = movieSize[2]

  baseWidget    = widget_base(TITLE='Animate plot')
  animateWidget = cw_animate(baseWidget, xwindowSize, ywindowSize, Ntime, $
			     TRACK=keyword_set(TRACK), $
			     CYCLE=keyword_set(CYCLE))
  widget_control, /REALIZE, baseWidget

  window, /FREE, XSIZE=xwindowSize, YSIZE=ywindowSize, /PIXMAP
  pixMap = !D.WINDOW

  IF (keyword_set(DISPLAY_MEAN)) THEN movie_mean = avg(movie, 1)

  FOR n_iter=0, Ntime-1 DO BEGIN
    plot, x, movie[*, n_iter], $
     YRANGE=yrange, _EXTRA=extra_keywords
    IF (keyword_set(DISPLAY_MEAN)) THEN $
     oplot, x, movie_mean, COLOR=meanColor, THICK=2
    cw_animate_load, animateWidget, FRAME=n_iter, WINDOW=pixMap

    IF (keyword_set(PNG_NAME)) THEN BEGIN
      wset, pixMap
      filename = png_name + '_' + $
       string(n_iter, FORMAT=((n_iter LT 10) ? '(I1)' : '(I2)')) + '.png'
      write_png, filename, tvrd(), r_curr, g_curr, b_curr
    ENDIF
  ENDFOR
  wdelete, pixMap

  cw_animate_run, animateWidget
  xmanager, "Xanim", baseWidget, EVENT_HANDLER="EHANDLER"
END
