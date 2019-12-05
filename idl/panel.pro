PRO panel, image, x, y, XPOS=xpos, YPOS=ypos, $
           LEFT=left, SCALETEXT=ScaleText, NOSCALE=noscale, $
           XTITLE=xtitle, YTITLE=ytitle, TITLE=title, NOIMAGE=noimage, $
           TOP_BAR=top_bar, ORDER=order, ORTHOSCOPIC=orthoscopic, $
           REVERSE_X=reverse_x, REVERSE_Y=reverse_y, $
           TABLE_SIZE=table_size, TABLE_OFFSET=table_offset, ZLOG=zlog, $
           SCALECHARSIZE=scaleCharsize, ERASE=erase, TICKLEN=ticklen, $
           SCALEWIDTH=scalewidth, _EXTRA=extra_keywords

;+
; NAME:
;	PANEL
;
; PURPOSE:
;	This procedure displays a two-dimensional image in a plot
;       window with x- and y-axes. It adds a bar showing the relation
;       between values in the image and their assigned (false) color.
;
; CATEGORY:
;	Image display
;
; CALLING SEQUENCE:
;       PANEL, Image, X, Y
;
; INPUTS:
;	Image:  The two-dimensional image array.
;
; OPTIONAL INPUTS:
;	X:      Scale for x-axis.
;	Y:      Scale for y-axis.
;	
; KEYWORD PARAMETERS:
;       XPOS:       X-coordinate (in device coordinates) of lower right corner
;                   of plot window
;
;       YPOS:       Y-coordinate (in device coordinates) of lower right corner
;                    of plot window
;
;       LEFT:       Set this keyword to place the color-scale bar to the left
;                    of the image. Default is to the right.
;
;       TITLE:      Main title of plot window
;
;       XTITLE:     Title for x-axis
;
;       YTITLE:     Title for y-axis
;
;       SCALETEXT:  Title for color-scale bar. By default it is placed along
;                    the vertical axis.
;
;       NOIMAGE:    Set this image to draw only the axes.
;       TOP:        Set this keyword to draw SCALETEXT on top of the
;                    color-scale bar.
;
;       ORDER:      Set this keyword to flip image vertically.
;
;       ORTHOSCOPIC: Use this keyword when the x- and/or y-axis are
;                     not equally spaced.
;
;      _EXTRA:      Extra keywords are passed to the plot command that
;                   draws the axes.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;       Uses routines PSTV and PSPOSITION. These require the image scale
;       to be set prior to calling by routine PSWINDOW
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Wed Oct  4 22:02:59 2006 --
;-

  IF (n_params(0) LT 1) THEN BEGIN
    print, "Usage: panel, image, x, y, XPOS=xpos, YPOS=ypos, "
    print, "LEFT=left, SCALETEXT=ScaleText, NOSCALE=noscale, "
    print, "XTITLE=xtitle, YTITLE=ytitle, TITLE=title, NOIMAGE=noimage, "
    print, "TOP_BAR=top_bar, ORDER=order, ORTHOSCOPIC=orthoscopic, "
    print, "ZLOG=zlog, SCALECHARSIZE=scaleCharsize, _EXTRA=extra_keywords"

    return
  ENDIF

  IF (keyword_set(ERASE)) THEN erase

  sz = size(image)
  Nx = sz(1)  &  Ny = sz(2)
  BlankTickname = strarr(30)+' '
  dx_scale = 20
  dy_scale = 27

  IF (NOT keyword_set(XPOS))  THEN xpos  = 25
  IF (NOT keyword_set(YPOS))  THEN ypos  = 25
  IF (n_elements(y) EQ 0)     THEN y     = [0, Ny]
  IF (n_elements(x) EQ 0)     THEN x     = [0, Nx]
  IF (NOT keyword_set(TITLE)) THEN title = ''
  IF (keyword_set(ORDER))     THEN order = 1  ELSE  order = 0
  IF (NOT keyword_set(ZLOG))  THEN zlog  = 0  ELSE  zlog = 1
  IF (NOT keyword_set(SCALECHARSIZE)) THEN scaleCharsize = 1.0
  IF (NOT keyword_set(TABLE_SIZE)) THEN table_size = !D.TABLE_SIZE
  IF (n_elements(TABLE_OFFSET) LE 0) THEN table_offset = 0B
  IF (NOT keyword_set(TICKLEN))  THEN ticklen = -0.02
  IF (NOT keyword_set(SCALEWIDTH)) THEN ScaleWidth = 15

  IF (keyword_set(TOP_BAR)) THEN BEGIN
    xpi = xpos
    xps = xpos
    ypi = ypos
    yps = ypos + Ny + dy_scale
  ENDIF ELSE BEGIN
    ypi = ypos
    yps = ypos
    IF (keyword_set(LEFT)) THEN BEGIN
      xps = xpos
      xpi = xpos + ScaleWidth + dx_scale
    ENDIF ELSE BEGIN
      xpi = xpos
      xps = xpos + Nx + dx_scale
    ENDELSE
  ENDELSE

  IF (ZLOG) THEN image = alog10(image)

  IF (keyword_set(ORTHOSCOPIC)) THEN BEGIN
    IF (n_params(0) LT 3) THEN BEGIN
      print, "panel: keyword ORTHOSCOPIC needs x and y grid values"
      return
    ENDIF
    Nxg = n_elements(x)
    tabinv, x, x(0) + (x(Nxg-1) - x(0)) * findgen(Nxg)/(Nxg-1), xeff
    Nyg = n_elements(y)
    tabinv, y, y(0) + (y(Nyg-1) - y(0)) * findgen(Nyg)/(Nyg-1), yeff

    xscale = Nx / float(Nxg-1)  &  yscale =  Ny / float(Nyg-1)
    image = interpolate(image, xscale * scaleline(xeff, Nx), $
                        yscale * scaleline(yeff, Ny), /GRID)
  ENDIF

  ;; displays with maximum detail, so scale between min and max

  imin = min(image, max=imax)
  IF (imin EQ imax) THEN BEGIN
    byte_image = bytarr(Nx, Ny) + 127
    IF (imin EQ 0.0) THEN BEGIN
      imin = -0.5  &  imax = +0.5
    ENDIF ELSE BEGIN
      tmp = imin
      imin = 0.9*tmp  &  imax = 1.1*tmp
    ENDELSE
  ENDIF ELSE BEGIN
    byte_image = table_offset + $
     byte((table_size - 1) * (image - imin)/float(imax - imin))
  ENDELSE

  IF (keyword_set(TOP_BAR)) THEN BEGIN
    scale = rebin(reform(table_offset + $
                         byte(indgen(Nx)/float(Nx-1)*(table_size-1)), $
                         Nx, 1), Nx, ScaleWidth)
  ENDIF ELSE BEGIN
    scale = rebin(reform(table_offset + $
                         byte(indgen(Ny)/float(Ny-1)*(table_size-1)), $
                         1, Ny), ScaleWidth, Ny)
  ENDELSE

  IF (ZLOG) THEN BEGIN
    imin = 10.0^imin
    imax = 10.0^imax
  ENDIF

  IF (NOT keyword_set(XTITLE)) THEN BEGIN
    xTickname = BlankTickname
    xtitle = ' '
  ENDIF ELSE BEGIN
    xTickName=''
    xtitle = xtitle
  ENDELSE
  IF (NOT keyword_set(YTITLE)) THEN BEGIN
    yTickname = BlankTickname
    ytitle = ' '
  ENDIF ELSE BEGIN
    yTickName=''
    ytitle = ytitle
  ENDELSE

  IF (NOT keyword_set(NOSCALE)) THEN BEGIN 
    IF (NOT keyword_set(SCALETEXT)) THEN ScaleText = ''

    IF (NOT keyword_set(NOIMAGE)) THEN pstv, scale, xps, yps

    IF (keyword_set(TOP_BAR)) THEN BEGIN 
      position = psposition([xps, yps, xps+Nx, yps+ScaleWidth])
      xstyle = 5
      ystyle = 1
      xticklen = 8.0*ticklen  &  yticklen = 0.0001
    ENDIF ELSE BEGIN
      position = psposition([xps, yps, xps+ScaleWidth, yps+Ny])
      xstyle = 1
      ystyle = 5
      yticklen = 12.0*ticklen  &  xticklen = 0.0001
    ENDELSE

    IF (keyword_set(TOP_BAR)) THEN scaletitle=ScaleText  ELSE  scaletitle=' '
    plot, /NODATA, /NOERASE, /DEVICE, YLOG=zlog, $
          imin + (imax - imin)*findgen(Ny)/(Ny-1), $
          imin + (imax - imin)*findgen(Ny)/(Ny-1), $
          YTICKLEN=yticklen, XTICKLEN=xticklen, $
          XSTYLE=xstyle, YSTYLE=ystyle, $
          XTICKNAME=BlankTickname, YTICKNAME=BlankTickname, $
          POSITION=position, CHARSIZE=scaleCharsize, $
          TITLE=scaletitle, _EXTRA=extra_keywords

    IF (keyword_set(TOP_BAR)) THEN BEGIN
        axis, XAXIS=0, XSTYLE=1, XTICKLEN=xticklen, XTITLE="", $
              XLOG=zlog, CHARSIZE=scaleCharsize, _EXTRA=extra_keywords
        axis, XAXIS=1, XSTYLE=1, XTICKLEN=0, $
              XLOG=zlog, XTICKNAME=BlankTickname, _EXTRA=extra_keywords
    ENDIF ELSE BEGIN
      IF (keyword_set(LEFT)) THEN BEGIN
        axis, YAXIS=0, YSTYLE=1, YTICKLEN=yticklen, $
              YTITLE=ScaleText, CHARSIZE=scaleCharsize, $
              YLOG=zlog, _EXTRA=extra_keywords
        axis, YAXIS=1, YSTYLE=1, YTICKLEN=xticklen, YLOG=zlog, $
              YTICKNAME=BlankTickname, _EXTRA=extra_keywords
      ENDIF ELSE BEGIN
        axis, YAXIS=0, YSTYLE=1, YTICKLEN=xticklen, $
              YTICKNAME=BlankTickname, YLOG=zlog, _EXTRA=extra_keywords
        axis, YAXIS=1, YTITLE=ScaleText, YSTYLE=1, YTICKLEN=yticklen, $
              YLOG=zlog, CHARSIZE=scaleCharsize, _EXTRA=extra_keywords
      ENDELSE
    ENDELSE
  ENDIF ELSE BEGIN
    xstyle = 1
    ystyle = 1
  ENDELSE

  IF (keyword_set(REVERSE_X)) THEN byte_image = reverse(byte_image, 1)
  IF (keyword_set(REVERSE_Y)) THEN byte_image = reverse(byte_image, 2)

  IF (NOT keyword_set(noimage)) THEN pstv, byte_image, xpi, ypi, ORDER=order
  plot, /NODATA, /NOERASE, /DEVICE, $
   [x(0), x(n_elements(x)-1)], [y(0), y(n_elements(y)-1)], TITLE=title, $
   XTICKNAME=xTickName, YTICKNAME=yTickName, $
   XTITLE=xtitle, XTICKLEN=ticklen, XSTYLE=xstyle, $
   YTITLE=ytitle, YTICKLEN=ticklen, YSTYLE=ystyle, $
   POSITION=psposition([xpi, ypi, xpi+Nx, ypi+Ny]), _EXTRA=extra_keywords

  IF (NOT keyword_set(NOSCALE)) THEN BEGIN
    IF (keyword_set(TOP_BAR)) THEN BEGIN
      axis, XAXIS=0, XTITLE=xtitle, XSTYLE=1, TICKLEN=ticklen, $
            XTICKNAME=xTickname, _EXTRA=extra_keywords
      axis, XAXIS=1, XSTYLE=1, TICKLEN=ticklen, XTICKNAME=BlankTickname, $
            _EXTRA=extra_keywords
    ENDIF ELSE BEGIN
      IF (keyword_set(LEFT)) THEN BEGIN
        axis, YAXIS=0, YSTYLE=1, TICKLEN=ticklen, YTICKNAME=BlankTickname, $
              _EXTRA=extra_keywords
        axis, YAXIS=1, YTITLE=ytitle, YSTYLE=1, TICKLEN=ticklen, $
              YTICKNAME=yTickname, _EXTRA=extra_keywords
      ENDIF ELSE BEGIN
        axis, YAXIS=0, YTITLE=ytitle, YSTYLE=1, TICKLEN=ticklen, $
              YTICKNAME=yTickname, _EXTRA=extra_keywords
        axis, YAXIS=1, YSTYLE=1, TICKLEN=ticklen, YTICKNAME=BlankTickname, $
              _EXTRA=extra_keywords
      ENDELSE
    ENDELSE
  ENDIF
END
