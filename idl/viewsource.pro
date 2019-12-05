; -------- file: -------------------------- viewsource.pro ----------- ;

; -------- begin -------------------------- setRatio.pro ------------- ;

FUNCTION setRatio, stash, ratio
  widget_control, stash, GET_UVALUE=state

  state.ratio = ratio
  widget_control, stash, SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- setRatio.pro ------------- ;

; -------- begin -------------------------- XViewSource_Event.pro ---- ;

PRO XViewSource_Event, Event

@geometry.common
@files.common

  COMMON screen_common, screenSize, scaleFactor

  Off = 0  &  On = 1
  screenSize  = fix([550, 400])

  ;; --- Main event handler --                          -------------- ;

  widget_control, Event.id, GET_UVALUE=Action
  stash = widget_info(Event.top, /CHILD)
  widget_control, stash, GET_UVALUE=state

  CASE Action OF
    'QUIT': widget_control, Event.top, /DESTROY

    'PRINT': BEGIN
      filename = 'viewSource-' + timeStamp() + '.ps'
      IF ((geometryType EQ "TWO_D_PLANE" OR $
           geometryType EQ "THREE_D_PLANE") AND state.displayType GT 0) THEN $
       font = -1 ELSE font = 0
      psopen, FILENAME=filename, /COLOR, font=font
      ScaleFactor = [!D.X_SIZE, !D.Y_SIZE] / float(ScreenSize)
      drawS, state
      psclose
      ScaleFactor = [1.0, 1.0]
      r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
    END

    'PNG': BEGIN
      widget_control, state.drawWidget, GET_VALUE=WindowNo
      rhwritepng, WindowNo, 'viewSource-'
    END

    'ORIENT': orient, state, EVENT_HANDLER='XViewSource_Event', $
     GROUP_LEADER=Event.top, TITLE='Orientation source function'
  
    'ORIENT_UPDATE': BEGIN
      OrientationUpdate, state
      drawS, state
    END

    'ORIENT_RESET': BEGIN
      OrientationUpdate, state, /SET
      drawS, state
    END

    'TRACK_SURFACE': drawS, state, /TRACK

    'SLIDER': drawS, setLambdaNo(stash, Event.value)

    'RAYMENU': drawS, setRayNo(stash, Event.value)

    'NEWOPACITYFILE': BEGIN
      IF (n_elements(opacunit) GT 0) THEN free_lun, opacunit
      IF (openOpacity("*.out")) THEN drawS, setLambdaNo(stash, 0)
    END

    'PANEL':     IF (Event.select) THEN drawS, setDisplayType(stash, 0)
    'WIREGRID':  IF (Event.select) THEN drawS, setDisplayType(stash, 1)
    'SHADESURF': IF (Event.select) THEN drawS, setDisplayType(stash, 2)
 
    'RATIO_OFF': IF (Event.select) THEN drawS, setRatio(stash, 0)
    'RATIO1_ON': IF (Event.select) THEN drawS, setRatio(stash, 1)
    'RATIO2_ON': IF (Event.select) THEN drawS, setRatio(stash, 2)

    'LOG_ON': IF (Event.select) THEN $
     drawS, setLog(stash, On)  ELSE  drawS, setLog(stash, Off)

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of monochromatic source function", $
             "", $
             "In 1-D mode the pair of three down arrows in the main panel", $
             "indicates the location of tau = 0.3, 1.0, 3.0. The up arrow", $
             "indicates the location of tau = 1.0 in the background opacity.",$
             "", $
             "In 2-D mode the locations of tau = 1.0 in total and", $
             "background opacity are indicated by the solid and dashed line", $
             "respectively.", $
             "", $
             "Version 1.0, Mar 11, 1999", $
             "Han Uitenbroek (HUitenbroek@cfa.harvard.edu)"])
  ELSE:
  ENDCASE

END
; -------- end ---------------------------- XViewSource_Event.pro ---- ;

; -------- begin -------------------------- SWidgetSetup.pro -- ------ ;

FUNCTION SWidgetSetup, lambda, lambdaDisplay

@geometry.common

  COMMON screen_common, screenSize, scaleFactor

  scaleFactor = [1.0, 1.0]
  screenSize  = fix([550, 400])

  IF (geometryType EQ "TWO_D_PLANE") THEN $
   rayNo = 0 $
  ELSE $
   rayNo = geometry.Nrays - 1

  state = {baseWidget: 0L, drawWidget: 0L, ray: rayNo, rayText: 0L, $
           xmuText: 0L, ymuText: 0L, zmuText: 0L, wmuText: 0L, $           
           lambdaText: 0L, lambdaSlider: 0L, xSlider: 0L, zSlider: 0L, $
           lambdaNo: long(lambdaDisplay), log: 1, displayType: 1, ratio: 0}

  state.baseWidget = widget_base(TITLE='XViewSource', /ROW, MBAR=menuBar, $
                                  RESOURCE_NAME='XViewSource')

  SourceBase = widget_base(state.baseWidget, /COLUMN)

  fileMenu  = widget_button(menuBar, VALUE='File', /MENU)
  openSpectrumButton = widget_button(fileMenu, VALUE='open opacity file', $
                                     UVALUE='NEWOPACITYFILE')
  quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                             RESOURCE_NAME='quitbutton')

  rayMenu = widget_button(menuBar, Value='Rays', $
                          UVALUE='RAYMENU', $
                          /MENU, EVENT_FUNC='rayMenu_Event_Func')
  IF (geometry.Nrays GT 20) THEN Nstep = 5 ELSE Nstep = 1 
  FOR i=0,geometry.Nrays-1, Nstep DO BEGIN
    menuItem = widget_button(rayMenu, UVALUE=i, $
                             VALUE=string(FORMAT='("ray ", I2)', i))
  ENDFOR

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    orientMenu = widget_button(menuBar, VALUE='Orientation', /MENU)
    orientButton = widget_button(orientMenu, VALUE='set orientation', $
                                  UVALUE='ORIENT')
    trackButton = widget_button(orientMenu, VALUE='track surface', $
                                  UVALUE='TRACK_SURFACE')
  ENDIF

  printMenu= widget_button(menuBar, VALUE='Print', /MENU)
  printButton = widget_button(printMenu, VALUE='PostScript', UVALUE='PRINT')
  pngButton   = widget_button(printMenu, VALUE='png', UVALUE='PNG')

  helpMenu   = widget_button(menuBar, VALUE='Help', /MENU, /HELP)
  infoButton = widget_button(helpMenu, VALUE='XViewSource', $
                              UVALUE='INFORMATION')

  ;; --- Draw frame base widget --                      ------------- ;;

  frame     = widget_base(sourceBase, /COLUMN, /FRAME)
  waveFrame = widget_base(frame, /ROW)
  waveLabel = widget_label(waveFrame, VALUE='Wavelength: ')
  state.lambdaText = widget_label(waveFrame, /FRAME, /DYNAMIC_RESIZE)
  waveLabel    = widget_label(waveFrame, VALUE='[nm]')

  rayFrame = widget_base(frame, /ROW)
  raylabel = widget_label(rayFrame, VALUE='Ray:')
  state.rayText = widget_label(rayFrame, /FRAME, $
                                VALUE=string(FORMAT='(I2)', state.ray))
  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    xmulabel      = widget_label(rayFrame, VALUE=' mu_x:')
    state.xmuText = widget_label(rayFrame, /FRAME, /DYNAMIC_RESIZE, $
                                 VALUE=string(FORMAT='(F6.3)', $
                                              geometry.xmu(state.ray)))
    ymulabel      = widget_label(rayFrame, VALUE=' mu_y:')
    state.ymuText = widget_label(rayFrame, /FRAME, /DYNAMIC_RESIZE, $
                                 VALUE=string(FORMAT='(F6.3)', $
                                               geometry.ymu(state.ray)))
    zmu = sqrt(1.0 - $
               (geometry.xmu(state.ray)^2 +  geometry.ymu(state.ray)^2))
    zmulabel      = widget_label(rayFrame, VALUE=' mu_z:')
    state.zmuText = widget_label(rayFrame, /FRAME, $
                               VALUE=string(FORMAT='(F6.3)', zmu))
  ENDIF ELSE BEGIN
    xmulabel      = widget_label(rayFrame, VALUE=' mu:')
    state.xmuText = widget_label(rayFrame, /FRAME, /DYNAMIC_RESIZE, $
                                 VALUE=string(FORMAT='(F6.3)', $
                                              geometry.xmu(state.ray)))
  ENDELSE
  wmulabel      = widget_label(rayFrame, VALUE=' weight:')
  state.wmuText = widget_label(rayFrame, /FRAME, /DYNAMIC_RESIZE, $
                                VALUE=string(FORMAT='(E10.4)', $
                                             geometry.wmu(state.ray)))

  drawFrame = widget_base(SourceBase, /FRAME)
  state.drawWidget = widget_draw(drawFrame, $
                                 XSIZE=screenSize[0], YSIZE=screenSize[1])

  makeupFrame  = widget_base(SourceBase, /ROW, /FRAME)
  linlogFrame  = widget_base(makeupFrame, /COLUMN, /EXCLUSIVE)
  linearButton = widget_button(linlogFrame, VALUE='Linear', $
                                UVALUE='LOG_OFF')
  logButton    = widget_button(linlogFrame, VALUE='Logarithmic', $
                                UVALUE='LOG_ON')
  widget_control, logButton, /SET_BUTTON
  
  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    shadeFrame   = widget_base(makeupFrame, /COLUMN, /EXCLUSIVE)
    panelButton  = widget_button(shadeFrame, VALUE='Panel', $
                                 UVALUE='PANEL')
    wireButton   = widget_button(shadeFrame, VALUE='Wire Grid', $
                                  UVALUE='WIREGRID')
    shadeButton  = widget_button(shadeFrame, VALUE='Shade Surface', $
                                  UVALUE='SHADESURF')
    widget_control, wireButton, /SET_BUTTON

    ratioFrame = widget_base(makeupFrame, /COLUMN, /EXCLUSIVE)
    SButton  = widget_button(ratioFrame, VALUE='S',   UVALUE='RATIO_OFF')
    JSButton = widget_button(ratioFrame, VALUE='S/J', UVALUE='RATIO1_ON')
    BSButton = widget_button(ratioFrame, VALUE='S/B', UVALUE='RATIO2_ON')
    widget_control, SButton, /SET_BUTTON
  ENDIF

  waveFrame = widget_base(makeupFrame)
  state.lambdaSlider = widget_slider(waveFrame, TITLE='Wavelength No:', $
                                     UVALUE='SLIDER', XSIZE=255, $
                                     MAXIMUM=n_elements(lambda) - 1, $
                                     VALUE=lambdaDisplay, SCROLL=10)

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- SWidgetSetup.pro --------- ;

; -------- begin -------------------------- oplotIntensity.pro ------- ;

PRO oplotIntensity, state

@spectrum.common
@geometry.common

  IF (!D.NAME EQ 'X') THEN BEGIN
    widget_control, state.drawWidget, GET_VALUE=WindowNo
    wset, WindowNo
  ENDIF

  oldxcrange = !X.CRANGE
  oldregion = !P.REGION
  !P.REGION = [0.4, 0.6, 0.8, 0.95]  &  !P.NOERASE=1

  IF (geometryType EQ "TWO_D_PLANE" AND state.displayType EQ 0) THEN $
   color = 255B $
  ELSE $
   color = !P.COLOR

  IF (spectrum.lambda[lambdaDisplay] GT 1.0E+03) THEN BEGIN
    xtitle = '!x[micron]'
    lambda = spectrum.lambda / 1.0E3
    range  = [-0.002, 0.002] / 3.0
  ENDIF ELSE BEGIN
    xtitle = '!x[nm]'
    lambda = spectrum.lambda
    range  = [-0.2, 0.2]
  ENDELSE

  plot, lambda, spectrum.I[*, state.ray], XSTYLE=8, $
   YSTYLE=4, XRANGE=lambda[lambdaDisplay] + range, $
   XTITLE='!7l' + xtitle, CHARSIZE=0.7, /YLOG, COLOR=color

  x = lambda[lambdaDisplay]
  y = spectrum.I[lambdaDisplay, state.ray]
  arrow_bottom = y / 10^(0.2 * (!Y.CRANGE[1] - !Y.CRANGE[0])) 
  arrow, /DATA, /SOLID, x, arrow_bottom, x, y, COLOR=200B

  !P.REGION = oldregion  &  !P.NOERASE=0
  !X.CRANGE = oldxcrange
END
; -------- end ---------------------------- oplotIntensity.pro ------- ;

; -------- begin -------------------------- drawS.pro ---------------- ;

PRO drawS, state, TRACK=track

@geometry.common
@atmos.common
@spectrum.common
@opacity.common

  COMMON screen_common, screenSize, scaleFactor

  KM_TO_M = 1.0E+03

  Jcolor = 200B  &  cColor = 80B  &  PlanckColor = 16B
  asColor = 160B

  IF (!D.NAME EQ 'X') THEN BEGIN
    widget_control, state.drawWidget, GET_VALUE=WindowNo
    wset, WindowNo
    thick = 2
  ENDIF ELSE $
   thick = 6

  readJ, state.lambdaNo
  readOpacity, state.lambdaNo, state.ray

  Bp = Planck(atmos.T, spectrum.lambda[state.lambdaNo], /HZ)
  S  = (eta_as + eta_c + J*scatt) / (chi_as + chi_c)

  CASE (state.ratio) OF 
    0: BEGIN
      sv = S
      IF (geometryType EQ "TWO_D_PLANE") THEN $
       title = 'Source Function [J m!U-2!N s!U-1!N Hz!U-1!N sr!U-1!N]' $
      ELSE $
       title = 'S, J [J m!U-2!N s!U-1!N Hz!U-1!N sr!U-1!N]'
    END
    1: BEGIN
      sv = S / J
      title = 'S/J'
    END
    2: BEGIN
      sv = S / Bp
      title = 'S/B'
    END
  ENDCASE

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    xkm = geometry.x / KM_TO_M
    zkm = geometry.z / KM_TO_M

    IF (keyword_set(TRACK)) THEN BEGIN
      wavelength = strtrim(string(spectrum.lambda[state.lambdaNo], $
                          FORMAT='(F9.3, " [nm]")'), 2)
      CASE (state.ratio) OF
        0: tracktitle = 'source function: S  at lambda = '   + wavelength
        1: tracktitle = 'source function: S/J  at lambda = ' + wavelength
        2: tracktitle = 'source function: S/B  at lambda = ' + wavelength
      END
      IF (state.log) THEN tracktitle = 'log of ' + tracktitle
      surf_track, sv, xkm, zkm, TITLE=tracktitle, ZLOG=state.log
    ENDIF ELSE BEGIN

      zmu = sqrt(1.0 - geometry.xmu[state.ray]^2 + geometry.ymu[state.ray]^2)

      xtauone_c = fltarr(geometry.Nx)   &  xtauone_tot = xtauone_c
      ztauone_c = fltarr(geometry.Nx)   &  ztauone_tot = ztauone_c

      s1_c = fltarr(geometry.Nx)        &  s1_tot = s1_c

      FOR l=1, geometry.Nx-1 DO BEGIN
        vis  = raytrace(geometry, state.ray, l, XRAY=xray, ZRAY=zray)
        path = (zray - zray[0]) / zmu

        peff_c   = tauone(path, rayinterpolate(chi_c, vis))
        peff_tot = tauone(path, rayinterpolate(chi_c + chi_as, vis))

        xtauone_c[l]   = interpolate(xray, peff_c)
        xtauone_tot[l] = interpolate(xray, peff_tot)
        ztauone_c[l]   = interpolate(zray, peff_c)
        ztauone_tot[l] = interpolate(zray, peff_tot)

        sray   = rayinterpolate(sv, vis)
        s1_c[l]   = interpolate(sray, peff_c)
        s1_tot[l] = interpolate(sray, peff_tot)

      ENDFOR

      CASE (state.displayType) OF 
      0: BEGIN
        erase
        panel, scaleimg_idl(sv, screenSize[0]-155, screenSize[1]-100),$
         XPOS=60, YPOS=50, xkm, zkm, /ORTHOSCOPIC, /ORDER, $
         XTITLE='x [km]', YTITLE='z [km]', ZLOG=state.log, SCALETEXT=title
        oplot, xtauone_c/KM_TO_M, ztauone_c/KM_TO_M, THICK=thick, COLOR=255B, PSYM=2
        as_non_zero = where(chi_as GT 0.0, count)
        IF (count GT 0) THEN $
         oplot, xtauone_tot/KM_TO_M, ztauone_tot/KM_TO_M, THICK=thick, COLOR=220B, PSYM=4

        IF (!D.NAME EQ 'PS') THEN rhannotate, xann(0.05), yann(0.925), $
         TEXT=string(FORMAT='(F9.3, " [nm]")', $
                     spectrum.lambda[state.lambdaNo]), CHARCOLOR=255B
      END
      1: surface, sv, xkm, zkm, /T3D, COLOR=!P.COLOR, $
         CHARSIZE=1.4, XTITLE='x', YTITLE='z', ZTITLE=title, $
         BOTTOM=200B, ZLOG=state.log
      2: shade_surf, sv, xkm, zkm, /T3D, $
         CHARSIZE=1.4, XTITLE='x', YTITLE='z', ZTITLE=title, ZLOG=state.log
      ENDCASE

      IF (state.displayType GT 0) THEN BEGIN
        plots, xtauone_c/KM_TO_M, ztauone_c/KM_TO_M, s1_c, /T3D, THICK=thick, $
         LINE=2, COLOR=PlanckColor
        as_non_zero = where(chi_as GT 0.0, count)
        IF (count GT 0) THEN $
         plots, xtauone_c/KM_TO_M, ztauone_c/KM_TO_M, s1_tot, $
         /T3D, THICK=thick, COLOR=PlanckColor
      ENDIF
    ENDELSE
  ENDIF ELSE BEGIN
    xmu = geometry.xmu[state.ray]

    IF (geometryType EQ "ONE_D_PLANE") THEN BEGIN
      tau   = getTau(geometry.height, (chi_as + chi_c) / xmu)
      tau_c = getTau(geometry.height, chi_c/xmu)
    ENDIF ELSE BEGIN
      tau   = getTau(geometry.r, (chi_as + chi_c) / xmu)
      tau_c = getTau(geometry.r, chi_c/xmu)
    ENDELSE

    as_non_zero = where(chi_as GT 0.0, count)
    IF (count GT 0) THEN $
     S_as = eta_as[as_non_zero] / chi_as[as_non_zero]
    S_c = (eta_c + scatt*J) / chi_c

    plot, geometry.cmass, S, THICK=thick, $
     XTITLE='Column Mass [kg m!U-2!N]', YTITLE=title, /XLOG, YLOG=state.log
    xann05 = xann(0.05)
    rhannotate, xann05, yann(0.95), THICK=thick, $
     TEXT='S!Dtotal!N', MARKTYPE=0

    rhannotate, xann05, yann(0.75), TEXT='J', MARKTYPE=0, $
     MARKCOLOR=JColor, THICK=thick
    oplot, geometry.cmass, J, THICK=thick, COLOR=JColor

    ;; --- Mark optical depth unity in total opacity -- ;

    tabinv, tau, [0.3, 1.0, 3.0], tau_unity
    x = linear(geometry.cmass, tau_unity)
    y = linear(S, tau_unity)
    IF (state.log) THEN $
     arrow_top = y * 10^(0.06 * (!Y.CRANGE[1] - !Y.CRANGE[0])) $
    ELSE $
     arrow_top = y + 0.06 * (!Y.CRANGE[1] - !Y.CRANGE[0])
    arrow, /DATA, /SOLID, THICK=2, x, arrow_top, x, y, COLOR=asColor

    oplot, geometry.cmass, Bp, THICK=thick, LINESTYLE=2, COLOR=PlanckColor
    rhannotate, xann05, yann(0.8), TEXT='B!DPlanck!N', $
     MARKTYPE=2, THICK=thick, MARKCOLOR=PlanckColor

    IF (count GT 0) THEN BEGIN 
      oplot, geometry.cmass[as_non_zero], S_as, COLOR=asColor, THICK=thick
      rhannotate, xann05, yann(0.9), TEXT='S!Dactive!N', $
       MARKTYPE=0, MARKCOLOR=asColor, THICK=thick
    ENDIF

    oplot, geometry.cmass, S_c, COLOR=cColor, THICK=thick
    rhannotate, xann05, yann(0.85), TEXT='S!Dbackgr!N', $
     MARKTYPE=0, MARKCOLOR=cColor, THICK=thick

    ;; --- Mark optical depth unity in background opacity        -- ;

    tabinv, tau_c, 1.0, tau_unity
    x = linear(geometry.cmass, tau_unity)
    y = linear(S_c, tau_unity)
    IF (state.log) THEN $
     arrow_bottom = y / 10^(0.06 * (!Y.CRANGE[1] - !Y.CRANGE[0])) $
    ELSE $
     arrow_bottom = y - 0.06 * (!Y.CRANGE[1] - !Y.CRANGE[0])
    arrow, /DATA, /SOLID, THICK=2, x, arrow_bottom, x, y, COLOR=cColor
    ENDELSE

   oplotIntensity, state
END
; -------- end ---------------------------- drawS.pro ---------------- ;

; -------- begin -------------------------- XViewSource.pro ---------- ;

PRO XViewSource, GROUP_LEADER=group_leader

;+
; NAME:
;	XVIEWSOURCE
;
; PURPOSE:
;
; CATEGORY:
;	Data reduction
;
; CALLING SEQUENCE:
;
; INPUTS:
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
; OPTIONAL OUTPUTS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; EXAMPLE:
; MODIFICATION HISTORY:
;
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Wed Jun 30 08:37:54 2010 --
;-

@atmos.common
@spectrum.common
@files.common

  IF (n_elements(atmos) EQ 0) THEN $
   IF (NOT readAtmos(atmosFile)) THEN return
  IF (n_elements(spectrum) EQ 0) THEN $
   IF (NOT readSpectrum(spectrumFile)) THEN return

  IF (NOT keyword_set(GROUP_LEADER)) THEN BEGIN
    group_leader  = 0
    surfr, AX=30, AZ=30
    lambdaDisplay = 0L
    JFile    = 'J.dat'
    opacFile = 'opacity.out'
    backgroundFile = 'background.dat'
    IF (NOT openJ(JFile)) THEN return
  ENDIF
  IF (NOT openOpacity(opacFile)) THEN return

  state = SWidgetSetup(spectrum.lambda, lambdaDisplay)
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader
  surfr, AX=30, AZ=30
  drawS, setLambdaNo(widget_info(state.baseWidget, /CHILD), lambdaDisplay)

  ;; --- Register with the XManager --               -------------- ;;

  xmanager, 'XViewSource', state.baseWidget, $
   EVENT_HANDLER='XViewSource_Event', GROUP_LEADER=group_leader

  IF (NOT keyword_set(GROUP_LEADER)) THEN BEGIN
    free_lun, Junit
    free_lun, opacunit
  ENDIF
END
; -------- end ---------------------------- XViewS.pro ------------- ;
