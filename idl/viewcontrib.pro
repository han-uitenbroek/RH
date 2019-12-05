; ----------------------------------------- viewcontrib.pro ---------- ;

; -------- begin -------------------------- XViewContrib_Event.pro --- ;

PRO XViewContrib_Event, Event

@files.common
@geometry.common

  COMMON contrib_conserve, contrib, Nlamb, x
  COMMON screen_common, screenSize, scaleFactor

  ;; --- Main event handler --                          -------------- ;

  widget_control, Event.id, GET_UVALUE=Action
  stash = widget_info(Event.top, /CHILD)
  widget_control, stash, GET_UVALUE=state

  CASE Action OF
    'QUIT': BEGIN

      ; --- Free memory and quit --                     -------------- ;

      contrib = 0.0  &  x = 0.0
      widget_control, Event.top, /DESTROY
    END

    'PRINT': BEGIN
      filename = '/tmp/viewContrib-' + timeStamp() + '.ps'
      IF (state.panel) THEN font = 0 ELSE font = -1
      psopen, FILENAME=filename, /COLOR, FONT=font
      pswindow, XSIZE=screenSize[0], YSIZE=screenSize[1]
      drawContrib, state
      psclose
      scaleFactor = [1.0, 1.0]
      r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
    END
    'PRINT_BW': BEGIN
      filename = '/tmp/viewContrib-' + timeStamp() + '.ps'
      IF (state.panel) THEN font = 0 ELSE font = -1
      psopen, FILENAME=filename, /COLOR, FONT=font
      pswindow, XSIZE=screenSize[0], YSIZE=screenSize[1]

      ;; Print B&W in reversed colors (for print publications).

      tvlct, reverse(indgen(256)), reverse(indgen(256)), reverse(indgen(256))
      !P.COLOR = 255B

      drawContrib, state
      psclose
      scaleFactor = [1.0, 1.0]
      r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
    END

    'PNG': BEGIN
      widget_control, state.drawWidget, GET_VALUE=WindowNo
      rhwritepng, WindowNo, '/tmp/viewContrib-'
    END

    'ORIENT': orient, state, EVENT_HANDLER='XViewContrib_Event', $
     GROUP_LEADER=Event.top, TITLE='Orientation Contribution Function'
  
    'ORIENT_UPDATE': BEGIN
      OrientationUpdate, state
      drawContrib, state
    END

    'ORIENT_RESET': BEGIN
      OrientationUpdate, state, /SET
      drawContrib, state
    END

    'SETREDBLUE': BEGIN
      widget_control, state.blueField, GET_VALUE=lambdablue
      widget_control, state.redField,  GET_VALUE=lambdared
      drawContrib, setBlue_n_Red(stash, lambdablue, lambdared)
    END

    'RAYMENU': drawContrib, setRayNo(stash, Event.value)

    'LOG_ON': IF (Event.select) THEN $
     drawContrib, setLog(stash, 1)  ELSE  drawContrib, setLog(stash, 0)

    'PANEL_ON': BEGIN
      state.panel = (Event.select) ? 1 : 0
      widget_control, state.orientButton, SENSITIVE=(Event.select) ? 0 : 1
      widget_control, stash, SET_UVALUE=state
      drawContrib, state
    END
    
    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of contribution function", $
             "", $
             "Version 1.0, Mar 29, 1996", $
             "Han Uitenbroek (HUitenbroek@cfa.harvard.edu)"])
  ELSE:
  ENDCASE

END
; -------- end ---------------------------- XViewContrib_Event.pro ---- ;

; -------- begin -------------------------- ContribWidgetSetup.pro ---- ;

FUNCTION ContribWidgetSetup, lambda, blue, red, type

@geometry.common

  COMMON screen_common, screenSize, scaleFactor

  scaleFactor = [1.0, 1.0]
  screenSize  = fix([600, 450])

  Nspect = n_elements(lambda)

  IF (geometryType EQ "TWO_D_PLANE") THEN $
   rayNo = 0 $
  ELSE $
   rayNo = geometry.Nrays - 1

  state = {baseWidget: 0L, drawWidget: 0L, $
           blueText: 0L, redText: 0L, rayText: 0L, $
           xmuText: 0L, ymuText: 0L, wmuText: 0L, $
           blueField: 0L, redField: 0L, xSlider: 0L, zSlider: 0L, $
           ray: rayNo, log: 0, panel: 0,  blue: blue,  red: red, $
           lambdablue: lambda[blue], lambdared: lambda[red], $
           orientButton: 0L, type: type}

  state.baseWidget = widget_base(TITLE='XViewContrib', /ROW, MBAR=menuBar, $
                                  RESOURCE_NAME='XViewContrib')

  spectrumBase = widget_base(state.baseWidget, /COLUMN)

  fileMenu  = widget_button(menuBar, VALUE='File', /MENU)
  quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                              RESOURCE_NAME='quitbutton')

  rayMenu = widget_button(menuBar, Value='Rays', $
                           UVALUE='RAYMENU', $
                           /MENU, EVENT_FUNC='rayMenu_Event_Func')
  IF (geometry.Nrays GT 20) THEN Nstep = 5 ELSE Nstep = 1 
  FOR i=0,geometry.Nrays-1, Nstep DO BEGIN
    menuItem = widget_button(rayMenu, UVALUE=i, $
                              VALUE=string(FORMAT='("ray ", I2)', i))
  endfor

  orientMenu = widget_button(menuBar, VALUE='Orientation', /MENU)
  state.orientButton = widget_button(orientMenu, VALUE='set orientation', $
                                     UVALUE='ORIENT')

  printMenu= widget_button(menuBar, VALUE='Print', /MENU)
  printButton = widget_button(printMenu, VALUE='PostScript', UVALUE='PRINT')
  printButton = widget_button(printMenu, VALUE='PostScript (B&W)', $
                              UVALUE='PRINT_BW')
  pngButton   = widget_button(printMenu, VALUE='png', UVALUE='PNG')


  helpMenu   = widget_button(menuBar, VALUE='Help', /MENU, /HELP)
  infoButton = widget_button(helpMenu, VALUE='XViewContrib', $
                              UVALUE='INFORMATION')

  drawFrame    = widget_base(spectrumBase, /FRAME, /COLUMN)

  waveFrame     = widget_base(drawFrame, /ROW)
  raylabel      = widget_label(waveFrame, VALUE='Ray:')
  state.rayText = widget_label(waveFrame, /FRAME, $
                                VALUE=string(FORMAT='(I2)', state.ray))

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    xmulabel      = widget_label(waveFrame, VALUE=' mu_x:')
    state.xmuText = widget_label(waveFrame, /FRAME, /DYNAMIC_RESIZE, $
                                  VALUE=string(FORMAT='(F6.3)', $
                                               geometry.xmu(state.ray)))
    ymulabel      = widget_label(waveFrame, VALUE=' mu_y:')
    state.ymuText = widget_label(waveFrame, /FRAME, /DYNAMIC_RESIZE, $
                                  VALUE=string(FORMAT='(F6.3)', $
                                               geometry.ymu(state.ray)))
  ENDIF ELSE BEGIN
    xmulabel      = widget_label(waveFrame, VALUE=' mu:')
    state.xmuText = widget_label(waveFrame, /FRAME, /DYNAMIC_RESIZE, $
                                  VALUE=string(FORMAT='(F6.3)', $
                                               geometry.xmu(state.ray)))
  ENDELSE
  wmulabel      = widget_label(waveFrame, VALUE=' weight:')
  state.wmuText = widget_label(waveFrame, /FRAME, /DYNAMIC_RESIZE, $
                                VALUE=string(FORMAT='(E10.4)', $
                                             geometry.wmu(state.ray)))

  waveFrame    = widget_base(drawFrame, /ROW)
  waveLabel    = widget_label(waveFrame, VALUE='Wavelength #: ')
  state.blueText = widget_label(waveFrame, /FRAME, /DYNAMIC_RESIZE, $
                VALUE=string(FORMAT='(I6)', blue))
  waveLabel    = widget_label(waveFrame, VALUE=' to ')
  state.redText = widget_label(waveFrame, /FRAME, /DYNAMIC_RESIZE, $
                 VALUE=string(FORMAT='(I6)', red))

  ;; --- Draw frame base widget --                      --------------- ;;

  state.drawWidget = widget_draw(drawFrame, $
                                 XSIZE=screenSize[0], YSIZE=screenSize[1])

  makeupFrame  = widget_base(spectrumBase, /ROW)
  linlogFrame  = widget_base(makeupFrame, /COLUMN, /EXCLUSIVE, /FRAME)
  linearButton = widget_button(linlogFrame, VALUE='Linear', $
                                UVALUE='LOG_OFF')
  logButton    = widget_button(linlogFrame, VALUE='Logarithmic', $
                                UVALUE='LOG_ON')
  widget_control, linearButton, /SET_BUTTON
  
  panelFrame   = widget_base(makeupFrame, /COLUMN, /EXCLUSIVE, /FRAME)
  wireButton   = widget_button(panelFrame, VALUE='Wire Grid', $
                                UVALUE='SHADE_OFF')
  shadeButton  = widget_button(panelFrame, VALUE='Panel', $
                                UVALUE='PANEL_ON')
  widget_control, wireButton, /SET_BUTTON

  waveFrame = widget_base(makeupFrame, /FRAME, /ROW)
  state.blueField = cw_field(waveFrame, TITLE='Wavelength', /FLOATING, $
                             UVALUE='SETREDBLUE', VALUE=state.lambdablue, $
                             XSIZE=12, /RETURN_EVENT)
  state.redField  = cw_field(waveFrame, TITLE='  to ', /FLOATING, $
                             UVALUE='SETREDBLUE', VALUE=state.lambdared, $
                             XSIZE=12, /RETURN_EVENT)
  waveLabel = widget_label(waveFrame, VALUE='[nm]')

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- ContribWidgetSetup.pro ---- ;

; -------- begin -------------------------- drawContrib.pro ----------- ;

PRO drawcontrib, state

@geometry.common
@spectrum.common
@opacity.common


  KM_TO_M = 1.0E3
  MIN_FRACTION = 1.0E-6

  IF (!D.NAME EQ 'X') THEN BEGIN
    widget_control, state.drawWidget, GET_VALUE=WindowNo
    wset, WindowNo
  ENDIF

  blue = state.blue
  red  = state.red

  widget_control, state.blueText, SET_VALUE=string(FORMAT='(I6)', blue)
  widget_control, state.redText, SET_VALUE=string(FORMAT='(I6)', red)

  x = geometry.height / geometry.xmu[state.ray]

  Nlamb = red - blue + 1
  S = fltarr(geometry.Ndep, Nlamb)
  tau = S

  dx = x - shift(x, 1)
  dx[0] = dx[1]

  FOR la=0, Nlamb-1 DO BEGIN
    readJ, blue + la
    readOpacity, blue + la, state.ray

    tau[*, la]  = getTau(x, chi_as + chi_c)
    dtau    = tau[*, la] - shift(tau[*, la], 1)
    dtau[0] = tau[0, la]

    CASE state.type OF
      0: S[*, la] = (eta_as + eta_c + J*scatt) / (chi_c + chi_as)
      1: S[*, la] = eta_as / (chi_c + chi_as)
      2: S[*, la] = (eta_c + J*scatt) / (chi_c + chi_as)
    ENDCASE

    S[*, la] = S[*, la] * (-dtau / dx)  
  ENDFOR
  contrib = transpose(S * exp(-(tau < 50.0))) * $
   KM_TO_M / geometry.xmu[state.ray]

  z = geometry.height / KM_TO_M

  ytitle = 'Height [km]'
  scaleText = $
   'Contribution function [J m!U-2!N s!U-1!N Hz!U-1!N sr!U-1!N km!U-1!N]'
  IF (state.panel) THEN BEGIN
    erase
    panel, scaleimg_idl(contrib, 400, 350) > (MIN_FRACTION * max(contrib)), $
     spectrum.lambda[blue:red], z, $
     XTITLE='Wavelength [nm]', YTITLE=ytitle, XPOS=75, YPOS=50, $
     SCALETEXT=scaleText, /ORTHOSCOPIC, /ORDER, ZLOG=state.log

    oplot, spectrum.lambda[state.blue:state.red], $
     spectrum.I[state.blue:state.red, state.ray] / $
     max(spectrum.I[state.blue:state.red, state.ray]) * $
     0.75 * (!Y.CRANGE[1] - !Y.CRANGE[0]) + $
     0.10 * (!Y.CRANGE[1] - !Y.CRANGE[0]), COLOR=255B
  ENDIF ELSE BEGIN
    surface, contrib > (MIN_FRACTION * max(contrib)), $
     spectrum.lambda[blue:red], z, ZLOG=state.log, $
     XTITLE='Wavelength [nm]', YTITLE=ytitle, $
     ZTITLE=scaleText, /T3D, CHARSIZE=1.4
  ENDELSE

END
; -------- end ---------------------------- drawContrib.pro ----------- ;

; -------- begin -------------------------- XViewContrib.pro ---------- ;

PRO XViewContrib, BLUE=blue, RED=red, GROUP_LEADER=group_leader, $
                  CONTINUUM=continuum, ACTIVE_SET=active_set

@geometry.common
@spectrum.common
@files.common

  IF (NOT keyword_set(GROUP_LEADER)) THEN group_leader=0

  IF (n_elements(atmos) EQ 0) THEN $
   IF (NOT readAtmos(atmosFile)) THEN return
  IF (n_elements(spectrum) EQ 0) THEN $
   IF (NOT readSpectrum(spectrumFile)) THEN return

  IF (NOT openOpacity(opacFile)) THEN return

  IF (NOT keyword_set(BLUE))  THEN blue = 0
  IF (NOT keyword_set(RED))   THEN red  = spectrum.Nspect - 1
  IF (keyword_set(ACTIVE_SET)) THEN $
   type = 1 $
  ELSE IF (keyword_set(CONTINUUM)) THEN $
   type = 2 $
  ELSE $
   type = 0

  state = ContribWidgetSetup(spectrum.lambda, blue, red, type)
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader

  drawContrib, state

  ;; --- Register with the XManager --                   ------------ ;;

  xmanager, 'XViewContrib', state.baseWidget, $
   EVENT_HANDLER='XViewContrib_Event', GROUP_LEADER=group_leader
END
; -------- end ---------------------------- XViewContrib.pro --------- ;
