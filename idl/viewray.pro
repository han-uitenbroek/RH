FUNCTION setQuantity, stash, qstring

@geometry.common
@atmos.common
@spectrum.common
@opacity.common

  widget_control, stash, GET_UVALUE=state
  IF (qstring EQ 'Source function') THEN $
   widget_control, state.waveslider, SENSITIVE=1 $
  ELSE $
   widget_control, state.waveslider, SENSITIVE=0

  state.qstring = qstring

  CASE (state.qstring) OF
    'Temperature': BEGIN
      state.f = atmos.T
      state.ztitle = 'Temperature [K]'
    END
    
    'Vlos': BEGIN
      mux = geometry.xmu[state.ray]
      muz = sqrt(1.0 - (mux^2 + geometry.ymu[state.ray]^2))
      KM_TO_M = 1.0E3
      state.f = (mux*geometry.vx + muz*geometry.vz) / KM_TO_M
      state.ztitle = 'LOS velocity [km/s]'
    END

    'Blos': BEGIN
      IF (tag_present(atmos, 'STOKES')) THEN BEGIN
        mux = geometry.xmu[state.ray]
        muy = geometry.ymu[state.ray]
        muz = sqrt(1.0 - (mux^2 +muy ^2))
        Bxe = cos(atmos.chi_B) * sin(atmos.gamma_B)
        Bye = sin(atmos.chi_B) * sin(atmos.gamma_B)
        Bze = cos(atmos.gamma_B)
        state.f = atmos.B * (mux*Bxe + muy*Bye + muz*Bze)
        state.ztitle = 'LOS magn. field [T]'
      ENDIF
    END

    'Field Angle': BEGIN
      IF (tag_present(atmos, 'STOKES')) THEN BEGIN
        mux = geometry.xmu[state.ray]
        muy = geometry.ymu[state.ray]
        muz = sqrt(1.0 - (mux^2 +muy ^2))
        Bxe = cos(atmos.chi_B) * sin(atmos.gamma_B)
        Bye = sin(atmos.chi_B) * sin(atmos.gamma_B)
        Bze = cos(atmos.gamma_B)
        state.f = mux*Bxe + muy*Bye + muz*Bze
        state.ztitle = 'cosine of field agngle'
      ENDIF
    END

    'Source function': BEGIN
      widget_control, state.waveslider, GET_VALUE=nspect
      readJ, nspect
      readOpacity, nspect, state.ray
      state.f = (eta_as + eta_c + J*scatt) / (chi_c + chi_as)
      state.ztitle = 'Source function [J m!U-2!N s!U-1!N Hz!U-1!N sr!U-1!N]'
      state.lambda = spectrum.lambda[nspect]
    END
  ENDCASE

  widget_control, stash, SET_UVALUE=state
  return, state
END


; -------- begin -------------------------- XViewRay_Event.pro ------- ;

PRO XViewRay_Event, Event

  ;; --- Main event handler --                          -------------- ;

  COMMON screen_common, screenSize, scaleFactor

  stash = widget_info(Event.handler, /CHILD)
  widget_control, stash, GET_UVALUE=state

  widget_control, Event.id, GET_UVALUE=Action
  CASE Action OF
    'QUIT': widget_control, Event.top, /DESTROY

    'PRINT': BEGIN
      filename = 'viewRay-' + timeStamp() + '.ps'
      PSopen, FILENAME=filename, /COLOR, XSIZE=7.0, YSIZE=6.5
      ScaleFactor = [!D.X_SIZE, !D.Y_SIZE] / float(ScreenSize)
      displayRay, state
      PSclose
      ScaleFactor = [1.0, 1.0]
      r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
    END

    'PNG': BEGIN
      widget_control, state.drawWidget, GET_VALUE=WindowNo
      rhwritepng, WindowNo, 'viewRay-'
    END

    'RAYMENU': BEGIN
      state = setRayNo(stash, Event.value)
      displayRay, setQuantity(stash, state.qstring)
    END

    'SLIDER': BEGIN
      displayRay, state
    END

    'WAVESLIDER': BEGIN
      displayRay, setQuantity(stash, state.qstring)
    END

    'VR_TEMPERATURE': displayRay, setQuantity(stash, "Temperature")

    'VR_VLOS': displayRay, setQuantity(stash, "Vlos")

    'VR_BLOS': displayRay, setQuantity(stash, "Blos")

    'VR_ANGLE': displayRay, setQuantity(stash, "Field Angle")

    'VR_SOURCE': displayRay, setQuantity(stash, "Source function")

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of quantities along ray path", $
             "", $
             "Han Uitenbroek (HUitenbroek@nso.edu)"])
  ELSE:
  ENDCASE
END
; -------- end ---------------------------- XViewRay_Event.pro ------- ;

PRO displayRay, state

  MM_TO_M = 1.0E6

  COMMON screen_common, screenSize, scaleFactor

@geometry.common
@spectrum.common
@opacity.common

  IF (!D.NAME EQ 'X') THEN BEGIN
    widget_control, state.drawWidget, GET_VALUE=WindowNo
    wset, WindowNo
    erase
  ENDIF
  widget_control, state.slider, GET_VALUE=x_index
  widget_control, state.xlabel, $
   SET_VALUE=string(geometry.x[x_index]/MM_TO_M, FORMAT='(F7.3)')
  widget_control, state.waveslider, GET_VALUE=w_index
  widget_control, state.wlabel, $
   SET_VALUE=string(spectrum.lambda[w_index], FORMAT='(F9.3)')

  vis = raytrace(geometry, state.ray, x_index, XRAY=xray, ZRAY=zray)

  IF (!D.NAME EQ 'PS') THEN CHARSIZE=0.6 ELSE CHARSIZE=1.0

  panel, scaleimg_idl(state.f, screenSize[0]-250, 300), /ORTHOSCOPIC, $
   geometry.x/MM_TO_M, geometry.z/MM_TO_M, SCALETEXT=state.ztitle, $
   XTITLE='x [Mm]', YTITLE='z [Mm]', XPOS=125, YPOS=300, /ORDER, $
   TITLE=string(FORMAT='("x = ", I4, ", mu = ", I2)', x_index, state.ray), $
         CHARSIZE=CHARSIZE
  oplot, xray/MM_TO_M, zray/MM_TO_M, PSYM=4, SYMSIZE=0.5, COLOR=255B

  fray = rayinterpolate(state.f, vis)
  plot, zray/MM_TO_M, fray, XTITLE='z [Mm]', YTITLE=state.ztitle, $
   /DEVICE, POSITION=PSposition([125, 50, screenSize[0]-125, 250]), /NOERASE, $
        CHARSIZE=CHARSIZE
  oplot, zray/MM_TO_M, fray, /PSYM, SYMSIZE=0.5, COLOR=200B

  IF (state.qstring EQ "Source function") THEN BEGIN
    mux = geometry.xmu[state.ray]
    muz = sqrt(1.0 - (mux^2 + geometry.ymu[state.ray]^2))
    chiray = rayinterpolate(chi_as + chi_c, vis)
    tau = gettau(zray/muz, chiray)

    dz = zray - shift(zray, 1)
    dz[0] = dz[1]
    dtau    = tau - shift(tau, 1)
    dtau[0] = tau[0]

    contrib = MM_TO_M * state.f*exp(-(tau < 50.0)) * (-dtau / dz)

    oplot, zray/MM_TO_M, 0.9*max(fray)/max(contrib) * contrib, $
     COLOR=100B, PSYM=4, SYMSIZE=0.5
    rhannotate, TEXT='Contribution function', xann(0.6), yann(0.9), $
     MARKTYPE=-4, MARKCOLOR=100B, SYMSIZE=0.5
    
    rhannotate, TEXT='lambda = ' + string(state.lambda, FORM='(F10.4)'), $
                xann(0.25), yann(0.9)

  ENDIF
END

; -------- begin -------------------------- rayWidgetSetup.pro ------- ;

FUNCTION rayWidgetSetup, f, rayNo, x_index, ztitle

@geometry.common
@atmos.common
@spectrum.common

  COMMON screen_common, screenSize, scaleFactor

  scaleFactor = [1.0, 1.0]
  screenSize  = fix([1300, 650])

  state = {baseWidget: 0L, drawWidget: 0L, log: 0, $
           slider: 0L, xlabel: 0L, f: f, ray: rayNo, rayText: 0L, $
           xmuText: 0L, ymuText: 0L, zmuText: 0L, wmuText: 0L, $
           ztitle:  ztitle, waveslider: 0L, wlabel: 0L, $
           qstring: "Temperature", lambda: 0.0}

  state.baseWidget = widget_base(TITLE='XViewRay', /COLUMN, $
				 RESOURCE_NAME='XViewRay', $
				 MBAR=menuBar)

  base = widget_base(state.baseWidget, /COLUMN)

  fileMenu = widget_button(menuBar, VALUE='File', /MENU)
  button = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                         RESOURCE_NAME='quitbutton')

  printMenu= widget_button(menuBar, VALUE='Print', /MENU)
  printButton = widget_button(printMenu, VALUE='PostScript', UVALUE='PRINT')
  pngButton   = widget_button(printMenu, VALUE='png', UVALUE='PNG')

  helpMenu = widget_button(menuBar, VALUE='Help', /HELP)
  button = widget_button(helpMenu, UVALUE='INFORMATION', $
                         VALUE='about XViewRay')

  rayMenu = widget_button(menuBar, Value='Rays', $
                           UVALUE='RAYMENU', $
                           /MENU, EVENT_FUNC='rayMenu_Event_Func')

  FOR i=0,geometry.Nrays-1 DO BEGIN
    menuItem = widget_button(rayMenu, UVALUE=i, $
                              VALUE=string(FORMAT='("ray ", I2)', i))
  ENDFOR

  quantityMenu = widget_button(menuBar, VALUE='Quantity', UVALUE='QUANTITY')
  vr_tempButton = widget_button(quantityMenu, VALUE='Temperature', $
                                UVALUE='VR_TEMPERATURE')
  vr_vlosButton = widget_button(quantityMenu, VALUE='LOS velocity', $
                                UVALUE='VR_VLOS')
  vr_blosButton = widget_button(quantityMenu, VALUE='LOS magn. field', $
                                UVALUE='VR_BLOS')
  vr_angleButton = widget_button(quantityMenu, VALUE='Field Angle', $
                                UVALUE='VR_ANGLE')
  IF (NOT tag_present(atmos, 'STOKES')) THEN $
   widget_control, vr_blosButton, SENSITIVE=0

  vr_contrButton = widget_button(quantityMenu, VALUE='Source function', $
                                UVALUE='VR_SOURCE')

  frame = widget_base(base, /ROW, /FRAME)

  ;; --- Draw frame --                                  ------------- ;;

  drawFrame = widget_base(base, /FRAME, /COLUMN)
  label = widget_label(drawFrame, VALUE="Display values along a ray", $
                       /ALIGN_CENTER)
  state.drawWidget = widget_draw(drawFrame, $
                                 XSIZE=screenSize[0], YSIZE=screenSize[1])

  slideFrame = widget_base(base, /FRAME, /ROW)
  state.Slider = widget_slider(slideFrame, TITLE='x index:', $
                               UVALUE='SLIDER', MAXIMUM=geometry.Nx - 1, $
                               VALUE=x_index, XSIZE=250)
  state.xlabel = widget_label(slideframe, /FRAME, $
                              VALUE=string(geometry.x[x_index]/1.0E6, $
                                           FORMAT='(F7.3)'))
  label = widget_label(slideframe, VALUE='[Mm]')

  state.waveslider = widget_slider(slideFrame, TITLE='wavelength:', $
                                   UVALUE='WAVESLIDER', $
                                   MAXIMUM=spectrum.Nspect - 1, $
                                   VALUE=0, XSIZE=250)
  state.wlabel = widget_label(slideframe, /FRAME, $
                              VALUE=string(spectrum.lambda[0], $
                                           FORMAT='(F9.3)'))
  label = widget_label(slideframe, VALUE='[nm]')
  widget_control, state.waveslider, SENSITIVE=0

  rayFrame = widget_base(frame, /ROW)
  raylabel = widget_label(rayFrame, VALUE='Ray:')
  state.rayText = widget_label(rayFrame, /FRAME, $
                                VALUE=string(FORMAT='(I2)', state.ray))

  xmulabel = widget_label(rayFrame, VALUE=' mu_x:')
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

  wmulabel      = widget_label(rayFrame, VALUE=' weight:')
  state.wmuText = widget_label(rayFrame, /FRAME, /DYNAMIC_RESIZE, $
                                VALUE=string(FORMAT='(E10.4)', $
                                             geometry.wmu(state.ray)))

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- rayWidgetSetup.pro ------- ;

PRO XViewRay, f, mu, x_index, GROUP_LEADER=group_leader, ZTITLE=ztitle

@spectrum.common
@files.common

  IF (NOT keyword_set(GROUP_LEADER)) THEN BEGIN
    group_leader = 0
    JFile    = 'J.dat'
    opacFile = 'opacity.out'
    backgroundFile = 'background.dat'
    IF (NOT openJ(JFile)) THEN return
  ENDIF

  IF (n_elements(spectrum) EQ 0) THEN $
   IF (NOT readSpectrum(spectrumFile)) THEN return
  IF (NOT openOpacity(opacFile)) THEN return
  
  IF (NOT keyword_set(ZTITLE)) THEN ztitle = "Function value"
  state = rayWidgetSetup(f, mu, x_index, ztitle)
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader
  displayRay, state

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewRay', state.baseWidget, $
   EVENT_HANDLER='XViewRay_Event', GROUP_LEADER=group_leader
END
