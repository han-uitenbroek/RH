; -------- file: -------------------------- viewavg.pro -------------- ;

; -------- begin -------------------------- XViewAvg_Event.pro ------- ;

PRO XViewAvg_Event, Event

@spectrum.common

  Off = 0  &  On = 1

  ;; --- Main event handler --                          -------------- ;

  stash = widget_info(Event.handler, /CHILD)
  widget_control, stash, GET_UVALUE=state

  widget_control, Event.id, GET_UVALUE=Action
  CASE Action OF
    'QUIT': widget_control, Event.top, /DESTROY

    'PRINT': BEGIN
      filename = '/tmp/viewAvg-' + timeStamp() + '.ps'
      psopen, FILENAME=filename, /COLOR
      displayAvg, state
      psclose
      r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
    END

    'PNG': BEGIN
      widget_control, state.drawWidget, GET_VALUE=WindowNo
      rhwritepng, WindowNo, '/tmp/viewAvg-'
    END

    'NEWSPECTRUMFILE': BEGIN
      IF (result = readSpectrum(dialog_pickfile(FILTER='*.out', $
                                                TITLE='Spectral data File', $
                                                /MUST_EXIST, /READ, $
                                                FILE=spectrumFile))) THEN $
       displayAvg, state $
      ELSE $
       widget_control, Event.top, /DESTROY
    END

    'RAYMENU': displayAvg, setRayNo(stash, Event.value, OPLOT=state.oplot)

    'MARKEDGES' : markEdges, state

    'SETREDBLUE': BEGIN
      widget_control, state.blueField, GET_VALUE=lambdablue
      widget_control, state.redField,  GET_VALUE=lambdared
      displayAvg, setBlue_n_Red(stash, lambdablue, lambdared)
    END

    'LOG_ON': IF (Event.select) THEN $
     displayAvg, setLog(stash, On)  ELSE  displayAvg, setLog(stash, Off)

    'OPLOT_ON': IF (Event.select) THEN $
      r = setOplot(stash, On)  ELSE  BEGIN
      r = setRayNo(stash, state.ray)
      displayAvg, setOplot(stash, Off)
    ENDELSE

    'ATLAS_TOGGLE': BEGIN
      IF (state.atlas_toggle) THEN BEGIN
        state.atlas_toggle = 0L
        state.atlas_scale  = 0.0
      ENDIF ELSE $
       state.atlas_toggle = 1L
      FOR n=0, n_elements(state.atlas_button)-1 DO $
       widget_control, state.atlas_button[n], SENSITIVE=state.atlas_toggle
      widget_control, stash, SET_UVALUE=state
    END
    'KPNO_ATLAS': BEGIN
      state.atlas = Action
      widget_control, stash, SET_UVALUE=state
      displayAvg, state
    END
    'BN_ATLAS': BEGIN
      state.atlas = Action
      widget_control, stash, SET_UVALUE=state
      displayAvg, state
    END
    'ATMOS_ATLAS':  BEGIN
      state.atlas = Action
      widget_control, stash, SET_UVALUE=state
      displayAvg, state
    END
    'KPIR_ATLAS':  BEGIN
      state.atlas = Action
      widget_control, stash, SET_UVALUE=state
      displayAvg, state
    END
    'HAWAII_ATLAS': BEGIN
      state.atlas = Action
      widget_control, stash, SET_UVALUE=state
      displayAvg, state
    END
    'SET_SCALE': BEGIN
      IF (Event.press AND Event.clicks EQ 2 AND state.atlas_toggle) THEN BEGIN
        state.x_press = Event.x
        state.y_press = Event.y
        displayAvg, state, /SET_ATLAS_SCALE
        widget_control, stash, SET_UVALUE=state
      ENDIF
    END

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of spatially averaged spectrum", $
             "", $
             "Version 1.0, Jun 26, 1995", $
             "Han Uitenbroek (HUitenbroek@cfa.harvard.edu)"])
    'WAVEINFO': BEGIN
      IF (spectrum.vacuum_to_air) THEN $
      result = dialog_message(/INFORMATION, $
       string(FORMAT='("Wavelengths at and above ", F8.3, ' + $
                     '" [nm] have been converted to air wavelengths")', $
              spectrum.air_limit)) $
      ELSE $
       result = dialog_message(/INFORMATION, "All wavelengths are in vacuum")
    END
  ELSE:
  ENDCASE

END
; -------- end ---------------------------- XViewAvg_Event.pro ------- ;

; -------- begin -------------------------- displayAvg.pro ----------- ;

PRO displayAvg, state, SET_ATLAS_SCALE=set_atlas_scale

@geometry.common
@spectrum.common

  COMMON oplot_save_theRays, theRays

  IF (!D.NAME EQ 'X') THEN BEGIN
    widget_control, state.drawWidget, GET_VALUE=WindowNo
    wset, WindowNo
  ENDIF

  IF (!D.NAME EQ 'PS') THEN thick = 4 ELSE thick = 1

  blue = state.blue
  red  = state.red
  Navg = red - blue + 1

  IF (NOT state.oplot) THEN $
   theRays = [state.ray] $
  ELSE BEGIN 
    idx = where(theRays EQ state.ray, count)
    IF (count EQ 0) THEN theRays = [theRays, state.ray]
  ENDELSE

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    center_spectrum = file_search('spectrum_0.00_1.00', COUNT=count)
  ENDIF ELSE BEGIN
    center_spectrum = file_search('spectrum_0.00_0.00', COUNT=count)
  ENDELSE

  FOR i=0, n_elements(theRays)-1 DO BEGIN
    avgSpectrum = dblarr(Navg)

    IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
      FOR la=0, Navg-1 DO $
       avgSpectrum[la] = $
       int_tabulated(geometry.x, $
                     cmprss(spectrum.I[blue+la, *, theRays[i]])) /$
       (geometry.x[geometry.Nx-1] - geometry.x[0])
    ENDIF

    IF (geometryType EQ "THREE_D_PLANE") THEN BEGIN
      FOR la=0, Navg-1 DO $
       avgSpectrum[la] = avg(spectrum.I[blue+la, *, *, theRays[i]])
    ENDIF

    IF (i EQ 0) THEN $
     plot, spectrum.lambda[blue:red], avgSpectrum, THICK=thick, $
     XTITLE='Wavelength [nm]', XMARGIN=[14, 2], YNOZERO=(NOT state.log), $
     YTITLE='Average Intensity [J m!U-2!N s!U-1!N Hz!U-1!N sr!U-1!N]', $
     XRANGE=[state.lambdablue, state.lambdared], XSTYLE=1, YLOG=state.log $
    ELSE BEGIN
      oplot, spectrum.lambda[blue:red], avgSpectrum, $
       COLOR=200B - i*32B, THICK=thick
      rhannotate, xann(0.85), yann(i*0.05), $
       TEXT=string(FORMAT='("ray ", I2)', theRays[i]), $
       MARKTYPE=0, MARKCOLOR=(200B - i*32B), THICK=thick
    ENDELSE

    oplot, spectrum.lambda[blue:red], avgSpectrum, $
     PSYM=5, SYMSIZE=0.4, COLOR=200B
  ENDFOR

  IF (count GT 0) THEN BEGIN
    ray = readray(center_spectrum[0])
    IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
      FOR la=0, Navg-1 DO $
       avgSpectrum[la] = int_tabulated(geometry.x, $
        cmprss(ray.I[*, blue+la])) / (geometry.x[geometry.Nx-1]-geometry.x[0])
     ENDIF ELSE $
      Iavg = avg(reform(ray.I[*, *, blue:red], Navg, $
                        geometry.Nx*geometry.Nz, /OVERWRITE), 1)
    oplot, spectrum.lambda[blue:red], avgSpectrum, COLOR=128B, THICK=thick
    rhannotate, xann(0.05), yann(0.1), TEXT='disk center', $
     MARKCOLOR=128B, MARKTYPE=0
  ENDIF
  IF (state.atlas_toggle) THEN BEGIN
    zmu = sqrt(1.0 - geometry.xmu[state.ray]^2 - geometry.ymu[state.ray]^2)
    IF (n_tags(ray) GT 0  OR  zmu EQ 1.0) THEN BEGIN
      atlasColor = 15B

      CASE (state.atlas) OF
        "KPNO_ATLAS": $
         KPNOatlas, state.lambdablue, state.lambdared, Iatlas, latlas
        "BN_ATLAS": BEGIN
          braultneckel, state.lambdablue, state.lambdared, Iatlas, latlas
          state.atlas_scale = 1.0
        END
        "ATMOS_ATLAS": $
         ATMOSatlas, state.lambdablue, state.lambdared, Iatlas, latlas
        "KPIR_ATLAS": $
         KPIRatlas, state.lambdablue, state.lambdared, Iatlas, latlas
        'HAWAII_ATLAS': BEGIN
          Hawaiiatlas, state.lambdablue, state.lambdared, Iatlas, latlas
          state.atlas_scale = 1.0
        END
      ENDCASE

      IF (state.atlas_scale EQ 0.0 AND $
          NOT keyword_set(SET_ATLAS_SCALE)) THEN BEGIN
        r = dialog_message(/INFORMATION, $
                           "Double click to mark reference level for atlas", $
                           DIALOG_PARENT=state.basewidget)
      ENDIF

      IF (keyword_set(SET_ATLAS_SCALE)) THEN BEGIN
        xc = convert_coord([state.x_press], [state.y_press], /DEVICE, /TO_DATA)
        tabinv, latlas, xc[0], x_eff
        yatlas = interpolate(Iatlas, x_eff)
        tabinv, spectrum.lambda[blue:red], xc[0], x_eff
        ydata = interpolate(avgSpectrum, x_eff, /CUBIC)
        state.atlas_scale = ydata[0] / yatlas[0]
      ENDIF

      IF (n_elements(latlas) GT 1 AND state.atlas_scale NE 0.0) THEN BEGIN
        oplot, latlas, Iatlas*state.atlas_scale, COLOR=atlasColor, THICK=thick
      ENDIF
    ENDIF ELSE $
     r = dialog_message("Need disk center intensities for this operation")
  ENDIF
END
; -------- end ---------------------------- displayAvg.pro ----------- ;

; -------- begin -------------------------- displayAvg.pro ----------- ;

PRO markEdges, state, T3D=t3d

@atmos.common
@spectrum.common

  labelColor = 255B
  IF (keyword_set(T3D)) THEN t3dset = 1  ELSE t3dset = 0

  widget_control, state.drawWidget, GET_VALUE=WindowNo
  wset, WindowNo
  usersym, [0,0], [-2,2], THICK=2, COLOR=100B

  edges = getEdges(LABELS=labels)
  index = where(edges GE state.lambdablue AND $
                edges LE state.lambdared, Nindex)

  IF (Nindex GT 0) THEN BEGIN
    edges  = edges[index]
    labels = labels[index]

    yann = !Y.CRANGE(0) + 0.7*(!Y.CRANGE(1) - !Y.CRANGE(0))
    dy = 0.02*(!Y.CRANGE(1) - !Y.CRANGE(0))
    IF (state.log) THEN BEGIN
      yann = 10^yann
      ylabel = yann*10^dy
    ENDIF ELSE ylabel = yann + dy
    oplot, edges, [replicate(yann, Nindex)], PSYM=8
    FOR n=0, Nindex-1 DO BEGIN
      xyouts, edges[n], ylabel, labels[n], $
       CHARSIZE=0.9, COLOR=labelColor, ORIENT=90.0
    ENDFOR
  ENDIF
END
; -------- end ---------------------------- displayAvg.pro ----------- ;

; -------- begin -------------------------- avgWidgetSetup.pro ------- ;

FUNCTION avgWidgetSetup, lambda, blue, red

@geometry.common
@spectrum.common

  Nspect = n_elements(lambda)

  state = {baseWidget: 0L, drawWidget: 0L, ray: 0, $
           rayText: 0L, log: 1, xmuText: 0L, ymuText: 0L, zmuText: 0L, $
	   wmuText: 0L, blueText: 0L, redText: 0L, $
	   blueField: 0L, redField: 0L, blue: long(blue),  red: long(red), $
           lambdared: lambda[red], lambdablue: lambda[blue], oplot: 0, $
           atlas_toggle: 0L, atlas: "KPNO", atlas_button: lonarr(5), $
           atlas_scale: 0.0, x_press: 0L, y_press: 0L}

  state.baseWidget = widget_base(TITLE='XViewAvg', /COLUMN, $
                                  RESOURCE_NAME='XViewAvg', $
                                  MBAR=menuBar)

  avgBase = widget_base(state.baseWidget, /COLUMN)

  fileMenu   = widget_button(menuBar, VALUE='File', /MENU)
  openAtomButton  = widget_button(fileMenu, VALUE='open Spectrum file', $
                                   UVALUE='NEWSPECTRUMFILE')
  quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                              RESOURCE_NAME='quitbutton')

  backgrMenu  = widget_button(menuBar, VALUE='Background', /MENU)
  edgesButton = widget_button(backgrMenu, VALUE='mark edges', $
                               UVALUE='MARKEDGES')

  rayMenu = widget_button(menuBar, Value='Rays', UVALUE='RAYMENU', $
                          /MENU, EVENT_FUNC='rayMenu_Event_Func')

  FOR i=0,geometry.Nrays-1 DO BEGIN
    menuItem = widget_button(rayMenu, UVALUE=i, $
                             VALUE=string(FORMAT='("Ray #", I2)', i))
  ENDFOR

  atlasMenu = widget_button(menuBar, VALUE='Atlas', /MENU)
  IF (strlen(getenv('RH_ATLAS_PATH')) GT 0) THEN BEGIN
    button = widget_button(atlasMenu, VALUE='toggle atlas', $
                           UVALUE='ATLAS_TOGGLE')
  ENDIF ELSE $
   widget_control, atlasMenu, SENSITIVE=0

  state.atlas_button[0] = $
   widget_button(atlasMenu, VALUE='KPNO atlas (optical)', $
                 UVALUE='KPNO_ATLAS', /SEPARATOR)
  state.atlas_button[1] = widget_button(atlasMenu, VALUE='Brault-Neckel', $
                                        UVALUE='BN_ATLAS')
  state.atlas_button[2] = widget_button(atlasMenu, VALUE='ATMOS (IR)', $
                                        UVALUE='ATMOS_ATLAS')
  state.atlas_button[3] = widget_button(atlasMenu, VALUE='KPIR', $
                                        UVALUE='KPIR_ATLAS')
  state.atlas_button[4] = widget_button(atlasMenu, VALUE='Hawaii (UV)', $
                                            UVALUE='HAWAII_ATLAS')
  FOR n=0, 4 DO widget_control, state.atlas_button[n], SENSITIVE=0

  printMenu= widget_button(menuBar, VALUE='Print', /MENU)
  printButton = widget_button(printMenu, VALUE='PostScript', UVALUE='PRINT')
  gifButton   = widget_button(printMenu, VALUE='png', UVALUE='PNG')

  helpMenu   = widget_button(menuBar, VALUE='Help', /MENU, /HELP)
  infoButton = widget_button(helpMenu, VALUE='XViewAvg', $
                              UVALUE='INFORMATION')
  waveButton = widget_button(helpMenu, VALUE='Wavelength info', $
                              UVALUE='WAVEINFO')

  infoFrame = widget_base(avgBase, /COLUMN, /FRAME)
  waveFrame = widget_base(infoFrame, /ROW)
  waveLabel    = widget_label(waveFrame, VALUE='Wavelength #: ')
  state.blueText = widget_label(waveFrame, /FRAME, /DYNAMIC_RESIZE, $
                VALUE=string(FORMAT='(I6)', blue))
  waveLabel    = widget_label(waveFrame, VALUE='to')
  state.redText = widget_label(waveFrame, /FRAME, /DYNAMIC_RESIZE, $
                 VALUE=string(FORMAT='(I6)', red))
  state.blueField = cw_field(waveFrame, TITLE=' Wavelength', /FLOATING, $
                             UVALUE='SETREDBLUE', VALUE=state.lambdablue, $
                             XSIZE=12, /RETURN_EVENT)
  state.redField  = cw_field(waveFrame, TITLE='to', /FLOATING, $
                             UVALUE='SETREDBLUE', VALUE=state.lambdared, $
                             XSIZE=12, /RETURN_EVENT)
  waveLabel    = widget_label(waveFrame, VALUE='[nm]')

  rayFrame = widget_base(infoFrame, /ROW)
  raylabel = widget_label(rayFrame, VALUE='Ray:')
  state.rayText = widget_label(rayFrame, /FRAME, $
                               VALUE=string(FORMAT='(I2)', state.ray))

  xmulabel      = widget_label(rayFrame, VALUE=' mu_x:')
  state.xmuText = widget_label(rayFrame, /FRAME, $
                               VALUE=string(FORMAT='(F6.3)', $
                                            geometry.xmu(state.ray)))
  ymulabel      = widget_label(rayFrame, VALUE=' mu_y:')
  state.ymuText = widget_label(rayFrame, /FRAME, $
                               VALUE=string(FORMAT='(F6.3)', $
                                            geometry.ymu(state.ray)))
  zmu = sqrt(1.0 - (geometry.xmu(state.ray)^2 +  geometry.ymu(state.ray)^2))
  zmulabel      = widget_label(rayFrame, VALUE=' mu_z:')
  state.zmuText = widget_label(rayFrame, /FRAME, $
                               VALUE=string(FORMAT='(F6.3)', zmu))

  wmulabel      = widget_label(rayFrame, VALUE=' weight:')
  state.wmuText = widget_label(rayFrame, /FRAME, $
                               VALUE=string(FORMAT='(E10.4)', $
                                            geometry.wmu(state.ray)))

  ;; --- Draw frame --                                  ------------- ;;

  drawFrame =  widget_base(avgBase, /FRAME, /COLUMN)
  state.drawWidget = widget_draw(drawFrame, XSIZE=600, YSIZE=450, $
                                 /BUTTON_EVENTS, UVALUE='SET_SCALE')

  fieldFrame = widget_base(avgBase, /FRAME, /ROW)

  linlogFrame  = widget_base(fieldFrame, /ROW, /EXCLUSIVE)
  linearButton = widget_button(linlogFrame, VALUE='Lin', $
                               UVALUE='LOG_OFF')
  logButton    = widget_button(linlogFrame, VALUE='Log', $
                               UVALUE='LOG_ON')
  widget_control, (state.log) ? logButton : linearButton, /SET_BUTTON

  oplotFrame = widget_base(fieldFrame, /ROW, /EXCLUSIVE)
  plotButton = widget_button(oplotFrame, VALUE='Plot', $
                             UVALUE='OPLOT_OFF')
  oplotButton = widget_button(oplotFrame, VALUE='Oplot', $
                              UVALUE='OPLOT_ON')
  widget_control, (state.oplot) ? oplotButton : plotButton, /SET_BUTTON

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- avgWidgetSetup.pro ------- ;

; -------- begin -------------------------- XViewAvg.pro ------------- ;

PRO XViewAvg, BLUE=blue, RED=red, GROUP_LEADER=group_leader

@geometry.common
@spectrum.common

  IF (NOT keyword_set(GROUP_LEADER)) THEN BEGIN
    group_leader = 0
    surfr, AX=30, AZ=30
  ENDIF

  IF (n_elements(geometry) EQ 0) THEN $
   IF (NOT readAtmos(geometryFile)) THEN return
  IF (n_elements(spectrum) EQ 0) THEN $
   IF (NOT readSpectrum(spectrumFile)) THEN return

  IF (NOT keyword_set(BLUE)) THEN  blue = 0
  IF (NOT keyword_set(RED))  THEN  red  = spectrum.Nspect - 1

  state = avgWidgetSetup(spectrum.lambda, blue, red)
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader
  displayAvg, state

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewAvg', state.baseWidget, $
   EVENT_HANDLER='XViewAvg_Event', GROUP_LEADER=group_leader
END
; -------- end ---------------------------- XViewAvg.pro ------------- ;
