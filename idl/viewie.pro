; ----------------------------------------- viewie.pro --------------- ;

FUNCTION rayMenu_Event_Func, event

  widget_control, event.id, GET_UVALUE=uvalue

  return, {ID:event.handler, TOP:event.top, HANDLER:0L, VALUE:uvalue}
END
; -------- begin -------------------------- setRayNo.pro ------------- ;

FUNCTION setRayNo, stash, rayNo, OPLOT=oplot

@geometry.common

  widget_control, stash, GET_UVALUE=state
  state.ray = rayNo

  IF (NOT keyword_set(OPLOT)) THEN BEGIN
    widget_control, state.rayText, SET_VALUE=string(FORMAT='(I2)', rayNo)
    widget_control, state.xmuText, SET_VALUE=string(FORMAT='(F6.3)', $
                                                    geometry.xmu(state.ray))

    IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
      widget_control, state.ymuText, SET_VALUE=string(FORMAT='(F6.3)', $
						      geometry.ymu(state.ray))
      zmu = sqrt(1.0 - $
                 (geometry.xmu(state.ray)^2 +  geometry.ymu(state.ray)^2))
      widget_control, state.zmuText, SET_VALUE=string(FORMAT='(F6.3)', zmu)
    ENDIF

    widget_control, state.wmuText, SET_VALUE=string(FORMAT='(E10.4)', $
                                                    geometry.wmu(state.ray))
  ENDIF
  widget_control, stash, SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- setRayNo.pro ------------- ;

; -------- begin -------------------------- setBlue_n_Red.pro -------- ;

FUNCTION setBlue_n_Red, stash, lambdablue, lambdared

@spectrum.common

  widget_control, stash, GET_UVALUE=state
  IF (lambdared LE lambdablue) THEN BEGIN
    result = dialog_message(/ERROR, $
                            "Minimum wavelength should be lower than maximum")
    return, state
  ENDIF

  tabinv, spectrum.lambda, [lambdablue, lambdared], leff
  blue = long(leff(0) + 1)
  red  = long(leff(1)) < (spectrum.Nspect-1)

  IF (red EQ blue) THEN BEGIN
    result = dialog_message(/ERROR, $
                            "Found no wavelengths in specified interval")
    return, state
  ENDIF

  widget_control, state.blueText, $
   SET_VALUE=string(FORMAT='(I6)', blue)
  widget_control, state.redText, $
   SET_VALUE=string(FORMAT='(I6)', red)
  state.blue = blue
  state.red  = red 

  state.lambdablue = lambdablue  &  state.lambdared = lambdared
  widget_control, stash, SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- setBlue_n_Red.pro -------- ;

; -------- begin -------------------------- XViewIe_Event.pro -------- ;

PRO XViewIe_Event, Event

@files.common
@spectrum.common
@geometry.common

  COMMON screen_common, screenSize, scaleFactor

  Off = 0  &  On = 1

  ;; --- Main event handler --                          -------------- ;

  widget_control, Event.id, GET_UVALUE=Action
  stash = widget_info(Event.top, /CHILD)
  widget_control, stash, GET_UVALUE=state

  CASE Action OF
    'QUIT': widget_control, Event.top, /DESTROY

    'PRINT': BEGIN
      filename = 'viewIe-' + timeStamp() + '.ps'
      psopen, FILENAME=filename, /COLOR, FONT=font
      scaleFactor = [!D.X_SIZE, !D.Y_SIZE] / float(screenSize)
      drawIe, state
      psclose
      ScaleFactor = [1.0, 1.0]
      r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
    END

    'PNG': BEGIN
      widget_control, state.drawWidget, GET_VALUE=WindowNo
      rhwritepng, WindowNo, '/tmp/viewIe-'
    END

    'ORIENT': orient, state, EVENT_HANDLER='XViewIe_Event', $
     GROUP_LEADER=Event.top, TITLE='Orientation Intensity'
  
    'ORIENT_UPDATE': BEGIN
      OrientationUpdate, state
      drawIe, state
    END

    'ORIENT_RESET': BEGIN
      OrientationUpdate, state, /SET
      drawIe, state
    END

    'TRACK_SURFACE': drawIe, state, /TRACK

    'NEWSPECTRUMFILE': BEGIN
      IF (result = $
          readSpectrum(dialog_pickfile(FILTER='*.out', $
                                       TITLE='Emergent intensity file', $
                                       /MUST_EXIST, /READ, $
                                       FILE=spectrumFile))) THEN $
        drawIe, setBlue_n_Red(stash, state.lambdablue, state.lambdared) $
      ELSE $
       widget_control, Event.top, /DESTROY
    END

    'SETREDBLUE': BEGIN
      widget_control, state.blueField, GET_VALUE=lambdablue
      widget_control, state.redField,  GET_VALUE=lambdared
      drawIe, setBlue_n_Red(stash, lambdablue, lambdared)
    END

    'RAYMENU': drawIe, setRayNo(stash, Event.value)

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

    'BN_ATLAS': BEGIN
      state.atlas = Action
      drawIe, state
      widget_control, stash, SET_UVALUE=state
    END
    'KPNO_ATLAS': BEGIN
      state.atlas = Action
      drawIe, state
      widget_control, stash, SET_UVALUE=state
    END
    'KPIR_ATLAS': BEGIN
      state.atlas = Action
      drawIe, state
      widget_control, stash, SET_UVALUE=state
    END
    'ATMOS_ATLAS':  BEGIN
      state.atlas = Action
      drawIe, state
      widget_control, stash, SET_UVALUE=state
    END
    'HAWAII_ATLAS': BEGIN
      state.atlas = Action
      drawIe, state
      widget_control, stash, SET_UVALUE=state
    END
    'KPK_ATLAS': BEGIN
      state.atlas = Action
      drawIe, state
      widget_control, stash, SET_UVALUE=state
    END
    'KPSPOT_ATLAS': BEGIN
      state.atlas = Action
      drawIe, state
      widget_control, stash, SET_UVALUE=state
    END
    'SET_SCALE': BEGIN
      IF (Event.press AND Event.clicks EQ 2 AND state.atlas_toggle) THEN BEGIN
        state.x_press = Event.x
        state.y_press = Event.y
        drawIe, state, /SET_ATLAS_SCALE
        widget_control, stash, SET_UVALUE=state
      ENDIF
    END

    'ACTIVE_CONTRIB':  XViewContrib, BLUE=state.blue, RED=state.red, $
                           GROUP_LEADER=Event.top, /ACTIVE_SET
    'CONTIN_CONTRIB':  XViewContrib, BLUE=state.blue, RED=state.red, $
                           GROUP_LEADER=Event.top, /CONTINUUM
    'TOTAL_CONTRIB':   XViewContrib, BLUE=state.blue, RED=state.red, $
                           GROUP_LEADER=Event.top


    'SPLITTING': BEGIN
      atom_files = file_search("atom.*.out", COUNT=count)
      FOR n=0, count-1 DO begin
        atom = readAtom(atom_files[n], /ESSENTIALS)
        FOR kr=0, atom.Nline-1 DO BEGIN
          IF (atom.transition[kr].lambda0 GE state.lambdablue AND $
              atom.transition[kr].lambda0 LE state.lambdared) THEN BEGIN
            lo = atom.transition[kr].i
            hi = atom.transition[kr].j
            xViewSplitting, atom.g[hi], atom.labels[hi], $
                            atom.g[lo], atom.labels[lo], GROUP_LEADER=Event.top
          ENDIF
        ENDFOR
      ENDFOR
    END

    'STOKES_Q': XViewStokes, BLUE=state.blue, RED=state.red, $
     GROUP_LEADER=Event.top, 'Q', state.ray

    'STOKES_U': XViewStokes, BLUE=state.blue, RED=state.red, $
     GROUP_LEADER=Event.top, 'U', state.ray

    'STOKES_V': XViewStokes, BLUE=state.blue, RED=state.red, $
     GROUP_LEADER=Event.top, 'V', state.ray
    
    'LOG_ON': IF (Event.select) THEN $
     drawIe, setLog(stash, On)  ELSE  drawIe, setLog(stash, Off)

    'WIRE': IF (Event.select) THEN $
     drawIe, setWire(stash, On)  ELSE  drawIe, setWire(stash, Off)

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of emergent intensity. Go to the Rays menu to", $
             "view intensities from other directions. Display can be", $
             "set to linear or logarithmic intensity scale and displayed", $
             "as wiregrid or paneled graphs. The minimum (blue)", $
             "and maximum (red) wavelength of the display can be set in", $
             "the lower right of the panel.", $
             "", $
             "Version 1.0, Jun 2, 1995", $
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
; -------- end ---------------------------- XViewIe_Event.pro --------- ;

; -------- begin -------------------------- IeWidgetSetup.pro --------- ;

FUNCTION IeWidgetSetup, lambda, blue, red

@geometry.common
@spectrum.common
@atmos.common
@files.common
@input.common

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
           xmuText: 0L, ymuText: 0L, zmuText: 0L, wmuText: 0L, $
           blueField: 0L, redField: 0L, xSlider: 0L, zSlider: 0L, $
           ray: rayNo, log: 0, wire: 1,  blue: long(blue),  red: long(red), $
           lambdablue: lambda[blue], lambdared: lambda[red], $
           atlas_toggle: 0L, atlas: "KPNO", atlas_button: lonarr(7), $
           atlas_scale: 0.0, x_press: 0L, y_press: 0L}

  state.baseWidget = widget_base(TITLE='XViewIe', /ROW, MBAR=menuBar, $
                                  RESOURCE_NAME='XViewIe')

  spectrumBase = widget_base(state.baseWidget, /COLUMN)

  fileMenu  = widget_button(menuBar, VALUE='File', /MENU)
  openSpectrumButton = widget_button(fileMenu, $
                                      VALUE='open Spectrum file', $
                                      UVALUE='NEWSPECTRUMFILE')
  quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                             RESOURCE_NAME='quitbutton')

  rayMenu = widget_button(menuBar, Value='Rays', $
                           UVALUE='RAYMENU', $
                           /MENU, EVENT_FUNC='rayMenu_Event_Func')

  FOR i=0,geometry.Nrays-1 DO BEGIN
    menuItem = widget_button(rayMenu, UVALUE=i, $
                              VALUE=string(FORMAT='("ray ", I2)', i))
  ENDFOR

  IF (geometryType EQ "ONE_D_PLANE") THEN BEGIN
    contribMenu = widget_button(menuBar, VALUE='Contribution', /MENU)
    button1 = widget_button(contribMenu, VALUE='total', $
                           UVALUE='TOTAL_CONTRIB')
    button2 = widget_button(contribMenu, VALUE='active set', $
                           UVALUE='ACTIVE_CONTRIB')
    button3 = widget_button(contribMenu, VALUE='continuum', $
                           UVALUE='CONTIN_CONTRIB')

    f = file_search("asrs.out", COUNT=count)
    IF (count EQ 0) THEN BEGIN
      widget_control, button1, SENSITIVE=0
      widget_control, button2, SENSITIVE=0
    ENDIF
  ENDIF

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    orientMenu = widget_button(menuBar, VALUE='Orientation', /MENU)
    orientButton = widget_button(orientMenu, VALUE='set orientation', $
                                  UVALUE='ORIENT')
    trackButton = widget_button(orientMenu, VALUE='track surface', $
                                  UVALUE='TRACK_SURFACE')
  ENDIF ELSE BEGIN
    atlasMenu  = widget_button(menuBar, VALUE='Atlas', /MENU)
    IF (strlen(getenv('RH_ATLAS_PATH')) GT 0) THEN BEGIN
      button = widget_button(atlasMenu, VALUE='toggle atlas', $
                             UVALUE='ATLAS_TOGGLE')

      state.atlas_button[0] = $
       widget_button(atlasMenu, VALUE='Kitt Peak atlas (optical)', $
                     UVALUE='KPNO_ATLAS', /SEPARATOR)

      state.atlas_button[1] = widget_button(atlasMenu, $
                                            VALUE='Brault-Neckel', $
                                            UVALUE='BN_ATLAS')
      state.atlas_button[2] = widget_button(atlasMenu, $
                                            VALUE='Kitt Peak (IR)', $
                                            UVALUE='KPIR_ATLAS')
      state.atlas_button[3] = widget_button(atlasMenu, VALUE='Hawaii (UV)', $
                                            UVALUE='HAWAII_ATLAS')
      state.atlas_button[4] = widget_button(atlasMenu, VALUE='Harvard (UV)', $
                                            UVALUE='KPK_ATLAS')
      state.atlas_button[5] = widget_button(atlasMenu, VALUE='ATMOS (IR)', $
                                            UVALUE='ATMOS_ATLAS')
      state.atlas_button[6] = widget_button(atlasMenu, VALUE='KP spot', $
                                            UVALUE='KPSPOT_ATLAS')
      FOR n=0, 6 DO widget_control, state.atlas_button[n], SENSITIVE=0
    ENDIF ELSE $
     widget_control, atlasMenu, SENSITIVE=0
  ENDELSE

  stokesMenu = widget_button(menuBar, VALUE='Polarization', /MENU)
  f = file_search("atom.*.out", COUNT=count_as)
  IF (tag_present(atmos, 'BACKGRFLAGS')) THEN $
   f = where(atmos.backgrflags.ispolarized EQ 1, count_c) $
  ELSE $
   count_c = 0 
  
  IF (tag_present(atmos, 'STOKES') OR inputData.backgr_pol) THEN BEGIN
    IF (count_as GT 0) THEN $
     button = widget_button(stokesMenu, VALUE='Line splittings', $
                            UVALUE='SPLITTING')

    button = widget_button(stokesMenu, VALUE='Stokes Q', UVALUE='STOKES_Q')
    button = widget_button(stokesMenu, VALUE='Stokes U', UVALUE='STOKES_U')
    button = widget_button(stokesMenu, VALUE='Stokes V', UVALUE='STOKES_V')
  ENDIF ELSE $
   widget_control, stokesMenu, SENSITIVE=0

  ;; --- Print and help menus --

  printMenu= widget_button(menuBar, VALUE='Print', /MENU)
  printButton = widget_button(printMenu, VALUE='PostScript', UVALUE='PRINT')
  pngButton   = widget_button(printMenu, VALUE='png', UVALUE='PNG')

  helpMenu   = widget_button(menuBar, VALUE='Help', /MENU, /HELP)
  infoButton = widget_button(helpMenu, VALUE='XViewIe', $
                              UVALUE='INFORMATION')
  waveButton = widget_button(helpMenu, VALUE='Wavelength info', $
                              UVALUE='WAVEINFO')

  infoFrame = widget_base(spectrumBase, /COLUMN, /FRAME)
  waveFrame = widget_base(infoFrame, /ROW)
  waveLabel    = widget_label(waveFrame, VALUE='Wavelength #: ')
  state.blueText = widget_label(waveFrame, /FRAME, /DYNAMIC_RESIZE, $
                VALUE=string(FORMAT='(I6)', blue))
  waveLabel = widget_label(waveFrame, VALUE='to')
  state.redText = widget_label(waveFrame, /FRAME, /DYNAMIC_RESIZE, $
                 VALUE=string(FORMAT='(I6)', red))

  state.blueField = cw_field(waveFrame, TITLE='  Wavelength', /FLOATING, $
                             UVALUE='SETREDBLUE', VALUE=state.lambdablue, $
                             XSIZE=12, /RETURN_EVENT)
  state.redField = cw_field(waveFrame, TITLE='to', /FLOATING, $
                            UVALUE='SETREDBLUE', VALUE=state.lambdared, $
                            XSIZE=12, /RETURN_EVENT)
  waveLabel = widget_label(waveFrame, VALUE='[nm]')

  rayFrame = widget_base(infoFrame, /ROW)
  raylabel = widget_label(rayFrame, VALUE='Ray:')
  state.rayText = widget_label(rayFrame, /FRAME, $
                                VALUE=string(FORMAT='(I2)', state.ray))

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
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

  ;; --- Draw frame base widget --                      ------------- ;;

  drawFrame = widget_base(spectrumBase, /FRAME)
  state.drawWidget = widget_draw(drawFrame,$
                                 XSIZE=screenSize[0], YSIZE=screenSize[1], $
                                 /BUTTON_EVENTS, UVALUE='SET_SCALE')

  makeupFrame = widget_base(spectrumBase, /ROW, /FRAME)
  linlogFrame = widget_base(makeupFrame, /ROW, /EXCLUSIVE)
  linearButton = widget_button(linlogFrame, VALUE='Linear', $
                               UVALUE='LOG_OFF')
  logButton    = widget_button(linlogFrame, VALUE='Logarithmic', $
                               UVALUE='LOG_ON')
  widget_control, (state.log) ? logButton : linearButton, /SET_BUTTON
  
  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    shadeFrame   = widget_base(makeupFrame, /ROW, /EXCLUSIVE)
    wireButton   = widget_button(shadeFrame, VALUE='Wire Grid', $
                                 UVALUE='WIRE')
    panelButton  = widget_button(shadeFrame, VALUE='Panel', $
                                 UVALUE='PANEL')
    widget_control, (state.wire) ? wireButton : panelButton, /SET_BUTTON
  ENDIF

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- IeWidgetSetup.pro ------ ;

; -------- begin -------------------------- drawIe.pro ------------- ;

PRO drawIe, state, SET_ATLAS_SCALE=set_atlas_scale, TRACK=track

@geometry.common
@spectrum.common

  IF (!D.NAME EQ 'X') THEN BEGIN
    widget_control, state.drawWidget, GET_VALUE=WindowNo
    wset, WindowNo
  ENDIF

  blue = state.blue
  red  = state.red

  pxsize = 415
  pysize = 375

  widget_control, state.blueText, $
   SET_VALUE=string(FORMAT='(I6)', blue)
  widget_control, state.redText, $
   SET_VALUE=string(FORMAT='(I6)', red)

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    IF (keyword_set(TRACK)) THEN BEGIN
      tracktitle = 'emergent intensity for ray ' + $
       strtrim(string(state.ray, FORMAT='(I3)'), 2)
      IF (state.log) THEN tracktitle = 'log of ' + tracktitle

      surf_track, spectrum.I[blue:red, *, state.ray], ZLOG=state.log, $
       spectrum.lambda[blue:red], geometry.x/1.0E3, TITLE=tracktitle
    ENDIF ELSE BEGIN
      ztitle = 'Intensity [J m!U-2!N s!U-1!N Hz!U-1!N sr!U-1!N]'    
      IF (state.log) THEN  zlog = 1  ELSE  zlog = 0
      IF (NOT state.wire) THEN BEGIN
        erase
        panel, scaleimg_idl(spectrum.I[blue:red, *, state.ray], $
                            pxsize, pysize), $
               spectrum.lambda[blue:red], geometry.x/1.0E3, $
               /ORTHOSCOPIC, XTITLE='lambda [nm]', YTITLE='x [km]', $
               ZLOG=zlog, XPOS=65, YPOS=50, SCALETEXT=ztitle, $
               XTHICK=3, YTHICK=3

        avgsp = avg(spectrum.I[blue:red, *, state.ray], 1)
        spmin = min(avgsp, MAX=spmax)
        yrange = !Y.CRANGE[1] - !Y.CRANGE[0]
        avgspectrum = 0.1*yrange + 0.8*yrange * (avgsp / spmax)
        oplot, spectrum.lambda[blue:red], avgspectrum, THICK=3, COLOR=255B
      ENDIF ELSE BEGIN
        surface, spectrum.I[blue:red, *, state.ray], $
         spectrum.lambda[blue:red], geometry.x/1.0E3, $
         /T3D, CHARSIZE=1.4, BOTTOM=175B, $
         XTITLE='lambda [nm]', YTITLE='x [km]', ZTITLE=ztitle, ZLOG=zlog, $
         XRANGE=[state.lambdablue, state.lambdared], XSTYLE=1
      ENDELSE
    ENDELSE
  ENDIF ELSE BEGIN
    center_spectrum = file_search('spectrum_1.00', COUNT=count)
    IF (count GT 0) THEN BEGIN
      ray = readray(center_spectrum[0])
      Imax = max([ray.I[blue:red], spectrum.I[blue:red, state.ray]], MIN=Imin) 
    ENDIF ELSE $
      Imax = max(spectrum.I[blue:red, state.ray], MIN=Imin)

    IF (!D.NAME EQ 'PS') THEN thick = 4 ELSE thick = 1

    IF (state.log) THEN  ylog = 1  ELSE  ylog = 0
    plot, /NODATA, /YNOZERO, spectrum.lambda[blue:red], $
     spectrum.I[blue:red, state.ray], YLOG=ylog, $
     XTITLE='Wavelength [nm]', XMARGIN=[14, 2], $
     YTITLE='Intensity [J m!U-2!N s!U-1!N Hz!U-1!N sr!U-1!N]', $
     XRANGE=[state.lambdablue, state.lambdared], XSTYLE=1, $
     YRANGE=[Imin, Imax]
    FOR mu=0, geometry.Nrays-1 DO $
     IF (mu NE state.ray) THEN oplot, spectrum.lambda[blue:red], $
     spectrum.I[blue:red, mu], COLOR=130B, THICK=thick
    oplot, spectrum.lambda[blue:red], spectrum.I[blue:red, state.ray]

    IF (count GT 0) THEN BEGIN
      oplot, spectrum.lambda[blue:red], ray.I[blue:red], COLOR=175B
      rhannotate, xann(0.05), yann(0.1), TEXT='Disk center', MARKTYPE=0, $
       MARKCOLOR=175B, THICK=thick
    ENDIF ELSE ray = 0

    IF (state.atlas_toggle) THEN BEGIN
      atlasColor = 15B

      CASE (state.atlas) OF
        'HAWAII_ATLAS': BEGIN
          Hawaiiatlas, state.lambdablue, state.lambdared, Iatlas, latlas
          state.atlas_scale = 1.0
        END
        'KPK_ATLAS': BEGIN
          KPKatlas, state.lambdablue, state.lambdared, Iatlas, latlas
          state.atlas_scale = 1.0
        END
        'BN_ATLAS': BEGIN
          braultneckel, state.lambdablue, state.lambdared, Iatlas, latlas
          state.atlas_scale = 1.0
        END
        ELSE: BEGIN
          IF (n_tags(ray) GT 0 $
              OR geometryType EQ "SPHERICAL_SYMMETRIC") THEN BEGIN
            IF (state.atlas EQ "KPNO_ATLAS") THEN $
             KPNOatlas, state.lambdablue, state.lambdared, Iatlas, latlas $
            ELSE IF (state.atlas EQ "KPIR_ATLAS") THEN $
             KPIRatlas, state.lambdablue, state.lambdared, Iatlas, latlas $
            ELSE IF (state.atlas EQ "ATMOS_ATLAS") THEN $
             ATMOSatlas, state.lambdablue, state.lambdared, Iatlas, latlas $
            ELSE $
             KPspotatlas, state.lambdablue, state.lambdared, Iatlas, latlas

            IF (state.atlas_scale EQ 0.0 AND $
                NOT keyword_set(SET_ATLAS_SCALE)) THEN BEGIN
              r = dialog_message(/INFORMATION, $
                   "Double click to mark reference level for atlas", $
                                 DIALOG_PARENT=state.basewidget)
            ENDIF

            IF (keyword_set(SET_ATLAS_SCALE)) THEN BEGIN
              xc = convert_coord([state.x_press], [state.y_press], $
                                 /DEVICE, /TO_DATA)
              tabinv, latlas, xc[0], x_eff
              yatlas = interpolate(Iatlas, x_eff)
              tabinv, spectrum.lambda[blue:red], xc[0], x_eff
              IF (geometryType EQ "SPHERICAL_SYMMETRIC") THEN $
               ydata = interpolate(spectrum.I[blue:red, geometry.Nrays-1], $
                                   x_eff, /CUBIC) $
              ELSE $
               ydata = interpolate(ray.I[blue:red], x_eff, /CUBIC)
              state.atlas_scale = ydata[0] / yatlas[0]
            ENDIF
          ENDIF ELSE $
           r = dialog_message("Need disk center intensities for" + $
                                "this operation")
        END 
      ENDCASE 

      IF (n_elements(latlas) GT 1 AND state.atlas_scale NE 0.0) THEN BEGIN
        oplot, latlas, Iatlas*state.atlas_scale, $
         COLOR=atlasColor, THICK=thick
      ENDIF
    ENDIF
  ENDELSE
END
; -------- end ---------------------------- drawIe.pro --------------- ;

; -------- begin -------------------------- XViewIe.pro -------------- ;

PRO XViewIe, BLUE=blue, RED=red, GROUP_LEADER=group_leader

;+
; NAME:
;	XVIEWIE
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
;   --- Last modified: Tue Jun  7 16:42:23 2011 --
;-

@geometry.common
@atmos.common
@spectrum.common
@files.common

  IF (NOT keyword_set(GROUP_LEADER)) THEN group_leader=0

  IF (n_elements(atmos) EQ 0) THEN $
   IF (NOT readAtmos(atmosFile)) THEN return
  IF (n_elements(spectrum) EQ 0) THEN $
   IF (NOT readSpectrum(spectrumFile)) THEN return

  IF (NOT keyword_set(BLUE)) THEN blue = 0
  IF (NOT keyword_set(RED))  THEN red  = spectrum.Nspect - 1

  state = IeWidgetSetup(spectrum.lambda, blue, red)
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader
  surfr, AX=30, AZ=30
  drawIe, state

  ;; --- Register with the XManager --                   ------------ ;;

  xmanager, 'XViewIe', state.baseWidget, $
   EVENT_HANDLER='XViewIe_Event', GROUP_LEADER=group_leader
END
; -------- end ---------------------------- XViewIe.pro -------------- ;
