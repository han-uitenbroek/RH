; ----------------------------------------- viewstokes.pro ----------- ;

; -------- begin -------------------------- XViewStokes_Event.pro ---- ;

PRO XViewStokes_Event, Event

@files.common
@spectrum.common

  Off = 0  &  On = 1

  ;; --- Main event handler --                          -------------- ;

  widget_control, Event.id, GET_UVALUE=Action
  stash = widget_info(Event.top, /CHILD)
  widget_control, stash, GET_UVALUE=state

  CASE Action OF
    'QUIT': widget_control, Event.top, /DESTROY

    'PRINT': BEGIN
common SCREEN_COMMON, ScreenSize, ScaleFactor

      filename = '/tmp/viewStokes-' + timeStamp() + '.ps'
      psopen, FILENAME=filename, /COLOR, /LANDSCAPE
      ScreenSize = [600, 450]
      ScaleFactor = [!D.X_SIZE, !D.Y_SIZE]/float(ScreenSize)

      drawStokes, state
      psclose
      r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
    END

    'PNG': BEGIN
      widget_control, state.drawWidget, GET_VALUE=WindowNo
      rhwritepng, WindowNo, '/tmp/viewStokes-'
    END

    'SETREDBLUE': BEGIN
      widget_control, state.blueField, GET_VALUE=lambdablue
      widget_control, state.redField,  GET_VALUE=lambdared
      drawStokes, setBlue_n_Red(stash, lambdablue, lambdared)
    END

    'RAYMENU': drawStokes, setRayNo(stash, Event.value)

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
    
    'WIRE': IF (Event.select) THEN $
     drawStokes, setWire(stash, On)  ELSE  drawStokes, setWire(stash, Off)

    'RELATIVE': IF (Event.select) THEN $
     drawStokes, setRelative(stash, On) $
    ELSE  $
     drawStokes, setRelative(stash, Off)

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of Stokes intensities. Go to the Rays menu to", $
             "view intensities from other directions. The minimum (blue)", $
             "and maximum (red) wavelength of the display can be set in", $
             "the lower right of the panel.", $
             "", $
             "Version 1.0, Apr 2 2000", $
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
; -------- end ---------------------------- XViewStokes_Event.pro --------- ;

; -------- begin -------------------------- StokesWidgetSetup.pro --------- ;

FUNCTION StokesWidgetSetup, lambda, blue, red, StokesType, rayNo

@geometry.common
@spectrum.common
@atmos.common
@files.common

  Nspect = n_elements(lambda)

  state = {baseWidget: 0L, drawWidget: 0L, $
           blueText: 0L, redText: 0L, rayText: 0L, $
           xmuText: 0L, ymuText: 0L, zmuText: 0L, wmuText: 0L, $
           blueField: 0L, redField: 0L, xSlider: 0L, zSlider: 0L, $
           wire: 1, ray: rayNo, shade: 0,  blue: long(blue),  red: long(red), $
           lambdablue: lambda[blue], lambdared: lambda[red], $
           StokesType: StokesType, relative: 0L}

  state.baseWidget = widget_base(TITLE='XViewStokes: ' + StokesType, $
                                 /ROW, MBAR=menuBar, $
                                 RESOURCE_NAME='XViewStokes')

  StokesBase = widget_base(state.baseWidget, /COLUMN)

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

  stokesMenu = widget_button(menuBar, VALUE='Polarization', /MENU)

  button = widget_button(stokesMenu, VALUE='Line splittings', $
                         UVALUE='SPLITTING')
  button = widget_button(stokesMenu, VALUE='Stokes Q', UVALUE='STOKES_Q')
  button = widget_button(stokesMenu, VALUE='Stokes U', UVALUE='STOKES_U')
  button = widget_button(stokesMenu, VALUE='Stokes V', UVALUE='STOKES_V')

  ;; --- Print and help menus --

  printMenu= widget_button(menuBar, VALUE='Print', /MENU)
  printButton = widget_button(printMenu, VALUE='PostScript', UVALUE='PRINT')
  pngButton   = widget_button(printMenu, VALUE='png', UVALUE='PNG')

  helpMenu   = widget_button(menuBar, VALUE='Help', /MENU, /HELP)
  infoButton = widget_button(helpMenu, VALUE='XViewStokes', $
                              UVALUE='INFORMATION')
  waveButton = widget_button(helpMenu, VALUE='Wavelength info', $
                              UVALUE='WAVEINFO')

  infoFrame = widget_base(StokesBase, /COLUMN, /FRAME)
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

  drawFrame = widget_base(StokesBase, /FRAME)
  state.drawWidget = widget_draw(drawFrame, XSIZE=600, YSIZE=450, $
                                 /BUTTON_EVENTS, UVALUE='SET_SCALE')

  makeupFrame = widget_base(StokesBase, /ROW, /FRAME)
  relativeFrame = widget_base(makeupFrame, /ROW, /EXCLUSIVE)
  absoluteButton = widget_button(relativeFrame, VALUE='Absolute', $
                                 UVALUE='ABSOLUTE')
  relativeButton = widget_button(relativeFrame, VALUE='Relative', $
                                 UVALUE='RELATIVE')
  widget_control, (state.relative) ? $
   relativeButton : absoluteButton, /SET_BUTTON

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
; -------- end ---------------------------- StokesWidgetSetup.pro ------ ;

; -------- begin -------------------------- drawStokes.pro ------------- ;

PRO drawStokes, state

@geometry.common
@spectrum.common

  IF (!D.NAME EQ 'X') THEN BEGIN
    widget_control, state.drawWidget, GET_VALUE=WindowNo
    wset, WindowNo
    thick = 1
  ENDIF ELSE $ 
   thick = 4

  blue = state.blue
  red  = state.red

  pxsize = 415
  pysize = 375

  widget_control, state.blueText, $
   SET_VALUE=string(FORMAT='(I6)', blue)
  widget_control, state.redText, $
   SET_VALUE=string(FORMAT='(I6)', red)

  center_spectrum = file_search('spectrum_*1.00', COUNT=count)
  IF (count GT 0) THEN ray = readray(center_spectrum[0])

   CASE (state.StokesType) OF
    'Q': BEGIN
      CASE (geometryType) OF
        "ONE_D_PLANE": BEGIN
          stokes = spectrum.Stokes_Q[blue:red, *]
          IF (count GT 0) THEN stokes_ray = ray.Stokes_Q[blue:red]
        END
        "TWO_D_PLANE": BEGIN
          stokes = spectrum.Stokes_Q[blue:red, *, *]
          IF (count GT 0) THEN stokes_ray = ray.Stokes_Q[*, blue:red]
        END
      ENDCASE
    END
    'U': BEGIN
      CASE (geometryType) OF
        "ONE_D_PLANE": BEGIN
          stokes = spectrum.Stokes_U[blue:red, *]
          IF (count GT 0) THEN stokes_ray = ray.Stokes_U[blue:red]
        END
        "TWO_D_PLANE": BEGIN
          stokes = spectrum.Stokes_U[blue:red, *, *]
          IF (count GT 0) THEN stokes_ray = ray.Stokes_U[*, blue:red]
        END
      ENDCASE
    END
    'V': BEGIN
      CASE (geometryType) OF
        "ONE_D_PLANE": BEGIN
          stokes = spectrum.Stokes_V[blue:red, *]
          IF (count GT 0) THEN stokes_ray = ray.Stokes_V[blue:red]
        END
        "TWO_D_PLANE": BEGIN
          stokes = spectrum.Stokes_V[blue:red, *, *]
          IF (count GT 0) THEN stokes_ray = ray.Stokes_V[*, blue:red]
        END
      ENDCASE
    END
  ENDCASE
  IF (state.relative) THEN BEGIN
    CASE (geometryType) OF
      "ONE_D_PLANE": BEGIN
        stokes /= spectrum.I[blue:red, *]
        IF (count GT 0) THEN stokes_ray /= ray.I[blue:red]
      END
      "TWO_D_PLANE": BEGIN
       stokes /= spectrum.I[blue:red, *, *]
        IF (count GT 0) THEN stokes_ray /= ray.I[*, blue:red]
      END
    ENDCASE
    ytitle = 'Stokes ' + state.StokesType + ' / I'
  ENDIF ELSE $
   ytitle = 'Stokes ' + state.StokesType + $
   ' [J m!U-2!N s!U-1!N Hz!U-1!N sr!U-1!N]'

  CASE (geometryType) OF
    "ONE_D_PLANE": BEGIN
      plot, /NODATA, /YNOZERO, $
      spectrum.lambda[blue:red], stokes[*, state.ray], $
       XTITLE='Wavelength [nm]', XMARGIN=[14, 2], YTITLE=ytitle, $
       XRANGE=[state.lambdablue, state.lambdared], XSTYLE=1
      FOR mu=0, geometry.Nrays-1 DO $
       IF (mu NE state.ray) THEN oplot, spectrum.lambda[blue:red], $
       stokes[*, mu], COLOR=100, THICK=thick
      oplot, spectrum.lambda[blue:red], stokes[*, state.ray]

      IF (count GT 0) THEN BEGIN
        oplot, spectrum.lambda[blue:red], stokes_ray, COLOR=200B
        rhannotate, xann(0.05), yann(0.1), TEXT='Disk center', MARKTYPE=0, $
         MARKCOLOR=200B, THICK=thick
      ENDIF
    END
    "TWO_D_PLANE": BEGIN
      IF (NOT state.wire) THEN BEGIN
        erase
        panel, scaleimg_idl(stokes[*, *, state.ray], $
                            pxsize, pysize), $
               spectrum.lambda[blue:red], geometry.x/1.0E3, $
               /ORTHOSCOPIC, XTITLE='lambda [nm]', YTITLE='x [km]', $
               XPOS=65, YPOS=50, SCALETEXT=ytitle, $
               XTHICK=3, YTHICK=3

        avgsp = avg(stokes[*, *, state.ray], 1)
        spmin = min(avgsp, MAX=spmax)
        yrange = !Y.CRANGE[1] - !Y.CRANGE[0]
        avgspectrum = !Y.CRANGE[0] + 0.1*yrange + $
         0.8*yrange * (avgsp - spmin)/ (spmax - spmin)
        oplot, spectrum.lambda[blue:red], avgspectrum, COLOR=255B
      ENDIF ELSE BEGIN
        surface, XTITLE='Wavelength [nm]', YTITLE='x [km]', $
         ZTITLE=ytitle, XRANGE=[state.lambdablue, state.lambdared], XSTYLE=1, $
         stokes[*, *, state.ray], spectrum.lambda[blue:red], geometry.x/1.0E3
      ENDELSE
    END
    ELSE:
  ENDCASE

END
; -------- end ---------------------------- drawStokes.pro --------------- ;

; -------- begin -------------------------- XViewStokes.pro -------------- ;

PRO XViewStokes, BLUE=blue, RED=red, GROUP_LEADER=group_leader, $
                 StokesType, ray

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
;   --- Last modified: Fri Apr 17 15:30:27 2009 --
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

  state = StokesWidgetSetup(spectrum.lambda, blue, red, StokesType, ray)
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader
  drawStokes, state

  ;; --- Register with the XManager --                   ------------ ;;

  xmanager, 'XViewStokes', state.baseWidget, $
   EVENT_HANDLER='XViewStokes_Event', GROUP_LEADER=group_leader
END
; -------- end ---------------------------- XViewStokes.pro ---------- ;
