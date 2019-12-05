; -------- file: -------------------------- viewflux.pro ------------- ;

; -------- begin -------------------------- XViewFluxEvent.pro ------- ;

PRO XViewFluxEvent, Event

@spectrum.common
@geometry.common

  Off = 0  &  On = 1

  ;; --- Main event handler --                          -------------- ;

  widget_control, Event.id, GET_UVALUE=Action
  stash = widget_info(Event.top, /CHILD)
  widget_control, stash, GET_UVALUE=state

  CASE Action OF
    'QUIT': widget_control, Event.top, /DESTROY

    'NEWFLUXFILE': BEGIN
      IF (result = readSpectrum(dialog_pickfile(FILTER='*.out', $
                                                TITLE='Flux data File', $
                                                /MUST_EXIST, /READ, $
                                                FILE=spectrumFile))) THEN $
       displayFlux, state $
      ELSE $
       widget_control, Event.top, /DESTROY
    END

    'PRINT': BEGIN
      filename = '/tmp/viewFlux-' + timeStamp() + '.ps'
      IF (geometryType EQ "TWO_D_PLANE" OR $
          geometryType EQ "THREE_D_PLANE") THEN font = -1 ELSE font = 0
      psopen, FILENAME=filename, /COLOR, FONT=font
      displayFlux, state
      psclose
      r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
    END

    'PNG': BEGIN
      widget_control, state.drawWidget, GET_VALUE=WindowNo
      rhwritepng, WindowNo, '/tmp/viewFlux-'
    END

    'MARKEDGES': BEGIN
      t3d, /YZEXCH
      markEdges, state, /T3D
      t3d, /YZEXCH
    END

    'TOGGLE_ATLAS': BEGIN
      state.atlas = (state.atlas) ? 0L : 1L
      widget_control, stash, SET_UVALUE=state
      displayFlux, state
    END

    'ORIENT': orient, state, EVENT_HANDLER='XViewFluxEvent', $
     GROUP_LEADER=Event.top, TITLE='Orientation Flux'
  
    'ORIENT_UPDATE': BEGIN
      OrientationUpdate, state
      displayFlux, state
    END

    'ORIENT_RESET': BEGIN
      OrientationUpdate, state, /SET
      displayFlux, state
    END

    'SETREDBLUE': BEGIN
      widget_control, state.blueField, GET_VALUE=lambdablue
      widget_control, state.redField,  GET_VALUE=lambdared
      displayFlux, setBlue_n_Red(stash, lambdablue, lambdared)
    END

    'LOG_ON': IF (Event.select) THEN $
     displayFlux, setLog(stash, On)  ELSE  displayFlux, setLog(stash, Off)

    'SHADE_ON': IF (Event.select) THEN $
     displayFlux, setShade(stash, On) ELSE displayFlux, setShade(stash, Off)

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of Flux spectrum", $
             "", $
             "Version 1.0, Dec, 1996", $
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
; -------- end ---------------------------- XViewFluxEvent.pro ------- ;

; -------- begin -------------------------- displayFlux.pro ---------- ;

PRO displayFlux, state

@geometry.common
@spectrum.common

  Off = 0  &  On = 1

  IF (!D.NAME EQ 'X') THEN BEGIN
    widget_control, state.drawWidget, GET_VALUE=WindowNo
    wset, WindowNo
  ENDIF
  IF (!D.NAME EQ 'PS') THEN thick = 4 ELSE thick = 1

  blue = state.blue
  red  = state.red
  ztitle = 'Flux [J m!U-2!N s!U-1!N Hz!U-1!N]'

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN  
    IF (state.shade) THEN BEGIN
      shade_surf, flux[blue:red, *], $
       spectrum.lambda[blue:red], geometry.x/1.0E3, $
       /T3D, CHARSIZE=1.2, XTITLE='lambda [nm]', YTITLE='x [km]', $
       ZTITLE=ztitle, ZLOG=state.log, $
       XRANGE=[state.lambdablue, state.lambdared], XSTYLE=1
    ENDIF ELSE BEGIN
      surface, flux[blue:red, *], $
       spectrum.lambda[blue:red], geometry.x/1.0E3, $
       /T3D, CHARSIZE=1.2, BOTTOM=200B, $
       XTITLE='lambda [nm]', YTITLE='x [km]', ZTITLE=ztitle, ZLOG=state.log, $
       XRANGE=[state.lambdablue, state.lambdared], XSTYLE=1
    ENDELSE
  ENDIF
  IF (geometryType EQ "ONE_D_PLANE" OR $
      geometryType EQ "SPHERICAL_SYMMETRIC") THEN BEGIN
    plot, /YNOZERO, THICK=thick, spectrum.lambda[blue:red], flux[blue:red],$
     YLOG=state.log, XTITLE='Wavelength [nm]', $
     XMARGIN=[14,2], YTITLE=ztitle, $
     XRANGE=[state.lambdablue, state.lambdared], XSTYLE=1
    oplot, spectrum.lambda[blue:red], flux[blue:red], THICK=thick, $
     PSYM=4, COLOR=15B, SYMSIZE=0.7

;    oplot, spectrum.lambda[blue:red], $
;     !PI*Planck(5770.0, spectrum.lambda[blue:red], /HZ), $
;     COLOR=15B, THICK=thick

    IF (state.atlas) THEN BEGIN
      SolFluxatlas, spectrum.lambda[blue], spectrum.lambda[red], $
       atlasFlux, lambda
      IF (n_elements(atlasFlux) GT 1) THEN $
       oplot, lambda, atlasFlux, COLOR=130B, THICK=thick
    ENDIF
  ENDIF
END
; -------- end ---------------------------- displayFlux.pro ---------- ;

; -------- begin -------------------------- avgWidgetSetup.pro ------- ;

FUNCTION fluxWidgetSetup, lambda, blue, red

@geometry.common
@spectrum.common

  Nspect = n_elements(lambda)

  state = {baseWidget: 0L, drawWidget: 0L, $
           blueText: 0L, redText: 0L, blueField: 0L, redField: 0L, $
           blue: long(blue),  red: long(red),  log: 0,  shade: 0, $
           xSlider: 0L, zSlider: 0L, atlas: 0L, $
           lambdared: lambda[red], lambdablue: lambda[blue]}

  state.baseWidget = widget_base(TITLE='XViewFlux', /COLUMN, $
                                 RESOURCE_NAME='XViewFlux', MBAR=menuBar)

  avgBase = widget_base(state.baseWidget, /COLUMN)

  fileMenu   = widget_button(menuBar, VALUE='File', /MENU)
  openAtomButton  = widget_button(fileMenu, VALUE='open Spectrum file', $
                                   UVALUE='NEWFLUXFILE')
  quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                              RESOURCE_NAME='quitbutton')

  backgrMenu  = widget_button(menuBar, VALUE='Background', /MENU)
  edgesButton = widget_button(backgrMenu, VALUE='mark edges', $
                               UVALUE='MARKEDGES')

  IF (geometryType EQ "TWO_D_PLANE" OR $
      geometryType EQ "THREE_D_PLANE") THEN BEGIN
    orientMenu = widget_button( menuBar, VALUE='Orientation', /MENU )
    orientButton = widget_button( orientMenu, VALUE='set orientation', $
                                  UVALUE='ORIENT' )
  ENDIF ELSE BEGIN
    atlasMenu = widget_button(menuBar, VALUE='Atlas', /MENU)
    IF (strlen(getenv('RH_ATLAS_PATH')) GT 0) THEN BEGIN
      atlasButton = widget_button(atlasMenu, VALUE='toggle atlas', $
                                  UVALUE='TOGGLE_ATLAS')
    ENDIF ELSE $
     widget_control, atlasMenu, SENSITIVE=0
  ENDELSE

  printMenu= widget_button(menuBar, VALUE='Print', /MENU)
  printButton = widget_button(printMenu, VALUE='PostScript', UVALUE='PRINT')
  pngButton   = widget_button(printMenu, VALUE='png', UVALUE='PNG')

  helpMenu   = widget_button(menuBar, VALUE='Help', /MENU, /HELP)
  infoButton = widget_button(helpMenu, VALUE='XViewFlux', $
                              UVALUE='INFORMATION')
  waveButton = widget_button(helpMenu, VALUE='Wavelength info', $
                              UVALUE='WAVEINFO')

  waveFrame    = widget_base(avgBase, /ROW, /FRAME)
  waveLabel    = widget_label(waveFrame, VALUE='Wavelength #: ')
  state.blueText = widget_label(waveFrame, /FRAME, /DYNAMIC_RESIZE, $
                VALUE=string(FORMAT='(I6)', blue))
  waveLabel    = widget_label(waveFrame, VALUE=' to ')
  state.redText = widget_label(waveFrame, /FRAME, /DYNAMIC_RESIZE, $
                 VALUE=string(FORMAT='(I6)', red))
  state.blueField = cw_field(waveFrame, TITLE=' Wavelength:', /FLOATING, $
                             UVALUE='SETREDBLUE', VALUE=state.lambdablue, $
                             XSIZE=12, /RETURN_EVENT)
  state.redField  = cw_field(waveFrame, TITLE=' to ', /FLOATING, $
                             UVALUE='SETREDBLUE', VALUE=state.lambdared, $
                             XSIZE=12, /RETURN_EVENT)
  waveLabel = widget_label(waveFrame, VALUE='[nm]')

  ;; --- Draw frame --                                  ------------- ;;

  drawFrame =  widget_base(avgBase, /FRAME, /COLUMN)
  state.drawWidget = widget_draw(drawFrame, XSIZE=600, YSIZE=450)

  fieldFrame = widget_base(avgBase, /FRAME, /ROW)
  linlogFrame  = widget_base(fieldFrame, /ROW, /EXCLUSIVE)
  linearButton = widget_button(linlogFrame, VALUE='Lin', $
                               UVALUE='LOG_OFF')
  logButton    = widget_button(linlogFrame, VALUE='Log', $
                               UVALUE='LOG_ON')
  widget_control, (state.log) ? logButton : linearButton, /SET_BUTTON

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    shadeFrame   = widget_base(fieldFrame, /ROW, /EXCLUSIVE)
    wireButton   = widget_button(shadeFrame, VALUE='Wire Grid', $
                                 UVALUE='SHADE_OFF')
    shadeButton  = widget_button(shadeFrame, VALUE='Shade Surface', $
                                 UVALUE='SHADE_ON')
    widget_control, (state.shade) ? shadeButton : wireButton, /SET_BUTTON
  ENDIF

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- fluxWidgetSetup.pro ------ ;

; -------- begin -------------------------- XViewFlux.pro ------------ ;

PRO XViewFlux, BLUE=blue, RED=red, GROUP_LEADER=group_leader

;+
; NAME:
;	XVIEWFLUX
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
;   --- Last modified: Thu Jan 25 14:56:27 2001 --
;-
@atmos.common
@spectrum.common

  IF (NOT keyword_set(GROUP_LEADER)) THEN BEGIN
    group_leader = 0
    surfr, AX=30, AZ=30
  ENDIF

  IF (n_elements(atmos) EQ 0) THEN $
   IF (NOT readAtmos(atmosFile)) THEN return
  IF (n_elements(spectrum) EQ 0) THEN $
   IF (NOT readSpectrum(spectrumFile)) THEN return
  IF (n_elements(flux) EQ 0) THEN $
   IF (NOT readFlux('flux.out')) THEN return

  IF (NOT keyword_set(BLUE)) THEN  blue = 0
  IF (NOT keyword_set(RED))  THEN  red  = spectrum.Nspect - 1

  state = fluxWidgetSetup(spectrum.lambda, blue, red)
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader
  displayFlux, state

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewFlux', state.baseWidget, $
   EVENT_HANDLER='XViewFluxEvent', GROUP_LEADER=group_leader
END
; -------- end ---------------------------- XViewFlux.pro ------------ ;
