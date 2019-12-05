; ----------------------------------------- viewj.pro ---------------- ;

; -------- begin -------------------------- setLambdaNo.pro ---------- ;

FUNCTION setLambdaNo, stash, index

@spectrum.common

  ;; --- This routine is shared by viewJ and viewsource --  ---------- ;

  widget_control, stash, GET_UVALUE=state

  lambdaDisplay  = index
  state.lambdaNo = index
  widget_control, state.lambdaText, $
   SET_VALUE=string(FORMAT='(F10.3)', spectrum.lambda[index])

  widget_control, stash, SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- setLambdaNo.pro ---------- ;

; -------- begin -------------------------- XViewJ_Event.pro --------- ;

PRO XViewJ_Event, Event

@files.common
@geometry.common
@spectrum.common

  Off = 0  &  On = 1

  ;; --- Main event handler --                          -------------- ;

  widget_control, Event.id, GET_UVALUE=Action
  stash = widget_info(Event.top, /CHILD)
  widget_control, stash, GET_UVALUE=state

  CASE Action OF
    'QUIT': widget_control, Event.top, /DESTROY

    'ORIENT': orient, state, EVENT_HANDLER='XViewJ_Event', $
     GROUP_LEADER=Event.top, TITLE='Orientation J'
  
    'ORIENT_UPDATE': BEGIN
      OrientationUpdate, state
      drawJ, state
    END

    'ORIENT_RESET': BEGIN
      OrientationUpdate, state, /SET
      drawJ, state
    END

    'SLICER': slicer3, ptr_new(reverse(J, 3))
    'SLICER_LOG': slicer3, ptr_new(alog10(reverse(J, 3)))

    'LAMBDA_SLIDER':  drawJ, setLambdaNo(stash, Event.value)

    'HEIGHT_SLIDER': BEGIN
      state.z = Event.value
      widget_control, stash, SET_UVALUE=state
      widget_control, state.heightText, $
       SET_VALUE=string(FORMAT='(F10.3)', geometry.z[state.z]/1.0E+03)
      drawJ, state
    END

    'LOG_ON': IF (Event.select) THEN $
     drawJ, setLog(stash, On)  ELSE  drawJ, setLog(stash, Off)

    'SHADE_ON': IF (Event.select) THEN $
     drawJ, setShade(stash, On)  ELSE  drawJ, setShade(stash, Off)

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of mean intensity", $
             "", $
             "", $
             "Version 1.0, Jun 2, 1995", $
             "Han Uitenbroek (HUitenbroek@cfa.harvard.edu)"])
  ELSE:
  ENDCASE

END
; -------- end ---------------------------- XViewJ_Event.pro --------- ;

; -------- begin -------------------------- WidgetSetup_J.pro -------- ;

FUNCTION JWidgetSetup, lambda, lambdaDisplay

@geometry.common

  state = {baseWidget: 0L, drawWidget: 0L, $
           lambdaText: 0L, lambdaSlider: 0L, xSlider: 0L, zSlider: 0L, $
           lambdaNo: long(lambdaDisplay), log: 1, shade: 0, $
           heightSlider: 0L, heightText: 0L, z: 0L}

  state.baseWidget = widget_base(TITLE='XViewJ', /ROW, MBAR=menuBar, $
                                 RESOURCE_NAME='XViewJ')

  JBase = widget_base(state.baseWidget, /COLUMN)

  fileMenu = widget_button(menuBar, VALUE='File', /MENU)
  openSpectrumButton = widget_button(fileMenu, VALUE='open J file', $
                                     UVALUE='NEWJFILE')
  quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                             RESOURCE_NAME='quitbutton')

  IF (geometryType EQ "THREE_D_PLANE") THEN BEGIN
    slicerMenu = widget_button(menuBar, VALUE='3-D Slicer')
    slicerButton = widget_button(slicerMenu, VALUE='slice J', UVALUE='SLICER')
    slicerButton = widget_button(slicerMenu, VALUE='slice log(J)', $
                                 UVALUE='SLICER_LOG')
  ENDIF

  IF (geometryType EQ "TWO_D_PLANE" OR $
      geometryType EQ "THREE_D_PLANE") THEN BEGIN
    orientMenu = widget_button(menuBar, VALUE='Orientation', /MENU)
    orientButton = widget_button(orientMenu, VALUE='set orientation', $
                                 UVALUE='ORIENT')
  ENDIF

  helpMenu   = widget_button(menuBar, VALUE='Help', /MENU, /HELP)
  infoButton = widget_button(helpMenu, VALUE='XViewJ', $
                              UVALUE='INFORMATION')

  ;; --- Draw frame base widget --                      ------------- ;;

  waveFrame = widget_base(Jbase, /ROW, /FRAME)
  state.lambdaSlider = widget_slider(waveFrame, TITLE='Wavelength No:', $
                                      UVALUE='LAMBDA_SLIDER', XSIZE=255, $
                                      MAXIMUM=n_elements(lambda) - 1, $
                                      VALUE=lambdaDisplay)
  waveLabel = widget_label(waveFrame, VALUE='  Wavelength: ')
  state.lambdaText = widget_label(waveFrame, /FRAME, /DYNAMIC_RESIZE)
  waveLabel = widget_label(waveFrame, VALUE='[nm]')

  drawFrame = widget_base(JBase, /FRAME)
  state.drawWidget = widget_draw(drawFrame, XSIZE=500, YSIZE=400)

  makeupFrame  = widget_base(JBase, /ROW, /FRAME)
  linlogFrame  = widget_base(makeupFrame, /ROW, /EXCLUSIVE)
  linearButton = widget_button(linlogFrame, VALUE='Linear', $
                                UVALUE='LOG_OFF')
  logButton    = widget_button(linlogFrame, VALUE='Logarithmic', $
                                UVALUE='LOG_ON')
  widget_control, (state.log) ? logButton : linearButton, /SET_BUTTON

  IF (geometryType EQ "TWO_D_PLANE" OR $
      geometryType EQ "THREE_D_PLANE") THEN BEGIN
    shadeFrame = widget_base(makeupFrame, /ROW, /EXCLUSIVE)
    wireButton = widget_button(shadeFrame, VALUE='Wire Grid', $
                                  UVALUE='SHADE_OFF')
    shadeButton = widget_button(shadeFrame, VALUE='Shade Surface', $
                                  UVALUE='SHADE_ON')
    widget_control, (state.shade) ? shadeButton : wireButton, /SET_BUTTON
  ENDIF
  IF (geometryType EQ "THREE_D_PLANE") THEN BEGIN
    heightFrame = widget_base(JBase, /ROW, /FRAME)
    state.heightSlider = widget_slider(heightFrame, TITLE='Depth index:', $
                                       UVALUE='HEIGHT_SLIDER', XSIZE=255, $
                                       MAXIMUM=geometry.Nz - 1, VALUE=state.z)

    label = widget_label(heightFrame, VALUE='  Height: ')
    state.heightText = widget_label(heightFrame, /FRAME, /DYNAMIC_RESIZE, $
                                    VALUE=string(geometry.z[state.z]/1.0E+3, $
                                                 FORMAT='(F10.3)'))
    label = widget_label(heightFrame, VALUE='[km]')
  ENDIF

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- WidgetSetup_J.pro ------ ;

; -------- begin -------------------------- drawJ.pro -------------- ;

PRO drawJ, state

@geometry.common
@atmos.common
@spectrum.common

  widget_control, state.drawWidget, GET_VALUE=WindowNo
  wset, WindowNo

  readJ, state.lambdaNo
  Jtitle = 'Mean Intensity [J m!U-2!N s!U-1!N Hz!U-1!N sr!U-1!N]'

  CASE geometryType OF
    'ONE_D_PLANE': BEGIN
      plot, geometry.cmass, J, /XLOG, YLOG=state.log, $
       XTITLE='Column Mass [kg m!U-2!N]', YTITLE=Jtitle
    END
    'TWO_D_PLANE': BEGIN
      IF (state.shade) THEN BEGIN
        shade_surf, J, geometry.x/1.0E3, geometry.z/1.0E3, $
         /T3D, CHARSIZE=1.4, ZLOG=state.log, $
         XTITLE='x [km]', YTITLE='z [km]', ZTITLE=Jtitle
      ENDIF ELSE BEGIN
        surface, J, geometry.x/1.0E3, geometry.z/1.0E3, $
         /T3D, CHARSIZE=1.4, BOTTOM=200B, ZLOG=state.log, $
         XTITLE='x [km]', YTITLE='z [km]', ZTITLE=ztitle
      ENDELSE
    END
    'THREE_D_PLANE': BEGIN
      IF (state.shade) THEN BEGIN
        shade_surf, J[*, *, state.z], $
         geometry.dx/1.0E3*findgen(geometry.Nx), $
         geometry.dy/1.0E3*findgen(geometry.Ny), $
         /T3D, CHARSIZE=1.4, ZLOG=state.log, $
         XTITLE='x [km]', YTITLE='y [km]', ZTITLE=Jtitle
      ENDIF ELSE BEGIN
        surface, J[*, *, state.z], $
         geometry.dx/1.0E3*findgen(geometry.Nx), $
         geometry.dy/1.0E3*findgen(geometry.Ny), $
         /T3D, CHARSIZE=1.4, BOTTOM=200B, ZLOG=state.log, $
         XTITLE='x [km]', YTITLE='y [km]', ZTITLE=ztitle
      ENDELSE
    END
    'SPHERICAL_SYMMETRIC': BEGIN
      plot, geometry.cmass, J, YTITLE=Jtitle, /XLOG, YLOG=state.log, $
       XTITLE='Column Mass [kg m!U-2!N]'
    END  
  ENDCASE
END
; -------- end ---------------------------- drawJ.pro ---------------- ;

; -------- begin -------------------------- openJ.pro ---------------- ;

FUNCTION openJ, fileName

@geometry.common
@spectrum.common
@files.common
@input.common

  WHILE (NOT existFile(fileName, UNIT=Junit, $
                       SWAP_ENDIAN=(inputdata.big_endian NE $
                                    is_big_endian()))) DO BEGIN
    answer = dialog_message(/QUESTION, $
                            ["File " + fileName + " does not exist", $
                             "Continue?"])
    IF (answer EQ 'Yes') THEN $
     fileName = dialog_pickfile(FILTER='*.dat', $
                                TITLE='Find mean intensity Data', $
                                /MUST_EXIST, /READ, FILE=fileName) $
    ELSE $
     return, 0
  ENDWHILE

  ;; --- Store name in common block Files_Common --     -------------- ;

  JFile = fileName

  CASE geometryType OF
    "ONE_D_PLANE":         J = dblarr(geometry.Ndep)
    "TWO_D_PLANE":         J = dblarr(geometry.Nx, geometry.Nz)
    "THREE_D_PLANE":       J = dblarr(geometry.Nx, geometry.Ny, geometry.Nz)
    "SPHERICAL_SYMMETRIC": J = dblarr(geometry.Nradius)
  ENDCASE
  return, 1
END
; -------- end ---------------------------- openJ.pro ---------------- ;

; -------- begin -------------------------- readJ.pro ---------------- ;

PRO readJ, lambda

@geometry.common
@spectrum.common
@files.common

  lambdaDisplay = lambda
  CASE geometryType OF
    "ONE_D_PLANE":         offset = 8UL*lambda * geometry.Ndep
    "TWO_D_PLANE":         offset = 8UL*long64(lambda) * $
                                    long64(geometry.Nx*geometry.Nz)
    "THREE_D_PLANE":       offset = 8UL*long64(lambda) * $
                                    (long64(geometry.Nx)*long64(geometry.Ny)* $
                                     long64(geometry.Nz))
    "SPHERICAL_SYMMETRIC": offset = 8UL*lambda * geometry.Nradius
  ENDCASE

  point_lun, Junit, offset
  readu, Junit, J
END
; -------- end ---------------------------- readJ.pro ---------------- ;

; -------- begin -------------------------- XViewJ.pro --------------- ;

PRO XViewJ, GROUP_LEADER=group_leader

;+
; NAME:
;	XVIEWJ
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
;   --- Last modified: Fri Jul  2 03:42:06 2010 --
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
    JFile = 'J.dat'
    IF (NOT openJ(JFile)) THEN return
  ENDIF

  state = JWidgetSetup(spectrum.lambda, lambdaDisplay)
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader
  surfr, AX=30, AZ=30

  drawJ, setLambdaNo(widget_info(state.baseWidget, /CHILD), lambdaDisplay)

  ;; --- Register with the XManager --               -------------- ;;

  xmanager, 'XViewJ', state.baseWidget, $
   EVENT_HANDLER='XViewJ_Event', GROUP_LEADER=group_leader

  IF (NOT keyword_set(GROUP_LEADER)) THEN  free_lun, Junit
END
; -------- end ---------------------------- XViewJ.pro --------------- ;
