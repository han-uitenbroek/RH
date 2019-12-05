; -------- file: -------------------------- viewangles.pro ----------- ;

; -------- begin -------------------------- XViewAngles_Event.pro ---- ;

PRO XViewAngles_Event, Event

@files.common

  Line = 0  &  Cont = 1

  ;; --- Main event handler --                          -------------- ;

  stash = widget_info( Event.handler, /CHILD )
  widget_control, stash, GET_UVALUE=state

  widget_control, Event.id, GET_UVALUE=Action
  CASE Action OF

    'QUIT': widget_control, Event.top, /DESTROY

    'NEWGEOMETRYFILE': BEGIN
      IF (result = readGeometry(dialog_pickfile(FILTER='*.out', $
                                                TITLE='Geometry File', $
                                                /MUST_EXIST, /READ, $
                                                FILE=geometryFile))) THEN $
       drawRays, state $
      ELSE $
       widget_control, Event.top, /DESTROY
    END

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of geometry data", $
             "", $
             "Version 1.0, Jun 6, 1995", $
             "Han Uitenbroek (HUitenbroek@cfa.harvard.edu)"])
  ELSE:
  ENDCASE
END
; -------- end ---------------------------- XViewAngles_Event.pro ---- ;

; -------- begin -------------------------- set_T3D.pro -------------- ;

PRO set_T3D

  translation = [0.5, 0.5, 0.5]
  t3d, /RESET

  !X.S = [1.2, 1.0] / 2.4  &  !Y.S = !X.S
  !Z.S = [0.4, 1.0] / 1.6

  t3d, TRANSLATE=-translation
  t3d, /YZEXCH
  t3d, OBLIQUE=[0.5, 60]
  t3d, TRANSLATE=translation
END
; -------- end ---------------------------- set_T3D.pro -------------- ;

; -------- begin -------------------------- drawAxes.pro ------------- ;

PRO drawAxes

  erase
  Ncurve = 100  &  curveColor = 200B
  XZ = 1
  !X.TYPE = 0  &  !Y.TYPE = 0

  plots, [-1.0, 1.0], [0, 0], [0, 0], /T3D, THICK=2.0
  xyouts, 1.0, 0.0, Z=0.0, 'x', /T3D, TEXT_AXES=XZ, SIZE=2.0
  plots, [0, 0], [-2.0, 2.0], [0, 0], /T3D, COLOR=curveColor
  plots, [0, 0], [-1.0, 1.0], [0, 0], /T3D, THICK=2.0
  xyouts, 0.0, 1.0, Z=0.0, 'y', /T3D, TEXT_AXES=XZ, SIZE=2.0
  plots, [0, 0], [0, 0], [-0.2, 1.0], /T3D, THICK=2.0
  xyouts, 0.0, 0.0, Z=1.0, 'z', /T3D, TEXT_AXES=XZ, SIZE=2.0

  x = 2.0*findgen(Ncurve)/(Ncurve - 1) - 1.0
  curve = sqrt(1.0 - x^2)
  plots, x, fltarr(Ncurve), curve, /T3D, /DATA, COLOR=curveColor
  plots, fltarr(Ncurve), x, curve, /T3D, /DATA, COLOR=curveColor
END
; -------- end ---------------------------- drawAxes.pro ------------- ;

; -------- begin -------------------------- drawRays.pro ------------- ;

PRO drawRays, state

@geometry.common

  shade = 128B  &  rayColor = 208B  &  White = 255B
  XZ = 1

  widget_control, state.drawWidget, GET_VALUE=WindowNo
  wset, WindowNo
  widget_control, state.NraysLabel, SET_VALUE=string(FORMAT='(I2)', $
                                                     geometry.Nrays)
  IF (geometryType EQ "TWO_D_PLANE"  OR $
      geometryType EQ "THREE_D_PLANE") THEN BEGIN
    CASE (geometry.angleSet) OF
      0: text = 'SET_VERTICAL'
      1: text = 'SET_GL'
      2: text = 'SET_A2'
      3: text = 'SET_A4'
      4: text = 'SET_A6'
      5: text = 'SET_A8'
      6: text = 'SET_B4'
      7: text = 'SET_B6'
      8: text = 'SET_B8'
    ENDCASE
    widget_control, state.angleSetLabel, SET_VALUE=text

    drawAxes
    FOR i=0, geometry.Nrays-1 DO BEGIN
      xmu = geometry.xmu(i)  &  ymu = geometry.ymu(i)
      muz = sqrt(1.0 - (xmu^2 + ymu^2))
      plots, [0.0, xmu], [0.0, ymu], [0.0, muz], /T3D, /DATA, COLOR=rayColor, $
       THICK=2.0
      plots, [0.0, xmu], [0.0, ymu], [0.0, 0.0], /T3D, COLOR=shade
      plots, [xmu, xmu], [ymu, ymu], [0.0, muz], /T3D, COLOR=shade
    ENDFOR
    FOR i=0, geometry.Nrays-1 DO BEGIN
      xmu = geometry.xmu(i)  &  ymu = geometry.ymu(i)
      muz = sqrt(1.0 - (xmu^2 + ymu^2))
      xyouts, xmu/2, ymu/2, Z=muz/2, /T3D, TEXT_AXES=XZ, SIZE=1.4, $
       string(FORMAT='(I2)', i), COLOR=White, ALIGNMENT=0.5
      xyouts, xmu, ymu, Z=muz, /T3D, TEXT_AXES=XZ, SIZE=1.4, $
       string(FORMAT='(F5.3)', geometry.wmu(i)), COLOR=White
    ENDFOR

    surfr, AX=30, AZ=30
  ENDIF ELSE BEGIN
    plot, [0,1], XMARGIN=[2,12]  &  erase
    plots, [0.0, 1.0], [0.0, 0.0], COLOR=White
    plots, [0.0, 1.0], [0.8, 0.8], COLOR=White
    plots, [0.2, 0.2], [0.0, 1.0], COLOR=shade

    FOR mu=0, geometry.Nrays-1 DO BEGIN
      sinTheta = sqrt(1.0 - geometry.xmu(mu)^2)
      tanTheta = sinTheta / geometry.xmu(mu)
      arrow, 0.2, 0.0, 0.2 + tanTheta, 1.0, COLOR=rayColor, /DATA
     
      xyouts, 0.2 + 0.5*sinTheta, 0.5*geometry.xmu(mu), $
       string(FORMAT='(I2)', mu), COLOR=White, CHARSIZE=1.4, ALIGNMENT=0.5
    ENDFOR
  ENDELSE
END
; -------- end ---------------------------- drawRays.pro ------------- ;

; -------- begin -------------------------- anglesWidgetSetup.pro ---- ;

FUNCTION anglesWidgetSetup

@geometry.common

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    xsize = 400  &  ysize = 300
  ENDIF ELSE BEGIN
    xsize = 400  &  ysize = 400
  ENDELSE

  state = {baseWidget: 0L, drawWidget: 0L, NraysLabel: 0L, angleSetLabel: 0L}

  state.baseWidget = widget_base( TITLE='XViewAngles', /COLUMN, $
                                  RESOURCE_NAME='XViewAngles', MBAR=menuBar )
  geometryBase = widget_base( state.baseWidget, /COLUMN )

  fileMenu   = widget_button( menuBar, VALUE='File', /MENU)
  openAtomButton  = widget_button( fileMenu, VALUE='open Geometry file', $
                                   UVALUE='NEWGEOMETRYFILE' )
  quitButton = widget_button( fileMenu, VALUE='Quit', UVALUE='QUIT', $
                              RESOURCE_NAME='quitbutton')

  helpMenu   = widget_button( menuBar, VALUE='Help', /MENU, /HELP )
  infoButton = widget_button( helpMenu, VALUE='XViewAngles', $
                              UVALUE='INFORMATION' )

  ;; --- Draw frame base widget --                      ------------- ;;

  drawFrame     = widget_base( geometryBase, /FRAME, /COLUMN )
  labelFrame    = widget_base( drawFrame, /ROW )

  IF (geometryType EQ "TWO_D_PLANE"  OR $
      geometryType EQ "THREE_D_PLANE") THEN BEGIN
    angleSetLabel = widget_label( labelFrame, VALUE='Angle Set:' )
    state.angleSetLabel = widget_label( labelFrame, VALUE='SET_A2', /FRAME )
  ENDIF
  nraysLabel       = widget_label( labelFrame, VALUE='  Nrays:' )
  state.nraysLabel = widget_label( labelFrame, /FRAME, /DYNAMIC_RESIZE )
  state.drawWidget = widget_draw( drawFrame, XSIZE=xsize, YSIZE=ysize )

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- anglesWidgetSetup.pro ---- ;

; -------- begin -------------------------- XViewAngles.pro ---------- ;

PRO XViewAngles, GROUP_LEADER=group_leader

;+
; NAME:
;	XVIEWANGLES
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
;   --- Last modified: Thu Jan 15 16:35:20 2004   --
;-

@geometry.common

  IF (NOT keyword_set(GROUP_LEADER)) THEN group_leader=0
  IF (n_elements(geometry) EQ 0) THEN $
   IF (NOT readGeometry(geometryFile)) THEN return

  state = anglesWidgetSetup(  )
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader

  set_T3D
  drawRays, state

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewAngles', state.baseWidget, $
   EVENT_HANDLER='XViewAngles_Event', GROUP_LEADER=group_leader
END

; -------- end ---------------------------- XViewAngles.pro ---------- ;
