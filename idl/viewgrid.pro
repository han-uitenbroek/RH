; ----------------------------------------- viewgrid.pro ------------- ;

; -------- begin -------------------------- XViewGrid_Event.pro ------ ;

PRO XViewGrid_Event, Event

  Off = 0  &  On = 1

  ;; --- Main event handler --                          -------------- ;

  stash = widget_info(Event.handler, /CHILD)
  widget_control, stash, GET_UVALUE=state

  widget_control, Event.id, GET_UVALUE=Action
  CASE Action OF
    'QUIT': widget_control, Event.top, /DESTROY

    'NEWATMOSFILE': BEGIN
      IF (result = readAtmos(dialog_pickfile(FILTER='*.dat', $
                                             TITLE='Atmospheric data File', $
                                             /MUST_EXIST, /READ, $
                                             FILE=atmosFile))) THEN $
       displayGrid, stash $
      ELSE $
       widget_control, Event.top, /DESTROY
    END

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of spatial grid", $
             "", $
             "Version 1.0, Jun 19, 1995", $
             "Han Uitenbroek (HUitenbroek@cfa.harvard.edu)"])
  ELSE:
  ENDCASE

END
; -------- end ---------------------------- XViewGrid_Event.pro ------ ;

; -------- begin -------------------------- displayGrid.pro ---------- ;

PRO displayGrid, stash

  COMMON screen_common, ScreenSize, ScaleFactor
@geometry.common
@atmos.common

  widget_control, stash, GET_UVALUE=state

  widget_control, state.NxText, $
   SET_VALUE=string(FORMAT='(I3)', geometry.Nx)
  widget_control, state.NzText, $
   SET_VALUE=string(FORMAT='(I3)', geometry.Nz)

  widget_control, state.drawWidget, GET_VALUE=WindowNo
  wset, WindowNo
  ScreenSize = state.ScreenSize  &  ScaleFactor = state.ScaleFactor

  panel, scaleimg_idl(atmos.T, 325, 275), xpos=60, ypos=50, $
   geometry.x/1.0E3, xtitle='x [km]', geometry.z/1.0E3, ytitle='z [km]', $
   scaleText='T [K]', /ORDER, /ORTHOSCOPIC
  
  for l=0, geometry.Nx-1 do $
   oplot, [1, 1]*geometry.x(l)/1.0E3,$
   [geometry.z(0), geometry.z(geometry.Nz-1)]/1.0E3, color=150
  for k=0, geometry.Nz-1 do $
   oplot, [geometry.x(0), geometry.x(geometry.Nx-1)]/1.0E3, $
   [1, 1]*geometry.z(k)/1.0E3, color=150
END
; -------- end ---------------------------- displayGrid.pro ---------- ;

; -------- begin -------------------------- gridWidgetSetup.pro ------ ;

FUNCTION gridWidgetSetup

  state = {baseWidget: 0L, drawWidget: 0L, $
           ScreenSize: fix([500, 400]), ScaleFactor: float([1.0, 1.0]), $
           NxText: 0L,  NzText: 0L}

  state.baseWidget = widget_base( TITLE='XViewGrid', /COLUMN, $
                                  RESOURCE_NAME='XViewGrid', $
                                  MBAR=menuBar )

  gridBase = widget_base( state.baseWidget, /COLUMN )

  fileMenu   = widget_button( menuBar, VALUE='File', /MENU)
  openAtomButton  = widget_button( fileMenu, VALUE='open Atmos file', $
                                   UVALUE='NEWATMOSFILE' )
  quitButton = widget_button( fileMenu, VALUE='Quit', UVALUE='QUIT', $
                              RESOURCE_NAME='quitbutton')

  helpMenu   = widget_button( menuBar, VALUE='Help', /MENU, /HELP )
  infoButton = widget_button( helpMenu, VALUE='XViewGrid', $
                              UVALUE='INFORMATION' )

  drawFrame =  widget_base( gridBase, /FRAME, /COLUMN )
  labelFrame = widget_base( drawFrame, /ROW )
  NxLabel      = widget_label( labelFrame, VALUE='  Nx:' )
  state.Nxtext = widget_label( labelFrame, /FRAME )
  NzLabel      = widget_label( labelFrame, VALUE='  Nz:' )
  state.Nztext = widget_label( labelFrame, /FRAME )

  ;; --- Draw frame --                                  ------------- ;;

  state.drawWidget = widget_draw( drawFrame, XSIZE=state.ScreenSize(0), $
                                  YSIZE=state.ScreenSize(1))

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- gridWidgetSetup.pro ------ ;

; -------- begin -------------------------- XViewGrid.pro ------------ ;

PRO XViewGrid, GROUP_LEADER=group_leader

;+
; NAME:
;	XVIEWGRID
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
;   --- Last modified: Thu Dec 16 17:44:26 1999 --
;-

@geometry.common
@atmos.common
@files.common

  IF (NOT keyword_set(GROUP_LEADER)) THEN group_leader=0

  IF (n_elements(atmos) EQ 0) THEN $
   IF ( (NOT readGeometry(geometryFile)) OR $
        (NOT readAtmos(atmosFile)) ) THEN return

  state = gridWidgetSetup( )
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader

  displayGrid, widget_info( state.baseWidget, /CHILD )

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewGrid', state.baseWidget, $
   EVENT_HANDLER='XViewGrid_Event', GROUP_LEADER=group_leader

END
; -------- end ---------------------------- XViewGrid.pro ------------ ;
